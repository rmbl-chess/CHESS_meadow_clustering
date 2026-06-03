#!/usr/bin/env python
"""
natureserve_fetch.py — gather candidate NatureServe IVC communities for the
meadow/shrub classes and cache them to JSON for offline crosswalk matching.

Strategy (see the crosswalk discussion in docs / commit history): our classes
are characterized by an IndVal diagnostic-species assemblage. We search
NatureServe Explorer for each reconciled diagnostic species, collect the
ECOSYSTEM (IVC vegetation-unit) hits, fetch each unique community's detail
record, and cache the fields needed for matching + enrichment. Matching itself
runs offline against this cache (code/meadow/20_ecosystem_crosswalk.R), so the
pipeline has no live API dependency and the snapshot is reproducible + citable.

Why this shape:
  - Diagnostic species drive the search (a community is named by its species).
  - communityCompositions is empty at Group level, so the diagnostic species
    are parsed from each community's scientific NAME (the IVC name encodes them).
  - Geography comes from elementNationals -> subnation codes (CO flag).
  - Names are reconciled first via species_natureserve_crosswalk.csv so no class
    silently searches a name NatureServe doesn't index (e.g. Veratrum
    tenuipetalum -> V. californicum).

API (confirmed empirically; the published docs mis-state the paths):
  POST /api/data/search                      body: quickSearch token, paging
  GET  /api/data/taxon/ELEMENT_GLOBAL.2.{id} detail record
No auth. Be respectful: User-Agent + --sleep between calls. Idempotent: cached
detail records are skipped unless --refresh.

Inputs:
  data/small_reference/species_natureserve_crosswalk.csv   (name reconciliation)
  data/derived/label_descriptions.csv                      (IndVal assemblages)
Output:
  data/derived/natureserve_cache.json   (fetch date, per-species hits, per-
                                          community detail; force-include in git)

Usage:
  python code/python/natureserve_fetch.py --sleep 0.2
  python code/python/natureserve_fetch.py --max-species 6   # quick test run
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import re
import time
import urllib.request
import urllib.error
from datetime import datetime, timezone
from pathlib import Path

logger = logging.getLogger("natureserve_fetch")

API = "https://explorer.natureserve.org/api/data"
UA = "CHESS-meadow-crosswalk/0.1 (RMBL; respectful cached pull)"
NONSP_PREFIX = ("Other", "NPV", "Bare")
# species binomial in an IVC community name: "Genus epithet" (epithet lowercase)
SPECIES_RE = re.compile(r"\b([A-Z][a-z]+ [a-z][a-z-]+)\b")


# ---------- API client -----------------------------------------------------

def _request(url: str, body: dict | None = None, retries: int = 3,
             sleep: float = 0.2):
    data = json.dumps(body).encode() if body is not None else None
    headers = {"Accept": "application/json", "User-Agent": UA}
    if data is not None:
        headers["Content-Type"] = "application/json"
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, data=data, headers=headers)
            with urllib.request.urlopen(req, timeout=45) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503) and attempt < retries - 1:
                time.sleep(2 ** attempt)         # backoff on transient errors
                continue
            logger.warning("HTTP %s for %s", e.code, url)
            return None
        except Exception as e:                   # noqa: BLE001 - network flake
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
                continue
            logger.warning("request failed (%s) for %s", e, url)
            return None
    return None


def search_ecosystems(token: str, n: int, sleep: float) -> list[dict]:
    """Return ECOSYSTEM (IVC vegetation-unit) hits for a quick-search token."""
    resp = _request(f"{API}/search", {
        "criteriaType": "combined",
        "textCriteria": [{"paramType": "quickSearch", "searchToken": token}],
        "pagingOptions": {"page": 0, "recordsPerPage": n},
    }, sleep=sleep)
    if not resp:
        return []
    out = []
    for r in resp.get("results", []):
        if r.get("recordType") == "ECOSYSTEM":
            out.append({"elementGlobalId": r.get("elementGlobalId"),
                        "scientificName": r.get("scientificName"),
                        "roundedGRank": r.get("roundedGRank"),
                        "nsxUrl": r.get("nsxUrl")})
    return out


def fetch_detail(eid, sleep: float) -> dict | None:
    d = _request(f"{API}/taxon/ELEMENT_GLOBAL.2.{eid}", sleep=sleep)
    if not d:
        return None
    # geography: subnation codes across all nations
    subnats = []
    for en in d.get("elementNationals") or []:
        iso = (en.get("nation") or {}).get("isoCode")
        for sn in en.get("elementSubnationals") or []:
            code = (sn.get("subnation") or {}).get("subnationCode")
            if code:
                subnats.append(f"{iso}-{code}" if iso else code)
    eg = d.get("ecosystemGlobal") or {}
    name = d.get("scientificName") or ""
    return {
        "elementGlobalId": d.get("elementGlobalId"),
        "scientificName": name,
        "primaryCommonName": d.get("primaryCommonName"),
        "elcode": d.get("elcode"),
        "roundedGRank": d.get("roundedGRank") or d.get("grank"),
        "classificationLevel": (d.get("classificationLevel") or {}).get(
            "classificationLevelNameEn"),
        "nsxUrl": d.get("nsxUrl"),
        "subnations": sorted(set(subnats)),
        "in_colorado": any(s.endswith("-CO") or s == "CO" for s in subnats),
        # diagnostic species parsed from the IVC community NAME (compositions
        # are empty at Group level; the name encodes the diagnostics).
        "name_species": sorted(set(SPECIES_RE.findall(name))),
        "conceptSentence": eg.get("conceptSentence"),
        "summary": eg.get("summary"),
        "macrogroup": ((eg.get("macrogroupHierarchy") or {})
                       .get("scientificName")),
    }


# ---------- search-token derivation ----------------------------------------

def load_crosswalk(path: Path) -> dict:
    m = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            m[row["our_name"]] = (row["natureserve_name"], row["match_level"])
    return m


def species_from_descriptions(path: Path) -> list[str]:
    """Unique binomials across the top-5 IndVal indicator + abundant columns."""
    names = set()
    with open(path) as f:
        for row in csv.DictReader(f):
            for col in ("indicators", "abundant"):
                for tok in (row.get(col) or "").split(";"):
                    nm = re.sub(r"\s*\(.*$", "", tok).strip()
                    if nm and not nm.startswith(NONSP_PREFIX) and " " in nm:
                        names.add(nm)
    return sorted(names)


def reconciled_tokens(species: list[str], crosswalk: dict) -> dict:
    """our_name -> search token (reconciled; identity if not in crosswalk)."""
    out = {}
    for s in species:
        if s in crosswalk:
            out[s] = crosswalk[s][0]          # natureserve_name (species|genus)
        else:
            out[s] = s
    return out


# ---------- main -----------------------------------------------------------

def run(args: argparse.Namespace) -> int:
    crosswalk = load_crosswalk(args.crosswalk)
    species = species_from_descriptions(args.species_source)
    tokens = reconciled_tokens(species, crosswalk)
    uniq_tokens = sorted(set(tokens.values()))
    if args.max_species:
        uniq_tokens = uniq_tokens[: args.max_species]
    logger.info("%d class species -> %d unique reconciled search tokens",
                len(species), len(uniq_tokens))

    out_path = Path(args.out)
    cache = {}
    if out_path.exists() and not args.refresh:
        cache = json.loads(out_path.read_text())
        logger.info("resuming from cache: %d communities already fetched",
                    len(cache.get("ecosystems", {})))
    ecosystems = cache.get("ecosystems", {})
    species_hits = cache.get("species_search", {})

    # 1. search each token, collect ecosystem ids
    for i, tok in enumerate(uniq_tokens, 1):
        if tok in species_hits and not args.refresh:
            continue
        hits = search_ecosystems(tok, args.records_per_species, args.sleep)
        species_hits[tok] = [h["elementGlobalId"] for h in hits]
        logger.info("[%d/%d] '%s' -> %d ecosystem hits",
                    i, len(uniq_tokens), tok, len(hits))
        time.sleep(args.sleep)

    # 2. fetch detail for each unique ecosystem id (skip already cached)
    all_ids = sorted({eid for ids in species_hits.values() for eid in ids
                      if eid is not None})
    todo = [e for e in all_ids if str(e) not in ecosystems]
    logger.info("%d unique communities (%d new to fetch)",
                len(all_ids), len(todo))
    for i, eid in enumerate(todo, 1):
        det = fetch_detail(eid, args.sleep)
        if det:
            ecosystems[str(eid)] = det
            if i % 20 == 0 or i == len(todo):
                logger.info("  fetched %d/%d details", i, len(todo))
        time.sleep(args.sleep)

    # 3. write cache
    out = {
        "fetch_date": args.today or datetime.now(timezone.utc).strftime("%Y-%m-%d"),
        "api_base": API,
        "citation": "NatureServe. NatureServe Explorer. https://explorer.natureserve.org",
        "n_search_tokens": len(uniq_tokens),
        "species_search": species_hits,
        "ecosystems": ecosystems,
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=1))
    co = sum(1 for e in ecosystems.values() if e.get("in_colorado"))
    logger.info("Wrote %s: %d communities (%d in Colorado)",
                out_path, len(ecosystems), co)
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--crosswalk", type=Path,
                   default=Path("data/small_reference/species_natureserve_crosswalk.csv"))
    p.add_argument("--species-source", type=Path,
                   default=Path("data/derived/label_descriptions.csv"))
    p.add_argument("--out", type=Path,
                   default=Path("data/derived/natureserve_cache.json"))
    p.add_argument("--records-per-species", type=int, default=25,
                   help="max search hits to keep per species token")
    p.add_argument("--max-species", type=int, default=None,
                   help="cap number of search tokens (quick test runs)")
    p.add_argument("--sleep", type=float, default=0.2,
                   help="seconds between API calls (be respectful)")
    p.add_argument("--refresh", action="store_true",
                   help="re-fetch even if already cached")
    p.add_argument("--today", default=None,
                   help="override fetch_date (YYYY-MM-DD) for reproducible runs")
    p.add_argument("--log-level", default="INFO")
    return p


def main() -> int:
    args = build_parser().parse_args()
    logging.basicConfig(level=args.log_level,
                        format="%(asctime)s  %(levelname)-7s  %(message)s",
                        datefmt="%H:%M:%S")
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
