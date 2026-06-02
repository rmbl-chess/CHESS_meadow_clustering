#!/usr/bin/env python
"""
mosaic_pc_maps.py — collapse per-tile PC COGs (output of
generate_aop_pc_maps.py) into a single multi-band COG per domain, so
the result can be moved off the server as one file per domain instead
of hundreds of tiles + a VRT.

Uses the per-domain VRT that generate_aop_pc_maps.py builds as the
input to gdal_translate — gdal_translate streams the VRT through the
COG driver, so memory use is bounded by the block size, not the full
raster size. Output is a strict COG tuned for GIS performance:
BAND-interleaved (single-PC views read one band, not all 20), ZSTD-
compressed (fast decode), with a deep internal overview pyramid.

If a per-domain VRT is missing (e.g., generate_aop_pc_maps.py was
interrupted before the final VRT pass), this script rebuilds it from
the per-tile COGs in <input-dir>/<DOMAIN>/pc_*.tif before mosaicking.

Usage:
    python mosaic_pc_maps.py \\
        --input-dir   data/derived/aop_pc_maps \\
        --output-dir  data/derived/aop_pc_maps_mosaic
        # --domains ALMO CRBU UPTA  (default: all VRTs in input-dir)

Required: gdal_translate on PATH (any modern GDAL install).
"""

from __future__ import annotations

import argparse
import logging
import subprocess
import sys
import time
from pathlib import Path


logger = logging.getLogger("mosaic_pc_maps")


def discover_domains(input_dir: Path,
                     requested: list[str] | None) -> list[str]:
    """Return the domain set to process. Look at existing VRTs and any
    per-tile subdirectories so we still find domains whose VRT hasn't
    been built yet."""
    domains: set[str] = set()
    for vrt in input_dir.glob("*_pc_mosaic.vrt"):
        domains.add(vrt.name.split("_pc_mosaic.vrt")[0])
    for sub in input_dir.iterdir():
        if sub.is_dir() and any(sub.glob("pc_*.tif")):
            domains.add(sub.name)
    if requested:
        domains = domains.intersection(requested)
    return sorted(domains)


def ensure_vrt(input_dir: Path, domain: str) -> Path:
    """Return the per-domain VRT path, building it from per-tile COGs if
    it doesn't already exist."""
    vrt_path = input_dir / f"{domain}_pc_mosaic.vrt"
    if vrt_path.exists():
        return vrt_path
    tile_dir = input_dir / domain
    tiles = sorted(tile_dir.glob("pc_*.tif"))
    if not tiles:
        raise FileNotFoundError(
            f"No per-tile COGs found under {tile_dir} — cannot build VRT."
        )
    logger.info("Building VRT for %s from %d tiles ...", domain, len(tiles))
    cmd = ["gdalbuildvrt", "-overwrite", str(vrt_path)] + [str(t) for t in tiles]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"gdalbuildvrt failed for {domain}: {res.stderr.strip()}")
    return vrt_path


def mosaic_domain(vrt: Path, out_tif: Path) -> tuple[bool, str]:
    """Convert one VRT mosaic into a single multi-band COG."""
    out_tif.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "gdal_translate",
        str(vrt), str(out_tif),
        "-of", "COG",
        # BAND interleave is the big GIS-performance win: a 20-band float
        # PC stack is almost always viewed one PC (or a 3-PC composite) at a
        # time. PIXEL interleave (the COG default) forces decompressing all
        # 20 bands per tile to render one; BAND lets GDAL touch only the
        # requested band(s) -> ~10-30x faster single-band reads in QGIS.
        "-co", "INTERLEAVE=BAND",
        # ZSTD decompresses far faster than DEFLATE (cost ~constant across
        # levels), so tile reads are quicker; LEVEL=13 keeps the ratio close
        # to DEFLATE-6 on this float data.
        "-co", "COMPRESS=ZSTD",
        "-co", "LEVEL=13",
        "-co", "PREDICTOR=YES",                 # floating-point predictor (=3)
        "-co", "BLOCKSIZE=512",
        # Force a fresh, deep overview pyramid so zoomed-out views cascade to
        # a tiny top level. IGNORE_EXISTING avoids reusing a shallow source
        # pyramid; COUNT=6 -> ~60 px top overview for these domains.
        "-co", "OVERVIEWS=IGNORE_EXISTING",
        "-co", "OVERVIEW_COUNT=6",
        "-co", "OVERVIEW_RESAMPLING=AVERAGE",
        "-co", "NUM_THREADS=ALL_CPUS",
        "-co", "BIGTIFF=IF_SAFER",
    ]
    logger.info("Running: %s", " ".join(cmd))
    t0 = time.time()
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        return False, res.stderr.strip()
    size_mb = out_tif.stat().st_size / 1e6
    return True, f"{size_mb:.1f} MB in {time.time() - t0:.1f}s"


def run(args: argparse.Namespace) -> int:
    input_dir  = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    if not input_dir.exists():
        logger.error("Input directory does not exist: %s", input_dir)
        return 1
    domains = discover_domains(input_dir, args.domains)
    if not domains:
        logger.error("No domains found under %s (looked for *_pc_mosaic.vrt "
                     "and <DOMAIN>/pc_*.tif)", input_dir)
        return 1
    logger.info("Domains to process: %s", domains)

    n_ok = n_fail = 0
    for domain in domains:
        out_tif = output_dir / f"{domain}_pc_mosaic.tif"
        if out_tif.exists() and not args.overwrite:
            logger.info("Skipping %s (exists)", out_tif.name)
            continue
        try:
            vrt = ensure_vrt(input_dir, domain)
        except (FileNotFoundError, RuntimeError) as e:
            logger.error("Cannot build VRT for %s: %s", domain, e)
            n_fail += 1
            continue
        ok, msg = mosaic_domain(vrt, out_tif)
        if ok:
            n_ok += 1
            logger.info("Wrote %s — %s", out_tif.name, msg)
        else:
            n_fail += 1
            logger.error("Failed %s: %s", out_tif.name, msg)

    logger.info("Done: %d ok, %d failed", n_ok, n_fail)
    return 0 if n_fail == 0 else 1


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input-dir",  type=Path,
                   default=Path("data/derived/aop_pc_maps"),
                   help="Directory containing {DOMAIN}_pc_mosaic.vrt files.")
    p.add_argument("--output-dir", type=Path,
                   default=Path("data/derived/aop_pc_maps_mosaic"),
                   help="Where to write the single-COG mosaics.")
    p.add_argument("--domains",    nargs="*", default=None,
                   help="Limit to these domain codes (default: all VRTs found).")
    p.add_argument("--overwrite",  action="store_true",
                   help="Re-mosaic even if the output COG already exists.")
    p.add_argument("--log-level",  default="INFO")
    return p


def main() -> int:
    args = build_parser().parse_args()
    logging.basicConfig(
        level=args.log_level,
        format="%(asctime)s  %(levelname)-7s  %(message)s",
        datefmt="%H:%M:%S",
    )
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
