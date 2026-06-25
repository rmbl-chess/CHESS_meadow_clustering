#!/usr/bin/env python
"""
extract_supplemental_spectra.py — extract per-pixel NEON AOP reflectance at the
2026 supplemental field-crown polygons, producing a CSV that matches the
delivered `site_extraction_spectra_2025 (1).csv` schema so the R meadow/shrub
pipeline can read it with no format changes.

Companion to data/raw/Supplemental_field_2026/ (69 plots, sites 2001-2069).
The polygons carry only a site label; this script derives `domain` (by
centroid-in-bbox) and `site_type` (Shrub iff a Named Species covers 100% of a
known shrub genus, else Meadow) in-script, so it is fully self-contained for
the Hub. Spectra come from the 2025 AOP (most recent flight); no radiometric
year-correction is applied (the 2025 basis is the reference) — absolute scale
is irrelevant anyway because 04_join_spectra.R L2-normalizes every pixel.

For each polygon:
  - Open the {ALMO,CRBU} 2025 icechunk virtual Zarr store.
  - Window the reflectance cube to the polygon bbox.
  - Keep every 1 m pixel whose CENTER falls inside the polygon (shapely
    contains); if a small polygon contains no center, fall back to the single
    nearest pixel to the centroid.
  - Drop the -9000 no-data sentinel pixels. Emit ONE ROW PER PIXEL.

Output columns (862, identical order to the 2025 delivery):
  site_number, domain, sampling_area, site_type, fid, row, col, x_utm, y_utm,
  shade, rfl_band_1..426, unc_band_1..426
  - shade is left blank (NEON's shade mask isn't reproducible here; the R
    aggregation keeps shade==1 OR NA, so blank == kept).
  - unc_band_* are left blank (carried for schema parity; never clustered on).

Required packages: icechunk xarray numpy pandas shapely  (NOT geopandas).

Usage (Hub, conda chess-hub):
    python code/python/extract_supplemental_spectra.py \\
        --polygons data/raw/Supplemental_field_2026/augment_polygons_2026_06_23_wgs_utm.geojson \\
        --cover    data/raw/Supplemental_field_2026/augment_cover_cleaned_2026_06_23.csv \\
        --year     2025 \\
        --output   "data/raw/ESS-DIVE-Spectra/site_extraction_spectra_2026.csv"
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from shapely.geometry import shape, Point

logger = logging.getLogger("extract_supplemental_spectra")

NO_DATA_SENTINEL = -9000.0
N_BANDS = 426

# Domain extents (EPSG:32613) from the 3 m CHM rasters. The three boxes are
# effectively disjoint (no x/y overlap between ALMO/CRBU/UPTA), so centroid-in-
# bbox assigns cleanly; ALMO is tested first as a tie-break safeguard.
DOMAIN_BBOX = {
    "ALMO": (337000.0, 356002.0, 4273998.0, 4299000.0),   # xmin,xmax,ymin,ymax
    "CRBU": (315000.0, 338001.0, 4297000.0, 4327000.0),
    "UPTA": (343000.0, 363001.0, 4302000.0, 4326000.0),
}
DOMAIN_ORDER = ("ALMO", "CRBU", "UPTA")

# Shrub genera (mirror code/shrub/01_load.R) — used only to label site_type.
SHRUB_GENERA = {
    "Salix", "Betula", "Alnus", "Cornus", "Lonicera", "Ribes", "Symphoricarpos",
    "Amelanchier", "Prunus", "Acer", "Artemisia", "Purshia", "Chrysothamnus",
    "Ericameria", "Dasiphora", "Potentilla", "Rosa", "Sambucus", "Shepherdia",
    "Juniperus", "Arctostaphylos", "Vaccinium", "Spiraea", "Rubus",
}


def open_aop_store(domain: str, year: int) -> xr.Dataset:
    """Open the icechunk virtual Zarr store for one (domain, year)."""
    import icechunk
    storage = icechunk.s3_storage(
        bucket="rmbl-chess-data",
        prefix=f"virtual/AOP/spectrometer/{domain}/{year}/",
        region="us-east-2",
        anonymous=True,
    )
    repo = icechunk.Repository.open(
        storage,
        authorize_virtual_chunk_access={"s3://rmbl-chess-data/": None},
    )
    return xr.open_zarr(repo.readonly_session(branch="main").store)


def assign_domain(cx: float, cy: float) -> str | None:
    """Domain whose bbox contains the polygon centroid (ALMO first)."""
    for dom in DOMAIN_ORDER:
        xmin, xmax, ymin, ymax = DOMAIN_BBOX[dom]
        if xmin <= cx <= xmax and ymin <= cy <= ymax:
            return dom
    # Fallback: nearest bbox center.
    best, best_d = None, np.inf
    for dom, (xmin, xmax, ymin, ymax) in DOMAIN_BBOX.items():
        d = np.hypot(cx - (xmin + xmax) / 2, cy - (ymin + ymax) / 2)
        if d < best_d:
            best, best_d = dom, d
    logger.warning("Centroid (%.1f, %.1f) in no bbox; nearest=%s", cx, cy, best)
    return best


def derive_site_types(cover: pd.DataFrame) -> dict[int, str]:
    """Shrub iff a Named Species covers 100% of a shrub genus, else Meadow."""
    named = cover[cover["Cover_Type"] == "Live Vegetation - Named Species"].copy()
    named["genus"] = (named["Cover_Class_Name"].fillna("")
                      .str.strip().str.split().str[0])
    is_shrub100 = ((named["Cover_Percent"] >= 100)
                   & named["genus"].isin(SHRUB_GENERA))
    shrub_sites = set(named.loc[is_shrub100, "Site_Number"].astype(int))
    return {int(s): ("Shrub" if int(s) in shrub_sites else "Meadow")
            for s in cover["Site_Number"].unique()}


def polygon_pixels(ds: xr.Dataset, geom) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (spectra[n_pix, bands], x_centers[n_pix], y_centers[n_pix]) for
    the 1 m pixels whose centers fall inside `geom`. Falls back to the nearest
    single pixel to the centroid when none are contained."""
    minx, miny, maxx, maxy = geom.bounds
    pad = 1.5
    cube = ds.reflectance.sel(
        northing=slice(maxy + pad, miny - pad),   # northing descends
        easting=slice(minx - pad, maxx + pad),
    )
    east = cube.easting.values
    north = cube.northing.values
    if east.size == 0 or north.size == 0:
        return np.empty((0, N_BANDS)), np.empty(0), np.empty(0)

    vals = cube.values                            # (bands, n_north, n_east)
    ex, ny = np.meshgrid(east, north)             # (n_north, n_east)
    flat_x = ex.ravel()
    flat_y = ny.ravel()
    spec = vals.reshape(vals.shape[0], -1).T      # (n_north*n_east, bands)

    inside = np.array([geom.contains(Point(x, y))
                       for x, y in zip(flat_x, flat_y)], dtype=bool)
    if not inside.any():
        c = geom.centroid
        j = int(np.argmin(np.hypot(flat_x - c.x, flat_y - c.y)))
        inside = np.zeros(flat_x.size, dtype=bool)
        inside[j] = True

    spec = spec[inside].astype(np.float64)
    spec = np.where(spec > NO_DATA_SENTINEL, spec, np.nan)
    valid = ~np.all(np.isnan(spec), axis=1)
    return spec[valid], flat_x[inside][valid], flat_y[inside][valid]


def quick_ndvi(spec: np.ndarray, wls: np.ndarray) -> float:
    nir = spec[int(np.argmin(np.abs(wls - 860)))]
    red = spec[int(np.argmin(np.abs(wls - 660)))]
    return (nir - red) / (nir + red) if (nir + red) > 0 else np.nan


def run(args: argparse.Namespace) -> int:
    cover = pd.read_csv(args.cover)
    cover["Site_Number"] = cover["Site_Number"].astype(int)
    site_type = derive_site_types(cover)

    gj = json.load(open(args.polygons))
    polys: dict[int, object] = {}
    for f in gj["features"]:
        sid = int(f["properties"]["Label_of_F"])
        polys[sid] = shape(f["geometry"])
    logger.info("Polygons: %d sites", len(polys))

    # Group sites by domain (centroid-in-bbox).
    by_domain: dict[str, list[int]] = {d: [] for d in DOMAIN_BBOX}
    for sid, geom in polys.items():
        c = geom.centroid
        by_domain[assign_domain(c.x, c.y)].append(sid)
    for d, sids in by_domain.items():
        logger.info("  domain %s: %d sites", d, len(sids))

    band_cols = [f"rfl_band_{i+1}" for i in range(N_BANDS)]
    unc_cols = [f"unc_band_{i+1}" for i in range(N_BANDS)]
    meta_cols = ["site_number", "domain", "sampling_area", "site_type",
                 "fid", "row", "col", "x_utm", "y_utm", "shade"]

    rows: list[dict] = []
    wls_ref: np.ndarray | None = None
    t_start = time.time()
    for dom, sids in by_domain.items():
        if not sids:
            continue
        logger.info("Opening %s %d store ...", dom, args.year)
        ds = open_aop_store(dom, args.year)
        wls = ds.wavelength.values
        if wls_ref is None:
            wls_ref = wls
        ok = empty = 0
        for sid in sids:
            spec, xs, ys = polygon_pixels(ds, polys[sid])
            if spec.shape[0] == 0:
                empty += 1
                logger.warning("site %d (%s): no valid pixels", sid, dom)
                continue
            ok += 1
            for k in range(spec.shape[0]):
                row = {
                    "site_number": sid, "domain": dom, "sampling_area": "",
                    "site_type": site_type.get(sid, "Meadow"),
                    "fid": f"AUG{sid}", "row": "", "col": "",
                    "x_utm": float(xs[k]), "y_utm": float(ys[k]), "shade": "",
                }
                row.update({c: float(v) for c, v in zip(band_cols, spec[k])})
                row.update({c: "" for c in unc_cols})
                rows.append(row)
            if ok % 20 == 0:
                logger.info("  %s: %d sites done (%.0fs)", dom, ok,
                            time.time() - t_start)
        logger.info("Finished %s: %d sites with pixels, %d empty", dom, ok, empty)

    df = pd.DataFrame(rows, columns=meta_cols + band_cols + unc_cols)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)
    logger.info("Wrote %s (%d pixel rows across %d sites, %.1f min)",
                args.output, len(df), df["site_number"].nunique(),
                (time.time() - t_start) / 60)

    # Per-site pixel count + NDVI sanity sidecar.
    if wls_ref is not None and len(df):
        bmat = df[band_cols].to_numpy(dtype=float)
        ndvi = [quick_ndvi(bmat[i], wls_ref) for i in range(len(df))]
        summ = (df.assign(ndvi=ndvi)
                  .groupby(["site_number", "domain", "site_type"])
                  .agg(n_pixels=("x_utm", "size"), ndvi=("ndvi", "mean"))
                  .reset_index())
        summ_path = Path(args.output).with_name("site_extraction_spectra_2026_summary.csv")
        summ.to_csv(summ_path, index=False)
        logger.info("Wrote %s (per-site pixel counts + mean NDVI)", summ_path)
        logger.info("Per-site pixels: min=%d median=%d max=%d",
                    summ.n_pixels.min(), int(summ.n_pixels.median()),
                    summ.n_pixels.max())
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--polygons", type=Path, required=True)
    p.add_argument("--cover", type=Path, required=True)
    p.add_argument("--year", type=int, default=2025)
    p.add_argument("--output", type=Path,
                   default=Path("data/raw/ESS-DIVE-Spectra/site_extraction_spectra_2026.csv"))
    p.add_argument("--log-level", default="INFO")
    return p


def main() -> int:
    args = build_parser().parse_args()
    logging.basicConfig(
        level=args.log_level,
        format="%(asctime)s  %(levelname)-7s  %(message)s", datefmt="%H:%M:%S")
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
