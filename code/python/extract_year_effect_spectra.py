#!/usr/bin/env python
"""
extract_year_effect_spectra.py — pull AOP reflectance spectra at the
same set of CRBU points from BOTH 2018 and 2025 imagery, for direct
year-to-year comparison.

Companion to code/joint/12_year_effect_points.R, which produces
data/derived/year_effect_points.csv (vegetated + non-vegetated points,
stratified by snow-free DOY).

For each point:
  - Open the CRBU icechunk virtual Zarr store for {2018, 2025}.
  - Read a 3x3 m window of 1 m reflectance centered on (x_utm, y_utm).
  - Drop the no-data sentinel pixels, then take the per-band mean
    (no L2 normalization, no PCA). Preserves raw radiometric +
    phenology information so the user can compare any band/index
    directly.
  - Compute a quick NDVI for filtering / sanity-checking.

Output:
  data/derived/year_effect_spectra.parquet
    long-ish format, one row per (point_id, year)
    columns: point_id, year, x_utm, y_utm, point_type, cover_class,
             doy_band, snow_free_doy, n_pixels, ndvi, rfl_band_1..N
  data/derived/year_effect_wavelengths.csv
    band index -> wavelength (nm)

Required packages:
    icechunk xarray numpy pandas pyarrow

Usage:
    python extract_year_effect_spectra.py \\
        --points  data/derived/year_effect_points.csv \\
        --output  data/derived/year_effect_spectra.parquet \\
        --years   2018 2025 \\
        --domain  CRBU \\
        --max-points 50          # optional, for local smoke testing
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


logger = logging.getLogger("extract_year_effect_spectra")

NO_DATA_SENTINEL = -9000.0   # NEON AOP reflectance no-data
PIXEL_WINDOW_M   = 3         # 3x3 m window around each point


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


def extract_pixel_window(ds: xr.Dataset,
                         x_utm: float, y_utm: float,
                         halfwin_m: float = PIXEL_WINDOW_M / 2.0
                        ) -> tuple[np.ndarray | None, int]:
    """Return (mean spectrum, n_valid_pixels) for a 3 m window around
    (x_utm, y_utm). Mean is over the no-data-masked pixels of the
    window; no normalization."""
    cube = ds.reflectance.sel(
        northing=slice(y_utm + halfwin_m, y_utm - halfwin_m),   # descending
        easting=slice(x_utm - halfwin_m, x_utm + halfwin_m),
    ).values
    if cube.size == 0 or cube.shape[1] == 0 or cube.shape[2] == 0:
        return None, 0
    bands = cube.shape[0]
    flat = cube.reshape(bands, -1).T.astype(np.float64)
    flat = np.where(flat > NO_DATA_SENTINEL, flat, np.nan)
    valid = ~np.all(np.isnan(flat), axis=1)
    flat = flat[valid]
    if flat.shape[0] == 0:
        return None, 0
    return np.nanmean(flat, axis=0), int(flat.shape[0])


def quick_ndvi(spec: np.ndarray, wls: np.ndarray) -> float:
    """Cheap NDVI (NIR ~860 nm, RED ~660 nm) for sanity checks."""
    if spec is None:
        return np.nan
    nir = spec[int(np.argmin(np.abs(wls - 860)))]
    red = spec[int(np.argmin(np.abs(wls - 660)))]
    return (nir - red) / (nir + red) if (nir + red) > 0 else np.nan


def run(args: argparse.Namespace) -> int:
    points = pd.read_csv(args.points)
    if args.max_points:
        points = points.head(args.max_points).copy()
    logger.info("Points: %d (%s)", len(points), args.points)

    # Open both stores once.
    stores: dict[int, xr.Dataset] = {}
    wls_ref: np.ndarray | None = None
    for yr in args.years:
        logger.info("Opening %s %d virtual store ...", args.domain, yr)
        ds = open_aop_store(args.domain, yr)
        stores[yr] = ds
        wls_this = ds.wavelength.values
        if wls_ref is None:
            wls_ref = wls_this
        else:
            max_drift = float(np.max(np.abs(wls_this - wls_ref)))
            if max_drift > 5:
                logger.warning(
                    "Wavelength drift between years is %.2f nm — bands "
                    "are not perfectly aligned; will index by wavelength "
                    "in downstream comparison rather than by band index.",
                    max_drift,
                )

    n_bands = len(wls_ref)
    band_cols = [f"rfl_band_{i + 1}" for i in range(n_bands)]

    # Output wavelengths sidecar.
    wl_path = Path(args.output).with_name("year_effect_wavelengths.csv")
    pd.DataFrame({
        "band_number":   np.arange(1, n_bands + 1),
        "wavelength_nm": wls_ref,
    }).to_csv(wl_path, index=False)
    logger.info("Wrote %s (%d bands)", wl_path, n_bands)

    # Extract spectra.
    rows: list[dict] = []
    t_start = time.time()
    for yr, ds in stores.items():
        ok = fail = 0
        t_dom = time.time()
        for _, p in points.iterrows():
            try:
                spec, n_pix = extract_pixel_window(ds, p.x_utm, p.y_utm)
            except Exception as e:
                logger.exception("point %d year %d: %s", p.point_id, yr, e)
                spec, n_pix = None, 0
            if spec is None:
                fail += 1
                continue
            ok += 1
            row = {
                "point_id":      int(p.point_id),
                "year":          yr,
                "x_utm":         p.x_utm,
                "y_utm":         p.y_utm,
                "point_type":    p.point_type,
                "cover_class":   getattr(p, "cover_class", p.point_type),
                "doy_band":      p.doy_band,
                "snow_free_doy": p.snow_free_doy,
                "n_pixels":      n_pix,
                "ndvi":          quick_ndvi(spec, wls_ref),
            }
            row.update({c: float(v) for c, v in zip(band_cols, spec)})
            rows.append(row)
            if ok % 200 == 0:
                logger.info("  %d: %d ok / %d fail (%.0fs)",
                            yr, ok, fail, time.time() - t_dom)
        logger.info("Finished %d: %d ok / %d fail in %.0fs",
                    yr, ok, fail, time.time() - t_dom)

    df = pd.DataFrame(rows)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(args.output, index=False)
    logger.info("Wrote %s (%d rows, %d cols, %.1f min total)",
                args.output, len(df), df.shape[1],
                (time.time() - t_start) / 60)
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--points", type=Path,
                   default=Path("data/derived/year_effect_points.csv"))
    p.add_argument("--output", type=Path,
                   default=Path("data/derived/year_effect_spectra.parquet"))
    p.add_argument("--years", nargs="+", type=int, default=[2018, 2025])
    p.add_argument("--domain", default="CRBU")
    p.add_argument("--max-points", type=int, default=None,
                   help="Limit to first N points (local smoke testing).")
    p.add_argument("--log-level", default="INFO")
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
