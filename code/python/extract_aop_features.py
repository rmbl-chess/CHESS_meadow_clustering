#!/usr/bin/env python
"""
Extract AOP pixel features at target meadow pixels.

For each target pixel from data/derived/target_pixels.csv:
  1. Read a 3x3 window of 1 m AOP reflectance centered on (x_utm, y_utm)
     from the appropriate domain's virtual Zarr store on S3.
  2. L2-normalize each 1 m pixel's spectrum.
  3. Mean across valid pixels in the window -> one site-equivalent spectrum.
  4. Apply the water-band mask (drop 1340-1450, 1800-1950, >2400 nm).
  5. Subtract the PCA center and project through saved PCA loadings ->
     PC01..PC20.
  6. Compute the 6 narrow-band indices from the masked spectrum.
  7. Carry snow_free_doy and domain/tile through.

Output: a Parquet file with one row per pixel and these columns:
  x_utm, y_utm, domain, tile, snow_free_doy, doy_band,
  spec_PC01 .. spec_PC20, ndvi, ndwi, pri, red_edge_slope, cai, ndli

The output schema is the feature set expected by the R-trained classifier.

Designed to be run on a remote server with the AOP S3 data. Tested
locally on a small subset by passing --max-pixels.

Required packages:
    icechunk xarray pandas numpy pyarrow dask[distributed]

Usage:
    python extract_aop_features.py \\
        --targets        data/derived/target_pixels.csv \\
        --pca-model      data/derived/aop_classifier_pca.csv \\
        --meta           data/derived/aop_classifier_meta.json \\
        --output         data/derived/target_pixel_features.parquet \\
        --year           2025 \\
        --workers        8 \\
        --max-pixels     50      # for local testing
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import xarray as xr


logger = logging.getLogger("extract_aop_features")


# ---------- Config ---------------------------------------------------------

@dataclass(frozen=True)
class PreprocConfig:
    water_band_ranges_nm: list[tuple[float, float]]
    no_data_sentinel: float
    pixel_window_m: int   # half-window in metres; 1 means 3x3 1m pixels


# ---------- AOP store access -----------------------------------------------

def open_aop_store(domain: str, year: int) -> xr.Dataset:
    """Open the icechunk-virtual Zarr store for one (domain, year)."""
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


# ---------- Preprocessing helpers -----------------------------------------

def water_mask(wls_nm: np.ndarray, ranges: list[tuple[float, float]]) -> np.ndarray:
    """True for bands to KEEP (i.e., not in any water-absorption range)."""
    keep = np.ones_like(wls_nm, dtype=bool)
    for lo, hi in ranges:
        keep &= ~((wls_nm >= lo) & (wls_nm <= hi))
    return keep


def extract_pixel_spectrum(
    ds: xr.Dataset, x_utm: float, y_utm: float, cfg: PreprocConfig
) -> np.ndarray | None:
    """Extract the 3x3 window centered on (x, y), L2-normalize each 1 m
    pixel, and mean across the 9 pixels. Returns a (n_bands,) array or
    None if the window is entirely no-data."""
    half = cfg.pixel_window_m / 2.0
    # northing axis is DESCENDING in the AOP store, so slice high -> low
    cube = ds.reflectance.sel(
        northing=slice(y_utm + half, y_utm - half),
        easting=slice(x_utm - half, x_utm + half),
    ).values  # (bands, ~3, ~3)
    if cube.size == 0 or cube.shape[1] == 0 or cube.shape[2] == 0:
        return None
    bands = cube.shape[0]
    # reshape to (n_pixels, bands)
    flat = cube.reshape(bands, -1).T.astype(np.float64)
    flat = np.where(flat > cfg.no_data_sentinel, flat, np.nan)
    # Drop pixels that are all-NaN
    valid_pixel = ~np.all(np.isnan(flat), axis=1)
    flat = flat[valid_pixel]
    if flat.shape[0] == 0:
        return None
    # L2-normalize per pixel (compute norm over valid bands only)
    norms = np.sqrt(np.nansum(flat ** 2, axis=1, keepdims=True))
    norms = np.where(norms == 0, np.nan, norms)
    flat_norm = flat / norms
    # Mean across valid pixels (NaN-aware)
    return np.nanmean(flat_norm, axis=0)


def compute_indices(spec_kept: np.ndarray, kept_wls: np.ndarray) -> dict[str, float]:
    """Compute the 6 narrow-band indices from the masked, normalized
    spectrum (length = n_kept_bands)."""
    def b(target_nm: float) -> float:
        return float(spec_kept[int(np.argmin(np.abs(kept_wls - target_nm)))])

    nir, red = b(860), b(660)
    ndvi = (nir - red) / (nir + red) if (nir + red) != 0 else np.nan

    swir = b(1240)
    ndwi = (nir - swir) / (nir + swir) if (nir + swir) != 0 else np.nan

    b531, b570 = b(531), b(570)
    pri = (b531 - b570) / (b531 + b570) if (b531 + b570) != 0 else np.nan

    b700, b750 = b(700), b(750)
    wl_700 = float(kept_wls[int(np.argmin(np.abs(kept_wls - 700)))])
    wl_750 = float(kept_wls[int(np.argmin(np.abs(kept_wls - 750)))])
    red_edge_slope = (b750 - b700) / (wl_750 - wl_700)

    cai = 0.5 * (b(2000) + b(2200)) - b(2100)

    b1680, b1754 = b(1680), b(1754)
    if b1680 > 0 and b1754 > 0:
        l1, l2 = np.log(1 / b1754), np.log(1 / b1680)
        ndli = (l1 - l2) / (l1 + l2) if (l1 + l2) != 0 else np.nan
    else:
        ndli = np.nan

    return dict(ndvi=ndvi, ndwi=ndwi, pri=pri,
                red_edge_slope=red_edge_slope, cai=cai, ndli=ndli)


def project_pcs(spec_kept: np.ndarray,
                pca_center: np.ndarray,
                pca_loadings: np.ndarray) -> np.ndarray:
    """Return PC scores: centered spectrum @ loadings -> shape (n_pcs,)."""
    return (spec_kept - pca_center) @ pca_loadings


# ---------- Main pipeline --------------------------------------------------

def feature_row(
    ds: xr.Dataset,
    row: pd.Series,
    cfg: PreprocConfig,
    keep_band_idx: np.ndarray,
    kept_wls: np.ndarray,
    pca_center: np.ndarray,
    pca_loadings: np.ndarray,
) -> dict | None:
    spec_full = extract_pixel_spectrum(ds, row.x_utm, row.y_utm, cfg)
    if spec_full is None or np.isnan(spec_full).all():
        return None
    spec_kept = spec_full[keep_band_idx]
    # If any kept band failed, we still might be able to interpolate; here we
    # just drop pixels with NaN in kept bands.
    if np.isnan(spec_kept).any():
        return None
    indices = compute_indices(spec_kept, kept_wls)
    pcs = project_pcs(spec_kept, pca_center, pca_loadings)
    out = {
        "x_utm":         row.x_utm,
        "y_utm":         row.y_utm,
        "domain":        row.domain,
        "tile":          row.tile,
        "snow_free_doy": row.snow_free_doy,
        "doy_band":      row.doy_band,
    }
    for i, v in enumerate(pcs, start=1):
        out[f"spec_PC{i:02d}"] = float(v)
    out.update(indices)
    return out


def run(args: argparse.Namespace) -> None:
    cfg_meta = json.loads(Path(args.meta).read_text())["preprocessing"]
    cfg = PreprocConfig(
        water_band_ranges_nm=[tuple(r) for r in cfg_meta["water_band_ranges_nm"]],
        no_data_sentinel=cfg_meta["no_data_sentinel"],
        pixel_window_m=cfg_meta["pixel_window_m"],
    )

    pca_df = pd.read_csv(args.pca_model)
    kept_wls    = pca_df["wavelength_nm"].to_numpy()
    pca_center  = pca_df["center"].to_numpy()
    pca_loadings = pca_df[[c for c in pca_df.columns if c.startswith("PC")]].to_numpy()
    logger.info("Loaded PCA model: %d retained bands x %d PCs",
                pca_loadings.shape[0], pca_loadings.shape[1])

    targets = pd.read_csv(args.targets)
    if args.max_pixels:
        targets = targets.head(args.max_pixels).copy()
        logger.info("Restricted to first %d pixels for local testing",
                    len(targets))
    logger.info("Targets: %d pixels across %d domains",
                len(targets), targets["domain"].nunique())

    # Optional Dask client. For a sequential pixel loop it's marginal
    # because each .values triggers its own chunk fetch; the bigger win
    # is parallelism across domains/tiles. Left here for the server.
    client = None
    if args.workers > 1:
        from dask.distributed import Client
        client = Client(n_workers=args.workers,
                        threads_per_worker=2,
                        memory_limit=args.memory_limit)
        logger.info("Dask client: %s", client)

    # Pre-compute the keep_band_idx once we open the first store. Match each
    # training wavelength to the nearest AOP-store band by wavelength. The R
    # training pipeline labeled bands as whole-nm centers (384, 389, ...) but
    # the AOP store carries the NEON actual values (383.88, 388.89, ...).
    # Matching by index ensures we pull the same physical band.
    keep_band_idx: np.ndarray | None = None

    results: list[dict] = []
    t_start = time.time()

    for domain, group in targets.groupby("domain", sort=False):
        logger.info("Opening %s virtual store (%d pixels)", domain, len(group))
        ds = open_aop_store(domain, args.year)

        if keep_band_idx is None:
            wls_full = ds.wavelength.values
            # Match each training wavelength to its nearest AOP-store band by
            # wavelength. Robust against label drift (R training rounded to
            # whole nm; AOP store carries true NEON centers).
            keep_band_idx = np.array(
                [int(np.argmin(np.abs(wls_full - w))) for w in kept_wls]
            )
            max_drift_nm = float(np.max(np.abs(wls_full[keep_band_idx] - kept_wls)))
            logger.info(
                "Matched %d training wavelengths to AOP-store bands (max drift %.2f nm)",
                len(keep_band_idx), max_drift_nm,
            )
            if max_drift_nm > 3:
                logger.warning("Some wavelength matches drift > 3 nm; check label compatibility.")

        ok = fail = 0
        t_dom = time.time()
        for _, row in group.iterrows():
            try:
                feat = feature_row(ds, row, cfg, keep_band_idx,
                                   kept_wls, pca_center, pca_loadings)
            except Exception as e:
                logger.exception("Failed pixel (%s, %s): %s",
                                 row.x_utm, row.y_utm, e)
                feat = None
            if feat is None:
                fail += 1
                continue
            ok += 1
            results.append(feat)
            if ok % 200 == 0:
                logger.info("  %s: %d ok / %d fail (%.1fs)",
                            domain, ok, fail, time.time() - t_dom)

        logger.info("Finished %s: %d ok / %d fail in %.1fs",
                    domain, ok, fail, time.time() - t_dom)

    elapsed = time.time() - t_start
    logger.info("Total: %d feature rows in %.1fs", len(results), elapsed)

    out_df = pd.DataFrame(results)
    out_df.to_parquet(args.output, index=False)
    logger.info("Wrote %s (%d rows, %d cols)",
                args.output, len(out_df), out_df.shape[1])

    if client is not None:
        client.close()


# ---------- CLI ------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--targets",     type=Path,
                   default=Path("data/derived/target_pixels.csv"))
    p.add_argument("--pca-model",   type=Path,
                   default=Path("data/derived/aop_classifier_pca.csv"))
    p.add_argument("--meta",        type=Path,
                   default=Path("data/derived/aop_classifier_meta.json"))
    p.add_argument("--output",      type=Path,
                   default=Path("data/derived/target_pixel_features.parquet"))
    p.add_argument("--year",        type=int, default=2025)
    p.add_argument("--workers",     type=int, default=1,
                   help="Dask workers; 1 = no Dask client.")
    p.add_argument("--memory-limit", default="4GB")
    p.add_argument("--max-pixels",  type=int, default=None,
                   help="Limit to first N pixels (local testing).")
    p.add_argument("--log-level",   default="INFO")
    return p


def main() -> int:
    args = build_parser().parse_args()
    logging.basicConfig(
        level=args.log_level,
        format="%(asctime)s  %(levelname)-7s  %(message)s",
        datefmt="%H:%M:%S",
    )
    run(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
