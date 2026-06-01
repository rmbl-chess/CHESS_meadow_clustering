#!/usr/bin/env python
"""
generate_aop_pc_maps.py — write 20-band PC COGs for AOP tiles.

For each AOP tile (1 km x 1 km), compute PC01..PC20 using the field-
plot-derived PCA model exported by 24_export_classifier_model.R, and write
a 20-band Cloud-Optimized GeoTIFF (one band per PC) at 3 m resolution. The
preprocessing recipe matches the classifier feature pipeline exactly:

    1. Per-1 m-pixel L2 normalize each spectrum.
    2. 3x3 block mean -> 3 m resolution.
    2b. (--year 2018 only) Add the NDVI-stratified 2018->2025 radiometric
        correction in the L2-normalized space, mirroring R's
        04_join_spectra.R::apply_year_correction so 2018 PCs land on the
        2025-fit basis. Per 3 m cell: NDVI (660/860 nm) -> NDVI bin ->
        add that bin's per-band delta. NaN delta or NaN NDVI -> no change.
    3. Apply the water-band mask (drop ~1340-1450, ~1800-1950, >2400 nm).
    4. Subtract PCA center, project through saved PCA loadings -> PC01..PC20.

Year coverage: 2025 AOP exists for ALMO/CRBU/UPTA; 2018 AOP is CRBU-only
(matching the field campaign). The 2018->2025 correction is applied only
for --year 2018; pass --no-year-correction to disable it.

Output layout:
    data/derived/aop_pc_maps/
      ALMO/
        pc_337000_4286000.tif      # 20-band COG, 3 m, EPSG:32613
        pc_338000_4286000.tif
        ...
      CRBU/
        ...
      UPTA/
        ...
      ALMO_pc_mosaic.vrt           # per-domain mosaic for QGIS
      CRBU_pc_mosaic.vrt
      UPTA_pc_mosaic.vrt

Designed for the cloud server with AOP S3 data. For local testing, pass
--tile-limit N to process only the first N tiles. Skips tiles whose output
already exists, so re-running fills in gaps.

Tile inventory is read directly from S3 via boto3 (anonymous access).

Required packages:
    icechunk xarray rasterio numpy pandas boto3

Usage:
    python generate_aop_pc_maps.py \\
        --pca-model      data/derived/aop_classifier_pca.csv \\
        --meta           data/derived/aop_classifier_meta.json \\
        --output-dir     data/derived/aop_pc_maps \\
        --year           2025 \\
        --tile-limit     2          # for local testing

    # 2018 (CRBU only) — applies the NDVI-stratified 2018->2025 correction:
    python generate_aop_pc_maps.py --year 2018 --domains CRBU \\
        --output-dir data/derived/aop_pc_maps_2018
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import subprocess
import sys
import time
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import xarray as xr
from rasterio.transform import from_origin


logger = logging.getLogger("generate_aop_pc_maps")

NO_DATA_OUT = -9999.0
TILE_SIZE_M = 1000          # AOP mosaic tile size
AGG_FACTOR  = 3             # 1m -> 3m
N_PC        = 20

# NDVI-stratified 2018->2025 correction bins. Mirrors R 04_join_spectra.R /
# 15_stratified_year_correction.R: breaks c(-Inf,0.2,0.4,0.6,0.8,Inf), cut
# with right=FALSE. NDVI_BREAKS holds only the internal cut points so that
# np.digitize(ndvi, NDVI_BREAKS) (right=False) returns 0..4 in label order.
CORRECTION_YEAR = 2018
NDVI_BREAKS = (0.20, 0.40, 0.60, 0.80)
NDVI_LABELS = ("ndvi_lt_0.20", "ndvi_0.20_0.40", "ndvi_0.40_0.60",
               "ndvi_0.60_0.80", "ndvi_ge_0.80")
NDVI_RED_NM = 660.0
NDVI_NIR_NM = 860.0


# ---------- Config ---------------------------------------------------------

@dataclass(frozen=True)
class PreprocConfig:
    water_band_ranges_nm: list[tuple[float, float]]
    no_data_sentinel: float


@dataclass(frozen=True)
class YearCorrection:
    """NDVI-stratified 2018->2025 per-band delta, aligned to the AOP cube's
    full band order. Applied in L2-normalized space at the 3 m stage.

    delta:    (n_bins, n_full_bands) float32, NaN deltas already -> 0.
    red_idx / nir_idx: cube band indices nearest 660 / 860 nm (for NDVI).
    breaks:   internal NDVI cut points for np.digitize (right=False).
    """
    delta: np.ndarray
    red_idx: int
    nir_idx: int
    breaks: np.ndarray


def load_year_correction_delta(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read the NDVI-stratified delta CSV into (delta, wavelength_nm).

    delta is (n_bins, n_bands) ordered by NDVI_LABELS x ascending band_number,
    with NaN -> 0 to match R's `ifelse(is.na(delta), 0, delta)`.
    """
    cor = pd.read_csv(path)
    bands = np.sort(cor["band_number"].unique())
    delta = np.full((len(NDVI_LABELS), len(bands)), np.nan, dtype=np.float64)
    wl = np.full(len(bands), np.nan, dtype=np.float64)
    for i, lab in enumerate(NDVI_LABELS):
        sub = cor[cor["ndvi_bin"] == lab].set_index("band_number").reindex(bands)
        delta[i] = sub["delta"].to_numpy()
        wl = sub["wavelength_nm"].to_numpy()
    delta = np.nan_to_num(delta, nan=0.0).astype(np.float32)
    return delta, wl


# ---------- AOP store access -----------------------------------------------

def open_aop_store(domain: str, year: int) -> xr.Dataset:
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


def match_kept_bands(wls_full: np.ndarray, kept_wls: np.ndarray) -> np.ndarray:
    """Match each training wavelength to its nearest AOP-store band index."""
    return np.array(
        [int(np.argmin(np.abs(wls_full - w))) for w in kept_wls]
    )


def list_tiles_s3(domains: list[str], year: int) -> pd.DataFrame:
    """List AOP tile filenames per domain via boto3 (anonymous).

    Returns a DataFrame with columns (domain, tile, easting, northing).
    Matches the inventory R 22_target_pixels.R writes to aop_coverage.gpkg.
    boto3 is used rather than the AWS CLI so this works on servers without
    the CLI installed (it is pulled in by icechunk's S3 backend).
    """
    import boto3
    from botocore import UNSIGNED
    from botocore.client import Config

    s3 = boto3.client(
        "s3", region_name="us-east-2",
        config=Config(signature_version=UNSIGNED),
    )
    rows = []
    pat = re.compile(r"([^/]+_rfl_(\d+)_(\d+)\.nc)$")
    for dom in domains:
        prefix = f"AOP/spectrometer/mosaic/{dom}/{year}/"
        paginator = s3.get_paginator("list_objects_v2")
        for page in paginator.paginate(Bucket="rmbl-chess-data", Prefix=prefix):
            for obj in page.get("Contents", []):
                m = pat.search(obj["Key"])
                if m:
                    rows.append({
                        "domain":   dom,
                        "tile":     m.group(1),
                        "easting":  int(m.group(2)),
                        "northing": int(m.group(3)),
                    })
    return pd.DataFrame(rows)


# ---------- Tile processing -----------------------------------------------

def process_tile(
    ds: xr.Dataset,
    easting: int,
    northing: int,
    cfg: PreprocConfig,
    keep_band_idx: np.ndarray,
    pca_center: np.ndarray,
    pca_loadings: np.ndarray,
    out_path: Path,
    correction: YearCorrection | None = None,
) -> tuple[bool, str]:
    """Compute PC01..PC20 for one tile and write a 20-band COG.

    Returns (success, message). On success the message is a short summary;
    on failure it explains why (empty tile, no valid pixels, etc.).
    """
    # AOP store has descending northing — slice high->low.
    cube = ds.reflectance.sel(
        easting=slice(easting, easting + TILE_SIZE_M),
        northing=slice(northing + TILE_SIZE_M, northing),
    )
    if cube.size == 0:
        return False, "empty slice"

    east_coords  = cube.easting.values
    north_coords = cube.northing.values
    arr = cube.values.astype(np.float32)        # (bands, ny, nx)
    arr = np.where(arr > cfg.no_data_sentinel, arr, np.nan)

    bands, ny, nx = arr.shape

    ny_t = (ny // AGG_FACTOR) * AGG_FACTOR
    nx_t = (nx // AGG_FACTOR) * AGG_FACTOR
    if ny_t == 0 or nx_t == 0:
        return False, "tile smaller than aggregation factor"
    arr = arr[:, :ny_t, :nx_t]

    # --- 1. Per-1m-pixel L2 normalize --------------------------------------
    flat = arr.reshape(bands, -1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        norms = np.sqrt(np.nansum(flat ** 2, axis=0, keepdims=True))
    norms = np.where(norms == 0, np.nan, norms)
    arr_norm = (flat / norms).reshape(bands, ny_t, nx_t)
    del flat, norms

    # --- 2. 3x3 block mean -> 3m -------------------------------------------
    ny3 = ny_t // AGG_FACTOR
    nx3 = nx_t // AGG_FACTOR
    blk = arr_norm.reshape(bands, ny3, AGG_FACTOR, nx3, AGG_FACTOR)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        arr_3m = np.nanmean(blk, axis=(2, 4))
    del arr_norm, blk

    # --- 2b. NDVI-stratified 2018->2025 correction (year 2018 only) --------
    # Applied here, on the L2-normalized 3 m mean spectrum and BEFORE the
    # water mask, to mirror R 04_join_spectra.R (which corrects the
    # L2-normalized per-site mean over the full 426-band grid). NDVI picks
    # the bin; that bin's per-band delta is added. NaN NDVI -> no change.
    if correction is not None:
        if correction.delta.shape[1] != arr_3m.shape[0]:
            return False, (f"correction has {correction.delta.shape[1]} bands "
                           f"but tile has {arr_3m.shape[0]}")
        red = arr_3m[correction.red_idx]
        nir = arr_3m[correction.nir_idx]
        denom = nir + red
        with np.errstate(invalid="ignore", divide="ignore"):
            ndvi = np.where(denom != 0, (nir - red) / denom, np.nan)
        bin_idx = np.digitize(ndvi, correction.breaks)   # right=False -> 0..4
        valid = ~np.isnan(ndvi)
        for k in range(correction.delta.shape[0]):
            m = valid & (bin_idx == k)
            if m.any():
                arr_3m[:, m] += correction.delta[k][:, None]

    # --- 3. Select kept bands ----------------------------------------------
    arr_kept = arr_3m[keep_band_idx, :, :]
    del arr_3m

    # --- 4. PCA project -----------------------------------------------------
    n_kept = arr_kept.shape[0]
    arr_kept_flat = arr_kept.reshape(n_kept, ny3 * nx3)
    arr_centered  = arr_kept_flat - pca_center[:, None]
    pixel_valid   = ~np.isnan(arr_centered).any(axis=0)

    pcs = np.full((N_PC, ny3 * nx3), NO_DATA_OUT, dtype=np.float32)
    if pixel_valid.any():
        pcs[:, pixel_valid] = (
            pca_loadings.T @ arr_centered[:, pixel_valid]
        ).astype(np.float32)
    pcs_3d = pcs.reshape(N_PC, ny3, nx3)

    # GeoTIFF rows go north -> south.
    if north_coords[0] < north_coords[-1]:
        pcs_3d = pcs_3d[:, ::-1, :]
    if east_coords[0] > east_coords[-1]:
        pcs_3d = pcs_3d[:, :, ::-1]

    transform = from_origin(easting, northing + TILE_SIZE_M,
                            AGG_FACTOR, AGG_FACTOR)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(
        out_path, "w",
        driver="COG",
        height=ny3, width=nx3, count=N_PC,
        dtype="float32",
        crs="EPSG:32613",
        transform=transform,
        nodata=NO_DATA_OUT,
        COMPRESS="DEFLATE",
        LEVEL=6,
        PREDICTOR="YES",
        BLOCKSIZE=256,
        OVERVIEW_RESAMPLING="AVERAGE",
    ) as dst:
        for i in range(N_PC):
            dst.write(pcs_3d[i], i + 1)
            dst.set_band_description(i + 1, f"spec_PC{i + 1:02d}")

    n_valid = int(pixel_valid.sum())
    return True, f"{ny3}x{nx3}px, {n_valid}/{ny3 * nx3} valid"


# ---------- Main pipeline --------------------------------------------------

def run(args: argparse.Namespace) -> None:
    cfg_meta = json.loads(Path(args.meta).read_text())["preprocessing"]
    cfg = PreprocConfig(
        water_band_ranges_nm=[tuple(r) for r in cfg_meta["water_band_ranges_nm"]],
        no_data_sentinel=cfg_meta["no_data_sentinel"],
    )

    pca_df = pd.read_csv(args.pca_model)
    kept_wls     = pca_df["wavelength_nm"].to_numpy()
    pca_center   = pca_df["center"].to_numpy()
    pca_loadings = pca_df[[c for c in pca_df.columns
                           if c.startswith("PC")]].to_numpy()
    if pca_loadings.shape[1] < N_PC:
        raise RuntimeError(
            f"PCA model has {pca_loadings.shape[1]} PCs but N_PC={N_PC}"
        )
    pca_loadings = pca_loadings[:, :N_PC]
    logger.info("Loaded PCA model: %d retained bands x %d PCs",
                pca_loadings.shape[0], N_PC)

    # 2018->2025 radiometric correction (loaded here; wired to the cube's
    # band grid lazily on the first tile, alongside keep_band_idx).
    correction_delta = correction_wl = None
    if args.year == CORRECTION_YEAR and not args.no_year_correction:
        if not args.correction_csv.exists():
            raise FileNotFoundError(
                f"--year {CORRECTION_YEAR} needs {args.correction_csv} for the "
                "2018->2025 correction. Pass --no-year-correction to project "
                "2018 reflectance uncorrected (NOT recommended — PCs will be "
                "off-basis).")
        correction_delta, correction_wl = load_year_correction_delta(args.correction_csv)
        logger.info("Loaded 2018->2025 NDVI-stratified correction: "
                    "%d bins x %d bands (%s)",
                    correction_delta.shape[0], correction_delta.shape[1],
                    args.correction_csv.name)
    elif args.year == CORRECTION_YEAR:
        logger.warning("--no-year-correction set: 2018 PCs computed WITHOUT the "
                       "2018->2025 correction; they will not match the 2025 basis.")

    domains = args.domains or ["ALMO", "CRBU", "UPTA"]
    logger.info("Listing tiles from S3 for %s ...", ", ".join(domains))
    tiles = list_tiles_s3(domains, args.year)
    if args.tile_limit:
        tiles = tiles.head(args.tile_limit)
    logger.info("Tiles to process: %d (%s)",
                len(tiles),
                ", ".join(f"{d}={n}" for d, n in
                          tiles["domain"].value_counts().items()))

    out_root = Path(args.output_dir)
    out_root.mkdir(parents=True, exist_ok=True)

    keep_band_idx: np.ndarray | None = None
    correction: YearCorrection | None = None
    ds_cache: dict[str, xr.Dataset] = {}
    t_start = time.time()
    n_ok = n_skip = n_fail = 0

    for _, t in tiles.iterrows():
        out_path = out_root / t.domain / f"pc_{int(t.easting)}_{int(t.northing)}.tif"
        if out_path.exists() and not args.overwrite:
            n_skip += 1
            continue

        if t.domain not in ds_cache:
            logger.info("Opening %s virtual store", t.domain)
            ds_cache[t.domain] = open_aop_store(t.domain, args.year)
        ds = ds_cache[t.domain]

        if keep_band_idx is None:
            wls_full = ds.wavelength.values
            keep_band_idx = match_kept_bands(wls_full, kept_wls)
            max_drift = float(np.max(np.abs(wls_full[keep_band_idx] - kept_wls)))
            logger.info(
                "Matched %d training wavelengths (max drift %.2f nm)",
                len(keep_band_idx), max_drift,
            )
            if correction_delta is not None:
                if correction_delta.shape[1] != len(wls_full):
                    raise RuntimeError(
                        f"correction table has {correction_delta.shape[1]} bands "
                        f"but AOP cube has {len(wls_full)}; cannot align by index")
                wl_drift = float(np.max(np.abs(correction_wl - wls_full)))
                red_idx = int(np.argmin(np.abs(wls_full - NDVI_RED_NM)))
                nir_idx = int(np.argmin(np.abs(wls_full - NDVI_NIR_NM)))
                correction = YearCorrection(
                    delta=correction_delta, red_idx=red_idx, nir_idx=nir_idx,
                    breaks=np.asarray(NDVI_BREAKS),
                )
                logger.info(
                    "Year correction wired: NDVI red=%.1fnm(idx %d) "
                    "nir=%.1fnm(idx %d); max correction-vs-cube band drift %.2f nm",
                    wls_full[red_idx], red_idx, wls_full[nir_idx], nir_idx, wl_drift,
                )

        t0 = time.time()
        try:
            ok, msg = process_tile(
                ds, int(t.easting), int(t.northing),
                cfg, keep_band_idx, pca_center, pca_loadings, out_path,
                correction=correction,
            )
        except Exception as e:
            logger.exception("Tile (%s, %d, %d) crashed: %s",
                             t.domain, t.easting, t.northing, e)
            ok, msg = False, str(e)
        dt = time.time() - t0

        if ok:
            n_ok += 1
            logger.info("  %s pc_%d_%d.tif  %.1fs  %s",
                        t.domain, t.easting, t.northing, dt, msg)
        else:
            n_fail += 1
            logger.warning("  %s pc_%d_%d.tif  FAILED  %s",
                           t.domain, t.easting, t.northing, msg)

    logger.info(
        "Done: %d ok, %d skipped (existing), %d failed in %.1f min",
        n_ok, n_skip, n_fail, (time.time() - t_start) / 60,
    )

    if args.skip_vrt:
        return
    for domain in sorted(tiles["domain"].unique()):
        tifs = sorted((out_root / domain).glob("pc_*.tif"))
        if not tifs:
            continue
        vrt_path = out_root / f"{domain}_pc_mosaic.vrt"
        cmd = ["gdalbuildvrt", "-overwrite", str(vrt_path)] + \
              [str(p) for p in tifs]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            logger.warning("gdalbuildvrt failed for %s: %s", domain, res.stderr)
        else:
            logger.info("Built %s (%d tiles)", vrt_path.name, len(tifs))


# ---------- CLI ------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--pca-model",   type=Path,
                   default=Path("data/derived/aop_classifier_pca.csv"))
    p.add_argument("--meta",        type=Path,
                   default=Path("data/derived/aop_classifier_meta.json"))
    p.add_argument("--output-dir",  type=Path,
                   default=Path("data/derived/aop_pc_maps"))
    p.add_argument("--year",        type=int, default=2025)
    p.add_argument("--domains",     nargs="*", default=None,
                   help="Restrict to these domain codes (e.g. ALMO CRBU UPTA)")
    p.add_argument("--tile-limit",  type=int, default=None,
                   help="Process only the first N tiles (local testing).")
    p.add_argument("--overwrite",   action="store_true",
                   help="Re-process tiles even if outputs already exist.")
    p.add_argument("--skip-vrt",    action="store_true",
                   help="Skip building per-domain VRT mosaics at the end.")
    p.add_argument("--correction-csv", type=Path,
                   default=Path("data/small_reference/"
                                "year_effect_correction_2018_to_2025_by_ndvi.csv"),
                   help="NDVI-stratified 2018->2025 per-band delta table; "
                        "applied only when --year 2018.")
    p.add_argument("--no-year-correction", action="store_true",
                   help="Disable the 2018->2025 radiometric correction even "
                        "for --year 2018 (PCs will be off the 2025 basis).")
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
