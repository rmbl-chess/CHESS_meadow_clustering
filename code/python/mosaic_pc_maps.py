#!/usr/bin/env python
"""
mosaic_pc_maps.py — collapse per-tile PC COGs (output of
generate_aop_pc_maps.py) into a single multi-band COG per domain, so
the result can be moved off the server as one file per domain instead
of hundreds of tiles + a VRT.

Uses the per-domain VRT that generate_aop_pc_maps.py already builds as
the input to gdal_translate — gdal_translate streams the VRT through
the COG driver, so memory use is bounded by the block size, not the
full raster size. Output is a strict COG with internal overviews.

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


def mosaic_domain(vrt: Path, out_tif: Path) -> tuple[bool, str]:
    """Convert one VRT mosaic into a single multi-band COG."""
    out_tif.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "gdal_translate",
        str(vrt), str(out_tif),
        "-of", "COG",
        "-co", "COMPRESS=DEFLATE",
        "-co", "LEVEL=6",
        "-co", "PREDICTOR=YES",
        "-co", "BLOCKSIZE=512",
        "-co", "OVERVIEW_RESAMPLING=AVERAGE",
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
    vrts = sorted(input_dir.glob("*_pc_mosaic.vrt"))
    if args.domains:
        vrts = [v for v in vrts
                if v.name.split("_pc_mosaic.vrt")[0] in args.domains]
    if not vrts:
        logger.error("No *_pc_mosaic.vrt files found under %s", input_dir)
        return 1

    logger.info("Found %d VRT mosaic(s): %s",
                len(vrts), [v.name for v in vrts])

    n_ok = n_fail = 0
    for vrt in vrts:
        domain  = vrt.name.split("_pc_mosaic.vrt")[0]
        out_tif = output_dir / f"{domain}_pc_mosaic.tif"
        if out_tif.exists() and not args.overwrite:
            logger.info("Skipping %s (exists)", out_tif.name)
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
