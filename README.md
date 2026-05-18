# CHESS_meadow_clustering

Cluster co-occurring plant species and NEON AOP spectra into training samples for a CHESS meadow map classifier.

## What this is

The CHESS map classifier needs training samples that are (a) spectrally distinguishable in NEON AOP imagery and (b) interpretable as known plant species compositions. This repo joins field-collected vegetation data from ESS-DIVE with NEON AOP spectra, then clusters the combined feature space to define those training samples. Audience is the CHESS team.

## Quick start

```bash
git clone <repo-url> CHESS_meadow_clustering
cd CHESS_meadow_clustering
R -e 'renv::restore()'
R -e 'source("code/01_load.R")'
```

## Layout

- `code/` — numbered R scripts (`01_load.R` → `04_figures.R`).
- `data/` — inputs; see `data/README.md` for provenance.
- `output/` — figures and tables (gitignored).
- `docs/` — methods notes; canonical figures.

## Data

- **CHESS 2025 field vegetation** — ESS-DIVE: species list, meadow cover, site metadata (in `data/ESS-DIVE-Vegetation-Field/`).
- **NEON AOP spectra** — to be added; document flight year, tile IDs, and QA in `data/README.md`.
- CRS for SDP-derived layers: `EPSG:32613` (UTM Zone 13N). Field data needs reprojection.

## Tests

No tests; this is an exploratory analysis pipeline.

## License

Internal-only; no LICENSE file.
