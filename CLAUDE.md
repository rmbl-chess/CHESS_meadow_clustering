# CHESS_meadow_clustering

Defines clusters of co-occurring plant species and NEON AOP spectra to seed training samples for a meadow map classifier in the CHESS project area. Goal: clusters that are spectrally distinct (and therefore mappable from AOP imagery) but also interpretable as recognizable plant species compositions. Audience is the CHESS team; downstream consumer is a supervised classifier over NEON AOP tiles. SDP CRS is `EPSG:32613`.

## Data pipeline

Two field sampling campaigns (2018 and 2025) are reconciled into a single combined cover table, then joined with NEON AOP spectra at matching plot footprints:

1. **Load** both campaigns' plot-level data (sites, cover, species lists).
2. **Reconcile taxonomy** across the two species-list tables. Most species names match; build a crosswalk for the ones that don't.
3. **Combine cover** into a single table keyed by `(SiteID, Year)`, with one column per harmonized species (`<Spp>_cover`).
4. **Join spectra** by site ID to produce the final vegetation–spectrum dataset that feeds clustering.

## Layout

- `code/` — numbered R scripts driving the pipeline.
  - `01_load.R` — read 2018 + 2025 field data.
  - `02_reconcile_taxonomy.R` — crosswalk species names between campaigns.
  - `03_combine_cover.R` — assemble the `(SiteID, Year, Spp*_cover)` table.
  - `04_join_spectra.R` — join NEON AOP spectra to the combined cover table.
  - Clustering / figure scripts will be added once the analysis approach is chosen.
- `data/` — inputs, with provenance in `data/README.md`. Three source dirs under `data/raw/` (gitignored — fetch from ESS-DIVE):
  - `data/raw/ESS-DIVE-Vegetation-Field-2018/` — 2018 fractional cover, species list, AOP crown polygons.
  - `data/raw/ESS-DIVE-Vegetation-Field-2025/` — 2025 cover, site metadata, species list.
  - `data/raw/ESS-DIVE-Spectra/` — extracted NEON AOP spectra for both years, with per-year wavelength tables and crown polygons.
  - `data/derived/`, `data/cache/` — gitignored.
  - `data/small_reference/` — committed; small canonical inputs only.
- `output/` — figures and tables (gitignored). Commit canonical figures to `docs/figures/` if they need to live in the repo.
- `docs/` — methods notes, manuscript drafts.
- `DESCRIPTION` — exists only to drive `renv`; this is **not** an R package.

## Conventions

- **R-based analysis.** `terra`, `sf`, `stars`, `rSDP`, `tidyverse`. Env managed by `renv`.
- **Number-prefixed scripts** so the pipeline order is obvious.
- **CRS-explicit.** Don't assume `EPSG:32613` — pass it through. Field data in WGS84 needs to be reprojected.
- **SDP access via `rSDP`** (not raw `s3://` URLs).
- **`output/` is gitignored.** Move keepers to `docs/figures/`.
- No tests; this is an exploratory pipeline.

## Common commands

```bash
R -e 'renv::restore()'       # install deps from renv.lock
R -e 'source("code/01_load.R")'
```

## Things to be careful about

- **Taxonomy reconciliation across campaigns** — 2018 and 2025 mostly share names but expect drift (synonyms, splits/lumps, typos, authority changes). The crosswalk in `02_reconcile_taxonomy.R` is the single source of truth; downstream cover columns should use the harmonized name.
- **SiteID stability across years** — confirm whether plots resampled in 2025 reuse the 2018 IDs or have a separate scheme. If they differ, build a site crosswalk before stacking.
- **Spectra are pre-extracted at crown footprints**, not raw AOP cubes — the join unit is a crown polygon (one per site × year). 2018 and 2025 each have their own wavelength table; always carry the year's wavelengths alongside the spectral matrix.
- **CRS check on crowns:** the polygons are GeoJSON, so they may be `EPSG:4326` (GeoJSON convention) or already `EPSG:32613`. Verify before spatial joins.
- **Cluster interpretability vs separability is a trade-off** — record the chosen `k` and any species-list / dimensionality-reduction choices so future-Claude can re-derive them.
- **`output/` is gitignored** — any figure that needs to be archived has to be moved to `docs/figures/` explicitly.

## Reference implementation / cross-refs

- `~/code/rSDP` — SDP catalog access.
- `~/code/CHESS_trait_upscaling` — sibling CHESS analysis; similar AOP integration patterns.
- `~/code/chess_workshop` — workshop materials for the broader CHESS project.

## Open questions

- Do the 2018 and 2025 wavelength sets actually match (same NEON AOP sensor / processing) or do they need to be resampled to a common grid?
- SiteID scheme — are 2018 and 2025 plots the same physical locations with the same IDs, or do they need a crosswalk?
- Clustering algorithm — k-means / GMM / spectral / hierarchical? Decide after EDA.
- How to weight species composition vs spectral distance in a joint clustering — concatenate features, fit separately and reconcile, or constrain one with the other?
