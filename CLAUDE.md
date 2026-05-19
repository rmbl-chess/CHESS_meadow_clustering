# CHESS_meadow_clustering

Joint meadow + shrub classifier for NEON AOP imagery in the CHESS project area (East River / Upper Taylor / Almont), plus a per-pixel sampling-priority deliverable that ranks where new field plots would do the most good. Audience is the CHESS team; downstream consumer is a supervised classifier over NEON AOP tiles. SDP CRS is `EPSG:32613`.

## Data pipeline

Two field sampling campaigns (2018 and 2025) are reconciled, joined to NEON AOP spectra at crown footprints, classified, and then used to score basin pixels for fieldwork priority. Three logical phases:

1. **Meadow classes (01 → 24).** Reconcile taxonomy → combine cover → join spectra → spectra-first Ward clustering on 2025 plots (PCs 2–12 + snow-free DOY, z-scaled) → infer 2018 labels by Hellinger composition + DOY similarity. 31 classes (26 spectral S01–S26 + 5 monotypic-species overrides M01–M05). Exports a PCA model + target pixel list for Python-side AOP inference.
2. **Shrub classes (30 → 35).** 2025 single-species crown table + 2018 fractional_cover rows where a shrub-dominated genus has cover = 100 %. Synonym reconciliation. Salix is collapsed 4-way (`wolfii`, `boothii`, `planifolia`, `other`); Ribes / Juniperus collapsed to genus. 16 classes. Pixel-level RF with NDVI ≥ 0.20 filter; site-level fold CV.
3. **Joint classifier + leverage (36 → 43).** CHM extraction at crown centroids → joint training set (858 sites × 47 classes; 28 features including DOY and CHM) → balanced RF CV for per-class recall AND unweighted RF predictions on 5,354 cloud-extracted basin pixels → Mahalanobis novelty + per-pixel leverage → class summary table + sampling-priority GeoPackage. See `docs/sampling_priority_guide.md` for the team-facing description.

## Layout

- `code/` — R scripts grouped by phase, numbered from 01 within each phase, plus `python/` and `examples/`:
  - `code/meadow/01_load.R` → `code/meadow/24_export_classifier_model.R` — meadow pipeline (load → cluster → describe → export classifier model + target pixels).
  - `code/shrub/01_load.R` → `code/shrub/06_pixel_training.R` — shrub pipeline (load → reconcile → join spectra → separability → label set → pixel-level RF).
  - `code/joint/01_canopy_height.R` → `code/joint/08_figures.R` — joint pipeline (CHM → joint RF → inference → novelty → leverage → summary + figures).
  - `code/python/extract_aop_features.py` — per-pixel feature extraction over the icechunk virtual Zarr (resume-capable, parquet checkpoints).
  - `code/python/generate_aop_pc_maps.py` — per-tile 20-band PC COGs + per-domain VRTs (boto3 S3 listing, no AWS CLI).
  - `code/python/mosaic_pc_maps.py` — single multi-band COG per domain for off-server transfer; rebuilds missing VRTs from per-tile COGs.
- `data/` — provenance in `data/README.md`. Three source dirs under `data/raw/` (gitignored — fetch from ESS-DIVE).
  - `data/derived/` — gitignored except force-included handoff artifacts (see `.gitignore`).
  - `data/small_reference/` — committed; small canonical inputs (taxonomy crosswalk, narratives, woody-taxa reference).
- `docs/` — committed: `sampling_priority_guide.md` (team-facing), `figures/joint_*.pdf` (vector figures from script 43).
- `output/` — gitignored throwaway plots; move keepers to `docs/figures/`.
- `DESCRIPTION` — drives `renv` only; this is NOT an R package.

## Conventions

- **R-based analysis.** `terra`, `sf`, `stars`, `rSDP`, `ranger`, `tidyverse`. Env managed by `renv`.
- **Python tools** in `code/python/` use icechunk, xarray, rasterio, boto3. The chess-hub conda env. `arrow` is NOT installed for R — parquet ↔ R goes through Python conversion to CSV.
- **Phase subdirectories**: `code/meadow/`, `code/shrub/`, `code/joint/`. Each script is numbered from 01 within its subdir; run in numeric order. The directory is the phase context, so script names drop redundant phase prefixes (e.g., `code/shrub/01_load.R`, not `01_shrub_load.R`).
- **CRS-explicit.** Don't assume `EPSG:32613` — pass it through.
- **SDP access via `rSDP`** (not raw `s3://` URLs). NEON CHM tifs are local Google Drive paths hard-coded at the top of `code/joint/01_canopy_height.R`.
- **Force-included `data/derived/` artifacts.** Handoff CSVs / GeoPackages are not gitignored (see `.gitignore` exceptions). Large parquets and intermediate RDS files are.
- **`output/` is gitignored.** Move keepers to `docs/figures/`.

## Common commands

```bash
R -e 'renv::restore()'                       # restore env
for f in code/meadow/*.R; do Rscript "$f"; done   # phase 1 — meadow
for f in code/shrub/*.R;  do Rscript "$f"; done   # phase 2 — shrub
for f in code/joint/*.R;  do Rscript "$f"; done   # phase 3 — joint
```

## Things to be careful about

- **Two RFs, two questions.** Use balanced-weighted RF for per-class recall (script 37, CV). Use UNWEIGHTED RF for inference predictions on basin pixels (script 38). Class-weighting at inference time turns the smallest class into a catch-all for uncertain pixels and over-predicts it ~20× — we hit this bug already (S26 wet-sedge predicted to 33 % of pixels before the fix).
- **CHM extraction strategy.** Use single-point reads at crown centroids, not polygon extraction. Polygon reads against cloud-mounted CHM tifs trigger one windowed HTTP request per crown — ~30 min total. Centroid point reads bulk through in < 1 s per domain.
- **Salix is at genus + 4 by design.** Don't try to recover species-level Salix separation — the centroid dendrogram puts all 12 Salix binomials in one cluster and within-Salix DOY adds < 2 pp accuracy. Treat `Salix other` as a real class, not a residual bin.
- **Salix planifolia exception**: late-melt (DOY 155 +); `planifolia` and `glauca` are the only Salix species DOY genuinely separates.
- **2018-only meadow clustering caveat**: 2018 AOP spectra have atmospheric-correction drift; 2018 sites get labels via Hellinger composition + DOY similarity to 2025 plots, with a `confidence` column.
- **Taxonomy synonyms baked in**: Pentaphylloides floribunda → Dasiphora fruticosa; Distegia involucrata → Lonicera involucrata; Psychrophila leptosepala → Caltha leptosepala. Check `code/meadow/02_reconcile_taxonomy.R` and `code/shrub/02_taxonomy.R`.
- **Shrub site_type filter**: 2025 spectra has `site_type ∈ {Meadow, Shrub, Tree}`. The shrub pipeline filters to `site_type == "Shrub"` before averaging — don't drop this filter or you contaminate the shrub training set with co-located meadow pixels.
- **Inference-pixel bias toward meadows.** The 5,354-pixel cloud-extracted set was filtered through R3D018 landcover class 3 + neighborhood-purity; shrub leverage is undercounted as a result. For shrub priorities, work from `class_summary_table.csv` directly.
- **NDVI filter for shrub pixel training**: 0.20 threshold rescues sagebrush (Artemisia tridentata) which sits in the 0.30–0.45 NDVI band. Min 3 pixels per site is a coupled knob — raising it drops small-crown sites entirely.
- **CRS check on crowns:** the 2018 GeoJSON is in `EPSG:4326`; the 2025 GeoJSON declares `urn:ogc:def:crs:EPSG::32613` in its header. Reproject in `code/meadow/01_load.R`.
- **`output/` is gitignored** — figures that need to live in the repo go in `docs/figures/` (committed via gitignore exception).

## Reference implementation / cross-refs

- `~/code/rSDP` — SDP catalog access.
- `~/code/CHESS_trait_upscaling` — sibling CHESS analysis; similar AOP integration patterns.
- `~/code/chess_workshop` — workshop materials for the broader CHESS project.

## Open questions

- A shrub-targeted inference run is planned (current inference is meadow-biased via R3D018 class 3 filter).
- ESS-DIVE DOIs + 2018↔2025 site crosswalk still pending in `data/README.md`.
- 2018 wavelength drift (~2 nm) is handled by nearest-index match; resample only if a band-specific artifact shows up.
