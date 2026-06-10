# CHESS_meadow_clustering

Joint meadow + shrub classifier for NEON AOP imagery in the CHESS project area (East River / Upper Taylor / Almont), plus a per-pixel sampling-priority deliverable that ranks where new field plots would do the most good. Audience is the CHESS team; downstream consumer is a supervised classifier over NEON AOP tiles. SDP CRS is `EPSG:32613`.

## Data pipeline

Two field sampling campaigns (2018 and 2025) are reconciled, joined to NEON AOP spectra at crown footprints, classified, and then used to score basin pixels for fieldwork priority. Three logical phases:

1. **Meadow classes (`code/meadow/` 01 → 20).** Reconcile taxonomy → combine cover → join AOP spectra, applying a per-band **NDVI-stratified 2018→2025 radiometric correction** so 2018/2025 cluster together (2018 is no longer dropped; the 2026 supplemental plots are extracted from 2025 AOP and need no correction) → spectra-first Ward clustering (variant_G: PCs 2–12 + snow-free DOY, z-scaled; PC1 dropped — a sweep showed it adds nuisance, not coherence) → **gated composition subclustering**: split each spectral cluster by species-Hellinger only where the split is RF-mappable (min CV recall ≥ 0.25 on the 22-feature inference space), splitting the *post-monotypic residual* membership; sub-classes below 0.40 recall are flagged `needs_ancillary` (need topographic-wetness/phenology covariates later — see `subclass_mappability.csv`). Monotypic-species overrides at ≥70% cover. **34 meadow classes** (spectral S## + sub-classes S06.a/b, S08.a/b/c + 5 monotypic species; the 2026 ALMO sagebrush campaign reshaped the partition so the prior S01/S09/S20 sub-splits no longer form and nothing is currently flagged `needs_ancillary`). Then `19_label_descriptions.R` (IndVal narratives) and `20_ecosystem_crosswalk.R` (NatureServe crosswalk, below). Exports the PCA model + target pixels for Python AOP inference.
2. **Shrub classes (`code/shrub/` 01 → 06).** 2025 single-species crown table + 2018 fractional_cover rows where a shrub-dominated genus has cover = 100 %. Synonym reconciliation. Salix is collapsed 4-way (`wolfii`, `boothii`, `planifolia`, `other`); Ribes / Juniperus collapsed to genus. 17 classes (the 2026 supplemental crowns reinforce Alnus/Cornus/Lonicera). Pixel-level RF with NDVI ≥ 0.20 filter; site-level fold CV.
3. **Joint classifier + maps + leverage (`code/joint/` 01 → 15).** CHM at crown centroids → joint training set (**927 sites × 51 classes** = 34 meadow + 17 shrub; **22 features** = 20 spectral PCs + snow-free DOY + CHM; the 6 narrow-band indices are dropped at inference) → balanced RF CV for per-class recall AND unweighted RF for inference. Two inference products: (a) **per-pixel sampling priority** on 5,064 cloud-extracted basin pixels with CHM coverage (Mahalanobis novelty + leverage → `sampling_priority.gpkg`, `class_summary_table.csv`; `docs/sampling_priority_guide.md`); (b) **wall-to-wall classified + confidence COGs** per domain (`09_inference.R` → `data/derived/aop_classified/{DOM}_{class,confidence,labeled}_3m_v1.tif` + QGIS `.qml`, with a 4-family physiognomy color scheme in `10_label_rasters.R`) — **note: the COGs are still on the prior 40/56-class system; re-inference (`09/10`) is deferred**. Inference reads PC mosaics regenerated on the *current* PCA basis (`data/derived/aop_pc_maps_mosaic/`), NOT the older JPL_delivered mosaics — see the basis caveat below.

**NatureServe IVC crosswalk.** `code/python/natureserve_fetch.py` queries the NatureServe Explorer API for candidate International Vegetation Classification communities per class's diagnostic species (reconciled via `data/small_reference/species_natureserve_crosswalk.csv`) and caches them (`data/derived/natureserve_cache.json`). `code/meadow/20_ecosystem_crosswalk.R` scores each class's IndVal assemblage against the cached Colorado catalog → `natureserve_candidates.csv` (top-3 draft for curation). The top-1 draft community is surfaced on the labeled maps' RAT for spatial review.

## Layout

- `code/` — R scripts grouped by phase, numbered from 01 within each phase, plus `python/` and `examples/`:
  - `code/meadow/01_load.R` → `code/meadow/25_diagnostics.R` — meadow pipeline (load → reconcile → join spectra+correction → cluster `10` → gated subcluster `11` → describe `19` → NatureServe crosswalk `20` → export classifier model `24` + target pixels `23` + diagnostics `25`).
  - `code/shrub/01_load.R` → `code/shrub/06_pixel_training.R` — shrub pipeline (load → reconcile → join spectra → separability → label set → pixel-level RF).
  - `code/joint/01_canopy_height.R` → `code/joint/15_...R` — joint pipeline (CHM `01` → joint RF `02` → sampling-priority pixels `03–06` → figures `08` → wall-to-wall inference COGs `09` → labeled rasters `10` → landcover mask `11`; `12–15` are the 2018→2025 year-effect diagnosis + radiometric-correction fit).
  - `code/python/extract_aop_features.py` — per-pixel feature extraction over the icechunk virtual Zarr (resume-capable, parquet checkpoints).
  - `code/python/generate_aop_pc_maps.py` — per-tile 20-band PC COGs + per-domain VRTs for 2025 AND 2018 (the 2018 run applies the NDVI-stratified correction). boto3 S3 listing, no AWS CLI.
  - `code/python/mosaic_pc_maps.py` — single multi-band COG per domain (BAND interleave + ZSTD + deep overview pyramid for fast GIS rendering); rebuilds missing VRTs from per-tile COGs.
  - `code/python/natureserve_fetch.py` — NatureServe Explorer API → cached candidate IVC communities for the crosswalk.
- `data/` — provenance in `data/README.md`. Three source dirs under `data/raw/` (gitignored — fetch from ESS-DIVE).
  - `data/derived/` — gitignored except force-included handoff artifacts (see `.gitignore`): training sets, sampling-priority gpkgs, `class_summary_table.csv`, NatureServe cache + candidates, `subclass_mappability.csv`, and the labeled-raster sidecars under `aop_classified/` (`.qml` + `class_lookup_labeled.csv`). The large class/confidence COGs themselves are gitignored — local only, moved off-server via the python mosaic tools.
  - `data/small_reference/` — committed canonical inputs: taxonomy crosswalk, `label_community_names.csv` (narratives), `class_categories.csv` (moisture × elevation, optional `category` color override), `species_natureserve_crosswalk.csv`, `year_effect_correction_2018_to_2025_by_ndvi.csv`.
- `docs/` — committed: `sampling_priority_guide.md` (team-facing), `figures/joint_*.pdf` (from `08_figures.R`).
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

- **Two RFs, two questions.** Use balanced-weighted RF for per-class recall (`02_training.R`, CV). Use UNWEIGHTED RF for inference predictions (`09_inference.R` refits unweighted on the 22 features). Class-weighting at inference time turns the smallest class into a catch-all for uncertain pixels and over-predicts it ~20× — we hit this bug already (S26 wet-sedge predicted to 33 % of pixels before the fix).
- **CHM extraction strategy.** Use single-point reads at crown centroids, not polygon extraction. Polygon reads against cloud-mounted CHM tifs trigger one windowed HTTP request per crown — ~30 min total. Centroid point reads bulk through in < 1 s per domain.
- **Salix is at genus + 4 by design.** Don't try to recover species-level Salix separation — the centroid dendrogram puts all 12 Salix binomials in one cluster and within-Salix DOY adds < 2 pp accuracy. Treat `Salix other` as a real class, not a residual bin.
- **Salix planifolia exception**: late-melt (DOY 155 +); `planifolia` and `glauca` are the only Salix species DOY genuinely separates.
- **2018 now clusters directly (was 2025-only).** The NDVI-stratified 2018→2025 radiometric correction (`year_effect_correction_2018_to_2025_by_ndvi.csv`, applied in `04_join_spectra.R`) removes the year drift, so both campaigns cluster together; only the rare spectra-less 2018 sites fall back to Hellinger+DOY nearest-cluster. Single-year classes are legitimate (elevation-stratified sampling: 2018 more alpine, 2025 more low sagebrush) — the subcluster gate does NOT require both years.
- **PCA basis must match train ↔ inference.** `09_inference.R` reads PC mosaics in `data/derived/aop_pc_maps_mosaic/`, regenerated by `generate_aop_pc_maps.py` from the CURRENT `aop_classifier_pca.csv`. The older `JPL_delivered/` mosaics were on a stale basis (sign-flipped PC1, rotated PC3/4) and produced spatially-coherent meadow↔shrub MISclassification. Regenerate the mosaics whenever the meadow PCA changes; verify with a co-located train-vs-mosaic PC correlation (should be strongly positive).
- **Gated subclustering runs on residual membership.** `11` splits a spectral cluster only where every sub-class is RF-mappable (`recall_bar` 0.25, try k=3 then k=2), AFTER monotypic-stand sites leave (split the class as it deploys). Sub-classes below `ancillary_recall` (0.40) are flagged `needs_ancillary` — kept for ecological fidelity, to be mapped later with topographic-wetness / phenology covariates.
- **NatureServe name reconciliation.** Diagnostic species must be mapped to NatureServe/IVC nomenclature before searching or a class silently gets zero candidates (Veratrum tenuipetalum → V. californicum + V. viride alias for the CO community; Bromopsis → Bromus; Vaccinium caespitosum spelling). See `species_natureserve_crosswalk.csv`; a genus-only token can't carry a match alone.
- **Taxonomy synonyms baked in**: Pentaphylloides floribunda → Dasiphora fruticosa; Distegia involucrata → Lonicera involucrata; Psychrophila leptosepala → Caltha leptosepala. Check `code/meadow/02_reconcile_taxonomy.R` and `code/shrub/02_taxonomy.R`.
- **Shrub site_type filter**: 2025 spectra has `site_type ∈ {Meadow, Shrub, Tree}`. The shrub pipeline filters to `site_type == "Shrub"` before averaging — don't drop this filter or you contaminate the shrub training set with co-located meadow pixels.
- **Inference-pixel bias toward meadows.** The cloud-extracted set (5,354 pixels; 5,064 with CHM coverage are scored) was filtered through R3D018 landcover class 3 + neighborhood-purity; shrub leverage is undercounted as a result. For shrub priorities, work from `class_summary_table.csv` directly.
- **NDVI filter for shrub pixel training**: 0.20 threshold rescues sagebrush (Artemisia tridentata) which sits in the 0.30–0.45 NDVI band. Min 3 pixels per site is a coupled knob — raising it drops small-crown sites entirely.
- **CRS check on crowns:** the 2018 GeoJSON is in `EPSG:4326`; the 2025 GeoJSON declares `urn:ogc:def:crs:EPSG::32613` in its header. Reproject in `code/meadow/01_load.R`.
- **`output/` is gitignored** — figures that need to live in the repo go in `docs/figures/` (committed via gitignore exception).

## Reference implementation / cross-refs

- `~/code/rSDP` — SDP catalog access.
- `~/code/CHESS_trait_upscaling` — sibling CHESS analysis; similar AOP integration patterns.
- `~/code/chess_workshop` — workshop materials for the broader CHESS project.

## Open questions

- **`needs_ancillary` mechanism.** Sub-class splits below 0.40 CV recall are kept for ecological fidelity but flagged for later resolution with ancillary covariates (topographic wetness index, phenology). After the 2026 ALMO sagebrush plots reshaped the partition, only S06/S08 split and **none are currently flagged** — but the gate still applies, so flagged sub-classes can reappear on future reclusters. See `subclass_mappability.csv`.
- **NatureServe crosswalk is a draft.** `natureserve_candidates.csv` holds top-3 per class; the curated `class_ecosystem_crosswalk.csv` + wiring the recognized community + G-rank into `19`/`06`/`10` is pending. `10` currently surfaces only the top-1 draft on the map for review.
- A shrub-targeted inference run is planned (sampling-priority pixels are meadow-biased via the R3D018 class-3 filter).
- ESS-DIVE DOIs + 2018↔2025 site crosswalk still pending in `data/README.md`.
- 2018 wavelength drift (~2 nm) is handled by nearest-index match; resample only if a band-specific artifact shows up.
