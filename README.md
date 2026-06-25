# CHESS_meadow_clustering

Build a joint meadow + shrub classifier for NEON AOP imagery in the East River / Upper Taylor / Almont basins, and use it to direct the next field-sampling campaign toward the highest-leverage gaps in the training data.

## What this is

The CHESS map classifier needs training samples that are (a) spectrally separable in 1 m NEON AOP imagery and (b) interpretable as recognizable plant communities. This repo:

1. Reconciles 2018 + 2025 + 2026 CHESS field vegetation campaigns into one harmonized cover table.
2. Clusters **2018 + 2025 + 2026** meadow plots (an NDVI-stratified 2018→2025 radiometric correction brings 2018/2025 onto one spectral basis; the 2026 supplemental plots are extracted directly from 2025 AOP, so they need no correction) into spectrally distinct community types, then **gated composition-subclustering** splits each spectral cluster further only where the split stays spectrally mappable. **35 meadow classes** (low-mappability sub-classes are flagged `needs_ancillary`).
3. Assembles a parallel shrub training set (species-level labels, with Salix collapsed to 4 spectrally distinguishable groups). 17 shrub classes.
4. Joins the two into a single classifier with **52 classes** (632 meadow + 323 shrub sites). Deployment features = 20 spectral PCs + snow-free DOY + canopy height (the 6 narrow-band indices are used in CV but dropped at inference).
5. Selects stratified target pixels, runs Python extraction + classification on a cloud server, and computes a per-pixel novelty + leverage score so the next field campaign targets the gaps.
6. Produces **wall-to-wall classified + confidence COGs** per AOP domain, and a **NatureServe IVC community crosswalk** mapping each class to a recognized vegetation community.

**Current state.** Joint 52-class RF reaches ~67 % CV accuracy (0.679 unweighted / 0.674 balanced) after the 2026 supplemental campaign (98 plots across ALMO / CRBU / UPTA — 84 meadow + 12 shrub in training) was folded in. Per-pixel sampling priority (5,064 basin pixels: predicted class + Mahalanobis novelty + leverage) and the wall-to-wall classified/confidence/labeled COGs are regenerated on the 52-class classifier. A NatureServe crosswalk draft (`natureserve_candidates.csv`) awaits curation. See `docs/sampling_priority_guide.md` for the team-facing deliverable.

> **Note:** the meadow PCA basis was refit and the classes resplit; inference now reads PC mosaics regenerated on the current basis (an older off-basis set caused meadow↔shrub misclassification — see Key choices).

## Pipeline

R scripts are grouped into one subdirectory per phase under `code/`; each phase is renumbered from 01 so the directory listing matches run order. `renv` env assumed restored.

### Phase 1 — Meadow classes (`code/meadow/`)

| Step | Script | Output |
|---|---|---|
| Load 2018/2025 field data | `01_load.R` | per-year cover + species + crowns + spectra RDS |
| Reconcile taxonomy | `02_reconcile_taxonomy.R` | `data/derived/taxonomy_crosswalk.csv` |
| Combine cover | `03_combine_cover.R` | `data/derived/cover_combined.rds` |
| Join AOP spectra + year correction | `04_join_spectra.R` | `data/derived/veg_spectra.rds` — L2-normalized per-site spectra; applies the NDVI-stratified 2018→2025 radiometric correction to 2018 so both years cluster together |
| Preprocess features | `05_preprocess_features.R` | brightness-normed spectra + Hellinger cover + PCA (fit on 2025, both years projected) |
| **Architecture B: spectra-first cluster** | `10_cluster_spectra.R` | Ward / variant_G = PCs 2–12 + snow-free DOY, z-scaled (PC1 dropped; K-selection exploration in `_attic/`) |
| **Gated composition sub-clusters** | `11_subcluster_composition.R` | splits each spectral cluster by species-Hellinger only where the split is RF-mappable (try k=3/2, min recall ≥ 0.25) on the post-monotypic residual; flags `needs_ancillary` (< 0.40). `final_clusters_B.rds`, `subclass_mappability.csv` |
| Environment join | `13_extract_environment.R` | snow-free DOY per site |
| Training-sample export | `17_export_training_samples.R` | `training_samples_*.{csv,gpkg}` |
| Class descriptions + narratives | `19_label_descriptions.R` | `data/derived/label_descriptions.csv` (IndVal + abundant species + physiognomy) AND `data/small_reference/label_community_names.csv` (auto-drafted starter narratives; preserves user-curated rows on re-run) |
| **NatureServe IVC crosswalk** | `20_ecosystem_crosswalk.R` | scores each class's IndVal assemblage against the cached Colorado NatureServe catalog (species-F1 + ecology) → `natureserve_candidates.csv` (top-3 draft per class for curation). Needs `code/python/natureserve_fetch.py` run first to build the cache |
| **Target pixel selection** | `23_target_pixels.R` | `target_pixels.{csv,gpkg}` (6,000 stratified meadow pixels) |
| Classifier-feature export | `24_export_classifier_model.R` | `aop_classifier_*.{csv,json}` — R → Python handoff |
| **Diagnostic suite** (last) | `25_diagnostics.R` | One file with seven sections: cluster figures (A), env coherence (B), PC1/PC2 (C), inference QC (D), size distribution (E), K-sweep (F), year-effect (G). All output PNGs + RDS artifacts |

### Phase 2 — Shrub classes (`code/shrub/`)

| Step | Script | Output |
|---|---|---|
| Load shrub records | `01_load.R` | 2025 from `chess_shrub_site_cleaned.csv`; 2018 from `fractional_cover` filtered to shrub-dominated genera at 100 % cover |
| Reconcile shrub taxonomy | `02_taxonomy.R` | `shrub_taxonomy_crosswalk.csv` (Pentaphylloides → Dasiphora; Distegia → Lonicera; …) |
| Join AOP spectra | `03_join_spectra.R` | `shrub_veg_spectra.rds` (site_type=Shrub for 2025) |
| Separability + dendrogram | `04_separability.R` | RF CV recall, centroid Ward dendrogram → identifies the Salix complex |
| Final shrub label set | `05_label_set.R` | 17 classes: Salix split 4-way, Ribes / Juniperus collapsed, small-N classes dropped or aggregated |
| Pixel-level shrub RF | `06_pixel_training.R` | Pixel-level training with NDVI ≥ 0.20 filter + DOY covariate; site-level fold CV. ~59 % accuracy on 17 classes |

### Phase 3 — Joint classifier + leverage analysis (`code/joint/`)

| Step | Script | Output |
|---|---|---|
| Canopy height extraction | `01_canopy_height.R` | `canopy_height.rds` — 1 m NEON CHM at every crown centroid (single-point extract, all three domains in < 1 s) |
| **Joint training set + RF** | `02_training.R` | `joint_training_set.{rds,csv}`, `punch_list.csv`. 955 sites × 52 classes; CV on 28 features but the deployed RF uses 22 (drops the 6 indices); balanced 5-fold CV ~ 67 % |
| Inference on basin pixels | `03_predict_inference_pixels.R` | `inference_predictions.csv` — unweighted RF predictions on the 5,064 extracted meadow pixels with CHM coverage; refreshes the punch list with `predicted_n_pixels` |
| Landscape novelty | `04_landscape_distance.R` | Mahalanobis distance from every inference pixel to every class centroid (pooled within-class covariance). `inference_pixel_distances.csv`, `novelty_by_class.csv`, `novelty_by_hex.gpkg` |
| **Per-pixel sampling priority** | `05_sampling_priority.R` | `sampling_priority.gpkg` — leverage = `nearest_d / sqrt(n_training_for_predicted_class)` per pixel + per-class top-10 candidate sites |
| Joint summary outputs | `06_summary_outputs.R` | `class_summary_table.csv` (one row per class with N, prevalence, recall, indicator + abundant taxa, description, median leverage) AND `joint_training.gpkg` (two layers, `training_sites_crowns` + `training_sites_points`, with class metadata for QGIS review) |
| Joint figures | `08_figures.R` | `docs/figures/joint_*.pdf` — recall + leverage scatter + feature space + confusion matrix |
| **Wall-to-wall inference** | `09_inference.R` | refits the unweighted 22-feature joint RF and predicts class + max-probability per 3 m pixel across each domain's PC mosaic; masks to CHM ≤ 4 m AND R3D018 meadow/shrub landcover. `aop_classified/{DOM}_{class,confidence}_3m_v1.tif` (COGs) + `class_lookup.csv` |
| **Labeled rasters** | `10_label_rasters.R` | attaches RAT + 4-family physiognomy color table + top-1 NatureServe draft community → `{DOM}_class_3m_v1_labeled.tif` + `.qml` + `class_lookup_labeled.csv` |
| Landcover mask | `11_landcover_mask.R` | standalone R3D018 meadow/shrub mask (now folded into `09`'s single pass) |
| Year-effect + correction | `12–15_*.R` | 2018→2025 spectral year-effect diagnosis and the NDVI-stratified per-band radiometric correction fit (consumed by `code/meadow/04_join_spectra.R`) |

### Python tools (`code/python/`)

| Script | Purpose |
|---|---|
| `extract_aop_features.py` | Extract 20 PCs + 6 narrow-band indices at each target pixel (3 m, per-pixel L2 norm → 3 × 3 mean → water-band mask → PCA project). Periodic parquet checkpointing + resume on re-run. |
| `generate_aop_pc_maps.py` | Write 20-band PC COGs per AOP tile + per-domain VRTs, for **2025 and 2018** (the `--year 2018` run applies the NDVI-stratified radiometric correction; 2018 is CRBU-only). Projects through the current `aop_classifier_pca.csv`. boto3 tile inventory. |
| `mosaic_pc_maps.py` | Collapse per-tile COGs into one multi-band COG per domain (**BAND interleave + ZSTD + deep overview pyramid** — ~10–30× faster single-band rendering in QGIS). Rebuilds the VRT if missing. |
| `natureserve_fetch.py` | Query the NatureServe Explorer API for candidate IVC communities per class's reconciled diagnostic species; cache to `natureserve_cache.json` (with fetch date + citation) for the offline crosswalk match. |

## Key choices

- **Architecture B (spectra-first clustering).** Meadow classes are spectral clusters first, composition second. A clustering that isn't spectrally separable can't be mapped.
- **Both years cluster together (2018→2025 correction).** 2018 AOP spectra had atmospheric-correction drift; an NDVI-stratified per-band correction (`year_effect_correction_2018_to_2025_by_ndvi.csv`, applied in `04_join_spectra.R`) removes it, so 2018 plots cluster directly alongside 2025 rather than being inferred. Single-year classes are legitimate (the elevation-stratified sampling design: 2018 more alpine, 2025 more low-elevation sagebrush).
- **Gated, mappability-aware subclustering.** Spectral clusters split into community sub-classes ONLY where each sub-class clears an RF-mappability bar (CV recall ≥ 0.25 on the 22-feature inference space), evaluated on the post-monotypic residual membership. Sub-classes below 0.40 are kept but flagged `needs_ancillary` — we preserve ecologically real distinctions even when 3 m spectra can't yet separate them, to be resolved later with topographic-wetness / phenology covariates.
- **Inference reads PC mosaics on the CURRENT basis.** Training and inference must share one PCA basis. The `JPL_delivered/` mosaics were on a stale basis (sign-flipped PC1, rotated higher PCs) and caused spatially-coherent meadow↔shrub misclassification; `09` now reads mosaics regenerated by `generate_aop_pc_maps.py` from the current `aop_classifier_pca.csv`. A co-located train-vs-mosaic PC correlation is the check.
- **Physiognomy color scheme.** Labeled rasters use 4 HCL hue families (shrubland blue-purple, grassland tan-brown, forb-meadow green, wetland pink-red); within-family lightness/hue ramps order by elevation→snowmelt, so hue = structure and the in-family gradient = community.
- **NatureServe crosswalk.** Each class is matched to a recognized IVC community by its IndVal diagnostic-species assemblage (reconciled to NatureServe nomenclature) against a cached Colorado catalog.
- **Salix 4-way split.** Centroid dendrogram from `33_shrub_separability.R` shows all 12 Salix binomials in one cluster. We keep the three highest-N species (`wolfii`, `boothii`, `planifolia`) and collapse the rest to `Salix other`. Within-Salix DOY barely helps — kept at genus-plus-four resolution.
- **Pixel-level shrub training with NDVI ≥ 0.20 filter + min 3 pixels/site.** Trades a slight overall-accuracy loss against rescuing small classes (Purshia 0 → 1.0, Betula 0 → 0.5, Prunus 0 → 0.67) by giving them more training signal per site.
- **CHM extracted at crown centroids (single point per crown).** 1 m NEON CHM is at the same scale as the 3 m AOP feature recipe — single-pixel reads are 100× faster than polygon extraction over the cloud-mounted rasters (< 1 s total vs ~30 min for polygons).
- **Two RFs for two questions.** Balanced-weighted RF (script 37, CV) measures per-class quality. Unweighted RF (script 38, inference) generates the realistic class-proportion map. Class-weighting at inference time turns the smallest class into a catch-all for uncertain pixels and over-predicts it dramatically (S26 with n = 7 was being assigned to 33 % of pixels before this was fixed).
- **Per-pixel leverage score.** `leverage = nearest_d / sqrt(n_training_for_predicted_class)` — high when a pixel is spectrally far from any class centroid AND its predicted class is undersampled. Drives the sampling-priority deliverable.
- **CRS-explicit.** SDP / AOP layers are `EPSG:32613` (UTM Zone 13N); field-data GeoJSONs are reprojected from WGS84 before any spatial join.

## Outputs

Canonical artifacts under `data/derived/` (force-included via `.gitignore` exceptions; the rest of `data/derived/` is gitignored):

| Path | What |
|---|---|
| `data/small_reference/label_community_names.csv` | Auto-drafted (curatable) narratives for the 35 meadow classes |
| `data/derived/training_samples_sites.csv` | Meadow training: one row per (plot, year) with `final_label`, confidence tier, env covariates |
| `data/derived/training_samples_crowns.gpkg` | Meadow crown polygons for QGIS review |
| `data/derived/shrub_training_set.csv` | Shrub training table (N=323) |
| `data/derived/shrub_label_crosswalk.csv` | canonical_binomial → final_label mapping |
| `data/derived/canopy_height.rds` | 1 m CHM at every crown centroid (1393 sites) |
| `data/derived/joint_training_set.csv` | Joint meadow+shrub training (955 sites × 52 classes) |
| `data/derived/subclass_mappability.csv` | Per split sub-class: gate CV recall + `needs_ancillary` flag |
| `data/derived/joint_training.gpkg` | Same as above with crown geometries + class metadata |
| `data/derived/punch_list.csv` | Class-level summary: training N, recall, predicted prevalence, top confusions, augmentation priority |
| `data/derived/class_summary_table.csv` | Per-class table with detailed IndVal indicator + abundant taxa strings, descriptions, median leverage |
| `data/derived/target_pixels.{csv,gpkg}` | 6,000 stratified meadow pixels for inference |
| `data/derived/aop_classifier_*.{csv,json}` | R → Python handoff (PCA loadings, preprocessing config, index formulas) |
| `data/derived/inference_predictions.csv` | Joint-RF predictions on the 5,064 cloud-extracted pixels (with CHM coverage) |
| `data/derived/inference_pixel_distances.csv` | Mahalanobis distance to every class centroid per pixel |
| `data/derived/novelty_by_class.csv`, `.../novelty_by_hex.gpkg` | Aggregated novelty for spatial review |
| **`data/derived/sampling_priority.gpkg`** | **Primary deliverable** — 5,064 candidate pixels ranked by leverage |
| `data/derived/sampling_priority_top.csv` | Top 10 per class (295 candidate sites for fieldwork planning) |
| `data/derived/aop_classified/{DOM}_{class,confidence}_3m_v1.tif` | Wall-to-wall classified + confidence COGs per domain (gitignored; local only) |
| `data/derived/aop_classified/{DOM}_class_3m_v1_labeled.tif` + `.qml` + `class_lookup_labeled.csv` | Labeled rasters with physiognomy colors + draft NatureServe community (sidecars committed) |
| `data/derived/natureserve_cache.json`, `natureserve_candidates.csv` | NatureServe IVC catalog snapshot + top-3 crosswalk draft per class |
| `docs/figures/joint_*.pdf` | Class recall + leverage scatter + feature space + confusion (vector PDFs) |
| `docs/sampling_priority_guide.md` | Team-facing one-pager: how to interpret + use `sampling_priority.gpkg` |

After running the Python tools on a cloud server:

| Path | What |
|---|---|
| `data/derived/target_pixel_features.parquet` | Per-target-pixel features ready for classifier inference |
| `data/derived/aop_pc_maps/{DOMAIN}/pc_{e}_{n}.tif` | 20-band PC COGs at 3 m, one per AOP tile |
| `data/derived/aop_pc_maps/{DOMAIN}_pc_mosaic.vrt` | Per-domain mosaic for QGIS |
| `data/derived/aop_pc_maps_mosaic/{DOMAIN}_pc_mosaic.tif` | Single multi-band COG per domain for off-server transfer |

## Layout

```
code/
  meadow/              Phase 1 — meadow classes (01_load.R → 24_export_classifier_model.R)
  shrub/               Phase 2 — shrub classes (01_load.R → 06_pixel_training.R)
  joint/               Phase 3 — joint training + leverage (01_canopy_height.R → 08_figures.R)
  python/              Cloud-side feature extraction + COG generation + mosaic
  examples/            Reference notebooks
data/raw/              ESS-DIVE working copies (gitignored, fetch from ESS-DIVE)
data/derived/          Pipeline outputs (gitignored except force-included handoffs)
data/small_reference/  Small canonical inputs (committed)
docs/                  sampling_priority_guide.md + figures/ (committed)
output/                Throwaway diagnostic plots (gitignored)
DESCRIPTION            Drives renv only — this is NOT an R package
renv.lock              R dependency lockfile
```

See `data/README.md` for ESS-DIVE source dirs and file-level provenance, and `docs/sampling_priority_guide.md` for the team-facing guide to the sampling-priority deliverable.

## Reproducing

```bash
git clone https://github.com/rmbl-chess/CHESS_meadow_clustering
cd CHESS_meadow_clustering
R -e 'renv::restore()'

# Fetch the ESS-DIVE source dirs into data/raw/ (see data/README.md):
#   ESS-DIVE-Vegetation-Field-2018/
#   ESS-DIVE-Vegetation-Field-2025/
#   ESS-DIVE-Spectra/
# Also: a local copy of R3D018 1 m landcover (for code/meadow/23_target_pixels.R
# and code/joint/09_inference.R), and the NEON CHM + JPL PC mosaic paths
# hard-coded at the top of code/joint/01_canopy_height.R / 09_inference.R;
# edit as needed for your environment.

# Phase 1: meadow classes + target pixels + R-side feature export
for f in code/meadow/*.R; do Rscript "$f"; done

# Phase 2: shrub classes
for f in code/shrub/*.R; do Rscript "$f"; done

# Phase 3: joint classifier + leverage analysis
#   (01 runs in < 1 s; 02 takes ~ 1 min; 03-08 each < 30 s)
for f in code/joint/*.R; do Rscript "$f"; done
```

Python pipeline (designed for a cloud server with AOP S3 access):

```bash
# Local smoke tests
python code/python/extract_aop_features.py --max-pixels 50 --output /tmp/test.parquet --workers 1
python code/python/generate_aop_pc_maps.py --tile-limit 1 --output-dir /tmp/pc_test
python code/python/mosaic_pc_maps.py --input-dir /tmp/pc_test --output-dir /tmp/pc_mosaic

# Full cloud run
python code/python/extract_aop_features.py --workers 16 --memory-limit 8GB
python code/python/generate_aop_pc_maps.py
python code/python/mosaic_pc_maps.py        # consolidates per-tile COGs into one COG per domain
```

The Python feature extractor checkpoints to parquet every `--checkpoint-every` pixels (default 200) and skips already-completed pixels on a re-run, so an interrupted multi-hour job picks up where it left off.

## Status / open questions

- **`needs_ancillary` sub-classes.** Sub-class splits below 0.40 CV recall are kept but flagged for later resolution with ancillary covariates (topographic wetness, phenology). In the current 35-class meadow set S01, S04, S18, and S25 split; **S18.a/b and S25.a/b are flagged `needs_ancillary`** (the long-broad S01 finally splits a/b at recall 0.76/0.41). The split set shifts with each reclustering. See `subclass_mappability.csv`.
- **NatureServe crosswalk is a draft.** `natureserve_candidates.csv` is top-3 per class; curating the final `class_ecosystem_crosswalk.csv` and wiring the recognized community + G-rank into the descriptions/labels is pending. `10` currently surfaces only the top-1 draft on the map.
- **Inference is meadow-biased.** The sampling-priority set is filtered to R3D018 meadow class with strict neighborhood-purity rules; shrub leverage is undercounted. A shrub-targeted inference run is planned to balance the picture.
- **Small-N classes need spatial diversity, not just more samples.** Prunus virginiana (n=3) and Purshia tridentata (n=4) remain concentrated in a single drainage — augmentation should prioritize new geographic locations, not just more N. The 2026 campaign rescued Alnus incana (4→9) and Cornus sericea (1→6); still-critical gaps are now Symphoricarpos rotundifolius, Lonicera involucrata, and the monotypic wet-forb meadows (Caltha, Osmorhiza) — see `punch_list.csv`.
- **ESS-DIVE DOIs and 2018↔2025 site crosswalk** still pending in `data/README.md`.
- **2018 wavelength drift** is implicit (nearest-index match, max ~2 nm). Resampling onto a common grid is on the list if a band-specific artifact shows up downstream.

## License

Internal CHESS-team use; no LICENSE file.
