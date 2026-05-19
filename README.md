# CHESS_meadow_clustering

Build a joint meadow + shrub classifier for NEON AOP imagery in the East River / Upper Taylor / Almont basins, and use it to direct the next field-sampling campaign toward the highest-leverage gaps in the training data.

## What this is

The CHESS map classifier needs training samples that are (a) spectrally separable in 1 m NEON AOP imagery and (b) interpretable as recognizable plant communities. This repo:

1. Reconciles 2018 + 2025 CHESS field vegetation campaigns into one harmonized cover table.
2. Clusters 2025 meadow plots into spectrally distinct community types and infers 2018 labels by composition + phenology similarity. 31 meadow classes.
3. Assembles a parallel shrub training set (species-level labels, with Salix collapsed to 4 spectrally distinguishable groups). 16 shrub classes.
4. Joins the two into a single classifier with 47 classes (548 meadow + 310 shrub sites). Features = 20 spectral PCs + 6 narrow-band indices + snow-free DOY + canopy height.
5. Selects stratified target pixels across the basin, runs Python extraction + classification on a cloud server, and computes a per-pixel novelty + leverage score so the next field campaign can be targeted at the gaps.

**Current state.** Joint 47-class RF reaches 63–65 % balanced CV accuracy. Inference (unweighted RF) has been run on 5,354 meadow pixels across the three AOP domains, with predicted class and Mahalanobis novelty scored for each. A sampling-priority GeoPackage + class summary table are ready for the team — see `docs/sampling_priority_guide.md`.

## Pipeline

R scripts are grouped into one subdirectory per phase under `code/`; each phase is renumbered from 01 so the directory listing matches run order. `renv` env assumed restored.

### Phase 1 — Meadow classes (`code/meadow/`)

| Step | Script | Output |
|---|---|---|
| Load 2018/2025 field data | `01_load.R` | per-year cover + species + crowns + spectra RDS |
| Reconcile taxonomy | `02_reconcile_taxonomy.R` | `data/derived/taxonomy_crosswalk.csv` |
| Combine cover | `03_combine_cover.R` | `data/derived/cover_combined.rds` |
| Join AOP spectra | `04_join_spectra.R` | `data/derived/veg_spectra.rds` |
| Preprocess features | `05_preprocess_features.R` | brightness-normed spectra + Hellinger cover + PCA |
| **Architecture B: spectra-first cluster** | `10_cluster_spectra.R` | Ward / PCs 2–12 + snow-free DOY, z-scaled (K-selection exploration in `_attic/`, see attic README) |
| Composition sub-clusters | `11_subcluster_composition.R` | within-spectral subclusters by composition |
| Environment join | `13_extract_environment.R` | snow-free DOY per site |
| Training-sample export | `17_export_training_samples.R` | `training_samples_*.{csv,gpkg}` |
| Class descriptions | `19_label_descriptions.R` | IndVal + abundant species tables |
| Narrative drafts | `20_label_narratives.R` | `data/small_reference/label_community_names.csv` |
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
| Final shrub label set | `05_label_set.R` | 16 classes: Salix split 4-way, Ribes / Juniperus collapsed, small-N classes dropped or aggregated |
| Pixel-level shrub RF | `06_pixel_training.R` | Pixel-level training with NDVI ≥ 0.20 filter + DOY covariate; site-level fold CV. 72 % accuracy on 16 classes |

### Phase 3 — Joint classifier + leverage analysis (`code/joint/`)

| Step | Script | Output |
|---|---|---|
| Canopy height extraction | `01_canopy_height.R` | `canopy_height.rds` — 1 m NEON CHM at every crown centroid (single-point extract, all three domains in < 1 s) |
| **Joint training set + RF** | `02_training.R` | `joint_training_set.{rds,csv}`, `punch_list.csv`. 858 sites × 47 classes; 28 features; balanced 5-fold CV ~ 64 % |
| Inference on basin pixels | `03_predict_inference_pixels.R` | `inference_predictions.csv` — unweighted RF predictions on the 5,354 extracted meadow pixels; refreshes the punch list with `predicted_n_pixels` |
| Landscape novelty | `04_landscape_distance.R` | Mahalanobis distance from every inference pixel to every class centroid (pooled within-class covariance). `inference_pixel_distances.csv`, `novelty_by_class.csv`, `novelty_by_hex.gpkg` |
| **Per-pixel sampling priority** | `05_sampling_priority.R` | `sampling_priority.gpkg` — leverage = `nearest_d / sqrt(n_training_for_predicted_class)` per pixel + per-class top-10 candidate sites |
| Class summary table | `06_class_summary_table.R` | `class_summary_table.csv` — one row per class with N (2018, 2025, total), basin prevalence, recall, indicator + abundant taxa (full IndVal cov/freq/IV strings for meadows), description, median leverage |
| Joint training GeoPackage | `07_training_gpkg.R` | `joint_training.gpkg` — two layers (`training_sites_crowns`, `training_sites_points`) with class metadata for QGIS review |
| Joint figures | `08_figures.R` | `docs/figures/joint_*.pdf` — recall + leverage scatter + feature space + confusion matrix |

### Python tools (`code/python/`)

| Script | Purpose |
|---|---|
| `extract_aop_features.py` | Extract 20 PCs + 6 narrow-band indices at each target pixel (3 m, per-pixel L2 norm → 3 × 3 mean → water-band mask → PCA project). Periodic parquet checkpointing + resume on re-run. |
| `generate_aop_pc_maps.py` | Write 20-band PC COGs for every 2025 AOP tile + per-domain VRT mosaics. boto3 tile inventory (no AWS CLI dependency). |
| `mosaic_pc_maps.py` | Collapse per-tile COGs + VRTs into one multi-band COG per domain (streams through `gdal_translate -of COG`). Rebuilds the VRT from tile COGs if missing. For off-server data migration. |

## Key choices

- **Architecture B (spectra-first clustering).** Meadow classes are spectral clusters first, composition second. A clustering that isn't spectrally separable can't be mapped.
- **2025-only meadow clustering.** 2018 AOP spectra have atmospheric-correction drift (see `17_year_effect_pcs.R`); 2018 sites are assigned the nearest 2025 cluster by composition (Hellinger) + snow-free DOY similarity, with a confidence tier.
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
| `data/small_reference/label_community_names.csv` | Hand-tuned narratives for the 31 meadow classes |
| `data/derived/training_samples_sites.csv` | Meadow training: one row per (plot, year) with `final_label`, confidence tier, env covariates |
| `data/derived/training_samples_crowns.gpkg` | Meadow crown polygons for QGIS review |
| `data/derived/shrub_training_set.csv` | Shrub training table (N=310) |
| `data/derived/shrub_label_crosswalk.csv` | canonical_binomial → final_label mapping |
| `data/derived/canopy_height.rds` | 1 m CHM at every crown centroid (1325 sites) |
| `data/derived/joint_training_set.csv` | Joint meadow+shrub training (858 sites × 47 classes) |
| `data/derived/joint_training.gpkg` | Same as above with crown geometries + class metadata |
| `data/derived/punch_list.csv` | Class-level summary: training N, recall, predicted prevalence, top confusions, augmentation priority |
| `data/derived/class_summary_table.csv` | Per-class table with detailed IndVal indicator + abundant taxa strings, descriptions, median leverage |
| `data/derived/target_pixels.{csv,gpkg}` | 6,000 stratified meadow pixels for inference |
| `data/derived/aop_classifier_*.{csv,json}` | R → Python handoff (PCA loadings, preprocessing config, index formulas) |
| `data/derived/inference_predictions.csv` | Joint-RF predictions on the 5,354 cloud-extracted pixels |
| `data/derived/inference_pixel_distances.csv` | Mahalanobis distance to every class centroid per pixel |
| `data/derived/novelty_by_class.csv`, `.../novelty_by_hex.gpkg` | Aggregated novelty for spatial review |
| **`data/derived/sampling_priority.gpkg`** | **Primary deliverable** — 5,354 candidate pixels ranked by leverage |
| `data/derived/sampling_priority_top.csv` | Top 10 per class (346 candidate sites for fieldwork planning) |
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
# Also: a local copy of R3D018 1 m landcover (for 22_target_pixels.R)
# and three NEON CHM tifs (for 36_canopy_height.R) - paths hard-coded
# at the top of those scripts; edit as needed for your environment.

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

- **Inference is meadow-biased.** The 5,354-pixel inference set is filtered to R3D018 meadow class with strict neighborhood-purity rules; shrub leverage is undercounted. A shrub-targeted inference run is planned to balance the picture.
- **Small-N classes need spatial diversity, not just more samples.** Prunus virginiana (n=3), Purshia tridentata (n=4), Alnus incana (n=4) are concentrated in a single drainage in the training data — augmentation should prioritize new geographic locations, not just more N.
- **ESS-DIVE DOIs and 2018↔2025 site crosswalk** still pending in `data/README.md`.
- **2018 wavelength drift** is implicit (nearest-index match, max ~2 nm). Resampling onto a common grid is on the list if a band-specific artifact shows up downstream.

## License

Internal CHESS-team use; no LICENSE file.
