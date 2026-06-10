# Sampling-Priority Guide

The 2026 CHESS sampling priority dataset ranks 5,064 candidate AOP pixels by how much each one would help the meadow / shrub classifier improve if a field crew went out and labeled it (the 5,354-pixel cloud extraction, minus 290 pixels without CHM coverage). This document explains what the labels mean, how the leverage score is built, and how to use the file in the field. Alongside it, the classifier now produces **wall-to-wall 3 m community maps** for the three AOP domains (ALMO / CRBU / UPTA) — see Files.

## Files for the team

Canonical outputs live in the analysis repo under `data/derived/`:

| File | What it is |
|---|---|
| `sampling_priority.gpkg` | **Primary fieldwork deliverable** — 5,064 candidate pixels with leverage scores, predicted class, and novelty metrics. |
| `class_summary_table.csv` | One row per class (**51 classes**: 34 meadow + 17 shrub) with 2018/2025/2026 support, basin prevalence, recall, indicator + abundant taxa, narratives, median leverage, and augmentation priority. |
| `joint_training.gpkg` | All training sites (604 meadow + 323 shrub) as crown polygons + centroids with class metadata — for inspecting where each class is currently sampled. |
| `aop_classified/{DOM}_class_3m_v1_labeled.tif` (+ `.qml`) | **Wall-to-wall 3 m community map** per AOP domain (ALMO / CRBU / UPTA). Each class is colored by physiognomy family (shrubland blue-purple, grassland tan-brown, forb-meadow green, wetland pink-red) and labeled with its code + a draft NatureServe community; `{DOM}_confidence_3m_v1.tif` holds per-pixel max class probability. Load the `_labeled.tif` — the sibling `.qml` auto-applies colors + labels. **⚠️ These COGs are still on the prior 40/56-class system — pending re-inference on the current 34/51-class classifier (`09/10` deferred).** |

Export dated snapshots to the team Google Drive (`.../Sampling_Priority_2026/datasets/`) for fieldwork. **Re-export from the current 51-class outputs (the 2026 supplemental ALMO sagebrush campaign reclustered the meadow set from 40 → 34 classes).**

## What's in the gpkg

| Column | Meaning |
|---|---|
| `predicted_label` | Class the joint RF classifier assigns to this pixel. |
| `class_description` | Human-readable name (meadow narrative or shrub binomial). |
| `nearest_class` | Closest training centroid in feature space — may differ from `predicted_label`. |
| `nearest_d`, `second_d`, `margin` | Mahalanobis distance to nearest + second-nearest class, and the gap. Big `nearest_d` = novel; small `margin` = classifier is on the fence. |
| `ood_flag` | `TRUE` when `nearest_d` exceeds the training 95th-percentile (most pixels are flagged — see caveats). |
| `n_total` | Number of training sites for the predicted class. |
| `balanced_recall` | Class-weighted RF CV recall for that class (from `02_training.R`). |
| `predicted_n_pixels` | How many of the 5,064 inference pixels the unweighted RF assigns to this class (from `03_predict_inference_pixels.R`). |
| `leverage` | Composite priority score (see below). |
| `augmentation_priority` | Bucket: `critical` / `high` / `medium` / `ok`. |
| `snow_free_doy`, `canopy_height_m` | Site covariates pulled at the pixel. |

## How training classes were built

**Meadow classes (34 total).** Spectra-first hierarchical Ward clustering on per-plot NEON AOP spectra (PCs 2–12 + snow-free DOY, z-scaled; PC1 dropped). **2018, 2025, and 2026** plots cluster together — an NDVI-stratified 2018→2025 radiometric correction removes the between-year drift (the 2026 supplemental plots were extracted directly from 2025 AOP, so they need no correction). Each spectral cluster is then split into community sub-classes by species composition, but only where the split stays **spectrally mappable** (an RF-recall gate); ecologically real but hard-to-map sub-classes are kept and flagged `needs_ancillary`. In the current set only S06 (k=2) and S08 (k=3) split, and none are flagged `needs_ancillary` — the 56 new ALMO sagebrush plots reshaped the partition so the prior S01/S09/S20 sub-splits no longer form. Monotypic-species overrides carve out near-pure stands (≥70% cover) as their own classes. Curatable narratives live in `data/small_reference/label_community_names.csv`.

**Shrub classes (17 total).** 2025 from the field shrub-crown table (one species per site by design). 2018 from `fractional_cover` filtered to rows where a shrub-dominated genus had cover = 100%. Synonyms reconciled (`Pentaphylloides floribunda → Dasiphora fruticosa`; `Distegia involucrata → Lonicera involucrata`). Salix species are not spectrally separable from each other; the three highest-N species (`Salix wolfii`, `boothii`, `planifolia`) are kept distinct and the remaining 9 binomials roll up to `Salix other`. Genera with <3 sites are dropped.

**Joint training set.** 927 sites × 51 classes (604 meadow + 323 shrub). Shrub spectra are projected onto the deployed meadow PCA basis (`aop_classifier_pca.csv`) so training and inference share one feature space. CV uses 28 features (20 PCs, 6 narrow-band indices — NDVI / NDWI / PRI / red-edge slope / CAI / NDLI — snow-free DOY, canopy height); the **deployed inference RF uses 22**, dropping the 6 indices, which aren't available in the PC-only inference mosaics.

**Two Random Forest fits, two purposes.** The same training data is used in two RFs that answer different questions:

- **Balanced-weighted RF (`02_training.R`, 5-fold site-level CV)** — inverse-frequency class weights so every class contributes equally to per-class recall. Source of `balanced_recall`. Overall CV accuracy ~67% (0.674 unweighted / 0.668 balanced).
- **Unweighted RF (`09_inference.R`, predicts the basin pixels AND the wall-to-wall map)** — no class weights, so predictions follow the natural class proportions. Source of `predicted_label` and `predicted_n_pixels`. Using weights here was tried first and produced a wildly biased map: 33% of pixels assigned to the smallest meadow class (S26, n=7), because inverse-frequency weighting amplified rare classes ~22× and turned them into catch-alls for any uncertain pixel.

The leverage score combines `balanced_recall`-style class quality with `predicted_n_pixels`-style prevalence, so both RFs feed into the prioritization.

## How leverage is computed

For every inference pixel:

```
leverage = nearest_d / sqrt(n_training_for_predicted_class)
```

- **`nearest_d`** is the Mahalanobis distance from the pixel's feature vector to the closest training class centroid, computed with pooled within-class covariance (z-scaled features). Higher = the pixel is further from anything in the training set.
- **`n_training`** is the number of training sites for the class the RF predicts. Lower = the classifier is leaning on a thin sample.

Multiplying these two signals captures the kind of pixel a new field plot helps the most with: spectrally novel *and* assigned to an under-trained class. The marginal value of one new sample falls as roughly `1/√n` for many learners, hence the square-root denominator. Pixels with high leverage are also where the classifier is most likely to be wrong on the final map.

## How to use it in QGIS / in the field

1. **Load** `Meadow_shrub_sampling_priority_2026_5_19.gpkg` in QGIS.
2. **Symbolize** the points: graduated colors on `leverage` (quantile bins). Highest-leverage pixels jump out visually.
3. **Filter** by `augmentation_priority IN ('critical', 'high')` to focus on the urgent classes first.
4. **Label** points by `class_description` so the predicted vegetation type shows on the map.
5. **Cross-check** with `Meadow_shrub_class_summary_table_2026_5_19.csv` for a class-level summary (training N, recall, predicted area, indicator taxa, top confusions) before heading out.
6. **See existing training** by overlaying `Meadow_shrub_joint_training_2026_05_19.gpkg` (the `training_sites_points` layer) so you can avoid re-sampling well-covered classes.

## What the priority buckets mean

- **`critical`** — class has fewer than 5 training sites or zero CV recall. Any new sample helps; geographically clustered training is a real risk.
- **`high`** — fewer than 10 sites or recall below 0.4, OR a class with predicted area ≥ 200 pixels but training N < 20 (high-leverage situations).
- **`medium`** — fewer than 20 sites or recall below 0.6.
- **`ok`** — well-trained classes with adequate recall.

## Caveats

- `predicted_label` reflects the *unweighted* RF (script 38) — it shows the realistic class proportions for the basin. `balanced_recall` reflects the *weighted* RF (script 37) — it shows per-class quality at training time. Both are correct for their own purpose; just don't expect a class with high `balanced_recall` to also have high `predicted_n_pixels`.
- The inference pixels were drawn from R3D018 landcover class 3 (meadow) with strict neighborhood-purity filters. **Shrub leverage is therefore undercounted** in this dataset because shrub crowns are rare in the meadow sample. For shrub-specific priorities, work from the class summary table directly (filter `class_type == "shrub"`) — a separate shrub-targeted inference run is planned.
- ~94% of inference pixels exceed the OOD threshold. Hand-picked field crowns are spectrally tighter than random basin pixels — the threshold is calibrated against training distribution, so OOD is genuinely common. Treat `nearest_d` as a continuous ranking rather than the binary `ood_flag`.
- Predicted classes for high-`nearest_d` pixels are extrapolations: the classifier picked the closest known class, but it might be a class that doesn't even exist in the training set. New samples in these regions sometimes warrant a new class.
- Salix is at genus-plus-4 granularity by design. Don't expect a `Salix drummondiana` candidate to recover species-level resolution within Salix — that's not currently mappable.
- Some classes (Prunus virginiana, Purshia tridentata, Alnus incana) have N ≤ 4 *and* are concentrated in a single sampling area in the training data. Marked `critical`; sampling elsewhere in the basin is the priority.

## Sources

| File | Role |
|---|---|
| `code/meadow/01_load.R` → `code/meadow/23_target_pixels.R` | Meadow class pipeline + inference target selection |
| `code/shrub/01_load.R` → `code/shrub/06_pixel_training.R` | Shrub class pipeline |
| `code/joint/01_canopy_height.R` | NEON 1 m CHM extracted at crown centroids |
| `code/joint/02_training.R` | Joint training set + 47-class RF |
| `code/joint/03_predict_inference_pixels.R` | RF predictions on the 5,064 inference pixels |
| `code/joint/04_landscape_distance.R` | Per-pixel Mahalanobis distance to class centroids |
| `code/joint/05_sampling_priority.R` | The leverage score that drives this gpkg |
