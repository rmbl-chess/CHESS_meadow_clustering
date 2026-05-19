# Sampling-Priority Guide

`data/derived/sampling_priority.gpkg` ranks 5,354 candidate AOP pixels by how much each one would help the meadow / shrub classifier improve if a field crew went out and labeled it. This document explains what the labels mean, how the leverage score is built, and how to use the file in the field.

## What's in the gpkg

| Column | Meaning |
|---|---|
| `predicted_label` | Class the joint RF classifier assigns to this pixel. |
| `class_description` | Human-readable name (meadow narrative or shrub binomial). |
| `nearest_class` | Closest training centroid in feature space — may differ from `predicted_label`. |
| `nearest_d`, `second_d`, `margin` | Mahalanobis distance to nearest + second-nearest class, and the gap. Big `nearest_d` = novel; small `margin` = classifier is on the fence. |
| `ood_flag` | `TRUE` when `nearest_d` exceeds the training 95th-percentile (most pixels are flagged — see caveats). |
| `n_total` | Number of training sites for the predicted class. |
| `balanced_recall` | Class-weighted RF CV recall for that class. |
| `predicted_n_pixels` | How many of the 5,354 inference pixels the classifier assigns to this class. |
| `leverage` | Composite priority score (see below). |
| `augmentation_priority` | Bucket: `critical` / `high` / `medium` / `ok`. |
| `snow_free_doy`, `canopy_height_m` | Site covariates pulled at the pixel. |

## How training classes were built

**Meadow classes (S01–S26 + M01–M05; 31 total).** Spectra-first hierarchical clustering on per-plot 2025 NEON AOP spectra, with PCs 2–12 and snow-free DOY z-scaled (PC1 dropped — too much overall-brightness leverage). 2018 plots were dropped from the clustering itself (atmospheric-correction drift between years) and labeled by inferring the nearest 2025 cluster from composition (Hellinger distance over the harmonized cover table) and phenology. A post-hoc override carves out five monotypic-species classes (M01–M05) wherever a single species exceeded 70% cover — these are spectrally distinctive enough to deserve their own labels. Curated narratives live in `data/small_reference/label_community_names.csv`.

**Shrub classes (16 total).** 2025 from the field shrub-crown table (one species per site by design). 2018 from `fractional_cover` filtered to rows where a shrub-dominated genus had cover = 100%. Synonyms reconciled (`Pentaphylloides floribunda → Dasiphora fruticosa`; `Distegia involucrata → Lonicera involucrata`). Salix species are not spectrally separable from each other; the three highest-N species (`Salix wolfii`, `boothii`, `planifolia`) are kept distinct and the remaining 9 binomials roll up to `Salix other`. Genera with <3 sites are dropped.

**Joint training set.** 858 sites × 47 classes (548 meadow + 310 shrub). Shrub spectra are projected onto the deployed meadow PCA basis (in `aop_classifier_pca.csv`) so training and inference share one feature space. The classifier sees 28 features: 20 PCs, 6 narrow-band indices (NDVI / NDWI / PRI / red-edge slope / CAI / NDLI), snow-free DOY, and canopy height. Random Forest with inverse-frequency class weights, 5-fold site-level CV: 63–65% overall accuracy.

## How leverage is computed

For every inference pixel:

```
leverage = nearest_d / sqrt(n_training_for_predicted_class)
```

- **`nearest_d`** is the Mahalanobis distance from the pixel's feature vector to the closest training class centroid, computed with pooled within-class covariance (z-scaled features). Higher = the pixel is further from anything in the training set.
- **`n_training`** is the number of training sites for the class the RF predicts. Lower = the classifier is leaning on a thin sample.

Multiplying these two signals captures the kind of pixel a new field plot helps the most with: spectrally novel *and* assigned to an under-trained class. The marginal value of one new sample falls as roughly `1/√n` for many learners, hence the square-root denominator. Pixels with high leverage are also where the classifier is most likely to be wrong on the final map.

## How to use it in QGIS / in the field

1. **Load** `sampling_priority.gpkg` in QGIS.
2. **Symbolize** the points: graduated colors on `leverage` (quantile bins). Highest-leverage pixels jump out visually.
3. **Filter** by `augmentation_priority IN ('critical', 'high')` to focus on the urgent classes first.
4. **Label** points by `class_description` so the predicted vegetation type shows on the map.
5. **Plan routes**: `data/derived/sampling_priority_top.csv` already holds the top 10 leverage pixels per predicted class (346 sites across 42 classes). Filter by `domain` and `augmentation_priority` to scope a single field day.
6. **Cross-check** with `data/derived/punch_list.csv` for a class-level summary (training N, recall, predicted area, top confusions) before heading out.

## What the priority buckets mean

- **`critical`** — class has fewer than 5 training sites or zero CV recall. Any new sample helps; geographically clustered training is a real risk.
- **`high`** — fewer than 10 sites or recall below 0.4, OR a class with predicted area ≥ 200 pixels but training N < 20 (high-leverage situations).
- **`medium`** — fewer than 20 sites or recall below 0.6.
- **`ok`** — well-trained classes with adequate recall.

## Caveats

- The 5,354 inference pixels were drawn from R3D018 landcover class 3 (meadow) with strict neighborhood-purity filters. **Shrub leverage is therefore undercounted** in this dataset because shrub crowns are rare in the meadow sample. For shrub-specific priorities, work from `punch_list.csv` directly or generate a shrub-targeted inference set.
- ~94% of inference pixels exceed the OOD threshold. Hand-picked field crowns are spectrally tighter than random basin pixels — the threshold is calibrated against training distribution, so OOD is genuinely common. Treat `nearest_d` as a continuous ranking rather than the binary `ood_flag`.
- Predicted classes for high-`nearest_d` pixels are extrapolations: the classifier picked the closest known class, but it might be a class that doesn't even exist in the training set. New samples in these regions sometimes warrant a new class.
- Salix is at genus-plus-4 granularity by design. Don't expect a `Salix drummondiana` candidate to recover species-level resolution within Salix — that's not currently mappable.
- Some classes (Prunus virginiana, Purshia tridentata, Alnus incana) have N ≤ 4 *and* are concentrated in a single sampling area in the training data. Marked `critical`; sampling elsewhere in the basin is the priority.

## Sources

| File | Role |
|---|---|
| `code/01_load.R` → `code/22_target_pixels.R` | Meadow class pipeline + inference target selection |
| `code/30_shrub_load.R` → `code/35_shrub_pixel_training.R` | Shrub class pipeline |
| `code/36_canopy_height.R` | NEON 1 m CHM extracted at crown centroids |
| `code/37_joint_training.R` | Joint training set + 47-class RF |
| `code/38_predict_inference_pixels.R` | RF predictions on the 5,354 inference pixels |
| `code/39_landscape_distance.R` | Per-pixel Mahalanobis distance to class centroids |
| `code/40_sampling_priority.R` | The leverage score that drives this gpkg |
