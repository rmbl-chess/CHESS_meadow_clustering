# 38_predict_inference_pixels.R — score the 5354-pixel Python-extracted
# inference set with the joint classifier, then enrich the punch list with
# predicted-prevalence info. The premise: classes the joint RF predicts
# OFTEN in the basin but were under-trained are the highest-leverage
# augmentation targets — every new field sample buys real map quality.
#
# Steps:
#   1. Read the extracted parquet (5354 pixels x 32 cols: 20 PCs + 6 indices
#      + snow_free_doy + meta). Drop the parquet pixels that have no CHM
#      coverage (rare edge tiles).
#   2. Extract CHM mean + p90 at each pixel from the local CHM rasters
#      using a 3x3 m window so the value matches the training-side recipe
#      (training-side CHM is averaged over each crown polygon).
#   3. Fit the joint RF on the full meadow + shrub training set (no CV —
#      we want one final model).
#   4. Predict class for every inference pixel; tally counts per class.
#   5. Augment punch_list.csv with `predicted_n_pixels` and a refined
#      `augmentation_priority` that boosts classes that are both
#      under-trained AND frequently predicted.
#
# Inputs:
#   data/derived/extract_meadow_spectra_6k_v1.csv      (parquet pre-converted
#                                                       in Python; arrow R
#                                                       package is not in
#                                                       the renv)
#   data/derived/joint_training_set.rds
#   data/derived/punch_list.csv
#   Local CHM rasters (Google Drive)
# Outputs:
#   data/derived/inference_predictions.parquet
#   data/derived/punch_list.csv (updated in place; previous version
#                                preserved in punch_list_v1.csv)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ranger)
  library(terra)
})

chm_dir <- "/Users/ian/Library/CloudStorage/GoogleDrive-ibreckhe@gmail.com/My Drive/BreckheimerLab2025/Projects/CHESS/Data/AOP_mosaics/NEON_delivered"
chm_paths <- list(
  ALMO = file.path(chm_dir, "CHESS25_ALMO_CHM_1m_v1.tif"),
  CRBU = file.path(chm_dir, "CHESS25_CRBU_CHM_1m_v2.tif"),
  UPTA = file.path(chm_dir, "CHESS25_UPTA_CHM_1m_v1.tif")
)

# --- 1. Read inference pixels --------------------------------------------
pix <- readr::read_csv("data/derived/extract_meadow_spectra_6k_v1.csv",
                       show_col_types = FALSE)
cat(sprintf("Inference pixels: %d (ALMO=%d, CRBU=%d, UPTA=%d)\n",
            nrow(pix),
            sum(pix$domain == "ALMO"),
            sum(pix$domain == "CRBU"),
            sum(pix$domain == "UPTA")))

# --- 1b. Refresh the 20 PCs from the CURRENT PC mosaics ------------------
# CRITICAL basis fix: the cached 6k extraction was projected on a STALE PCA
# basis (sign-flipped vs the current aop_classifier_pca.csv the training PCs +
# inference mosaics use; co-located 6k-vs-mosaic PC correlations were ~ -1.0).
# Re-read the 20 PCs straight from the domain mosaics so the sampling-priority
# pixels share the classifier's basis. The 6 narrow-band indices are computed
# from raw reflectance (basis-independent) and kept as-is. Without this, the
# RF predicts a different basis from the one it was trained on (~8% agreement).
pc_mosaics <- c(ALMO = "ALMO_pc_mosaic.tif", CRBU = "CRBU_pc_mosaic_2025.tif",
                UPTA = "UPTA_pc_mosaic.tif")
pc_cols20 <- sprintf("spec_PC%02d", 1:20)
for (dom in names(pc_mosaics)) {
  idx <- which(pix$domain == dom); if (!length(idx)) next
  mp <- file.path("data/derived/aop_pc_maps_mosaic", pc_mosaics[[dom]])
  if (!file.exists(mp)) stop("missing PC mosaic for basis refresh: ", mp)
  v <- terra::extract(terra::rast(mp), cbind(pix$x_utm[idx], pix$y_utm[idx]))
  for (k in 1:20) pix[[pc_cols20[k]]][idx] <- v[[k]]
  cat(sprintf("  %s: refreshed 20 PCs from mosaic at %d pixels\n", dom, length(idx)))
}
na_pc <- !stats::complete.cases(pix[, pc_cols20])
if (any(na_pc)) {
  cat(sprintf("  dropping %d pixels outside mosaic PC coverage\n", sum(na_pc)))
  pix <- pix[!na_pc, ]
}

# --- 2. CHM at each pixel ------------------------------------------------
# Single-point extraction matches the 1 m CHM resolution. The inference
# pixels are 3 m cell centers, so a single-pixel read is the natural CHM
# value at that center (per-pixel; classifier sees per-3m).
terra::terraOptions(progress = 0)
extract_chm_points <- function(domain, df) {
  chm  <- terra::rast(chm_paths[[domain]])
  vect <- terra::vect(df, geom = c("x_utm", "y_utm"), crs = "EPSG:32613")
  cat(sprintf("  %s: extracting CHM at %d pixels ... ", domain, nrow(df)))
  t0 <- Sys.time()
  vals <- terra::extract(chm, vect, ID = FALSE)
  cat(sprintf("done (%.1fs)\n",
              as.numeric(Sys.time() - t0, units = "secs")))
  df |> dplyr::mutate(canopy_height_m = vals[[1]])
}

pix_full <- purrr::map_dfr(c("ALMO", "CRBU", "UPTA"), function(dom) {
  sub <- dplyr::filter(pix, domain == dom)
  if (nrow(sub) == 0) return(NULL)
  extract_chm_points(dom, sub)
}) |>
  dplyr::filter(!is.na(canopy_height_m))
cat(sprintf("Pixels with CHM coverage: %d / %d\n", nrow(pix_full), nrow(pix)))

# --- 3. Fit joint RF on full training set -------------------------------
js <- readRDS("data/derived/joint_training_set.rds")
train <- js$training
features <- js$feature_cols
X_tr <- as.matrix(train[, features])
y_tr <- factor(train$final_label)
cat(sprintf("\nFitting joint RF on %d training sites, %d features, %d classes\n",
            nrow(X_tr), ncol(X_tr), nlevels(y_tr)))

# Important: fit UNWEIGHTED for inference. Balanced class weights
# (1/freq) inflate small classes ~20x; in CV that buys per-class recall
# but in inference it turns the smallest classes into catch-all
# predictions for any pixel the model is unsure about. The balanced
# recall metric is still tracked in 37 (the punch list); 38 is only
# for the realistic class-proportion map.
fit_joint <- ranger::ranger(
  x = X_tr, y = y_tr, num.trees = 1000,
  classification = TRUE, seed = 42
)
cat(sprintf("RF OOB prediction error (unweighted): %.1f%%\n",
            100 * fit_joint$prediction.error))

# --- 4. Predict for inference pixels ------------------------------------
X_inf <- as.matrix(pix_full[, features])
pred  <- predict(fit_joint, X_inf)$predictions
pix_full$predicted_label <- as.character(pred)
cat(sprintf("\nPredicted class distribution across %d inference pixels:\n",
            nrow(pix_full)))
print(as.data.frame(pix_full |>
        dplyr::count(predicted_label) |>
        dplyr::arrange(dplyr::desc(n))))

# --- 5. Refresh punch list with predicted-prevalence info ----------------
punch_old <- readr::read_csv("data/derived/punch_list.csv",
                             show_col_types = FALSE)
if (!file.exists("data/derived/punch_list_v1.csv")) {
  readr::write_csv(punch_old, "data/derived/punch_list_v1.csv")
  cat("Snapshotted previous punch list as punch_list_v1.csv\n")
}
pred_counts <- pix_full |>
  dplyr::count(predicted_label, name = "predicted_n_pixels") |>
  dplyr::rename(final_label = predicted_label)

# Drop any old prediction columns so a re-run doesn't create
# predicted_n_pixels.x / .y after the join. The new column comes from
# pred_counts below.
punch_old <- punch_old |>
  dplyr::select(-dplyr::any_of(c("predicted_n_pixels", "pixels_per_site")))

# Pixels per training site = a rough proxy for how much real-estate each
# new field sample would buy in the predicted map.
punch <- punch_old |>
  dplyr::left_join(pred_counts, by = "final_label") |>
  dplyr::mutate(
    predicted_n_pixels = dplyr::coalesce(predicted_n_pixels, 0L),
    pixels_per_site    = predicted_n_pixels / pmax(n_total, 1L),
    augmentation_priority = dplyr::case_when(
      n_total < 5 | balanced_recall == 0                           ~ "critical",
      predicted_n_pixels >= 200 & n_total < 20                     ~ "high",
      n_total < 10 | balanced_recall < 0.4                         ~ "high",
      n_total < 20 | balanced_recall < 0.6                         ~ "medium",
      TRUE                                                          ~ "ok"
    )
  ) |>
  dplyr::arrange(factor(augmentation_priority,
                        levels = c("critical", "high", "medium", "ok")),
                 dplyr::desc(pixels_per_site))

cat("\n=== Updated punch list ===\n")
print(as.data.frame(punch))

readr::write_csv(punch, "data/derived/punch_list.csv")
readr::write_csv(pix_full, "data/derived/inference_predictions.csv")
cat("\nWrote data/derived/punch_list.csv\n")
cat("Wrote data/derived/inference_predictions.csv\n")
