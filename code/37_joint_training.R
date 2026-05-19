# 37_joint_training.R — combine the meadow + shrub training sets, project
# shrub spectra onto the canonical meadow PCA basis, add CHM as a covariate,
# fit a joint Random-Forest classifier, and emit a punch list of classes
# whose training data is the most likely bottleneck on final-map quality.
#
# The deployed AOP-tile classifier (extract_aop_features.py +
# generate_aop_pc_maps.py) is hard-wired to the meadow PCA basis already
# saved in aop_classifier_pca.csv. To keep training and inference on the
# same basis, shrub spectra are projected through that same basis here.
#
# Inputs:
#   data/derived/training_samples_sites.csv     (548 meadow sites x features)
#   data/derived/shrub_veg_spectra.rds          (310 shrub sites, raw rfl bands)
#   data/derived/shrub_training_set.rds         (canonical -> final_label crosswalk)
#   data/derived/spectral_features.rds          (meadow PCA: rotation, center,
#                                                 keep_cols, keep_wl)
#   data/derived/environment.rds                (snow_free_doy per site)
#   data/derived/canopy_height.rds              (CHM mean + p90 per site)
# Outputs:
#   data/derived/joint_training_set.rds     (training table, recall, confusion)
#   data/derived/joint_training_set.csv     committed for review/handoff
#   data/derived/punch_list.csv             committed; augmentation targets
#   output/joint_confusion.png

suppressPackageStartupMessages({
  library(tidyverse)
  library(ranger)
  library(ggplot2)
})

# --- 1. Load meadow training ---------------------------------------------
meadow <- readr::read_csv("data/derived/training_samples_sites.csv",
                          show_col_types = FALSE)
cat(sprintf("Meadow training: %d sites, %d classes\n",
            nrow(meadow), dplyr::n_distinct(meadow$final_label)))

# --- 2. Load and project shrub spectra onto meadow PCA basis -------------
sf  <- readRDS("data/derived/spectral_features.rds")
pca <- sf$pca                       # rotation (348 x 348), center
keep_wl <- sf$keep_wl               # wavelengths kept after water mask (nm)

shrub_vs    <- readRDS("data/derived/shrub_veg_spectra.rds")
shrub_train <- readRDS("data/derived/shrub_training_set.rds")$training |>
  dplyr::select(site_number, Year, canonical_binomial, final_label)

# Inner-join: only shrub sites that survived the crosswalk (have a final_label).
shrub_features_raw <- shrub_vs$joined |>
  dplyr::inner_join(shrub_train, by = c("site_number", "Year")) |>
  dplyr::filter(!is.na(final_label))

# Reconstruct keep_cols: pick the rfl_band_* columns whose wavelengths match
# the meadow PCA's keep_wl (within tolerance). Order matches keep_wl.
rfl_cols  <- grep("^rfl_band_", names(shrub_features_raw), value = TRUE)
band_nums <- as.integer(stringr::str_extract(rfl_cols, "\\d+$"))
band_wl_all <- shrub_vs$wavelengths$center_wavelength_nm[
                 match(band_nums, shrub_vs$wavelengths$band_number)]
keep_idx  <- vapply(keep_wl,
                    function(w) which.min(abs(band_wl_all - w)),
                    integer(1))
keep_cols <- rfl_cols[keep_idx]
stopifnot(length(keep_cols) == length(keep_wl))
cat(sprintf("Matched %d shrub bands to meadow PCA basis (max drift %.2f nm)\n",
            length(keep_cols),
            max(abs(band_wl_all[keep_idx] - keep_wl))))

spec_mat <- as.matrix(shrub_features_raw[, keep_cols])
spec_ctr <- sweep(spec_mat, 2, pca$center, FUN = "-")
spec_pcs <- spec_ctr %*% pca$rotation[, seq_len(20)]
colnames(spec_pcs) <- sprintf("spec_PC%02d", seq_len(20))

# --- 3. Compute the same six narrow-band indices for shrub spectra -------
b_at <- function(target_nm) {
  spec_mat[, which.min(abs(keep_wl - target_nm))]
}
wl_at <- function(target_nm) keep_wl[which.min(abs(keep_wl - target_nm))]

shrub_indices <- tibble::tibble(
  ndvi  = (b_at(860) - b_at(660)) / (b_at(860) + b_at(660)),
  ndwi  = (b_at(860) - b_at(1240)) / (b_at(860) + b_at(1240)),
  pri   = (b_at(531) - b_at(570)) / (b_at(531) + b_at(570)),
  red_edge_slope = (b_at(750) - b_at(700)) / (wl_at(750) - wl_at(700)),
  cai   = 0.5 * (b_at(2000) + b_at(2200)) - b_at(2100),
  ndli  = (log(1 / b_at(1754)) - log(1 / b_at(1680))) /
          (log(1 / b_at(1754)) + log(1 / b_at(1680)))
)

# --- 4. Pull DOY for shrub sites ----------------------------------------
env <- readRDS("data/derived/environment.rds")

shrub <- shrub_features_raw |>
  dplyr::select(site_number, Year, sampling_area, final_label) |>
  dplyr::bind_cols(as.data.frame(spec_pcs), shrub_indices) |>
  dplyr::left_join(env, by = c("site_number", "Year"))

cat(sprintf("Shrub training (projected onto meadow PCA): %d sites, %d classes\n",
            nrow(shrub), dplyr::n_distinct(shrub$final_label)))

# --- 5. Combine meadow + shrub into one table ----------------------------
common_cols <- c("site_number", "Year", "final_label",
                 paste0("spec_PC", sprintf("%02d", 1:20)),
                 "ndvi", "ndwi", "pri", "red_edge_slope", "cai", "ndli",
                 "snow_free_doy")
meadow_slim <- meadow |>
  dplyr::mutate(class_type = "meadow") |>
  dplyr::select(dplyr::any_of(c(common_cols, "class_type")))
shrub_slim <- shrub |>
  dplyr::mutate(class_type = "shrub") |>
  dplyr::select(dplyr::any_of(c(common_cols, "class_type")))
joint <- dplyr::bind_rows(meadow_slim, shrub_slim)
cat(sprintf("\nJoint training set: %d sites (%d meadow + %d shrub), %d classes\n",
            nrow(joint), sum(joint$class_type == "meadow"),
            sum(joint$class_type == "shrub"),
            dplyr::n_distinct(joint$final_label)))

# --- 6. Add canopy height covariate (CHM) -------------------------------
chm <- readRDS("data/derived/canopy_height.rds")
joint <- joint |>
  dplyr::left_join(chm, by = c("site_number", "Year"))
cat(sprintf("Sites with CHM: %d / %d (drop sites with no canopy data)\n",
            sum(!is.na(joint$canopy_height_m)), nrow(joint)))
joint_full <- joint |> dplyr::filter(!is.na(canopy_height_m),
                                     !is.na(snow_free_doy))
cat(sprintf("Joint training (complete-case): %d sites\n", nrow(joint_full)))

# --- 7. RF stratified site-level CV --------------------------------------
feature_cols <- c(
  paste0("spec_PC", sprintf("%02d", 1:20)),
  "ndvi", "ndwi", "pri", "red_edge_slope", "cai", "ndli",
  "snow_free_doy", "canopy_height_m"
)
X <- as.matrix(joint_full[, feature_cols])
y <- factor(joint_full$final_label)
n_folds <- 5

set.seed(42)
fold <- integer(nrow(X))
for (lvl in levels(y)) {
  idx <- which(y == lvl)
  fold[idx] <- ((sample(seq_along(idx)) - 1L) %% n_folds) + 1L
}

run_cv <- function(use_weights) {
  preds <- factor(rep(NA_character_, nrow(X)), levels = levels(y))
  for (f in seq_len(n_folds)) {
    tr <- which(fold != f); te <- which(fold == f)
    if (length(tr) == 0 || length(te) == 0) next
    cw <- NULL
    if (use_weights) {
      tab <- table(y[tr])
      cw  <- as.numeric(sum(tab) / (length(tab) * tab))
      names(cw) <- names(tab)
    }
    fit <- ranger::ranger(
      x = X[tr, , drop = FALSE], y = y[tr],
      num.trees = 800, classification = TRUE,
      class.weights = cw, seed = 42 + f
    )
    preds[te] <- predict(fit, X[te, , drop = FALSE])$predictions
  }
  preds
}

cat("\nRunning joint RF CV (this may take a minute) ...\n")
preds_unw <- run_cv(FALSE)
preds_w   <- run_cv(TRUE)

summarise_recall <- function(preds, name) {
  tibble::tibble(truth = y, pred = preds) |>
    dplyr::filter(!is.na(pred)) |>
    dplyr::group_by(truth) |>
    dplyr::summarise(n_test = dplyr::n(),
                     recall = mean(pred == truth), .groups = "drop") |>
    dplyr::rename(final_label = truth) |>
    dplyr::mutate(model = name)
}
recall_df <- dplyr::bind_rows(
  summarise_recall(preds_unw, "unweighted"),
  summarise_recall(preds_w,   "balanced")
) |>
  dplyr::left_join(joint_full |> dplyr::count(final_label, name = "n"),
                   by = "final_label") |>
  dplyr::mutate(final_label = as.character(final_label)) |>
  tidyr::pivot_wider(id_cols = c(final_label, n),
                     names_from = model,
                     values_from = recall,
                     names_glue = "{model}_recall")

overall_unw <- mean(preds_unw == y, na.rm = TRUE)
overall_w   <- mean(preds_w   == y, na.rm = TRUE)
cat(sprintf("\nJoint RF CV overall accuracy:  unweighted = %.1f%%   balanced = %.1f%%\n",
            100 * overall_unw, 100 * overall_w))

# --- 8. Confusion analysis: top off-diagonal pairs per class -------------
conf <- tibble::tibble(truth = y, pred = preds_w) |>
  dplyr::filter(!is.na(pred), truth != pred) |>
  dplyr::count(truth, pred, name = "n_confused")
top_confusions <- conf |>
  dplyr::group_by(truth) |>
  dplyr::arrange(dplyr::desc(n_confused), .by_group = TRUE) |>
  dplyr::slice_head(n = 2) |>
  dplyr::summarise(top_confusions = paste(sprintf("%s (%d)", pred, n_confused),
                                          collapse = "; "),
                   .groups = "drop") |>
  dplyr::rename(final_label = truth) |>
  dplyr::mutate(final_label = as.character(final_label))

# --- 9. Build punch list ------------------------------------------------
# Per-class spatial diversity: how many distinct sampling areas does it
# come from? A class with N=10 all from one drainage is far more fragile
# than a class with N=10 from 10 separate areas.
sampling_area_diversity <- function() {
  # Meadow side: sampling_area in training_samples_sites.csv -> not present
  # in the current export. Best proxy is unique site_number per Year combos,
  # which is just N (no signal). For meadow we fall back to a stub of 'NA'
  # and rely on the shrub-side detail. Future improvement: export sampling
  # area in the meadow training set.
  shrub_div <- shrub |>
    dplyr::group_by(final_label) |>
    dplyr::summarise(n_sampling_areas = dplyr::n_distinct(sampling_area),
                     .groups = "drop")
  shrub_div
}
div <- sampling_area_diversity()

punch <- joint_full |>
  dplyr::count(final_label, class_type, name = "n_total") |>
  dplyr::left_join(
    joint_full |> dplyr::count(final_label, Year) |>
      tidyr::pivot_wider(names_from = Year, values_from = n,
                         names_prefix = "n_", values_fill = 0L),
    by = "final_label"
  ) |>
  dplyr::left_join(div, by = "final_label") |>
  dplyr::left_join(recall_df, by = "final_label") |>
  dplyr::left_join(top_confusions, by = "final_label") |>
  dplyr::mutate(
    augmentation_priority = dplyr::case_when(
      n_total < 5 | balanced_recall == 0                        ~ "critical",
      n_total < 10 | (balanced_recall < 0.4)                    ~ "high",
      n_total < 20 | (balanced_recall < 0.6)                    ~ "medium",
      TRUE                                                      ~ "ok"
    )
  ) |>
  dplyr::arrange(factor(augmentation_priority,
                        levels = c("critical", "high", "medium", "ok")),
                 dplyr::desc(n_total))

cat("\n=== Punch list (top of list = most urgent) ===\n")
print(as.data.frame(punch))

# --- 10. Persist --------------------------------------------------------
dir.create("output", showWarnings = FALSE)
readr::write_csv(punch, "data/derived/punch_list.csv")
readr::write_csv(
  joint_full |> dplyr::select(site_number, Year, class_type, final_label,
                              snow_free_doy, canopy_height_m),
  "data/derived/joint_training_set.csv"
)
saveRDS(list(
  training        = joint_full,
  feature_cols    = feature_cols,
  recall          = recall_df,
  confusion       = conf,
  top_confusions  = top_confusions,
  punch_list      = punch,
  overall_acc     = list(unweighted = overall_unw, balanced = overall_w)
), "data/derived/joint_training_set.rds")
cat("\nWrote data/derived/punch_list.csv\n")
cat("Wrote data/derived/joint_training_set.{rds,csv}\n")
