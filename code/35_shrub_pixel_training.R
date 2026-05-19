# 35_shrub_pixel_training.R — train RF on shrub-crown pixels (not the
# site-averaged spectrum), then aggregate pixel-level class probabilities
# back to a per-site prediction. Goal: see whether the extra training
# signal from ~20-30 pixels/site rescues the weak minority classes that
# stayed at zero recall in 34.
#
# Design points:
#   - Folds are assigned at the SITE level (all pixels of one site share
#     a fold) to prevent the obvious within-site leakage that would make
#     a pixel-level CV look great but tell us nothing.
#   - Class weights are inverse-frequency at the PIXEL level (so the Salix
#     mega-class doesn't dominate any worse than at site level).
#   - PCA is re-fit on the pixel-level data — pixel-scale PC structure
#     can differ from site-mean PC structure.
#   - Site prediction = argmax of mean class probability across the site's
#     pixels. Compared to site-averaged baseline from 34.
#
# Inputs:
#   data/derived/shrub_records_canonical.rds
#   data/derived/shrub_label_crosswalk.csv
#   data/derived/spectra_2018.rds, data/derived/spectra_2025.rds
#   data/derived/shrub_training_set.rds   (site-averaged baseline)
# Outputs:
#   data/derived/shrub_pixel_training.rds
#   output/shrub_pixel_vs_site_recall.png
#   output/shrub_pixel_confusion.png

suppressPackageStartupMessages({
  library(tidyverse)
  library(ranger)
  library(ggplot2)
})

records   <- readRDS("data/derived/shrub_records_canonical.rds")
crosswalk <- readr::read_csv("data/derived/shrub_label_crosswalk.csv",
                             show_col_types = FALSE)
sp_2018   <- readRDS("data/derived/spectra_2018.rds")
sp_2025   <- readRDS("data/derived/spectra_2025.rds")
env       <- readRDS("data/derived/environment.rds")           # snow-free DOY
baseline  <- readRDS("data/derived/shrub_training_set.rds")  # site-averaged

# Filter knobs. Pixels are dropped if NDVI < ndvi_min (filters shadows and
# mixed/non-vegetation edges). A site is dropped if too few pixels survive.
ndvi_min <- 0.30
min_pixels_per_site <- 5L

# --- 1. Build per-pixel feature table ------------------------------------
# Map each pixel to its site's final_label via canonical binomial.
labels <- records |>
  dplyr::left_join(crosswalk |> dplyr::select(canonical_binomial, final_label),
                   by = "canonical_binomial") |>
  dplyr::filter(!is.na(final_label)) |>
  dplyr::select(site_number, Year, final_label)

# Filter 2025 to Shrub site_type pixels at shrub sites; 2018 to all pixels
# at shrub-dominated sites (no site_type column in 2018 spectra).
pix_2025 <- sp_2025$spectra |>
  dplyr::filter(site_type == "Shrub",
                shade == 1 | is.na(shade)) |>
  dplyr::semi_join(dplyr::filter(labels, Year == 2025L), by = "site_number") |>
  dplyr::mutate(Year = 2025L)
pix_2018 <- sp_2018$spectra |>
  dplyr::filter(shade == 1 | is.na(shade)) |>
  dplyr::semi_join(dplyr::filter(labels, Year == 2018L), by = "site_number") |>
  dplyr::mutate(Year = 2018L)
cat(sprintf("Pixels: %d (2025 Shrub) + %d (2018 shrub-dominated)\n",
            nrow(pix_2025), nrow(pix_2018)))

pixels <- dplyr::bind_rows(pix_2025, pix_2018)
rfl_cols  <- grep("^rfl_band_", names(pixels), value = TRUE)
band_nums <- as.integer(stringr::str_extract(rfl_cols, "\\d+$"))
band_wl   <- sp_2025$wavelengths$center_wavelength[
                match(band_nums, sp_2025$wavelengths$band_number)] * 1000
water_mask <- (band_wl >= 1340 & band_wl <= 1450) |
              (band_wl >= 1800 & band_wl <= 1950) |
              (band_wl >  2400)
keep_cols <- rfl_cols[!water_mask]
cat(sprintf("Bands retained after water mask: %d (of %d)\n",
            length(keep_cols), length(rfl_cols)))

# --- Per-pixel NDVI filter ------------------------------------------------
# NDVI on raw (non-normalized) reflectance: NIR ~860 nm, RED ~660 nm. The
# ratio is invariant to L2 normalization, so we can compute it on raw and
# apply the threshold before any other processing.
idx_red <- rfl_cols[which.min(abs(band_wl - 660))]
idx_nir <- rfl_cols[which.min(abs(band_wl - 860))]
ndvi <- (pixels[[idx_nir]] - pixels[[idx_red]]) /
        (pixels[[idx_nir]] + pixels[[idx_red]])
cat(sprintf("Pre-NDVI pixels: %d  NDVI distribution: min=%.2f med=%.2f max=%.2f\n",
            nrow(pixels), min(ndvi, na.rm=TRUE),
            median(ndvi, na.rm=TRUE), max(ndvi, na.rm=TRUE)))
pass_ndvi <- !is.na(ndvi) & ndvi >= ndvi_min
cat(sprintf("After NDVI >= %.2f filter: %d pixels (%.1f%% retained)\n",
            ndvi_min, sum(pass_ndvi), 100 * mean(pass_ndvi)))
pixels <- pixels[pass_ndvi, ]

# Per-pixel L2 normalize, then drop rows that are NA or zero-norm.
spec_mat <- as.matrix(pixels[, keep_cols])
norms    <- sqrt(rowSums(spec_mat^2, na.rm = TRUE))
good     <- is.finite(norms) & norms > 0 & complete.cases(spec_mat)
spec_mat <- spec_mat[good, , drop = FALSE]
pixels   <- pixels[good, ]
spec_mat <- spec_mat / sqrt(rowSums(spec_mat^2))
cat(sprintf("Pixels after QC: %d\n", nrow(spec_mat)))

# Drop sites whose surviving pixel count is below the threshold.
site_pix <- pixels |>
  dplyr::group_by(site_number, Year) |>
  dplyr::summarise(n_pix = dplyr::n(), .groups = "drop")
good_sites <- site_pix |> dplyr::filter(n_pix >= min_pixels_per_site)
cat(sprintf("Sites with >= %d surviving pixels: %d / %d\n",
            min_pixels_per_site, nrow(good_sites), nrow(site_pix)))
keep_pix <- paste(pixels$site_number, pixels$Year, sep = "_") %in%
              paste(good_sites$site_number, good_sites$Year, sep = "_")
pixels   <- pixels[keep_pix, ]
spec_mat <- spec_mat[keep_pix, , drop = FALSE]

# Join the site → final_label map AND snow-free DOY onto each pixel.
pixels <- pixels |>
  dplyr::left_join(labels, by = c("site_number", "Year")) |>
  dplyr::left_join(env,    by = c("site_number", "Year")) |>
  dplyr::filter(!is.na(final_label), !is.na(snow_free_doy))
spec_mat <- spec_mat[seq_len(nrow(pixels)), , drop = FALSE]   # alignment guard

cat(sprintf("Pixels with assigned label: %d across %d sites\n",
            nrow(pixels),
            dplyr::n_distinct(pixels[, c("site_number", "Year")])))
cat("\n=== Pixel counts per class ===\n")
print(pixels |> dplyr::count(final_label, name = "n_pixels") |>
        dplyr::arrange(dplyr::desc(n_pixels)) |> as.data.frame())

# --- 2. PCA on pixel-level data -----------------------------------------
pca <- prcomp(spec_mat, center = TRUE, scale. = FALSE)
n_pc <- 20
PCs <- pca$x[, seq_len(n_pc), drop = FALSE]
colnames(PCs) <- sprintf("PC%02d", seq_len(n_pc))
cat(sprintf("\nPixel-PCA: PC1-%d explain %.1f%% of variance\n",
            n_pc, 100 * sum(pca$sdev[seq_len(n_pc)]^2) / sum(pca$sdev^2)))

# Feature matrix = pixel PCs concatenated with the site-level snow-free DOY
# (z-scaled across all pixels — all pixels at the same site share a value).
doy_z <- as.numeric(scale(pixels$snow_free_doy))
X <- cbind(PCs, snow_free_doy_z = doy_z)
y <- factor(pixels$final_label)
site_key <- paste(pixels$site_number, pixels$Year, sep = "_")

# --- 3. Site-level fold assignment --------------------------------------
sites_per_label <- pixels |>
  dplyr::distinct(site_number, Year, final_label) |>
  dplyr::count(final_label, name = "n_sites")
cat("\n=== Sites per class ===\n")
print(sites_per_label |> dplyr::arrange(dplyr::desc(n_sites)) |>
        as.data.frame())

set.seed(42)
site_table <- pixels |>
  dplyr::distinct(site_number, Year, final_label) |>
  dplyr::group_by(final_label) |>
  dplyr::mutate(fold = ((sample(seq_along(site_number)) - 1L) %% 5L) + 1L) |>
  dplyr::ungroup() |>
  dplyr::mutate(site_key = paste(site_number, Year, sep = "_"))

fold <- site_table$fold[match(site_key, site_table$site_key)]
n_folds <- 5

# --- 4. RF CV (probability) with class weights --------------------------
preds_prob <- matrix(NA_real_, nrow = nrow(X), ncol = nlevels(y),
                     dimnames = list(NULL, levels(y)))
for (f in seq_len(n_folds)) {
  tr <- which(fold != f); te <- which(fold == f)
  if (length(tr) == 0 || length(te) == 0) next
  tab <- table(y[tr])
  cw  <- as.numeric(sum(tab) / (length(tab) * tab))
  names(cw) <- names(tab)
  fit <- ranger::ranger(
    x = X[tr, , drop = FALSE], y = y[tr],
    num.trees = 500,
    probability = TRUE,
    class.weights = cw,
    seed = 42 + f
  )
  pp <- predict(fit, X[te, , drop = FALSE])$predictions
  preds_prob[te, colnames(pp)] <- pp
  cat(sprintf("  fold %d: train %d, test %d pixels (%d sites)\n",
              f, length(tr), length(te),
              dplyr::n_distinct(site_key[te])))
}

# Pixel-level RF prediction (argmax probability)
preds_pixel <- factor(levels(y)[apply(preds_prob, 1, which.max)],
                      levels = levels(y))

# --- 5. Aggregate to site level: mean prob across pixels -> argmax -------
site_pred <- pixels |>
  dplyr::mutate(.row = dplyr::row_number()) |>
  dplyr::select(site_number, Year, final_label, .row) |>
  dplyr::group_by(site_number, Year, final_label) |>
  dplyr::summarise(.rows = list(.row), .groups = "drop") |>
  dplyr::mutate(
    mean_probs = purrr::map(.rows, \(rr)
      colMeans(preds_prob[rr, , drop = FALSE], na.rm = TRUE)
    ),
    site_pred = purrr::map_chr(mean_probs, \(p)
      levels(y)[which.max(replace(p, is.na(p), -Inf))]
    )
  ) |>
  dplyr::mutate(correct = site_pred == final_label)

site_pred <- site_pred |>
  dplyr::mutate(site_pred = factor(site_pred, levels = levels(y)),
                truth     = factor(final_label, levels = levels(y)))

# --- 6. Site-level recall vs site-averaged baseline ---------------------
pixel_recall <- site_pred |>
  dplyr::group_by(truth) |>
  dplyr::summarise(n_test = dplyr::n(),
                   pixel_recall = mean(correct),
                   .groups = "drop") |>
  dplyr::rename(final_label = truth) |>
  dplyr::mutate(final_label = as.character(final_label))

# baseline recall from 34 (balanced model on site-averaged spectra).
baseline_recall <- baseline$recall  # has unweighted/balanced cols
compare <- pixel_recall |>
  dplyr::left_join(
    baseline_recall |>
      dplyr::select(final_label, balanced_recall, n) |>
      dplyr::rename(site_avg_recall = balanced_recall),
    by = "final_label"
  ) |>
  dplyr::arrange(dplyr::desc(n))

cat("\n=== Site-level recall: pixel-trained vs site-averaged baseline ===\n")
print(as.data.frame(compare))

overall_pixel <- mean(site_pred$correct)
overall_site  <- baseline$overall_acc
cat(sprintf(
  "\nOverall accuracy:  pixel-trained = %.1f%%   site-averaged = %.1f%%\n",
  100 * overall_pixel, 100 * overall_site
))

# --- 7. Plots -----------------------------------------------------------
dir.create("output", showWarnings = FALSE)
plot_df <- compare |>
  tidyr::pivot_longer(c(pixel_recall, site_avg_recall),
                      names_to = "model", values_to = "recall")
p_rec <- ggplot(plot_df,
                aes(x = forcats::fct_reorder(final_label, n),
                    y = recall, fill = model)) +
  geom_col(position = "dodge") +
  coord_flip() +
  scale_fill_manual(values = c(pixel_recall = "steelblue",
                               site_avg_recall = "tomato")) +
  labs(x = NULL, y = "Recall",
       title = "Shrub class recall: pixel-trained vs site-averaged",
       subtitle = sprintf("Pixel-trained overall = %.1f%%   site-averaged overall = %.1f%%",
                          100 * overall_pixel, 100 * overall_site)) +
  theme_minimal(base_size = 11)
ggsave("output/shrub_pixel_vs_site_recall.png", p_rec,
       width = 8, height = 5, dpi = 150)

conf <- tibble::tibble(truth = site_pred$truth,
                      pred  = site_pred$site_pred) |>
  dplyr::count(truth, pred, name = "n") |>
  tidyr::complete(truth, pred, fill = list(n = 0L))
label_order <- compare$final_label
p_conf <- ggplot(conf,
                 aes(x = factor(pred, levels = label_order),
                     y = factor(truth, levels = rev(label_order)),
                     fill = n)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = ifelse(n > 0, n, "")), size = 3) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(x = "Predicted (site-level argmax)", y = "Truth",
       title = "Pixel-trained shrub classifier — site-level confusion",
       subtitle = sprintf("Overall site accuracy = %.1f%%",
                          100 * overall_pixel)) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  coord_fixed()
ggsave("output/shrub_pixel_confusion.png", p_conf,
       width = 9, height = 9, dpi = 150)

# --- 8. Persist ---------------------------------------------------------
saveRDS(list(
  pixel_predictions = tibble::tibble(
    site_number = pixels$site_number,
    Year        = pixels$Year,
    final_label = pixels$final_label,
    preds_pixel = preds_pixel
  ),
  site_predictions  = dplyr::select(site_pred, site_number, Year, truth,
                                    site_pred, correct),
  recall_compare    = compare,
  overall           = list(pixel_trained = overall_pixel,
                           site_averaged = overall_site)
), "data/derived/shrub_pixel_training.rds")
cat("\nWrote data/derived/shrub_pixel_training.rds\n")
cat("Wrote output/shrub_pixel_vs_site_recall.png\n")
cat("Wrote output/shrub_pixel_confusion.png\n")
