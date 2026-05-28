# 15_stratified_year_correction.R — prototype an NDVI-stratified
# 2018->2025 radiometric/phenology correction.
#
# The non-vegetated correction in script 13 (single global per-band delta
# fit on rock/bare points only) collapsed ~75% of the residual at
# non-vegetated control points but, per the pair diagnostic in script 14,
# left a Cohen's d up to -8 in the BLUE region for vegetated communities
# clustering apart by year. That suggests vegetation-specific phenology /
# chlorophyll absorption that the rock baseline can't carry.
#
# Approach: bin the year_effect points by NDVI (so they cover the entire
# bare→dense-vegetation gradient instead of just bare), and for each bin
# compute a per-band mean delta. Apply by per-site NDVI lookup.
#
# Inputs:
#   data/derived/year_effect_spectra.parquet
#   data/raw/ESS-DIVE-Spectra/wavelengths_2025.csv
# Outputs:
#   data/derived/year_effect_correction_by_ndvi.csv
#   data/small_reference/year_effect_correction_2018_to_2025_by_ndvi.csv
#   docs/figures/year_effect_stratified_correction.pdf
#   docs/figures/year_effect_stratified_validation.pdf

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(nanoparquet)
})
dir.create("docs/figures", showWarnings = FALSE, recursive = TRUE)

spec <- nanoparquet::read_parquet("data/derived/year_effect_spectra.parquet")
wls  <- readr::read_csv("data/raw/ESS-DIVE-Spectra/wavelengths_2025.csv",
                         show_col_types = FALSE) |>
  dplyr::transmute(band_number, wavelength_nm = center_wavelength * 1000)
band_cols <- sprintf("rfl_band_%d", seq_len(426))

# L2-normalize per row.
mat <- as.matrix(spec[, band_cols])
norms <- sqrt(rowSums(mat^2, na.rm = TRUE))
norms[norms == 0 | !is.finite(norms)] <- NA_real_
spec[, band_cols] <- as.data.frame(mat / norms)

# Pair by point_id.
s18 <- spec |> dplyr::filter(year == 2018L)
s25 <- spec |> dplyr::filter(year == 2025L)
paired <- dplyr::inner_join(
  s18 |> dplyr::select(point_id, point_type, doy_band, snow_free_doy,
                       ndvi_18 = ndvi, dplyr::all_of(band_cols)),
  s25 |> dplyr::select(point_id,
                       ndvi_25 = ndvi, dplyr::all_of(band_cols)),
  by = "point_id",
  suffix = c("_18", "_25")
)
cat(sprintf("Paired points: %d (NDVI[2018] median %.2f, NDVI[2025] median %.2f)\n",
            nrow(paired),
            stats::median(paired$ndvi_18, na.rm = TRUE),
            stats::median(paired$ndvi_25, na.rm = TRUE)))

# Per-pair mean NDVI determines the bin (so a point that browned/greened
# between years lands somewhere in between).
paired$ndvi_avg <- (paired$ndvi_18 + paired$ndvi_25) / 2
ndvi_breaks <- c(-Inf, 0.20, 0.40, 0.60, 0.80, Inf)
ndvi_labels <- c("ndvi_lt_0.20", "ndvi_0.20_0.40",
                 "ndvi_0.40_0.60", "ndvi_0.60_0.80", "ndvi_ge_0.80")
paired$ndvi_bin <- cut(paired$ndvi_avg, breaks = ndvi_breaks,
                       labels = ndvi_labels, right = FALSE,
                       include.lowest = TRUE)
cat("\nPaired-point counts per NDVI bin:\n")
print(paired |> dplyr::count(ndvi_bin) |> as.data.frame())

mat_18 <- as.matrix(paired[, paste0(band_cols, "_18")])
mat_25 <- as.matrix(paired[, paste0(band_cols, "_25")])
diff_mat <- mat_25 - mat_18

# --- Per-bin per-band delta ---------------------------------------------
per_bin_delta <- function(idx) {
  if (length(idx) == 0) {
    return(matrix(NA_real_, nrow = 1, ncol = 426))
  }
  d <- diff_mat[idx, , drop = FALSE]
  matrix(colMeans(d, na.rm = TRUE), nrow = 1)
}
bin_idx <- split(seq_len(nrow(paired)), paired$ndvi_bin)

deltas <- purrr::map_dfr(ndvi_labels, function(b) {
  d <- per_bin_delta(bin_idx[[b]])
  tibble::tibble(
    ndvi_bin      = b,
    band_number   = seq_len(426),
    wavelength_nm = wls$wavelength_nm,
    delta         = as.numeric(d),
    n_pairs       = length(bin_idx[[b]] %||% integer(0))
  )
})
readr::write_csv(deltas, "data/derived/year_effect_correction_by_ndvi.csv")
cat(sprintf("\nWrote data/derived/year_effect_correction_by_ndvi.csv (%d rows)\n",
            nrow(deltas)))

# Save the small_reference copy (canonical, committed).
dir.create("data/small_reference", showWarnings = FALSE, recursive = TRUE)
readr::write_csv(
  deltas,
  "data/small_reference/year_effect_correction_2018_to_2025_by_ndvi.csv"
)

# --- 50/50 held-out validation, stratified vs global ---------------------
set.seed(42)
fold <- sample(c("fit", "val"), nrow(paired),
                replace = TRUE, prob = c(0.5, 0.5))

# Strategy 1: global non-vegetated-only correction (current production).
nv_fit <- which(fold == "fit" & paired$point_type == "non_vegetated")
delta_nv <- colMeans(diff_mat[nv_fit, , drop = FALSE], na.rm = TRUE)

# Strategy 2: NDVI-stratified delta from the fit set (vegetated +
# non-vegetated pooled per bin).
delta_strat_fit <- vapply(ndvi_labels, function(b) {
  idx <- which(fold == "fit" & paired$ndvi_bin == b)
  if (length(idx) == 0) return(rep(NA_real_, 426))
  colMeans(diff_mat[idx, , drop = FALSE], na.rm = TRUE)
}, numeric(426))
# delta_strat_fit[band, bin]

# Validation set
val_idx     <- which(fold == "val")
val_diff    <- diff_mat[val_idx, , drop = FALSE]
val_ndvi_bin <- paired$ndvi_bin[val_idx]
val_type    <- paired$point_type[val_idx]

corrected_global <- sweep(val_diff, 2, delta_nv, FUN = "-")
corrected_strat  <- val_diff
for (k in seq_along(val_idx)) {
  b <- as.character(val_ndvi_bin[k])
  if (is.na(b) || !(b %in% colnames(delta_strat_fit))) next
  d <- delta_strat_fit[, b]
  if (all(is.na(d))) next
  corrected_strat[k, ] <- val_diff[k, ] - d
}

scoring <- function(label, point_type_filter = NULL) {
  if (is.null(point_type_filter)) {
    sel <- seq_len(nrow(val_diff))
  } else {
    sel <- which(val_type == point_type_filter)
  }
  if (length(sel) == 0) return(NULL)
  pooled_sd <- apply(rbind(mat_18[val_idx[sel], ],
                            mat_25[val_idx[sel], ]),
                     2, sd, na.rm = TRUE)
  raw_mean       <- colMeans(val_diff[sel, ],      na.rm = TRUE)
  global_mean    <- colMeans(corrected_global[sel, ], na.rm = TRUE)
  strat_mean     <- colMeans(corrected_strat[sel, ],  na.rm = TRUE)
  c(label              = label,
    n                  = length(sel),
    mean_abs_raw       = sprintf("%.5f", mean(abs(raw_mean))),
    mean_abs_global    = sprintf("%.5f", mean(abs(global_mean))),
    mean_abs_strat     = sprintf("%.5f", mean(abs(strat_mean))),
    mean_d_raw         = sprintf("%.3f", mean(abs(raw_mean    / pmax(pooled_sd, 1e-9)))),
    mean_d_global      = sprintf("%.3f", mean(abs(global_mean / pmax(pooled_sd, 1e-9)))),
    mean_d_strat       = sprintf("%.3f", mean(abs(strat_mean  / pmax(pooled_sd, 1e-9))))
  )
}

cat("\n=== Held-out residual: raw vs global (non-veg only) vs NDVI-stratified ===\n")
tbl <- dplyr::bind_rows(
  as.list(scoring("ALL")),
  as.list(scoring("non_vegetated", "non_vegetated")),
  as.list(scoring("vegetated",     "vegetated"))
)
print(as.data.frame(tbl), row.names = FALSE)

# --- Plot: per-bin delta spectrum ---------------------------------------
water_bands <- tibble::tibble(
  xmin = c(1340, 1800, 2400),
  xmax = c(1450, 1950, 2510)
)

p_strat <- ggplot(deltas |> dplyr::filter(!is.na(delta)),
                   aes(x = wavelength_nm, y = delta, colour = ndvi_bin)) +
  geom_rect(data = water_bands, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey90", alpha = 0.6) +
  geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
  geom_line(linewidth = 0.6) +
  scale_colour_brewer(palette = "RdYlGn", name = "NDVI bin",
                      direction = 1) +
  labs(x = "Wavelength (nm)",
       y = "Mean rfl difference (2025 − 2018), per NDVI bin",
       title = "NDVI-stratified 2018->2025 correction",
       subtitle = "Drift signature differs sharply with NDVI — single global delta misses the vegetated story") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        plot.title.position = "plot")
ggsave("docs/figures/year_effect_stratified_correction.pdf", p_strat,
       width = 11, height = 6, device = cairo_pdf)

# --- Plot: residual reduction comparison --------------------------------
val_long <- tibble::tibble(
  wavelength_nm = wls$wavelength_nm,
  raw      = colMeans(val_diff,         na.rm = TRUE),
  global   = colMeans(corrected_global, na.rm = TRUE),
  stratified = colMeans(corrected_strat, na.rm = TRUE)
) |>
  tidyr::pivot_longer(c(raw, global, stratified),
                      names_to = "method", values_to = "mean_resid") |>
  dplyr::mutate(method = factor(method,
                                 levels = c("raw", "global", "stratified")))
p_val <- ggplot(val_long,
                 aes(x = wavelength_nm, y = mean_resid, colour = method)) +
  geom_rect(data = water_bands, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey90", alpha = 0.6) +
  geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = c(raw = "#d6604d",
                                  global = "#fdb863",
                                  stratified = "#4393c3")) +
  labs(x = "Wavelength (nm)",
       y = "Mean per-band residual (2025 − 2018), held-out half",
       title = "NDVI-stratified vs global non-vegetated correction",
       subtitle = "Flatter line at 0 = better generalization") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        plot.title.position = "plot")
ggsave("docs/figures/year_effect_stratified_validation.pdf", p_val,
       width = 11, height = 5, device = cairo_pdf)

cat("\nWrote docs/figures/year_effect_stratified_correction.pdf\n")
cat("Wrote docs/figures/year_effect_stratified_validation.pdf\n")
