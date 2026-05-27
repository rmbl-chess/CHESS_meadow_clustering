# 13_year_effect_analysis.R — paired 2018-vs-2025 comparison at the
# CRBU year-effect points (script 12). Goal: separate two kinds of
# year-to-year spectral drift:
#
#   - phenology shift     (vegetated points only) — same site sampled
#                         in two different growing seasons may differ
#                         genuinely; expected in NDVI, red-edge slope,
#                         and the NIR plateau.
#   - radiometric drift   (non-vegetated points) — rock / bare / sparse
#                         veg should be spectrally stable year-to-year;
#                         any systematic shift in per-band reflectance
#                         is instrument / processing.
#
# Method:
#   - Pair extracted spectra by point_id. Keep only points present in
#     BOTH years (inner join).
#   - Per band: diff = rfl(2025) - rfl(2018); paired t-test; effect
#     size = mean_diff / pooled_sd_per_band.
#   - Per point: per-pair NDVI, NDWI, brightness (mean rfl across the
#     non-water-band wavelengths).
#   - Stratify by point_type and (optionally) doy_band.
#
# Inputs:
#   data/derived/year_effect_spectra.parquet
#   data/raw/ESS-DIVE-Spectra/wavelengths_2025.csv  (band -> wavelength)
# Outputs:
#   data/derived/year_effect_summary.csv             per-band stats
#   data/derived/year_effect_per_point.csv           per-pair indices
#   docs/figures/year_effect_mean_diff_spectrum.pdf
#   docs/figures/year_effect_ndvi_scatter.pdf
#   docs/figures/year_effect_brightness_scatter.pdf
#   docs/figures/year_effect_per_band_effect.pdf

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  if (!requireNamespace("nanoparquet", quietly = TRUE)) {
    renv::install("nanoparquet")
  }
  library(nanoparquet)
})

dir.create("docs/figures", showWarnings = FALSE, recursive = TRUE)

# --- Load ---------------------------------------------------------------
spec <- nanoparquet::read_parquet("data/derived/year_effect_spectra.parquet")
cat(sprintf("Spectra: %d rows x %d cols\n", nrow(spec), ncol(spec)))

wls <- readr::read_csv("data/raw/ESS-DIVE-Spectra/wavelengths_2025.csv",
                       show_col_types = FALSE) |>
  dplyr::transmute(band_number, wavelength_nm = center_wavelength * 1000)
stopifnot(nrow(wls) == 426)

band_cols <- sprintf("rfl_band_%d", seq_len(426))

# --- L2-normalize each row over all 426 bands ----------------------------
# Matches the classifier preprocessing recipe: per-spectrum L2 unit
# vector so the comparison is on spectral SHAPE, not on absolute
# brightness. Raw values are kept on the side (`spec_raw`) for the
# brightness scatter — L2 norm by definition removes overall brightness,
# which is exactly the radiometric-calibration signal we want to see for
# non-vegetated pixels.
spec_raw <- spec
spec_l2  <- spec
mat <- as.matrix(spec_l2[, band_cols])
norms <- sqrt(rowSums(mat^2, na.rm = TRUE))
norms[norms == 0 | !is.finite(norms)] <- NA_real_
mat <- mat / norms
spec_l2[, band_cols] <- as.data.frame(mat)
dropped <- sum(is.na(norms))
cat(sprintf("L2-normalized %d rows (%d dropped: zero / invalid norm)\n",
            sum(!is.na(norms)), dropped))

# --- Pair by point_id (inner join 2018 + 2025) on L2-normed spectra -----
pair_join <- function(df) {
  s18 <- df |> dplyr::filter(year == 2018L)
  s25 <- df |> dplyr::filter(year == 2025L)
  dplyr::inner_join(
    s18 |> dplyr::select(point_id, point_type, doy_band, snow_free_doy,
                         dplyr::all_of(band_cols)),
    s25 |> dplyr::select(point_id, dplyr::all_of(band_cols)),
    by = "point_id",
    suffix = c("_18", "_25")
  )
}
paired      <- pair_join(spec_l2)    # L2-normed pairs (primary)
paired_raw  <- pair_join(spec_raw)   # raw pairs (for brightness only)
cat(sprintf("Paired points (in both years): %d\n", nrow(paired)))
cat("By point_type x doy_band:\n")
print(paired |> dplyr::count(point_type, doy_band) |>
        tidyr::pivot_wider(names_from = doy_band, values_from = n,
                           values_fill = 0L) |> as.data.frame())

# --- Difference matrix (2025 - 2018) per band ---------------------------
mat_18 <- as.matrix(paired[, paste0(band_cols, "_18")])
mat_25 <- as.matrix(paired[, paste0(band_cols, "_25")])
diff_mat <- mat_25 - mat_18

# Per (point_type, band) stats
per_stratum_band <- function(strat) {
  idx <- which(paired$point_type == strat)
  d   <- diff_mat[idx, , drop = FALSE]
  v18 <- mat_18  [idx, , drop = FALSE]
  v25 <- mat_25  [idx, , drop = FALSE]
  tibble::tibble(
    point_type     = strat,
    band_number    = seq_len(426),
    wavelength_nm  = wls$wavelength_nm,
    mean_18        = colMeans(v18, na.rm = TRUE),
    mean_25        = colMeans(v25, na.rm = TRUE),
    mean_diff      = colMeans(d, na.rm = TRUE),
    sd_diff        = apply(d, 2, sd, na.rm = TRUE),
    n_pairs        = colSums(!is.na(d)),
    pooled_sd      = sqrt(0.5 * (apply(v18, 2, var, na.rm = TRUE) +
                                  apply(v25, 2, var, na.rm = TRUE))),
    t_stat         = mean_diff / (sd_diff / sqrt(n_pairs)),
    cohens_d       = mean_diff / pmax(pooled_sd, 1e-9),
    ci95_lo        = mean_diff - 1.96 * sd_diff / sqrt(n_pairs),
    ci95_hi        = mean_diff + 1.96 * sd_diff / sqrt(n_pairs)
  )
}
summary_band <- purrr::map_dfr(c("vegetated", "non_vegetated"),
                                per_stratum_band)
readr::write_csv(summary_band, "data/derived/year_effect_summary.csv")
cat(sprintf("\nWrote data/derived/year_effect_summary.csv (%d rows)\n",
            nrow(summary_band)))

# --- Per-point indices for paired comparison ----------------------------
# Pick the nearest band to each target wavelength (in nm).
b_at <- function(target_nm) which.min(abs(wls$wavelength_nm - target_nm))

idx_red <- b_at(660); idx_nir <- b_at(860)
idx_swir <- b_at(1240)
idx_blue <- b_at(490); idx_green <- b_at(560)

ndvi   <- function(spec) (spec[, idx_nir] - spec[, idx_red]) /
                          (spec[, idx_nir] + spec[, idx_red])
ndwi   <- function(spec) (spec[, idx_nir] - spec[, idx_swir]) /
                          (spec[, idx_nir] + spec[, idx_swir])
# Brightness on visible-NIR (drop water bands implicitly): mean reflectance
# across bands < 1340 nm.
vis_idx <- which(wls$wavelength_nm < 1340)
brightness <- function(spec) rowMeans(spec[, vis_idx], na.rm = TRUE)

# Raw matrices for brightness only (paired_raw has the same point_id
# order as `paired` — match by point_id to be safe).
mat_raw_18 <- as.matrix(paired_raw[match(paired$point_id,
                                          paired_raw$point_id),
                                    paste0(band_cols, "_18")])
mat_raw_25 <- as.matrix(paired_raw[match(paired$point_id,
                                          paired_raw$point_id),
                                    paste0(band_cols, "_25")])

per_point <- paired |>
  dplyr::select(point_id, point_type, doy_band, snow_free_doy) |>
  dplyr::mutate(
    # NDVI / NDWI are band-ratio invariants — same on L2-normed or raw.
    ndvi_18       = ndvi(mat_18),
    ndvi_25       = ndvi(mat_25),
    ndvi_diff     = ndvi_25 - ndvi_18,
    ndwi_18       = ndwi(mat_18),
    ndwi_25       = ndwi(mat_25),
    ndwi_diff     = ndwi_25 - ndwi_18,
    # Brightness on RAW reflectance — L2 norm would zero this signal.
    brightness_18 = brightness(mat_raw_18),
    brightness_25 = brightness(mat_raw_25),
    brightness_diff = brightness_25 - brightness_18
  )
readr::write_csv(per_point, "data/derived/year_effect_per_point.csv")
cat("Wrote data/derived/year_effect_per_point.csv\n")

cat("\nPer-(point_type) summary of per-pair indices:\n")
print(per_point |>
        dplyr::group_by(point_type) |>
        dplyr::summarise(
          n               = dplyr::n(),
          mean_ndvi_diff  = mean(ndvi_diff,       na.rm = TRUE),
          sd_ndvi_diff    = sd(ndvi_diff,         na.rm = TRUE),
          mean_brightness_diff = mean(brightness_diff, na.rm = TRUE),
          sd_brightness_diff   = sd(brightness_diff,   na.rm = TRUE),
          .groups = "drop"
        ) |> as.data.frame())

# --- Plot 1: mean difference spectrum with 95% CI ribbon ----------------
water_bands <- tibble::tibble(
  xmin = c(1340, 1800, 2400),
  xmax = c(1450, 1950, 2510)
)

p1 <- ggplot(summary_band,
             aes(x = wavelength_nm, y = mean_diff, colour = point_type,
                 fill = point_type)) +
  geom_rect(data = water_bands, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey90", alpha = 0.6) +
  geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
  geom_ribbon(aes(ymin = ci95_lo, ymax = ci95_hi), alpha = 0.2,
              colour = NA) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = c(vegetated = "#1b7837",
                                 non_vegetated = "#8c510a")) +
  scale_fill_manual(values   = c(vegetated = "#1b7837",
                                 non_vegetated = "#8c510a")) +
  labs(x = "Wavelength (nm)",
       y = "L2-normalized reflectance difference (2025 − 2018)",
       title = "Paired year-effect on AOP spectral shape (CRBU)",
       subtitle = sprintf("Mean ± 95%% CI per band on L2-normalized spectra; %d vegetated + %d non-vegetated paired points",
                          sum(paired$point_type == "vegetated"),
                          sum(paired$point_type == "non_vegetated"))) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        plot.title.position = "plot")
ggsave("docs/figures/year_effect_mean_diff_spectrum.pdf", p1,
       width = 11, height = 6, device = cairo_pdf)

# --- Plot 2: per-band Cohen's d (standardized effect) -------------------
p2 <- ggplot(summary_band,
             aes(x = wavelength_nm, y = cohens_d, colour = point_type)) +
  geom_rect(data = water_bands, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey90", alpha = 0.6) +
  geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
  geom_hline(yintercept = c(-0.2, 0.2),
             colour = "grey70", linetype = "dotted") +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c(vegetated = "#1b7837",
                                 non_vegetated = "#8c510a")) +
  labs(x = "Wavelength (nm)",
       y = "Cohen's d (mean_diff / pooled SD)",
       title = "Year-effect magnitude by wavelength (L2-normalized; standardized)",
       subtitle = "Dotted lines at |d| = 0.2 mark the conventional small-effect band") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        plot.title.position = "plot")
ggsave("docs/figures/year_effect_per_band_effect.pdf", p2,
       width = 11, height = 5, device = cairo_pdf)

# --- Plot 3: NDVI 2018 vs 2025 scatter ----------------------------------
p3 <- ggplot(per_point, aes(x = ndvi_18, y = ndvi_25, colour = doy_band)) +
  geom_abline(slope = 1, intercept = 0,
              colour = "grey50", linetype = "dashed") +
  geom_point(alpha = 0.4, size = 0.7) +
  facet_wrap(~ point_type) +
  scale_colour_brewer(palette = "Set1") +
  labs(x = "NDVI (2018)", y = "NDVI (2025)",
       title = "NDVI 2018 vs 2025 at common CRBU points",
       subtitle = "Above-diagonal = greener in 2025; below-diagonal = greener in 2018") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        plot.title.position = "plot",
        strip.text = element_text(face = "bold")) +
  coord_fixed()
ggsave("docs/figures/year_effect_ndvi_scatter.pdf", p3,
       width = 10, height = 6, device = cairo_pdf)

# --- Plot 4: brightness 2018 vs 2025 scatter ----------------------------
p4 <- ggplot(per_point,
             aes(x = brightness_18, y = brightness_25, colour = doy_band)) +
  geom_abline(slope = 1, intercept = 0,
              colour = "grey50", linetype = "dashed") +
  geom_point(alpha = 0.4, size = 0.7) +
  facet_wrap(~ point_type, scales = "free") +
  scale_colour_brewer(palette = "Set1") +
  labs(x = "Visible-NIR brightness (2018) — raw reflectance",
       y = "Visible-NIR brightness (2025) — raw reflectance",
       title = "Per-point visible-NIR brightness, paired CRBU points (raw)",
       subtitle = "Computed on raw reflectance (NOT L2-normed) so radiometric drift survives.\nNon-vegetated drift off the 1:1 line = strongest calibration signal.") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        plot.title.position = "plot",
        strip.text = element_text(face = "bold"))
ggsave("docs/figures/year_effect_brightness_scatter.pdf", p4,
       width = 10, height = 6, device = cairo_pdf)

cat("\nWrote 4 PDFs in docs/figures/year_effect_*\n")

# --- Tabular bottom-line ------------------------------------------------
cat("\n=== Bands with largest standardized effect (|Cohen's d| top 8 per type) ===\n")
print(summary_band |>
        dplyr::group_by(point_type) |>
        dplyr::slice_max(abs(cohens_d), n = 8) |>
        dplyr::ungroup() |>
        dplyr::select(point_type, wavelength_nm, mean_diff,
                       cohens_d, t_stat, n_pairs) |>
        dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 3))) |>
        dplyr::arrange(point_type, dplyr::desc(abs(cohens_d))) |>
        as.data.frame())
