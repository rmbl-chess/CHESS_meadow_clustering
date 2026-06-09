# 16_year_mean_spectra.R — mean brightness-normalized (L2) spectrum per year.
#
# Companion to the year-effect DELTA figures (script 15): instead of the
# 2025 − 2018 difference, plot the mean L2-normalized reflectance spectrum
# for each campaign year, faceted by vegetated / non-vegetated point type,
# so the actual spectral shapes (and the year offset within each type) are
# visible. Same L2-normalization (per row) and water-band shading as the
# correction fit.
#
# Inputs:
#   data/derived/year_effect_spectra.parquet
#   data/raw/ESS-DIVE-Spectra/wavelengths_2025.csv
# Outputs:
#   docs/figures/year_mean_spectra_2018_2025.pdf

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

# Drop the latest snowmelt stratum (doy_band == "late", snow-free DOY ~155+):
# late-melting bare sites can still hold snow at AOP imaging time, which would
# contaminate the non-vegetated spectra with snow (bright, flat-ish).
drop_doy <- "late"
n_before <- nrow(spec)
spec <- dplyr::filter(spec, !doy_band %in% drop_doy)
cat(sprintf("Excluded doy_band {%s}: %d -> %d points\n",
            paste(drop_doy, collapse = ", "), n_before, nrow(spec)))

# Keep only MATCHED PAIRS: point_ids with spectra present in BOTH years, so
# the 2018 vs 2025 means compare the same physical locations.
paired_ids <- spec |> dplyr::distinct(point_id, year) |>
  dplyr::count(point_id, name = "n_yr") |>
  dplyr::filter(n_yr == 2L) |> dplyr::pull(point_id)
n_pre <- nrow(spec)
spec <- dplyr::filter(spec, point_id %in% paired_ids)
cat(sprintf("Matched pairs (both years): %d point_ids -> %d rows (from %d)\n",
            length(paired_ids), nrow(spec), n_pre))

# L2-normalize each spectrum (unit brightness), same as 04/15.
mat   <- as.matrix(spec[, band_cols])
norms <- sqrt(rowSums(mat^2, na.rm = TRUE))
norms[norms == 0 | !is.finite(norms)] <- NA_real_
matn  <- mat / norms

# Mean + SD per band, split by year AND vegetated / non-vegetated point type.
combos <- tidyr::expand_grid(year = c(2018L, 2025L),
                             ptype = c("vegetated", "non_vegetated"))
per_year <- purrr::pmap_dfr(combos, function(year, ptype) {
  m <- matn[spec$year == year & spec$point_type == ptype, , drop = FALSE]
  tibble::tibble(year = year, point_type = ptype, band_number = seq_len(426),
                 mean_rfl = colMeans(m, na.rm = TRUE),
                 sd_rfl   = apply(m, 2, stats::sd, na.rm = TRUE))
}) |>
  dplyr::left_join(wls, by = "band_number") |>
  dplyr::mutate(year = factor(year),
                point_type = factor(point_type,
                  levels = c("vegetated", "non_vegetated"),
                  labels = c("Vegetated", "Non-vegetated")))

# facet strip labels carry per-facet sample sizes
ns <- table(spec$point_type, spec$year)
fac_labels <- c(
  "Vegetated"     = sprintf("Vegetated (2018 n=%d, 2025 n=%d)",
                            ns["vegetated", "2018"], ns["vegetated", "2025"]),
  "Non-vegetated" = sprintf("Non-vegetated (2018 n=%d, 2025 n=%d)",
                            ns["non_vegetated", "2018"], ns["non_vegetated", "2025"]))

# Water-absorption bands (noisy after atmospheric correction) — shade them.
water_bands <- tibble::tibble(xmin = c(1340, 1800, 2400),
                              xmax = c(1450, 1950, 2510))

yr_cols <- c("2018" = "#d6604d", "2025" = "#4393c3")

p <- ggplot(per_year, aes(x = wavelength_nm, colour = year, fill = year)) +
  geom_rect(data = water_bands, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey90", alpha = 0.6) +
  geom_ribbon(aes(ymin = mean_rfl - sd_rfl, ymax = mean_rfl + sd_rfl),
              colour = NA, alpha = 0.18) +
  geom_line(aes(y = mean_rfl), linewidth = 0.7) +
  facet_wrap(~ point_type, ncol = 1,
             labeller = ggplot2::as_labeller(fac_labels)) +
  scale_colour_manual(values = yr_cols, name = "Year") +
  scale_fill_manual(values = yr_cols, guide = "none") +
  labs(x = "Wavelength (nm)",
       y = "Mean L2-normalized reflectance",
       title = "Mean brightness-normalized spectra, 2018 vs 2025",
       subtitle = paste("Matched pairs (same points in both years), by point",
                        "type; ribbon = ±1 SD. Late-snowmelt stratum excluded",
                        "(snow risk). Grey = masked water bands.")) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", plot.title.position = "plot",
        strip.text = element_text(face = "bold"))

ggsave("docs/figures/year_mean_spectra_2018_2025.pdf", p,
       width = 11, height = 8.5, device = cairo_pdf)
cat("Wrote docs/figures/year_mean_spectra_2018_2025.pdf\n")
