# 14_year_separated_pairs.R — diagnose the residual year-effect by
# spectral feature.
#
# Strategy: pick pairs of clusters that we suspect are the SAME ecological
# community split by year (same dominant indicator species but one cluster
# is 2018-only and the other is 2025-only). For each pair, compare mean
# spectra band-by-band on the corrected, L2-normalized, mean-per-site
# spectra that the cluster sees. Where do they actually disagree, and
# what does the disagreement look like — phenology (red-edge, NIR), or
# something more global?
#
# Pairs (chosen from code/joint/13_year_effect_analysis.R follow-up):
#   Veratrum tall-forb wet meadow      S02 (2018-only) vs S21 (2025-only)
#   Ligusticum mesic tall-forb         S06 (2018-only) vs S19 (2025-only)
#   Vaccinium dwarf-shrub meadow       S05 (2018-only) vs S18 (2025-only)
#   Carex aquatilis wet sedge          S11 (2018-only) vs S22 (2025-only)
#
# Inputs:
#   data/derived/veg_spectra.rds              corrected per-site spectra
#   data/derived/final_clusters_B.rds         final cluster assignments
#   data/raw/ESS-DIVE-Spectra/wavelengths_2025.csv
# Outputs:
#   data/derived/year_separated_pair_summary.csv  per-pair per-band stats
#   docs/figures/year_separated_pair_diffs.pdf    4-panel diff spectra
#   docs/figures/year_separated_pair_means.pdf    overlaid mean spectra

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})
dir.create("docs/figures", showWarnings = FALSE, recursive = TRUE)

vs <- readRDS("data/derived/veg_spectra.rds")$joined
fc <- readRDS("data/derived/final_clusters_B.rds")$assignments
wls <- readr::read_csv("data/raw/ESS-DIVE-Spectra/wavelengths_2025.csv",
                       show_col_types = FALSE) |>
  dplyr::transmute(band_number,
                   wavelength_nm = center_wavelength * 1000)

band_cols <- sprintf("rfl_band_%d", seq_len(426))

# Pairs to investigate: (community, label_2018, label_2025)
pairs <- tibble::tribble(
  ~community,                       ~label_2018, ~label_2025,
  "Veratrum tall-forb wet meadow",  "S02",       "S21",
  "Ligusticum mesic tall-forb",     "S06",       "S19",
  "Vaccinium dwarf-shrub",          "S05",       "S18",
  "Carex aquatilis wet sedge",      "S11",       "S22"
)

# Join cluster labels onto veg_spectra.
spec <- vs |>
  dplyr::select(site_number, Year, dplyr::all_of(band_cols)) |>
  dplyr::inner_join(fc |> dplyr::select(site_number, Year, final_label),
                     by = c("site_number", "Year"))

water_bands <- tibble::tibble(
  xmin = c(1340, 1800, 2400),
  xmax = c(1450, 1950, 2510)
)

# Per-pair stats (band-by-band).
per_pair <- purrr::map_dfr(seq_len(nrow(pairs)), function(i) {
  comm <- pairs$community[i]
  lab18 <- pairs$label_2018[i]
  lab25 <- pairs$label_2025[i]
  m18 <- spec |> dplyr::filter(final_label == lab18) |>
    dplyr::select(dplyr::all_of(band_cols)) |> as.matrix()
  m25 <- spec |> dplyr::filter(final_label == lab25) |>
    dplyr::select(dplyr::all_of(band_cols)) |> as.matrix()
  if (nrow(m18) == 0 || nrow(m25) == 0) return(NULL)
  mean_18 <- colMeans(m18, na.rm = TRUE)
  mean_25 <- colMeans(m25, na.rm = TRUE)
  pooled_sd <- sqrt(0.5 * (apply(m18, 2, var, na.rm = TRUE) +
                            apply(m25, 2, var, na.rm = TRUE)))
  tibble::tibble(
    community     = comm,
    pair_label    = sprintf("%s: %s (n=%d) vs %s (n=%d)",
                             comm, lab18, nrow(m18), lab25, nrow(m25)),
    band_number   = seq_len(426),
    wavelength_nm = wls$wavelength_nm,
    n_2018        = nrow(m18),
    n_2025        = nrow(m25),
    mean_18       = mean_18,
    mean_25       = mean_25,
    mean_diff     = mean_18 - mean_25,
    pooled_sd     = pooled_sd,
    cohens_d      = mean_diff / pmax(pooled_sd, 1e-9)
  )
})

readr::write_csv(per_pair, "data/derived/year_separated_pair_summary.csv")
cat(sprintf("Wrote data/derived/year_separated_pair_summary.csv (%d rows)\n",
            nrow(per_pair)))

# Top differing bands per pair.
cat("\n=== Top 6 differing bands per pair (by |Cohen's d|) ===\n")
print(per_pair |>
        dplyr::group_by(pair_label) |>
        dplyr::slice_max(abs(cohens_d), n = 6) |>
        dplyr::ungroup() |>
        dplyr::select(pair_label, wavelength_nm, mean_18, mean_25,
                       mean_diff, cohens_d) |>
        dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 3))) |>
        dplyr::arrange(pair_label, dplyr::desc(abs(cohens_d))) |>
        as.data.frame(), row.names = FALSE)

# Plot 1: difference spectrum per pair.
p_diff <- ggplot(per_pair,
                  aes(x = wavelength_nm, y = mean_diff)) +
  geom_rect(data = water_bands, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey90", alpha = 0.6) +
  geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
  geom_line(linewidth = 0.6, colour = "#d6604d") +
  facet_wrap(~ pair_label, ncol = 1, scales = "free_y") +
  labs(x = "Wavelength (nm)",
       y = "Mean reflectance difference (2018-cluster − 2025-cluster)",
       title = "Where the residual year-effect lives (L2-normalized, corrected spectra)",
       subtitle = "Same ecological community split into 2018-only vs 2025-only clusters") +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(face = "bold"),
        plot.title.position = "plot")
ggsave("docs/figures/year_separated_pair_diffs.pdf", p_diff,
       width = 11, height = 9, device = cairo_pdf)

# Plot 2: overlaid mean spectra per pair.
overlaid <- per_pair |>
  dplyr::select(pair_label, wavelength_nm, mean_18, mean_25) |>
  tidyr::pivot_longer(c(mean_18, mean_25),
                      names_to = "year", values_to = "rfl",
                      names_prefix = "mean_") |>
  dplyr::mutate(year = dplyr::recode(year,
                                      "18" = "2018-cluster",
                                      "25" = "2025-cluster"))
p_mean <- ggplot(overlaid,
                  aes(x = wavelength_nm, y = rfl, colour = year)) +
  geom_rect(data = water_bands, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey90", alpha = 0.6) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c("2018-cluster" = "#1b7837",
                                  "2025-cluster" = "#762a83")) +
  facet_wrap(~ pair_label, ncol = 1) +
  labs(x = "Wavelength (nm)",
       y = "L2-normalized reflectance (cluster mean)",
       title = "Cluster mean spectra for each year-separated pair",
       subtitle = "Lines overlap = no residual year-effect; visible offsets = where it lives") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold"),
        plot.title.position = "plot")
ggsave("docs/figures/year_separated_pair_means.pdf", p_mean,
       width = 11, height = 9, device = cairo_pdf)

cat("\nWrote docs/figures/year_separated_pair_diffs.pdf\n")
cat("Wrote docs/figures/year_separated_pair_means.pdf\n")
