# 17_year_effect_pcs.R — characterize systematic 2018-vs-2025 spectral
# differences for sites that are otherwise similar.
#
# Method: for each 2018 site, find the closest 2025 site by Euclidean
# distance on the full Hellinger composition (named species + the 6
# non-species categories — so "similar composition" includes similar
# bare/NPV fractions). Compute per-feature differences (2025 - 2018) for
# spectral PCs and narrow-band indices on the matched pairs. Test which
# features differ systematically.
#
# Caveats:
#   - Nearest-neighbor matching is asymmetric: a 2025 site can be matched
#     to multiple 2018 sites. That's fine for diagnosing systematic shifts;
#     it's not for unbiased causal estimation.
#   - We filter to "good" matches (Hellinger distance below the median) to
#     avoid the long tail of pairs that look nothing alike.
#
# Inputs:  data/derived/{veg_spectra, composition_species, spectral_features}.rds
# Outputs: output/figures/year_effect_pcs.png
#          stdout: per-feature mean diff + t-stat table

suppressPackageStartupMessages({
  library(tidyverse)
})

cs        <- readRDS("data/derived/composition_species.rds")$hellinger
spec_feat <- readRDS("data/derived/spectral_features.rds")$features
spec_meta <- readRDS("data/derived/spectral_features.rds")

hell_cols <- setdiff(names(cs), c("site_number", "Year"))
H <- as.matrix(cs[, hell_cols])

idx_18 <- cs$Year == 2018L
idx_25 <- cs$Year == 2025L
H_18   <- H[idx_18, , drop = FALSE]
H_25   <- H[idx_25, , drop = FALSE]
keys_18 <- cs |> dplyr::filter(Year == 2018L) |>
  dplyr::select(site_number, Year)
keys_25 <- cs |> dplyr::filter(Year == 2025L) |>
  dplyr::select(site_number, Year)

# Cross-distance matrix (rows=2018, cols=2025). Use the identity
# ||a - b||² = ||a||² + ||b||² - 2 a·b to vectorize.
sq_norms_18 <- rowSums(H_18^2)
sq_norms_25 <- rowSums(H_25^2)
cross_d2    <- outer(sq_norms_18, sq_norms_25, "+") - 2 * H_18 %*% t(H_25)
cross_d2[cross_d2 < 0] <- 0
cross_d <- sqrt(cross_d2)

best_idx <- apply(cross_d, 1, which.min)
best_d   <- cross_d[cbind(seq_along(best_idx), best_idx)]

matches <- tibble::tibble(
  site_18   = keys_18$site_number,
  site_25   = keys_25$site_number[best_idx],
  hell_dist = best_d
)

threshold <- stats::median(matches$hell_dist)
matches_close <- matches |> dplyr::filter(hell_dist < threshold)
cat(sprintf("All pairs: %d   Close pairs (Hell dist < %.3f, median): %d\n",
            nrow(matches), threshold, nrow(matches_close)))

# Pull per-pair spectral features.
spec_18 <- spec_feat |>
  dplyr::filter(Year == 2018L) |>
  dplyr::select(site_number, dplyr::starts_with("spec_PC"),
                ndvi, ndwi, pri, red_edge_slope, cai, ndli)
spec_25 <- spec_feat |>
  dplyr::filter(Year == 2025L) |>
  dplyr::select(site_number, dplyr::starts_with("spec_PC"),
                ndvi, ndwi, pri, red_edge_slope, cai, ndli)

pairs <- matches_close |>
  dplyr::inner_join(spec_18, by = c("site_18" = "site_number")) |>
  dplyr::inner_join(spec_25, by = c("site_25" = "site_number"),
                    suffix = c("_18", "_25"))

all_features <- c(
  grep("^spec_PC", names(spec_feat), value = TRUE),
  intersect(c("ndvi", "ndwi", "pri", "red_edge_slope", "cai", "ndli"),
            names(spec_feat))
)
feature_cols <- all_features[paste0(all_features, "_18") %in% names(pairs) &
                             paste0(all_features, "_25") %in% names(pairs)]

diffs <- purrr::map_dfc(feature_cols, function(f) {
  v18 <- pairs[[paste0(f, "_18")]]
  v25 <- pairs[[paste0(f, "_25")]]
  d   <- v25 - v18
  tibble::tibble(!!f := d)
})

# Per-feature summary: mean, sd, paired t-stat, effect size standardized
# to the feature's own SD across all sites.
feat_sds <- vapply(feature_cols,
                   function(f) sd(spec_feat[[f]], na.rm = TRUE),
                   numeric(1))

per_feat <- purrr::map_dfr(feature_cols, function(f) {
  d <- diffs[[f]]
  if (all(is.na(d))) return(NULL)
  tt <- stats::t.test(d)
  tibble::tibble(
    feature       = f,
    mean_diff     = mean(d, na.rm = TRUE),
    sd_diff       = sd(d,   na.rm = TRUE),
    t_stat        = unname(tt$statistic),
    p_value       = tt$p.value,
    feature_sd    = feat_sds[f],
    abs_effect_sd = abs(mean(d, na.rm = TRUE)) / feat_sds[f]
  )
}) |>
  dplyr::arrange(dplyr::desc(abs_effect_sd))

cat("\n=== Per-feature systematic shift (2025 - 2018), sorted by |effect/SD| ===\n")
print(per_feat |> dplyr::mutate(dplyr::across(where(is.numeric),
                                              ~ signif(.x, 3))),
      n = Inf, width = Inf)

# Plot top differing features.
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

top_features <- per_feat |> dplyr::slice_head(n = 12) |> dplyr::pull(feature)

plot_df <- diffs |>
  dplyr::select(dplyr::all_of(top_features)) |>
  tidyr::pivot_longer(dplyr::everything(), names_to = "feature", values_to = "diff") |>
  dplyr::mutate(feature = factor(feature, levels = top_features))

p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = feature, y = diff)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_violin(fill = "steelblue", alpha = 0.4, scale = "width",
                       draw_quantiles = c(0.25, 0.5, 0.75)) +
  ggplot2::geom_jitter(width = 0.15, height = 0, alpha = 0.3, size = 0.7) +
  ggplot2::labs(
    x = NULL, y = "Value (2025) − Value (2018)",
    title = "Year effect on spectral features (matched pairs by composition)",
    subtitle = sprintf("%d matched pairs (Hellinger dist < median).  Mean ≠ 0 ⇒ systematic year shift.",
                       nrow(pairs))
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

ggplot2::ggsave("output/figures/year_effect_pcs.png", p,
                width = 12, height = 6, dpi = 150)
cat("\nWrote output/figures/year_effect_pcs.png\n")

# --- Loadings for top-shifting PCs (so we can see WHICH wavelengths) -------
top_pcs <- per_feat |>
  dplyr::filter(grepl("^spec_PC", feature)) |>
  dplyr::slice_head(n = 6) |>
  dplyr::pull(feature)

pca   <- spec_meta$pca
wl_nm <- spec_meta$keep_wl   # already in nm (from 05_preprocess_features.R)

loadings_long <- tibble::tibble(
  wavelength_nm = wl_nm
)
for (pc_name in top_pcs) {
  pc_idx <- as.integer(stringr::str_extract(pc_name, "\\d+"))
  loadings_long[[pc_name]] <- pca$rotation[, pc_idx]
}
loadings_long <- loadings_long |>
  tidyr::pivot_longer(-wavelength_nm, names_to = "PC", values_to = "loading") |>
  dplyr::mutate(PC = factor(PC, levels = top_pcs))

water_bands <- tibble::tibble(
  xmin = c(1340, 1800, 2400),
  xmax = c(1450, 1950, 2510)
)

p_load <- ggplot2::ggplot(loadings_long,
                          ggplot2::aes(x = wavelength_nm, y = loading, colour = PC)) +
  ggplot2::geom_rect(data = water_bands, inherit.aes = FALSE,
                     ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                     fill = "grey90", alpha = 0.6) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_line(linewidth = 0.6) +
  ggplot2::facet_wrap(~ PC, ncol = 2, scales = "free_y") +
  ggplot2::labs(x = "Wavelength (nm)", y = "PCA loading",
                title = "Wavelength loadings for top year-shifting PCs",
                subtitle = "Loadings show which spectral bands each PC weights; helps interpret the year-effect signature") +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(legend.position = "none",
                 strip.text = ggplot2::element_text(face = "bold"))

ggplot2::ggsave("output/figures/year_effect_pc_loadings.png", p_load,
                width = 12, height = 7, dpi = 150)
cat("Wrote output/figures/year_effect_pc_loadings.png\n")

# Save the raw differences and per-feature stats for follow-up.
saveRDS(list(matches = matches, matches_close = matches_close,
             diffs = diffs, per_feature = per_feat,
             threshold = threshold),
        "data/derived/year_effect.rds")
