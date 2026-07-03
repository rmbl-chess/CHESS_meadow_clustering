# 19_band_importance_by_year.R — per-wavelength importance, split by campaign
# year (2018 / 2025 / both), to see whether the two years' imagery leans on
# different wavelengths for the SAME communities.
#
# Method mirrors 18_band_importance.R (permutation importance of the 20 spectral
# PCs in the deployed 22-feature RF -> reprojected to bands via squared
# loadings), but:
#   - class labels are FIXED to the existing joint-training clusters
#     (final_label) so nothing is re-clustered per year;
#   - to keep the three profiles comparable we CONTROL for community
#     composition: all three RFs use only the classes present with >= MIN_N
#     sites in BOTH 2018 and 2025 (a shared task). Single-year classes (2025's
#     low-sagebrush, 2018's alpine) are dropped so wavelength differences
#     reflect the YEAR, not different communities.
#   - the PCA basis (2025-fit) and loadings are shared across all three, so
#     the reprojection is directly comparable.
#
# Inputs:  data/derived/joint_training_set.rds, spectral_features.rds
# Outputs: data/derived/band_importance_by_year.csv
#          docs/figures/band_importance_by_year.pdf

suppressPackageStartupMessages({
  library(tidyverse)
  library(ranger)
})
dir.create("docs/figures", showWarnings = FALSE, recursive = TRUE)

MIN_N <- 2L
SEEDS <- 1:10

jt   <- readRDS("data/derived/joint_training_set.rds")$training
sf   <- readRDS("data/derived/spectral_features.rds")
n_pc <- sf$n_pc
pc_cols <- sprintf("spec_PC%02d", seq_len(n_pc))
feat <- c(pc_cols, "snow_free_doy", "canopy_height_m")
V2   <- sf$pca$rotation[, seq_len(n_pc), drop = FALSE]^2   # 348 x 20

# --- Shared class set: >= MIN_N sites in BOTH 2018 and 2025 -----------------
yr <- jt |> dplyr::filter(Year %in% c(2018L, 2025L)) |>
  dplyr::count(final_label, Year) |>
  tidyr::pivot_wider(names_from = Year, values_from = n, values_fill = 0,
                     names_prefix = "y")
shared <- yr$final_label[yr$y2018 >= MIN_N & yr$y2025 >= MIN_N]
cat(sprintf("Shared classes (>= %d sites in both years): %d\n", MIN_N, length(shared)))

base <- jt |>
  dplyr::filter(Year %in% c(2018L, 2025L), final_label %in% shared) |>
  dplyr::select(final_label, Year, dplyr::all_of(feat)) |>
  tidyr::drop_na()

# --- Per-band importance for one row subset --------------------------------
band_importance <- function(df) {
  X <- as.matrix(df[, feat]); y <- factor(df$final_label)
  imp_mat <- vapply(SEEDS, function(s) {
    fit <- ranger::ranger(x = X, y = y, num.trees = 1000,
                          importance = "permutation", classification = TRUE,
                          seed = s, num.threads = 0)
    fit$variable.importance[pc_cols]
  }, numeric(n_pc))
  imp_pc <- pmax(0, rowMeans(imp_mat))
  band_imp <- as.numeric(V2 %*% imp_pc)
  100 * band_imp / sum(band_imp)             # share of spectral importance (%)
}

subsets <- list(
  "2018"       = base |> dplyr::filter(Year == 2018L),
  "2025"       = base |> dplyr::filter(Year == 2025L),
  "Both years" = base
)
for (nm in names(subsets))
  cat(sprintf("  %-11s %d sites, %d classes\n", nm, nrow(subsets[[nm]]),
              dplyr::n_distinct(subsets[[nm]]$final_label)))

prof <- purrr::imap_dfr(subsets, function(df, nm) {
  tibble::tibble(subset = nm, band = seq_along(sf$keep_wl),
                 wavelength_nm = sf$keep_wl, importance_pct = band_importance(df))
}) |>
  dplyr::mutate(subset = factor(subset, levels = c("2018", "2025", "Both years")))

# Wide by unique band index (keep_wl has a couple of duplicate nm values).
w <- prof |>
  tidyr::pivot_wider(id_cols = c(band, wavelength_nm),
                     names_from = subset, values_from = importance_pct)
readr::write_csv(w, "data/derived/band_importance_by_year.csv")

# How similar are the 2018 and 2025 profiles? (residual year fingerprint)
cat(sprintf("\n2018 vs 2025 per-band profile correlation: r = %.3f\n",
            stats::cor(w$`2018`, w$`2025`)))

# --- Plot: three profiles overlaid -----------------------------------------
prof <- prof |> dplyr::mutate(segment = findInterval(wavelength_nm, c(1395, 1875)))
water_bands <- tibble::tibble(xmin = c(1340, 1800, 2400), xmax = c(1450, 1950, 2510))
yr_cols <- c("2018" = "#1b9e77", "2025" = "#7570b3", "Both years" = "grey20")

p <- ggplot(prof, aes(wavelength_nm, importance_pct, colour = subset)) +
  geom_rect(data = water_bands, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey90", alpha = 0.6) +
  geom_line(aes(group = interaction(subset, segment)), linewidth = 0.55) +
  scale_colour_manual(values = yr_cols, name = "Training data") +
  labs(x = "Wavelength (nm)", y = "Share of spectral importance (%)",
       title = "Per-band classification importance by campaign year",
       subtitle = paste0("Permutation importance of 20 PCs -> bands (squared loadings), ",
                         "fixed clusters, ", length(shared), " shared classes. ",
                         "Grey = water-masked.")) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank(),
        legend.position = "top")
ggsave("docs/figures/band_importance_by_year.pdf", p, width = 10.5, height = 6,
       device = cairo_pdf)
cat("Wrote docs/figures/band_importance_by_year.pdf\n")
cat("Wrote data/derived/band_importance_by_year.csv\n")
