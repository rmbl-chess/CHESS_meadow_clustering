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

# --- Per-band importance for one subset, with CV-fold profiles -------------
# Repeated stratified K-fold CV: each fold-model gives a per-band % profile; the
# point estimate is the mean and the uncertainty is the SD across fold-models.
cv_band_profiles <- function(df, K = 5L, R = 4L) {
  X <- as.matrix(df[, feat]); y <- factor(df$final_label); n <- nrow(X)
  out <- list()
  for (rep in seq_len(R)) {
    set.seed(2000L + rep)
    fold <- integer(n)
    for (lvl in levels(y)) { ii <- which(y == lvl)
      fold[ii] <- ((sample(seq_along(ii)) - 1L) %% K) + 1L }
    for (k in seq_len(K)) {
      tr <- which(fold != k); yt <- droplevels(y[tr])
      if (nlevels(yt) < 2) next
      fit <- ranger::ranger(x = X[tr, , drop = FALSE], y = yt, num.trees = 800,
                            importance = "permutation", classification = TRUE,
                            seed = 2000L + rep * 10L + k, num.threads = 0)
      ipc <- pmax(0, fit$variable.importance[pc_cols])
      bi  <- as.numeric(V2 %*% ipc)
      out[[length(out) + 1L]] <- 100 * bi / sum(bi)
    }
  }
  do.call(cbind, out)                                    # 348 x (K*R)
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
  M <- cv_band_profiles(df)
  tibble::tibble(subset = nm, band = seq_along(sf$keep_wl), wavelength_nm = sf$keep_wl,
                 importance_pct = rowMeans(M), sd_pct = apply(M, 1, stats::sd))
}) |>
  dplyr::mutate(subset = factor(subset, levels = c("2018", "2025", "Both years")))

# Wide by unique band index (keep_wl has a couple of duplicate nm values).
w <- prof |>
  tidyr::pivot_wider(id_cols = c(band, wavelength_nm),
                     names_from = subset, values_from = c(importance_pct, sd_pct))
readr::write_csv(w, "data/derived/band_importance_by_year.csv")

# How similar are the 2018 and 2025 profiles? (residual year fingerprint)
cat(sprintf("\n2018 vs 2025 per-band profile correlation: r = %.3f\n",
            stats::cor(w$`importance_pct_2018`, w$`importance_pct_2025`)))

# --- Companion panel: mean L2-normalized reflectance per year --------------
mean_spectrum_by_year <- function(keep_wl) {
  vs  <- readRDS("data/derived/veg_spectra.rds")
  wln <- vs$wavelengths$center_wavelength_nm
  rc  <- sprintf("rfl_band_%d", seq_along(wln))
  wm  <- (wln >= 1340 & wln <= 1450) | (wln >= 1800 & wln <= 1950) | (wln > 2400)
  keep <- rc[!wm]; stopifnot(length(keep) == length(keep_wl))
  ss <- readRDS("data/derived/shrub_veg_spectra.rds")$joined
  r  <- dplyr::bind_rows(vs$joined[, c("site_number", "Year", keep)],
                         ss[, c("site_number", "Year", keep)]) |>
    dplyr::filter(Year %in% c(2018L, 2025L))
  # 2018 is CRBU-only, so restrict 2025 to CRBU too — domain-matched comparison.
  crbu25 <- readRDS("data/derived/spectra_2025.rds")$spectra |>
    dplyr::filter(domain == "CRBU") |> dplyr::pull(site_number) |> unique()
  r <- dplyr::filter(r, Year == 2018L | (Year == 2025L & site_number %in% crbu25))
  purrr::map_dfr(c(2018L, 2025L), function(y)
    tibble::tibble(series = as.character(y), wavelength_nm = keep_wl,
                   reflectance = colMeans(as.matrix(r[r$Year == y, keep]), na.rm = TRUE)))
}
refl <- mean_spectrum_by_year(sf$keep_wl)

# --- Two-panel plot: mean spectrum (top) + importance by year (bottom) ------
water_bands <- tibble::tibble(xmin = c(1340, 1800, 2400), xmax = c(1450, 1950, 2510))
yr_cols <- c("2018" = "#1b9e77", "2025" = "#7570b3", "Both years" = "grey20")
p_ref <- "Mean reflectance (L2-normalized)"
p_imp <- "Band importance (% of spectral)"
plotdf <- dplyr::bind_rows(
  refl |> dplyr::transmute(panel = p_ref, wavelength_nm, series, value = reflectance,
                           ymin = NA_real_, ymax = NA_real_),
  prof |> dplyr::transmute(panel = p_imp, wavelength_nm, series = as.character(subset),
                           value = importance_pct,
                           ymin = pmax(0, importance_pct - sd_pct),
                           ymax = importance_pct + sd_pct)
) |>
  dplyr::mutate(panel = factor(panel, levels = c(p_ref, p_imp)),
                series = factor(series, levels = c("2018", "2025", "Both years")),
                segment = findInterval(wavelength_nm, c(1395, 1875)))

p <- ggplot(plotdf, aes(wavelength_nm, value, colour = series)) +
  geom_rect(data = water_bands, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey90", alpha = 0.6) +
  geom_ribbon(data = function(d) dplyr::filter(d, panel == p_imp),
              aes(x = wavelength_nm, ymin = ymin, ymax = ymax, fill = series,
                  group = interaction(series, segment)),
              inherit.aes = FALSE, alpha = 0.18) +
  geom_line(aes(group = interaction(series, segment)), linewidth = 0.55) +
  facet_wrap(~ panel, ncol = 1, scales = "free_y", strip.position = "left") +
  scale_colour_manual(values = yr_cols, name = "Training data") +
  scale_fill_manual(values = yr_cols, guide = "none") +
  labs(x = "Wavelength (nm)", y = NULL,
       title = "Per-band classification importance by campaign year",
       subtitle = paste0("Permutation importance of 20 PCs -> bands (squared loadings), ",
                         "fixed clusters, ", length(shared), " shared classes; ",
                         "top = mean spectra (2018 & 2025, both CRBU-only); ribbon = ±1 SD across ",
                         "5-fold CV. Grey = water-masked.")) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank(),
        strip.placement = "outside", strip.text.y.left = element_text(angle = 90),
        legend.position = "top")
ggsave("docs/figures/band_importance_by_year.pdf", p, width = 10.5, height = 7.5,
       device = cairo_pdf)
cat("Wrote docs/figures/band_importance_by_year.pdf\n")
cat("Wrote data/derived/band_importance_by_year.csv\n")
