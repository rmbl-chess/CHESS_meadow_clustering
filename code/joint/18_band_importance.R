# 18_band_importance.R — which wavelengths drive the joint classification?
#
# Strategy (see also the accompanying notes):
#   1. Permutation importance of each spectral PC in the DEPLOYED joint RF
#      (22 features = 20 PC + snow-free DOY + CHM; unweighted, as in 09).
#      Permutation (not impurity) because it is ~scale-invariant, and the PCs
#      are orthogonal so importance is not diluted by collinearity the way a
#      raw-band importance would be. Averaged over several seeds for stability.
#   2. Reproject PC importance onto the 348 water-masked bands via SQUARED
#      loadings: band_imp[b] = sum_j imp_j * v_j[b]^2. Loadings are orthonormal
#      (sum_b v_j[b]^2 = 1), so this preserves the total importance and gives a
#      per-band share; sign is dropped (magnitude only).
#   3. Summarise per band (+ wavelength) and per spectral region; plot.
#
# Because the deployed model drops the 6 narrow-band indices, the 20 PCs carry
# ALL spectral signal — the reprojection fully accounts for spectral importance.
#
# Inputs:
#   data/derived/joint_training_set.rds     (training, feature_cols)
#   data/derived/spectral_features.rds      (pca$rotation, keep_wl, n_pc)
# Outputs:
#   data/derived/band_importance.csv        (per band: wavelength, importance, %)
#   data/derived/band_importance_regions.csv
#   docs/figures/band_importance.pdf

suppressPackageStartupMessages({
  library(tidyverse)
  library(ranger)
})
dir.create("docs/figures", showWarnings = FALSE, recursive = TRUE)

jt   <- readRDS("data/derived/joint_training_set.rds")
sf   <- readRDS("data/derived/spectral_features.rds")
n_pc <- sf$n_pc                                  # 20
pc_cols <- sprintf("spec_PC%02d", seq_len(n_pc))
# Deployed feature set: PCs + DOY + CHM (the 6 indices are dropped at inference).
feat <- c(pc_cols, "snow_free_doy", "canopy_height_m")

train <- jt$training |>
  dplyr::select(final_label, dplyr::all_of(feat)) |>
  tidyr::drop_na()
X <- as.matrix(train[, feat]); y <- factor(train$final_label)
cat(sprintf("Permutation importance: %d sites, %d features, %d classes\n",
            nrow(X), ncol(X), nlevels(y)))

# --- Step 1: permutation importance, averaged over seeds --------------------
seeds <- 1:10
imp_mat <- vapply(seeds, function(s) {
  fit <- ranger::ranger(x = X, y = y, num.trees = 1000, importance = "permutation",
                        classification = TRUE, seed = s, num.threads = 0)
  fit$variable.importance[feat]
}, numeric(length(feat)))
rownames(imp_mat) <- feat
imp <- rowMeans(imp_mat)
imp[imp < 0] <- 0                                # tiny negatives -> 0 (noise)
imp_pc <- imp[pc_cols]

cat("\nNon-spectral feature importance (context):\n")
cat(sprintf("  snow_free_doy = %.4f (%.1f%% of total)\n",
            imp["snow_free_doy"], 100 * imp["snow_free_doy"] / sum(imp)))
cat(sprintf("  canopy_height_m = %.4f (%.1f%% of total)\n",
            imp["canopy_height_m"], 100 * imp["canopy_height_m"] / sum(imp)))
cat(sprintf("  20 spectral PCs together = %.1f%% of total importance\n",
            100 * sum(imp_pc) / sum(imp)))
cat("\nTop PCs by permutation importance:\n")
print(round(sort(imp_pc, decreasing = TRUE)[1:8], 4))

# --- Step 2: reproject to bands via squared loadings ------------------------
V  <- sf$pca$rotation[, seq_len(n_pc), drop = FALSE]   # 348 x 20, orthonormal cols
V2 <- V^2                                              # fraction of each PC per band
band_imp <- as.numeric(V2 %*% imp_pc)                  # 348-vector
communality <- rowSums(V2)                             # sum_{j<=20} v_j[b]^2 (<=1)
stopifnot(abs(sum(band_imp) - sum(imp_pc)) < 1e-8)     # total preserved

bands <- tibble::tibble(
  band          = seq_along(sf$keep_wl),
  wavelength_nm = sf$keep_wl,
  importance    = band_imp,
  importance_pct = 100 * band_imp / sum(band_imp),
  communality   = communality
) |> dplyr::arrange(wavelength_nm)

# --- Bootstrap uncertainty on the per-band importance -----------------------
# Resample sites with replacement (B times), refit, permutation importance ->
# reproject -> per-band % share. The SD across resamples is the sampling
# uncertainty of the importance -- wider and more honest than a CV-fold SD,
# whose 80%-overlapping training folds are highly correlated and understate it.
boot_band_profiles <- function(B = 40L) {
  n <- nrow(X); out <- list()
  for (b in seq_len(B)) {
    set.seed(3000L + b)
    ii <- sample(n, n, replace = TRUE); yt <- droplevels(y[ii])
    if (nlevels(yt) < 2) next
    fit <- ranger::ranger(x = X[ii, , drop = FALSE], y = yt, num.trees = 800,
                          importance = "permutation", classification = TRUE,
                          seed = 3000L + b, num.threads = 0)
    ipc <- pmax(0, fit$variable.importance[pc_cols])
    bi  <- as.numeric(V2 %*% ipc)
    out[[length(out) + 1L]] <- 100 * bi / sum(bi)
  }
  do.call(cbind, out)                                    # 348 x B
}
bootM <- boot_band_profiles()
bands$sd_pct <- apply(bootM, 1, stats::sd)[bands$band]   # bands arranged by wl
cat(sprintf("Bootstrap uncertainty from %d resamples; boot-mean vs point corr = %.3f\n",
            ncol(bootM), stats::cor(rowMeans(bootM)[bands$band], bands$importance_pct)))

readr::write_csv(bands, "data/derived/band_importance.csv")
cat(sprintf("\nWrote data/derived/band_importance.csv (%d bands)\n", nrow(bands)))

# --- Step 3: region summary + plot ------------------------------------------
regions <- tibble::tribble(
  ~region,       ~lo,   ~hi,
  "VIS",          380,   700,
  "Red edge",     700,   750,
  "NIR",          750,  1300,
  "SWIR1",       1300,  1900,
  "SWIR2",       1900,  2401)
bands$region <- purrr::map_chr(bands$wavelength_nm, function(nm) {
  hit <- regions$region[nm >= regions$lo & nm < regions$hi]
  if (length(hit)) hit[1] else NA_character_
})
region_summary <- bands |>
  dplyr::group_by(region) |>
  dplyr::summarise(n_bands = dplyr::n(),
                   importance_pct = sum(importance_pct),
                   .groups = "drop") |>
  dplyr::arrange(dplyr::desc(importance_pct))
cat("\n=== Importance by spectral region ===\n")
print(as.data.frame(region_summary))
readr::write_csv(region_summary, "data/derived/band_importance_regions.csv")

cat("\nTop 10 individual wavelengths:\n")
print(bands |> dplyr::arrange(dplyr::desc(importance_pct)) |>
        dplyr::transmute(wavelength_nm = round(wavelength_nm), importance_pct = round(importance_pct, 2)) |>
        utils::head(10) |> as.data.frame())

# --- Companion panel: mean L2-normalized reflectance per year --------------
# Mean training spectrum (meadow + shrub) on the same 348 water-masked bands,
# split 2018 / 2025, so importance peaks can be read against spectral features.
mean_spectrum_by_year <- function(keep_wl) {
  vs  <- readRDS("data/derived/veg_spectra.rds")
  wln <- vs$wavelengths$center_wavelength_nm
  rc  <- sprintf("rfl_band_%d", seq_along(wln))
  wm  <- (wln >= 1340 & wln <= 1450) | (wln >= 1800 & wln <= 1950) | (wln > 2400)
  keep <- rc[!wm]; stopifnot(length(keep) == length(keep_wl))
  ss <- readRDS("data/derived/shrub_veg_spectra.rds")$joined
  refl <- dplyr::bind_rows(vs$joined[, c("site_number", "Year", keep)],
                           ss[, c("site_number", "Year", keep)]) |>
    dplyr::filter(Year %in% c(2018L, 2025L))
  # 2018 is CRBU-only, so restrict 2025 to CRBU too — the reflectance comparison
  # is then domain-matched (year effect, not the ALMO/UPTA composition mix).
  crbu25 <- readRDS("data/derived/spectra_2025.rds")$spectra |>
    dplyr::filter(domain == "CRBU") |> dplyr::pull(site_number) |> unique()
  refl <- dplyr::filter(refl, Year == 2018L |
                          (Year == 2025L & site_number %in% crbu25))
  purrr::map_dfr(c(2018L, 2025L), function(y)
    tibble::tibble(series = as.character(y), wavelength_nm = keep_wl,
                   reflectance = colMeans(as.matrix(refl[refl$Year == y, keep]),
                                          na.rm = TRUE)))
}
refl <- mean_spectrum_by_year(sf$keep_wl)

# --- Two-panel figure: mean spectrum (top) + band importance (bottom) -------
water_bands <- tibble::tibble(xmin = c(1340, 1800, 2400), xmax = c(1450, 1950, 2510))
p_ref <- "Mean reflectance (L2-normalized)"
p_imp <- "Band importance (% of spectral)"
plotdf <- dplyr::bind_rows(
  refl  |> dplyr::transmute(panel = p_ref, wavelength_nm, series, value = reflectance,
                            ymin = NA_real_, ymax = NA_real_),
  bands |> dplyr::transmute(panel = p_imp, wavelength_nm, series = "All samples",
                            value = importance_pct,
                            ymin = pmax(0, importance_pct - sd_pct),
                            ymax = importance_pct + sd_pct)
) |>
  dplyr::mutate(panel = factor(panel, levels = c(p_ref, p_imp)),
                series = factor(series, levels = c("2018", "2025", "All samples")),
                segment = findInterval(wavelength_nm, c(1395, 1875)))
cols <- c("2018" = "#1b9e77", "2025" = "#7570b3", "All samples" = "#2c7fb8")

p <- ggplot(plotdf, aes(wavelength_nm, value, colour = series)) +
  geom_rect(data = water_bands, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey90", alpha = 0.6) +
  geom_ribbon(data = function(d) dplyr::filter(d, panel == p_imp),
              aes(x = wavelength_nm, ymin = ymin, ymax = ymax, fill = series,
                  group = interaction(series, segment)),
              inherit.aes = FALSE, alpha = 0.35) +
  geom_line(aes(group = interaction(series, segment)), linewidth = 0.4) +
  facet_wrap(~ panel, ncol = 1, scales = "free_y", strip.position = "left") +
  scale_colour_manual(values = cols, name = NULL) +
  scale_fill_manual(values = cols, guide = "none") +
  labs(x = "Wavelength (nm)", y = NULL,
       title = "Per-band contribution to the joint classification",
       subtitle = paste0("RF permutation importance of 20 PCs reprojected via squared loadings; ",
                         "top = mean training spectra (2018 & 2025, both CRBU-only). ",
                         "Ribbon = ±1 SD across 40 bootstrap resamples. Grey = water-masked.")) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank(),
        strip.placement = "outside", strip.text.y.left = element_text(angle = 90),
        legend.position = "top")
ggsave("docs/figures/band_importance.pdf", p, width = 10, height = 7.5, device = cairo_pdf)
cat("\nWrote docs/figures/band_importance.pdf\n")
