# 05_preprocess_features.R — derive feature blocks for clustering.
#
# Spectral block (built from already-L2-normalized mean spectra in 04):
#   - Mask water absorption bands (1340-1450, 1800-1950, >2400 nm).
#   - Top N PCA components (broad-shape variance).
#   - Narrow-band indices: NDVI, NDWI, PRI, red-edge slope, CAI, NDLI.
#
# Composition block — built at two granularity levels:
#   - Species-level: one column per canonical "Genus species".
#   - Genus-level:   species summed within genus (Salix x6 -> Salix, etc.).
# Both include the 6 non-species cover categories (NPV, Bare, Other_Forb,
# Other_Graminoid, Other_Moss_Lichen, Other_Deciduous_Shrub).
# Both are rare-trimmed (drop features present in <5% of sites) and
# Hellinger-transformed (sqrt of relative abundance per site).
#
# Inputs:  data/derived/veg_spectra.rds
# Outputs:
#   data/derived/spectral_features.rds
#   data/derived/composition_species.rds
#   data/derived/composition_genus.rds

library(tidyverse)

vs <- readRDS("data/derived/veg_spectra.rds")
joined <- vs$joined
wl     <- vs$wavelengths

# ============================================================================
# SPECTRAL FEATURES
# ============================================================================

rfl_cols  <- grep("^rfl_band_", names(joined), value = TRUE)
band_nums <- as.integer(stringr::str_extract(rfl_cols, "\\d+$"))
band_wl   <- wl$center_wavelength_nm[match(band_nums, wl$band_number)]

water_mask <- (band_wl >= 1340 & band_wl <= 1450) |
              (band_wl >= 1800 & band_wl <= 1950) |
              (band_wl >  2400)
keep_cols <- rfl_cols[!water_mask]
keep_wl   <- band_wl[!water_mask]
message(sprintf("Spectral: %d bands retained after water mask (of %d total).",
                length(keep_cols), length(rfl_cols)))

spec_mat <- as.matrix(joined[, keep_cols])

# PCA decorrelates the high-correlation reflectance bands. Spectra are
# already L2-normalized in 04; here we just center.
#
# Fit the basis on the 2025 rows ONLY so it matches JPL's pre-computed
# 3 m PC mosaics (which were generated on the 2025-only basis the
# original aop_classifier_pca.csv was exported on). 2018 + 2026 spectra
# (2018 year-corrected in 04; 2026 extracted directly from 2025 AOP) are
# PROJECTED onto this basis so they contribute to clustering and the joint
# training set without dragging the basis itself. This keeps the inference
# COGs (data/derived/aop_classified/) coherent with the new RF.
is_2025 <- joined$Year == 2025L
pca   <- prcomp(spec_mat[is_2025, , drop = FALSE],
                 center = TRUE, scale. = FALSE)
# Project ALL rows (every year) onto the 2025-fit basis.
spec_centered <- sweep(spec_mat, 2, pca$center, FUN = "-")
n_pc  <- 20
spec_pcs <- (spec_centered %*% pca$rotation)[, seq_len(n_pc), drop = FALSE]
colnames(spec_pcs) <- sprintf("spec_PC%02d", seq_len(n_pc))

# Cumulative variance explained on the FIT data (2025 only).
varexp <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
message(sprintf(
  "Spectral PCA fit on %d 2025 rows; top %d PCs explain %.1f%% of fit variance.",
  sum(is_2025), n_pc, 100 * varexp[n_pc]
))
message(sprintf("Projected %d sites onto basis (%d 2018 + %d 2025 + %d 2026).",
                nrow(spec_mat), sum(joined$Year == 2018L), sum(is_2025),
                sum(joined$Year == 2026L)))

# Narrow-band indices. Pick the band whose center is closest to each target.
band_of <- function(target_nm) which.min(abs(keep_wl - target_nm))
b <- function(nm) spec_mat[, band_of(nm)]

ndvi <- (b(860) - b(660)) / (b(860) + b(660))
ndwi <- (b(860) - b(1240)) / (b(860) + b(1240))
pri  <- (b(531) - b(570)) / (b(531) + b(570))
red_edge_slope <- (b(750) - b(700)) / (keep_wl[band_of(750)] - keep_wl[band_of(700)])
cai  <- 0.5 * (b(2000) + b(2200)) - b(2100)          # Cellulose Absorption Index
ndli <- (log(1 / b(1754)) - log(1 / b(1680))) /
        (log(1 / b(1754)) + log(1 / b(1680)))         # Normalized Lignin Index

spectral_features <- bind_cols(
  joined |> select(site_number, Year),
  as_tibble(spec_pcs),
  tibble(ndvi = ndvi, ndwi = ndwi, pri = pri,
         red_edge_slope = red_edge_slope, cai = cai, ndli = ndli)
)

saveRDS(list(features    = spectral_features,
             pca         = pca,
             keep_wl     = keep_wl,
             n_pc        = n_pc,
             var_explained = varexp),
        "data/derived/spectral_features.rds")

# ============================================================================
# COMPOSITION FEATURES
# ============================================================================

cover_cols <- grep("_cover$", names(joined), value = TRUE)
nonsp_cats <- c("Other_Forb", "Other_Graminoid", "NPV", "Bare",
                "Other_Moss_Lichen", "Other_Deciduous_Shrub")

long_comp <- joined |>
  select(site_number, Year, all_of(cover_cols)) |>
  pivot_longer(all_of(cover_cols), names_to = "feature", values_to = "cover") |>
  filter(cover > 0) |>
  mutate(
    name  = stringr::str_replace(feature, "_cover$", ""),
    is_species = !name %in% nonsp_cats,
    # For species, take first underscore-separated token as genus (e.g.,
    # "Salix_boothii" -> "Salix"). For non-species, the genus key is the
    # category name so the column survives genus aggregation intact.
    genus = if_else(is_species, stringr::word(name, 1, sep = "_"), name)
  )

prep_composition <- function(long_df, key, min_prev = 0.05) {
  total_sites <- length(unique(long_df$site_number))
  prev <- long_df |>
    group_by(.feat = .data[[key]]) |>
    summarise(n_sites = n_distinct(site_number), .groups = "drop") |>
    mutate(prev = n_sites / total_sites) |>
    rename(!!key := .feat)

  keep <- prev[[key]][prev$prev >= min_prev]
  long_kept <- long_df |> filter(.data[[key]] %in% keep)

  agg <- long_kept |>
    group_by(site_number, Year, .feat = .data[[key]]) |>
    summarise(cover = sum(cover), .groups = "drop") |>
    rename(!!key := .feat)

  wide <- agg |>
    mutate(col = paste0(.data[[key]], "_cover")) |>
    select(-all_of(key)) |>
    pivot_wider(names_from = col, values_from = cover, values_fill = 0)

  feat_cols <- setdiff(names(wide), c("site_number", "Year"))
  rel <- as.matrix(wide[, feat_cols])
  totals <- rowSums(rel)
  totals[totals == 0] <- NA_real_
  hellinger <- sqrt(rel / totals)
  hell_df <- bind_cols(wide |> select(site_number, Year), as_tibble(hellinger))

  # Drop rows where total cover was 0 (shouldn't happen after prevalence trim,
  # but guard anyway).
  hell_df <- hell_df |> filter(rowSums(across(all_of(feat_cols))) > 0)

  list(wide = wide, hellinger = hell_df,
       feature_cols = feat_cols, prevalence = prev)
}

comp_species <- prep_composition(long_comp, "name")
comp_genus   <- prep_composition(long_comp, "genus")

message(sprintf(
  "Composition after rare-drop (>=5%% of sites): species=%d features, genus=%d features.",
  length(comp_species$feature_cols),
  length(comp_genus$feature_cols)))

saveRDS(comp_species, "data/derived/composition_species.rds")
saveRDS(comp_genus,   "data/derived/composition_genus.rds")
