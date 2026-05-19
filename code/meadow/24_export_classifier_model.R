# 24_export_classifier_model.R — package everything Python needs to extract
# AOP pixel features and predict classes:
#
#   aop_classifier_pca.csv          per-band PCA loadings + centering for
#                                    the 348 retained (non-water) bands
#   aop_classifier_indices.csv       target wavelengths for the 6 narrow-
#                                    band indices (NDVI, NDWI, PRI,
#                                    red_edge_slope, CAI, NDLI)
#   aop_classifier_meta.json         scalar config (water-band ranges,
#                                    preprocessing notes)
#
# The actual classifier RF is left to Python-side training from
# training_samples_sites.csv -- it's small and Python's RF will be fine.

suppressPackageStartupMessages({
  library(tidyverse)
  library(jsonlite)
})

spec_meta <- readRDS("data/derived/spectral_features.rds")
pca       <- spec_meta$pca
keep_wl   <- spec_meta$keep_wl   # nm, length 348

stopifnot(length(keep_wl) == nrow(pca$rotation),
          length(pca$center) == length(keep_wl))

n_pc <- 20

# --- 1. PCA model: loadings + center per retained band -----------------------
pca_df <- tibble::tibble(
  band_idx       = seq_along(keep_wl),
  wavelength_nm  = keep_wl,
  center         = as.numeric(pca$center)
)
loadings <- as.data.frame(pca$rotation[, seq_len(n_pc)])
colnames(loadings) <- sprintf("PC%02d", seq_len(n_pc))
pca_df <- dplyr::bind_cols(pca_df, loadings)

readr::write_csv(pca_df, "data/derived/aop_classifier_pca.csv")
cat(sprintf("Wrote aop_classifier_pca.csv: %d retained bands x %d PCs\n",
            nrow(pca_df), n_pc))

# --- 2. Index definitions --------------------------------------------------
indices_df <- tibble::tibble(
  index = c("ndvi", "ndwi", "pri", "red_edge_slope", "cai", "ndli"),
  formula = c(
    "(b860 - b660) / (b860 + b660)",
    "(b860 - b1240) / (b860 + b1240)",
    "(b531 - b570) / (b531 + b570)",
    "(b750 - b700) / (wl750 - wl700)",
    "0.5*(b2000 + b2200) - b2100",
    "(log(1/b1754) - log(1/b1680)) / (log(1/b1754) + log(1/b1680))"
  ),
  target_wavelengths_nm = c(
    "660,860", "860,1240", "531,570", "700,750",
    "2000,2100,2200", "1680,1754"
  )
)
readr::write_csv(indices_df, "data/derived/aop_classifier_indices.csv")
cat("Wrote aop_classifier_indices.csv\n")

# --- 3. Metadata / preprocessing config ------------------------------------
meta <- list(
  preprocessing = list(
    description = paste(
      "Per-pixel L2 normalization, then average across the 9 1m pixels in",
      "the 3x3 window centered on each 3m target. Apply water-band mask",
      "(see water_band_ranges_nm). Subtract PCA center, multiply by PCA",
      "rotation to get PC1..PC20."
    ),
    water_band_ranges_nm = list(c(1340, 1450), c(1800, 1950), c(2400, 9999)),
    no_data_sentinel     = -9000,
    pixel_window_m       = 3,   # 3x3 1m AOP pixels
    bands_total          = 426,
    bands_retained       = nrow(pca_df),
    n_pcs                = n_pc
  ),
  training_label_levels = sort(unique(
    readRDS("data/derived/final_clusters_B.rds")$final_summary$final_label
  )),
  notes = paste(
    "Classifier trained on training_samples_sites.csv with final_label as",
    "target. Features: spec_PC01..PC20, ndvi, ndwi, pri, red_edge_slope,",
    "cai, ndli, snow_free_doy."
  )
)
writeLines(jsonlite::toJSON(meta, pretty = TRUE, auto_unbox = TRUE),
           "data/derived/aop_classifier_meta.json")
cat("Wrote aop_classifier_meta.json\n")

cat("\nSummary:\n")
cat(sprintf("  Retained wavelength range: %.1f - %.1f nm\n",
            min(keep_wl), max(keep_wl)))
cat(sprintf("  PC1-%d explain %.1f%% of variance\n",
            n_pc, 100 * spec_meta$var_explained[n_pc]))
