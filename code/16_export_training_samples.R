# 16_export_training_samples.R — package the final clusters as training
# samples for downstream AOP classifier training.
#
# Three artefacts in data/derived/:
#
#   training_samples_sites.csv    one row per (site_number, Year). Contains:
#                                   - keys: site_number, Year
#                                   - labels: final_label, spec_cluster,
#                                             sub_cluster, indicator_species,
#                                             indicator_genus
#                                   - quality: recall, tier (strong /
#                                              marginal / weak)
#                                   - features: 20 spec PCs, 6 indices,
#                                               snow_free_doy
#
#   training_samples_crowns.gpkg  one row per crown polygon (1327 total
#                                  before filtering); same attributes as
#                                  the sites CSV plus geometry in
#                                  EPSG:32613. Use this to extract AOP
#                                  pixels and train.
#
#   training_labels_summary.csv   one row per final_label. Documents the
#                                  class: size, year breakdown, recall,
#                                  tier, indicator species/genus, top
#                                  features (composition profile), and
#                                  snow-free DOY range.
#
# Tiering scheme (recall = 5-fold CV with PCs + indices + snow_free_doy):
#   strong    recall ≥ 0.80   — classifier reliably finds these classes;
#                               canonical training pool.
#   marginal  0.50 ≤ recall   — usable but classifier confuses with
#                               neighbours; weight down or use carefully.
#   weak      recall < 0.50   — not recommended as a separate class.

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
})

fc        <- readRDS("data/derived/final_clusters_B.rds")
spec_feat <- readRDS("data/derived/spectral_features.rds")$features
env       <- readRDS("data/derived/environment.rds")
sc        <- readRDS("data/derived/spectral_clusters.rds")

# Pull the variant used by 10 (currently variant_F).
variant_used <- "variant_F"
features_in_clustering <- sc[[variant_used]]$pc_cols  # e.g. PC02..PC12 + DOY

# Feature columns to export. Keep the full PC set (1-20) plus indices and
# DOY — even though only a subset was used in clustering, the classifier
# may want the full feature space for prediction.
pc_cols   <- grep("^spec_PC", names(spec_feat), value = TRUE)
idx_cols  <- intersect(c("ndvi", "ndwi", "pri", "red_edge_slope",
                         "cai", "ndli"), names(spec_feat))

# --- Per-site attribute table -----------------------------------------------
tier_of <- function(r) dplyr::case_when(
  is.na(r) | r < 0.50 ~ "weak",
  r        < 0.80     ~ "marginal",
  TRUE                ~ "strong"
)

# Curated narrative labels (from data/small_reference/label_community_names.csv).
# Loaded here so each per-site row carries the narrative — that flows into the
# crowns GeoPackage and supports spatial review in QGIS / similar.
narr_path <- "data/small_reference/label_community_names.csv"
narratives <- if (file.exists(narr_path)) {
  readr::read_csv(narr_path, show_col_types = FALSE) |>
    dplyr::select(final_label,
                  dplyr::any_of(c("narrative_draft", "narrative_curated", "notes")))
} else {
  tibble::tibble(final_label = character())
}

label_meta <- fc$final_summary |>
  dplyr::transmute(
    final_label,
    indicator_species_label = indicator_species,
    label_top_features      = top_features,
    label_recall            = as.numeric(recall),
    label_n_sites           = n,
    label_n_2018            = n_2018,
    label_n_2025            = n_2025,
    tier                    = tier_of(label_recall)
  ) |>
  dplyr::left_join(narratives, by = "final_label")

asg_cols <- intersect(c("site_number", "Year", "spec_cluster", "sub_cluster",
                        "final_label", "source",
                        "inference_distance", "inference_gap"),
                      names(fc$assignments))
sites <- fc$assignments |>
  dplyr::select(dplyr::all_of(asg_cols)) |>
  dplyr::left_join(label_meta, by = "final_label") |>
  dplyr::inner_join(spec_feat |>
                      dplyr::select(site_number, Year,
                                    dplyr::all_of(pc_cols),
                                    dplyr::any_of(idx_cols)),
                    by = c("site_number", "Year")) |>
  dplyr::inner_join(env, by = c("site_number", "Year")) |>
  dplyr::arrange(final_label, site_number, Year)

readr::write_csv(sites, "data/derived/training_samples_sites.csv")
cat(sprintf("Wrote training_samples_sites.csv: %d sites, %d cols\n",
            nrow(sites), ncol(sites)))

# --- Per-label summary ------------------------------------------------------
env_per_label <- sites |>
  dplyr::group_by(final_label) |>
  dplyr::summarise(
    snow_free_doy_mean = mean(snow_free_doy, na.rm = TRUE),
    snow_free_doy_sd   = sd(snow_free_doy,   na.rm = TRUE),
    snow_free_doy_min  = min(snow_free_doy,  na.rm = TRUE),
    snow_free_doy_max  = max(snow_free_doy,  na.rm = TRUE),
    .groups = "drop"
  )

# Rich label descriptions (from 18_label_descriptions.R) — top indicator
# species, top abundant species, physiognomic profile. Optional: only join
# if the file is present.
desc_path <- "data/derived/label_descriptions.csv"
descriptions <- if (file.exists(desc_path)) {
  readr::read_csv(desc_path, show_col_types = FALSE)
} else {
  tibble::tibble(final_label = character())
}

# (Narratives already loaded above and joined into the per-site table.)

label_summary <- fc$final_summary |>
  dplyr::transmute(
    final_label,
    spec_cluster,
    n_sites = n,
    n_2018, n_2025,
    recall = as.numeric(recall),
    tier = tier_of(recall),
    indicator_species,
    top_features
  ) |>
  dplyr::left_join(env_per_label, by = "final_label") |>
  dplyr::left_join(descriptions, by = "final_label") |>
  dplyr::left_join(narratives,   by = "final_label") |>
  dplyr::arrange(dplyr::desc(recall), dplyr::desc(n_sites))

readr::write_csv(label_summary, "data/derived/training_labels_summary.csv")
cat(sprintf("Wrote training_labels_summary.csv: %d labels\n",
            nrow(label_summary)))

# --- Per-crown GeoPackage ---------------------------------------------------
crowns_2018 <- sf::st_read("data/derived/crowns_2018.gpkg", quiet = TRUE) |>
  dplyr::mutate(Year = 2018L) |>
  dplyr::select(site_number, Year)
crowns_2025 <- sf::st_read("data/derived/crowns_2025.gpkg", quiet = TRUE) |>
  dplyr::mutate(Year = 2025L) |>
  dplyr::select(site_number, Year)
crowns_all <- dplyr::bind_rows(crowns_2018, crowns_2025) |>
  sf::st_transform(32613)

crowns_labeled <- crowns_all |>
  dplyr::inner_join(sites, by = c("site_number", "Year"))

sf::st_write(crowns_labeled,
             "data/derived/training_samples_crowns.gpkg",
             delete_dsn = TRUE, quiet = TRUE)
cat(sprintf("Wrote training_samples_crowns.gpkg: %d crown polygons, %d cols\n",
            nrow(crowns_labeled), ncol(crowns_labeled)))

# --- Per-crown CENTROID point file (lighter-weight for display in QGIS) ----
crowns_points <- crowns_labeled |>
  sf::st_centroid()
sf::st_write(crowns_points,
             "data/derived/training_samples_points.gpkg",
             delete_dsn = TRUE, quiet = TRUE)
cat(sprintf("Wrote training_samples_points.gpkg: %d centroids, %d cols\n",
            nrow(crowns_points), ncol(crowns_points)))

# --- Console summary by tier ------------------------------------------------
cat("\n=== Tier breakdown ===\n")
tier_summary <- label_summary |>
  dplyr::group_by(tier) |>
  dplyr::summarise(n_labels = dplyr::n(),
                   n_sites  = sum(n_sites),
                   recall_min = min(recall),
                   recall_max = max(recall),
                   .groups = "drop") |>
  dplyr::mutate(tier = factor(tier, levels = c("strong","marginal","weak"))) |>
  dplyr::arrange(tier)
print(tier_summary)

cat("\nFiles in data/derived/:\n")
for (f in c("training_samples_sites.csv",
            "training_samples_crowns.gpkg",
            "training_samples_points.gpkg",
            "training_labels_summary.csv")) {
  p <- file.path("data/derived", f)
  if (file.exists(p)) {
    cat(sprintf("  %-40s  %.1f KB\n", f, file.size(p) / 1024))
  }
}
