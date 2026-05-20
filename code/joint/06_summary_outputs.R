# 06_summary_outputs.R — joint phase deliverable bundle.
#
# Two artifacts in one script (they share inputs and one feeds the other):
#
#   data/derived/class_summary_table.csv   one row per class with N,
#                                          prevalence, recall, indicator
#                                          + abundant taxa, description,
#                                          median leverage, augmentation
#                                          priority.
#   data/derived/joint_training.gpkg       two layers (training_sites_crowns
#                                          / _points) with class metadata
#                                          for QGIS review.
#
# Inputs:
#   data/derived/punch_list.csv
#   data/derived/sampling_priority.gpkg
#   data/derived/joint_training_set.rds
#   data/derived/crowns_2018.gpkg, crowns_2025.gpkg
#   data/derived/label_descriptions.csv          (meadow IndVal indicators)
#   data/derived/shrub_training_set.rds          (per-record canonical info)
#   data/small_reference/label_community_names.csv (meadow narratives)

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
})

# ============================================================================
# Section 1 — class_summary_table.csv
# ============================================================================
punch <- readr::read_csv("data/derived/punch_list.csv",
                         show_col_types = FALSE)
spri  <- sf::st_read("data/derived/sampling_priority.gpkg", quiet = TRUE) |>
  sf::st_drop_geometry()
narrs <- readr::read_csv("data/small_reference/label_community_names.csv",
                         show_col_types = FALSE) |>
  dplyr::transmute(final_label,
                   description = dplyr::coalesce(narrative_curated,
                                                 narrative_draft),
                   short_label)

# Meadow side: full IndVal indicator + abundant + physiognomy strings.
meadow_desc <- readr::read_csv("data/derived/label_descriptions.csv",
                               show_col_types = FALSE) |>
  dplyr::transmute(final_label,
                   indicator_taxa = indicators,
                   abundant_taxa  = abundant,
                   physiognomy    = physiognomy)

# Shrub side: "binomial (n=X, %)" listing built from canonical binomials.
shrub_train <- readRDS("data/derived/shrub_training_set.rds")$training
shrub_desc <- shrub_train |>
  dplyr::count(final_label, canonical_binomial, name = "n_sites") |>
  dplyr::group_by(final_label) |>
  dplyr::mutate(pct = round(100 * n_sites / sum(n_sites), 1)) |>
  dplyr::arrange(dplyr::desc(n_sites), .by_group = TRUE) |>
  dplyr::summarise(
    indicator_taxa = paste(
      sprintf("%s (n=%d, %0.1f%%)", canonical_binomial, n_sites, pct),
      collapse = "; "),
    .groups = "drop"
  ) |>
  dplyr::mutate(abundant_taxa = indicator_taxa,
                physiognomy   = "shrub crown")

class_extra <- dplyr::bind_rows(meadow_desc, shrub_desc) |>
  dplyr::left_join(narrs, by = "final_label") |>
  dplyr::mutate(
    description = dplyr::coalesce(description, final_label),
    # Shrub classes don't go through label_community_names.csv, so build
    # their short label here. Format: "Shrub - {binomial}"; convert the
    # collapsed-genus labels (e.g., "Salix sp.") to the same shape.
    short_label = dplyr::case_when(
      !is.na(short_label)              ~ short_label,
      TRUE                              ~ paste("Shrub -", final_label)
    )
  )

med_lev <- spri |>
  dplyr::group_by(predicted_label) |>
  dplyr::summarise(median_leverage = stats::median(leverage, na.rm = TRUE),
                   .groups = "drop") |>
  dplyr::rename(final_label = predicted_label)

summary_tbl <- punch |>
  dplyr::transmute(
    final_label, class_type, n_2018, n_2025, n_total,
    predicted_n_pixels,
    pct_of_inference = round(100 * predicted_n_pixels /
                              sum(predicted_n_pixels, na.rm = TRUE), 2),
    balanced_recall  = round(balanced_recall, 2),
    augmentation_priority, top_confusions
  ) |>
  dplyr::left_join(med_lev,     by = "final_label") |>
  dplyr::left_join(class_extra, by = "final_label") |>
  dplyr::select(final_label, short_label, class_type, description,
                indicator_taxa, abundant_taxa, physiognomy,
                n_2018, n_2025, n_total,
                predicted_n_pixels, pct_of_inference,
                balanced_recall, median_leverage,
                augmentation_priority, top_confusions) |>
  dplyr::arrange(
    factor(augmentation_priority,
           levels = c("critical", "high", "medium", "ok")),
    dplyr::desc(median_leverage)
  ) |>
  dplyr::mutate(median_leverage = round(median_leverage, 2))

readr::write_csv(summary_tbl, "data/derived/class_summary_table.csv")
cat(sprintf("Wrote data/derived/class_summary_table.csv (%d rows)\n",
            nrow(summary_tbl)))

# ============================================================================
# Section 2 — joint_training.gpkg (two layers)
# ============================================================================
js    <- readRDS("data/derived/joint_training_set.rds")
train <- js$training

# Keep only the columns we want in the gpkg join (the summary_tbl already
# has class_type; drop it from train to avoid class_type.x / .y).
class_summ <- summary_tbl |>
  dplyr::select(final_label, class_type,
                class_description = description,
                indicator_taxa, abundant_taxa,
                n_total, predicted_n_pixels, pct_of_inference,
                balanced_recall, median_leverage, augmentation_priority)
train <- train |> dplyr::select(-class_type)
train_enriched <- train |> dplyr::left_join(class_summ, by = "final_label")

crowns_2018 <- sf::st_read("data/derived/crowns_2018.gpkg", quiet = TRUE) |>
  dplyr::mutate(Year = 2018L) |>
  dplyr::select(site_number, Year)
crowns_2025 <- sf::st_read("data/derived/crowns_2025.gpkg", quiet = TRUE) |>
  dplyr::mutate(Year = 2025L) |>
  dplyr::select(site_number, Year)
crowns <- dplyr::bind_rows(crowns_2018, crowns_2025) |>
  sf::st_transform(32613)

# Union multi-crown sites so the polygon layer is one row per (site, Year).
crowns_one_per_site <- crowns |>
  dplyr::group_by(site_number, Year) |>
  dplyr::summarise(.groups = "drop")

crowns_layer <- crowns_one_per_site |>
  dplyr::inner_join(train_enriched, by = c("site_number", "Year"))
cat(sprintf("training_sites_crowns layer: %d polygons\n", nrow(crowns_layer)))

crown_pts <- sf::st_centroid(crowns_layer)
missing_geom <- train_enriched |>
  dplyr::anti_join(sf::st_drop_geometry(crowns_layer),
                   by = c("site_number", "Year"))
if (nrow(missing_geom) > 0) {
  cat(sprintf("Note: %d training rows have no crown polygon; dropped from gpkg\n",
              nrow(missing_geom)))
}
points_layer <- crown_pts
cat(sprintf("training_sites_points layer: %d points\n", nrow(points_layer)))

sf::st_write(crowns_layer, "data/derived/joint_training.gpkg",
             layer = "training_sites_crowns",
             delete_dsn = TRUE, quiet = TRUE)
sf::st_write(points_layer, "data/derived/joint_training.gpkg",
             layer = "training_sites_points",
             append = TRUE, quiet = TRUE)
cat("\nWrote data/derived/joint_training.gpkg (2 layers)\n")

cat("\n=== Class breakdown in gpkg ===\n")
print(as.data.frame(sf::st_drop_geometry(crowns_layer) |>
        dplyr::count(class_type, augmentation_priority, name = "n_sites") |>
        tidyr::pivot_wider(names_from = augmentation_priority,
                           values_from = n_sites, values_fill = 0L)))
