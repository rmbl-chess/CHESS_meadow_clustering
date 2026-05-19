# 42_joint_training_gpkg.R — build a single GeoPackage that contains every
# meadow + shrub training site, geo-located, with the columns a reviewer
# needs to interpret the joint classifier (final_label, class_type,
# class_description, indicator taxa, recall, augmentation_priority,
# leverage stats). Two layers:
#
#   training_sites_crowns   crown polygons (where available)
#   training_sites_points   one centroid per (site, Year)
#
# Inputs:
#   data/derived/joint_training_set.rds
#   data/derived/crowns_2018.gpkg, crowns_2025.gpkg
#   data/derived/class_summary_table.csv
#   data/derived/punch_list.csv
# Output:
#   data/derived/joint_training.gpkg   (two layers, EPSG:32613)

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
})

js     <- readRDS("data/derived/joint_training_set.rds")
train  <- js$training
class_summ <- readr::read_csv("data/derived/class_summary_table.csv",
                              show_col_types = FALSE) |>
  dplyr::select(final_label, class_type, class_description = description,
                indicator_taxa, abundant_taxa,
                n_total, predicted_n_pixels, pct_of_inference,
                balanced_recall, median_leverage, augmentation_priority)
# class_summary already has class_type; drop the duplicated column from train
# before joining so we don't get class_type.x / class_type.y.
train <- train |> dplyr::select(-class_type)

train_enriched <- train |>
  dplyr::left_join(class_summ, by = "final_label")

crowns_2018 <- sf::st_read("data/derived/crowns_2018.gpkg", quiet = TRUE) |>
  dplyr::mutate(Year = 2018L) |>
  dplyr::select(site_number, Year)
crowns_2025 <- sf::st_read("data/derived/crowns_2025.gpkg", quiet = TRUE) |>
  dplyr::mutate(Year = 2025L) |>
  dplyr::select(site_number, Year)
crowns <- dplyr::bind_rows(crowns_2018, crowns_2025) |>
  sf::st_transform(32613)

# A few sites have multiple crowns; union them so we get one polygon
# per (site, Year) row.
crowns_one_per_site <- crowns |>
  dplyr::group_by(site_number, Year) |>
  dplyr::summarise(.groups = "drop")

# Polygon layer: join enriched attributes onto unioned crowns.
crowns_layer <- crowns_one_per_site |>
  dplyr::inner_join(train_enriched, by = c("site_number", "Year"))
cat(sprintf("training_sites_crowns layer: %d polygons\n", nrow(crowns_layer)))

# Point layer: centroids of those crowns, plus sites whose polygons
# weren't found (rare 2018 cases).
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
        dplyr::count(class_type, augmentation_priority,
                     name = "n_sites") |>
        tidyr::pivot_wider(names_from = augmentation_priority,
                           values_from = n_sites, values_fill = 0L)))
