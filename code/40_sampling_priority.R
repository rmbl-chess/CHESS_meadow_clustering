# 40_sampling_priority.R — per-pixel sampling priority for next field
# campaigns. Combines two existing signals: (a) Mahalanobis novelty
# (nearest_d, from 39) and (b) training-N of the pixel's predicted class
# (from punch_list / 38), into a single leverage score per pixel:
#
#     leverage = nearest_d / sqrt(n_training_predicted_class)
#
# Rationale:
#   - nearest_d in the numerator: a pixel far from any training centroid
#     buys MORE information than one that's already well within an
#     existing class envelope (Mahalanobis is in z-scaled feature units).
#   - 1/sqrt(n_train) in the denominator: marginal value of one more
#     sample falls roughly as 1/sqrt(n) for many learners, so a pixel
#     predicted as a low-N class gives a bigger expected gain.
#   - Multiplicative not additive: a pixel that's both novel AND
#     predicted-as-undersampled scores well above either alone.
#
# Outputs a point GeoPackage and a per-class top-K candidate site list
# so a field crew can drop the gpkg in QGIS, pick top candidates per
# class, and plan a sampling route.
#
# Inputs:
#   data/derived/inference_pixel_distances.csv
#   data/derived/punch_list.csv
#   data/derived/inference_predictions.csv
# Outputs:
#   data/derived/sampling_priority.gpkg     (5354 points, EPSG:32613)
#   data/derived/sampling_priority_top.csv  (top-K per class)

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
})

dist <- readr::read_csv("data/derived/inference_pixel_distances.csv",
                        show_col_types = FALSE)
punch <- readr::read_csv("data/derived/punch_list.csv",
                         show_col_types = FALSE)
pred  <- readr::read_csv("data/derived/inference_predictions.csv",
                         show_col_types = FALSE)

# nearest_class from 39 and predicted_label from 38 don't have to match
# (39 uses centroid distance; 38 uses RF). Use the RF prediction as the
# canonical "predicted class" for leverage scoring — that's the class
# that would actually appear on the final map.
joined <- dist |>
  dplyr::select(x_utm, y_utm, domain, nearest_class, nearest_d,
                second_class, second_d, margin, ood_flag) |>
  dplyr::bind_cols(
    pred |> dplyr::select(predicted_label, snow_free_doy, canopy_height_m)
  ) |>
  dplyr::left_join(
    punch |> dplyr::select(final_label, n_total, predicted_n_pixels,
                            balanced_recall, augmentation_priority),
    by = c("predicted_label" = "final_label")
  )

# --- Leverage score ------------------------------------------------------
joined <- joined |>
  dplyr::mutate(
    leverage = nearest_d / sqrt(pmax(n_total, 1L)),
    novelty_rank          = dplyr::min_rank(dplyr::desc(nearest_d)),
    class_scarcity_rank   = dplyr::min_rank(n_total),
    leverage_rank         = dplyr::min_rank(dplyr::desc(leverage))
  )

cat(sprintf("Inference pixels: %d\n", nrow(joined)))
cat("\nLeverage distribution:\n")
print(summary(joined$leverage))

cat("\nTop 20 highest-leverage pixels (= best next field sites):\n")
print(joined |>
        dplyr::arrange(dplyr::desc(leverage)) |>
        dplyr::select(x_utm, y_utm, domain, predicted_label, n_total,
                      balanced_recall, nearest_d, margin, leverage,
                      augmentation_priority) |>
        utils::head(20) |> as.data.frame())

# --- Top-K candidates per predicted class -------------------------------
top_per_class <- 10L
top_candidates <- joined |>
  dplyr::group_by(predicted_label) |>
  dplyr::arrange(dplyr::desc(leverage), .by_group = TRUE) |>
  dplyr::slice_head(n = top_per_class) |>
  dplyr::ungroup() |>
  dplyr::arrange(augmentation_priority, dplyr::desc(leverage))
cat(sprintf("\nTop-%d candidates per predicted class: %d rows across %d classes\n",
            top_per_class, nrow(top_candidates),
            dplyr::n_distinct(top_candidates$predicted_label)))

# --- Persist ------------------------------------------------------------
joined_sf <- sf::st_as_sf(joined, coords = c("x_utm", "y_utm"),
                          crs = 32613, remove = FALSE)
sf::st_write(joined_sf, "data/derived/sampling_priority.gpkg",
             delete_dsn = TRUE, quiet = TRUE)
readr::write_csv(top_candidates, "data/derived/sampling_priority_top.csv")
cat("\nWrote data/derived/sampling_priority.gpkg\n")
cat("Wrote data/derived/sampling_priority_top.csv\n")
