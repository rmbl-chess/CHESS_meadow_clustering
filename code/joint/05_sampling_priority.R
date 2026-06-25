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
# Vegetation-cover gate: raw leverage rewards spectral novelty that is often
# sparse bare / rock / dry pixels (leverage ~ -0.5 correlated with NDVI), which
# crews won't sample (target is >= cover_min_pct live cover). An NDVI->cover
# calibration from the training plots gives `pct_cover_est`; `leverage_gated`
# zeroes sub-threshold pixels so the fieldwork ranking and per-class candidate
# list only surface vegetated targets. Raw `leverage` is kept for context.
#
# Outputs a point GeoPackage and a per-class top-K candidate site list
# so a field crew can drop the gpkg in QGIS, pick top candidates per
# class, and plan a sampling route.
#
# Inputs:
#   data/derived/inference_pixel_distances.csv
#   data/derived/punch_list.csv
#   data/derived/inference_predictions.csv   (has per-pixel ndvi)
#   data/derived/cover_combined.rds, spectral_features.rds  (NDVI->cover fit)
# Outputs:
#   data/derived/sampling_priority.gpkg     (EPSG:32613; adds pct_cover_est,
#                                            meets_cover_min, leverage_gated)
#   data/derived/sampling_priority_top.csv  (top-K per class, cover-gated)

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
})

# --- Vegetation-cover gate -------------------------------------------------
# Raw leverage rewards spectral novelty, which is often sparse bare / rock /
# dry pixels (leverage ~ -0.5 correlated with NDVI) — not valid field targets
# (crews only sample sites with >= cover_min_pct live plant cover). Calibrate
# NDVI -> live-cover from the training plots (which carry measured % cover) and
# gate priority by an estimated cover. CAVEAT: training plots are all well-
# vegetated (lowest ~35% cover), so the 25% boundary is an EXTRAPOLATION
# (R^2 ~ 0.63) — `pct_cover_est` is approximate; tune `cover_min_pct` and
# inspect before trusting the exact cutoff.
cover_min_pct <- 25
cc <- readRDS("data/derived/cover_combined.rds")
.cov_cols <- grep("_cover$", names(cc), value = TRUE)
.tot  <- rowSums(cc[, .cov_cols], na.rm = TRUE)
.live <- rowSums(cc[, setdiff(.cov_cols, c("NPV_cover", "Bare_cover"))], na.rm = TRUE)
cover_cal <- tibble::tibble(site_number = cc$site_number, Year = cc$Year,
                            live_frac = pmin(100, 100 * .live / pmax(.tot, 1e-9))) |>
  dplyr::inner_join(readRDS("data/derived/spectral_features.rds")$features |>
                      dplyr::select(site_number, Year, ndvi),
                    by = c("site_number", "Year"))
cover_fit <- stats::lm(live_frac ~ ndvi, data = cover_cal)
ndvi_at_min <- unname((cover_min_pct - coef(cover_fit)[1]) / coef(cover_fit)[2])
cat(sprintf("NDVI->cover fit: live_frac = %.1f + %.1f*NDVI (R2=%.2f); %d%% cover ~ NDVI %.3f\n",
            coef(cover_fit)[1], coef(cover_fit)[2], summary(cover_fit)$r.squared,
            cover_min_pct, ndvi_at_min))

dist <- readr::read_csv("data/derived/inference_pixel_distances.csv",
                        show_col_types = FALSE)
punch <- readr::read_csv("data/derived/punch_list.csv",
                         show_col_types = FALSE)
pred  <- readr::read_csv("data/derived/inference_predictions.csv",
                         show_col_types = FALSE)

# Meadow class narratives (S01..S26, M01..M05) live in small_reference.
# Shrub labels are species binomials and need no extra description.
narrs <- readr::read_csv("data/small_reference/label_community_names.csv",
                         show_col_types = FALSE) |>
  dplyr::transmute(predicted_label = final_label,
                   class_description = dplyr::coalesce(narrative_curated,
                                                       narrative_draft))

# nearest_class from 39 and predicted_label from 38 don't have to match
# (39 uses centroid distance; 38 uses RF). Use the RF prediction as the
# canonical "predicted class" for leverage scoring — that's the class
# that would actually appear on the final map.
joined <- dist |>
  dplyr::select(x_utm, y_utm, domain, nearest_class, nearest_d,
                second_class, second_d, margin, ood_flag) |>
  dplyr::bind_cols(
    pred |> dplyr::select(predicted_label, snow_free_doy, canopy_height_m, ndvi)
  ) |>
  dplyr::left_join(
    punch |> dplyr::select(final_label, n_total, predicted_n_pixels,
                            balanced_recall, augmentation_priority),
    by = c("predicted_label" = "final_label")
  ) |>
  dplyr::left_join(narrs, by = "predicted_label") |>
  dplyr::mutate(
    # For shrub labels (species binomials), use the label itself as the
    # description so every row has a human-readable name.
    class_description = dplyr::coalesce(class_description, predicted_label)
  )

# --- Leverage score + vegetation-cover gate ------------------------------
joined <- joined |>
  dplyr::mutate(
    leverage = nearest_d / sqrt(pmax(n_total, 1L)),
    # Estimated live plant cover from the NDVI calibration above.
    pct_cover_est = pmin(100, pmax(0, as.numeric(
      stats::predict(cover_fit, newdata = data.frame(ndvi = ndvi))))),
    meets_cover_min = !is.na(pct_cover_est) & pct_cover_est >= cover_min_pct,
    # Gated leverage: sites below the cover minimum are not field targets, so
    # they drop to the bottom of the ranking (raw `leverage` kept for context).
    leverage_gated = dplyr::if_else(meets_cover_min, leverage, 0),
    novelty_rank          = dplyr::min_rank(dplyr::desc(nearest_d)),
    class_scarcity_rank   = dplyr::min_rank(n_total),
    leverage_rank         = dplyr::min_rank(dplyr::desc(leverage)),
    leverage_gated_rank   = dplyr::min_rank(dplyr::desc(leverage_gated))
  )

cat(sprintf("Inference pixels: %d\n", nrow(joined)))
cat(sprintf("Below %d%% est. live cover (gated out of fieldwork ranking): %d / %d (%.1f%%)\n",
            cover_min_pct, sum(!joined$meets_cover_min), nrow(joined),
            100 * mean(!joined$meets_cover_min)))
cat("\nGated-leverage distribution (vegetated targets only):\n")
print(summary(joined$leverage_gated[joined$meets_cover_min]))

cat("\nTop 20 highest GATED-leverage pixels (>= cover min; best next field sites):\n")
print(joined |>
        dplyr::filter(meets_cover_min) |>
        dplyr::arrange(dplyr::desc(leverage_gated)) |>
        dplyr::select(x_utm, y_utm, domain, predicted_label, n_total,
                      balanced_recall, nearest_d, pct_cover_est, leverage,
                      augmentation_priority) |>
        utils::head(20) |> as.data.frame())

# --- Top-K candidates per predicted class (cover-gated) -----------------
# Only vegetated sites (>= cover_min_pct) are real field targets.
top_per_class <- 10L
top_candidates <- joined |>
  dplyr::filter(meets_cover_min) |>
  dplyr::group_by(predicted_label) |>
  dplyr::arrange(dplyr::desc(leverage_gated), .by_group = TRUE) |>
  dplyr::slice_head(n = top_per_class) |>
  dplyr::ungroup() |>
  dplyr::arrange(augmentation_priority, dplyr::desc(leverage_gated))
cat(sprintf("\nTop-%d cover-gated candidates per predicted class: %d rows across %d classes\n",
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
