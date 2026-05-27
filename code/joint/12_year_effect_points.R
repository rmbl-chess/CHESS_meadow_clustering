# 12_year_effect_points.R — stratified random point set in the CRBU domain
# for evaluating 2018-vs-2025 AOP year effects via direct re-extraction
# at the same locations.
#
# Two point types:
#   vegetated      R3D018 class 3 (meadow) + class 10 (deciduous shrub)
#                  with >= 80% purity in the 3 m cell. Use this set to
#                  test phenology differences between years (NDVI, red
#                  edge, NIR plateau, etc.).
#   non_vegetated  R3D018 class 6 (rock / bare soil / sparse veg) with
#                  >= 80% purity. Surfaces that don't grow back season-
#                  to-season — any spectral shift between 2018 and 2025
#                  is most likely instrument / radiometric calibration.
#
# Both strata are further binned by snow-free DOY (early / mid / late)
# so phenology / elevation isn't confounded with year. Sample size is
# ~n_per_stratum per (point_type x doy_band) cell.
#
# Inputs:
#   data/raw/SDP/UG_landcover_1m_v4.tif       (R3D018, 1 m)
#   data/derived/aop_chm_3m/CRBU_chm_max_3m.tif   (defines CRBU extent + 3 m grid)
#   R4D061 via rSDP                                (snow-free DOY)
# Outputs:
#   data/derived/year_effect_points.csv
#   data/derived/year_effect_points.gpkg

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(terra)
  library(rSDP)
})
terra::terraOptions(progress = 0)

domain_code      <- "CRBU"
purity_threshold <- 0.80
n_per_stratum    <- 500L
doy_breaks <- tibble::tribble(
  ~band,  ~min_doy, ~max_doy,
  "early",       0,    130,
  "mid",       130,    155,
  "late",      155,    400
)
set.seed(42)

# --- 1. Crop landcover to CRBU extent and reclassify -----------------------
lc_1m <- terra::rast("data/raw/SDP/UG_landcover_1m_v4.tif")
ref_3m <- terra::rast(
  file.path("data/derived/aop_chm_3m",
            sprintf("%s_chm_max_3m.tif", domain_code))
)
crbu_ext <- terra::ext(ref_3m)

cat(sprintf("Cropping R3D018 to %s extent ...\n", domain_code))
t0 <- Sys.time()
lc_dom <- terra::crop(lc_1m, crbu_ext)
cat(sprintf("  %d x %d px (%.1fs)\n",
            terra::nrow(lc_dom), terra::ncol(lc_dom),
            as.numeric(Sys.time() - t0, units = "secs")))

# Binary masks per stratum (1 where in target, 0 elsewhere).
veg_1m <- terra::classify(
  lc_dom,
  matrix(c(3, 3, 1, 10, 10, 1), ncol = 3, byrow = TRUE),
  others = 0L, include.lowest = TRUE, right = NA,
  datatype = "INT1U"
)
bare_1m <- terra::classify(
  lc_dom,
  matrix(c(6, 6, 1), ncol = 3, byrow = TRUE),
  others = 0L, include.lowest = TRUE, right = NA,
  datatype = "INT1U"
)

# Aggregate to 3 m by mean fraction in each 3x3 block, then threshold.
veg_frac  <- terra::aggregate(veg_1m,  fact = 3, fun = "mean", na.rm = TRUE)
bare_frac <- terra::aggregate(bare_1m, fact = 3, fun = "mean", na.rm = TRUE)
veg_3m  <- veg_frac  >= purity_threshold
bare_3m <- bare_frac >= purity_threshold

# Snap both to the CRBU CHM grid.
veg_3m  <- terra::resample(veg_3m,  ref_3m, method = "near")
bare_3m <- terra::resample(bare_3m, ref_3m, method = "near")

# --- 2. Candidate cell centers per stratum ---------------------------------
candidates <- function(mask, type_label) {
  idx <- which(terra::values(mask, mat = FALSE) == 1)
  if (length(idx) == 0) return(NULL)
  xy <- terra::xyFromCell(mask, idx)
  tibble::tibble(x_utm = xy[, 1], y_utm = xy[, 2],
                 point_type = type_label)
}
veg_pts  <- candidates(veg_3m,  "vegetated")
bare_pts <- candidates(bare_3m, "non_vegetated")
cat(sprintf("Candidates: %s vegetated, %s non-vegetated\n",
            format(nrow(veg_pts),  big.mark = ","),
            format(nrow(bare_pts), big.mark = ",")))

# Pre-sample to manageable size before the DOY extract.
prefilter_n <- 50000L
prefilter <- function(df, n) {
  if (nrow(df) > n) dplyr::slice_sample(df, n = n, replace = FALSE) else df
}
all_pts <- dplyr::bind_rows(
  prefilter(veg_pts,  prefilter_n),
  prefilter(bare_pts, prefilter_n)
)

# --- 3. Snow-free DOY at each candidate ------------------------------------
cat("Extracting snow-free DOY (R4D061) ... ")
t0 <- Sys.time()
r4d061 <- rSDP::sdp_get_raster("R4D061")
if (terra::nlyr(r4d061) > 1) r4d061 <- r4d061[[1]]
vals <- terra::extract(
  r4d061,
  terra::vect(all_pts, geom = c("x_utm", "y_utm"), crs = "EPSG:32613"),
  ID = FALSE
)[[1]]
cat(sprintf("done (%.1fs)\n",
            as.numeric(Sys.time() - t0, units = "secs")))

all_pts <- all_pts |>
  dplyr::mutate(snow_free_doy = vals) |>
  dplyr::filter(!is.na(snow_free_doy)) |>
  dplyr::mutate(doy_band = cut(snow_free_doy,
                                breaks = c(doy_breaks$min_doy,
                                           max(doy_breaks$max_doy)),
                                labels = doy_breaks$band,
                                include.lowest = TRUE))

cat("Candidate counts by point_type x doy_band:\n")
print(all_pts |> dplyr::count(point_type, doy_band) |>
        tidyr::pivot_wider(names_from = doy_band, values_from = n,
                           values_fill = 0L) |> as.data.frame())

# --- 4. Stratified random sample ------------------------------------------
sampled <- all_pts |>
  dplyr::group_by(point_type, doy_band) |>
  dplyr::group_split() |>
  purrr::map_dfr(function(df) {
    dplyr::slice_sample(df, n = min(n_per_stratum, nrow(df)),
                         replace = FALSE)
  })

cat(sprintf("\nFinal: %d points\n", nrow(sampled)))
print(sampled |> dplyr::count(point_type, doy_band) |>
        tidyr::pivot_wider(names_from = doy_band, values_from = n,
                           values_fill = 0L) |> as.data.frame())

out <- sampled |>
  dplyr::transmute(
    point_id  = seq_len(dplyr::n()),
    domain    = domain_code,
    point_type, doy_band, snow_free_doy, x_utm, y_utm
  )

readr::write_csv(out, "data/derived/year_effect_points.csv")
sf::st_write(
  sf::st_as_sf(out, coords = c("x_utm", "y_utm"),
                crs = 32613, remove = FALSE),
  "data/derived/year_effect_points.gpkg",
  delete_dsn = TRUE, quiet = TRUE
)
cat("\nWrote data/derived/year_effect_points.csv\n")
cat("Wrote data/derived/year_effect_points.gpkg\n")
