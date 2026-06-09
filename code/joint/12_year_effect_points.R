# 12_year_effect_points.R — stratified random point set in the CRBU domain
# for evaluating 2018-vs-2025 AOP year effects via direct re-extraction
# at the same locations.
#
# Strata (point_type = coarse veg/non-veg; cover_class = lifeform):
#   meadow_shrub    R3D018 class 3 (meadow) + class 10 (deciduous shrub)
#                   [point_type vegetated]. Phenology differences between
#                   years (NDVI, red edge, NIR plateau, etc.).
#   tree_deciduous  R3D018 class 2 (deciduous trees > 2 m) [vegetated].
#   tree_evergreen  R3D018 class 1 (evergreen trees & shrubs) [vegetated].
#                   The two tree strata let us break the vegetated year-
#                   effect comparison out by canopy lifeform.
#   bare            R3D018 class 6 (rock / bare soil / sparse veg)
#                   [point_type non_vegetated]. Surfaces that don't grow
#                   back season-to-season — any spectral shift between
#                   2018 and 2025 is most likely instrument / radiometric
#                   calibration.
# All strata require >= 80% purity in the 3 m cell.
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

# One binary 3 m purity mask per landcover stratum (>= purity_threshold of
# the 9 underlying 1 m cells are the target class), snapped to the CHM grid.
purity_mask <- function(target_classes) {
  rcl <- cbind(target_classes, target_classes, 1)
  m1  <- terra::classify(lc_dom, rcl, others = 0L,
                         include.lowest = TRUE, right = NA, datatype = "INT1U")
  frac <- terra::aggregate(m1, fact = 3, fun = "mean", na.rm = TRUE)
  terra::resample(frac >= purity_threshold, ref_3m, method = "near")
}

# point_type is the coarse veg/non-veg split scripts 13/15 rely on;
# cover_class is the lifeform breakout this set adds.
strata <- tibble::tribble(
  ~cover_class,     ~classes,   ~point_type,
  "meadow_shrub",   c(3L, 10L), "vegetated",
  "tree_deciduous", 2L,         "vegetated",
  "tree_evergreen", 1L,         "vegetated",
  "bare",           6L,         "non_vegetated"
)

# --- 2. Candidate cell centers per stratum ---------------------------------
candidates <- function(mask, cover_label, type_label) {
  idx <- which(terra::values(mask, mat = FALSE) == 1)
  if (length(idx) == 0) return(NULL)
  xy <- terra::xyFromCell(mask, idx)
  tibble::tibble(x_utm = xy[, 1], y_utm = xy[, 2],
                 cover_class = cover_label, point_type = type_label)
}
cand_list <- purrr::pmap(strata, function(cover_class, classes, point_type) {
  candidates(purity_mask(classes), cover_class, point_type)
})
cat("Candidates per cover_class:\n")
for (d in cand_list) if (!is.null(d))
  cat(sprintf("  %-15s %s\n", d$cover_class[1],
              format(nrow(d), big.mark = ",")))

# Pre-sample each stratum to manageable size before the DOY extract.
prefilter_n <- 50000L
prefilter <- function(df, n) {
  if (is.null(df) || nrow(df) == 0) return(df)
  if (nrow(df) > n) dplyr::slice_sample(df, n = n, replace = FALSE) else df
}
all_pts <- purrr::map_dfr(cand_list, prefilter, n = prefilter_n)

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

cat("Candidate counts by cover_class x doy_band:\n")
print(all_pts |> dplyr::count(cover_class, doy_band) |>
        tidyr::pivot_wider(names_from = doy_band, values_from = n,
                           values_fill = 0L) |> as.data.frame())

# --- 4. Stratified random sample (n_per_stratum per cover_class x DOY) ------
sampled <- all_pts |>
  dplyr::group_by(cover_class, doy_band) |>
  dplyr::group_split() |>
  purrr::map_dfr(function(df) {
    dplyr::slice_sample(df, n = min(n_per_stratum, nrow(df)),
                         replace = FALSE)
  })

cat(sprintf("\nFinal: %d points\n", nrow(sampled)))
print(sampled |> dplyr::count(cover_class, doy_band) |>
        tidyr::pivot_wider(names_from = doy_band, values_from = n,
                           values_fill = 0L) |> as.data.frame())

out <- sampled |>
  dplyr::transmute(
    point_id  = seq_len(dplyr::n()),
    domain    = domain_code,
    point_type, cover_class, doy_band, snow_free_doy, x_utm, y_utm
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
