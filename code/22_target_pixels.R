# 22_target_pixels.R — build a target meadow-pixel list for AOP extraction.
#
# Strategy: iterate over each 2025 AOP tile (1 km x 1 km blocks listed on
# S3 by filename), and for each tile:
#   1. Read a 1 km windowed view of the LOCAL R3D018 landcover
#   2. Aggregate by factor 3 -> 3 m pixels with meadow-class fraction
#   3. Threshold at >= 0.80 (3x3 80% neighborhood criterion)
#   4. Accumulate the cell centers of qualifying 3 m pixels
#
# Then once:
#   5. Sample snow-free DOY (R4D061) at each candidate
#   6. Stratified random sample by DOY band (early/mid/late)
#   7. Tag each pixel with its parent tile name (for downstream extraction)
#
# Outputs:
#   data/derived/target_pixels.csv     one row per target pixel
#   data/derived/target_pixels.gpkg    same as POINT geometry (EPSG:32613)
#   data/derived/aop_coverage.gpkg     1 polygon per 2025 tile (1 km square)

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(terra)
  library(rSDP)
})

# --- Config -----------------------------------------------------------------
domains            <- c("ALMO", "CRBU", "UPTA")
year               <- 2025
landcover_path     <- "data/raw/SDP/UG_landcover_1m_v4.tif"
meadow_class       <- 3L
# Two-stage filter:
#   1. Within-pixel: each 3m pixel must have >= meadow_pixel_threshold of its
#      9 underlying 1m landcover pixels classified as meadow.
#   2. Neighborhood: each candidate 3m pixel must have >= neighborhood_min
#      of its 9 3x3 neighbors (including itself) also pass filter 1.
# Combined, target pixels are spectrally clean (interior to 9m x 9m meadow).
meadow_pixel_threshold     <- 0.80
neighborhood_min_neighbors <- 6L    # 6 of 9 = >=67%
agg_factor                 <- 3L
doy_bands <- tibble::tribble(
  ~band,    ~min_doy, ~max_doy,
  "early",        0,    130,
  "mid",        130,    155,
  "late",       155,    400
)
n_per_band <- 2000
set.seed(42)

# Point terra temp at a non-root volume in case any op spills to disk.
tmpdir <- file.path("data", "raw", "SDP_tmp")
dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
terra::terraOptions(tempdir = tmpdir)

# --- 1. S3 tile inventory ---------------------------------------------------
cat("Inventorying 2025 AOP tiles on S3 ...\n")
tile_inventory <- purrr::map_dfr(domains, function(dom) {
  raw <- system(sprintf(
    "aws s3 ls --no-sign-request s3://rmbl-chess-data/AOP/spectrometer/mosaic/%s/%d/",
    dom, year), intern = TRUE)
  fns <- stringr::str_extract(raw, "[^ ]+_rfl_\\d+_\\d+\\.nc$")
  fns <- fns[!is.na(fns)]
  if (length(fns) == 0) return(NULL)
  coords <- stringr::str_match(fns,
    "_rfl_(\\d+)_(\\d+)\\.nc$")[, 2:3, drop = FALSE]
  tibble::tibble(
    domain   = dom,
    tile     = fns,
    easting  = as.integer(coords[, 1]),
    northing = as.integer(coords[, 2])
  )
})
cat(sprintf("Total tiles: %d (ALMO=%d CRBU=%d UPTA=%d)\n",
            nrow(tile_inventory),
            sum(tile_inventory$domain == "ALMO"),
            sum(tile_inventory$domain == "CRBU"),
            sum(tile_inventory$domain == "UPTA")))

# Tile polygons: 1 km square; easting/northing = SW corner (NEON / GDAL UTM convention).
tile_polys <- tile_inventory |>
  dplyr::mutate(geometry = purrr::map2(easting, northing, function(e, n) {
    sf::st_polygon(list(matrix(c(
      e, n, e + 1000, n, e + 1000, n + 1000, e, n + 1000, e, n
    ), ncol = 2, byrow = TRUE)))
  })) |>
  sf::st_as_sf(sf_column_name = "geometry", crs = 32613)
sf::st_write(tile_polys, "data/derived/aop_coverage.gpkg",
             delete_dsn = TRUE, quiet = TRUE)
cat("Wrote data/derived/aop_coverage.gpkg\n")

# --- 2-4. Tile-by-tile meadow-pixel extraction ------------------------------
cat(sprintf("Scanning %d tiles for meadow pixels ...\n", nrow(tile_inventory)))
lc <- terra::rast(landcover_path)
t0 <- Sys.time()
candidates_list <- vector("list", nrow(tile_inventory))
report_every <- 100

for (i in seq_len(nrow(tile_inventory))) {
  e <- tile_inventory$easting[i]
  n <- tile_inventory$northing[i]
  win_ext <- terra::ext(e, e + 1000, n, n + 1000)
  # Skip tiles outside the landcover extent
  if (!terra::relate(win_ext, terra::ext(lc), "intersects")) next

  win <- terra::crop(lc, win_ext)
  if (terra::ncell(win) == 0) next

  # 1 m -> 3 m: fraction of class 3 in each 3x3 block (within-pixel purity)
  win_meadow <- (win == meadow_class)
  agg <- terra::aggregate(win_meadow, fact = agg_factor,
                          fun = "mean", na.rm = TRUE)
  pass1 <- terra::classify(agg,
    matrix(c(-Inf, meadow_pixel_threshold, 0,
             meadow_pixel_threshold, Inf, 1), ncol = 3, byrow = TRUE),
    right = FALSE)

  # 3x3 focal: count of pass1 neighbors (incl. self). Neighborhood filter
  # requires >= neighborhood_min_neighbors of 9.
  focal_count <- terra::focal(pass1, w = matrix(1, 3, 3),
                              fun = "sum", na.policy = "all", fillvalue = 0)
  pass2_vals <- terra::values(focal_count, mat = FALSE) >=
                  neighborhood_min_neighbors
  pass1_vals <- terra::values(pass1, mat = FALSE) == 1
  keep_idx <- which(pass1_vals & pass2_vals & !is.na(pass1_vals) &
                                              !is.na(pass2_vals))
  if (length(keep_idx) == 0) next

  xy <- terra::xyFromCell(agg, keep_idx)
  candidates_list[[i]] <- tibble::tibble(
    x = xy[, 1], y = xy[, 2],
    domain = tile_inventory$domain[i],
    tile   = tile_inventory$tile[i]
  )

  if (i %% report_every == 0) {
    cat(sprintf("  %4d / %d tiles  %.1fs elapsed  %s candidates\n",
                i, nrow(tile_inventory),
                as.numeric(Sys.time() - t0, units = "secs"),
                format(sum(purrr::map_int(candidates_list, NROW)),
                       big.mark = ",")))
  }
}
candidates <- dplyr::bind_rows(candidates_list)
cat(sprintf("Tile pass done: %.1fs   total candidates: %s\n",
            as.numeric(Sys.time() - t0, units = "secs"),
            format(nrow(candidates), big.mark = ",")))

# Pre-sample: 37M+ candidates is too many for the DOY extract step. Sample
# uniformly to a manageable size first; we'll stratify by DOY after the
# extract. The pre-sample is uniform across tiles so it doesn't favor large
# domains too much.
prefilter_n <- 200000L
if (nrow(candidates) > prefilter_n) {
  candidates <- candidates |>
    dplyr::slice_sample(n = prefilter_n, replace = FALSE)
  cat(sprintf("Pre-sampled to %s candidates for the DOY extract step.\n",
              format(nrow(candidates), big.mark = ",")))
}

# --- 5. Snow-free DOY at each candidate -------------------------------------
cat("Extracting snow-free DOY (R4D061) at candidates ...\n")
sf_doy <- rSDP::sdp_get_raster("R4D061")
if (terra::nlyr(sf_doy) > 1) sf_doy <- sf_doy[[1]]
candidates$snow_free_doy <- terra::extract(
  sf_doy,
  terra::vect(candidates, geom = c("x", "y"), crs = "EPSG:32613"),
  ID = FALSE
)[[1]]
candidates <- candidates |>
  dplyr::filter(!is.na(snow_free_doy)) |>
  dplyr::mutate(doy_band = cut(snow_free_doy,
    breaks = c(doy_bands$min_doy, max(doy_bands$max_doy)),
    labels = doy_bands$band, include.lowest = TRUE))

cat("DOY band counts among candidates:\n")
print(candidates |> dplyr::count(doy_band, domain) |>
  tidyr::pivot_wider(names_from = domain, values_from = n, values_fill = 0L))

# --- 6. Stratified random sample -------------------------------------------
cat(sprintf("Sampling %d per DOY band (or all if fewer) ...\n", n_per_band))
sampled <- candidates |>
  dplyr::group_split(doy_band) |>
  purrr::map_dfr(function(df) {
    n <- min(n_per_band, nrow(df))
    dplyr::slice_sample(df, n = n, replace = FALSE)
  })
cat(sprintf("Final target pixels: %s\n", format(nrow(sampled), big.mark = ",")))
print(sampled |> dplyr::count(doy_band, domain) |>
  tidyr::pivot_wider(names_from = domain, values_from = n, values_fill = 0L))

# --- 7. Save -----------------------------------------------------------------
out_csv <- sampled |>
  dplyr::transmute(x_utm = x, y_utm = y, snow_free_doy, doy_band, domain, tile)
readr::write_csv(out_csv, "data/derived/target_pixels.csv")
cat("Wrote data/derived/target_pixels.csv\n")

out_sf <- sf::st_as_sf(out_csv, coords = c("x_utm", "y_utm"),
                       crs = 32613, remove = FALSE)
sf::st_write(out_sf, "data/derived/target_pixels.gpkg",
             delete_dsn = TRUE, quiet = TRUE)
cat("Wrote data/derived/target_pixels.gpkg\n")

unlink(tmpdir, recursive = TRUE)
