# 12_extract_environment.R — extract snow-free DOY climatology from the SDP
# at each site's crown centroid. Snow-free date integrates elevation,
# aspect, and microsite drainage; in RMBL it's the single covariate that
# explains the most ecological variation at meadow plot scale. 27m raster
# (R4D061) is much smaller than 1-3m topographic layers, so extraction is
# fast.
#
#   R4D061  Snowpack Persistence DOY Mean (1993-2022)  -> snow_free_doy
#
# (If diagnosis suggests we need more env covariates, add them here.)
#
# Sites with multiple crown polygons (a few 2018 sites) get the mean of
# their crown centroids' env values.
#
# Inputs:  data/derived/crowns_2018.gpkg, .../crowns_2025.gpkg
# Outputs: data/derived/environment.rds  (site_number, Year, snow_free_doy)

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(terra)
  library(rSDP)
})

# Load crowns and reduce to (site_number, Year) centroids in EPSG:32613.
crowns_2018 <- sf::st_read("data/derived/crowns_2018.gpkg", quiet = TRUE) |>
  dplyr::mutate(Year = 2018L) |>
  dplyr::select(site_number, Year)
crowns_2025 <- sf::st_read("data/derived/crowns_2025.gpkg", quiet = TRUE) |>
  dplyr::mutate(Year = 2025L) |>
  dplyr::select(site_number, Year)

crowns <- dplyr::bind_rows(crowns_2018, crowns_2025) |>
  sf::st_centroid() |>
  sf::st_transform(32613)

cat(sprintf("Total crown centroids: %d (2018=%d, 2025=%d)\n",
            nrow(crowns),
            sum(crowns$Year == 2018L), sum(crowns$Year == 2025L)))

# Convert to terra SpatVector for raster extraction.
crowns_vect <- terra::vect(crowns)

# --- Pull raster and extract values ----------------------------------------
cat("Loading R4D061 (snow-free DOY) ... ")
t0 <- Sys.time()
r  <- rSDP::sdp_get_raster("R4D061")
if (terra::nlyr(r) > 1) r <- r[[1]]
cat(sprintf("done (%.1fs)\n", as.numeric(Sys.time() - t0, units = "secs")))

cat("Extracting at crown centroids ... ")
t1 <- Sys.time()
pts  <- terra::project(crowns_vect, terra::crs(r))
vals <- terra::extract(r, pts, ID = FALSE)[[1]]
cat(sprintf("done (%.1fs)\n", as.numeric(Sys.time() - t1, units = "secs")))

env_df <- dplyr::bind_cols(
  crowns |> sf::st_drop_geometry() |> dplyr::select(site_number, Year),
  tibble::tibble(snow_free_doy = vals)
)

# Average across crown centroids per (site_number, Year).
env_site <- env_df |>
  dplyr::group_by(site_number, Year) |>
  dplyr::summarise(dplyr::across(dplyr::everything(),
                                 ~ mean(.x, na.rm = TRUE)),
                   .groups = "drop")

cat("\nEnvironmental summary (per-site means):\n")
print(env_site |> dplyr::summarise(dplyr::across(
  -c(site_number, Year),
  list(min = ~ min(.x, na.rm = TRUE),
       med = ~ stats::median(.x, na.rm = TRUE),
       max = ~ max(.x, na.rm = TRUE)),
  .names = "{.col}__{.fn}")) |> tidyr::pivot_longer(everything()) |>
  tidyr::separate(name, c("var", "stat"), sep = "__"))

saveRDS(env_site, "data/derived/environment.rds")
cat(sprintf("\nWrote data/derived/environment.rds (%d sites).\n",
            nrow(env_site)))
