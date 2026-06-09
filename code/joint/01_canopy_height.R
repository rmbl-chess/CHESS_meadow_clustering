# 01_canopy_height.R — extract NEON canopy height at every meadow + shrub
# crown centroid, using a 3 m maximum-statistic aggregation of the 1 m CHM.
#
# Step A: aggregate each 1 m domain CHM to 3 m by max (one-time pre-pass;
#         the 3 m COG is written under data/derived/aop_chm_3m/ and reused
#         by the inference script 09_inference.R, both as the
#         canopy_height_m feature AND as the tree mask).
# Step B: point-extract the 3 m max value at every crown centroid for the
#         training table.
#
# Using max (rather than mean) is what makes the same raster usable for
# tree masking: a 3 m cell containing any tree pixel registers as tall,
# even if most of the cell is grass.
#
# 2018 plots are all in CRBU and reuse the 2025 CRBU CHM (canopy structure
# is stable over 7 years for established communities).
#
# Inputs:
#   data/derived/crowns_2018.gpkg, .../crowns_2025.gpkg
#   /Users/ian/Library/CloudStorage/.../CHESS25_{ALMO,CRBU,UPTA}_CHM_1m_v*.tif
# Outputs:
#   data/derived/aop_chm_3m/{ALMO,CRBU,UPTA}_chm_max_3m.tif   (COGs)
#   data/derived/canopy_height.rds                            (per-site)

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(terra)
})

terra::terraOptions(progress = 0)   # progress bar slows scattered I/O

chm_dir <- "/Users/ian/Library/CloudStorage/GoogleDrive-ibreckhe@gmail.com/My Drive/BreckheimerLab2025/Projects/CHESS/Data/AOP_mosaics/NEON_delivered"
chm_paths_1m <- list(
  ALMO = file.path(chm_dir, "CHESS25_ALMO_CHM_1m_v1.tif"),
  CRBU = file.path(chm_dir, "CHESS25_CRBU_CHM_1m_v2.tif"),
  UPTA = file.path(chm_dir, "CHESS25_UPTA_CHM_1m_v1.tif")
)
stopifnot(all(file.exists(unlist(chm_paths_1m))))

agg_dir <- "data/derived/aop_chm_3m"
dir.create(agg_dir, showWarnings = FALSE, recursive = TRUE)
chm_paths_3m <- setNames(
  file.path(agg_dir, sprintf("%s_chm_max_3m.tif", names(chm_paths_1m))),
  names(chm_paths_1m)
)

# --- Step A: 1m -> 3m max aggregation (one-time per domain) --------------
for (dom in names(chm_paths_1m)) {
  out <- chm_paths_3m[[dom]]
  if (file.exists(out)) {
    cat(sprintf("%s: 3m max CHM already present (%s)\n", dom, basename(out)))
    next
  }
  cat(sprintf("%s: aggregating 1m -> 3m max ... ", dom))
  t0 <- Sys.time()
  chm_1m <- terra::rast(chm_paths_1m[[dom]])
  terra::aggregate(
    chm_1m, fact = 3, fun = "max", na.rm = TRUE,
    filename = out,
    overwrite = FALSE,
    wopt = list(
      filetype = "GTiff",
      gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2",
               "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256",
               "BIGTIFF=IF_SAFER")
    )
  )
  cat(sprintf("done (%.1fs)\n",
              as.numeric(Sys.time() - t0, units = "secs")))
}

# --- Step B: centroids by domain + point extract from 3m rasters ---------
crowns_2018 <- sf::st_read("data/derived/crowns_2018.gpkg", quiet = TRUE) |>
  dplyr::mutate(Year = 2018L, domain = "CRBU") |>
  dplyr::select(site_number, Year, domain)
crowns_2025 <- sf::st_read("data/derived/crowns_2025.gpkg", quiet = TRUE) |>
  dplyr::mutate(Year = 2025L) |>
  dplyr::select(site_number, Year, domain)

# 2026 supplemental crowns (ALMO + CRBU); domain backfilled in 01_load.R from
# the extracted spectra. Reuses the 2025 CHM for its domain (canopy structure),
# same convention as 2018 reusing the CRBU CHM.
crowns_2026_path <- "data/derived/crowns_2026.gpkg"
crowns_2026 <- if (file.exists(crowns_2026_path)) {
  c26 <- sf::st_read(crowns_2026_path, quiet = TRUE)
  if ("domain" %in% names(c26)) {
    c26 |> dplyr::mutate(Year = 2026L) |>
      dplyr::select(site_number, Year, domain)
  } else NULL
} else NULL

centroids <- dplyr::bind_rows(crowns_2018, crowns_2025, crowns_2026) |>
  sf::st_transform(32613) |>
  sf::st_centroid()
cat(sprintf("\nCrown centroids: %d (2018=%d, 2025=%d, 2026=%d)\n",
            nrow(centroids), sum(centroids$Year == 2018L),
            sum(centroids$Year == 2025L), sum(centroids$Year == 2026L)))

extract_chm_for_domain <- function(centroids_dom, chm_path_3m) {
  cat(sprintf("%s: extracting 3m max CHM at %d centroids ... ",
              basename(chm_path_3m), nrow(centroids_dom)))
  t0   <- Sys.time()
  chm  <- terra::rast(chm_path_3m)
  vect <- terra::vect(centroids_dom)
  vals <- terra::extract(chm, vect, ID = FALSE)
  cat(sprintf("done (%.1fs)\n",
              as.numeric(Sys.time() - t0, units = "secs")))
  centroids_dom |>
    sf::st_drop_geometry() |>
    dplyr::mutate(canopy_height_m = vals[[1]])
}

per_crown <- purrr::map_dfr(names(chm_paths_3m), function(dom) {
  sub <- centroids |>
    dplyr::filter(domain == dom |
                  (Year == 2018L & dom == "CRBU"))
  if (nrow(sub) == 0) return(NULL)
  extract_chm_for_domain(sub, chm_paths_3m[[dom]])
})

# --- Aggregate to per-site (a few sites have multiple crowns) ------------
chm_per_site <- per_crown |>
  dplyr::filter(!is.na(canopy_height_m)) |>
  dplyr::group_by(site_number, Year) |>
  dplyr::summarise(
    n_crowns        = dplyr::n(),
    canopy_height_m = mean(canopy_height_m, na.rm = TRUE),
    .groups         = "drop"
  )

cat(sprintf("\nPer-site CHM table: %d sites (2018=%d, 2025=%d, 2026=%d)\n",
            nrow(chm_per_site),
            sum(chm_per_site$Year == 2018L),
            sum(chm_per_site$Year == 2025L),
            sum(chm_per_site$Year == 2026L)))
cat("\n3m max canopy height summary (m):\n")
print(summary(chm_per_site$canopy_height_m))

saveRDS(chm_per_site, "data/derived/canopy_height.rds")
cat(sprintf("\nWrote data/derived/canopy_height.rds (%d sites)\n",
            nrow(chm_per_site)))
