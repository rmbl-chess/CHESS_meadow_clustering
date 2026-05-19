# 36_canopy_height.R — extract NEON 1 m canopy height model (CHM) at every
# meadow and shrub crown centroid. Output is a (site_number, Year) table
# with canopy height in meters. Used as a covariate in joint meadow+shrub
# training (script 37).
#
# Implementation: per-domain SINGLE-PASS point extract at crown centroids
# (no per-polygon buffer iteration). Two earlier passes were killed
# because terra was triggering an HTTP windowed read for each polygon —
# slow against the Google Drive-mounted CHM tifs. Point extraction is one
# pixel read per point and ~100x faster.
#
# Canopy-variability info (p90 / max within a small buffer) is dropped;
# the CHM is 1 m and the AOP-feature recipe uses a 3x3 mean anyway, so a
# single-pixel centroid value tracks the same thing the classifier sees.
#
# 2018 plots are all in CRBU and reuse the 2025 CRBU CHM (canopy structure
# is stable over 7 years for established communities).
#
# Inputs:
#   data/derived/crowns_2018.gpkg, .../crowns_2025.gpkg
#   /Users/ian/Library/CloudStorage/.../CHESS25_{ALMO,CRBU,UPTA}_CHM_1m_v*.tif
# Outputs:
#   data/derived/canopy_height.rds
#     site_number, Year, n_crowns, canopy_height_m

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(terra)
})

terra::terraOptions(progress = 0)   # progress bar slows scattered I/O

chm_dir <- "/Users/ian/Library/CloudStorage/GoogleDrive-ibreckhe@gmail.com/My Drive/BreckheimerLab2025/Projects/CHESS/Data/AOP_mosaics/NEON_delivered"
chm_paths <- list(
  ALMO = file.path(chm_dir, "CHESS25_ALMO_CHM_1m_v1.tif"),
  CRBU = file.path(chm_dir, "CHESS25_CRBU_CHM_1m_v2.tif"),
  UPTA = file.path(chm_dir, "CHESS25_UPTA_CHM_1m_v1.tif")
)
stopifnot(all(file.exists(unlist(chm_paths))))

# --- Centroids by domain --------------------------------------------------
crowns_2018 <- sf::st_read("data/derived/crowns_2018.gpkg", quiet = TRUE) |>
  dplyr::mutate(Year = 2018L, domain = "CRBU") |>
  dplyr::select(site_number, Year, domain)
crowns_2025 <- sf::st_read("data/derived/crowns_2025.gpkg", quiet = TRUE) |>
  dplyr::mutate(Year = 2025L) |>
  dplyr::select(site_number, Year, domain)

centroids <- dplyr::bind_rows(crowns_2018, crowns_2025) |>
  sf::st_transform(32613) |>
  sf::st_centroid()
cat(sprintf("Crown centroids: %d (2018=%d, 2025=%d)\n",
            nrow(centroids), sum(centroids$Year == 2018L),
            sum(centroids$Year == 2025L)))

# --- Per-domain bulk point extract ----------------------------------------
extract_chm_for_domain <- function(centroids_dom, chm_path) {
  cat(sprintf("%s: opening CHM ... ", basename(chm_path)))
  chm  <- terra::rast(chm_path)
  vect <- terra::vect(centroids_dom)
  cat(sprintf("extracting at %d centroids ... ", length(vect)))
  t0 <- Sys.time()
  vals <- terra::extract(chm, vect, ID = FALSE)
  cat(sprintf("done (%.1fs)\n",
              as.numeric(Sys.time() - t0, units = "secs")))
  centroids_dom |>
    sf::st_drop_geometry() |>
    dplyr::mutate(canopy_height_m = vals[[1]])
}

per_crown <- purrr::map_dfr(names(chm_paths), function(dom) {
  sub <- centroids |>
    dplyr::filter(domain == dom |
                  (Year == 2018L & dom == "CRBU"))
  if (nrow(sub) == 0) return(NULL)
  extract_chm_for_domain(sub, chm_paths[[dom]])
})

# --- Aggregate to per-site -----------------------------------------------
chm_per_site <- per_crown |>
  dplyr::filter(!is.na(canopy_height_m)) |>
  dplyr::group_by(site_number, Year) |>
  dplyr::summarise(
    n_crowns        = dplyr::n(),
    canopy_height_m = mean(canopy_height_m, na.rm = TRUE),
    .groups         = "drop"
  )

cat(sprintf("\nPer-site CHM table: %d sites (2018=%d, 2025=%d)\n",
            nrow(chm_per_site),
            sum(chm_per_site$Year == 2018L),
            sum(chm_per_site$Year == 2025L)))
cat("\nCanopy height summary (m):\n")
print(summary(chm_per_site$canopy_height_m))

saveRDS(chm_per_site, "data/derived/canopy_height.rds")
cat(sprintf("\nWrote data/derived/canopy_height.rds (%d sites)\n",
            nrow(chm_per_site)))
