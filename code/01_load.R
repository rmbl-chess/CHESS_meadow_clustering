# 01_load.R — read 2018 and 2025 field vegetation + extracted-spectra inputs.
#
# Inputs (see data/README.md):
#   data/raw/ESS-DIVE-Vegetation-Field-2018/
#       fractional_cover (1).csv
#       species_list (1).csv
#       metadata_column_key (1).csv
#       CRBU2018_AOP_Crowns.geojson
#   data/raw/ESS-DIVE-Vegetation-Field-2025/
#       chess_meadow_cover_cleaned.csv
#       chess_meadow_site_cleaned.csv
#       chess_species_list_cleaned.csv
#   data/raw/ESS-DIVE-Spectra/
#       site_extraction_spectra_2018 (1).csv
#       site_extraction_spectra_2025 (1).csv
#       wavelengths_2018.csv
#       wavelengths_2025.csv
#       CHESS_2025_crowns (1).geojson
#
# Outputs (data/derived/):
#   veg_2018.rds, veg_2025.rds          # cover + species list per year
#   crowns_2018.gpkg, crowns_2025.gpkg  # reprojected to EPSG:32613
#   spectra_2018.rds, spectra_2025.rds  # spectral matrix + wavelengths per year

library(tidyverse)
library(sf)

raw_2018 <- "data/raw/ESS-DIVE-Vegetation-Field-2018"
raw_2025 <- "data/raw/ESS-DIVE-Vegetation-Field-2025"
raw_spec <- "data/raw/ESS-DIVE-Spectra"

dir.create("data/derived", showWarnings = FALSE, recursive = TRUE)

# --- 2018 vegetation -------------------------------------------------------
cover_2018   <- readr::read_csv(file.path(raw_2018, "fractional_cover (1).csv"))
species_2018 <- readr::read_csv(file.path(raw_2018, "species_list (1).csv"))
colkey_2018  <- readr::read_csv(file.path(raw_2018, "metadata_column_key (1).csv"))

# --- 2025 vegetation -------------------------------------------------------
cover_2025   <- readr::read_csv(file.path(raw_2025, "chess_meadow_cover_cleaned.csv"))
sites_2025   <- readr::read_csv(file.path(raw_2025, "chess_meadow_site_cleaned.csv"))
species_2025 <- readr::read_csv(file.path(raw_2025, "chess_species_list_cleaned.csv"))

# --- Crown polygons (spatial join unit) ------------------------------------
# Reproject everything to EPSG:32613 up front so downstream joins assume one CRS.
crowns_2018 <- sf::st_read(file.path(raw_2018, "CRBU2018_AOP_Crowns.geojson"), quiet = TRUE) |>
  sf::st_transform(32613)
crowns_2025 <- sf::st_read(file.path(raw_spec, "CHESS_2025_crowns (1).geojson"), quiet = TRUE) |>
  sf::st_transform(32613)

# --- Spectra (pre-extracted at crown footprints) ---------------------------
spectra_2018 <- readr::read_csv(file.path(raw_spec, "site_extraction_spectra_2018 (1).csv"))
spectra_2025 <- readr::read_csv(file.path(raw_spec, "site_extraction_spectra_2025 (1).csv"))
wl_2018      <- readr::read_csv(file.path(raw_spec, "wavelengths_2018.csv"))
wl_2025      <- readr::read_csv(file.path(raw_spec, "wavelengths_2025.csv"))

# --- Persist ---------------------------------------------------------------
saveRDS(list(cover = cover_2018, species = species_2018, colkey = colkey_2018),
        "data/derived/veg_2018.rds")
saveRDS(list(cover = cover_2025, species = species_2025, sites = sites_2025),
        "data/derived/veg_2025.rds")
sf::st_write(crowns_2018, "data/derived/crowns_2018.gpkg", delete_dsn = TRUE, quiet = TRUE)
sf::st_write(crowns_2025, "data/derived/crowns_2025.gpkg", delete_dsn = TRUE, quiet = TRUE)
saveRDS(list(spectra = spectra_2018, wavelengths = wl_2018),
        "data/derived/spectra_2018.rds")
saveRDS(list(spectra = spectra_2025, wavelengths = wl_2025),
        "data/derived/spectra_2025.rds")
