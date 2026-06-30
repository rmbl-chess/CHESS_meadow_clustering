# 01_load.R — read 2018 and 2025 field vegetation + extracted-spectra inputs.
#
# Inputs (see data/README.md):
#   data/raw/ESS-DIVE-Vegetation-Field-2018/
#       fractional_cover (1).csv      SampleSiteCode, SamplingArea, CollectionDate,
#                                      CoverCode, FractionalCover, Note
#       species_list (1).csv          CoverCode, Family, Genus, Species,
#                                      AltFieldCode, Notes
#       metadata_column_key (1).csv
#       CRBU2018_AOP_Crowns.geojson   id, SiteCode  (EPSG:32613)
#   data/raw/ESS-DIVE-Vegetation-Field-2025/
#       chess_meadow_cover_cleaned.csv   Site_Number, Sampling_Area, Notes,
#                                         Cover_Type, Cover_Class_Name, Cover_Percent
#       chess_meadow_site_cleaned.csv    Location_Type, Sampling_Area, Site_Number,
#                                         Collection_Date, ...
#       chess_species_list_cleaned.csv   Taxon_family, Taxon_binomial, Taxon_label,
#                                         Taxon_full_name, GBIF_Taxon_ID, ...
#   data/raw/ESS-DIVE-Spectra/
#       site_extraction_spectra_2018 (1).csv  site_number, fid, row, col, x_utm,
#                                              y_utm, shade, rfl_band_1..426,
#                                              unc_band_1..426   (859 cols)
#       site_extraction_spectra_2025 (1).csv  site_number, domain, sampling_area,
#                                              site_type, fid, row, col, x_utm,
#                                              y_utm, shade, rfl_band_1..426,
#                                              unc_band_1..426   (862 cols)
#       wavelengths_{2018,2025}.csv   band_number, center_wavelength, fwhm
#       CHESS_2025_crowns (1).geojson  domain, sampling_area, site_type,
#                                       site_number, ...        (EPSG:32613)
#
# Outputs (data/derived/):
#   veg_2018.rds, veg_2025.rds            cover + species list per year
#   crowns_2018.gpkg, crowns_2025.gpkg    standardized to EPSG:32613
#   spectra_2018.rds, spectra_2025.rds    spectral matrix + wavelengths per year

library(tidyverse)
library(sf)

raw_2018 <- "data/raw/ESS-DIVE-Vegetation-Field-2018"
raw_2025 <- "data/raw/ESS-DIVE-Vegetation-Field-2025"
raw_2026 <- "data/raw/Supplemental_field_2026"
raw_spec <- "data/raw/ESS-DIVE-Spectra"

dir.create("data/derived", showWarnings = FALSE, recursive = TRUE)

# --- 2018 vegetation -------------------------------------------------------
cover_2018   <- readr::read_csv(file.path(raw_2018, "fractional_cover (1).csv"),
                                show_col_types = FALSE)
species_2018 <- readr::read_csv(file.path(raw_2018, "species_list (1).csv"),
                                show_col_types = FALSE)
colkey_2018  <- readr::read_csv(file.path(raw_2018, "metadata_column_key (1).csv"),
                                show_col_types = FALSE)

# 2018 SampleSiteCode is e.g. "325-ER18"; the leading integer matches
# crowns_2018$id and spectra_2018$site_number. Make that explicit here.
cover_2018 <- cover_2018 |>
  dplyr::mutate(site_number = as.integer(stringr::str_extract(SampleSiteCode, "^\\d+")))

# --- 2025 vegetation -------------------------------------------------------
cover_2025   <- readr::read_csv(file.path(raw_2025, "chess_meadow_cover_cleaned.csv"),
                                show_col_types = FALSE)
sites_2025   <- readr::read_csv(file.path(raw_2025, "chess_meadow_site_cleaned.csv"),
                                show_col_types = FALSE)
species_2025 <- readr::read_csv(file.path(raw_2025, "chess_species_list_cleaned.csv"),
                                show_col_types = FALSE)

cover_2025 <- cover_2025 |>
  dplyr::rename(site_number = Site_Number) |>
  dplyr::mutate(site_number = as.integer(site_number))

# --- 2026 supplemental vegetation ------------------------------------------
# Same cover schema as 2025 (Cover_Type / Cover_Class_Name / Cover_Percent);
# no separate species_list (Cover_Class_Name is already the binomial). Spectra
# are extracted from 2025 AOP by code/python/extract_supplemental_spectra.py.
cover_2026 <- readr::read_csv(
  file.path(raw_2026, "augment_cover_cleaned_2026_06_29.csv"),
  show_col_types = FALSE) |>
  dplyr::rename(site_number = Site_Number) |>
  dplyr::mutate(site_number = as.integer(site_number))

# --- Crown polygons --------------------------------------------------------
# Both GeoJSONs are already EPSG:32613; transform is defensive (no-op when matched).
crowns_2018 <- sf::st_read(file.path(raw_2018, "CRBU2018_AOP_Crowns.geojson"),
                           quiet = TRUE) |>
  sf::st_transform(32613) |>
  dplyr::rename(site_number = id, sample_site_code = SiteCode)

crowns_2025 <- sf::st_read(file.path(raw_spec, "CHESS_2025_crowns (1).geojson"),
                           quiet = TRUE) |>
  sf::st_transform(32613) |>
  dplyr::mutate(site_number = as.integer(site_number))

# 2026 crowns: only carry a site label (Label_of_F). domain + site_type are
# backfilled below from the extracted spectra (the python extractor is the
# authoritative source for both); force to MULTIPOLYGON to match 2025.
crowns_2026 <- sf::st_read(
  file.path(raw_2026, "augment_polygons_2026_06_29_wgs_utm.geojson"),
  quiet = TRUE) |>
  sf::st_zm(drop = TRUE) |>
  sf::st_transform(32613) |>
  dplyr::transmute(site_number = as.integer(Label_of_F))

# --- Spectra (pre-extracted at crown footprints) ---------------------------
spectra_2018 <- readr::read_csv(file.path(raw_spec, "site_extraction_spectra_2018 (1).csv"),
                                show_col_types = FALSE) |>
  dplyr::mutate(site_number = as.integer(site_number))
spectra_2025 <- readr::read_csv(file.path(raw_spec, "site_extraction_spectra_2025 (1).csv"),
                                show_col_types = FALSE) |>
  dplyr::mutate(site_number = as.integer(site_number))
wl_2018      <- readr::read_csv(file.path(raw_spec, "wavelengths_2018.csv"),
                                show_col_types = FALSE)
wl_2025      <- readr::read_csv(file.path(raw_spec, "wavelengths_2025.csv"),
                                show_col_types = FALSE)

# 2026 spectra (extracted from 2025 AOP on the Hub). Same schema + wavelength
# grid as 2025. Guarded so the rest of the pipeline runs before the Hub
# extraction lands; domain/site_type from the extractor are authoritative.
spectra_2026_path <- file.path(raw_spec, "site_extraction_spectra_2026_06_29.csv")
have_2026_spectra <- file.exists(spectra_2026_path)
if (have_2026_spectra) {
  spectra_2026 <- readr::read_csv(spectra_2026_path, show_col_types = FALSE) |>
    dplyr::mutate(site_number = as.integer(site_number))
  # Backfill 2026 crowns' domain + site_type from the spectra (one row per site).
  site_meta_2026 <- spectra_2026 |>
    dplyr::distinct(site_number, domain, site_type)
  crowns_2026 <- crowns_2026 |>
    dplyr::left_join(site_meta_2026, by = "site_number")
} else {
  message("No 2026 spectra at ", spectra_2026_path,
          " — run code/python/extract_supplemental_spectra.py on the Hub first. ",
          "Skipping spectra_2026.rds; 2026 sites will not enter clustering yet.")
}

# --- Persist ---------------------------------------------------------------
saveRDS(list(cover = cover_2018, species = species_2018, colkey = colkey_2018),
        "data/derived/veg_2018.rds")
saveRDS(list(cover = cover_2025, species = species_2025, sites = sites_2025),
        "data/derived/veg_2025.rds")
# 2026 has no species_list; Cover_Class_Name is already the binomial.
saveRDS(list(cover = cover_2026, species = NULL, sites = NULL),
        "data/derived/veg_2026.rds")
sf::st_write(crowns_2018, "data/derived/crowns_2018.gpkg",
             delete_dsn = TRUE, quiet = TRUE)
sf::st_write(crowns_2025, "data/derived/crowns_2025.gpkg",
             delete_dsn = TRUE, quiet = TRUE)
sf::st_write(crowns_2026, "data/derived/crowns_2026.gpkg",
             delete_dsn = TRUE, quiet = TRUE)
saveRDS(list(spectra = spectra_2018, wavelengths = wl_2018),
        "data/derived/spectra_2018.rds")
saveRDS(list(spectra = spectra_2025, wavelengths = wl_2025),
        "data/derived/spectra_2025.rds")
# 2026 spectra share the 2025 wavelength grid (extracted from 2025 AOP).
if (have_2026_spectra) {
  saveRDS(list(spectra = spectra_2026, wavelengths = wl_2025),
          "data/derived/spectra_2026.rds")
}
