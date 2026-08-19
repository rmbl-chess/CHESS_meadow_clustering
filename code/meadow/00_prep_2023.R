# 00_prep_2023.R — standardize the 2023 UER vegmap campaign into the shapes
# the rest of the pipeline already consumes (the 2026 supplemental campaign is
# the template: Cover_Class_Name IS the canonical binomial, no species list).
#
# The 2023 campaign (238 ~1 m2 plots, all CRBU) arrives with its own schema:
#   vegmap_2023_cover.csv                    Date, Site ("BV-1"), Species
#                                            (6-letter code), Cover (%)
#   vegmap_2023_wgs_utm_pixel_select.geojson AOP-pixel-snapped plot polygons
#                                            (EPSG:32613; Name = site label;
#                                            pix_sel = reviewed flag)
#   data/small_reference/taxonomy_crosswalk_2023.csv
#                                            collaborator's 3-campaign code ->
#                                            canonical crosswalk (2023 rows used)
#
# This script:
#   1. Normalizes site names on both sides (LT4 -> LT-4, SN-01 -> SN-1,
#      uu-1 -> UU-1) and reports cover<->polygon mismatches (known open issue:
#      SH-1..4 have cover but no polygon; JF-13..17 the reverse).
#   2. Assigns stable integer site_numbers (3001+) via a committed lookup
#      (data/small_reference/site_numbers_2023.csv) — existing assignments are
#      reused, new names appended, so numbers never shift across re-runs.
#   3. Maps species codes -> canonical binomials (two overrides below align
#      2023 with project conventions) and morphotype codes -> the 2025-style
#      non-species tokens, emitting a 2026-schema cover table.
#
# Outputs (data/derived/):
#   vegmap_2023_cover_std.csv      Site_Number, Site_Name, Collection_Date,
#                                  Cover_Type, Cover_Class_Name, Cover_Percent
#   vegmap_2023_polygons_std.geojson  site_number, site_name, geometry
#                                  (EPSG:32613; input to 01_load.R AND to
#                                  extract_supplemental_spectra.py
#                                  --site-field site_number on the Hub)
# Plus the (committed) site-number lookup in data/small_reference/.

library(tidyverse)
library(sf)

raw_dir <- "data/raw/Vegmap_2023"
dir.create("data/derived", showWarnings = FALSE, recursive = TRUE)

# --- Canonical-name overrides ----------------------------------------------
# GEUROS: collaborator crosswalk has misspelled "Geum rossi", and the 2018
#   campaign's canonical for the same taxon is Acomastylis rossii.
# VERCAL: project convention treats all Veratrum as V. tenuipetalum
#   (V. californicum is a synonym; see taxonomy_crosswalk_manual_edits.csv).
canonical_overrides <- c(
  "GEUROS" = "Acomastylis rossii",
  "VERCAL" = "Veratrum tenuipetalum"
)

# Morphotype codes -> the 2025/2026 non-species tokens, so the campaign-2023
# rows in nonspecies_category_map.csv mirror the 2026 block verbatim.
morphotype_map <- c(
  "GRAM"   = "Other Graminoid",
  "FORB"   = "Other Forb",
  "LITTER" = "Non-Photosynthetic Vegetation",
  "BARE"   = "Nonvegetated Dirt",
  "ROCK"   = "Nonvegetated Rock",
  "SHRUB"  = "Other Deciduous Shrub"
)

# --- Site-name normalization -----------------------------------------------
norm_site <- function(x) {
  x <- stringr::str_to_upper(stringr::str_squish(x))
  stringr::str_replace(x, "^([A-Z]+)-?0*([0-9]+)$", "\\1-\\2")
}

# --- Load ------------------------------------------------------------------
cover_raw <- readr::read_csv(file.path(raw_dir, "vegmap_2023_cover.csv"),
                             show_col_types = FALSE) |>
  dplyr::rename_with(~ stringr::str_remove(.x, "^\\ufeff"))  # BOM on col 1

polys_raw <- sf::st_read(
  file.path(raw_dir, "vegmap_2023_wgs_utm_pixel_select.geojson"),
  quiet = TRUE) |>
  sf::st_transform(32613)

crosswalk_2023 <- readr::read_csv(
  "data/small_reference/taxonomy_crosswalk_2023.csv",
  show_col_types = FALSE) |>
  dplyr::filter(campaign == "2023")

cover_raw <- cover_raw |> dplyr::mutate(site_name = norm_site(Site))
polys     <- polys_raw  |> dplyr::mutate(site_name = norm_site(Name))

stopifnot(!anyDuplicated(polys$site_name))

# --- Report cover <-> polygon mismatches (do not drop either side) ---------
cov_sites  <- unique(cover_raw$site_name)
poly_sites <- unique(polys$site_name)
no_poly  <- setdiff(cov_sites, poly_sites)
no_cover <- setdiff(poly_sites, cov_sites)
if (length(no_poly))
  message("Cover without polygon (awaiting collaborator resolution): ",
          paste(sort(no_poly), collapse = ", "))
if (length(no_cover))
  message("Polygon without cover (awaiting collaborator resolution): ",
          paste(sort(no_cover), collapse = ", "))

# --- Stable site-number assignment (3001+) ---------------------------------
lookup_path <- "data/small_reference/site_numbers_2023.csv"
lookup <- if (file.exists(lookup_path)) {
  readr::read_csv(lookup_path, show_col_types = FALSE) |>
    dplyr::mutate(site_number = as.integer(site_number))
} else {
  tibble::tibble(site_name = character(), site_number = integer())
}
new_names <- setdiff(sort(union(cov_sites, poly_sites)), lookup$site_name)
if (length(new_names)) {
  start <- max(3000L, suppressWarnings(max(lookup$site_number)))
  lookup <- dplyr::bind_rows(
    lookup,
    tibble::tibble(site_name = new_names,
                   site_number = start + seq_along(new_names))
  )
  readr::write_csv(dplyr::arrange(lookup, site_number), lookup_path)
  message(sprintf("Assigned %d new site numbers (now %d total) -> %s",
                  length(new_names), nrow(lookup), lookup_path))
}

# --- Standardized cover table (2026 schema) --------------------------------
unmapped <- setdiff(unique(cover_raw$Species),
                    c(crosswalk_2023$raw_name, names(morphotype_map)))
if (length(unmapped))
  stop("2023 species codes missing from taxonomy_crosswalk_2023.csv: ",
       paste(unmapped, collapse = ", "))

code_to_canonical <- crosswalk_2023 |>
  dplyr::filter(!is.na(canonical_name), canonical_name != "") |>
  dplyr::select(raw_name, canonical_name) |>
  tibble::deframe()
code_to_canonical[names(canonical_overrides)] <- canonical_overrides

cover_std <- cover_raw |>
  dplyr::inner_join(lookup, by = "site_name") |>
  dplyr::transmute(
    Site_Number      = site_number,
    Site_Name        = site_name,
    Collection_Date  = lubridate::mdy(Date),
    Cover_Type       = dplyr::if_else(Species %in% names(morphotype_map),
                                      "Non-Species Cover",
                                      "Live Vegetation - Named Species"),
    Cover_Class_Name = dplyr::coalesce(morphotype_map[Species],
                                       code_to_canonical[Species]),
    Cover_Percent    = as.numeric(Cover)
  )
stopifnot(!any(is.na(cover_std$Cover_Class_Name)))

readr::write_csv(cover_std, "data/derived/vegmap_2023_cover_std.csv")

# --- Standardized polygons -------------------------------------------------
polys_std <- polys |>
  dplyr::inner_join(lookup, by = "site_name") |>
  dplyr::select(site_number, site_name)

out_gj <- "data/derived/vegmap_2023_polygons_std.geojson"
if (file.exists(out_gj)) file.remove(out_gj)
# RFC7946=NO keeps coordinates in EPSG:32613 (with a crs member), which the
# python extractor requires — it reads geometry coords as UTM directly.
sf::st_write(polys_std, out_gj, quiet = TRUE,
             layer_options = c("RFC7946=NO", "WRITE_BBOX=NO"))

# --- Summary ---------------------------------------------------------------
totals <- cover_std |>
  dplyr::summarise(total = sum(Cover_Percent), .by = Site_Number)
message(sprintf(
  "2023 standardized: %d cover sites (totals min=%.0f max=%.0f), %d polygons, %d taxa + %d morphotypes.",
  dplyr::n_distinct(cover_std$Site_Number), min(totals$total), max(totals$total),
  nrow(polys_std),
  dplyr::n_distinct(cover_std$Cover_Class_Name[
    cover_std$Cover_Type == "Live Vegetation - Named Species"]),
  dplyr::n_distinct(cover_std$Cover_Class_Name[
    cover_std$Cover_Type != "Live Vegetation - Named Species"])
))
