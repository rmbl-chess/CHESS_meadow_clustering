# 03_combine_cover.R — assemble the combined cover table.
#
# Produces one row per (SiteID, Year) with one column per canonical species,
# named `<Spp>_cover`. Missing species in a year fill as 0 (absent) — flag any
# species recorded in only one campaign so they can be excluded if needed.
#
# Inputs:
#   data/derived/veg_2018.rds, data/derived/veg_2025.rds
#   data/derived/taxonomy_crosswalk.csv
# Outputs:
#   data/derived/cover_combined.rds  (SiteID, Year, <Spp1>_cover, <Spp2>_cover, ...)
#   data/derived/cover_combined.csv  (same, for inspection)

library(tidyverse)

veg_2018  <- readRDS("data/derived/veg_2018.rds")
veg_2025  <- readRDS("data/derived/veg_2025.rds")
crosswalk <- readr::read_csv("data/derived/taxonomy_crosswalk.csv")

# TODO: pivot each campaign's cover table to long form (SiteID, raw_name, cover),
#       then join canonical_name from the crosswalk, drop NAs / unresolved names.
long_2018 <- tibble::tibble(SiteID = character(), raw_name = character(), cover = double()) |>
  dplyr::mutate(Year = 2018L)
long_2025 <- tibble::tibble(SiteID = character(), raw_name = character(), cover = double()) |>
  dplyr::mutate(Year = 2025L)

# TODO: confirm SiteID scheme is consistent across years (or apply site
#       crosswalk before stacking).

cover_long <- dplyr::bind_rows(long_2018, long_2025) |>
  dplyr::left_join(crosswalk |> dplyr::select(raw_name, canonical_name),
                   by = "raw_name") |>
  dplyr::filter(!is.na(canonical_name))

# Sum within (SiteID, Year, canonical_name) in case the crosswalk lumps species.
cover_long <- cover_long |>
  dplyr::group_by(SiteID, Year, canonical_name) |>
  dplyr::summarise(cover = sum(cover, na.rm = TRUE), .groups = "drop")

cover_wide <- cover_long |>
  dplyr::mutate(species_col = paste0(canonical_name, "_cover")) |>
  dplyr::select(-canonical_name) |>
  tidyr::pivot_wider(names_from = species_col, values_from = cover, values_fill = 0)

saveRDS(cover_wide, "data/derived/cover_combined.rds")
readr::write_csv(cover_wide, "data/derived/cover_combined.csv")
