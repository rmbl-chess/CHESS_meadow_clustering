# 03_combine_cover.R — assemble the combined cover table.
#
# One row per (site_number, Year) with one column per canonical species,
# named `<Genus_species>_cover`. Missing species in a year fill as 0 (treat as
# absent). Non-species cover (Litter, Bare, "Other Forb", etc.) is excluded —
# only named species participate in clustering features.
#
# Inputs:
#   data/derived/veg_2018.rds, data/derived/veg_2025.rds
#   data/derived/taxonomy_crosswalk.csv
# Outputs:
#   data/derived/cover_combined.rds  (site_number, Year, <Spp>_cover, ...)
#   data/derived/cover_combined.csv  (same, for inspection)

library(tidyverse)

veg_2018  <- readRDS("data/derived/veg_2018.rds")
veg_2025  <- readRDS("data/derived/veg_2025.rds")
crosswalk <- readr::read_csv("data/derived/taxonomy_crosswalk.csv",
                             show_col_types = FALSE)
woody     <- readr::read_csv("data/small_reference/woody_taxa.csv",
                             show_col_types = FALSE)

cw <- crosswalk |>
  dplyr::filter(!is.na(canonical_name)) |>
  dplyr::select(campaign, raw_name, canonical_name)

# Woody-pure-pixel filter: 2018 records cover up to 100% (often a single woody
# crown), 2025 does not. Drop 2018 records where (cover == 100) AND the
# canonical name matches a woody entry — either genus alone (species blank in
# the woody list) or genus+species (e.g., Artemisia tridentata only). Woody
# records at lower cover are kept.
is_woody <- function(canonical_names, woody_tbl) {
  parts <- stringr::str_split_fixed(canonical_names, " ", 2)
  g <- parts[, 1]; s <- parts[, 2]
  genus_only <- woody_tbl |> dplyr::filter(is.na(species)) |> dplyr::pull(genus)
  genus_sp   <- woody_tbl |> dplyr::filter(!is.na(species)) |>
                  dplyr::transmute(key = paste(genus, species)) |> dplyr::pull(key)
  g %in% genus_only | paste(g, s) %in% genus_sp
}

# --- 2018: long-form cover keyed by site_number and CoverCode --------------
long_2018 <- veg_2018$cover |>
  dplyr::transmute(
    site_number, raw_name = CoverCode,
    cover = as.numeric(FractionalCover)
  ) |>
  dplyr::inner_join(cw |> dplyr::filter(campaign == "2018") |>
                      dplyr::select(raw_name, canonical_name),
                    by = "raw_name") |>
  dplyr::mutate(Year = 2018L)

n_before <- nrow(long_2018)
long_2018 <- long_2018 |>
  dplyr::filter(!(cover == 100 & is_woody(canonical_name, woody)))
message(sprintf("Dropped %d 2018 woody-pure-pixel records (cover==100 & woody).",
                n_before - nrow(long_2018)))

# --- 2025: filter to named species, then join on Cover_Class_Name ----------
long_2025 <- veg_2025$cover |>
  dplyr::filter(Cover_Type == "Live Vegetation - Named Species") |>
  dplyr::transmute(
    site_number, raw_name = Cover_Class_Name,
    cover = as.numeric(Cover_Percent)
  ) |>
  dplyr::inner_join(cw |> dplyr::filter(campaign == "2025") |>
                      dplyr::select(raw_name, canonical_name),
                    by = "raw_name") |>
  dplyr::mutate(Year = 2025L)

cover_long <- dplyr::bind_rows(long_2018, long_2025) |>
  dplyr::group_by(site_number, Year, canonical_name) |>
  dplyr::summarise(cover = sum(cover, na.rm = TRUE), .groups = "drop")

cover_wide <- cover_long |>
  dplyr::mutate(species_col = paste0(stringr::str_replace_all(canonical_name, "\\s+", "_"),
                                     "_cover")) |>
  dplyr::select(-canonical_name) |>
  tidyr::pivot_wider(names_from = species_col, values_from = cover, values_fill = 0)

saveRDS(cover_wide, "data/derived/cover_combined.rds")
readr::write_csv(cover_wide, "data/derived/cover_combined.csv")
