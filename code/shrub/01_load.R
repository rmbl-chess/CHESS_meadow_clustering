# 30_shrub_load.R — assemble a single-species shrub training dataset from
# the 2018 and 2025 CHESS field campaigns.
#
# Design:
#   2025: chess_shrub_site_cleaned.csv is one row per shrub site, each
#         tagged with a single Vegetation_Species in "Family-Genus species"
#         form. Site_Number links to the 2025 crowns (site_type == "Shrub")
#         and per-pixel spectra.
#   2018: shrub records are mixed into fractional_cover (1).csv. Per user
#         direction we keep only sites where ONE shrub-dominated genus has
#         FractionalCover == 100 — i.e., shrub-dominated single-species
#         crowns. Tree genera (Picea, Pinus, Abies, Pseudotsuga, Populus,
#         Larix, Quercus) are excluded entirely. Genus list is below.
#
# Inputs:
#   data/raw/ESS-DIVE-Vegetation-Field-2018/chess_shrub_site_cleaned.csv
#   data/derived/veg_2018.rds  (cover + species_list, written by 01_load.R)
#
# Output:
#   data/derived/shrub_records.rds   list(records, source_summary)
#   data/derived/shrub_records.csv   the same `records` table, for review
#
# The records table is one row per (site_number, Year) with the parsed
# species, intentionally mirroring the meadow training-sample layout so
# downstream scripts can reuse 04_join_spectra-style joins.

suppressPackageStartupMessages({
  library(tidyverse)
})

shrub_genera <- c(
  "Salix", "Acer", "Ribes", "Juniperus", "Alnus", "Betula",
  "Pentaphylloides", "Dasiphora", "Sambucus", "Symphoricarpos",
  "Distegia", "Lonicera", "Vaccinium", "Rosa", "Mahonia", "Berberis",
  "Cornus", "Shepherdia", "Prunus", "Crataegus", "Amelanchier",
  "Holodiscus", "Physocarpus", "Rubus", "Sorbus"
)
# Artemisia is herbaceous except for A. tridentata (big sagebrush), which
# is woody. Handled as a special-case CoverCode below.

# --- 2025 ------------------------------------------------------------------
shrub_2025_raw <- readr::read_csv(
  "data/raw/ESS-DIVE-Vegetation-Field-2018/chess_shrub_site_cleaned.csv",
  show_col_types = FALSE
)
shrub_2025 <- shrub_2025_raw |>
  dplyr::transmute(
    site_number    = as.integer(Site_Number),
    Year           = 2025L,
    sampling_area  = Sampling_Area,
    raw_species    = Vegetation_Species,
    vegetation_height_cm = Vegetation_Height
  ) |>
  tidyr::separate(raw_species, into = c("family", "binomial"),
                  sep = "-", remove = FALSE, extra = "merge", fill = "right") |>
  tidyr::separate(binomial, into = c("genus", "species"),
                  sep = " ", remove = FALSE, extra = "merge", fill = "right") |>
  dplyr::mutate(
    family   = stringr::str_squish(family),
    genus    = stringr::str_squish(genus),
    species  = stringr::str_squish(species),
    binomial = stringr::str_squish(binomial)
  )

cat(sprintf("2025 shrub sites loaded: %d\n", nrow(shrub_2025)))

# --- 2018 ------------------------------------------------------------------
veg_2018 <- readRDS("data/derived/veg_2018.rds")
cover_2018   <- veg_2018$cover
species_2018 <- veg_2018$species

species_2018_clean <- species_2018 |>
  dplyr::transmute(
    CoverCode,
    family   = stringr::str_squish(Family),
    genus    = stringr::str_squish(Genus),
    species  = stringr::str_squish(Species),
    binomial = stringr::str_squish(paste(genus, species))
  )

# Keep 2018 rows that are (a) shrub-dominated genus or A. tridentata
# AND (b) cover == 100 (single-species crown).
shrub_2018 <- cover_2018 |>
  dplyr::filter(FractionalCover == 100) |>
  dplyr::inner_join(species_2018_clean, by = "CoverCode") |>
  dplyr::filter(genus %in% shrub_genera |
                (genus == "Artemisia" & species == "tridentata")) |>
  dplyr::transmute(
    site_number,
    Year           = 2018L,
    sampling_area  = SamplingArea,
    raw_species    = CoverCode,
    family, genus, species, binomial,
    vegetation_height_cm = NA_real_
  )

cat(sprintf("2018 shrub records (single-species, cover==100): %d\n",
            nrow(shrub_2018)))

# --- 2026 supplemental -----------------------------------------------------
# Same cover==100 single-species-crown rule as 2018, but the 2026 cover names
# species directly (Cover_Class_Name is the binomial). No species_list, so
# family is left NA. Spectra come from 2025 AOP (site_type == "Shrub").
veg_2026 <- readRDS("data/derived/veg_2026.rds")
shrub_2026 <- veg_2026$cover |>
  dplyr::filter(Cover_Type == "Live Vegetation - Named Species",
                Cover_Percent == 100) |>
  tidyr::separate(Cover_Class_Name, into = c("genus", "species"),
                  sep = " ", remove = FALSE, extra = "merge", fill = "right") |>
  dplyr::mutate(genus = stringr::str_squish(genus),
                species = stringr::str_squish(species)) |>
  dplyr::filter(genus %in% shrub_genera |
                (genus == "Artemisia" & species == "tridentata")) |>
  dplyr::transmute(
    site_number    = as.integer(site_number),
    Year           = 2026L,
    sampling_area  = Sampling_Area,
    raw_species    = Cover_Class_Name,
    family         = NA_character_,
    genus, species,
    binomial       = stringr::str_squish(Cover_Class_Name),
    vegetation_height_cm = NA_real_
  )

cat(sprintf("2026 shrub records (single-species, cover==100): %d\n",
            nrow(shrub_2026)))

# --- Combine + summarize --------------------------------------------------
records <- dplyr::bind_rows(shrub_2018, shrub_2025, shrub_2026) |>
  dplyr::arrange(Year, site_number)

source_summary <- records |>
  dplyr::count(Year, genus, binomial, name = "n_records") |>
  dplyr::arrange(Year, dplyr::desc(n_records))

cat("\n--- Records per genus per year ---\n")
print(records |> dplyr::count(Year, genus) |>
        tidyr::pivot_wider(names_from = Year, values_from = n,
                           values_fill = 0L))

cat("\n--- Records per binomial per year (top 25) ---\n")
print(source_summary |> dplyr::group_by(Year) |>
        dplyr::slice_head(n = 25))

# --- Persist --------------------------------------------------------------
dir.create("data/derived", showWarnings = FALSE, recursive = TRUE)
saveRDS(list(records = records, source_summary = source_summary),
        "data/derived/shrub_records.rds")
readr::write_csv(records, "data/derived/shrub_records.csv")
cat("\nWrote data/derived/shrub_records.{rds,csv}\n")
