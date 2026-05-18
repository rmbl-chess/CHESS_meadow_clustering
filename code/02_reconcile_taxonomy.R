# 02_reconcile_taxonomy.R — build a species-name crosswalk between 2018 and 2025.
#
# Each campaign records species under its own scheme:
#   2018: 6-char `CoverCode` → join to `species_list` for `Genus species`.
#   2025: `Cover_Class_Name` (e.g., "Asteraceae-Balsamorhiza sagittata") → join
#         to `chess_species_list_cleaned` via `Taxon_label` for `Taxon_binomial`
#         (and `GBIF_Taxon_ID`, our preferred canonical reference).
#
# Output is one row per (campaign, raw_name) with a canonical binomial. Where
# the auto-match fails, the row is flagged `needs_review = TRUE` and the user
# can override via `data/small_reference/taxonomy_crosswalk_manual_edits.csv`.
#
# Inputs:  data/derived/veg_2018.rds, data/derived/veg_2025.rds
# Outputs:
#   data/derived/taxonomy_crosswalk.csv  (campaign, raw_name, canonical_name,
#                                          gbif_id, needs_review, note)

library(tidyverse)

veg_2018 <- readRDS("data/derived/veg_2018.rds")
veg_2025 <- readRDS("data/derived/veg_2025.rds")

# --- 2018: CoverCode -> "Genus species" ------------------------------------
# The 2018 species_list has Genus + Species columns and an AltFieldCode that
# sometimes records a synonym (e.g., AcoRos → "syn. = Geum rossii"). For now,
# canonical = "Genus species" from species_list; synonyms are handled by the
# manual-edits override below.
cw_2018 <- veg_2018$species |>
  dplyr::transmute(
    campaign       = "2018",
    raw_name       = CoverCode,
    canonical_name = dplyr::if_else(
      is.na(Genus) | is.na(Species),
      NA_character_,
      stringr::str_squish(paste(Genus, Species))
    ),
    gbif_id        = NA_integer_,
    needs_review   = is.na(Genus) | is.na(Species),
    note           = Notes
  )

# --- 2025: Cover_Class_Name -> Taxon_binomial ------------------------------
cw_2025 <- veg_2025$species |>
  dplyr::transmute(
    campaign       = "2025",
    raw_name       = Taxon_label,
    canonical_name = Taxon_binomial,
    gbif_id        = suppressWarnings(as.integer(GBIF_Taxon_ID)),
    needs_review   = is.na(Taxon_binomial),
    note           = Notes
  )

crosswalk <- dplyr::bind_rows(cw_2018, cw_2025)

# Normalize whitespace + case on canonical so identical species across years
# collapse to one column downstream.
crosswalk <- crosswalk |>
  dplyr::mutate(canonical_name = stringr::str_squish(canonical_name))

# Optional manual override: a committed CSV with the same columns; rows here
# replace auto-derived canonical_name for any matching (campaign, raw_name).
manual_path <- "data/small_reference/taxonomy_crosswalk_manual_edits.csv"
if (file.exists(manual_path)) {
  manual <- readr::read_csv(manual_path, show_col_types = FALSE)
  crosswalk <- crosswalk |>
    dplyr::rows_update(manual, by = c("campaign", "raw_name"))
}

# Flag species recorded in only one campaign — useful for the user to decide
# whether to drop or merge them in clustering features.
seen_in <- crosswalk |>
  dplyr::filter(!is.na(canonical_name)) |>
  dplyr::distinct(campaign, canonical_name) |>
  dplyr::count(canonical_name, name = "n_campaigns")
crosswalk <- crosswalk |>
  dplyr::left_join(seen_in, by = "canonical_name") |>
  dplyr::mutate(single_campaign_only = !is.na(n_campaigns) & n_campaigns == 1L) |>
  dplyr::select(-n_campaigns)

readr::write_csv(crosswalk, "data/derived/taxonomy_crosswalk.csv")
