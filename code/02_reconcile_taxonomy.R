# 02_reconcile_taxonomy.R — build a species-name crosswalk between 2018 and 2025.
#
# The two campaigns' species lists mostly overlap, but expect drift: synonyms,
# splits/lumps, typos, authority changes. Output is a canonical name per
# (campaign, raw_name) so downstream cover columns can use one harmonized name.
#
# Inputs:  data/derived/veg_2018.rds, data/derived/veg_2025.rds
# Outputs:
#   data/derived/taxonomy_crosswalk.csv  (campaign, raw_name, canonical_name, note)
#   data/small_reference/taxonomy_crosswalk_manual_edits.csv  (committed; hand
#       edits that override auto-matches — load and join over auto crosswalk)

library(tidyverse)

veg_2018 <- readRDS("data/derived/veg_2018.rds")
veg_2025 <- readRDS("data/derived/veg_2025.rds")

# TODO: identify the column holding the species name in each list (likely
# `taxon` / `scientific_name` / `species`); standardize to `raw_name`.
names_2018 <- veg_2018$species |> dplyr::transmute(raw_name = NA_character_)  # FIXME
names_2025 <- veg_2025$species |> dplyr::transmute(raw_name = NA_character_)  # FIXME

# Strategy:
#   1. Normalize whitespace, case, and authority strings.
#   2. Exact match on the normalized name → canonical.
#   3. Flag unmatched names for manual review (these go into the manual-edits CSV).
#   4. Final canonical name = manual override if present, else exact match,
#      else raw name with a `needs_review = TRUE` flag.

# TODO: implement normalization + match.

crosswalk <- tibble::tibble(
  campaign       = character(),
  raw_name       = character(),
  canonical_name = character(),
  note           = character()
)

readr::write_csv(crosswalk, "data/derived/taxonomy_crosswalk.csv")
