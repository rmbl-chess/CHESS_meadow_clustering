# 31_shrub_taxonomy.R — reconcile shrub species names between the 2018 and
# 2025 campaigns. Output is a per-record `canonical_binomial` that fully
# replaces year-specific raw names downstream.
#
# Known synonyms (preferred = modern accepted name on the right):
#   Pentaphylloides floribunda  -> Dasiphora fruticosa
#   Distegia involucrata        -> Lonicera involucrata
#
# Sambucus microbotrys (2018) is treated as a synonym of Sambucus racemosa
# (2025) — Flora of Colorado lumps them as S. racemosa subsp. microbotrys.
# Flagged via `note` so it can be unrolled if the user disagrees.
#
# Salix "spp" entries (n=4, 2018) lack a species epithet — kept under the
# generic "Salix sp." canonical so they can either be dropped from
# classifier training or rolled into a genus-level class downstream.
#
# Inputs:
#   data/derived/shrub_records.rds        (one row per (site, year))
#
# Outputs:
#   data/derived/shrub_records_canonical.rds  records with canonical_binomial
#   data/derived/shrub_taxonomy_crosswalk.csv per-(raw_binomial, year) crosswalk

suppressPackageStartupMessages({
  library(tidyverse)
})

records <- readRDS("data/derived/shrub_records.rds")$records

# --- Crosswalk rules -------------------------------------------------------
# Map raw_binomial (Year-specific, as parsed in 30) -> canonical_binomial.
# Anything not listed here is kept as-is via the default branch below.
rename_rules <- tibble::tribble(
  ~raw_binomial,                    ~canonical_binomial,         ~note,
  "Pentaphylloides floribunda",     "Dasiphora fruticosa",       "Pentaphylloides is a junior synonym of Dasiphora",
  "Distegia involucrata",           "Lonicera involucrata",      "Distegia is a junior synonym of Lonicera",
  "Sambucus microbotrys",           "Sambucus racemosa",         "S. microbotrys treated as S. racemosa ssp. microbotrys",
  "Salix spp",                      "Salix sp.",                 "2018 generic Salix — no species epithet"
)

records_canon <- records |>
  dplyr::mutate(raw_binomial = binomial) |>
  dplyr::left_join(rename_rules, by = "raw_binomial") |>
  dplyr::mutate(
    canonical_binomial = dplyr::coalesce(canonical_binomial, raw_binomial),
    canonical_genus    = stringr::word(canonical_binomial, 1),
    note               = dplyr::if_else(is.na(note), "", note)
  )

# --- Crosswalk: one row per (Year, raw_binomial) for review ---------------
crosswalk <- records_canon |>
  dplyr::count(Year, raw_binomial, canonical_binomial, note,
               name = "n_records") |>
  dplyr::arrange(canonical_binomial, Year)

# --- Class summary on canonical labels ------------------------------------
class_summary <- records_canon |>
  dplyr::count(canonical_genus, canonical_binomial, Year,
               name = "n_records") |>
  tidyr::pivot_wider(names_from = Year, values_from = n_records,
                     names_prefix = "n_", values_fill = 0L) |>
  dplyr::mutate(n_total = rowSums(dplyr::across(dplyr::starts_with("n_")))) |>
  dplyr::arrange(canonical_genus, dplyr::desc(n_total))

cat("=== Canonical labels (after rename rules) ===\n")
print(as.data.frame(class_summary))
cat(sprintf("\nTotal records: %d  Unique canonical binomials: %d  Genera: %d\n",
            nrow(records_canon), dplyr::n_distinct(records_canon$canonical_binomial),
            dplyr::n_distinct(records_canon$canonical_genus)))

# --- Persist --------------------------------------------------------------
saveRDS(records_canon, "data/derived/shrub_records_canonical.rds")
readr::write_csv(crosswalk, "data/derived/shrub_taxonomy_crosswalk.csv")
cat("\nWrote data/derived/shrub_records_canonical.rds\n")
cat("Wrote data/derived/shrub_taxonomy_crosswalk.csv\n")
