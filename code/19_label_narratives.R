# 19_label_narratives.R — generate a starter CSV for curated community-type
# narratives. Combines label_descriptions.csv + training_labels_summary.csv,
# auto-drafts a one-line narrative for each cluster using simple rules
# (elevation cue from snow_free_doy + habitat cue from physiognomy +
# dominant indicator), and writes to data/small_reference/.
#
# If data/small_reference/label_community_names.csv already exists, this
# script preserves the user-curated rows and only updates the auto-drafted
# columns. The narrative_curated column is the canonical name once filled.
#
# Inputs:
#   data/derived/label_descriptions.csv
#   data/derived/training_labels_summary.csv
# Output:
#   data/small_reference/label_community_names.csv  (committed; user-curated)

suppressPackageStartupMessages({
  library(tidyverse)
})

desc <- readr::read_csv("data/derived/label_descriptions.csv", show_col_types = FALSE)
summ <- readr::read_csv("data/derived/training_labels_summary.csv", show_col_types = FALSE)

# Helper: strip the "(cov=..., freq=..., IV=...)" tails to get clean species name.
first_species <- function(top_string) {
  if (is.na(top_string)) return(NA_character_)
  first <- stringr::str_split(top_string, ";\\s*", simplify = TRUE)[1]
  stringr::str_extract(first, "^[^(]+") |> stringr::str_trim()
}

# Helper: parse pct value out of a "Bare (37.6%)" type string.
extract_pct <- function(physiognomy, name) {
  m <- stringr::str_match(physiognomy,
                          sprintf("%s \\(([0-9.]+)%%\\)", name))[, 2]
  as.numeric(m)
}

# Build auto-draft narratives.
auto <- summ |>
  dplyr::select(final_label, n_sites, n_2018, n_2025, recall, tier,
                snow_free_doy_mean) |>
  dplyr::left_join(desc, by = "final_label") |>
  dplyr::mutate(
    top_indicator      = purrr::map_chr(indicators, first_species),
    top_abundant       = purrr::map_chr(abundant,   first_species),
    bare_pct           = extract_pct(physiognomy, "Bare"),
    npv_pct            = extract_pct(physiognomy, "NPV"),
    elevation_cue      = dplyr::case_when(
      snow_free_doy_mean < 110           ~ "Low-elevation",
      snow_free_doy_mean < 130           ~ "Mid-elevation montane",
      snow_free_doy_mean < 150           ~ "Subalpine",
      snow_free_doy_mean < 165           ~ "Upper subalpine",
      TRUE                                ~ "Late-melt alpine / subalpine"
    ),
    habitat_cue        = dplyr::case_when(
      !is.na(bare_pct) & bare_pct > 35  ~ "rocky shrub-steppe",
      !is.na(bare_pct) & bare_pct > 20  ~ "sparse meadow / shrubland",
      !is.na(npv_pct)  & npv_pct  > 25  ~ "dry meadow",
      TRUE                              ~ "meadow"
    ),
    narrative_draft = sprintf("%s %s, %s-dominated",
                              elevation_cue, habitat_cue, top_indicator)
  )

# Preserve user-edited columns (narrative_draft, narrative_curated, notes).
# narrative_draft is auto-generated on first run but treated as user-editable
# afterwards. Re-running this script never clobbers a non-NA hand-tuned value.
# To force regeneration of all narrative_draft values, set the env var
# CHESS_REGEN_NARRATIVES=1 before running.
out_path <- "data/small_reference/label_community_names.csv"
existing <- if (file.exists(out_path)) {
  readr::read_csv(out_path, show_col_types = FALSE) |>
    dplyr::select(final_label,
                  dplyr::any_of(c("narrative_draft", "narrative_curated", "notes"))) |>
    dplyr::rename_with(~ paste0(.x, "_existing"), -final_label)
} else {
  tibble::tibble(final_label = character())
}

force_regen <- isTRUE(nzchar(Sys.getenv("CHESS_REGEN_NARRATIVES", "")))

out <- auto |>
  dplyr::select(final_label, n_sites, n_2018, n_2025, recall, tier,
                snow_free_doy_mean, top_indicator, top_abundant,
                bare_pct, npv_pct, narrative_draft) |>
  dplyr::left_join(existing, by = "final_label")
# Ensure the _existing columns are present (left_join doesn't add them if
# `existing` is empty).
for (col in c("narrative_draft_existing", "narrative_curated_existing",
              "notes_existing")) {
  if (!col %in% names(out)) out[[col]] <- NA_character_
}
out <- out |>
  dplyr::mutate(
    narrative_draft = if (force_regen) {
      narrative_draft
    } else {
      dplyr::coalesce(narrative_draft_existing, narrative_draft)
    },
    narrative_curated = narrative_curated_existing,
    notes             = notes_existing
  ) |>
  dplyr::select(-narrative_draft_existing, -narrative_curated_existing,
                -notes_existing) |>
  dplyr::arrange(snow_free_doy_mean)

dir.create("data/small_reference", showWarnings = FALSE, recursive = TRUE)
readr::write_csv(out, out_path)
cat(sprintf("Wrote %s (%d labels)\n", out_path, nrow(out)))

cat("\n=== Auto-drafted narratives (sort by snow-free DOY) ===\n\n")
for (i in seq_len(nrow(out))) {
  r <- out[i, ]
  cat(sprintf("%s  n=%d  doy=%.0f  recall=%.2f  [%s]\n  -> %s\n",
              r$final_label, r$n_sites, r$snow_free_doy_mean,
              r$recall, r$tier, r$narrative_draft))
}
