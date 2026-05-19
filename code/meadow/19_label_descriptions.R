# 19_label_descriptions.R — per-label IndVal description + auto-drafted
# narrative. Combines the old 18 + 19 into one script: they share inputs
# (final_clusters_B.rds, environment.rds, veg_spectra.rds) and the second
# half consumes the first half's output directly.
#
# Section 1 — Per-label IndVal stats
#   mean_cover_in_cluster, freq_in_cluster, specificity, indval per
#   species and per non-species category; top-5 indicators + abundant +
#   top-3 physiognomy strings.
#   Output: data/derived/label_descriptions.csv
#
# Section 2 — Auto-drafted starter narratives
#   "{elevation cue} {habitat cue}, {top-indicator}-dominated", with cues
#   derived from snow-free DOY + bare/NPV pct. Preserves curated rows in
#   the existing CSV; CHESS_REGEN_NARRATIVES=1 forces a full regenerate.
#   Output: data/small_reference/label_community_names.csv

suppressPackageStartupMessages({
  library(tidyverse)
})

fc  <- readRDS("data/derived/final_clusters_B.rds")
vs  <- readRDS("data/derived/veg_spectra.rds")$joined
env <- readRDS("data/derived/environment.rds")

# ============================================================================
# Section 1 — label_descriptions.csv
# ============================================================================
nonsp <- c("Other_Forb", "Other_Graminoid", "NPV", "Bare",
           "Other_Moss_Lichen", "Other_Deciduous_Shrub")

cover_cols <- grep("_cover$", names(vs), value = TRUE)
cover_long <- vs |>
  dplyr::select(site_number, Year, dplyr::all_of(cover_cols)) |>
  tidyr::pivot_longer(dplyr::all_of(cover_cols),
                      names_to = "feature", values_to = "cover") |>
  dplyr::mutate(
    feature  = stringr::str_replace(feature, "_cover$", ""),
    is_named = !feature %in% nonsp
  ) |>
  # Only clustered (2025) sites characterize a cluster — inferred 2018
  # sites were assigned by similarity to its centroid, so they're
  # circular for defining the cluster.
  dplyr::inner_join(
    fc$assignments |>
      dplyr::filter(is.na(source) | source == "clustered_2025") |>
      dplyr::select(site_number, Year, final_label),
    by = c("site_number", "Year")
  )

per_cluster <- cover_long |>
  dplyr::group_by(final_label, feature, is_named) |>
  dplyr::summarise(
    mean_cover_in_cluster = mean(cover, na.rm = TRUE),
    freq_in_cluster       = mean(cover > 0, na.rm = TRUE),
    .groups = "drop"
  )

overall <- cover_long |>
  dplyr::group_by(feature) |>
  dplyr::summarise(mean_cover_overall = mean(cover, na.rm = TRUE),
                   .groups = "drop")

stats <- per_cluster |>
  dplyr::left_join(overall, by = "feature") |>
  dplyr::mutate(
    specificity = dplyr::if_else(
      mean_cover_overall > 0,
      mean_cover_in_cluster / mean_cover_overall,
      0
    ),
    indval = specificity * freq_in_cluster
  )

top_indicators <- stats |>
  dplyr::filter(is_named, mean_cover_in_cluster > 0) |>
  dplyr::group_by(final_label) |>
  dplyr::arrange(dplyr::desc(indval), .by_group = TRUE) |>
  dplyr::slice_head(n = 5) |>
  dplyr::summarise(
    indicators = paste(
      sprintf("%s (cov=%.1f%%, freq=%.0f%%, IV=%.1f)",
              stringr::str_replace_all(feature, "_", " "),
              mean_cover_in_cluster, 100 * freq_in_cluster, indval),
      collapse = "; "
    ),
    .groups = "drop"
  )

top_abundant <- stats |>
  dplyr::filter(is_named, mean_cover_in_cluster > 0) |>
  dplyr::group_by(final_label) |>
  dplyr::arrange(dplyr::desc(mean_cover_in_cluster), .by_group = TRUE) |>
  dplyr::slice_head(n = 5) |>
  dplyr::summarise(
    abundant = paste(
      sprintf("%s (%.1f%%)",
              stringr::str_replace_all(feature, "_", " "),
              mean_cover_in_cluster),
      collapse = "; "
    ),
    .groups = "drop"
  )

top_phys <- stats |>
  dplyr::filter(!is_named) |>
  dplyr::group_by(final_label) |>
  dplyr::arrange(dplyr::desc(mean_cover_in_cluster), .by_group = TRUE) |>
  dplyr::slice_head(n = 3) |>
  dplyr::summarise(
    physiognomy = paste(
      sprintf("%s (%.1f%%)", feature, mean_cover_in_cluster),
      collapse = "; "
    ),
    .groups = "drop"
  )

descriptions <- top_indicators |>
  dplyr::left_join(top_abundant, by = "final_label") |>
  dplyr::left_join(top_phys,     by = "final_label")

readr::write_csv(descriptions, "data/derived/label_descriptions.csv")
cat(sprintf("[1] Wrote data/derived/label_descriptions.csv (%d labels)\n",
            nrow(descriptions)))

# Pretty print for human review (sort by recall desc).
recall_lookup <- fc$final_summary |>
  dplyr::select(final_label, recall, n, n_2018, n_2025)
human_table <- descriptions |>
  dplyr::left_join(recall_lookup, by = "final_label") |>
  dplyr::arrange(dplyr::desc(recall))
for (i in seq_len(nrow(human_table))) {
  r <- human_table[i, ]
  cat(sprintf("\n%s   n=%d  (%d/%d 2018/2025)  recall=%.2f\n",
              r$final_label, r$n, r$n_2018, r$n_2025, r$recall))
  cat("  Indicators:  ", r$indicators,  "\n", sep = "")
  cat("  Abundant:    ", r$abundant,    "\n", sep = "")
  cat("  Physiognomy: ", r$physiognomy, "\n", sep = "")
}

# ============================================================================
# Section 2 — label_community_names.csv (auto-drafted starter narratives)
# ============================================================================
env_per_label <- fc$assignments |>
  dplyr::inner_join(env, by = c("site_number", "Year")) |>
  dplyr::group_by(final_label) |>
  dplyr::summarise(snow_free_doy_mean = mean(snow_free_doy, na.rm = TRUE),
                   .groups = "drop")

tier_of <- function(r) dplyr::case_when(
  is.na(r) | r < 0.50 ~ "weak",
  r        < 0.80     ~ "marginal",
  TRUE                ~ "strong"
)

summ <- fc$final_summary |>
  dplyr::transmute(final_label,
                   n_sites = n,
                   n_2018, n_2025,
                   recall = as.numeric(recall),
                   tier   = tier_of(recall)) |>
  dplyr::left_join(env_per_label, by = "final_label")

first_species <- function(top_string) {
  if (is.na(top_string)) return(NA_character_)
  first <- stringr::str_split(top_string, ";\\s*", simplify = TRUE)[1]
  stringr::str_extract(first, "^[^(]+") |> stringr::str_trim()
}
extract_pct <- function(physiognomy, name) {
  m <- stringr::str_match(physiognomy,
                          sprintf("%s \\(([0-9.]+)%%\\)", name))[, 2]
  as.numeric(m)
}

auto <- summ |>
  dplyr::select(final_label, n_sites, n_2018, n_2025, recall, tier,
                snow_free_doy_mean) |>
  dplyr::left_join(descriptions, by = "final_label") |>
  dplyr::mutate(
    top_indicator = purrr::map_chr(indicators, first_species),
    top_abundant  = purrr::map_chr(abundant,   first_species),
    bare_pct      = extract_pct(physiognomy, "Bare"),
    npv_pct       = extract_pct(physiognomy, "NPV"),
    elevation_cue = dplyr::case_when(
      snow_free_doy_mean < 110 ~ "Low-elevation",
      snow_free_doy_mean < 130 ~ "Mid-elevation montane",
      snow_free_doy_mean < 150 ~ "Subalpine",
      snow_free_doy_mean < 165 ~ "Upper subalpine",
      TRUE                     ~ "Late-melt alpine / subalpine"
    ),
    habitat_cue = dplyr::case_when(
      !is.na(bare_pct) & bare_pct > 35 ~ "rocky shrub-steppe",
      !is.na(bare_pct) & bare_pct > 20 ~ "sparse meadow / shrubland",
      !is.na(npv_pct)  & npv_pct  > 25 ~ "dry meadow",
      TRUE                              ~ "meadow"
    ),
    narrative_draft = sprintf("%s %s, %s-dominated",
                              elevation_cue, habitat_cue, top_indicator)
  )

# Preserve user-edited rows. Set CHESS_REGEN_NARRATIVES=1 to overwrite.
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
for (col in c("narrative_draft_existing", "narrative_curated_existing",
              "notes_existing")) {
  if (!col %in% names(out)) out[[col]] <- NA_character_
}
out <- out |>
  dplyr::mutate(
    narrative_draft = if (force_regen) narrative_draft
                       else dplyr::coalesce(narrative_draft_existing,
                                            narrative_draft),
    narrative_curated = narrative_curated_existing,
    notes             = notes_existing
  ) |>
  dplyr::select(-narrative_draft_existing, -narrative_curated_existing,
                -notes_existing) |>
  dplyr::arrange(snow_free_doy_mean)

dir.create("data/small_reference", showWarnings = FALSE, recursive = TRUE)
readr::write_csv(out, out_path)
cat(sprintf("\n[2] Wrote %s (%d labels)\n", out_path, nrow(out)))

cat("\n=== Auto-drafted narratives (sort by snow-free DOY) ===\n\n")
for (i in seq_len(nrow(out))) {
  r <- out[i, ]
  cat(sprintf("%s  n=%d  doy=%.0f  recall=%.2f  [%s]\n  -> %s\n",
              r$final_label, r$n_sites, r$snow_free_doy_mean,
              r$recall, r$tier, r$narrative_draft))
}
