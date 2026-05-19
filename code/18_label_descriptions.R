# 18_label_descriptions.R — richer per-label cluster descriptions.
#
# For each final_label, computes per-species and per-non-species:
#   mean_cover_in_cluster   mean cover (%) among sites in this cluster
#   freq_in_cluster         proportion of sites with cover > 0
#   mean_cover_overall      mean cover across ALL sites
#   specificity             mean_in_cluster / mean_overall  (>1 = enriched)
#   indval                  specificity * freq  (Dufrêne-Legendre style)
#
# Then produces three top-N summary columns per cluster:
#   indicators   top 5 named species by indval (what's characteristic)
#   abundant     top 5 named species by mean cover (what dominates by mass)
#   physiognomy  top 3 non-species categories by mean cover
#
# Inputs:  data/derived/final_clusters_B.rds, .../veg_spectra.rds
# Outputs: data/derived/label_descriptions.csv (one row per final_label)

suppressPackageStartupMessages({
  library(tidyverse)
})

fc <- readRDS("data/derived/final_clusters_B.rds")
vs <- readRDS("data/derived/veg_spectra.rds")$joined  # has raw cover columns

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
  # Use only clustered (2025) sites for cluster characterization. Inferred
  # 2018 sites should not feed the definition of what a cluster IS, since
  # they were assigned by composition similarity to that cluster's
  # centroid — that's circular.
  dplyr::inner_join(fc$assignments |>
                      dplyr::filter(is.na(source) | source == "clustered_2025") |>
                      dplyr::select(site_number, Year, final_label),
                    by = c("site_number", "Year"))

# Per-(cluster, feature) stats: mean cover in cluster, frequency.
per_cluster <- cover_long |>
  dplyr::group_by(final_label, feature, is_named) |>
  dplyr::summarise(
    mean_cover_in_cluster = mean(cover, na.rm = TRUE),
    freq_in_cluster       = mean(cover > 0, na.rm = TRUE),
    .groups = "drop"
  )

# Overall mean cover across all sites for each feature.
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

# --- Top N indicator species per cluster (named species, by indval) --------
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

# --- Top N abundant species per cluster (by mean cover) --------------------
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

# --- Top non-species (physiognomic) ----------------------------------------
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
cat(sprintf("Wrote data/derived/label_descriptions.csv (%d labels)\n\n",
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
  cat("  Indicators:  ", r$indicators,   "\n", sep = "")
  cat("  Abundant:    ", r$abundant,     "\n", sep = "")
  cat("  Physiognomy: ", r$physiognomy,  "\n", sep = "")
}
