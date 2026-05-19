# 41_class_summary_table.R — single-row-per-class summary that pulls
# together training support, basin prevalence (RF predictions on the
# 5,354-pixel inference set), recall, indicator taxa, curated narrative
# names, and median per-pixel leverage. Intended as the canonical
# class-level reference table for sharing with the team.
#
# Inputs:
#   data/derived/punch_list.csv
#   data/derived/sampling_priority.gpkg          (per-pixel leverage)
#   data/small_reference/label_community_names.csv (meadow narratives)
#   data/derived/label_descriptions.csv          (meadow IndVal indicators)
#   data/derived/shrub_training_set.rds          (per-record canonical info)
# Output:
#   data/derived/class_summary_table.csv         the single deliverable

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
})

punch <- readr::read_csv("data/derived/punch_list.csv",
                         show_col_types = FALSE)
spri  <- sf::st_read("data/derived/sampling_priority.gpkg", quiet = TRUE) |>
  sf::st_drop_geometry()
narrs <- readr::read_csv("data/small_reference/label_community_names.csv",
                         show_col_types = FALSE) |>
  dplyr::mutate(
    # M-classes are post-hoc monotypic-species overrides; the species IS
    # the indicator and lives in the label itself (e.g.
    # "M05_Osmorhiza_occidentalis"). Fall back to that when IndVal didn't
    # produce a top_indicator.
    fallback_indicator = dplyr::if_else(
      stringr::str_starts(final_label, "M\\d"),
      stringr::str_replace(stringr::str_replace(final_label,
                                                 "^M\\d+_", ""),
                           "_", " "),
      NA_character_
    )
  ) |>
  dplyr::transmute(final_label,
                   indicator_taxa = dplyr::coalesce(top_indicator,
                                                    fallback_indicator),
                   abundant_taxa  = dplyr::coalesce(top_abundant,
                                                    fallback_indicator),
                   description    = dplyr::coalesce(narrative_curated,
                                                    narrative_draft))

# Shrub side: indicator taxon is just the canonical binomial (each shrub
# site = one species; collapsed labels like "Salix sp." aggregate
# several binomials, list those for the indicator column).
shrub_train <- readRDS("data/derived/shrub_training_set.rds")$training
shrub_inds <- shrub_train |>
  dplyr::distinct(final_label, canonical_binomial) |>
  dplyr::group_by(final_label) |>
  dplyr::summarise(indicator_taxa = paste(sort(unique(canonical_binomial)),
                                          collapse = "; "),
                   .groups = "drop") |>
  dplyr::mutate(abundant_taxa = indicator_taxa,
                description   = final_label)

descriptions <- dplyr::bind_rows(narrs, shrub_inds)

# Per-class median leverage from the per-pixel scores.
med_lev <- spri |>
  dplyr::group_by(predicted_label) |>
  dplyr::summarise(median_leverage = stats::median(leverage, na.rm = TRUE),
                   .groups = "drop") |>
  dplyr::rename(final_label = predicted_label)

summary_tbl <- punch |>
  dplyr::transmute(
    final_label,
    class_type,
    n_2018,
    n_2025,
    n_total,
    predicted_n_pixels,
    pct_of_inference     = round(100 * predicted_n_pixels /
                                   sum(predicted_n_pixels, na.rm = TRUE), 2),
    balanced_recall      = round(balanced_recall, 2),
    augmentation_priority,
    top_confusions
  ) |>
  dplyr::left_join(med_lev,      by = "final_label") |>
  dplyr::left_join(descriptions, by = "final_label") |>
  dplyr::select(final_label, class_type, description,
                indicator_taxa, abundant_taxa,
                n_2018, n_2025, n_total,
                predicted_n_pixels, pct_of_inference,
                balanced_recall, median_leverage,
                augmentation_priority, top_confusions) |>
  dplyr::arrange(
    factor(augmentation_priority,
           levels = c("critical", "high", "medium", "ok")),
    dplyr::desc(median_leverage)
  ) |>
  dplyr::mutate(median_leverage = round(median_leverage, 2))

readr::write_csv(summary_tbl, "data/derived/class_summary_table.csv")
cat(sprintf("Wrote data/derived/class_summary_table.csv (%d rows)\n",
            nrow(summary_tbl)))
cat("\nPreview:\n")
print(as.data.frame(summary_tbl), max = 250)
