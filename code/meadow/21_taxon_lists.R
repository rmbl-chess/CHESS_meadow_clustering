# 21_taxon_lists.R — per-class taxon lists ranked by aggregate abundance.
#
# For every meadow class (by class code), every taxon observed in its plots,
# ranked by mean cover across the class's plots (zeros included, so the rank
# reflects aggregate abundance, not just presence). Frequency and per-plot
# mean-when-present are carried so reviewers can tell constant low-cover taxa
# from patchy dominants. Follows 19's convention: clustered sites only (the
# few inferred 2018 fallback sites were assigned BY composition, so including
# them would be circular). Shrub classes are species-defined and not listed.
#
# Inputs:  data/derived/final_clusters_B.rds, data/derived/cover_combined.rds
# Output:  data/derived/class_taxon_lists.csv
#          (final_label, internal_label, rank, taxon, mean_cover_pct,
#           freq_pct, mean_cover_when_present_pct, n_plots)

suppressPackageStartupMessages(library(tidyverse))

fc  <- readRDS("data/derived/final_clusters_B.rds")
cov <- readRDS("data/derived/cover_combined.rds")

asg <- fc$assignments |>
  dplyr::filter(source == "clustered") |>
  dplyr::select(site_number, Year, final_label, internal_label)

nonsp <- paste0(c("Other_Forb", "Other_Graminoid", "NPV", "Bare",
                  "Other_Moss_Lichen", "Other_Deciduous_Shrub"), "_cover")
sp_cols <- setdiff(grep("_cover$", names(cov), value = TRUE), nonsp)

long <- asg |>
  dplyr::inner_join(cov, by = c("site_number", "Year")) |>
  tidyr::pivot_longer(dplyr::all_of(sp_cols),
                      names_to = "taxon", values_to = "cover") |>
  dplyr::mutate(taxon = stringr::str_replace_all(
                          stringr::str_remove(taxon, "_cover$"), "_", " "),
                cover = tidyr::replace_na(cover, 0))

lists <- long |>
  dplyr::group_by(final_label, internal_label, taxon) |>
  dplyr::summarise(
    n_plots        = sum(cover > 0),
    mean_cover_pct = mean(cover),
    freq_pct       = 100 * mean(cover > 0),
    mean_cover_when_present_pct = ifelse(any(cover > 0),
                                         mean(cover[cover > 0]), NA_real_),
    .groups = "drop") |>
  dplyr::filter(n_plots > 0) |>
  dplyr::group_by(final_label) |>
  dplyr::arrange(dplyr::desc(mean_cover_pct), dplyr::desc(freq_pct),
                 taxon, .by_group = TRUE) |>
  dplyr::mutate(rank = dplyr::row_number()) |>
  dplyr::ungroup() |>
  dplyr::mutate(dplyr::across(dplyr::ends_with("_pct"), ~ round(.x, 2))) |>
  dplyr::select(final_label, internal_label, rank, taxon, mean_cover_pct,
                freq_pct, mean_cover_when_present_pct, n_plots)

readr::write_csv(lists, "data/derived/class_taxon_lists.csv")
cat(sprintf("Wrote class_taxon_lists.csv: %d taxon rows across %d classes (median %d taxa/class)\n",
            nrow(lists), dplyr::n_distinct(lists$final_label),
            round(median(table(lists$final_label)))))
