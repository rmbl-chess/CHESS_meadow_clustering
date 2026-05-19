# 21_size_distribution.R — cluster-size distribution diagnostic.
#
# Outputs: output/figures/size_distribution.png  (stacked bar of cluster
#          sizes by source, ordered by snow-free DOY, colored by tier).
# Plus stdout summaries by tier and size band.

suppressPackageStartupMessages({
  library(tidyverse)
})

fc   <- readRDS("data/derived/final_clusters_B.rds")
desc <- readr::read_csv("data/small_reference/label_community_names.csv",
                        show_col_types = FALSE)

sizes <- fc$assignments |>
  dplyr::count(final_label, source) |>
  tidyr::pivot_wider(names_from = source, values_from = n, values_fill = 0L) |>
  dplyr::rename(n_anchor = clustered_2025, n_inferred = inferred_2018) |>
  dplyr::mutate(n_total = n_anchor + n_inferred)

joined <- desc |>
  dplyr::select(final_label, recall, tier, snow_free_doy_mean) |>
  dplyr::left_join(sizes, by = "final_label") |>
  dplyr::mutate(tier = factor(tier, levels = c("strong", "marginal", "weak")))

# --- Stacked bar plot -------------------------------------------------------
plot_df <- joined |>
  dplyr::arrange(snow_free_doy_mean) |>
  dplyr::mutate(final_label = factor(final_label, levels = final_label)) |>
  tidyr::pivot_longer(c(n_anchor, n_inferred),
                      names_to = "source", values_to = "n") |>
  dplyr::mutate(source = dplyr::recode(source,
                                       n_anchor = "Clustered (2025)",
                                       n_inferred = "Inferred (2018)"))

tier_colors <- c(strong = "#1b7837", marginal = "#dfc27d", weak = "#bababa")

p <- ggplot2::ggplot(plot_df,
                     ggplot2::aes(x = final_label, y = n, fill = source)) +
  ggplot2::geom_col(width = 0.75) +
  ggplot2::geom_text(
    data = joined |> dplyr::arrange(snow_free_doy_mean) |>
      dplyr::mutate(final_label = factor(final_label, levels = final_label)),
    ggplot2::aes(y = n_total + 1.5, x = final_label,
                 label = sprintf("%.2f", recall), fill = NULL),
    size = 2.9, inherit.aes = FALSE
  ) +
  ggplot2::scale_fill_manual(values = c("Clustered (2025)" = "#3a8fb7",
                                        "Inferred (2018)"  = "#dca066"),
                             name = "Source") +
  ggplot2::labs(
    x = NULL, y = "Site count",
    title = "Cluster size distribution (26-class set), ordered by snow-free DOY",
    subtitle = "Number above bar = CV recall (clustered-2025 only).  Source colors split anchor vs inferred."
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1),
                 legend.position = "top",
                 panel.grid.major.x = ggplot2::element_blank())

dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)
ggplot2::ggsave("output/figures/size_distribution.png", p,
                width = 11, height = 6, dpi = 150)

cat("Wrote output/figures/size_distribution.png\n\n")

cat("=== Size summary by tier ===\n")
print(joined |>
  dplyr::group_by(tier) |>
  dplyr::summarise(
    n_clusters = dplyr::n(),
    anchor_total = sum(n_anchor),
    inferred_total = sum(n_inferred),
    grand_total = sum(n_total),
    min_size = min(n_total),
    median_size = stats::median(n_total),
    max_size = max(n_total),
    .groups = "drop"
  ) |>
  dplyr::arrange(tier))

cat("\n=== Size band distribution ===\n")
print(joined |>
  dplyr::mutate(size_band = dplyr::case_when(
    n_total <  5 ~ "tiny (<5)",
    n_total < 15 ~ "small (5-14)",
    n_total < 30 ~ "moderate (15-29)",
    n_total < 50 ~ "large (30-49)",
    TRUE         ~ "very large (50+)"
  ),
  size_band = factor(size_band,
                     levels = c("tiny (<5)", "small (5-14)",
                                "moderate (15-29)",
                                "large (30-49)", "very large (50+)"))) |>
  dplyr::count(size_band, .drop = FALSE))
