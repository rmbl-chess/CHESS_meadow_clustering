# 43_joint_figures.R — produce summary figures for the joint meadow+shrub
# classifier. Outputs go straight to docs/figures/ since they are
# intended deliverables (not gitignored output/).
#
# Figures (vector PDFs):
#   joint_class_recall.pdf       Per-class N + balanced CV recall
#   joint_leverage_scatter.pdf   Training N vs basin prevalence, leverage
#   joint_feature_space.pdf      DOY x CHM training-site scatter by type
#   joint_confusion.pdf          47x47 row-normalized RF CV confusion
#
# Inputs:
#   data/derived/joint_training_set.rds
#   data/derived/class_summary_table.csv

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ranger)
})

dir.create("docs/figures", recursive = TRUE, showWarnings = FALSE)

js <- readRDS("data/derived/joint_training_set.rds")
train <- js$training
features <- js$feature_cols
summ <- readr::read_csv("data/derived/class_summary_table.csv",
                        show_col_types = FALSE) |>
  dplyr::mutate(
    augmentation_priority = factor(augmentation_priority,
                                   levels = c("critical","high","medium","ok"))
  )

pri_colors <- c(critical = "#d7301f", high = "#fc8d59",
                medium = "#fdcc8a", ok = "#878787")
type_colors <- c(meadow = "#1b7837", shrub = "#762a83")

# ===========================================================================
# 1. Per-class N + recall
# ===========================================================================
long <- summ |>
  dplyr::transmute(final_label, class_type, augmentation_priority,
                   N = n_total, Recall = balanced_recall) |>
  tidyr::pivot_longer(c(N, Recall), names_to = "metric", values_to = "value")

# Order by N (descending) within each class_type so the figure groups
# meadow / shrub blocks clearly.
order_by <- summ |>
  dplyr::arrange(class_type, dplyr::desc(n_total)) |>
  dplyr::pull(final_label)
long <- long |>
  dplyr::mutate(final_label = factor(final_label, levels = rev(order_by)))

p1 <- ggplot(long, aes(x = final_label, y = value,
                        fill = augmentation_priority)) +
  geom_col() +
  geom_text(aes(label = ifelse(metric == "N", as.character(round(value)),
                                sprintf("%.2f", value))),
            hjust = -0.1, size = 2.6) +
  coord_flip() +
  facet_grid(class_type ~ metric, scales = "free", space = "free_y") +
  scale_fill_manual(values = pri_colors, name = "Augmentation\npriority") +
  labs(x = NULL, y = NULL,
       title = "Joint meadow+shrub classifier — training N and CV recall per class",
       subtitle = "Balanced 5-fold site-level CV; 858 sites, 47 classes") +
  theme_minimal(base_size = 9) +
  theme(strip.text = element_text(face = "bold"),
        panel.grid.major.y = element_blank(),
        plot.title.position = "plot")
ggsave("docs/figures/joint_class_recall.pdf", p1,
       width = 11, height = 11, device = cairo_pdf)

# ===========================================================================
# 2. Leverage scatter: training N vs predicted basin pixels
# ===========================================================================
lev <- summ |>
  dplyr::filter(predicted_n_pixels > 0) |>
  dplyr::mutate(
    label_for_plot = dplyr::case_when(
      augmentation_priority %in% c("critical", "high") ~ final_label,
      median_leverage > stats::quantile(median_leverage, 0.85,
                                        na.rm = TRUE) ~ final_label,
      TRUE                                              ~ NA_character_
    )
  )

p2 <- ggplot(lev, aes(x = n_total, y = predicted_n_pixels,
                       color = augmentation_priority,
                       shape = class_type)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_text(aes(label = label_for_plot),
            size = 3, hjust = -0.15, vjust = -0.4,
            check_overlap = TRUE, show.legend = FALSE) +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = pri_colors, name = "Augmentation\npriority") +
  scale_shape_manual(values = c(meadow = 16, shrub = 17), name = "Class type") +
  labs(x = "Training sites (log scale)",
       y = "Predicted pixels in basin (log scale)",
       title = "Leverage: under-trained classes that map a lot of the basin",
       subtitle = "Upper-left = small training set + lots of basin predicted ⇒ highest leverage") +
  theme_minimal(base_size = 11) +
  theme(plot.title.position = "plot")
ggsave("docs/figures/joint_leverage_scatter.pdf", p2,
       width = 9, height = 7, device = cairo_pdf)

# ===========================================================================
# 3. Feature space: DOY × CHM by class_type
# ===========================================================================
feat_plot <- train |>
  dplyr::mutate(class_type = factor(class_type, levels = c("meadow","shrub"))) |>
  dplyr::filter(!is.na(snow_free_doy), !is.na(canopy_height_m))

p3 <- ggplot(feat_plot, aes(x = snow_free_doy, y = canopy_height_m,
                             color = class_type)) +
  geom_point(alpha = 0.55, size = 1.6) +
  scale_color_manual(values = type_colors) +
  scale_y_sqrt(breaks = c(0, 0.5, 1, 2, 5, 10, 20)) +
  labs(x = "Snow-free DOY",
       y = "Canopy height (m, √-scale)",
       color = "Class type",
       title = "Training sites in the DOY × canopy-height covariate space",
       subtitle = "Meadows (low CHM, broad DOY range) vs. shrubs (taller, mostly mid-DOY)") +
  theme_minimal(base_size = 11) +
  theme(plot.title.position = "plot")
ggsave("docs/figures/joint_feature_space.pdf", p3,
       width = 9, height = 6, device = cairo_pdf)

# ===========================================================================
# 4. Joint confusion matrix (re-derived from 5-fold CV)
# ===========================================================================
# Re-run a quick CV so we have predictions per site to build the matrix;
# matches what 37 does but on the same training table that's in js.
X <- as.matrix(train[, features])
y <- factor(train$final_label)
set.seed(42)
fold <- integer(nrow(X))
for (lvl in levels(y)) {
  idx <- which(y == lvl)
  fold[idx] <- ((sample(seq_along(idx)) - 1L) %% 5L) + 1L
}
preds <- factor(rep(NA_character_, nrow(X)), levels = levels(y))
for (f in 1:5) {
  tr <- which(fold != f); te <- which(fold == f)
  tab <- table(y[tr]); cw <- as.numeric(sum(tab) / (length(tab) * tab))
  names(cw) <- names(tab)
  fit <- ranger::ranger(x = X[tr, ], y = y[tr], num.trees = 500,
                        classification = TRUE, class.weights = cw,
                        seed = 42 + f)
  preds[te] <- predict(fit, X[te, ])$predictions
}

# Order classes by class_type then by training-N within type.
class_order <- summ |>
  dplyr::arrange(class_type, dplyr::desc(n_total)) |>
  dplyr::pull(final_label)

conf <- tibble::tibble(truth = y, pred = preds) |>
  dplyr::filter(!is.na(pred)) |>
  dplyr::count(truth, pred, name = "n") |>
  dplyr::group_by(truth) |>
  dplyr::mutate(frac = n / sum(n)) |>
  dplyr::ungroup() |>
  tidyr::complete(truth, pred, fill = list(n = 0L, frac = 0))

p4 <- ggplot(conf,
             aes(x = factor(pred,  levels = class_order),
                 y = factor(truth, levels = rev(class_order)),
                 fill = frac)) +
  geom_tile(color = "grey92") +
  scale_fill_gradient(low = "white", high = "#08306b",
                      name = "Row\nfraction",
                      limits = c(0, 1)) +
  labs(x = "Predicted (5-fold CV)", y = "Truth",
       title = "Joint meadow+shrub RF — row-normalized confusion matrix",
       subtitle = "Classes ordered by class_type, then training-N within type") +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 75, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        panel.grid = element_blank(),
        plot.title.position = "plot") +
  coord_fixed()
ggsave("docs/figures/joint_confusion.pdf", p4,
       width = 12, height = 12, device = cairo_pdf)

cat("Wrote:\n",
    "  docs/figures/joint_class_recall.pdf\n",
    "  docs/figures/joint_leverage_scatter.pdf\n",
    "  docs/figures/joint_feature_space.pdf\n",
    "  docs/figures/joint_confusion.pdf\n", sep = "")
