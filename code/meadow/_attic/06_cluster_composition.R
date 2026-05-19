# 06_cluster_composition.R — hierarchical Ward clustering on Hellinger-
# transformed composition, run side-by-side at species and genus levels.
# Euclidean distance on Hellinger-transformed proportions IS the Hellinger
# distance (Legendre & Gallagher 2001).
#
# Cluster assignments are produced at multiple k so the spectral-separability
# step (Phase 2) can pick the right resolution. The "right" k will likely be
# different at species vs genus level — that's the point of running both.
#
# Inputs:  data/derived/composition_species.rds, .../composition_genus.rds
# Outputs: data/derived/clusters_species.rds, .../clusters_genus.rds
#          plus stdout diagnostics (cluster sizes + top features per cluster).

library(tidyverse)

cs <- readRDS("data/derived/composition_species.rds")
cg <- readRDS("data/derived/composition_genus.rds")

ks <- c(5, 10, 15, 20, 25, 30)

cluster_one <- function(comp, label) {
  hell <- comp$hellinger
  feat_cols <- setdiff(names(hell), c("site_number", "Year"))
  d <- dist(as.matrix(hell[, feat_cols]), method = "euclidean")
  hc <- hclust(d, method = "ward.D2")
  cuts <- as_tibble(setNames(
    lapply(ks, function(k) cutree(hc, k = k)),
    sprintf("k%02d", ks)
  ))
  assignments <- bind_cols(hell |> select(site_number, Year), cuts)
  message(sprintf("[%s] n_sites=%d, n_features=%d, hclust height range %.2f–%.2f",
                  label, nrow(hell), length(feat_cols),
                  min(hc$height), max(hc$height)))
  list(hclust = hc, dist = d, assignments = assignments,
       feature_cols = feat_cols, label = label)
}

cl_species <- cluster_one(cs, "species")
cl_genus   <- cluster_one(cg, "genus")

saveRDS(cl_species, "data/derived/clusters_species.rds")
saveRDS(cl_genus,   "data/derived/clusters_genus.rds")

# --- Diagnostics: cluster sizes + top features per cluster ---------------
diag_cluster <- function(cl, comp, k, top_n = 4) {
  k_col <- sprintf("k%02d", k)
  asg <- cl$assignments |>
    select(site_number, Year, cluster = all_of(k_col))
  feat_cols <- cl$feature_cols

  hell_long <- comp$hellinger |>
    select(site_number, Year, all_of(feat_cols)) |>
    pivot_longer(all_of(feat_cols), names_to = "feature", values_to = "h") |>
    inner_join(asg, by = c("site_number", "Year"))

  top <- hell_long |>
    group_by(cluster, feature) |>
    summarise(mean_h = mean(h), .groups = "drop") |>
    group_by(cluster) |>
    slice_max(mean_h, n = top_n, with_ties = FALSE) |>
    summarise(top_features = paste(stringr::str_replace(feature, "_cover$", ""),
                                   collapse = ", "), .groups = "drop")

  sizes <- asg |>
    count(cluster, name = "n") |>
    left_join(asg |> count(cluster, Year) |>
                pivot_wider(names_from = Year, values_from = n,
                            names_prefix = "n_", values_fill = 0L),
              by = "cluster")
  sizes |> left_join(top, by = "cluster") |> arrange(desc(n))
}

for (k in c(10, 15, 20)) {
  cat(sprintf("\n=== Species clustering, k=%d ===\n", k))
  print(diag_cluster(cl_species, cs, k), n = Inf, width = Inf)
  cat(sprintf("\n=== Genus clustering, k=%d ===\n", k))
  print(diag_cluster(cl_genus, cg, k), n = Inf, width = Inf)
}
