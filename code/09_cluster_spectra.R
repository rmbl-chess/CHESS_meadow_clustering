# 09_cluster_spectra.R — Architecture B step 1: Ward hierarchical clustering
# on the top 12 spectral PCs, swept across k. Each spectral cluster is
# characterized by:
#   - dominant indicator genus (top non-NPV/Bare feature in Hellinger centroid)
#   - dominance score      = max(named_centroid) / sum(named_centroid)
#   - heterogeneity score  = mean Hellinger Euclidean distance from each site
#                            to the cluster's composition centroid
#   - year breakdown       = n_2018, n_2025
#
# Spectra-first means every cluster is mappable by construction. The
# characterization tells us _what_ those mappable clusters resolve to in
# ecological terms.
#
# Inputs:  data/derived/spectral_features.rds, .../composition_genus.rds
# Outputs: data/derived/spectral_clusters.rds  (hclust object + per-k cuts +
#                                                per-k characterization tables)

suppressPackageStartupMessages({
  library(tidyverse)
})

spec_feat  <- readRDS("data/derived/spectral_features.rds")$features
comp_genus <- readRDS("data/derived/composition_genus.rds")

# Run two parallel clusterings:
#   variant_A  PCs 1-12 (full)        — full brightness + greenness signal
#   variant_B  PCs 3-12 (drop 1,2)    — PC1/PC2 typically encode brightness
#                                       and greenness; dropping them isolates
#                                       subtler spectral-shape variation.
ks <- c(4, 5, 6, 7, 8, 10, 12)

run_clustering <- function(pc_idx, label) {
  pc_cols <- sprintf("spec_PC%02d", pc_idx)
  mat     <- as.matrix(spec_feat[, pc_cols])
  d       <- dist(mat, method = "euclidean")
  hc      <- hclust(d, method = "ward.D2")
  cuts <- as_tibble(setNames(
    lapply(ks, function(k) cutree(hc, k = k)),
    sprintf("k%02d", ks)
  ))
  assignments <- bind_cols(spec_feat |> dplyr::select(site_number, Year), cuts)
  list(hclust = hc, dist = d, assignments = assignments,
       pc_cols = pc_cols, label = label)
}

variant_A <- run_clustering(seq_len(12), "PCs_1to12")
variant_B <- run_clustering(3:12,        "PCs_3to12")

# --- Per-cluster characterization helper -----------------------------------
nonsp_set <- c("Other_Forb", "Other_Graminoid", "NPV", "Bare",
               "Other_Moss_Lichen", "Other_Deciduous_Shrub")

hell <- comp_genus$hellinger
hell_feat_cols <- setdiff(names(hell), c("site_number", "Year"))
hell_mat <- as.matrix(hell[, hell_feat_cols])
rownames(hell_mat) <- paste(hell$site_number, hell$Year, sep = "_")

characterize_k <- function(assignments, k_col) {
  asg <- assignments |>
    dplyr::select(site_number, Year, cluster = dplyr::all_of(k_col)) |>
    dplyr::mutate(key = paste(site_number, Year, sep = "_"))

  hell_long <- hell |>
    pivot_longer(-c(site_number, Year), names_to = "feature", values_to = "h") |>
    mutate(feature = stringr::str_replace(feature, "_cover$", "")) |>
    inner_join(asg |> dplyr::select(site_number, Year, cluster),
               by = c("site_number", "Year"))

  centroids <- hell_long |>
    group_by(cluster, feature) |>
    summarise(mean_h = mean(h), .groups = "drop")

  top <- centroids |>
    group_by(cluster) |>
    arrange(desc(mean_h), .by_group = TRUE) |>
    summarise(
      indicator_genus = {
        sp <- feature[!feature %in% nonsp_set]
        if (length(sp) == 0) "(no named genus)" else sp[1]
      },
      top_features = paste(head(feature, 5), collapse = ", "),
      .groups = "drop"
    )

  dominance <- centroids |>
    filter(!feature %in% nonsp_set) |>
    group_by(cluster) |>
    summarise(
      dominance = if (sum(mean_h) > 0) max(mean_h) / sum(mean_h) else 0,
      .groups = "drop"
    )

  # Heterogeneity: mean Euclidean distance from each site (Hellinger vector)
  # to its cluster centroid. Some sites may be missing from hell_mat (e.g.,
  # sites whose entire cover was in rare genera dropped during 05's
  # prevalence trim) — filter to keys present in hell_mat.
  available <- rownames(hell_mat)
  het <- purrr::map_dfr(unique(asg$cluster), function(cl) {
    keys <- intersect(asg$key[asg$cluster == cl], available)
    if (length(keys) < 2) return(tibble(cluster = cl, heterogeneity = NA_real_))
    mat   <- hell_mat[keys, , drop = FALSE]
    cent  <- colMeans(mat)
    dists <- sqrt(rowSums(sweep(mat, 2, cent)^2))
    tibble(cluster = cl, heterogeneity = mean(dists))
  })

  sizes <- asg |> count(cluster, name = "n_sites")
  year_b <- asg |> count(cluster, Year) |>
    pivot_wider(names_from = Year, values_from = n, names_prefix = "n_",
                values_fill = 0L)

  sizes |>
    left_join(year_b, by = "cluster") |>
    left_join(dominance, by = "cluster") |>
    left_join(het, by = "cluster") |>
    left_join(top, by = "cluster") |>
    arrange(desc(n_sites))
}

characterize_variant <- function(variant) {
  setNames(
    lapply(ks, function(k) characterize_k(variant$assignments, sprintf("k%02d", k))),
    sprintf("k%02d", ks)
  )
}

char_A <- characterize_variant(variant_A)
char_B <- characterize_variant(variant_B)

# How much variance do we drop with PC1+PC2? Useful context.
spec_meta <- readRDS("data/derived/spectral_features.rds")
varexp    <- spec_meta$var_explained
ve_1_2    <- varexp[2]
ve_3_12   <- varexp[12] - varexp[2]
cat(sprintf("\nVariance check: PC1+PC2 = %.2f%%; PCs 3-12 = %.2f%%; total PCs 1-12 = %.2f%%\n",
            100 * ve_1_2, 100 * ve_3_12, 100 * varexp[12]))

for (k in ks) {
  cat(sprintf("\n========= k=%d =========\n", k))
  cat("\n[variant A: PCs 1-12]\n")
  print(char_A[[sprintf("k%02d", k)]], n = Inf, width = Inf)
  cat("\n[variant B: PCs 3-12, drop brightness/greenness axes]\n")
  print(char_B[[sprintf("k%02d", k)]], n = Inf, width = Inf)
}

saveRDS(list(
  variant_A = list(hclust = variant_A$hclust, dist = variant_A$dist,
                   assignments = variant_A$assignments,
                   characterizations = char_A, pc_cols = variant_A$pc_cols),
  variant_B = list(hclust = variant_B$hclust, dist = variant_B$dist,
                   assignments = variant_B$assignments,
                   characterizations = char_B, pc_cols = variant_B$pc_cols),
  ks = ks
), "data/derived/spectral_clusters.rds")
