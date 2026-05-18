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

spec_feat   <- readRDS("data/derived/spectral_features.rds")$features
comp_genus  <- readRDS("data/derived/composition_genus.rds")
comp_species <- readRDS("data/derived/composition_species.rds")

# Run three parallel clusterings:
#   variant_A  PCs 1-12 (full)       — full brightness + greenness signal.
#                                      Diagnostic: clusters tend to stratify
#                                      by bare-soil ↔ vegetation brightness,
#                                      which isn't the ecological axis we
#                                      want to map.
#   variant_B  PCs 3-12 (drop 1,2)   — Variant B isolates the ~1.5% of
#                                      variance after brightness+greenness;
#                                      dominated by year acquisition
#                                      artifacts (clusters split by year).
#   variant_C  PCs 2-12 (drop 1)     — Drop only the brightness axis; keep
#                                      greenness (PC2) and beyond. Should
#                                      cluster on shape variation rather
#                                      than albedo.
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
variant_C <- run_clustering(2:12,        "PCs_2to12")

# --- Per-cluster characterization helper -----------------------------------
# Characterization uses SPECIES-level composition: more granular ecological
# indicators (e.g. "Veratrum tenuipetalum" not just "Veratrum"). Genus-level
# indicator shown as a secondary column for cross-reference. Dominance and
# heterogeneity are computed on species-level Hellinger.
nonsp_set <- c("Other_Forb", "Other_Graminoid", "NPV", "Bare",
               "Other_Moss_Lichen", "Other_Deciduous_Shrub")

hell_sp <- comp_species$hellinger
hell_sp_feat <- setdiff(names(hell_sp), c("site_number", "Year"))
hell_sp_mat  <- as.matrix(hell_sp[, hell_sp_feat])
rownames(hell_sp_mat) <- paste(hell_sp$site_number, hell_sp$Year, sep = "_")

hell_g <- comp_genus$hellinger

characterize_k <- function(assignments, k_col) {
  asg <- assignments |>
    dplyr::select(site_number, Year, cluster = dplyr::all_of(k_col)) |>
    dplyr::mutate(key = paste(site_number, Year, sep = "_"))

  sp_long <- hell_sp |>
    pivot_longer(-c(site_number, Year), names_to = "feature", values_to = "h") |>
    mutate(feature = stringr::str_replace(feature, "_cover$", "")) |>
    inner_join(asg |> dplyr::select(site_number, Year, cluster),
               by = c("site_number", "Year"))
  sp_centroids <- sp_long |>
    group_by(cluster, feature) |>
    summarise(mean_h = mean(h), .groups = "drop")

  top_sp <- sp_centroids |>
    group_by(cluster) |>
    arrange(desc(mean_h), .by_group = TRUE) |>
    summarise(
      indicator_species = {
        sp <- feature[!feature %in% nonsp_set]
        if (length(sp) == 0) "(no named species)" else sp[1]
      },
      top_features_sp = paste(head(feature, 5), collapse = ", "),
      .groups = "drop"
    )

  g_long <- hell_g |>
    pivot_longer(-c(site_number, Year), names_to = "feature", values_to = "h") |>
    mutate(feature = stringr::str_replace(feature, "_cover$", "")) |>
    inner_join(asg |> dplyr::select(site_number, Year, cluster),
               by = c("site_number", "Year"))
  top_g <- g_long |>
    group_by(cluster, feature) |>
    summarise(mean_h = mean(h), .groups = "drop") |>
    group_by(cluster) |>
    arrange(desc(mean_h), .by_group = TRUE) |>
    summarise(
      indicator_genus = {
        gn <- feature[!feature %in% nonsp_set]
        if (length(gn) == 0) "(no named genus)" else gn[1]
      },
      .groups = "drop"
    )

  dominance <- sp_centroids |>
    filter(!feature %in% nonsp_set) |>
    group_by(cluster) |>
    summarise(
      dominance = if (sum(mean_h) > 0) max(mean_h) / sum(mean_h) else 0,
      .groups = "drop"
    )

  # Heterogeneity uses species-level Hellinger.
  available <- rownames(hell_sp_mat)
  het <- purrr::map_dfr(unique(asg$cluster), function(cl) {
    keys <- intersect(asg$key[asg$cluster == cl], available)
    if (length(keys) < 2) return(tibble(cluster = cl, heterogeneity = NA_real_))
    mat   <- hell_sp_mat[keys, , drop = FALSE]
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
    left_join(top_sp, by = "cluster") |>
    left_join(top_g,  by = "cluster") |>
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
char_C <- characterize_variant(variant_C)

# Variance breakdown for context.
spec_meta <- readRDS("data/derived/spectral_features.rds")
varexp    <- spec_meta$var_explained
ve_PC1    <- varexp[1]
ve_PC2    <- varexp[2] - varexp[1]
ve_3_12   <- varexp[12] - varexp[2]
cat(sprintf("\nVariance check: PC1=%.2f%%  PC2=%.2f%%  PCs 3-12=%.2f%%  total PCs 1-12=%.2f%%\n",
            100 * ve_PC1, 100 * ve_PC2, 100 * ve_3_12, 100 * varexp[12]))

for (k in ks) {
  cat(sprintf("\n========= k=%d =========\n", k))
  cat("\n[variant A: PCs 1-12]\n")
  print(char_A[[sprintf("k%02d", k)]], n = Inf, width = Inf)
  cat("\n[variant C: PCs 2-12, drop brightness only]\n")
  print(char_C[[sprintf("k%02d", k)]], n = Inf, width = Inf)
  cat("\n[variant B: PCs 3-12, drop brightness AND greenness]\n")
  print(char_B[[sprintf("k%02d", k)]], n = Inf, width = Inf)
}

saveRDS(list(
  variant_A = list(hclust = variant_A$hclust, dist = variant_A$dist,
                   assignments = variant_A$assignments,
                   characterizations = char_A, pc_cols = variant_A$pc_cols),
  variant_B = list(hclust = variant_B$hclust, dist = variant_B$dist,
                   assignments = variant_B$assignments,
                   characterizations = char_B, pc_cols = variant_B$pc_cols),
  variant_C = list(hclust = variant_C$hclust, dist = variant_C$dist,
                   assignments = variant_C$assignments,
                   characterizations = char_C, pc_cols = variant_C$pc_cols),
  ks = ks
), "data/derived/spectral_clusters.rds")
