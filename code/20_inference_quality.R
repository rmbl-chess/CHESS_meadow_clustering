# 20_inference_quality.R — diagnostic for the 2018 nearest-centroid
# inference. Surfaces per-cluster and per-site quality so we can catch
# mis-assigned 2018 sites going forward.
#
# Inputs:  data/derived/final_clusters_B.rds, .../environment.rds
# Outputs: stdout report (per-cluster summary + worst-fit sites)

suppressPackageStartupMessages({
  library(tidyverse)
})

fc  <- readRDS("data/derived/final_clusters_B.rds")
env <- readRDS("data/derived/environment.rds")
asg <- fc$assignments |> dplyr::inner_join(env, by = c("site_number", "Year"))

clusters_with_inf <- asg |>
  dplyr::filter(source == "inferred_2018") |>
  dplyr::distinct(final_label) |>
  dplyr::pull(final_label)

# --- Per-cluster summary ----------------------------------------------------
per_cluster <- purrr::map_dfr(clusters_with_inf, function(cl) {
  rows <- asg |> dplyr::filter(final_label == cl)
  a <- rows |> dplyr::filter(source == "clustered_2025")
  i <- rows |> dplyr::filter(source == "inferred_2018")
  tibble::tibble(
    label  = cl,
    n_anch = nrow(a), n_inf = nrow(i),
    anch_doy = round(mean(a$snow_free_doy), 1),
    inf_doy  = round(mean(i$snow_free_doy), 1),
    doy_shift = round(mean(i$snow_free_doy) - mean(a$snow_free_doy), 1),
    med_doy_diff = round(stats::median(abs(i$inference_doy_diff_days)), 1),
    p90_doy_diff = round(stats::quantile(abs(i$inference_doy_diff_days), 0.9), 1),
    med_hell     = round(stats::median(i$inference_hell_distance), 2),
    p90_hell     = round(stats::quantile(i$inference_hell_distance, 0.9), 2)
  )
}) |> dplyr::arrange(dplyr::desc(p90_doy_diff))

cat("=== Per-cluster inference quality ===\n")
cat("Sorted by p90 DOY mismatch (worst at top).\n\n")
print(per_cluster, n = Inf, width = Inf)

# --- Inferred-site confidence breakdown -------------------------------------
conf_tally <- asg |>
  dplyr::filter(source == "inferred_2018") |>
  dplyr::count(inference_confidence) |>
  dplyr::mutate(pct = round(100 * n / sum(n), 1))
cat("\n=== Inference confidence breakdown ===\n")
print(conf_tally)

# --- Worst-fit sites (by combined inference_distance) -----------------------
cat("\n=== Top 10 worst-fit inferred sites (highest combined distance) ===\n")
cat("These are 2018 sites whose composition (and/or DOY) does not closely\n")
cat("match any 2025 cluster centroid -- likely 2018-only ecological niches.\n\n")
worst <- asg |>
  dplyr::filter(source == "inferred_2018") |>
  dplyr::arrange(dplyr::desc(inference_distance)) |>
  dplyr::slice_head(n = 10) |>
  dplyr::select(site_number, final_label, snow_free_doy,
                inference_doy_diff_days, inference_hell_distance,
                inference_distance, inference_gap, inference_confidence)
print(worst, n = Inf, width = Inf)

cat("\nFilter on inference_confidence in QGIS to inspect spatially:\n")
cat("  high   ~half of inferred sites; close match in composition and DOY\n")
cat("  medium reasonable assignment, watch for ecological coherence\n")
cat("  low    poor composition match -- treat as low-confidence label\n")
