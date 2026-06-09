# 32_shrub_join_spectra.R — join AOP spectra to the canonical shrub records.
# Mirrors 04_join_spectra.R: per-pixel L2 normalize -> filter to sunlit
# (shade==1) -> mean across pixels per site -> one row per (site, Year).
#
# 2025-specific: spectra rows have `site_type` in {Meadow, Shrub, Tree}.
# We filter to site_type == "Shrub" so a shrub site_number that also
# appears in a co-located meadow crown isn't contaminated.
#
# Inputs:
#   data/derived/shrub_records_canonical.rds
#   data/derived/spectra_2018.rds, data/derived/spectra_2025.rds
# Outputs:
#   data/derived/shrub_veg_spectra.rds  list(joined, wavelengths)
#       joined: site_number, Year, canonical_binomial, canonical_genus,
#               sampling_area, n_pixels, rfl_band_1..426

suppressPackageStartupMessages({
  library(tidyverse)
})

records  <- readRDS("data/derived/shrub_records_canonical.rds")
sp_2018  <- readRDS("data/derived/spectra_2018.rds")
sp_2025  <- readRDS("data/derived/spectra_2025.rds")
sp_2026_path <- "data/derived/spectra_2026.rds"
sp_2026  <- if (file.exists(sp_2026_path)) readRDS(sp_2026_path) else NULL

# Wavelength sanity check (same as 04_join_spectra.R; 2025 is canonical).
wl_2018 <- sp_2018$wavelengths
wl_2025 <- sp_2025$wavelengths
stopifnot(nrow(wl_2018) == nrow(wl_2025))
max_center_diff_nm <- max(abs(wl_2018$center_wavelength -
                                wl_2025$center_wavelength)) * 1000
if (max_center_diff_nm > 20) {
  stop(sprintf("2018/2025 wavelengths differ by up to %.1f nm — resample.",
               max_center_diff_nm))
}
wavelengths <- wl_2025 |>
  dplyr::transmute(band_number,
                   center_wavelength_nm = center_wavelength * 1000,
                   fwhm_nm              = fwhm * 1000)

# --- Aggregate spectra to one mean L2-normalized vector per shrub site ---
agg_spectra <- function(df, year, keep_sunlit = TRUE) {
  rfl_cols <- grep("^rfl_band_", names(df), value = TRUE)
  if (keep_sunlit && "shade" %in% names(df)) {
    df <- dplyr::filter(df, shade == 1 | is.na(shade))
  }
  spec_mat <- as.matrix(df[, rfl_cols])
  norms <- sqrt(rowSums(spec_mat^2, na.rm = TRUE))
  norms[norms == 0] <- NA_real_
  df[, rfl_cols] <- spec_mat / norms
  df |>
    dplyr::group_by(site_number) |>
    dplyr::summarise(n_pixels = dplyr::n(),
                     dplyr::across(dplyr::all_of(rfl_cols),
                                   ~ mean(.x, na.rm = TRUE)),
                     .groups = "drop") |>
    dplyr::mutate(Year = year)
}

# 2025: restrict to site_type == "Shrub" before aggregation.
sp_2025_shrub <- dplyr::filter(sp_2025$spectra, site_type == "Shrub")
cat(sprintf("2025 spectra rows: %d total, %d at Shrub crowns\n",
            nrow(sp_2025$spectra), nrow(sp_2025_shrub)))

# 2018 spectra carries no site_type; the records-side filter (cover==100
# at a shrub-dominated genus) is the join key, so any pixels at a matching
# site_number are by definition from a shrub-dominated crown.
sp_2018_shrub <- sp_2018$spectra
cat(sprintf("2018 spectra rows: %d total (no site_type column)\n",
            nrow(sp_2018_shrub)))

spectra_2025 <- agg_spectra(sp_2025_shrub, 2025L)
spectra_2018 <- agg_spectra(sp_2018_shrub, 2018L)

# 2026: shrub crowns extracted from 2025 AOP (site_type == "Shrub"); no year
# correction (2025 basis). NULL when the Hub extraction hasn't run yet.
spectra_2026 <- if (!is.null(sp_2026)) {
  sp_2026_shrub <- dplyr::filter(sp_2026$spectra, site_type == "Shrub")
  cat(sprintf("2026 spectra rows: %d total, %d at Shrub crowns\n",
              nrow(sp_2026$spectra), nrow(sp_2026_shrub)))
  agg_spectra(sp_2026_shrub, 2026L)
} else NULL

# 2018 -> 2025 NDVI-stratified correction (mirror of meadow 04 logic;
# see code/joint/15_stratified_year_correction.R for fit + validation).
NDVI_BREAKS <- c(-Inf, 0.20, 0.40, 0.60, 0.80, Inf)
NDVI_LABELS <- c("ndvi_lt_0.20", "ndvi_0.20_0.40",
                 "ndvi_0.40_0.60", "ndvi_0.60_0.80", "ndvi_ge_0.80")
apply_year_correction <- function(df, correction_path, wavelengths) {
  if (!file.exists(correction_path)) return(df)
  cor <- readr::read_csv(correction_path, show_col_types = FALSE)
  bands <- sort(unique(cor$band_number))
  rfl_cols <- sprintf("rfl_band_%d", bands)
  stopifnot(all(rfl_cols %in% names(df)))
  band_at <- function(nm) which.min(abs(wavelengths$center_wavelength_nm - nm))
  bn_red <- wavelengths$band_number[band_at(660)]
  bn_nir <- wavelengths$band_number[band_at(860)]
  red <- df[[sprintf("rfl_band_%d", bn_red)]]
  nir <- df[[sprintf("rfl_band_%d", bn_nir)]]
  ndvi_site <- (nir - red) / (nir + red)
  bin_site  <- as.character(cut(ndvi_site, breaks = NDVI_BREAKS,
                                 labels = NDVI_LABELS, right = FALSE,
                                 include.lowest = TRUE))
  wide <- cor |>
    dplyr::select(band_number, ndvi_bin, delta) |>
    tidyr::pivot_wider(names_from = ndvi_bin, values_from = delta)
  band_order <- match(bands, wide$band_number)
  delta_mat <- matrix(NA_real_, nrow = nrow(df), ncol = length(bands))
  for (i in seq_len(nrow(df))) {
    b <- bin_site[i]
    if (is.na(b) || !(b %in% names(wide))) next
    delta_mat[i, ] <- wide[[b]][band_order]
  }
  df[, rfl_cols] <- as.matrix(df[, rfl_cols]) +
                    ifelse(is.na(delta_mat), 0, delta_mat)
  n_corr <- sum(!is.na(rowSums(delta_mat)))
  message(sprintf(
    "Applied NDVI-stratified 2018->2025 correction to %d / %d shrub sites.",
    n_corr, nrow(df)
  ))
  df
}
spectra_2018 <- apply_year_correction(
  spectra_2018,
  "data/small_reference/year_effect_correction_2018_to_2025_by_ndvi.csv",
  wavelengths
)

spectra_combined <- dplyr::bind_rows(spectra_2018, spectra_2025, spectra_2026)

# --- Join to records ------------------------------------------------------
joined <- records |>
  dplyr::select(site_number, Year, sampling_area,
                canonical_binomial, canonical_genus,
                raw_binomial, vegetation_height_cm) |>
  dplyr::inner_join(spectra_combined, by = c("site_number", "Year"))

cat(sprintf(
  "\nJoined shrub records to spectra: %d sites (%d 2018, %d 2025, %d 2026)\n",
  nrow(joined), sum(joined$Year == 2018L), sum(joined$Year == 2025L),
  sum(joined$Year == 2026L)
))

# Report records that fell out of the join (no matching pixels at the site).
unjoined <- records |>
  dplyr::anti_join(joined, by = c("site_number", "Year"))
cat(sprintf("Records with no matching spectra: %d\n", nrow(unjoined)))
if (nrow(unjoined) > 0) {
  cat("Top unjoined classes:\n")
  print(unjoined |> dplyr::count(Year, canonical_binomial, name = "n") |>
          dplyr::arrange(dplyr::desc(n)) |> utils::head(15))
}

# --- Class counts on the joined set --------------------------------------
cat("\n=== Joined-class N per canonical_binomial per year ===\n")
print(joined |>
        dplyr::count(canonical_binomial, Year) |>
        tidyr::pivot_wider(names_from = Year, values_from = n,
                           names_prefix = "n_", values_fill = 0L) |>
        dplyr::mutate(n_total = rowSums(dplyr::across(dplyr::starts_with("n_")))) |>
        dplyr::arrange(dplyr::desc(n_total)) |>
        as.data.frame())

# --- Persist --------------------------------------------------------------
saveRDS(list(joined = joined, wavelengths = wavelengths,
             unjoined = unjoined),
        "data/derived/shrub_veg_spectra.rds")
cat("\nWrote data/derived/shrub_veg_spectra.rds\n")
