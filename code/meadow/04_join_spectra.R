# 04_join_spectra.R — join the combined cover table to extracted spectra by
# (site_number, Year) to produce the vegetation–spectrum dataset for clustering.
#
# Notes on the spectra inputs (see 01_load.R header):
#   site_number       integer site ID (matches cover_combined)
#   fid               crown polygon identifier (one site can have multiple)
#   row, col          pixel coords within the crown extraction
#   x_utm, y_utm      UTM coords (EPSG:32613)
#   shade             shade flag
#   rfl_band_1..426   reflectance per band
#   unc_band_1..426   per-band uncertainty (carry but don't cluster on)
#   2025-only: domain, sampling_area, site_type
#
# A crown / site has many pixels. For a first pass we aggregate to one mean
# reflectance vector per site_number (filtering out shaded pixels). Tweak this
# in clustering once the analysis approach is chosen.
#
# Inputs:
#   data/derived/cover_combined.rds
#   data/derived/spectra_2018.rds, data/derived/spectra_2025.rds
# Outputs:
#   data/derived/veg_spectra.rds   list with:
#       joined        (site_number, Year, <Spp>_cover, ..., rfl_band_1..K)
#       wavelengths   (band_number, center_wavelength_nm, fwhm_nm)

library(tidyverse)

cover_combined <- readRDS("data/derived/cover_combined.rds")
sp_2018        <- readRDS("data/derived/spectra_2018.rds")
sp_2025        <- readRDS("data/derived/spectra_2025.rds")
# 2026 supplemental spectra (extracted from 2025 AOP); present only after the
# Hub extraction has run. Optional so the pipeline still runs without it.
sp_2026_path   <- "data/derived/spectra_2026.rds"
sp_2026        <- if (file.exists(sp_2026_path)) readRDS(sp_2026_path) else NULL

# Wavelength grids: 2018 and 2025 are nearly identical (426 bands, 0.384-2.510
# μm) but drift up to 10 nm in bands 81-135 (NIR plateau) — likely a NEON
# processing change. 2025 is canonical (newer processing). Offset is within
# band spacing (~5 nm) so we use band index directly for now; if clustering
# shows artifacts, resample 2018 onto the 2025 grid with approx().
wl_2018 <- sp_2018$wavelengths
wl_2025 <- sp_2025$wavelengths
stopifnot(nrow(wl_2018) == nrow(wl_2025))
max_center_diff_nm <- max(abs(wl_2018$center_wavelength - wl_2025$center_wavelength)) * 1000
if (max_center_diff_nm > 20) {
  stop(sprintf("2018/2025 center wavelengths differ by up to %.1f nm — resample required.",
               max_center_diff_nm))
}
message(sprintf("Wavelength grids: max center offset %.1f nm; using 2025 as canonical.",
                max_center_diff_nm))
wavelengths <- wl_2025 |>
  dplyr::transmute(
    band_number,
    center_wavelength_nm = center_wavelength * 1000,
    fwhm_nm              = fwhm * 1000
  )

# Aggregate spectra to one mean *brightness-normalized* reflectance vector per
# (site_number, Year). Order matters: averaging raw spectra mixes pixels with
# different illumination, blurring the shape; brightness-normalize each pixel
# (L2 unit-vector) first, then mean.
#
# `shade == 1` flags sunlit pixels (the column polarity is inverted relative
# to its name; verified empirically: 2025 meadow pixels are 1935:14 sunlit:
# shaded; 2018 (all meadow) is 13524:2624 same direction). Default keeps
# sunlit pixels; flip `keep_sunlit = FALSE` to retain all.
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
    dplyr::summarise(dplyr::across(dplyr::all_of(rfl_cols),
                                   ~ mean(.x, na.rm = TRUE)),
                     .groups = "drop") |>
    dplyr::mutate(Year = year)
}

# 2018 -> 2025 NDVI-stratified correction. Per-band delta fit separately
# for 5 NDVI bins on the 2115 paired CRBU points from script 13's
# extraction (see code/joint/15_stratified_year_correction.R for the
# fit + held-out validation: mean |Cohen's d| dropped from 0.36 -> 0.06
# for vegetated points). Each per-site spectrum gets its own bin by
# computing NDVI on the L2-normalized + mean-per-site reflectance, then
# the bin's per-band delta is added.
NDVI_BREAKS <- c(-Inf, 0.20, 0.40, 0.60, 0.80, Inf)
NDVI_LABELS <- c("ndvi_lt_0.20", "ndvi_0.20_0.40",
                 "ndvi_0.40_0.60", "ndvi_0.60_0.80", "ndvi_ge_0.80")

apply_year_correction <- function(df, year, correction_path, wavelengths) {
  if (year != 2018L) return(df)
  if (!file.exists(correction_path)) {
    message(sprintf("No correction CSV at %s; skipping year correction.",
                    correction_path))
    return(df)
  }
  cor <- readr::read_csv(correction_path, show_col_types = FALSE)
  bands <- sort(unique(cor$band_number))
  rfl_cols <- sprintf("rfl_band_%d", bands)
  stopifnot(all(rfl_cols %in% names(df)))

  # NDVI per site on the L2-normalized averaged spectrum (NIR ~860,
  # RED ~660). wavelengths is the per-band centers in nm.
  band_at <- function(nm) which.min(abs(wavelengths$center_wavelength_nm - nm))
  bn_red <- wavelengths$band_number[band_at(660)]
  bn_nir <- wavelengths$band_number[band_at(860)]
  red <- df[[sprintf("rfl_band_%d", bn_red)]]
  nir <- df[[sprintf("rfl_band_%d", bn_nir)]]
  ndvi_site <- (nir - red) / (nir + red)
  bin_site  <- as.character(cut(ndvi_site, breaks = NDVI_BREAKS,
                                 labels = NDVI_LABELS, right = FALSE,
                                 include.lowest = TRUE))

  # Wide: rows = bands, columns = ndvi bins
  wide <- cor |>
    dplyr::select(band_number, ndvi_bin, delta) |>
    tidyr::pivot_wider(names_from = ndvi_bin, values_from = delta)
  band_order <- match(bands, wide$band_number)

  delta_mat <- matrix(NA_real_, nrow = nrow(df), ncol = length(bands))
  n_per_bin <- integer(length(NDVI_LABELS))
  names(n_per_bin) <- NDVI_LABELS
  for (i in seq_len(nrow(df))) {
    b <- bin_site[i]
    if (is.na(b) || !(b %in% names(wide))) next
    delta_mat[i, ] <- wide[[b]][band_order]
    n_per_bin[b] <- n_per_bin[b] + 1L
  }
  df[, rfl_cols] <- as.matrix(df[, rfl_cols]) +
                    ifelse(is.na(delta_mat), 0, delta_mat)
  n_corr <- sum(!is.na(rowSums(delta_mat)))
  message(sprintf(
    "Applied NDVI-stratified 2018->2025 correction to %d / %d sites; per-bin n:",
    n_corr, nrow(df)
  ))
  message(paste(sprintf("  %s: %d", names(n_per_bin), n_per_bin),
                collapse = "\n"))
  df
}

correction_path <- "data/small_reference/year_effect_correction_2018_to_2025_by_ndvi.csv"
spectra_combined <- dplyr::bind_rows(
  apply_year_correction(agg_spectra(sp_2018$spectra, 2018L),
                         2018L, correction_path, wavelengths),
  agg_spectra(sp_2025$spectra, 2025L),
  # 2026 extracted from 2025 AOP -> no year correction (same basis as 2025).
  if (!is.null(sp_2026)) agg_spectra(sp_2026$spectra, 2026L) else NULL
)

# Inner join on (site_number, Year) so only sites with both cover and spectra
# survive into the clustering input.
joined <- dplyr::inner_join(
  cover_combined, spectra_combined,
  by = c("site_number", "Year")
)

message(sprintf("veg_spectra: %d sites (%d 2018, %d 2025, %d 2026), %d cover cols, %d bands",
                nrow(joined),
                sum(joined$Year == 2018L),
                sum(joined$Year == 2025L),
                sum(joined$Year == 2026L),
                sum(grepl("_cover$", names(joined))),
                sum(grepl("^rfl_band_", names(joined)))))

saveRDS(list(joined = joined, wavelengths = wavelengths),
        "data/derived/veg_spectra.rds")
