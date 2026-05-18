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

spectra_combined <- dplyr::bind_rows(
  agg_spectra(sp_2018$spectra, 2018L),
  agg_spectra(sp_2025$spectra, 2025L)
)

# Inner join on (site_number, Year) so only sites with both cover and spectra
# survive into the clustering input.
joined <- dplyr::inner_join(
  cover_combined, spectra_combined,
  by = c("site_number", "Year")
)

message(sprintf("veg_spectra: %d sites (%d 2018, %d 2025), %d cover cols, %d bands",
                nrow(joined),
                sum(joined$Year == 2018L),
                sum(joined$Year == 2025L),
                sum(grepl("_cover$", names(joined))),
                sum(grepl("^rfl_band_", names(joined)))))

saveRDS(list(joined = joined, wavelengths = wavelengths),
        "data/derived/veg_spectra.rds")
