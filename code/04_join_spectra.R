# 04_join_spectra.R — join the combined cover table to extracted spectra by
# (SiteID, Year) to produce the single vegetation–spectrum dataset that feeds
# clustering.
#
# Inputs:
#   data/derived/cover_combined.rds
#   data/derived/spectra_2018.rds, data/derived/spectra_2025.rds
# Outputs:
#   data/derived/veg_spectra.rds   list with:
#       cover         (SiteID, Year, <Spp>_cover, ...)
#       spectra       (SiteID, Year, b1, b2, ...) — possibly resampled to a
#                                                   common wavelength grid
#       wavelengths   (band, wavelength_nm)        — the chosen common grid

library(tidyverse)

cover_combined <- readRDS("data/derived/cover_combined.rds")
sp_2018        <- readRDS("data/derived/spectra_2018.rds")
sp_2025        <- readRDS("data/derived/spectra_2025.rds")

# TODO: identify the SiteID column in each spectra table; align with the
#       cover-table SiteID scheme.

# Decide whether wavelengths can be stacked directly or need resampling:
wl_equal <- isTRUE(all.equal(sp_2018$wavelengths, sp_2025$wavelengths))
if (!wl_equal) {
  # TODO: resample one (or both) to a common wavelength grid before stacking.
  stop("2018 and 2025 wavelength tables differ — implement resampling.")
}

spectra_combined <- dplyr::bind_rows(
  sp_2018$spectra |> dplyr::mutate(Year = 2018L),
  sp_2025$spectra |> dplyr::mutate(Year = 2025L)
)

veg_spectra <- list(
  cover       = cover_combined,
  spectra     = spectra_combined,
  wavelengths = sp_2018$wavelengths
)

# Inner join on (SiteID, Year) ensures only sites with both cover and spectra
# survive — that's the clustering input.
joined <- dplyr::inner_join(
  veg_spectra$cover,
  veg_spectra$spectra,
  by = c("SiteID", "Year")
)
saveRDS(list(joined = joined, wavelengths = veg_spectra$wavelengths),
        "data/derived/veg_spectra.rds")
