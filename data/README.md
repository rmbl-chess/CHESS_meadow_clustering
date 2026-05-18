# Data provenance

Three ESS-DIVE source dirs live under `data/raw/` (gitignored). Fetch from ESS-DIVE (DOIs TBD); the files below are working copies.

## `raw/ESS-DIVE-Vegetation-Field-2018/`

CHESS 2018 field campaign — fractional cover at AOP-aligned crown footprints.

| File | Purpose |
|---|---|
| `fractional_cover (1).csv` | Plot × species cover values. |
| `species_list (1).csv` | 2018 species list (canonical names for the campaign). |
| `metadata_column_key (1).csv` | Column dictionary for `fractional_cover`. |
| `CRBU2018_AOP_Crowns.geojson` | Crown polygons — the spatial join unit for 2018. |

## `raw/ESS-DIVE-Vegetation-Field-2025/`

CHESS 2025 field campaign — cover + site metadata, parallel structure to 2018.

| File | Purpose |
|---|---|
| `chess_meadow_cover_cleaned.csv` | Plot × species cover values. |
| `chess_meadow_site_cleaned.csv` | Site metadata (coords, attributes). |
| `chess_species_list_cleaned.csv` | 2025 species list. |
| `dd.csv`, `flmd.csv` | ESS-DIVE data dictionary + file-level metadata. |
| `CHESS_2025_Field_collected_vegetation_attributes.xml` | FGDC metadata. |

## `raw/ESS-DIVE-Spectra/`

NEON AOP spectra extracted at crown footprints for both campaign years. Spectral matrices and wavelengths are year-specific.

| File | Purpose |
|---|---|
| `site_extraction_spectra_2018 (1).csv` | Per-crown extracted spectra, 2018. |
| `site_extraction_spectra_2025 (1).csv` | Per-crown extracted spectra, 2025. |
| `wavelengths_2018.csv` | Band → wavelength table for 2018 spectra. |
| `wavelengths_2025.csv` | Band → wavelength table for 2025 spectra. |
| `CHESS_2025_crowns (1).geojson` | 2025 crown polygons (matches the 2018 GeoJSON in the veg-2018 dir). |
| `dd (2).csv`, `flmd (2).csv` | ESS-DIVE data dictionary + file-level metadata. |
| `CHESS_2025_Crown_polygons_and_extracted.xml` | FGDC metadata. |

## Subdirectories

- `raw/` — gitignored; ESS-DIVE working copies.
- `derived/` — gitignored; outputs from `code/` scripts (combined cover table, joined vegetation–spectrum dataset).
- `small_reference/` — committed; small canonical inputs (taxonomy crosswalk, AOI polygons).

## To do

- [ ] Confirm CRS of each crown GeoJSON (`EPSG:4326` vs `EPSG:32613`).
- [ ] Verify 2018 vs 2025 wavelength sets are on the same grid (same NEON AOP sensor / pipeline).
- [ ] Confirm SiteID scheme: 2018 vs 2025 IDs the same physical plots, or crosswalk needed?
- [ ] Add ESS-DIVE DOIs and any future updates.
