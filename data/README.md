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

## `raw/Vegmap_2023/`

2023 UER vegmap campaign (collaborator) — 238 ~1 m² plots, all within CRBU, sampled Jun–Aug 2023. Working copies from
`Google Drive: BreckheimerLab2025/Projects/CHESS/Projects/Vegmap_2023/data`. Standardized into the pipeline by `code/meadow/00_prep_2023.R`; spectra come from the **2025 AOP** (no 2023 CRBU flight — assumes 2023→2025 composition stability).

| File | Purpose |
|---|---|
| `vegmap_2023_cover.csv` | Plot × species cover (6-letter codes; totals = 100%). 2026-08-19 update added her authoritative `site_number` column (2500–2736). |
| `vegmap_2023_wgs_utm_pixel_select.geojson` | AOP-pixel-snapped plot polygons, EPSG:32613; `pix_sel` = reviewed. |
| `uer_vegmap_locations_2023.csv` | Raw GNSS corner points (provenance only). |

The code → canonical crosswalk is committed at `small_reference/taxonomy_crosswalk_2023.csv` (collaborator's 3-campaign file, dated 2026-08-06); the adopted site-number mapping is mirrored at `small_reference/site_numbers_2023.csv`. The 9 ambiguous plots (cover-only SH-01..04, polygon-only JF-13..17) are dropped by `00_prep_2023.R` — 233 clean sites enter the pipeline.

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
  - `class_ecosystem_crosswalk.csv` — the NatureServe **curation sheet** for the review round: one row per meadow class, pre-filled with the top-1 draft community; reviewers fill `decision` (accept / replace / reject), the `curated_*` columns, and `notes`. Curated results get wired into `19`/`06`/`10` once the round closes.

## To do

- [ ] Confirm CRS of each crown GeoJSON (`EPSG:4326` vs `EPSG:32613`).
- [ ] Verify 2018 vs 2025 wavelength sets are on the same grid (same NEON AOP sensor / pipeline).
- [ ] Confirm SiteID scheme: 2018 vs 2025 IDs the same physical plots, or crosswalk needed?
- [ ] Add ESS-DIVE DOIs and any future updates.
