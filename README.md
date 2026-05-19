# CHESS_meadow_clustering

Cluster co-occurring plant species and NEON AOP spectra into training samples for a CHESS meadow map classifier across the East River / Upper Taylor / Almont basins.

## What this is

The CHESS map classifier needs training samples that are (a) spectrally separable in 1 m NEON AOP imagery and (b) interpretable as recognizable plant species compositions. This repo:

1. Reconciles 2018 + 2025 CHESS field vegetation campaigns into one harmonized cover table.
2. Joins per-crown NEON AOP spectra and clusters 2025 plots into spectrally distinct community types.
3. Infers 2018 labels by composition + phenology similarity (2018 atmospheric correction was unreliable, so 2018 sites are excluded from the clustering itself).
4. Selects a stratified set of meadow target pixels across the 2025 AOP mosaic.
5. Provides Python tooling to extract classifier features and full-tile PC maps on a cloud server.

**Current state:** 31 final classes (26 spectral + 5 monotypic-species overrides), labeled with uniform "Life zone + taxa + physiognomy" narratives. Target pixel set is 6,000 (stratified by snow-free DOY × domain). The classifier itself is fit Python-side from the exported features.

## Pipeline

R scripts under `code/` are numbered in run order. Targets (`renv` env assumed restored):

| Step | Script | Output |
|---|---|---|
| Load 2018/2025 field data | `01_load.R` | cover + sites + species lists, in-memory |
| Reconcile taxonomy | `02_reconcile_taxonomy.R` | `data/derived/taxonomy_crosswalk.csv` |
| Combine cover | `03_combine_cover.R` | `data/derived/cover_combined.rds` |
| Join AOP spectra | `04_join_spectra.R` | `data/derived/veg_spectra.rds` |
| Preprocess features | `05_preprocess_features.R` | brightness-normed spectra + Hellinger cover |
| EDA: composition clusters | `06_cluster_composition.R` | initial composition-only k-means |
| RF spectral separability | `07b_separability_rf.R` | recall + tier classification |
| Iterative merging | `08_iterative_merge.R` | merged class set, threshold = 0.4 |
| **Architecture B: spectra-first cluster** | `09_cluster_spectra.R` | Ward / PCs 2–12 + snow-free DOY, z-scaled |
| Composition sub-clusters | `10_subcluster_composition.R` | within-spectral subclusters by composition |
| Diagnostics | `11_visualizations.R`, `13_diagnose_env.R`, `15_diagnose_brightness.R` | PNGs in `output/` |
| Environment join | `12_extract_environment.R` | snow-free DOY, elevation, etc. |
| Sweep over K | `14_sweep_kspec.R` | RF accuracy vs spectral K |
| Training-sample export | `16_export_training_samples.R` | `training_samples_*.{csv,gpkg}` |
| Year-effect check | `17_year_effect_pcs.R` | confirms 2018 spectra drift; clustering uses 2025 only |
| Class descriptions | `18_label_descriptions.R` | IndVal + abundant species tables |
| Narrative drafts | `19_label_narratives.R` | `data/small_reference/label_community_names.csv` |
| Inference QC | `20_inference_quality.R` | confidence tier per 2018 site |
| Size distribution | `21_size_distribution.R` | class sizes for review |
| **Target pixel selection** | `22_target_pixels.R` | `target_pixels.{csv,gpkg}` |
| Classifier-feature export | `24_export_classifier_model.R` | `aop_classifier_*.{csv,json}` (R → Python handoff) |

Python tools under `code/python/` are designed to run on a cloud server with the AOP S3 data:

| Script | Purpose |
|---|---|
| `extract_aop_features.py` | Extract 20 PCs + 6 narrow-band indices at each target pixel (3 m, L2 norm → 3×3 mean → water mask → PCA project). |
| `generate_aop_pc_maps.py` | Write 20-band PC COGs for every 2025 AOP tile + per-domain VRT mosaics. Same preprocessing as the feature extractor, so PC maps and target-pixel features are on the same scale. |

## Key choices

- **Architecture B (spectra-first).** Spectral clusters are the primary classes; composition is used post-hoc to interpret and sub-cluster. Reasoning: a clustering that's spectrally inseparable can't be mapped no matter how ecologically interpretable.
- **2025-only training.** 2018 AOP spectra show atmospheric-correction drift (see `17_year_effect_pcs.R`). 2018 sites are assigned labels by composition (Hellinger) + snow-free DOY similarity to 2025 plots and tagged with a confidence tier (high / medium / low).
- **PCs 2–12 + snow-free DOY**, all z-scaled, Ward linkage at K = 15 (variant G). PC1 dropped because it loaded heavily on overall brightness and dominated the clustering.
- **Monotypic-species overrides.** Five species (e.g., *Veratrum californicum*, *Dasiphora fruticosa*) are split off as their own classes wherever cover ≥ 70 %; they are spectrally distinctive enough to deserve dedicated classes even at small-N.
- **Target pixel filter.** Each 3 m pixel must be ≥ 80 % meadow (within-pixel, from R3D018 1 m landcover) **and** have ≥ 6 of 9 neighbors that also pass — so target pixels sit interior to ≥ 9 m × 9 m meadow patches and are unlikely to be tree-edge contaminated.
- **CRS-explicit.** SDP / AOP layers are `EPSG:32613` (UTM Zone 13N); field-data GeoJSONs are reprojected from WGS84 before any spatial join.

## Outputs

Most outputs are gitignored. Canonical artifacts are in `data/small_reference/` (committed) and `data/derived/` (rebuild from scripts):

| Path | What |
|---|---|
| `data/small_reference/label_community_names.csv` | Hand-tuned narratives for the 31 classes |
| `data/derived/training_samples_sites.csv` | One row per plot × year, with `final_label`, confidence tier, env covariates |
| `data/derived/training_samples_crowns.gpkg` | Crown polygons for QGIS review |
| `data/derived/target_pixels.{csv,gpkg}` | 6,000 stratified meadow pixels, with domain + parent AOP tile filename |
| `data/derived/aop_classifier_pca.csv` | PCA loadings + center for the 348 retained bands × 20 PCs |
| `data/derived/aop_classifier_indices.csv` | Narrow-band index formulas (NDVI / NDWI / PRI / red-edge slope / CAI / NDLI) |
| `data/derived/aop_classifier_meta.json` | Preprocessing config (water-band ranges, no-data sentinel) |

After running the Python tools on a cloud server:

| Path | What |
|---|---|
| `data/derived/target_pixel_features.parquet` | Per-target-pixel features ready for classifier inference |
| `data/derived/aop_pc_maps/{DOMAIN}/pc_{e}_{n}.tif` | 20-band PC COGs at 3 m, one per AOP tile |
| `data/derived/aop_pc_maps/{DOMAIN}_pc_mosaic.vrt` | Per-domain mosaic for QGIS |

## Layout

```
code/                  Numbered R scripts (01 → 24) + python/ + examples/
data/raw/              ESS-DIVE working copies (gitignored, fetch from ESS-DIVE)
data/derived/          Pipeline outputs (gitignored)
data/small_reference/  Small canonical inputs (committed)
docs/                  Methods notes, archived figures
output/                Figures + tables (gitignored)
DESCRIPTION            Drives renv only — this is NOT an R package
renv.lock              R dependency lockfile
```

See `data/README.md` for ESS-DIVE source dirs and file-level provenance.

## Reproducing

```bash
git clone https://github.com/rmbl-chess/CHESS_meadow_clustering
cd CHESS_meadow_clustering
R -e 'renv::restore()'
# Fetch the three ESS-DIVE source dirs into data/raw/ (see data/README.md)
for f in code/[0-9]*.R; do Rscript "$f"; done
```

Python tooling for the AOP extraction stage:

```bash
# Local dry run
python code/python/extract_aop_features.py --max-pixels 50 --output /tmp/test.parquet --workers 1
python code/python/generate_aop_pc_maps.py --tile-limit 1 --output-dir /tmp/pc_test

# Full cloud run
python code/python/extract_aop_features.py --workers 16 --memory-limit 8GB
python code/python/generate_aop_pc_maps.py
```

## Status / open questions

- Classifier itself (RF / GBM in Python) is the next step — features and target pixels are ready.
- ESS-DIVE DOIs and SiteID crosswalk between 2018 and 2025 plots still need to be locked down in `data/README.md`.
- 2018 wavelength set is close to 2025 (max drift 2 nm) but resampling is implicit via nearest-index match.

## License

Internal CHESS-team use; no LICENSE file.
