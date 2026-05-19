# meadow/_attic/

K-selection and merging-threshold exploration scripts that fed into the
final Architecture-B choice in `10_cluster_spectra.R`. They are run-once
methodology decisions, not part of the production pipeline. Nothing
downstream reads their outputs.

| Script | What it did |
|---|---|
| `06_cluster_composition.R` | Initial composition-only k-means at species + genus granularity. Compared LDA / RF separability before the spectra-first switch. |
| `07_spectral_separability.R` | LDA separability per K (composition-defined clusters). |
| `08_separability_rf.R` | Same as 07 but with Random Forest CV. Confirmed PC + env feature space was non-trivially better than composition alone. |
| `09_iterative_merge.R` | Iterative cluster-merging at LDA threshold = 0.4 — the precursor to spectra-first clustering. Established that aggregating across composition clusters until they were RF-separable was the right move. |

Re-run via `Rscript code/meadow/_attic/<file>` if you want to revisit
K-selection from scratch; outputs land in `data/derived/clusters_*.rds`
and `data/derived/separability_*.rds`.
