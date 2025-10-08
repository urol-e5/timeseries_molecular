`urol-e5/timeseries_molecular/M-multi-species/output/13.00-multiomics-barnacle/barnacle_factors`

## Directory overview

This directory is produced after running the Barnacle decomposition and contains the factor matrices, metadata and visualizations used for downstream interpretation.

---

### Files (alphabetical)

- `component_weights.csv` — component weights / importance values used to rank extracted components.
- `figures/` — directory containing visualization PNGs (component weights, timepoint loadings, sample heatmaps, PCA plots, top ortholog barplots, etc.).
- `gene_factors.csv` — gene loadings per component (rows = ortholog `group_id`, columns = components).
- `metadata.json` — JSON summary of analysis parameters and tensor metadata (tensor shape, rank, lambdas, convergence status, timepoint order, etc.).
- `sample_factors.csv` — combined sample (species_sample) loadings per component, includes `Species` and `Sample_ID` metadata columns for interpretation.
- `sample_mapping.csv` — mapping from combined sample index to `species` and `sample_id` (useful to link sample_factors rows back to species-specific samples).
- `time_factors.csv` — timepoint loadings per component (rows correspond to ordered timepoints).