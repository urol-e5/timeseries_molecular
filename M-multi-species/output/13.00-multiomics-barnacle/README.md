`urol-e5/timeseries_molecular/M-multi-species/output/13.00-multiomics-barnacle`

## Directory overview

This folder contains outputs produced by `scripts/13.00-multiomics-barnacle.Rmd` which prepares ortholog-mapped expression matrices for three species, normalizes the data, builds a 3D tensor for multiomics analysis, and runs the Barnacle sparse CP decomposition.

---

### Files (alphabetical)

- `apul_normalized_expression.csv` — Normalized Acropora pulchra expression matrix (sctransform::vst or log2(CPM+1) fallback).
- `apul_ortholog_expression.csv` — Acropora pulchra expression matrix aligned to ortholog `group_id`.
- `multiomics_tensor.npy` — NumPy binary file storing the 3D tensor used for decomposition. Tensor shape: (genes, combined_species_samples, timepoints).
- `normalized_merged_ortholog_expression.csv` — Merged ortholog expression matrix combining all three species with sample columns in format `<SPECIES_PREFIX>-<SAMPLE_ID>-<TP#>` (e.g., `ACR-139-TP1`).
- `peve_normalized_expression.csv` — Normalized Pocillopora verrucosa expression matrix (sctransform::vst or log2(CPM+1) fallback).
- `peve_ortholog_expression.csv` — Pocillopora verrucosa expression matrix aligned to ortholog `group_id`.
- `ptua_normalized_expression.csv` — Normalized Pocillopora meandrina expression matrix (sctransform::vst or log2(CPM+1) fallback).
- `ptua_ortholog_expression.csv` — Pocillopora meandrina expression matrix aligned to ortholog `group_id`.

