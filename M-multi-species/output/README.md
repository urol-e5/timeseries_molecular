`timeseries_molecular/M-multi-species/output/`

## Directory overview

This folder contains the output folders from scripts in the `M-multi-species` analysis pipeline. Each subdirectory should correspond to a script in the `scripts/` directory with the same name.

## Barnacle Analysis Reports

For comprehensive documentation of all Barnacle sparse tensor decomposition analyses:

- **[BARNACLE_ANALYSIS_REPORT.md](BARNACLE_ANALYSIS_REPORT.md)** - Detailed report covering all 6 major analysis streams, systematic convergence testing (94 parameter combinations), multi-rank comparison, findings, and recommendations
- **[BARNACLE_QUICK_REFERENCE.md](BARNACLE_QUICK_REFERENCE.md)** - Quick reference guide with best parameter sets, command examples, and result locations

### Barnacle Analysis Directories

- `13.00-multiomics-barnacle/` - Initial baseline analysis with normalized data
- `14-barnacle/` - Refined approach (Rank 5)
- `14.1-barnacle/` - Single run analysis
- `14.1-barnacle-convergence-tests/` - Systematic testing of 94 parameter combinations
- `14.5-barnacle/` - Advanced configuration testing
- `15-barnacle/` - Multi-rank comparison (Ranks 3, 8, 10, 12)
- `15.5-barnacle/` - Additional rank variations

### Key Finding

None of the 94 parameter combinations tested achieved convergence, but the best-performing parameter set was identified:
- **Best for Rank 5**: λ_gene=0.01, λ_sample=0.1, λ_time=0.05 → Final loss ~840,915
- **Best overall**: Rank 12 with standard parameters → Final loss 587,541