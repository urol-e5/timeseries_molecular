# Compare Multiway Interactions Across 3 Species

## Overview

This script analyzes and compares multiway interactions across three species by:
1. Downloading or loading CSV files containing multi-way interaction data
2. Counting occurrences of "CpG", "lncRNA", and "miRNA" features in each row
3. Generating summary statistics
4. Creating visualizations comparing the three species

## Script Location

`M-multi-species/scripts/30-compare-multiway-interactions.py`

## Output Location

All outputs are saved to: `M-multi-species/output/30-compare-multiway-interactions/`

## Data Sources

The script processes data from three species with the following URLs:
- Species 1: `https://gannet.fish.washington.edu/v1_web/owlshell/bu-github/ConTra/output/context_dependent_analysis_20251108_151255/tables/multi_way_interactions.csv`
- Species 2: `https://gannet.fish.washington.edu/v1_web/owlshell/bu-github/ConTra/output/context_dependent_analysis_20251108_144034/tables/multi_way_interactions.csv`
- Species 3: `https://gannet.fish.washington.edu/v1_web/owlshell/bu-github/ConTra/output/context_dependent_analysis_20251108_140602/tables/multi_way_interactions.csv`

## Usage

### Basic Usage

```bash
python M-multi-species/scripts/30-compare-multiway-interactions.py
```

### Using Local Files

If the URLs are not accessible, you can manually download the CSV files and place them in the output directory with these names:
- `species_1_multi_way_interactions.csv`
- `species_2_multi_way_interactions.csv`
- `species_3_multi_way_interactions.csv`

The script will automatically use existing local files if available.

## Outputs

### Data Files

1. **Raw Downloaded Data**:
   - `species_1_multi_way_interactions.csv`
   - `species_2_multi_way_interactions.csv`
   - `species_3_multi_way_interactions.csv`

2. **Annotated Data** (with count columns):
   - `species_1_annotated.csv`
   - `species_2_annotated.csv`
   - `species_3_annotated.csv`

3. **Summary Statistics**:
   - `summary_statistics.csv` - Contains:
     - Total rows per species
     - Total occurrences of each feature type (CpG, lncRNA, miRNA)
     - Mean occurrences per row
     - Number of rows with at least one occurrence of each feature

### Visualizations

1. **`total_occurrences_by_species.png`**
   - Bar chart showing total count of each feature type across all rows for each species

2. **`mean_occurrences_by_species.png`**
   - Bar chart showing the mean number of occurrences per row for each feature type

3. **`percentage_rows_with_occurrence.png`**
   - Bar chart showing the percentage of rows containing at least one occurrence of each feature

4. **`count_heatmaps.png`**
   - Heatmaps showing the distribution of feature counts across the first 50 rows for each species

## Feature Detection

The script performs case-insensitive detection of the following features across all columns in each row:
- **CpG**: Searches for "cpg" (case-insensitive)
- **lncRNA**: Searches for "lncrna" (case-insensitive)
- **miRNA**: Searches for "mirna" (case-insensitive)

## Dependencies

- pandas
- numpy
- matplotlib
- seaborn

## Example Output

Summary statistics example:
```
species_1:
  Total rows: 100
  CpG - Total: 150, Mean: 1.50, Rows with occurrence: 85
  lncRNA - Total: 120, Mean: 1.20, Rows with occurrence: 75
  miRNA - Total: 95, Mean: 0.95, Rows with occurrence: 60
```

## Notes

- The script handles network errors gracefully and will use local files if downloads fail
- All visualizations are saved as high-resolution PNG files (300 dpi)
- The counting is case-insensitive to ensure all variants are detected
