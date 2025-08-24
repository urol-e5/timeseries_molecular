# DNA Methylation Bedgraph Files

This directory contains species-level DNA methylation bedgraph files generated from filtered WGBS CpG count data for the timeseries molecular analysis project.

## Generated Files

### D-Apul_methylation.bedgraph
- **Species**: *Acropora pulchra* 
- **CpG sites**: 96,785
- **Samples**: 39 (across 4 timepoints)
- **Methylation range**: 0.00% - 93.59%
- **Source data**: `D-Apul/data/count-matrices-01/Apul-filtered-WGBS-CpG-counts.csv`

### F-Ptua_methylation.bedgraph  
- **Species**: *Pocillopora tuahiniensis*
- **CpG sites**: 137,292
- **Samples**: 32 (across multiple timepoints)
- **Methylation range**: 0.31% - 96.41%
- **Source data**: `F-Ptua/output/05-Ptua-bismark-CG/filtered-WGBS-CpG-counts.csv`

### E-Peve_methylation.bedgraph
- **Status**: Not generated - missing source data
- **Species**: *Porites evermanni*
- **Issue**: The required processed WGBS coverage files (*.CpG_report.merged_CpG_evidence.cov.gz) are not present in the expected locations
- **Next steps**: Need to complete the coverage2cytosine processing step or locate existing processed files

## File Format

The bedgraph files follow the standard UCSC bedGraph format:
```
track type=bedGraph name="[Species] DNA Methylation" description="Mean CpG methylation levels (0-100%)"
chromosome    start    end    methylation_percentage
```

Where:
- **chromosome**: Chromosome/scaffold identifier
- **start**: 0-based start position of CpG site
- **end**: 1-based end position of CpG site (start + 1)
- **methylation_percentage**: Mean methylation level across all samples (0.00-100.00%)

## Processing Details

### Data Processing Steps
1. **Input**: Filtered WGBS CpG count CSV files with sample-level methylation percentages
2. **CpG ID parsing**: Extract chromosome and position from CpG identifiers (format: `CpG_chromosome_position`)
3. **Mean calculation**: Calculate mean methylation percentage across all samples for each CpG site
4. **Sorting**: Sort by chromosome and genomic position
5. **Output**: Write in bedGraph format with track header

### Quality Control
- All CpG sites with valid genomic coordinates were included
- Sites with all NaN values across samples were excluded
- Mean methylation calculated excluding NaN values
- Coordinates converted to proper 0-based start, 1-based end format

## Usage

These bedgraph files can be:
- Loaded into genome browsers (UCSC Genome Browser, IGV, etc.)
- Used for comparative methylation analysis across species
- Integrated with other genomic annotation tracks
- Analyzed for methylation patterns and regulatory regions

## Generation Script

Files were generated using `generate_methylation_bedgraphs.py` which:
- Processes existing filtered WGBS CpG count CSV files
- Converts sample-level data to species-level averages
- Outputs properly formatted bedGraph files
- Provides summary statistics for each species

## Notes

- Methylation levels represent mean values across all available samples per species
- Different species have different numbers of samples and timepoints
- CpG coordinates are relative to species-specific genome assemblies
- Files are sorted for efficient genome browser loading

For questions about the data processing or file formats, refer to the generation script or the original WGBS analysis workflows in each species directory.