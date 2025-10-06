# GFF File Evaluation

## Overview

This directory contains code to evaluate all GFF and GFF3 files in the timeseries_molecular repository. The analysis identifies file properties, validates GFF3 compliance, and provides detailed statistics for each file.

## Files

- `evaluate-gff-files.sh` - Main bash script that analyzes all GFF/GFF3 files in the repository

## Outputs

The script generates four output files in `M-multi-species/output/`:

1. **gff-file-evaluation-summary.txt** - Summary table with key properties of all GFF files
2. **gff-file-detailed-analysis.txt** - Detailed analysis of each file including:
   - File size and line counts
   - Feature types present
   - GFF3 compliance status
   - Sample attributes
   - Source information
3. **gff-feature-type-matrix.txt** - Cross-file feature matrix showing:
   - All unique feature types found across all files
   - Count of each feature type in each GFF file (matrix format)
   - Summary statistics for each feature type
4. **gff-feature-type-by-file.txt** - Feature counts organized by file:
   - Each GFF file listed with all its feature types and counts
   - Total feature count per file
   - Easy to read, file-by-file format

## Usage

To run the analysis:

```bash
cd /home/runner/work/timeseries_molecular/timeseries_molecular
bash M-multi-species/code/evaluate-gff-files.sh
```

The script will:
1. Find all `.gff` and `.gff3` files in the repository
2. Analyze each file for:
   - Total size and line count
   - Number of feature lines (non-comment, non-empty)
   - Field count validation (GFF3 requires exactly 9 tab-separated fields)
   - Attribute count in the 9th column
   - Feature types present in the file
   - Source information
   - Presence of GFF3 header (`##gff-version 3`)
3. Generate summary statistics including counts of compliant vs. non-compliant files

## GFF3 Specification

A valid GFF3 file should have:
- Exactly 9 tab-separated fields per line (excluding comments)
- Fields: seqid, source, type, start, end, score, strand, phase, attributes
- Optional header line: `##gff-version 3`
- Attributes in field 9 should be semicolon-separated key=value pairs
- Common attributes include: ID, Name, Parent, etc.

## Results Summary

As of the last analysis (2025-10-06):
- **Total files analyzed:** 30
- **GFF3 compliant (9 fields):** 29
- **Non-compliant:** 1
- **Unique feature types found:** 20

### Feature Type Distribution

The analysis identified 20 unique feature types across all GFF files:
- **CDS** - Present in 7 files, total: 1,484,549 features
- **exon** - Present in 5 files, total: 1,082,562 features
- **mRNA** - Present in 8 files, total: 687,451 features
- **gene** - Present in 6 files, total: 233,200 features
- **Unknown_sRNA_locus** - Present in 3 files, total: 41,472 features
- **3prime_UTR** - Present in 2 files, total: 73,806 features
- **five_prime_UTR** - Present in 2 files, total: 73,605 features
- **UTR** - Present in 2 files, total: 30,196 features
- **tRNA** - Present in 2 files, total: 15,848 features
- **miRNA_binding** - Present in 3 files, total: 8,791 features
- **sRNA** - Present in 3 files, total: 6,031 features
- **nucleotide_motif** - Present in 1 file, total: 195,480 features
- And 8 additional feature types (MIRNA_hairpin, mature_miRNA, siRNA loci, etc.)

### Non-compliant File

The file `./D-Apul/output/06-Apul-miRNA-mRNA-RNAhybrid/Apulcra-genome-mRNA_only_MAX1000.gff` contains some lines with 10 fields instead of 9, making it non-compliant with GFF3 specifications. A corrected version exists as `Apulcra-genome-mRNA_only_MAX1000_formatted.gff`.

## File Categories

The repository contains GFF files for multiple species and purposes:

### Species-specific genome annotations:
- **D-Apul** (Acropora pulchra): Main genome GFF and various derived files
- **E-Peve** (Porites evermanni): Genome annotations and validated GFF3
- **F-Ptua** (Pocillopora tuahiniensis/meandrina): Gene annotations and validated GFF3

### Analysis-derived GFF files:
- **UTR annotations** (3' and 5' UTRs)
- **miRNA-mRNA interactions** (RNAhybrid results)
- **sRNA discovery** (ShortStack results)
- **CpG annotations**
- **Gene-only subsets**

## Attribute Counts by File Type

The analysis shows varying attribute complexity:
- **1 attribute**: Simple gene/feature annotations (ID only)
- **2 attributes**: CDS/exon files (ID + Parent)
- **3 attributes**: UTR and miRNA files (ID + Parent + product/other)
- **5 attributes**: Complex annotations (E-Peve original annotation)

## Dependencies

The script uses standard Unix tools:
- `bash` (v4.0+)
- `find`
- `awk`
- `grep`
- `stat`
- `sort`

No additional software installation required.
