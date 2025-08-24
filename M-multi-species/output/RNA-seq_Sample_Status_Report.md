# RNA-seq Sample Status Report

This report analyzes which RNA-seq samples were dropped from count matrices and identifies the reasons for any exclusions.

## Summary

**All submitted samples successfully made it through to the final count matrices:**

| Species | Submitted | Aligned | In Count Matrix | Success Rate |
|---------|-----------|---------|-----------------|--------------|
| *Acropora pulchra* | 40 | 40 | 40 | 100% |
| *Porites evermanni* | 38 | 38 | 38 | 100% |
| *Pocillopora tuahiniensis* | 39 | 39 | 39 | 100% |
| **Total** | **117** | **117** | **117** | **100%** |

## Key Findings

### No Samples Were Dropped

Contrary to the initial question about which samples were dropped, **all 117 submitted RNA-seq samples successfully completed the pipeline** and are included in the final count matrices:

- All samples passed trimming and quality control
- All samples successfully aligned to their respective genomes  
- All samples are included in the final count matrices

### Quality Issues Identified

While no samples were dropped, several samples had suboptimal quality metrics:

#### *Acropora pulchra* (6 samples with quality issues):
- **ACR-225-TP2**: Low total amount (936 ng)
- **ACR-229-TP2**: Low total amount (957 ng)
- **ACR-265-TP3**: Low total amount (905 ng) + **documented as outlier in sRNAseq analysis**
- **ACR-145-TP3**: Low total amount (905 ng)
- **ACR-173-TP2**: Low concentration (8.78 ng/μL) + Low total amount (790 ng)
- **ACR-173-TP4**: Low total amount (980 ng)

#### *Porites evermanni* (3 samples with quality issues):
- **POR-73-TP2**: **Zero concentration (0 ng/μL)** - sample likely failed extraction
- **POR-74-TP2**: Low concentration (9.55 ng/μL) + Low total amount (774 ng)
- **POR-83-TP1**: Low total amount (936 ng)

#### *Pocillopora tuahiniensis* (3 samples with quality issues):
- **POC-52-TP1**: Low total amount (871 ng)
- **POC-53-TP1**: Low concentration (8.23 ng/μL) + Low total amount (716 ng)
- **POC-57-TP1**: Low total amount (940 ng)

## Understanding the "40 samples" Reference

The issue mentions "40 samples" but this likely refers to:
- 40 *Acropora pulchra* samples (which is correct)
- NOT a total of 40 samples across all species

The actual total is 117 samples across the three species.

## Technical Notes

### Count Matrix Structure
- **Acropora**: 80 columns in count matrix = 40 unique samples × 2 (paired-end reads)
- **Porites**: 38 columns = 38 samples (single column per sample)
- **Pocillopora**: 39 columns = 39 samples (single column per sample)

### Quality Thresholds Used
- Low concentration: < 10 ng/μL
- Low total amount: < 1000 ng

### Known Outliers
- **ACR-265-TP3**: Documented in analysis files as an outlier based on PCA analysis of sRNAseq data (but still included in count matrices)
- **POR-73-TP2**: Had zero RNA concentration but still produced sequencing data

## Conclusion

**No RNA-seq samples were actually dropped from the count matrices.** All 117 submitted samples are represented in the final datasets. While some samples had suboptimal starting material quality (low concentration or total amount), the sequencing pipeline was robust enough to generate usable data from all samples.

The documented "outlier" sample (ACR-265-TP3) was identified during downstream analysis but was not removed from the count matrices, allowing researchers to make their own decisions about inclusion/exclusion in specific analyses.