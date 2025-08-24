#!/bin/bash
# Script to complete E-Peve methylation data processing and generate bedgraph file
# 
# This script should be run when the appropriate Bismark tools are available
# and sufficient computational resources are allocated.

echo "E-Peve Methylation Data Processing Script"
echo "========================================"

# Set working directory
WORK_DIR="/home/runner/work/timeseries_molecular/timeseries_molecular/E-Peve"
BISMARK_OUT_DIR="$WORK_DIR/output/03-Peve-bismark"
OUTPUT_DIR="/home/runner/work/timeseries_molecular/timeseries_molecular/methylation_bedgraphs"

echo "Working directory: $WORK_DIR"
echo "Bismark output directory: $BISMARK_OUT_DIR"

# Check if BAM files exist
echo ""
echo "Checking for BAM files..."
BAM_COUNT=$(find "$BISMARK_OUT_DIR" -name "*.bam" | wc -l)
echo "Found $BAM_COUNT BAM files"

if [ $BAM_COUNT -eq 0 ]; then
    echo "ERROR: No BAM files found. Bismark alignment may not be complete."
    echo "Expected files: POR-*_R1_001.fastp-trim_bismark_bt2_pe.bam"
    exit 1
fi

# List some example BAM files
echo "Example BAM files:"
find "$BISMARK_OUT_DIR" -name "*.bam" | head -3

echo ""
echo "REQUIRED STEPS TO COMPLETE E-PEVE PROCESSING:"
echo "============================================"

echo ""
echo "1. Run Bismark methylation extractor on all BAM files:"
echo "   This step extracts methylation calls from the aligned BAM files"
echo ""
echo "   Example command (adjust paths as needed):"
echo "   find $BISMARK_OUT_DIR -name '*_pe.bam' | \\"
echo "   xargs -n 1 -I{} bismark_methylation_extractor \\"
echo "   --bedGraph --counts --comprehensive --merge_non_CpG \\"
echo "   --multicore 24 --buffer_size 75% --output $BISMARK_OUT_DIR {}"

echo ""
echo "2. Run coverage2cytosine to generate merged coverage files:"
echo "   This step creates genome-wide methylation coverage files"
echo ""
echo "   Example command (adjust genome folder path):"
echo "   find $BISMARK_OUT_DIR -name '*_pe.deduplicated.bismark.cov.gz' | \\"
echo "   xargs -n 1 basename -s _pe.deduplicated.bismark.cov.gz | \\"
echo "   parallel -j 24 coverage2cytosine \\"
echo "   --genome_folder $WORK_DIR/data/ \\"
echo "   -o $BISMARK_OUT_DIR/{} \\"
echo "   --merge_CpG --zero_based \\"
echo "   $BISMARK_OUT_DIR/{}_pe.deduplicated.bismark.cov.gz"

echo ""
echo "3. Process coverage files to bedgraph format:"
echo "   This step filters for 10x coverage and creates processed files"
echo ""
echo "   # Filter for 10x coverage"
echo "   cd $BISMARK_OUT_DIR"
echo "   for f in *.CpG_report.merged_CpG_evidence.cov.gz; do"
echo "     STEM=\$(basename \"\${f}\") "
echo "     STEM=\"\${STEM%*.CpG_report.merged_CpG_evidence.cov.gz}\""
echo "     zcat \"\${f}\" | awk -F \$'\\t' 'BEGIN {OFS = FS} {if (\$5+\$6 >= 10) {print \$1, \$2, \$3, \$4}}' > \"\${STEM}\"_10x.bedgraph"
echo "   done"

echo ""
echo "   # Create processed tables"
echo "   for file in *10x.bedgraph; do"
echo "     awk -F\"\\t\" -v fname=\"\${file%_10x*}\" 'BEGIN {print \"CpG\\t\" fname}{print \"CpG_\"_\$1\"_\"\$2\"\\t\"\$4}' \"\$file\" > \"\${file%.bedgraph}_processed.txt\""
echo "   done"

echo ""
echo "4. Merge processed files:"
echo "   Use the provided Python script to merge sample-level data"
echo ""
echo "   cd $BISMARK_OUT_DIR"
echo "   python $WORK_DIR/code/11.2-merge_processed_txt.py"

echo ""
echo "5. Generate final bedgraph file:"
echo "   Use the main processing script to create species-level bedgraph"
echo ""
echo "   cd /home/runner/work/timeseries_molecular/timeseries_molecular"
echo "   python generate_methylation_bedgraphs.py"

echo ""
echo "SYSTEM REQUIREMENTS:"
echo "==================="
echo "- Bismark methylation extractor"
echo "- coverage2cytosine tool"
echo "- GNU parallel"
echo "- Python with pandas"
echo "- Sufficient disk space for intermediate files"
echo "- 24+ CPU cores recommended for parallel processing"

echo ""
echo "EXPECTED OUTPUT:"
echo "==============="
echo "- E-Peve_methylation.bedgraph in $OUTPUT_DIR"
echo "- Intermediate files in $BISMARK_OUT_DIR"

echo ""
echo "NOTE: This processing may take several hours to complete depending on"
echo "      data size and system resources."

# Check if the final output already exists
if [ -f "$OUTPUT_DIR/E-Peve_methylation.bedgraph" ]; then
    echo ""
    echo "SUCCESS: E-Peve_methylation.bedgraph already exists!"
    echo "File: $OUTPUT_DIR/E-Peve_methylation.bedgraph"
    wc -l "$OUTPUT_DIR/E-Peve_methylation.bedgraph"
else
    echo ""
    echo "STATUS: E-Peve_methylation.bedgraph not yet generated"
    echo "Run the steps above to complete processing"
fi