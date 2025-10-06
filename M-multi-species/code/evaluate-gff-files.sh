#!/bin/bash
# Script to evaluate all GFF and GFF3 files in the repository
# Author: GitHub Copilot
# Date: 2025
# Description: Inspects all .gff and .gff3 files and generates a summary report
#              including file properties and GFF3 specification compliance

# Set up output directory and files
REPO_DIR="/home/runner/work/timeseries_molecular/timeseries_molecular"
OUTPUT_DIR="${REPO_DIR}/M-multi-species/output"
OUTPUT_FILE="${OUTPUT_DIR}/gff-file-evaluation-summary.txt"
DETAILED_FILE="${OUTPUT_DIR}/gff-file-detailed-analysis.txt"

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# Initialize output files
echo "GFF/GFF3 File Evaluation Summary" > "${OUTPUT_FILE}"
echo "Generated on: $(date)" >> "${OUTPUT_FILE}"
echo "Repository: urol-e5/timeseries_molecular" >> "${OUTPUT_FILE}"
echo "========================================" >> "${OUTPUT_FILE}"
echo "" >> "${OUTPUT_FILE}"

echo "GFF/GFF3 File Detailed Analysis" > "${DETAILED_FILE}"
echo "Generated on: $(date)" >> "${DETAILED_FILE}"
echo "========================================" >> "${DETAILED_FILE}"
echo "" >> "${DETAILED_FILE}"

# Find all GFF and GFF3 files
cd "${REPO_DIR}"
FILES=$(find . -type f \( -name "*.gff" -o -name "*.gff3" \) | sort)

# Count total files
TOTAL_FILES=$(echo "${FILES}" | wc -l)
echo "Total GFF/GFF3 files found: ${TOTAL_FILES}" >> "${OUTPUT_FILE}"
echo "" >> "${OUTPUT_FILE}"

# Table header
echo "File Properties Table:" >> "${OUTPUT_FILE}"
echo "----------------------" >> "${OUTPUT_FILE}"
printf "%-80s %-12s %-10s %-10s %-15s %-10s\n" "File Path" "Size (bytes)" "Lines" "Features" "Attribute Count" "GFF3 Valid" >> "${OUTPUT_FILE}"
printf "%-80s %-12s %-10s %-10s %-15s %-10s\n" "--------" "-----------" "-----" "--------" "---------------" "----------" >> "${OUTPUT_FILE}"

# Counter for valid/invalid files
VALID_COUNT=0
INVALID_COUNT=0

# Process each file
while IFS= read -r file; do
    if [ -z "$file" ]; then
        continue
    fi
    
    echo "Processing: ${file}" >&2
    
    # Get file size
    FILE_SIZE=$(stat -c%s "${file}")
    
    # Count lines (excluding empty lines)
    LINE_COUNT=$(grep -cv '^$' "${file}" 2>/dev/null || echo "0")
    
    # Count feature lines (non-comment, non-empty lines)
    FEATURE_COUNT=$(grep -v '^#' "${file}" 2>/dev/null | grep -cv '^$' || echo "0")
    
    # Check field count (GFF3 should have exactly 9 tab-separated fields)
    FIELD_COUNTS=$(grep -v '^#' "${file}" 2>/dev/null | grep -v '^$' | awk -F'\t' '{print NF}' | sort -u)
    
    # Check if all non-comment lines have 9 fields
    GFF3_VALID="Yes"
    if echo "${FIELD_COUNTS}" | grep -qv '^9$'; then
        GFF3_VALID="No"
        INVALID_COUNT=$((INVALID_COUNT + 1))
    else
        VALID_COUNT=$((VALID_COUNT + 1))
    fi
    
    # Get unique feature types
    FEATURE_TYPES=$(grep -v '^#' "${file}" 2>/dev/null | grep -v '^$' | awk -F'\t' '{print $3}' | sort -u | paste -sd ',' - || echo "N/A")
    
    # Count attributes in the 9th column
    # Extract the attributes column and count unique attribute keys
    ATTRIBUTE_SAMPLE=$(grep -v '^#' "${file}" 2>/dev/null | grep -v '^$' | head -1 | awk -F'\t' '{print $9}')
    ATTRIBUTE_COUNT=$(echo "${ATTRIBUTE_SAMPLE}" | grep -o '[A-Za-z_][A-Za-z0-9_]*=' | wc -l)
    
    # Print summary line
    printf "%-80s %-12s %-10s %-10s %-15s %-10s\n" "${file}" "${FILE_SIZE}" "${LINE_COUNT}" "${FEATURE_COUNT}" "${ATTRIBUTE_COUNT}" "${GFF3_VALID}" >> "${OUTPUT_FILE}"
    
    # Write detailed analysis
    echo "========================================" >> "${DETAILED_FILE}"
    echo "File: ${file}" >> "${DETAILED_FILE}"
    echo "----------------------------------------" >> "${DETAILED_FILE}"
    echo "Size: ${FILE_SIZE} bytes" >> "${DETAILED_FILE}"
    echo "Total Lines: ${LINE_COUNT}" >> "${DETAILED_FILE}"
    echo "Feature Lines: ${FEATURE_COUNT}" >> "${DETAILED_FILE}"
    echo "GFF3 Compliant (9 fields): ${GFF3_VALID}" >> "${DETAILED_FILE}"
    echo "Field counts found: ${FIELD_COUNTS}" >> "${DETAILED_FILE}"
    echo "Attribute count (first line): ${ATTRIBUTE_COUNT}" >> "${DETAILED_FILE}"
    echo "" >> "${DETAILED_FILE}"
    echo "Feature Types:" >> "${DETAILED_FILE}"
    echo "${FEATURE_TYPES}" | tr ',' '\n' | sed 's/^/  - /' >> "${DETAILED_FILE}"
    echo "" >> "${DETAILED_FILE}"
    
    # Count each feature type
    echo "Feature Type Counts:" >> "${DETAILED_FILE}"
    grep -v '^#' "${file}" 2>/dev/null | grep -v '^$' | awk -F'\t' '{print $3}' | sort | uniq -c | awk '{printf "  %-30s %s\n", $2":", $1}' >> "${DETAILED_FILE}"
    echo "" >> "${DETAILED_FILE}"
    
    # Check for GFF3 header
    GFF3_HEADER=$(head -1 "${file}" | grep -c '##gff-version 3' || echo "0")
    if [ "${GFF3_HEADER}" -eq 1 ]; then
        echo "Has GFF3 header (##gff-version 3): Yes" >> "${DETAILED_FILE}"
    else
        echo "Has GFF3 header (##gff-version 3): No" >> "${DETAILED_FILE}"
    fi
    echo "" >> "${DETAILED_FILE}"
    
    # Show sample attributes from first few lines
    echo "Sample Attributes (from first 3 feature lines):" >> "${DETAILED_FILE}"
    grep -v '^#' "${file}" 2>/dev/null | grep -v '^$' | head -3 | awk -F'\t' '{print "  " $9}' >> "${DETAILED_FILE}"
    echo "" >> "${DETAILED_FILE}"
    
    # Get source column info
    SOURCES=$(grep -v '^#' "${file}" 2>/dev/null | grep -v '^$' | awk -F'\t' '{print $2}' | sort -u | paste -sd ',' - || echo "N/A")
    echo "Sources: ${SOURCES}" >> "${DETAILED_FILE}"
    echo "" >> "${DETAILED_FILE}"
    
done <<< "${FILES}"

# Add summary statistics
echo "" >> "${OUTPUT_FILE}"
echo "========================================" >> "${OUTPUT_FILE}"
echo "Summary Statistics:" >> "${OUTPUT_FILE}"
echo "-------------------" >> "${OUTPUT_FILE}"
echo "Total files analyzed: ${TOTAL_FILES}" >> "${OUTPUT_FILE}"
echo "GFF3 compliant (9 fields): ${VALID_COUNT}" >> "${OUTPUT_FILE}"
echo "Non-compliant: ${INVALID_COUNT}" >> "${OUTPUT_FILE}"
echo "" >> "${OUTPUT_FILE}"

# Add GFF3 specification notes
echo "GFF3 Specification Notes:" >> "${OUTPUT_FILE}"
echo "-------------------------" >> "${OUTPUT_FILE}"
echo "A valid GFF3 file should have:" >> "${OUTPUT_FILE}"
echo "1. Exactly 9 tab-separated fields per line (excluding comments)" >> "${OUTPUT_FILE}"
echo "2. Fields: seqid, source, type, start, end, score, strand, phase, attributes" >> "${OUTPUT_FILE}"
echo "3. Optional header line: ##gff-version 3" >> "${OUTPUT_FILE}"
echo "4. Attributes in field 9 should be semicolon-separated key=value pairs" >> "${OUTPUT_FILE}"
echo "5. Common attributes include: ID, Name, Parent, etc." >> "${OUTPUT_FILE}"
echo "" >> "${OUTPUT_FILE}"

# Create feature matrix file showing all features across all files
MATRIX_FILE="${OUTPUT_DIR}/gff-feature-type-matrix.txt"
echo "Creating feature type matrix..." >&2

echo "GFF Feature Type Matrix" > "${MATRIX_FILE}"
echo "Generated on: $(date)" >> "${MATRIX_FILE}"
echo "========================================" >> "${MATRIX_FILE}"
echo "" >> "${MATRIX_FILE}"
echo "This matrix shows the count of each feature type across all GFF files." >> "${MATRIX_FILE}"
echo "" >> "${MATRIX_FILE}"

# First pass: collect all unique feature types across all files
echo "Collecting all unique feature types..." >&2
ALL_FEATURES=$(find . -type f \( -name "*.gff" -o -name "*.gff3" \) -exec grep -v '^#' {} \; 2>/dev/null | grep -v '^$' | awk -F'\t' '{print $3}' | sort -u)

# Create temporary directory for feature counts
TEMP_DIR=$(mktemp -d)

# Second pass: for each file, get counts for each feature type
echo "Building feature matrix..." >&2
while IFS= read -r file; do
    if [ -z "$file" ]; then
        continue
    fi
    
    # Create a temporary file with feature counts for this GFF file
    TEMP_FILE="${TEMP_DIR}/$(echo ${file} | sed 's/[\/.]/_/g').counts"
    
    # Get counts for all features in this file
    grep -v '^#' "${file}" 2>/dev/null | grep -v '^$' | awk -F'\t' '{print $3}' | sort | uniq -c | awk '{print $2 "\t" $1}' > "${TEMP_FILE}"
    
done <<< "${FILES}"

# Build the matrix
echo "Feature Type Counts Across All GFF Files" >> "${MATRIX_FILE}"
echo "=========================================" >> "${MATRIX_FILE}"
echo "" >> "${MATRIX_FILE}"

# Write header
printf "%-40s" "Feature Type" >> "${MATRIX_FILE}"
while IFS= read -r file; do
    if [ -z "$file" ]; then
        continue
    fi
    # Shorten file path for header
    SHORT_NAME=$(basename "${file}")
    printf " | %-12s" "${SHORT_NAME:0:12}" >> "${MATRIX_FILE}"
done <<< "${FILES}"
echo "" >> "${MATRIX_FILE}"

# Write separator
printf "%-40s" "----------------------------------------" >> "${MATRIX_FILE}"
while IFS= read -r file; do
    if [ -z "$file" ]; then
        continue
    fi
    printf " | %-12s" "------------" >> "${MATRIX_FILE}"
done <<< "${FILES}"
echo "" >> "${MATRIX_FILE}"

# For each feature type, write a row with counts from each file
while IFS= read -r feature; do
    if [ -z "$feature" ]; then
        continue
    fi
    
    printf "%-40s" "${feature}" >> "${MATRIX_FILE}"
    
    while IFS= read -r file; do
        if [ -z "$file" ]; then
            continue
        fi
        
        TEMP_FILE="${TEMP_DIR}/$(echo ${file} | sed 's/[\/.]/_/g').counts"
        
        # Get count for this feature in this file
        COUNT=$(grep "^${feature}	" "${TEMP_FILE}" 2>/dev/null | awk '{print $2}')
        if [ -z "$COUNT" ]; then
            COUNT=0
        fi
        
        printf " | %12s" "${COUNT}" >> "${MATRIX_FILE}"
    done <<< "${FILES}"
    echo "" >> "${MATRIX_FILE}"
done <<< "${ALL_FEATURES}"

# Add summary section
echo "" >> "${MATRIX_FILE}"
echo "========================================" >> "${MATRIX_FILE}"
echo "Feature Type Summary:" >> "${MATRIX_FILE}"
echo "--------------------" >> "${MATRIX_FILE}"

# Count how many files each feature appears in
while IFS= read -r feature; do
    if [ -z "$feature" ]; then
        continue
    fi
    
    FILE_COUNT=0
    TOTAL_COUNT=0
    
    while IFS= read -r file; do
        if [ -z "$file" ]; then
            continue
        fi
        
        TEMP_FILE="${TEMP_DIR}/$(echo ${file} | sed 's/[\/.]/_/g').counts"
        COUNT=$(grep "^${feature}	" "${TEMP_FILE}" 2>/dev/null | awk '{print $2}')
        
        if [ -n "$COUNT" ] && [ "$COUNT" -gt 0 ]; then
            FILE_COUNT=$((FILE_COUNT + 1))
            TOTAL_COUNT=$((TOTAL_COUNT + COUNT))
        fi
    done <<< "${FILES}"
    
    printf "  %-30s: present in %2d files, total count: %8d\n" "${feature}" "${FILE_COUNT}" "${TOTAL_COUNT}" >> "${MATRIX_FILE}"
done <<< "${ALL_FEATURES}"

# Clean up temporary directory
rm -rf "${TEMP_DIR}"

# Create a transposed version showing files as rows and features as columns
TRANSPOSED_FILE="${OUTPUT_DIR}/gff-feature-type-by-file.txt"
echo "Creating transposed feature matrix..." >&2

echo "GFF Files with Feature Type Counts" > "${TRANSPOSED_FILE}"
echo "Generated on: $(date)" >> "${TRANSPOSED_FILE}"
echo "========================================" >> "${TRANSPOSED_FILE}"
echo "" >> "${TRANSPOSED_FILE}"
echo "Each file is listed with counts for all feature types present." >> "${TRANSPOSED_FILE}"
echo "" >> "${TRANSPOSED_FILE}"

# For each file, show all its feature counts
while IFS= read -r file; do
    if [ -z "$file" ]; then
        continue
    fi
    
    echo "========================================" >> "${TRANSPOSED_FILE}"
    echo "File: ${file}" >> "${TRANSPOSED_FILE}"
    echo "----------------------------------------" >> "${TRANSPOSED_FILE}"
    
    # Get all feature types and counts for this file
    grep -v '^#' "${file}" 2>/dev/null | grep -v '^$' | awk -F'\t' '{print $3}' | sort | uniq -c | awk '{printf "  %-30s: %8d\n", $2, $1}' >> "${TRANSPOSED_FILE}"
    
    # Calculate total features
    TOTAL=$(grep -v '^#' "${file}" 2>/dev/null | grep -v '^$' | wc -l)
    echo "  ----------------------------------------" >> "${TRANSPOSED_FILE}"
    printf "  %-30s: %8d\n" "TOTAL FEATURES" "${TOTAL}" >> "${TRANSPOSED_FILE}"
    echo "" >> "${TRANSPOSED_FILE}"
    
done <<< "${FILES}"

echo "Analysis complete. Results saved to:" >&2
echo "  Summary: ${OUTPUT_FILE}" >&2
echo "  Detailed: ${DETAILED_FILE}" >&2
echo "  Feature Matrix: ${MATRIX_FILE}" >&2
echo "  Feature by File: ${TRANSPOSED_FILE}" >&2
