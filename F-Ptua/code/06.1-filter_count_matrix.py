#!/usr/bin/env python3
"""
Filter count matrix to remove rows where more than 50% of samples have counts < 10
"""

import sys
import os

def filter_count_matrix(input_file, output_file, min_count=10, max_low_count_ratio=0.5):
    """
    Filter count matrix based on count thresholds
    
    Args:
        input_file: Path to input count matrix
        output_file: Path to output filtered matrix
        min_count: Minimum count threshold (default: 10)
        max_low_count_ratio: Maximum ratio of samples below threshold (default: 0.5)
    """
    
    print(f"Filtering {input_file}...")
    print(f"Keeping rows where at least {(1-max_low_count_ratio)*100:.1f}% of samples have counts >= {min_count}")
    
    # Count total rows for progress tracking
    total_lines = sum(1 for _ in open(input_file, 'r'))
    print(f"Total lines in file: {total_lines}")
    
    kept_rows = 0
    filtered_rows = 0
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        line_num = 0
        
        for line in infile:
            line_num += 1
            
            # Write header lines (starting with #) and column header
            if line.startswith('#') or line.startswith('Geneid'):
                outfile.write(line)
                continue
            
            # Skip empty lines
            if line.strip() == '':
                continue
            
            # Process data lines
            parts = line.strip().split('\t')
            
            # Skip if line doesn't have enough columns
            if len(parts) < 7:
                outfile.write(line)
                continue
            
            # Extract count columns (starting from column 7, index 6)
            try:
                count_values = [float(parts[i]) for i in range(6, len(parts))]
            except ValueError:
                # If we can't parse counts, keep the line
                outfile.write(line)
                continue
            
            # Count how many samples are below threshold
            low_count_samples = sum(1 for count in count_values if count < min_count)
            total_samples = len(count_values)
            
            # Calculate ratio of low-count samples
            low_count_ratio = low_count_samples / total_samples if total_samples > 0 else 0
            
            # Keep row if ratio of low-count samples is below threshold
            if low_count_ratio <= max_low_count_ratio:
                outfile.write(line)
                kept_rows += 1
            else:
                filtered_rows += 1
            
            # Progress update every 1000 lines
            if line_num % 1000 == 0:
                print(f"Processed {line_num}/{total_lines} lines...")
    
    print(f"\nFiltering complete!")
    print(f"Kept rows: {kept_rows}")
    print(f"Filtered rows: {filtered_rows}")
    print(f"Output saved to: {output_file}")

def main():
    if len(sys.argv) != 3:
        print("Usage: python filter_count_matrix.py <input_file> <output_file>")
        print("Example: python filter_count_matrix.py Apul_lncRNA_counts.txt Apul_lncRNA_counts_filtered.txt")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found")
        sys.exit(1)
    
    # Filter the matrix
    filter_count_matrix(input_file, output_file)

if __name__ == "__main__":
    main()
