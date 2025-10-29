#!/usr/bin/env python3
"""
Convert annotation_summary.md to CSV format
"""

import csv
import re

# Input and output file paths
input_file = "/home/runner/work/timeseries_molecular/timeseries_molecular/M-multi-species/output/22-Visualizing-Rank-outs/annotation_summary.md"
output_file = "/home/runner/work/timeseries_molecular/timeseries_molecular/M-multi-species/output/22-Visualizing-Rank-outs/annotation_summary.csv"

# Read the markdown file
with open(input_file, 'r') as f:
    lines = f.readlines()

# Find the header row and data rows
header = None
data_rows = []

for i, line in enumerate(lines):
    line = line.strip()
    
    # Skip empty lines and title
    if not line or line.startswith('#'):
        continue
    
    # Skip the "Total files processed" line
    if line.startswith('Total files processed'):
        continue
    
    # Check if this is a table row (contains pipe characters)
    if '|' in line:
        # Split by pipe and clean up
        cells = [cell.strip() for cell in line.split('|')]
        # Remove empty first and last elements if they exist
        cells = [cell for cell in cells if cell]
        
        # Check if this is the separator line (contains dashes)
        if all(re.match(r'^-+$', cell) for cell in cells):
            continue
        
        # If we haven't found the header yet, this is it
        if header is None:
            header = cells
        else:
            data_rows.append(cells)

# Write to CSV
with open(output_file, 'w', newline='') as f:
    writer = csv.writer(f)
    
    # Write header
    if header:
        writer.writerow(header)
    
    # Write data rows
    for row in data_rows:
        writer.writerow(row)

print(f"Converted {len(data_rows)} rows from {input_file} to {output_file}")
