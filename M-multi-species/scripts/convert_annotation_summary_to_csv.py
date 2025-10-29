#!/usr/bin/env python3
"""
Convert annotation_summary.md to CSV format
"""

import csv
import re
import os
from pathlib import Path

# Get the script directory and construct relative paths
script_dir = Path(__file__).parent
repo_root = script_dir.parent.parent
output_dir = repo_root / "M-multi-species" / "output" / "22-Visualizing-Rank-outs"

# Input and output file paths
input_file = output_dir / "annotation_summary.md"
output_file = output_dir / "annotation_summary.csv"

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
        # Remove only the leading and trailing empty elements from pipe splitting
        if cells and cells[0] == '':
            cells = cells[1:]
        if cells and cells[-1] == '':
            cells = cells[:-1]
        
        # Check if this is the separator line (contains dashes and possibly colons/spaces)
        if all(re.match(r'^[\s:-]+$', cell) for cell in cells):
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
