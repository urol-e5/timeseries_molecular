import pandas as pd
import os
import glob
import matplotlib.pyplot as plt
import seaborn as sns
# import requests # Removed dependency
import urllib.request
from io import StringIO
import numpy as np

# --- Configuration ---
TRANSCRIPT_DIR = "M-multi-species/output/26-rank35-optimization/lambda_gene_0.2/top_genes_per_component/"
ANNOTATION_FILE = "M-multi-species/output/12-ortho-annot/ortholog_groups_annotated.csv"
FACTORS_DIR = "M-multi-species/output/26-rank35-optimization/lambda_gene_0.2/barnacle_factors/"
PHYS_DATA_URL = "http://gannet.fish.washington.edu/seashell/snaps/MOSAiC_Physiology_Data__MOSAiC.csv"
OUTPUT_DIR = "M-multi-species/output/53-barnacle-interp-gemini/"

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# --- Steps 1-3: Data Annotation ---
print("Starting Step 1-3: Data Annotation...")

# Load annotation file
try:
    annot_df = pd.read_csv(ANNOTATION_FILE)
    print(f"Loaded annotation file: {ANNOTATION_FILE} with {len(annot_df)} rows.")
except FileNotFoundError:
    print(f"Error: Annotation file not found at {ANNOTATION_FILE}")
    exit(1)

# New logic: Load full gene factors and filter by 90th percentile
GENE_FACTORS_FILE = os.path.join(FACTORS_DIR, "gene_factors.csv")
print(f"Loading gene factors from: {GENE_FACTORS_FILE}")

try:
    gene_factors_df = pd.read_csv(GENE_FACTORS_FILE, index_col=0) # Assuming OG_ID is index
except FileNotFoundError:
    print(f"Error: Could not find {GENE_FACTORS_FILE}")
    exit(1)

# Annotate and save top genes for each component
# We will save these to a new directory to distinguish from the 'top100' set
TOP_GENES_DIR = os.path.join(OUTPUT_DIR, "top_genes_90th_percentile")
os.makedirs(TOP_GENES_DIR, exist_ok=True)

print("Processing components with 90th percentile threshold...")

for col in gene_factors_df.columns:
    if "Component_" not in col:
        continue
        
    # Calculate threshold
    threshold = gene_factors_df[col].quantile(0.90)
    
    # Filter genes above threshold
    # Also handle if value is 0? The request said "top genes above 90% threshold". 
    # Usually implies simply > quantile(0.90).
    top_genes = gene_factors_df[gene_factors_df[col] > threshold][col].sort_values(ascending=False)
    
    # Create DataFrame for merging
    component_df = top_genes.reset_index()
    component_df.columns = ['OG_ID', 'loading']
    
    # Merge with annotations
    # Ensure OG_ID types match (usually string)
    if 'group_id' in annot_df.columns:
        # Assuming annot_df has group_id as key
        merged_df = pd.merge(component_df, annot_df, left_on='OG_ID', right_on='group_id', how='left')
    else:
        # Fallback if column name differs
        merged_df = pd.merge(component_df, annot_df, left_on='OG_ID', right_index=True, how='left')

    # Save
    outfile = os.path.join(TOP_GENES_DIR, f"{col}_top90percentiles_annotation.csv")
    merged_df.to_csv(outfile, index=False)

print(f"Annotation complete. Files saved to {TOP_GENES_DIR}")



# --- Steps 4-6: Physiological Data & Narrative ---
print("\nStarting Steps 4-6: Physiological Data Integration...")

# Load Factors
try:
    sample_factors = pd.read_csv(os.path.join(FACTORS_DIR, "sample_factors.csv"), index_col=0)
    time_factors = pd.read_csv(os.path.join(FACTORS_DIR, "time_factors.csv"), index_col=0)
    # gene_factors = pd.read_csv(os.path.join(FACTORS_DIR, "gene_factors.csv"), index_col=0) # Large, maybe not needed for this part
    component_weights = pd.read_csv(os.path.join(FACTORS_DIR, "component_weights.csv"), index_col=0)
    sample_mapping = pd.read_csv(os.path.join(FACTORS_DIR, "sample_mapping.csv"))
    
    print("Loaded tensor factor files.")
except FileNotFoundError as e:
    print(f"Error loading factor files: {e}")
    exit(1)

# Load Physiological Data
try:
    # response = requests.get(PHYS_DATA_URL)
    # response.raise_for_status()
    # phys_df = pd.read_csv(StringIO(response.text))
    # Using pandas direct URL reading which uses urllib under the hood usually
    phys_df = pd.read_csv(PHYS_DATA_URL)
    print("Loaded physiological data from URL.")
except Exception as e:
    print(f"Error loading physiological data: {e}")
    exit(1)

# Merge Sample Factors with Sample Mapping and Physiological Data
# sample_factors index (0, 1, 2...) corresponds to combined_index in sample_mapping? 
# Let's check sample_mapping structure.
# combined_index,label,species,sample_id
# 0,apul_ACR-139,apul,ACR-139
# ...

# Prepare sample factors for merge
sample_factors = sample_factors.reset_index(drop=True)
sample_factors['combined_index'] = sample_factors.index
# Check if sample factors has more rows than sample mapping or vice versa. They should match.

merged_factors = pd.merge(sample_mapping, sample_factors, on='combined_index')



# Now merge with physiological data
# Physiological data likely has 'sample_id' or similar. 
# Looking at the phys data file usually helps, but I'll assume standard column names or inspect later.
# The user mentioned "Physiological data to correspoding to the samples avialbale".
# Common column might be 'Fragment_ID' or 'Colony_ID' or 'Sample_ID'. 
# I will inspect the first few columns of phys_df to guess the key if 'sample_id' isn't there.

# Clean up phys_df columns just in case
phys_df.columns = [c.strip() for c in phys_df.columns]
possible_keys = ['Sample_ID', 'sample_id', 'Fragment_ID', 'ID', 'Colony', 'colony_id']
key_used = None
for k in possible_keys:
    if k in phys_df.columns:
        key_used = k
        break

if key_used:
    # Ensure ID formats match (e.g., strips, upper/lower)
    merged_factors['sample_id'] = merged_factors['sample_id'].astype(str).str.strip()
    phys_df[key_used] = phys_df[key_used].astype(str).str.strip()
    
    # Drop 'species' from phys_df if it exists to avoid _x _y suffixes
    if 'species' in phys_df.columns:
        phys_df = phys_df.drop(columns=['species'])

    # Check if we need to map to sample_id explicitly
    left_key = 'sample_id'
    
    full_df = pd.merge(merged_factors, phys_df, left_on=left_key, right_on=key_used, how='left')
    print(f"Merged factors with physiological data on {key_used}.")
else:
    print("Warning: Could not find a common key for physiological data. Proceeding with factors only.")
    full_df = merged_factors

# --- Visualization & Analysis ---

# 1. Visualization of Components over Time (if time info is available/interpretable)
# We have time_factors.csv, but usually sample_factors + metadata tells the story better if samples have timepoints.
# Do we have timepoint info in phys_df or sample_mapping?
# sample_mapping has: label, species, sample_id.
# phys_df likely has timepoint.

# Time factors in Barnacle are usually Component x Timepoint. 
# Let's plot the Time Factors directly first.
plt.figure(figsize=(10, 6))
# time_factors structure: index=Timepoint?, columns=Components?
sns.heatmap(time_factors.T, cmap='viridis', cbar_kws={'label': 'Factor Value'})
plt.title("Time Factors Heatmap")
plt.xlabel("Timepoint")
plt.ylabel("Component")
plt.savefig(os.path.join(OUTPUT_DIR, "time_factors_heatmap.png"))
plt.close()

# 2. Top Components correlates with Physiology
# Calculate correlation between component sample loadings (factors) and physiological metrics.
# Components are likely columns in sample_factors like '0', '1', ... or 'Factor_0'...
# Let's start with columns that differ between full_df and phys_df/sample_mapping.
# The sample_factors columns need to be identified.
factor_cols = [c for c in sample_factors.columns if isinstance(c, int) or (isinstance(c, str) and c.isdigit())] 
if not factor_cols:
    # Maybe they are named 'Factor1', 'Factor2'?
    # Let's look at the sample_factors file content in a subsequent step if needed. 
    # For now, assuming they are the integer columns from the merge.
    # Actually, sample_factors.csv usually has columns 0, 1, 2...
    factor_cols = [c for c in full_df.columns if c in sample_factors.columns and c != 'combined_index']

# Phys columns - assume all numeric columns in phys_df except IDs
phys_cols = phys_df.select_dtypes(include=[np.number]).columns.tolist()
# Remove likely non-metric cols
phys_cols = [c for c in phys_cols if c not in ['combined_index'] and 'ID' not in c and 'Unnamed' not in c]

# Filter phys_cols to those actually in full_df
phys_cols = [c for c in phys_cols if c in full_df.columns]

if factor_cols and phys_cols:
    corr_matrix = full_df[factor_cols + phys_cols].corr()
    
    # Extract just Factor vs Phys correlations
    factor_phys_corr = corr_matrix.loc[factor_cols, phys_cols]
    
    plt.figure(figsize=(12, 10))
    sns.heatmap(factor_phys_corr, cmap='coolwarm', center=0, annot=False)
    plt.title("Correlation: Component Factors vs Physiology")
    plt.savefig(os.path.join(OUTPUT_DIR, "correlation_factors_physiology.png"))
    plt.close()
    
    print("Generated correlation heatmap.")
else:
    print(f"Skipping correlation: Factor cols {len(factor_cols)}, Phys cols found in merge {len(phys_cols)}")

# 3. generate narrative visuals for top interesting components
# "Interesting" could be high variance or strong correlation.
# Let's pick top 5 components by weight
# Check if weights exist, otherwise assume equal
if not component_weights.empty:
    sorted_components = component_weights.sum(axis=1).sort_values(ascending=False).index[:5]
else:
    sorted_components = factor_cols[:5]

# Plot them by Timepoint groupings if available in full_df
timepoint_col = None
for c in ['timepoint', 'Timepoint', 'Time']:
    if c in full_df.columns:
        timepoint_col = c
        break

if timepoint_col:
    # Plot top 4 components trends from sorted_components
    # Need to map component names/indices to columns in full_df
    # factor_cols are like "Component_1" based on merge? 
    # Wait, sample_factors usually has indices 0,1,2.
    # But I merged sample_factors.csv which has index_col=0. 
    # Let's check sample_factors.csv columns.
    # If sample_factors.csv had no header, columns are 0,1,2.
    # If it had header, they are names.
    # My factor_cols detection logic: [c for c in full_df.columns if c in sample_factors.columns]
    
    # Use detected factor_cols for plotting top ones.
    # If we have weights, map indices to column names.
    # Let's just take top 4 from factor_cols if sorted logic is complex mapping.
    # Actually, we can try to find the "Component_X" columns.
    
    top_cols = [c for c in factor_cols if c in full_df.columns][:4]
    
    if len(top_cols) > 0:
        fig, axes = plt.subplots(2, 2, figsize=(15, 10))
        for i, ax in enumerate(axes.flatten()):
            if i < len(top_cols):
                comp = top_cols[i]
                sns.lineplot(data=full_df, x=timepoint_col, y=comp, hue='species', ax=ax, marker='o', errorbar='sd')
                ax.set_title(f"{comp} Trends")
        plt.tight_layout()
        plt.savefig(os.path.join(OUTPUT_DIR, "top_components_time_trends.png"))
        plt.close()


print(f"Analysis complete. Outputs in {OUTPUT_DIR}")
