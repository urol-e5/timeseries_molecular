import pandas as pd
import glob

# Define the directory and file pattern
file_pattern = '*_processed.txt'  # Adjust for your file type and location

# Get a list of all matching file paths
file_paths = glob.glob(file_pattern)

# Read each file into a DataFrame and store them in a list
list_of_dfs = [pd.read_csv(file_path, sep ='\t') for file_path in file_paths]

# merge dfs together
if list_of_dfs:
    merged_data = list_of_dfs[0]
    for i in range(1, len(list_of_dfs)):
        merged_data = pd.merge(merged_data, list_of_dfs[i], on="CpG", how="outer")
else:
    merged_data = pd.DataFrame() # Or handle the empty list case as needed


# Save the merged DataFrame to a CSV file
if not merged_data.empty: # Only save if the DataFrame is not empty
    merged_data.to_csv('merged-WGBS-CpG-counts.csv', index=False)

# Remove rows with NA values
if not merged_data.empty:
    filtered_data = merged_data.dropna()  # This removes any row with at least one NA value

# Save the filtered DataFrame to a new CSV file
filtered_data.to_csv('merged-WGBS-CpG-counts_filtered.csv', index=False)
