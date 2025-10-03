
#!/usr/bin/env Rscript
#
# To run this script from command line, run the following:
# Rscript path/to/10-format-miRNA-counts.R \
#         path/to/raw/counts.txt \
#         path/to/metadata.csv \
#         path/to/output/file.txt
#
# UPDATED to rename using the following format: ColonyID-Timepoint (e.g. ACR-145-TP4)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

# Read arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript format_counts.R <counts_file> <metadata_file> <output_file>")
}

counts_file <- args[1]
metadata_file <- args[2]
output_file <- args[3]

# Read files
counts_df <- read.delim(counts_file, check.names = FALSE)  # for .txt
metadata_df <- read.csv(metadata_file, check.names = FALSE)  # for .csv

# Only keep rows with "Y" in MIRNA column
counts_df <- counts_df %>% filter(MIRNA == "Y")

# Remove columns 1 and 3 (Coords and MIRNA)
counts_df <- counts_df %>% select(-1, -3)

# Make column 2 the rownames
counts_df <- counts_df %>% column_to_rownames(var = colnames(counts_df)[1])

# Only keep characters before the first "-"

colnames(counts_df) <- sapply(colnames(counts_df), function(x) {
  sub("-.*", "", x) 
})

# Rename columns by matching to metadata (AzentaSampleName → ColonyID-Timepoint)
#  Create mapping

name_map <- metadata_df %>%
  select(AzentaSampleName, ColonyID, Timepoint) %>%
  mutate(new_name = paste0(ColonyID, "-", Timepoint)) %>%
  select(AzentaSampleName, new_name)

# Save old-named df to check mapping
old_counts_df <- counts_df

# Apply mapping where matches exist
current_names <- colnames(counts_df)

updated_names <- sapply(current_names, function(x) {
  if (x %in% name_map$AzentaSampleName) {
    name_map$new_name[name_map$AzentaSampleName == x]
  } else {
    x
  }
})

colnames(counts_df) <- updated_names

# Confirm mapping is accurate
ifelse(identical(unname(old_counts_df[1,]), unname(counts_df[1,])), print("accurate mapping"), print("mapping failed"))

# Save output
write_tsv(counts_df %>% rownames_to_column(var = "Name"), output_file)
