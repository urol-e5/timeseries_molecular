#!/usr/bin/env Rscript

# Script to normalize and merge ortholog expression data from three species
# Uses conservative library size normalization and log transformation

# No external packages required - uses base R functions only

# Set working directory
setwd("/Users/sr320/GitHub/timeseries_molecular/M-multi-species/output/13.00-multiomics-barnacle")

# Function to read and process expression data
read_expression_data <- function(file_path, species_name) {
  cat(sprintf("Reading %s data from %s...\n", species_name, file_path))

  # Read CSV file
  data <- read.csv(file_path, row.names = 1, check.names = FALSE)

  # Extract sample information from column names
  sample_cols <- colnames(data)
  sample_info <- data.frame(
    sample = sample_cols,
    species = species_name,
    stringsAsFactors = FALSE
  )

  # Add species prefix to ortholog IDs for clarity
  rownames(data) <- paste0(species_name, "_", rownames(data))

  return(list(
    counts = data,
    samples = sample_info
  ))
}

# Function to perform conservative normalization using library size factors
normalize_conservative <- function(count_data, species_name) {
  cat(sprintf("Normalizing %s data with library size normalization...\n", species_name))

  # Filter out genes with very low counts (conservative approach)
  keep <- rowSums(count_data >= 10) >= 3  # At least 3 samples with >= 10 counts
  filtered_counts <- count_data[keep, , drop = FALSE]

  cat(sprintf("  Retained %d genes after filtering (from %d)\n",
              nrow(filtered_counts), nrow(count_data)))

  # Calculate library sizes (total counts per sample)
  lib_sizes <- colSums(filtered_counts)

  # Calculate size factors (relative to median library size)
  median_lib_size <- median(lib_sizes)
  size_factors <- lib_sizes / median_lib_size

  # Normalize counts by size factors
  normalized_counts <- t(t(filtered_counts) / size_factors)

  # Apply log2 transformation with pseudocount for stability
  # Using log2(x + 1) which is more conservative than log2(x) for low counts
  normalized_log <- log2(normalized_counts + 1)

  return(normalized_log)
}

# Main execution
main <- function() {
  cat("Starting ortholog expression normalization and merging...\n\n")

  # File paths
  files <- c(
    "apul" = "apul_ortholog_expression.csv",
    "peve" = "peve_ortholog_expression.csv",
    "ptua" = "ptua_ortholog_expression.csv"
  )

  # Read all datasets
  datasets <- list()
  sample_info_all <- data.frame()

  for (species in names(files)) {
    result <- read_expression_data(files[species], species)
    datasets[[species]] <- result$counts
    sample_info_all <- rbind(sample_info_all, result$samples)
  }

  cat(sprintf("\nRead data for %d samples across 3 species\n", nrow(sample_info_all)))

  # Normalize each dataset
  normalized_data <- list()

  for (species in names(datasets)) {
    cat(sprintf("\nProcessing %s...\n", species))
    counts <- datasets[[species]]

    cat(sprintf("  Raw data: %d genes × %d samples\n", nrow(counts), ncol(counts)))

    # Normalize with conservative library size normalization
    normalized <- normalize_conservative(counts, species)

    cat(sprintf("  Normalized data: %d genes × %d samples\n", nrow(normalized), ncol(normalized)))

    normalized_data[[species]] <- normalized
  }

  # Find common orthologs across all species
  cat("\nFinding common orthologs...\n")

  # Extract ortholog IDs from row names
  get_ortholog_ids <- function(data, species) {
    gsub(paste0("^", species, "_"), "", rownames(data))
  }

  common_orthologs <- Reduce(intersect, list(
    get_ortholog_ids(normalized_data$apul, "apul"),
    get_ortholog_ids(normalized_data$peve, "peve"),
    get_ortholog_ids(normalized_data$ptua, "ptua")
  ))

  cat(sprintf("Found %d orthologs present in all three species\n", length(common_orthologs)))

  # Filter datasets to common orthologs
  filtered_data <- list()

  for (species in names(normalized_data)) {
    species_orthologs <- paste0(species, "_", common_orthologs)
    filtered_data[[species]] <- normalized_data[[species]][species_orthologs, ]
    rownames(filtered_data[[species]]) <- common_orthologs
  }

  # Merge datasets
  cat("\nMerging datasets...\n")

  # Combine all normalized data
  merged_data <- cbind(
    data.frame(ortholog_id = common_orthologs),
    filtered_data$apul,
    filtered_data$peve,
    filtered_data$ptua
  )

  # Add species prefixes to column names for clarity
  apul_cols <- paste0("apul_", colnames(filtered_data$apul))
  peve_cols <- paste0("peve_", colnames(filtered_data$peve))
  ptua_cols <- paste0("ptua_", colnames(filtered_data$ptua))

  colnames(merged_data) <- c("ortholog_id", apul_cols, peve_cols, ptua_cols)

  cat(sprintf("Merged data: %d orthologs × %d samples\n",
              nrow(merged_data), ncol(merged_data) - 1))

  # Note: Samples don't have direct correspondence across species
  cat("\nNote: Sample IDs are species-specific and don't have direct correspondence across species.\n")
  cat("Including all samples for each species in the final dataset.\n")

  # Check for any missing values in the merged data
  missing_count <- sum(!is.finite(as.matrix(merged_data[, -1])))
  cat(sprintf("Found %d missing/infinite values in the dataset\n", missing_count))

  if (missing_count > 0) {
    cat("Removing orthologs with missing values...\n")
    # Keep only orthologs where all samples have finite values
    finite_rows <- apply(is.finite(as.matrix(merged_data[, -1])), 1, all)
    final_data <- merged_data[finite_rows, ]
    cat(sprintf("Retained %d orthologs with complete data across all samples\n",
                nrow(final_data)))
  } else {
    final_data <- merged_data
    cat("All values are finite - no filtering needed\n")
  }

  cat(sprintf("Final dataset: %d orthologs × %d samples\n",
              nrow(final_data), ncol(final_data) - 1))

  # Save results
  output_file <- "normalized_merged_ortholog_expression.csv"

  cat(sprintf("\nSaving results to %s...\n", output_file))

  write.csv(final_data, output_file, row.names = FALSE)

  cat("Normalization and merging completed successfully!\n")

  # Print summary
  cat("\nSummary:\n")
  cat(sprintf("  - Input orthologs per species: ~%d (after filtering low-count genes)\n",
              nrow(normalized_data$apul)))
  cat(sprintf("  - Common orthologs: %d\n", length(common_orthologs)))
  cat(sprintf("  - Total samples: %d (species-specific samples included)\n",
              ncol(final_data) - 1))
  cat(sprintf("  - Final dataset: %d orthologs × %d samples\n",
              nrow(final_data), ncol(final_data) - 1))
}

# Run main function
if (!interactive()) {
  main()
}
