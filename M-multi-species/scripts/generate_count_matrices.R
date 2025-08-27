#!/usr/bin/env Rscript

# Generate count matrices for physiology data (lipidomics and metabolomics)
# This script creates CSV count matrices in the same format as CpG count matrices

cat("Loading lipidomics data...\n")

# Load lipidomics data
lipid_data <- readRDS("M-multi-species/data/lipidomics/lipids.rds")

# Create sample names in expected format (colony-timepoint)
lipid_data$sample_name <- paste(lipid_data$colony, lipid_data$timepoint, sep = "-")

cat("Processing lipidomics data...\n")
cat("Species found:", as.character(unique(lipid_data$species)), "\n")
cat("Number of samples:", length(unique(lipid_data$sample_name)), "\n")
cat("Number of lipids:", length(unique(lipid_data$lipid)), "\n")

# Convert to wide format manually
unique_lipids <- unique(lipid_data$lipid)
unique_samples <- sort(unique(lipid_data$sample_name))

# Create empty matrix
lipid_matrix <- matrix(0, nrow = length(unique_lipids), ncol = length(unique_samples))
rownames(lipid_matrix) <- unique_lipids
colnames(lipid_matrix) <- unique_samples

# Fill the matrix
cat("Creating lipidomics count matrix...\n")
for (i in 1:nrow(lipid_data)) {
  lipid <- lipid_data$lipid[i]
  sample <- lipid_data$sample_name[i]
  value <- lipid_data$nmol.g_plasma.ug_protein[i]
  
  if (!is.na(value)) {
    lipid_matrix[lipid, sample] <- value
  }
}

# Save combined lipidomics matrix
cat("Saving combined lipidomics count matrix...\n")
write.csv(lipid_matrix, "M-multi-species/output/12-generate-physiology-count-matrices/combined-lipidomics-count-matrix.csv", quote = FALSE)

# Create species-specific matrices
species_list <- list(
  "Acropora" = "D-Apul",
  "Porites" = "E-Peve", 
  "Pocillopora" = "F-Ptua"
)

for (species in names(species_list)) {
  prefix <- species_list[[species]]
  cat("Creating", species, "lipidomics matrix...\n")
  
  # Get samples for this species
  species_data <- lipid_data[lipid_data$species == species, ]
  species_samples <- sort(unique(species_data$sample_name))
  
  # Create species matrix
  species_matrix <- lipid_matrix[, species_samples, drop = FALSE]
  
  # Save species matrix
  filename <- paste0("M-multi-species/output/12-generate-physiology-count-matrices/", prefix, "-lipidomics-count-matrix.csv")
  write.csv(species_matrix, filename, quote = FALSE)
  
  cat("Saved", species, "matrix with", nrow(species_matrix), "lipids and", ncol(species_matrix), "samples\n")
}

cat("\nLipidomics matrices completed.\n\n")

# Now process metabolomics data
cat("Loading processed metabolomics data...\n")

# Load the processed CSV file
metab_data <- read.csv("M-multi-species/data/metabolomics/processed_metabolomics.csv", stringsAsFactors = FALSE)

cat("Processing metabolomics data...\n")
cat("Number of compounds:", nrow(metab_data), "\n")
cat("Number of columns:", ncol(metab_data), "\n")

# Convert to long format
metab_long <- data.frame()
for (i in 2:ncol(metab_data)) {
  sample_name <- colnames(metab_data)[i]
  for (j in 1:nrow(metab_data)) {
    compound <- metab_data$compound[j]
    value <- metab_data[j, i]
    
    # Skip if value is NA or if sample is PBQC
    if (!is.na(value) && !grepl("PBQC", sample_name)) {
      # Extract colony and timepoint
      colony <- gsub("_.*", "", sample_name)
      timepoint <- gsub(".*_", "", sample_name)
      
      # Determine species
      species <- if (grepl("ACR", colony)) "Acropora" 
                else if (grepl("POR", colony)) "Porites"
                else if (grepl("POC", colony)) "Pocillopora"
                else NA
      
      if (!is.na(species)) {
        sample_formatted <- paste(colony, timepoint, sep = "-")
        
        metab_long <- rbind(metab_long, data.frame(
          compound = compound,
          sample_name = sample_formatted,
          rel.quant = value,
          species = species,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

cat("Created long format data with", nrow(metab_long), "observations\n")
cat("Species found:", unique(metab_long$species), "\n")

# Convert to wide format manually
unique_compounds <- unique(metab_long$compound)
unique_metab_samples <- sort(unique(metab_long$sample_name))

cat("Number of compounds:", length(unique_compounds), "\n")
cat("Number of samples:", length(unique_metab_samples), "\n")

# Create empty matrix
metab_matrix <- matrix(0, nrow = length(unique_compounds), ncol = length(unique_metab_samples))
rownames(metab_matrix) <- unique_compounds
colnames(metab_matrix) <- unique_metab_samples

# Fill the matrix
cat("Creating metabolomics count matrix...\n")
for (i in 1:nrow(metab_long)) {
  compound <- metab_long$compound[i]
  sample <- metab_long$sample_name[i]
  value <- metab_long$rel.quant[i]
  
  if (!is.na(value)) {
    metab_matrix[compound, sample] <- value
  }
}

# Fix column names to ensure they use hyphens not periods
colnames(metab_matrix) <- gsub("\\.", "-", colnames(metab_matrix))

# Save combined metabolomics matrix
cat("Saving combined metabolomics count matrix...\n")
write.csv(metab_matrix, "M-multi-species/output/12-generate-physiology-count-matrices/combined-metabolomics-count-matrix.csv", quote = FALSE)

# Create species-specific metabolomics matrices
for (species in names(species_list)) {
  prefix <- species_list[[species]]
  cat("Creating", species, "metabolomics matrix...\n")
  
  # Get samples for this species from metabolomics data
  species_metab_data <- metab_long[metab_long$species == species, ]
  species_samples <- sort(unique(species_metab_data$sample_name))
  
  # Get columns that exist in the metabolomics matrix
  existing_samples <- intersect(species_samples, colnames(metab_matrix))
  
  if (length(existing_samples) > 0) {
    species_matrix <- metab_matrix[, existing_samples, drop = FALSE]
    
    # Save species matrix
    filename <- paste0("M-multi-species/output/12-generate-physiology-count-matrices/", prefix, "-metabolomics-count-matrix.csv")
    write.csv(species_matrix, filename, quote = FALSE)
    
    cat("Saved", species, "matrix with", nrow(species_matrix), "compounds and", ncol(species_matrix), "samples\n")
  }
}

cat("\nAll count matrices generated successfully!\n")
cat("Files saved in: M-multi-species/output/12-generate-physiology-count-matrices/\n")