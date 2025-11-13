#!/usr/bin/env Rscript
#
# Compare multiway interactions across 3 species
#
# This script reads 3 CSV files containing multi-way interactions from different species,
# counts occurrences of "CpG", "lncRNA", and "miRNA" in each row,
# generates summary statistics, and creates visualizations.
#
# URLs for the 3 CSV files:
# - Species 1: https://gannet.fish.washington.edu/v1_web/owlshell/bu-github/ConTra/output/context_dependent_analysis_20251108_151255/tables/multi_way_interactions.csv
# - Species 2: https://gannet.fish.washington.edu/v1_web/owlshell/bu-github/ConTra/output/context_dependent_analysis_20251108_144034/tables/multi_way_interactions.csv
# - Species 3: https://gannet.fish.washington.edu/v1_web/owlshell/bu-github/ConTra/output/context_dependent_analysis_20251108_140602/tables/multi_way_interactions.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
})

# Create output directory
output_dir <- "M-multi-species/output/30-compare-multiway-interactions"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# URLs for the three CSV files
urls <- c(
  species_1 = "https://gannet.fish.washington.edu/v1_web/owlshell/bu-github/ConTra/output/context_dependent_analysis_20251108_151255/tables/multi_way_interactions.csv",
  species_2 = "https://gannet.fish.washington.edu/v1_web/owlshell/bu-github/ConTra/output/context_dependent_analysis_20251108_144034/tables/multi_way_interactions.csv",
  species_3 = "https://gannet.fish.washington.edu/v1_web/owlshell/bu-github/ConTra/output/context_dependent_analysis_20251108_140602/tables/multi_way_interactions.csv"
)

cat("=== Comparing Multiway Interactions Across 3 Species ===\n\n")

# Function to count feature occurrences in each row
count_features <- function(df) {
  # Count occurrences in each row across all columns
  df$CpG_count <- apply(df, 1, function(row) {
    sum(grepl("cpg", tolower(as.character(row))))
  })
  
  df$lncRNA_count <- apply(df, 1, function(row) {
    sum(grepl("lncrna", tolower(as.character(row))))
  })
  
  df$miRNA_count <- apply(df, 1, function(row) {
    sum(grepl("mirna", tolower(as.character(row))))
  })
  
  return(df)
}

# Read and process each file
all_data <- list()
summary_stats <- data.frame()

for (species_name in names(urls)) {
  url <- urls[species_name]
  local_file <- file.path(output_dir, paste0(species_name, "_multi_way_interactions.csv"))
  
  cat(sprintf("Processing %s...\n", species_name))
  
  # Try to download or read local file
  tryCatch({
    if (file.exists(local_file)) {
      cat(sprintf("  Reading existing file: %s\n", local_file))
      df <- read_csv(local_file, show_col_types = FALSE)
    } else {
      cat(sprintf("  Downloading from: %s\n", url))
      df <- read_csv(url, show_col_types = FALSE)
      write_csv(df, local_file)
      cat(sprintf("  Saved to: %s\n", local_file))
    }
    
    cat(sprintf("  Loaded %d rows with %d columns\n", nrow(df), ncol(df)))
    
    # Count feature occurrences
    df <- count_features(df)
    
    # Save annotated data
    annotated_file <- file.path(output_dir, paste0(species_name, "_annotated.csv"))
    write_csv(df, annotated_file)
    cat(sprintf("  Saved annotated data to: %s\n", annotated_file))
    
    all_data[[species_name]] <- df
    
    # Calculate summary statistics
    stats <- data.frame(
      species = species_name,
      total_rows = nrow(df),
      cpg_total = sum(df$CpG_count),
      lncrna_total = sum(df$lncRNA_count),
      mirna_total = sum(df$miRNA_count),
      cpg_mean = mean(df$CpG_count),
      lncrna_mean = mean(df$lncRNA_count),
      mirna_mean = mean(df$miRNA_count),
      cpg_rows_with_occurrence = sum(df$CpG_count > 0),
      lncrna_rows_with_occurrence = sum(df$lncRNA_count > 0),
      mirna_rows_with_occurrence = sum(df$miRNA_count > 0)
    )
    
    summary_stats <- rbind(summary_stats, stats)
    
    cat("\n")
    
  }, error = function(e) {
    cat(sprintf("  Error processing %s: %s\n\n", species_name, e$message))
  })
}

if (nrow(summary_stats) == 0) {
  cat("\nError: No data files could be loaded.\n")
  cat("To use this script, ensure network access to gannet.fish.washington.edu\n")
  cat("or manually download the CSV files to the output directory.\n")
  quit(status = 1)
}

# Save summary statistics
cat("\n=== Summary Statistics ===\n\n")
print(summary_stats)

summary_file <- file.path(output_dir, "summary_statistics.csv")
write_csv(summary_stats, summary_file)
cat(sprintf("\nSummary statistics saved to: %s\n", summary_file))

# Create visualizations
cat("\n=== Creating Visualizations ===\n\n")

# Prepare data for plotting
plot_data <- summary_stats %>%
  select(species, cpg_total, lncrna_total, mirna_total) %>%
  pivot_longer(cols = c(cpg_total, lncrna_total, mirna_total),
               names_to = "feature", values_to = "count") %>%
  mutate(feature = recode(feature,
                         cpg_total = "CpG",
                         lncrna_total = "lncRNA",
                         mirna_total = "miRNA"))

# 1. Total occurrences by species
p1 <- ggplot(plot_data, aes(x = species, y = count, fill = feature)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  labs(title = "Total Feature Occurrences by Species",
       x = "Species",
       y = "Total Occurrences",
       fill = "Feature") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave(file.path(output_dir, "total_occurrences_by_species.png"), 
       p1, width = 10, height = 6, dpi = 300)
cat(sprintf("Saved: %s\n", file.path(output_dir, "total_occurrences_by_species.png")))

# 2. Mean occurrences per row
mean_data <- summary_stats %>%
  select(species, cpg_mean, lncrna_mean, mirna_mean) %>%
  pivot_longer(cols = c(cpg_mean, lncrna_mean, mirna_mean),
               names_to = "feature", values_to = "mean") %>%
  mutate(feature = recode(feature,
                         cpg_mean = "CpG",
                         lncrna_mean = "lncRNA",
                         mirna_mean = "miRNA"))

p2 <- ggplot(mean_data, aes(x = species, y = mean, fill = feature)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  labs(title = "Mean Feature Occurrences per Row by Species",
       x = "Species",
       y = "Mean Occurrences per Row",
       fill = "Feature") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave(file.path(output_dir, "mean_occurrences_by_species.png"), 
       p2, width = 10, height = 6, dpi = 300)
cat(sprintf("Saved: %s\n", file.path(output_dir, "mean_occurrences_by_species.png")))

# 3. Percentage of rows with occurrence
pct_data <- summary_stats %>%
  mutate(
    cpg_pct = (cpg_rows_with_occurrence / total_rows) * 100,
    lncrna_pct = (lncrna_rows_with_occurrence / total_rows) * 100,
    mirna_pct = (mirna_rows_with_occurrence / total_rows) * 100
  ) %>%
  select(species, cpg_pct, lncrna_pct, mirna_pct) %>%
  pivot_longer(cols = c(cpg_pct, lncrna_pct, mirna_pct),
               names_to = "feature", values_to = "percentage") %>%
  mutate(feature = recode(feature,
                         cpg_pct = "CpG",
                         lncrna_pct = "lncRNA",
                         mirna_pct = "miRNA"))

p3 <- ggplot(pct_data, aes(x = species, y = percentage, fill = feature)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  labs(title = "Percentage of Rows with Feature Occurrence by Species",
       x = "Species",
       y = "Percentage of Rows (%)",
       fill = "Feature") +
  ylim(0, 100) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave(file.path(output_dir, "percentage_rows_with_occurrence.png"), 
       p3, width = 10, height = 6, dpi = 300)
cat(sprintf("Saved: %s\n", file.path(output_dir, "percentage_rows_with_occurrence.png")))

cat("\n=== Analysis Complete ===\n")
cat(sprintf("All results saved to: %s\n", output_dir))
