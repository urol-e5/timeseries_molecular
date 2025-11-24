
cat("loading libraries", "\n")
library(tidyverse)
library(ggplot2)
library(DESeq2)
library(igraph)
library(psych)
library(tidygraph)
library(ggraph)
library(WGCNA)
library(edgeR)
library(reshape2)
library(ggcorrplot)
library(corrplot)
library(rvest)
library(purrr)
library(pheatmap)
library(glmnet)
library(caret)
library(factoextra)
library(vegan)
library(ggfortify)
library(genefilter)
library(scales)
library(purrr)


# set seed
set.seed(703)

# --------- Robust command-line parsing + defaults ----------
cat("loading libraries\n")

# (your library() calls here) ...

# set seed
set.seed(703)

# read args
args <- commandArgs(trailingOnly = TRUE)
cat("Raw args (length ", length(args), "):\n")
print(args)

# Expected positional args:
# 1 genes_file
# 2 miRNA_file
# 3 lncRNA_file
# 4 WGBS_file
# 5 metadata_file
# 6 output_dir
# 7 excluded_samples_raw  (comma-separated string like "s1,s2")
# 8 alpha (optional; default 0.5)
# 9 bootstrap1_reps (optional; default 50)
#10 bootstrap2_reps (optional; default 50)
#11 r2_threshold (optional; default 0.5)

min_expected <- 6
if (length(args) < min_expected) {
  stop(sprintf("Usage: Rscript script.R <genes_file> <miRNA_file> <lncRNA_file> <WGBS_file> <metadata_file> <output_dir> [excluded_samples_raw] [alpha_value] [bootstrap1_reps] [bootstrap2_reps] [r2_threshold]\nYou provided %d args.", length(args)))
}

genes_file     <- args[1]
miRNA_file     <- args[2]
lncRNA_file    <- args[3]
WGBS_file      <- args[4]
metadata_file  <- args[5]
output_dir     <- args[6]

# -------- Optional parameters ----------
excluded_samples_raw <- if (length(args) >= 7) args[7] else NA_character_
alpha_value          <- if (length(args) >= 8) as.numeric(args[8]) else 0.5
bootstrap1_reps      <- if (length(args) >= 9) as.integer(args[9]) else 50
bootstrap2_reps      <- if (length(args) >= 10) as.integer(args[10]) else 50
r2_threshold         <- if (length(args) >= 11) as.numeric(args[11]) else 0.5

# -------- Parse excluded samples (always string → vector) --------
if (is.na(excluded_samples_raw) || excluded_samples_raw == "" || excluded_samples_raw == "NA") {
  excluded_samples <- character(0)
} else {
  excluded_samples <- unlist(strsplit(excluded_samples_raw, "[,;\\s]+"))
  excluded_samples <- excluded_samples[excluded_samples != ""]
}

# diagnostics
cat("Running EN model with:\n")
cat(" Response genes: ", genes_file, "\n")
cat(" Predictor miRNA: ", miRNA_file, "\n")
cat(" Predictor lncRNA: ", lncRNA_file, "\n")
cat(" Predictor WGBS: ", WGBS_file, "\n")
cat(" Metadata: ", metadata_file, "\n")
cat(" Output dir: ", output_dir, "\n")
cat(" Excluded samples (parsed): ", paste(excluded_samples, collapse = ", "), "\n")
cat(" alpha: ", alpha_value, "\n")
cat(" bootstrap1_reps: ", bootstrap1_reps, "\n")
cat(" bootstrap2_reps: ", bootstrap2_reps, "\n")
cat(" r2_threshold: ", r2_threshold, "\n")

# make output dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# nice to have for debugging: sessionInfo and time
cat("Session info:\n")
print(sessionInfo())
cat("Start time: ", Sys.time(), "\n")

# --------- Defining model functions ----------
cat("Defining model functions")

# Define model functions
train_models <- function(response_features, predictor_features) {
  models <- list()
  
  for (feature in colnames(response_features)) {
    y <- response_features[[feature]]  # Gene expression
    X <- as.matrix(predictor_features)  # miRNA/lncRNA/methylation as predictors
    
    # Train elastic net model (alpha = 0.5 for mix of LASSO & Ridge)
    # Reducing nfolds from 10 (default) to 5, since I'm using smaller sample sizes of just 10 samples per timepoint
    model <- cv.glmnet(X, y, alpha = alpha_value, nfolds = 3)
    
    models[[feature]] <- model
  }
  
  return(models)
}

get_feature_importance <- function(models) {
  importance_list <- lapply(models, function(model) {
    coefs <- as.matrix(coef(model, s = "lambda.min"))[-1, , drop = FALSE]  # Convert to regular matrix & remove intercept
    
    # Convert to data frame
    coefs_df <- data.frame(Feature = rownames(coefs), Importance = as.numeric(coefs))
    
    return(coefs_df)
  })
  
  # Combine feature importance across all predicted genes
  importance_df <- bind_rows(importance_list) %>%
    group_by(Feature) %>%
    summarize(MeanImportance = mean(abs(Importance)), .groups = "drop") %>%
    arrange(desc(MeanImportance))
  
  return(importance_df)
}

evaluate_model_performance <- function(models, response_features, predictor_features) {
  results <- data.frame(Feature = colnames(response_features), R2 = NA)
  
  for (feature in colnames(response_features)) {
    y <- response_features[[feature]]
    X <- as.matrix(predictor_features)
    
    model <- models[[feature]]
    preds <- predict(model, X, s = "lambda.min")
    
    R2 <- cor(y, preds)^2  # R-squared metric
    results[results$Feature == feature, "R2"] <- R2
  }
  
  return(results)
}

# Define training and testing function
train_models_split <- function(response_features, predictor_features, train_frac = 0.8, alpha = alpha_value, nfolds = 10) {
  
  # Set up storage
  models <- list()
  performance <- data.frame(Feature = colnames(response_features), R2 = NA)
  
  # Select fraction of data
  sample_idx <- sample(seq_len(nrow(predictor_features)), size = floor(train_frac * nrow(predictor_features)))
  
  # Split predictor data
  X_train <- predictor_features[sample_idx, ]
  X_test  <- predictor_features[-sample_idx, ]
  
  for (feature in colnames(response_features)) {
    # Split response data
    y_train <- response_features[sample_idx, feature]
    y_test  <- response_features[-sample_idx, feature]
    
    
    # Check if y_train is constant
    if (length(unique(y_train)) == 1) {
      warning(sprintf("Skipping feature '%s': y_train is constant", feature))
      next
    }
    
    # Run model on training data
    model <- cv.glmnet(as.matrix(X_train), y_train, alpha = alpha_value, nfolds = nfolds)
    models[[feature]] <- model
    
    # Now apply the trained model to test data
    preds <- predict(model, newx = as.matrix(X_test), s = "lambda.min")
    R2 <- cor(y_test, preds)^2
    performance[performance$Feature == feature, "R2"] <- R2
  }
  
  return(list(models = models, performance = performance))
}


## Load and format data 
### mRNA ###
# raw gene counts data (will filter and variance stabilize)
genes <- as.data.frame(read_csv(genes_file))
# format gene IDs as rownames (instead of a column)
rownames(genes) <- genes$gene_id
genes <- genes%>%select(!gene_id)

### miRNA ###
# raw miRNA counts (will filter and variance stabilize)
miRNA <- read.table(file = miRNA_file, header = TRUE, sep = "\t", check.names = FALSE)
# format miRNA IDs as rownames (instead of a column)
rownames(miRNA) <- miRNA$Name
miRNA <- miRNA%>%select(!Name)

### lncRNA ###
# raw lncRNA counts (will filter and variance stabilize)
lncRNA_full <- read.table(lncRNA_file, header = TRUE, sep = "\t", check.names = FALSE)
# Remove info on genomic location, set lncRNA IDs as rownames
rownames(lncRNA_full) <- lncRNA_full$Geneid
lncRNA <- lncRNA_full %>% select(-Geneid, -Chr, -Start, -End, -Strand, -Length)
# Format sample names to match standard (SPEC-Sample#-TP#)
colnames(lncRNA) <- gsub("\\.", "-", colnames(lncRNA))

### WGBS data ###
WGBS <- read.csv(WGBS_file)
# Format sample names to match standard (SPEC-Sample#-TP#)
colnames(WGBS) <- gsub("\\.", "-", colnames(WGBS))

### load and format metadata ###
metadata <- read_csv("../../M-multi-species/data/rna_metadata.csv")%>%select(AzentaSampleName, ColonyID, Timepoint) %>%
  filter(grepl("ACR", ColonyID))
metadata$Sample <- paste0(metadata$ColonyID, "-", metadata$Timepoint)
rownames(metadata) <- metadata$Sample
colonies <- unique(metadata$ColonyID)


# --------- Data preprocessing ----------

### Filter data sets 
genes <- genes %>% select(-any_of(excluded_samples))
miRNA <- miRNA %>% select(-any_of(excluded_samples))
lncRNA <- lncRNA %>% select(-any_of(excluded_samples))
WGBS <- WGBS %>% select(-any_of(excluded_samples))
metadata <-  metadata[!(rownames(metadata) %in% excluded_samples), ]
rownames(metadata) <- metadata$Sample

#### WGBS
# Only keep CpGs that have a non-zero value in all samples. 
WGBS_filt <- WGBS %>% filter(if_all(-X, ~ .x > 0))
# Ensure it's formatted as a data frame
WGBS_filt <- as.data.frame(WGBS_filt)
# Set CpG IDs to rownames
rownames(WGBS_filt) <- WGBS_filt$X
WGBS_filt <- WGBS_filt %>% select(-X)

cat("Number of raw WGBS sites: ", nrow(WGBS), "\n")
cat("Number of WGBS sites retained after filtering: ", nrow(WGBS_filt), "\n")

#### RNA
# Ensure we only have genes, miRNA, and lncRNA that are present in at least one sample
# genes
genes_red <- genes[rowSums(genes) != 0, ]
# miRNA
miRNA_red <- miRNA[rowSums(miRNA) != 0, ]
# lncRNA
lncRNA_red <- lncRNA[rowSums(lncRNA) != 0, ]

cat("Retained ", nrow(genes_red), " of ", nrow(genes), "genes; ",
    nrow(miRNA_red), " of ", nrow(miRNA), " miRNA; and ", 
    nrow(lncRNA_red), " of ", nrow(lncRNA), " lncRNA through presence filtering", "\n")

# *pOverA*: 
# Specifying the minimum count for a proportion of samples for each gene. Setting 3/40 = 0.08. This would retain genes that are only expressed in a single season in a couple of the colonies. Additionally, setting the minimum count so that the minimum number of samples must have a gene count above a certain threshold. 
filtering_prop <- 3/ncol(genes_red)
#genes:
filt <- filterfun(pOverA(filtering_prop, 5))
#create filter for the counts data
gfilt <- genefilter(genes_red, filt)
#identify genes to keep by count filter
gkeep <- genes_red[gfilt,]
#identify gene lists
gn.keep <- rownames(gkeep)
#gene count data filtered in PoverA, P percent of the samples have counts over A
genes_filt <- as.data.frame(genes_red[which(rownames(genes_red) %in% gn.keep),])

# miRNA:
mifilt <- filterfun(pOverA(filtering_prop, 5))
#create filter for the counts data
mifilt <- genefilter(miRNA_red, mifilt)
#identify miRNA to keep by count filter
mikeep <- miRNA_red[mifilt,]
#identify miRNA to keep by count filter
mikeep <- miRNA_red[mifilt,]
#identify miRNA lists
mi.keep <- rownames(mikeep)
#miRNA count data filtered in PoverA, P percent of the samples have counts over A
miRNA_filt <- as.data.frame(miRNA_red[which(rownames(miRNA_red) %in% mi.keep),])

# lncRNA:
lncfilt <- filterfun(pOverA(filtering_prop, 5))
#create filter for the counts data
lncfilt <- genefilter(lncRNA_red, lncfilt)
#identify lncRNA to keep by count filter
lnckeep <- lncRNA_red[lncfilt,]
#identify lncRNA to keep by count filter
lnckeep <- lncRNA_red[lncfilt,]
#identify lncRNA lists
lnc.keep <- rownames(lnckeep)
#lncRNA count data filtered in PoverA, P percent of the samples have counts over A
lncRNA_filt <- as.data.frame(lncRNA_red[which(rownames(lncRNA_red) %in% lnc.keep),])

cat("Retained ", nrow(genes_filt), " of ", nrow(genes_red), "genes; ",
    nrow(miRNA_filt), " of ", nrow(miRNA_red), " miRNA; and ", 
    nrow(lncRNA_filt), " of ", nrow(lncRNA_red), " lncRNA through pOverA filtering", "\n")


### Transform data 
#Set the order of genes, miRNA, lncRNA, wgbs, and metadata to all be the same. 
# Ensure rownames of metadata are used as the desired column order
desired_order <- rownames(metadata)

# Reorder data frame columns
genes_filt <- genes_filt[, desired_order]
miRNA_filt <- miRNA_filt[, desired_order]
lncRNA_filt <- lncRNA_filt[, desired_order]
WGBS_filt <- WGBS_filt[, desired_order]

# Use a variance stabilized transformation for all the RNA data sets. Variance stabilization essentially tries to make variance independent of the mean, and we'll implement via the differential expression analysis package `DESeq2`
### genes ###
dds_Apul <- DESeqDataSetFromMatrix(countData = genes_filt, 
                              colData = metadata, 
                              design = ~Timepoint+ColonyID)
# Variance Stabilizing Transformation
vsd_genes <- varianceStabilizingTransformation(dds_Apul, blind = TRUE)  # Must use varianceStabilizingTransformation() instead of vst() due to few input genes
vsd_genes <- assay(vsd_genes)

### miRNA ###
dds_miRNA <- DESeqDataSetFromMatrix(countData = miRNA_filt, 
                              colData = metadata, 
                              design = ~Timepoint+ColonyID)
# Variance Stabilizing Transformation
vsd_miRNA <- varianceStabilizingTransformation(dds_miRNA, blind=TRUE) # Must use varianceStabilizingTransformation() instead of vst() due to few input genes
vsd_miRNA <- assay(vsd_miRNA)

### lncRNA ###
# All values must be integers
lncRNA_filt <- lncRNA_filt %>% mutate(across(where(is.numeric), round))

dds_lncRNA <- DESeqDataSetFromMatrix(countData = lncRNA_filt, 
                              colData = metadata, 
                              design = ~Timepoint+ColonyID)
# Variance Stabilizing Transformation
vsd_lncRNA <- assay(vst(dds_lncRNA, blind = TRUE))


# For the methylation, we'll convert from beta values to M-values
### WGBS ###
# Convert beta values to M-values: M = log2(beta / (1 - beta))
# Convert percent to proportion (0-100 -> 0-1)
WGBS_prop <- WGBS_filt/100
# Protect against zeros or ones by bounding beta values (as 0s and 1s have disproportionate impact)
WGBS_prop[WGBS_prop == 0] <- 1e-6
WGBS_prop[WGBS_prop == 1] <- 1 - 1e-6
# Again, this is not a vst transformation, I'm just retaining the vsd name prefix for consistency
vsd_WGBS <- log2(WGBS_prop / (1 - WGBS_prop))

### Scale standardization
#  Elastic Net, like many regression models, requires inputs to be on the same scale
# Scale each predictor block (miRNA, lncRNA, WGBS)
vsd_miRNA_scaled <- scale(vsd_miRNA)
vsd_lncRNA_scaled <- scale(vsd_lncRNA)
vsd_WGBS_scaled <- scale(vsd_WGBS)
# Scale response (genes)
vsd_genes_scaled <- scale(vsd_genes)


### Merge predictor features
# Create a merged dataset that contains the variance stabilized counts for all miRNA, lncRNA, and CpGs.
# Triple check that all three data frames have sample names in the same order
identical(colnames(vsd_lncRNA_scaled), colnames(vsd_miRNA_scaled))
identical(colnames(vsd_lncRNA_scaled), colnames(vsd_WGBS_scaled))

# Bind (stack dataframes vertically, so that they match by column/sample)
pred_counts <- rbind(vsd_lncRNA_scaled, vsd_miRNA_scaled, vsd_WGBS_scaled)

# Transform so that samples are on rows and features are in columns
pred_counts <- t(pred_counts)
cat("Predictor set dimensions: ", dim(pred_counts), "\n")


### Format
# Transform the counts matrices so that samples are on the rows and gene IDs on the columns
vsd_genes_t <- t(vsd_genes_scaled)

# Ensure both are formatted as data frames
pred_counts <- as.data.frame(pred_counts)
vsd_genes_t <- as.data.frame(vsd_genes_t)

# Ensure sample matching between gene and epigenetic dfs
common_samples <- intersect(rownames(vsd_genes_t), rownames(pred_counts))

pred_counts <- pred_counts[common_samples, ]
vsd_genes_t <- vsd_genes_t[common_samples,]


# --------- Running the Elastic Net ----------

## Split training and test
## Bootstrapping
TM_list <- list()
MP_list <- list()
FI_list <- list()

# Perform bootstrapping
for (i in 1:bootstrap1_reps) {
  cat("Beginning replicate", i)
  cat("\n")
  result <- train_models_split(vsd_genes_t, pred_counts)
  
  cat("Extracting results for replicate", i)
  cat("\n")
  trained_models <- result$models
  model_performance <- result$performance
  feature_importance <- get_feature_importance(trained_models)
  
  TM_list[[i]] <- trained_models
  TM_list[[i]]$Replicate <- i
  
  MP_list[[i]] <- model_performance
  MP_list[[i]]$Replicate <- i
  
  FI_list[[i]] <- feature_importance
  FI_list[[i]]$Replicate <- i
}

# Combine all results
all_TM <- do.call(rbind, TM_list)
all_MP <- do.call(rbind, MP_list)
all_FI <- do.call(rbind, FI_list)

# Count how often each feature exceeds the R² threshold
MP_summary <- all_MP %>%
  group_by(Feature) %>%
  summarize(
    Mean_R2 = mean(R2, na.rm = TRUE),
    SD_R2 = sd(R2, na.rm = TRUE),
    SE_R2 = SD_R2/sqrt(n()),
    High_Perf_Count = sum(R2 >= r2_threshold),
    .groups = 'drop'
  ) %>%
  arrange(desc(High_Perf_Count))

# How many genes are consistently well-predicted?
cat("On first round of bootstrapping, ", sum(!is.na(MP_summary$Feature) & MP_summary$Mean_R2 > 0.5, na.rm = TRUE), " genes are consistently well-predicted", "\n")


# Sort genes by average R^2 across all replicates, then plot R^2 for all replicates
# Merge mean R2 info into full replicate-level data
all_MP <- merge(all_MP, MP_summary[, c("Feature", "Mean_R2", "SD_R2")], by = "Feature")

# Reorder Feature factor by Mean_R2
all_MP$Feature <- factor(all_MP$Feature,
                              levels = MP_summary$Feature[order(MP_summary$Mean_R2)])

# Plot mean R² with error bars
# Choose SD to describe variability in model performance
p <- ggplot(all_MP, aes(x = Feature, y = Mean_R2)) +
  geom_errorbar(aes(ymin = Mean_R2 - SD_R2, ymax = Mean_R2 + SD_R2), width = 0.4, color = "gray40", linewidth = 0.25) +
  geom_point(size = 2, color = "black") +
  geom_hline(yintercept = 0.5,
             linetype = "dashed", color = "blue") +
  theme_minimal() +
  labs(title = "Apul: Mean R² with Error Bars Across Gene Features",
       x = "Gene Expression Feature (Ordered by Mean R²)",
       y = "Mean R² ± SD") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5)) +
  ylim(-0.2, 1.2) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))
ggsave(paste0(output_dir, "/bootstrap1_performance_errorbars.png"), p)

# Plot with features sorted by mean R²
p <- ggplot(all_MP, aes(x = Feature, y = R2)) +
  geom_point(aes(color = as.factor(Replicate)), size = 1) +
  geom_point(aes(y = Mean_R2), color = "black", shape = 18, size = 3) + 
  geom_hline(yintercept = mean(all_MP$R2, na.rm = TRUE),
             linetype = "dashed", color = "blue") +
  theme_minimal() +
  labs(title = "Apul: Model Performance Across Gene Expression",
       x = "Gene Expression Feature (Ordered by Mean R²)",
       y = "R² (Variance Explained)",
       color = "Replicate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 1)
ggsave(paste0(output_dir, "/bootstrap1_performance.png"), p)

# Average importance of all predictive features and plot those with highest average importance across all replicates
# Average importance across replicates
mean_importance_df <- all_FI %>%
  group_by(Feature) %>%
  summarise(MeanImportance = mean(MeanImportance, na.rm = TRUE)) %>%
  arrange(desc(MeanImportance))

# Select top 50 features
top_features <- mean_importance_df %>% top_n(50, MeanImportance)

# Plot
p <- ggplot(top_features, aes(x = reorder(Feature, MeanImportance), y = MeanImportance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # Flip for readability
  theme_minimal() +
  labs(title = "Apul: Top Predictive Features by Mean Importance",
       x = "Feature",
       y = "Mean Importance")
ggsave(paste0(output_dir, "/bootstrap1_top50_features.png"), p)

## Boostrapping on predictable genes
# Isolate well-predicted genes
well_predicted <- MP_summary[MP_summary$Mean_R2 > r2_threshold,]

vsd_high_perf_t <- vsd_genes_t[,colnames(vsd_genes_t) %in% well_predicted$Feature]
dim(vsd_high_perf_t)

# Run bootstrapping on isolated genes
TM_list <- list()
MP_list <- list()
FI_list <- list()

# Perform bootstrapping
for (i in 1:bootstrap2_reps) {
  cat("Beginning replicate", i)
  cat("\n")
  result <- train_models_split(vsd_high_perf_t, pred_counts)
  
  cat("Extracting results for replicate", i)
  cat("\n")
  trained_models <-result$models
  model_performance <- result$performance
  feature_importance <- get_feature_importance(trained_models)
  
  TM_list[[i]] <- trained_models
  TM_list[[i]]$Replicate <- i
  
  MP_list[[i]] <- model_performance
  MP_list[[i]]$Replicate <- i
  
  FI_list[[i]] <- feature_importance
  FI_list[[i]]$Replicate <- i
}

# Combine all results
all_TM_highperf <- do.call(rbind, TM_list)
all_MP_highperf <- do.call(rbind, MP_list)
all_FI_highperf <- do.call(rbind, FI_list)

# Count how often each feature exceeds the R² threshold
MP_summary_highperf <- all_MP_highperf %>%
  group_by(Feature) %>%
  summarize(
    Mean_R2 = mean(R2, na.rm = TRUE),
    SD_R2 = sd(R2, na.rm = TRUE),
    SE_R2 = SD_R2 / sqrt(n()),
    High_Perf_Count = sum(R2 >= r2_threshold),
    .groups = 'drop'
  ) %>%
  arrange(desc(High_Perf_Count))

# Sort genes by average R^2 across all replicates, then plot R^2 for all replicates
# Merge mean R2 info into full replicate-level data
all_MP_highperf <- merge(all_MP_highperf, MP_summary_highperf[, c("Feature", "Mean_R2", "SD_R2")], by = "Feature")

# Reorder Feature factor by Mean_R2
all_MP_highperf$Feature <- factor(all_MP_highperf$Feature,
                                          levels = MP_summary_highperf$Feature[order(MP_summary_highperf$Mean_R2)])

# Plot mean R² with error bars
# Choose SD to describe variability in model performance
p <- ggplot(all_MP_highperf, aes(x = Feature, y = Mean_R2)) +
  geom_errorbar(aes(ymin = Mean_R2 - SD_R2, ymax = Mean_R2 + SD_R2), width = 0.3, color = "gray40") +
  geom_point(size = 2, color = "black") +
  geom_hline(yintercept = 0.5,
             linetype = "dashed", color = "blue") +
  theme_minimal() +
  labs(title = "Apul: Mean R² with Error Bars Across Gene Features",
       x = "Gene Expression Feature (Ordered by Mean R²)",
       y = "Mean R² ± SD") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(-0.2, 1.2) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))
ggsave(paste0(output_dir, "/bootstrap2_performance_errorbars.png"), p)

# Plot all replicate R2s, with features sorted by mean R²
p <- ggplot(all_MP_highperf, aes(x = Feature, y = R2)) +
  geom_point(aes(color = as.factor(Replicate)), size = 1.5) +
  geom_point(aes(y = Mean_R2), color = "black", shape = 18, size = 3) + 
  geom_hline(yintercept = mean(all_MP_highperf$R2, na.rm = TRUE),
             linetype = "dashed", color = "blue") +
  theme_minimal() +
  labs(title = "Apul: Model Performance Across Gene Expression",
       x = "Gene Expression Feature (Ordered by Mean R²)",
       y = "R² (Variance Explained)",
       color = "Replicate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 1)
ggsave(paste0(output_dir, "/bootstrap2_performance.png"), p)

# Average importance across replicates
mean_importance_df_highperf <- all_FI_highperf %>%
  group_by(Feature) %>%
  summarise(MeanImportance = mean(MeanImportance, na.rm = TRUE)) %>%
  arrange(desc(MeanImportance))

# Select top 50 features
top_features_highperf <- mean_importance_df_highperf %>% top_n(50, MeanImportance)

# Plot
p <- ggplot(top_features_highperf, aes(x = reorder(Feature, MeanImportance), y = MeanImportance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # Flip for readability
  theme_minimal() +
  labs(title = "Apul: Top Predictive Features by Mean Importance",
       x = "Feature",
       y = "Mean Importance")
ggsave(paste0(output_dir, "/bootstrap2_top50_features.png"), p)

# How many genes are consistently well-predicted?
cat("On second round of bootstrapping, ", nrow(MP_summary_highperf[MP_summary_highperf$Mean_R2 > r2_threshold,]), " genes are again consistently well-predicted", "\n")

# Look at most predictive epigenetic features for these well-predicted genes
# Define color assignments for predictor types
color_palette <- c(miRNA = "#E69F00", lncRNA = "#0072B2", CpG = "#009E73")

#  Define well-predicted features
all_features_highR2 <- all_MP_highperf %>%
  filter(Mean_R2 > 0.5) %>%
  pull(Feature) %>%
  unique()

# Function to extract importance from a single model
get_feature_importance_for_feature <- function(model) {
  coefs <- as.matrix(coef(model, s = "lambda.min"))[-1, , drop = FALSE]  # remove intercept
  data.frame(Predictor = rownames(coefs), Importance = abs(as.numeric(coefs)))
}

# Function to assign predictor type
get_predictor_type <- function(predictor_name) {
  if (startsWith(predictor_name, "Cluster")) return("miRNA")
  if (startsWith(predictor_name, "lncRNA")) return("lncRNA")
  if (startsWith(predictor_name, "CpG")) return("CpG")
  return("Other")
}

# Initialize list to save top predictors
top_predictors <- list()

# Loop over features, aggregate importance across replicates, and plot
for (target_feature in all_features_highR2) {
  
  # Extract and combine importance scores across all replicates
  importance_all_reps <- lapply(TM_list, function(rep_entry) {
    model <- rep_entry[[target_feature]]
    if (!is.null(model)) {
      get_feature_importance_for_feature(model)
    } else {
      NULL
    }
  }) %>% bind_rows()

  # If no importance data was found, skip
  if (nrow(importance_all_reps) == 0) next

  # Compute average importance across replicates
  mean_importance <- importance_all_reps %>%
    group_by(Predictor) %>%
    summarise(MeanImportance = mean(Importance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(MeanImportance)) %>%
    slice_head(n = 20) %>%
    mutate(Type = sapply(Predictor, get_predictor_type))


  # Plot once per target feature
  p <- ggplot(mean_importance, aes(x = reorder(Predictor, MeanImportance), y = MeanImportance, fill=Type)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = color_palette) +
    coord_flip() +
    theme_minimal() +
    labs(title = paste("Top Predictors for", target_feature),
         x = "Predictor",
         y = "Mean Importance (across replicates)")

  # Save top 20 most important predictors for later use
  top20 <- mean_importance %>% 
    mutate(Feature = target_feature)
  
  top_predictors[[target_feature]] <- top20
  ggsave(paste0(output_dir, "/", target_feature,"_predictors.png"), plot=p)
}



## Plot expression of genes and their predictors
# Put predictor expression in format with samples in columns
pred_counts_t <- t(pred_counts)

# Loop over each feature in the list
for (feature in names(top_predictors)) {
  
  # Get top 5 predictors for this feature
  top_df <- top_predictors[[feature]] %>%
    slice_head(n = 5)
  predictors <- top_df$Predictor
  
  # Combine feature and predictors
  all_genes <- c(feature, predictors)
  
  # Check which genes are actually present
  present_genes <- all_genes[all_genes %in% rownames(pred_counts_t) | all_genes %in% rownames(vsd_genes_scaled)]
  
  if (!(feature %in% rownames(vsd_genes_scaled))) {
    warning(paste("Feature", feature, "not found in vsd_genes_scaled Skipping..."))
    next
  }

  # Extract predictor expression (only those present)
  predictor_expr <- pred_counts_t[intersect(predictors, rownames(pred_counts_t)), , drop = FALSE] %>% as.data.frame()

  # Add the feature expression
  feature_expr <- vsd_genes_scaled[feature, , drop = FALSE] %>% as.data.frame()
  
  # Combine into one data frame
  combined_expr <- rbind(predictor_expr, feature_expr)
  combined_expr$Gene <- rownames(combined_expr)
  
  # Convert to long format
  expr_long <- combined_expr %>%
    pivot_longer(-Gene, names_to = "Sample", values_to = "Expression")
  
  # Join with sample metadata
  expr_long <- expr_long %>%
    left_join(metadata, by = "Sample")
  
  # Plot, faceted by colony
  p <- ggplot(expr_long, aes(x = Timepoint, y = Expression, color = Gene, group = Gene)) +
  geom_point(alpha = 0.7) +
  geom_smooth(se = FALSE, method = "loess", linewidth = 0.5) +
  facet_wrap(~ColonyID) +
  labs(title = paste("Expression of", feature, "and Top 5 Predictors by Colony"),
       x = "Timepoint",
       y = "Expression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(paste0(output_dir, "/", feature,"_expression.png"), plot=p)
}

# Save top predictors and metrics of each well-predicted gene
# Combine all sub-tibbles into a single data frame with the parent list name as an ID
top_predictors_df <- imap_dfr(top_predictors, ~ {
  .x
})

write.csv(top_predictors_df, paste0(output_dir, "/top_predictors.csv"), row.names = FALSE)



