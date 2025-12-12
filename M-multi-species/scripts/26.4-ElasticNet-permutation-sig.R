#!/usr/bin/env Rscript

################################################################################
# Elastic Net Regression with Permutation-Based Predictor Significance Testing
# WITH PARALLEL PROCESSING FOR BOOTSTRAPPING AND PERMUTATION STEPS
################################################################################
#
# OVERVIEW:
# This script extends the original Elastic Net approach by adding permutation-based
# significance testing for each predictor coefficient. This addresses a key limitation
# of penalized regression methods: they provide coefficient estimates but not p-values.
#
# PARALLEL PROCESSING:
# -------------------
# This version incorporates parallel processing for:
#   1. Bootstrap replicates (each replicate runs independently)
#   2. Permutation testing across genes (each gene's permutation test is independent)
#   3. Permutations within each gene (each permutation is independent)
#
# This can dramatically reduce runtime when multiple cores are available.
# On a 4-core machine, expect ~3-4x speedup for bootstrap and permutation steps.
#
# STATISTICAL RATIONALE:
# ----------------------
# Elastic Net (and LASSO/Ridge) regression performs variable selection through
# regularization, but the resulting coefficients lack traditional inference statistics.
# This is because:
#   1. The penalization introduces bias into coefficient estimates
#   2. Standard asymptotic theory doesn't apply to penalized estimators
#   3. The selection process invalidates classical hypothesis testing
#
# PERMUTATION APPROACH:
# --------------------
# To obtain valid p-values, we use a permutation test framework:
#   1. For each gene, fit the real Elastic Net model and record predictor coefficients
#   2. Then, permute the response variable (gene expression) many times
#   3. For each permutation, refit the model and record coefficients
#   4. The p-value for each predictor is the proportion of permuted coefficients
#      that are >= the observed coefficient (in absolute value)
#
# REFERENCES:
# ----------
# - Phipson & Smyth (2010). "Permutation P-values Should Never Be Zero"
#   Statistical Applications in Genetics and Molecular Biology. DOI: 10.2202/1544-6115.1585
# - Ojala & Garriga (2010). "Permutation Tests for Studying Classifier Performance"
#   Journal of Machine Learning Research. DOI: 10.5555/1756006.1859920
# - Altmann et al. (2010). "Permutation importance: a corrected feature importance measure"
#   Bioinformatics. DOI: 10.1093/bioinformatics/btq134
#
################################################################################

cat("=== Elastic Net with Permutation Significance Testing (Parallel) ===\n")
cat("Loading libraries\n")

suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(DESeq2)))
suppressMessages(suppressWarnings(library(igraph)))
suppressMessages(suppressWarnings(library(psych)))
suppressMessages(suppressWarnings(library(tidygraph)))
suppressMessages(suppressWarnings(library(ggraph)))
suppressMessages(suppressWarnings(library(WGCNA)))
suppressMessages(suppressWarnings(library(edgeR)))
suppressMessages(suppressWarnings(library(reshape2)))
suppressMessages(suppressWarnings(library(ggcorrplot)))
suppressMessages(suppressWarnings(library(corrplot)))
suppressMessages(suppressWarnings(library(rvest)))
suppressMessages(suppressWarnings(library(purrr)))
suppressMessages(suppressWarnings(library(pheatmap)))
suppressMessages(suppressWarnings(library(glmnet)))
suppressMessages(suppressWarnings(library(caret)))
suppressMessages(suppressWarnings(library(factoextra)))
suppressMessages(suppressWarnings(library(vegan)))
suppressMessages(suppressWarnings(library(ggfortify)))
suppressMessages(suppressWarnings(library(genefilter)))
suppressMessages(suppressWarnings(library(scales)))
suppressMessages(suppressWarnings(library(parallel)))  # For parallel processing
suppressMessages(suppressWarnings(library(doParallel))) # Additional parallel support

# Set seed for reproducibility
set.seed(703)

################################################################################
# LOGGING SYSTEM SETUP
################################################################################
# 
# This logging system writes all output to both:
#   1. The terminal (stdout/stderr) for real-time monitoring
#   2. A log file in the output directory for record-keeping
#
# We achieve this by:
#   1. Creating a log file connection
#   2. Using a custom tee_output() function that writes to both destinations
#   3. Wrapping key output functions to use tee behavior
#

# Initialize logging variables (will be set up after output_dir is known)
log_file_path <- NULL
log_file_con <- NULL

#' Initialize the logging system
#' Call this after output_dir is created
#' 
#' @param output_dir Directory where log file will be saved
#' @param prefix Prefix for log filename
#' @return Path to the log file
init_logging <- function(output_dir, prefix = "analysis") {
  # Create timestamped log filename
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  log_filename <- sprintf("%s_log_%s.txt", prefix, timestamp)
  log_path <- file.path(output_dir, log_filename)
  
  # Open log file connection for appending
  log_file_con <<- file(log_path, open = "wt")
  log_file_path <<- log_path
  
  # Write header to log file
  writeLines(sprintf("=== Log initialized: %s ===", log_path), log_file_con)
  writeLines(sprintf("Timestamp: %s\n", Sys.time()), log_file_con)
  flush(log_file_con)
  
  # Print to console too
  cat(sprintf("=== Log initialized: %s ===\n", log_path))
  cat(sprintf("Timestamp: %s\n\n", Sys.time()))
  
  return(log_path)
}

#' Write message to both console and log file
#' Use this instead of cat() for important messages
#' 
#' @param ... Arguments to pass to cat()
log_message <- function(...) {
  msg <- paste0(...)
  
  # Write to console
  cat(msg)
  
  
  # Write to log file if it's open
  if (!is.null(log_file_con) && isOpen(log_file_con)) {
    writeLines(gsub("\n$", "", msg), log_file_con)
    flush(log_file_con)
  }
}

#' Capture and log the output of print() calls
#' 
#' @param x Object to print
log_print <- function(x) {
  # Capture print output
  output <- capture.output(print(x))
  
  # Write to console
  cat(paste(output, collapse = "\n"), "\n")
  
  # Write to log file
  if (!is.null(log_file_con) && isOpen(log_file_con)) {
    writeLines(output, log_file_con)
    flush(log_file_con)
  }
}

#' Close the logging system
#' Call this at the end of the script
close_logging <- function() {
  if (!is.null(log_file_con)) {
    tryCatch({
      if (isOpen(log_file_con)) {
        writeLines(sprintf("\n=== Log closed: %s ===", Sys.time()), log_file_con)
        flush(log_file_con)
        close(log_file_con)
      }
    }, error = function(e) {
      # Silently ignore errors during cleanup
    })
    log_file_con <<- NULL
  }
}

#' Wrapper for cat that also logs to file
#' Override the default cat to enable dual output
tee_cat <- function(..., file = "", sep = " ", fill = FALSE, labels = NULL, append = FALSE) {
  # Call original cat for console output
  base::cat(..., file = file, sep = sep, fill = fill, labels = labels, append = append)
  
  # Also write to log file if logging is active and not writing to another file
  if (file == "" && !is.null(log_file_con) && isOpen(log_file_con)) {
    msg <- paste(..., sep = sep)
    writeLines(gsub("\n$", "", msg), log_file_con)
    flush(log_file_con)
  }
}

#' Wrapper for print that also logs to file  
tee_print <- function(x, ...) {
  # Capture print output
  output <- capture.output(base::print(x, ...))
  
  # Print to console
  base::cat(paste(output, collapse = "\n"), "\n")
  
  # Write to log file
  if (!is.null(log_file_con) && isOpen(log_file_con)) {
    writeLines(output, log_file_con)
    flush(log_file_con)
  }
  
  invisible(x)
}

# Store original functions before overriding
original_cat <- base::cat
original_print <- base::print

################################################################################
# PROGRESS TRACKING UTILITIES
################################################################################
#
# These functions provide progress reporting for long-running operations.
# 
# For parallel operations, real-time progress is challenging in R because
# mclapply child processes can't communicate back to the parent. We use two
# strategies:
#   1. If pbapply package is available, use pbmclapply for built-in progress
#   2. Otherwise, use a shared file counter with a monitoring thread
#

# Check if pbapply is available for nicer parallel progress bars
has_pbapply <- requireNamespace("pbapply", quietly = TRUE)
if (has_pbapply) {
  suppressMessages(library(pbapply))
  cat("Using pbapply for progress bars\n")
} else {
  cat("Note: Install 'pbapply' package for better parallel progress bars\n")
  cat("      install.packages('pbapply')\n")
}

#' Create a text-based progress bar string
#' 
#' @param current Current progress count
#' @param total Total count
#' @param width Width of the progress bar in characters
#' @param prefix Text to show before the bar
#' @return Formatted progress string
format_progress_bar <- function(current, total, width = 40, prefix = "Progress") {
  if (total == 0) return(sprintf("%s: 0/0 (0%%)", prefix))
  pct <- current / total
  filled <- round(pct * width)
  empty <- width - filled
  bar <- paste0("[", paste(rep("=", filled), collapse = ""), 
                ifelse(filled < width, ">", ""),
                paste(rep(" ", max(0, empty - 1)), collapse = ""), "]")
  sprintf("\r%s: %s %d/%d (%.1f%%)", prefix, bar, current, total, pct * 100)
}

#' Print progress update (overwrites current line for clean display)
#' 
#' @param current Current progress count
#' @param total Total count
#' @param prefix Text to show before the bar
#' @param newline Whether to add newline at end (use TRUE when complete)
print_progress <- function(current, total, prefix = "Progress", newline = FALSE) {
  msg <- format_progress_bar(current, total, prefix = prefix)
  if (newline) {
    base::cat(msg, "\n")
  } else {
    base::cat(msg)
    flush.console()  # Force immediate display
  }
  # Also write to log (without carriage return, as a simple status line)
  if (!is.null(log_file_con) && isOpen(log_file_con)) {
    if (current == total || current %% max(1, floor(total / 10)) == 0) {
      # Only log every 10% or at completion to avoid huge log files
      writeLines(sprintf("%s: %d/%d (%.1f%%)", prefix, current, total, current/total*100), 
                 log_file_con)
      flush(log_file_con)
    }
  }
}

#' Log a progress milestone to file only
#' 
#' @param current Current count
#' @param total Total count
#' @param prefix Prefix for the message
log_progress_milestone <- function(current, total, prefix = "Progress") {
  if (!is.null(log_file_con) && isOpen(log_file_con)) {
    writeLines(sprintf("%s: %d/%d (%.1f%%)", prefix, current, total, 
                       ifelse(total > 0, current/total*100, 0)), 
               log_file_con)
    flush(log_file_con)
  }
}

################################################################################
# COMMAND LINE ARGUMENTS
################################################################################

args <- commandArgs(trailingOnly = TRUE)
cat("Raw args (length ", length(args), "):\n")
print(args)

# Expected positional args:
# 1  genes_file
# 2  miRNA_file
# 3  lncRNA_file
# 4  WGBS_file
# 5  metadata_file
# 6  output_dir
# 7  excluded_samples_raw  (comma-separated string like "s1,s2")
# 8  alpha (optional; default 0.5)
# 9  bootstrap1_reps (optional; default 50)
# 10 bootstrap2_reps (optional; default 50)
# 11 r2_threshold (optional; default 0.5)
# 12 n_permutations (optional; default 1000)
# 13 n_cores (optional; default 1) -- for parallel processing

min_expected <- 6
if (length(args) < min_expected) {
  stop(sprintf("Usage: Rscript script.R <genes_file> <miRNA_file> <lncRNA_file> <WGBS_file> <metadata_file> <output_dir> [excluded_samples_raw] [alpha_value] [bootstrap1_reps] [bootstrap2_reps] [r2_threshold] [n_permutations] [n_cores]\nYou provided %d args.", length(args)))
}

genes_file     <- args[1]
miRNA_file     <- args[2]
lncRNA_file    <- args[3]
WGBS_file      <- args[4]
metadata_file  <- args[5]
output_dir     <- args[6]

# -------- Optional parameters with safe parsing ----------
# Helper functions to safely parse arguments
safe_numeric <- function(x, default) {
  if (is.na(x) || x == "" || x == "NA" || trimws(x) == "") return(default)
  val <- suppressWarnings(as.numeric(trimws(x)))
  if (is.na(val)) return(default)
  return(val)
}

safe_integer <- function(x, default) {
  if (is.na(x) || x == "" || x == "NA" || trimws(x) == "") return(default)
  val <- suppressWarnings(as.integer(trimws(x)))
  if (is.na(val)) return(default)
  return(val)
}

excluded_samples_raw <- if (length(args) >= 7)  args[7]                          else NA_character_
alpha_value          <- if (length(args) >= 8)  safe_numeric(args[8], 0.5)       else 0.5
bootstrap1_reps      <- if (length(args) >= 9)  safe_integer(args[9], 50)        else 50
bootstrap2_reps      <- if (length(args) >= 10) safe_integer(args[10], 50)       else 50
r2_threshold         <- if (length(args) >= 11) safe_numeric(args[11], 0.5)      else 0.5
n_permutations       <- if (length(args) >= 12) safe_integer(args[12], 1000)     else 1000
n_cores              <- if (length(args) >= 13) safe_integer(args[13], 1)        else 1

# Validate critical parameters
if (is.na(n_permutations) || n_permutations < 1) {
  cat("Warning: Invalid n_permutations value, using default of 1000\n")
  n_permutations <- 1000
}
if (is.na(n_cores) || n_cores < 1) {
  cat("Warning: Invalid n_cores value, using default of 1\n")
  n_cores <- 1
}

# Cap n_cores at available cores
max_cores <- detectCores()
if (n_cores > max_cores) {
  cat(sprintf("Warning: Requested %d cores but only %d available. Using %d cores.\n", 
              n_cores, max_cores, max_cores))
  n_cores <- max_cores
}

# Parse excluded samples
if (is.na(excluded_samples_raw) || excluded_samples_raw == "" || excluded_samples_raw == "NA") {
  excluded_samples <- character(0)
} else {
  excluded_samples <- unlist(strsplit(excluded_samples_raw, "[,;\\s]+"))
  excluded_samples <- excluded_samples[excluded_samples != ""]
}

cat("List of samples to exclude:", excluded_samples, "\n")

# Diagnostics
cat("\nRunning EN model with permutation significance testing:\n")
cat(" Response genes:     ", genes_file, "\n")
cat(" Predictor miRNA:    ", miRNA_file, "\n")
cat(" Predictor lncRNA:   ", lncRNA_file, "\n")
cat(" Predictor WGBS:     ", WGBS_file, "\n")
cat(" Metadata:           ", metadata_file, "\n")
cat(" Output dir:         ", output_dir, "\n")
cat(" Excluded samples:   ", paste(excluded_samples, collapse = ", "), "\n")
cat(" Alpha:              ", alpha_value, "\n")
cat(" Bootstrap1 reps:    ", bootstrap1_reps, "\n")
cat(" Bootstrap2 reps:    ", bootstrap2_reps, "\n")
cat(" R2 threshold:       ", r2_threshold, "\n")
cat(" N permutations:     ", n_permutations, "\n")
cat(" N cores:            ", n_cores, "\n")
cat(" Max available cores:", max_cores, "\n")

# Create output directories
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

################################################################################
# INITIALIZE LOGGING
################################################################################
# Now that we have the output directory, initialize logging to file + console
log_file_path <- init_logging(output_dir, prefix = "ElasticNet_permutation")

# Override cat and print to use tee versions (dual output to console + log)
# This ensures all subsequent cat() and print() calls go to both destinations
cat <- tee_cat
print <- tee_print

cat("\nSession info:\n")
print(sessionInfo())
cat("Start time:", as.character(Sys.time()), "\n\n")

################################################################################
# FUNCTION DEFINITIONS
################################################################################

cat("Defining model functions\n")

#' Train Elastic Net Models for All Response Features
#' 
#' @param response_features Data frame with response variables (genes) in columns
#' @param predictor_features Data frame with predictor variables in columns
#' @return List of cv.glmnet model objects, one per gene
train_models <- function(response_features, predictor_features, alpha = alpha_value) {
  models <- list()
  
  for (feature in colnames(response_features)) {
    y <- response_features[[feature]]
    X <- as.matrix(predictor_features)
    
    model <- cv.glmnet(X, y, alpha = alpha, nfolds = 3)
    models[[feature]] <- model
  }
  
  return(models)
}

#' Extract Feature Importance (Coefficient Magnitudes) from Models
#' 
#' @param models List of cv.glmnet models
#' @return Data frame with Feature, MeanImportance columns
get_feature_importance <- function(models) {
  importance_list <- lapply(models, function(model) {
    coefs <- as.matrix(coef(model, s = "lambda.min"))[-1, , drop = FALSE]
    coefs_df <- data.frame(Feature = rownames(coefs), Importance = as.numeric(coefs))
    return(coefs_df)
  })
  
  importance_df <- bind_rows(importance_list) %>%
    group_by(Feature) %>%
    summarize(MeanImportance = mean(abs(Importance)), .groups = "drop") %>%
    arrange(desc(MeanImportance))
  
  return(importance_df)
}

#' Evaluate Model Performance (R-squared)
#' 
#' @param models List of cv.glmnet models
#' @param response_features Data frame with response variables
#' @param predictor_features Data frame with predictor variables
#' @return Data frame with Feature, R2 columns
evaluate_model_performance <- function(models, response_features, predictor_features) {
  results <- data.frame(Feature = colnames(response_features), R2 = NA)
  
  for (feature in colnames(response_features)) {
    y <- response_features[[feature]]
    X <- as.matrix(predictor_features)
    
    model <- models[[feature]]
    preds <- predict(model, X, s = "lambda.min")
    
    R2 <- cor(y, preds)^2
    results[results$Feature == feature, "R2"] <- R2
  }
  
  return(results)
}

#' Train Models with Train/Test Split (for a single bootstrap replicate)
#' 
#' @param response_features Data frame with response variables
#' @param predictor_features Data frame with predictor variables
#' @param train_frac Fraction of data for training
#' @param alpha Elastic net mixing parameter
#' @param nfolds Number of CV folds
#' @param seed Optional seed for this replicate
#' @return List with 'models' and 'performance' components
train_models_split <- function(response_features, predictor_features, 
                               train_frac = 0.8, alpha = alpha_value, 
                               nfolds = 10, seed = NULL) {
  
  # Set seed if provided (important for parallel reproducibility)
  if (!is.null(seed)) set.seed(seed)
  
  models <- list()
  performance <- data.frame(Feature = colnames(response_features), R2 = NA)
  
  # Sample training indices
  sample_idx <- sample(seq_len(nrow(predictor_features)), 
                       size = floor(train_frac * nrow(predictor_features)))
  
  X_train <- predictor_features[sample_idx, ]
  X_test  <- predictor_features[-sample_idx, ]
  
  for (feature in colnames(response_features)) {
    y_train <- response_features[sample_idx, feature]
    y_test  <- response_features[-sample_idx, feature]
    
    # Skip constant response variables
    if (length(unique(y_train)) == 1) {
      warning(sprintf("Skipping feature '%s': y_train is constant", feature))
      next
    }
    
    model <- tryCatch({
      cv.glmnet(as.matrix(X_train), y_train, alpha = alpha, nfolds = nfolds)
    }, error = function(e) {
      warning(sprintf("Error fitting model for %s: %s", feature, e$message))
      return(NULL)
    })
    
    if (!is.null(model)) {
      models[[feature]] <- model
      preds <- predict(model, newx = as.matrix(X_test), s = "lambda.min")
      R2 <- cor(y_test, preds)^2
      performance[performance$Feature == feature, "R2"] <- R2
    }
  }
  
  return(list(models = models, performance = performance))
}


################################################################################
# PARALLEL BOOTSTRAP FUNCTION (NEW)
################################################################################

#' Run Bootstrap Replicates in Parallel
#' 
#' This function distributes bootstrap replicates across multiple cores.
#' Each replicate gets a unique seed derived from the base seed + replicate number
#' to ensure reproducibility while maintaining independence.
#'
#' @param response_features Data frame with response variables
#' @param predictor_features Data frame with predictor variables
#' @param n_reps Number of bootstrap replicates
#' @param n_cores Number of cores to use
#' @param base_seed Base seed for reproducibility
#' @param alpha Elastic net alpha parameter
#' @param progress_prefix Prefix for progress display
#' @return List with TM_list, MP_list, FI_list
#' 
run_bootstrap_parallel <- function(response_features, predictor_features, 
                                   n_reps, n_cores, base_seed = 703, 
                                   alpha = alpha_value,
                                   progress_prefix = "Bootstrap") {
  
  cat(sprintf("Running %d bootstrap replicates using %d core(s)\n", n_reps, n_cores))
  
  # Function to run a single bootstrap replicate
  run_single_bootstrap <- function(i) {
    # Set unique seed for this replicate
    rep_seed <- base_seed + i
    
    result <- train_models_split(response_features, predictor_features, 
                                 alpha = alpha, seed = rep_seed)
    
    trained_models <- result$models
    model_performance <- result$performance
    feature_importance <- get_feature_importance(trained_models)
    
    # Add replicate ID
    model_performance$Replicate <- i
    feature_importance$Replicate <- i
    
    return(list(
      models = trained_models,
      performance = model_performance,
      importance = feature_importance,
      replicate = i
    ))
  }
  
  # Run in parallel or sequentially
  if (n_cores > 1) {
    cat(sprintf("%s: Running in parallel mode...\n", progress_prefix))
    log_progress_milestone(0, n_reps, progress_prefix)
    
    if (has_pbapply) {
      # Use pbapply for nice progress bars
      pboptions(type = "txt", style = 3, char = "=")
      results <- pblapply(1:n_reps, run_single_bootstrap, cl = n_cores)
    } else {
      # Fall back to mclapply without real-time progress
      # Show start message and completion
      cat(sprintf("%s: Processing %d replicates (progress shown at completion)...\n", 
                  progress_prefix, n_reps))
      results <- mclapply(1:n_reps, run_single_bootstrap, mc.cores = n_cores)
      cat(sprintf("%s: Completed %d/%d replicates\n", progress_prefix, n_reps, n_reps))
    }
    
    log_progress_milestone(n_reps, n_reps, progress_prefix)
    
  } else {
    cat(sprintf("%s: Running in sequential mode...\n", progress_prefix))
    results <- vector("list", n_reps)
    
    for (i in 1:n_reps) {
      print_progress(i - 1, n_reps, prefix = progress_prefix)
      results[[i]] <- run_single_bootstrap(i)
    }
    print_progress(n_reps, n_reps, prefix = progress_prefix, newline = TRUE)
  }
  
  # Extract and organize results
  TM_list <- lapply(results, function(x) {
    if (!is.null(x)) {
      x$models$Replicate <- x$replicate
      x$models
    } else {
      NULL
    }
  })
  
  MP_list <- lapply(results, function(x) if (!is.null(x)) x$performance else NULL)
  FI_list <- lapply(results, function(x) if (!is.null(x)) x$importance else NULL)
  
  # Remove NULLs
  TM_list <- Filter(Negate(is.null), TM_list)
  MP_list <- Filter(Negate(is.null), MP_list)
  FI_list <- Filter(Negate(is.null), FI_list)
  
  return(list(
    TM_list = TM_list,
    MP_list = MP_list,
    FI_list = FI_list
  ))
}


################################################################################
# PERMUTATION SIGNIFICANCE FUNCTIONS (WITH PARALLEL SUPPORT)
################################################################################

#' Calculate Permutation-Based P-values for Elastic Net Coefficients
#' 
#' This function implements a permutation test to assess the statistical significance
#' of each predictor's coefficient in an Elastic Net model.
#' 
#' IMPORTANT: Only predictors with non-zero observed coefficients are tested.
#' Predictors with zero coefficients are assigned p-value = NA because:
#'   1. Testing whether 0 differs from a null distribution is not meaningful
#'   2. Elastic Net has already decided these predictors are not important
#'   3. This dramatically reduces the multiple testing burden
#'
#' @param y Response vector (gene expression)
#' @param X Predictor matrix
#' @param n_perm Number of permutations
#' @param alpha Elastic net mixing parameter
#' @param n_cores Number of cores for parallel permutations
#' @param seed Random seed for reproducibility
#' @param test_only_nonzero If TRUE (default), only test predictors with non-zero coefficients
#' @param zero_threshold Threshold below which coefficients are considered zero
#' @return Data frame with columns: Predictor, Observed_Coef, Perm_Mean, Perm_SD, P_value
#' 
calculate_permutation_pvalues <- function(y, X, n_perm = 1000, alpha = 0.5, 
                                          n_cores = 1, seed = NULL,
                                          test_only_nonzero = TRUE,
                                          zero_threshold = 1e-10) {
  
  if (!is.null(seed)) set.seed(seed)
  
  X <- as.matrix(X)
  n_predictors <- ncol(X)
  predictor_names <- colnames(X)
  
  # Step 1: Fit the observed (real) model
  observed_model <- tryCatch({
    cv.glmnet(X, y, alpha = alpha, nfolds = 3)
  }, error = function(e) {
    warning("Error fitting observed model: ", e$message)
    return(NULL)
  })
  
  if (is.null(observed_model)) {
    return(data.frame(
      Predictor = predictor_names,
      Observed_Coef = NA,
      Abs_Observed_Coef = NA,
      Perm_Mean = NA,
      Perm_SD = NA,
      P_value = NA,
      Tested = FALSE,
      stringsAsFactors = FALSE
    ))
  }
  
  # Extract observed coefficients (excluding intercept)
  observed_coefs <- as.numeric(coef(observed_model, s = "lambda.min")[-1])
  names(observed_coefs) <- predictor_names
  
  # Store the lambda used for the observed model (for potential fixed-lambda approach)
  observed_lambda <- observed_model$lambda.min
  
  # Identify which predictors to test (non-zero coefficients only)
  nonzero_mask <- abs(observed_coefs) > zero_threshold
  n_nonzero <- sum(nonzero_mask)
  
  # If no non-zero coefficients, return early
  if (n_nonzero == 0) {
    return(data.frame(
      Predictor = predictor_names,
      Observed_Coef = observed_coefs,
      Abs_Observed_Coef = abs(observed_coefs),
      Perm_Mean = NA,
      Perm_SD = NA,
      P_value = NA,
      Tested = FALSE,
      stringsAsFactors = FALSE
    ))
  }
  
  # Get indices and names of predictors to test
  if (test_only_nonzero) {
    test_indices <- which(nonzero_mask)
  } else {
    test_indices <- seq_len(n_predictors)
  }
  
  # Step 2: Perform permutations
  # Function for a single permutation - using FIXED lambda from observed model
  # This ensures we're comparing apples to apples
  run_single_permutation <- function(perm_idx) {
    # Generate unique seed for this permutation
    if (!is.null(seed)) set.seed(seed + perm_idx)
    
    y_perm <- sample(y)
    
    perm_model <- tryCatch({
      # Use glmnet with fixed lambda (not cv.glmnet) for consistency
      glmnet(X, y_perm, alpha = alpha, lambda = observed_lambda)
    }, error = function(e) {
      return(NULL)
    })
    
    if (!is.null(perm_model)) {
      # Extract coefficients at the fixed lambda
      perm_coefs <- as.numeric(coef(perm_model, s = observed_lambda)[-1])
      return(perm_coefs[test_indices])
    } else {
      return(rep(NA, length(test_indices)))
    }
  }
  
  # Run permutations
  if (n_cores > 1 && n_perm >= n_cores) {
    perm_results <- mclapply(1:n_perm, run_single_permutation, mc.cores = n_cores)
  } else {
    perm_results <- lapply(1:n_perm, run_single_permutation)
  }
  
  # Convert to matrix (rows = permutations, cols = tested predictors)
  perm_coefs <- do.call(rbind, perm_results)
  
  # Step 3: Calculate p-values for tested predictors only
  p_values <- rep(NA, n_predictors)
  perm_means <- rep(NA, n_predictors)
  perm_sds <- rep(NA, n_predictors)
  tested <- rep(FALSE, n_predictors)
  
  for (idx in seq_along(test_indices)) {
    j <- test_indices[idx]
    obs_coef <- abs(observed_coefs[j])
    perm_coef_vec <- abs(perm_coefs[, idx])
    perm_coef_vec <- perm_coef_vec[!is.na(perm_coef_vec)]
    
    tested[j] <- TRUE
    
    if (length(perm_coef_vec) > 0) {
      n_extreme <- sum(perm_coef_vec >= obs_coef)
      p_values[j] <- (n_extreme + 1) / (length(perm_coef_vec) + 1)
      perm_means[j] <- mean(perm_coef_vec)
      perm_sds[j] <- sd(perm_coef_vec)
    }
  }
  
  results <- data.frame(
    Predictor = predictor_names,
    Observed_Coef = observed_coefs,
    Abs_Observed_Coef = abs(observed_coefs),
    Perm_Mean = perm_means,
    Perm_SD = perm_sds,
    P_value = p_values,
    Tested = tested,
    stringsAsFactors = FALSE
  )
  
  return(results)
}


#' Wrapper to Calculate Permutation P-values for All Genes (Parallelized)
#' 
#' This function distributes genes across cores for parallel processing.
#' Each gene's permutation test runs independently.
#' 
#' Only predictors with non-zero coefficients are tested for significance.
#' This dramatically reduces the multiple testing burden and focuses on
#' biologically relevant associations.
#'
#' @param response_features Data frame with genes in columns
#' @param predictor_features Data frame with predictors in columns  
#' @param n_perm Number of permutations per gene
#' @param alpha Elastic net alpha parameter
#' @param n_cores Number of cores for parallel processing
#' @param base_seed Base seed for reproducibility
#' @param progress_prefix Prefix for progress display
#' @return Data frame with Gene, Predictor, Observed_Coef, P_value, etc.
#' 
calculate_all_permutation_pvalues <- function(response_features, 
                                              predictor_features,
                                              n_perm = 1000,
                                              alpha = 0.5,
                                              n_cores = 1,
                                              base_seed = 703,
                                              progress_prefix = "Permutation") {
  
  gene_names <- colnames(response_features)
  n_genes <- length(gene_names)
  n_predictors <- ncol(predictor_features)
  X <- as.matrix(predictor_features)
  
  cat(sprintf("\nCalculating permutation p-values for %d genes with %d permutations each\n", 
              n_genes, n_perm))
  cat(sprintf("Total predictors per gene: %d\n", n_predictors))
  cat(sprintf("Using %d core(s) for gene-level parallelization\n", n_cores))
  cat(sprintf("NOTE: Only non-zero coefficients will be tested (reduces multiple testing burden)\n"))
  
  # Function to process a single gene
  process_gene <- function(gene_idx) {
    gene_name <- gene_names[gene_idx]
    y <- response_features[[gene_name]]
    
    # Skip if response is constant
    if (length(unique(y)) <= 1) {
      return(NULL)
    }
    
    # Unique seed per gene for reproducibility
    gene_seed <- base_seed + gene_idx * 1000
    
    # Note: we use n_cores=1 here because we're already parallelizing across genes
    # This avoids nested parallelization which can cause issues
    result <- calculate_permutation_pvalues(y, X, n_perm = n_perm, 
                                            alpha = alpha, n_cores = 1, 
                                            seed = gene_seed,
                                            test_only_nonzero = TRUE)
    result$Gene <- gene_name
    
    return(result)
  }
  
  # Process genes in parallel or sequentially
  if (n_cores > 1) {
    cat(sprintf("%s: Running in parallel mode across %d genes...\n", progress_prefix, n_genes))
    log_progress_milestone(0, n_genes, progress_prefix)
    
    if (has_pbapply) {
      # Use pbapply for nice progress bars
      pboptions(type = "txt", style = 3, char = "=")
      results_list <- pblapply(1:n_genes, process_gene, cl = n_cores)
    } else {
      # Fall back to mclapply without real-time progress
      cat(sprintf("%s: Processing %d genes (progress shown at completion)...\n", 
                  progress_prefix, n_genes))
      results_list <- mclapply(1:n_genes, process_gene, mc.cores = n_cores)
      cat(sprintf("%s: Completed %d/%d genes\n", progress_prefix, n_genes, n_genes))
    }
    
    log_progress_milestone(n_genes, n_genes, progress_prefix)
    
  } else {
    cat(sprintf("%s: Running in sequential mode...\n", progress_prefix))
    results_list <- vector("list", n_genes)
    
    for (i in seq_along(gene_names)) {
      # Update progress bar
      print_progress(i - 1, n_genes, prefix = progress_prefix)
      results_list[[i]] <- process_gene(i)
    }
    print_progress(n_genes, n_genes, prefix = progress_prefix, newline = TRUE)
  }
  
  # Combine results
  all_results <- bind_rows(results_list)
  
  # Print diagnostic summary
  n_tested <- sum(all_results$Tested, na.rm = TRUE)
  n_total <- nrow(all_results)
  n_nonzero <- sum(all_results$Abs_Observed_Coef > 1e-10, na.rm = TRUE)
  
  cat(sprintf("\n=== Permutation Test Summary ===\n"))
  cat(sprintf("Total predictor-gene combinations: %d\n", n_total))
  cat(sprintf("Non-zero coefficients: %d (%.2f%%)\n", n_nonzero, 100*n_nonzero/n_total))
  cat(sprintf("Predictors tested: %d\n", n_tested))
  cat(sprintf("Predictors skipped (zero coef): %d\n", n_total - n_tested))
  
  # Reorder columns
  all_results <- all_results %>%
    select(Gene, Predictor, Observed_Coef, Abs_Observed_Coef, 
           Perm_Mean, Perm_SD, P_value, Tested)
  
  return(all_results)
}


#' Apply Multiple Testing Correction (FDR)
#' 
#' FDR correction is applied ONLY to tested predictors (those with non-zero coefficients).
#' This is the correct approach because:
#'   1. We only performed hypothesis tests on non-zero coefficients
#'   2. Including untested predictors would artificially inflate the correction
#'   3. The multiple testing burden is already greatly reduced by testing only non-zero coefficients
#' 
#' @param perm_results Data frame from calculate_all_permutation_pvalues
#' @return Data frame with additional P_adj column
#' 
apply_fdr_correction <- function(perm_results) {
  
  # Apply FDR correction only to tested predictors, within each gene
  perm_results <- perm_results %>%
    group_by(Gene) %>%
    mutate(
      # Only adjust p-values for predictors that were actually tested
      P_adj = ifelse(Tested, 
                     p.adjust(P_value[Tested], method = "BH")[cumsum(Tested)],
                     NA)
    ) %>%
    ungroup()
  
  # Actually, the above is tricky. Let's do it more explicitly:
  perm_results <- perm_results %>%
    group_by(Gene) %>%
    mutate(
      P_adj = {
        # Get indices of tested predictors
        tested_idx <- which(Tested)
        # Initialize with NA
        padj <- rep(NA_real_, n())
        # Apply FDR only to tested predictors
        if (length(tested_idx) > 0) {
          padj[tested_idx] <- p.adjust(P_value[tested_idx], method = "BH")
        }
        padj
      }
    ) %>%
    ungroup()
  
  # Report correction statistics
  n_tested <- sum(perm_results$Tested, na.rm = TRUE)
  n_genes <- n_distinct(perm_results$Gene)
  avg_tested_per_gene <- n_tested / n_genes
  
  cat(sprintf("\nFDR correction applied to %d tests across %d genes\n", n_tested, n_genes))
  cat(sprintf("Average tests per gene: %.1f (vs. %d total predictors)\n", 
              avg_tested_per_gene, n_distinct(perm_results$Predictor)))
  
  return(perm_results)
}


#' Summarize Significant Predictors
#' 
#' @param perm_results Data frame with permutation test results
#' @param fdr_threshold FDR threshold for significance
#' @return List with summary statistics and significant predictor table
#' 
summarize_significant_predictors <- function(perm_results, fdr_threshold = 0.05) {
  
  # Only consider tested predictors for significance
  significant <- perm_results %>%
    filter(Tested, P_adj < fdr_threshold) %>%
    arrange(Gene, P_adj)
  
  gene_summary <- significant %>%
    group_by(Gene) %>%
    summarize(
      N_significant = n(),
      Top_predictor = Predictor[which.min(P_adj)],
      Min_padj = min(P_adj),
      .groups = "drop"
    )
  
  predictor_summary <- significant %>%
    mutate(
      Predictor_Type = case_when(
        grepl("^Cluster", Predictor) ~ "miRNA",
        grepl("^lncRNA", Predictor) ~ "lncRNA",
        grepl("^CpG", Predictor) ~ "CpG",
        TRUE ~ "Other"
      )
    ) %>%
    group_by(Predictor_Type) %>%
    summarize(
      N_significant = n(),
      N_unique_predictors = n_distinct(Predictor),
      N_genes_affected = n_distinct(Gene),
      .groups = "drop"
    )
  
  # Also provide summary of tested (non-zero) predictors
  tested_summary <- perm_results %>%
    filter(Tested) %>%
    mutate(
      Predictor_Type = case_when(
        grepl("^Cluster", Predictor) ~ "miRNA",
        grepl("^lncRNA", Predictor) ~ "lncRNA",
        grepl("^CpG", Predictor) ~ "CpG",
        TRUE ~ "Other"
      )
    ) %>%
    group_by(Predictor_Type) %>%
    summarize(
      N_tested = n(),
      N_significant = sum(P_adj < fdr_threshold, na.rm = TRUE),
      Pct_significant = 100 * N_significant / N_tested,
      .groups = "drop"
    )
  
  return(list(
    significant_predictors = significant,
    gene_summary = gene_summary,
    predictor_type_summary = predictor_summary,
    tested_summary = tested_summary
  ))
}


################################################################################
# VISUALIZATION FUNCTIONS
################################################################################

#' Plot Permutation Test Results for a Single Gene
plot_permutation_results <- function(gene_name, perm_results, top_n = 20, 
                                     output_dir, species_code) {
  
  gene_data <- perm_results %>%
    filter(Gene == gene_name) %>%
    arrange(P_adj) %>%
    head(top_n)
  
  if (nrow(gene_data) == 0) return(NULL)
  
  gene_data <- gene_data %>%
    mutate(
      Significant = P_adj < 0.05,
      Label = sprintf("%s\n(p=%.3f)", Predictor, P_adj)
    )
  
  p <- ggplot(gene_data, aes(x = reorder(Predictor, Abs_Observed_Coef), 
                             y = Abs_Observed_Coef)) +
    geom_point(aes(color = Significant), size = 3) +
    geom_errorbar(aes(ymin = Perm_Mean - Perm_SD, 
                      ymax = Perm_Mean + Perm_SD),
                  width = 0.2, color = "gray50", alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
    scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "red")) +
    coord_flip() +
    labs(
      title = sprintf("Permutation Test Results: %s", gene_name),
      subtitle = sprintf("Top %d predictors by coefficient magnitude", top_n),
      x = "Predictor",
      y = "Absolute Coefficient",
      color = "Significant\n(FDR < 0.05)",
      caption = "Error bars show mean ± SD of permuted coefficients"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  outfile <- file.path(output_dir, paste0(species_code, "_gene_plots/"),
                       paste0(gene_name, "_permutation_test.png"))
  ggsave(outfile, plot = p, width = 8, height = 8)
  
  return(p)
}


#' Plot Summary of Permutation Results Across All Genes
plot_permutation_summary <- function(perm_results, output_dir, species_code) {
  
  # 1. Distribution of p-values
  p1 <- ggplot(perm_results, aes(x = P_value)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "white") +
    geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
    labs(
      title = "Distribution of Permutation P-values",
      subtitle = "Red line indicates p = 0.05",
      x = "P-value",
      y = "Count"
    ) +
    theme_minimal()
  
  outfile <- file.path(output_dir, paste0(species_code, "_pvalue_distribution.png"))
  ggsave(outfile, plot = p1, width = 7, height = 5)
  
  # 2. Volcano plot
  p2 <- ggplot(perm_results, aes(x = Observed_Coef, y = -log10(P_adj))) +
    geom_point(aes(color = P_adj < 0.05), alpha = 0.3, size = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "red")) +
    labs(
      title = "Volcano Plot: Effect Size vs. Significance",
      x = "Observed Coefficient",
      y = "-log10(FDR-adjusted p-value)",
      color = "Significant\n(FDR < 0.05)"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  outfile <- file.path(output_dir, paste0(species_code, "_volcano_plot.png"))
  ggsave(outfile, plot = p2, width = 8, height = 6)
  
  # 3. Summary by predictor type
  predictor_summary <- perm_results %>%
    mutate(
      Predictor_Type = case_when(
        grepl("^Cluster", Predictor) ~ "miRNA",
        grepl("^lncRNA", Predictor) ~ "lncRNA",
        grepl("^CpG", Predictor) ~ "CpG",
        TRUE ~ "Other"
      ),
      Significant = P_adj < 0.05
    ) %>%
    group_by(Predictor_Type) %>%
    summarize(
      N_total = n(),
      N_significant = sum(Significant, na.rm = TRUE),
      Pct_significant = N_significant / N_total * 100,
      .groups = "drop"
    )
  
  p3 <- ggplot(predictor_summary, aes(x = Predictor_Type, y = Pct_significant, 
                                      fill = Predictor_Type)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)", 
                                  Pct_significant, N_significant, N_total)),
              vjust = -0.5, size = 3) +
    labs(
      title = "Percentage of Significant Associations by Predictor Type",
      x = "Predictor Type",
      y = "% Significant (FDR < 0.05)"
    ) +
    theme_minimal() +
    theme(legend.position = "none") +
    ylim(0, max(predictor_summary$Pct_significant) * 1.2)
  
  outfile <- file.path(output_dir, paste0(species_code, "_predictor_type_summary.png"))
  ggsave(outfile, plot = p3, width = 6, height = 5)
  
  return(list(p1 = p1, p2 = p2, p3 = p3, predictor_summary = predictor_summary))
}


################################################################################
# DATA LOADING AND PREPROCESSING
################################################################################

### mRNA ###
cat("Loading gene counts\n")
genes <- as.data.frame(read.csv(genes_file)) %>% select(-gene_id)
colnames(genes) <- gsub("\\.", "-", colnames(genes))
rownames(genes) <- genes[,1]
genes <- genes %>% select(-1)

### miRNA ###
cat("Loading miRNA counts\n")
miRNA <- read.table(file = miRNA_file, header = TRUE, sep = "\t", check.names = FALSE)
colnames(miRNA) <- gsub("\\.", "-", colnames(miRNA))
rownames(miRNA) <- miRNA[,1]
miRNA <- miRNA %>% select(-1)

### lncRNA ###
cat("Loading lncRNA counts\n")
lncRNA_full <- read.table(lncRNA_file, header = TRUE, sep = "\t", check.names = FALSE)
rownames(lncRNA_full) <- lncRNA_full$Geneid
lncRNA <- lncRNA_full %>% select(-Geneid, -Chr, -Start, -End, -Strand, -Length)
colnames(lncRNA) <- gsub("\\.", "-", colnames(lncRNA))

### WGBS data ###
cat("Loading WGBS data\n")
WGBS <- read.table(WGBS_file, header = TRUE)
colnames(WGBS) <- gsub("\\.", "-", colnames(WGBS))

### Load and format metadata ###
cat("Loading metadata table\n")

# Detect species prefix
if (any(grepl("ACR", colnames(miRNA)))) {
  species_prefix <- "ACR"
  species_code   <- "Apul"
} else if (any(grepl("POR", colnames(miRNA)))) {
  species_prefix <- "POR"
  species_code   <- "Peve"
} else if (any(grepl("POC", colnames(miRNA)))) {
  species_prefix <- "POC"
  species_code   <- "Ptuh"
} else {
  stop("No recognized species prefix found in column names of miRNA data")
}
cat("Species prefix and code:", species_prefix, ",", species_code, "\n")

metadata <- read_csv("../../M-multi-species/data/rna_metadata.csv", show_col_types = FALSE) %>%
  select(AzentaSampleName, ColonyID, Timepoint) %>%
  filter(grepl(species_prefix, ColonyID)) %>%
  as.data.frame()
metadata$Sample <- paste0(metadata$ColonyID, "-", metadata$Timepoint)
rownames(metadata) <- metadata$Sample
colonies <- unique(metadata$ColonyID)

metadata$Timepoint <- factor(metadata$Timepoint)
metadata$ColonyID  <- factor(metadata$ColonyID)

# Create output subdirectory for gene plots
dir.create(paste0(output_dir, "/", species_code, "_gene_plots"), showWarnings = FALSE)

# --------- Data preprocessing ----------

### Filter data sets 
genes <- genes %>% select(-any_of(excluded_samples))
miRNA <- miRNA %>% select(-any_of(excluded_samples))
lncRNA <- lncRNA %>% select(-any_of(excluded_samples))
WGBS <- WGBS %>% select(-any_of(excluded_samples))
metadata <- metadata[!(rownames(metadata) %in% excluded_samples), ]
rownames(metadata) <- metadata$Sample

#### WGBS
WGBS_filt <- WGBS %>% filter(if_all(-1, ~ .x > 0))
WGBS_filt <- as.data.frame(WGBS_filt)
rownames(WGBS_filt) <- WGBS_filt[,1]
WGBS_filt <- WGBS_filt %>% select(-1)

cat("Number of raw WGBS sites:", nrow(WGBS), "\n")
cat("Number of WGBS sites retained after filtering:", nrow(WGBS_filt), "\n")

#### RNA
genes_red <- genes[rowSums(genes) != 0, ]
miRNA_red <- miRNA[rowSums(miRNA) != 0, ]
lncRNA_red <- lncRNA[rowSums(lncRNA) != 0, ]

cat("Retained", nrow(genes_red), "of", nrow(genes), "genes;",
    nrow(miRNA_red), "of", nrow(miRNA), "miRNA; and",
    nrow(lncRNA_red), "of", nrow(lncRNA), "lncRNA through presence filtering\n")

# pOverA filtering
filtering_prop <- 3/ncol(genes_red)

# genes:
filt <- filterfun(pOverA(filtering_prop, 5))
gfilt <- genefilter(genes_red, filt)
gkeep <- genes_red[gfilt,]
gn.keep <- rownames(gkeep)
genes_filt <- as.data.frame(genes_red[which(rownames(genes_red) %in% gn.keep),])

# miRNA:
mifilt <- filterfun(pOverA(filtering_prop, 5))
mifilt <- genefilter(miRNA_red, mifilt)
mikeep <- miRNA_red[mifilt,]
mi.keep <- rownames(mikeep)
miRNA_filt <- as.data.frame(miRNA_red[which(rownames(miRNA_red) %in% mi.keep),])

# lncRNA:
lncfilt <- filterfun(pOverA(filtering_prop, 5))
lncfilt <- genefilter(lncRNA_red, lncfilt)
lnckeep <- lncRNA_red[lncfilt,]
lnc.keep <- rownames(lnckeep)
lncRNA_filt <- as.data.frame(lncRNA_red[which(rownames(lncRNA_red) %in% lnc.keep),])

cat("Retained", nrow(genes_filt), "of", nrow(genes_red), "genes;",
    nrow(miRNA_filt), "of", nrow(miRNA_red), "miRNA; and",
    nrow(lncRNA_filt), "of", nrow(lncRNA_red), "lncRNA through pOverA filtering\n")


### Transform data 
desired_order <- rownames(metadata)
genes_filt <- genes_filt[, desired_order]
miRNA_filt <- miRNA_filt[, desired_order]
lncRNA_filt <- lncRNA_filt[, desired_order]
WGBS_filt <- WGBS_filt[, desired_order]

# Variance stabilizing transformation
### genes ###
dds_genes <- suppressMessages(DESeqDataSetFromMatrix(countData = genes_filt,
                                                     colData = metadata,
                                                     design = ~Timepoint+ColonyID))
vsd_genes <- suppressMessages(varianceStabilizingTransformation(dds_genes, blind = TRUE))
vsd_genes <- assay(vsd_genes)

### miRNA ###
dds_miRNA <- suppressMessages(DESeqDataSetFromMatrix(countData = miRNA_filt,
                                                     colData = metadata,
                                                     design = ~Timepoint+ColonyID))
vsd_miRNA <- suppressMessages(varianceStabilizingTransformation(dds_miRNA, blind = TRUE))
vsd_miRNA <- assay(vsd_miRNA)

### lncRNA ###
lncRNA_filt <- lncRNA_filt %>% mutate(across(where(is.numeric), round))
dds_lncRNA <- suppressMessages(DESeqDataSetFromMatrix(countData = lncRNA_filt,
                                                      colData = metadata,
                                                      design = ~Timepoint+ColonyID))
vsd_lncRNA <- assay(vst(dds_lncRNA, blind = TRUE))

### WGBS ###
WGBS_prop <- WGBS_filt/100
WGBS_prop[WGBS_prop == 0] <- 1e-6
WGBS_prop[WGBS_prop == 1] <- 1 - 1e-6
vsd_WGBS <- log2(WGBS_prop / (1 - WGBS_prop))

### Scale standardization
vsd_miRNA_scaled <- scale(vsd_miRNA)
vsd_lncRNA_scaled <- scale(vsd_lncRNA)
vsd_WGBS_scaled <- scale(vsd_WGBS)
vsd_genes_scaled <- scale(vsd_genes)

### Merge predictor features
identical(colnames(vsd_lncRNA_scaled), colnames(vsd_miRNA_scaled))
identical(colnames(vsd_lncRNA_scaled), colnames(vsd_WGBS_scaled))

pred_counts <- rbind(vsd_lncRNA_scaled, vsd_miRNA_scaled, vsd_WGBS_scaled)
pred_counts <- t(pred_counts)
cat("Predictor set dimensions:", dim(pred_counts), "\n")

### Format
vsd_genes_t <- t(vsd_genes_scaled)
cat("Gene set dimensions:", dim(vsd_genes_t), "\n")

pred_counts <- as.data.frame(pred_counts)
vsd_genes_t <- as.data.frame(vsd_genes_t)

common_samples <- intersect(rownames(vsd_genes_t), rownames(pred_counts))
pred_counts <- pred_counts[common_samples, ]
vsd_genes_t <- vsd_genes_t[common_samples, ]


################################################################################
# PART 1: BOOTSTRAPPED ELASTIC NET (WITH PARALLEL PROCESSING)
################################################################################

cat("\n=== PART 1: Bootstrapped Elastic Net Training ===\n")

# Record start time for bootstrap
bootstrap_start_time <- Sys.time()

## First round of bootstrapping (PARALLEL)
cat("\n--- First Bootstrap Round ---\n")
bootstrap1_results <- run_bootstrap_parallel(
  response_features = vsd_genes_t,
  predictor_features = pred_counts,
  n_reps = bootstrap1_reps,
  n_cores = n_cores,
  base_seed = 703,
  alpha = alpha_value,
  progress_prefix = "Bootstrap1"
)

TM_list <- bootstrap1_results$TM_list
MP_list <- bootstrap1_results$MP_list
FI_list <- bootstrap1_results$FI_list

cat("Finished first round of bootstrapping\n")

# Combine results
all_TM <- do.call(rbind, TM_list)
all_MP <- do.call(rbind, MP_list)
all_FI <- do.call(rbind, FI_list)

# Summary statistics
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

cat("On first round of bootstrapping,", 
    sum(!is.na(MP_summary$Feature) & MP_summary$Mean_R2 > 0.5, na.rm = TRUE), 
    "genes are consistently well-predicted\n")

# Plotting
all_MP <- merge(all_MP, MP_summary[, c("Feature", "Mean_R2", "SD_R2")], by = "Feature")
all_MP$Feature <- factor(all_MP$Feature, levels = MP_summary$Feature[order(MP_summary$Mean_R2)])

p <- ggplot(all_MP, aes(x = Feature, y = Mean_R2)) +
  geom_errorbar(aes(ymin = Mean_R2 - SD_R2, ymax = Mean_R2 + SD_R2), 
                width = 0.4, color = "gray40", linewidth = 0.25) +
  geom_point(size = 2, color = "black") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "blue") +
  theme_minimal() +
  labs(title = paste0(species_code, ": Bootstrap1, Mean R^2 with Error Bars"),
       x = "Gene (Ordered by Mean R^2)",
       y = "Mean R^2 ± SD") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(-0.2, 1.2))

outfile <- file.path(output_dir, paste0(species_code, "_bootstrap1_performance_errorbars.png"))
ggsave(outfile, plot = p, width = 7, height = 7)


## Second round: bootstrapping on well-predicted genes (PARALLEL)
cat("\n--- Second Bootstrap Round ---\n")
well_predicted <- MP_summary[MP_summary$Mean_R2 > r2_threshold, ]
vsd_high_perf_t <- vsd_genes_t[, colnames(vsd_genes_t) %in% well_predicted$Feature]
cat("Number of well-predicted genes for second bootstrap:", ncol(vsd_high_perf_t), "\n")

bootstrap2_results <- run_bootstrap_parallel(
  response_features = vsd_high_perf_t,
  predictor_features = pred_counts,
  n_reps = bootstrap2_reps,
  n_cores = n_cores,
  base_seed = 1703,  # Different base seed for second round
  alpha = alpha_value,
  progress_prefix = "Bootstrap2"
)

TM_list <- bootstrap2_results$TM_list
MP_list <- bootstrap2_results$MP_list
FI_list <- bootstrap2_results$FI_list

all_TM_highperf <- do.call(rbind, TM_list)
all_MP_highperf <- do.call(rbind, MP_list)
all_FI_highperf <- do.call(rbind, FI_list)

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

cat("On second round of bootstrapping,", 
    nrow(MP_summary_highperf[MP_summary_highperf$Mean_R2 > r2_threshold, ]), 
    "genes are again consistently well-predicted\n")

# Report bootstrap timing
bootstrap_end_time <- Sys.time()
bootstrap_duration <- difftime(bootstrap_end_time, bootstrap_start_time, units = "mins")
cat(sprintf("\nBootstrapping completed in %.2f minutes\n", bootstrap_duration))


################################################################################
# PART 2: PERMUTATION SIGNIFICANCE TESTING (WITH PARALLEL PROCESSING)
################################################################################

cat("\n=== PART 2: Permutation-Based Significance Testing ===\n")
cat("This step assesses statistical significance of each predictor coefficient\n")
cat("by comparing observed coefficients to a null distribution from permuted data.\n\n")

genes_for_permutation <- colnames(vsd_high_perf_t)
cat(sprintf("Running permutation tests for %d well-predicted genes\n", 
            length(genes_for_permutation)))

# Record start time
perm_start_time <- Sys.time()

# Calculate permutation p-values (PARALLEL)
perm_results <- calculate_all_permutation_pvalues(
  response_features = vsd_high_perf_t,
  predictor_features = pred_counts,
  n_perm = n_permutations,
  alpha = alpha_value,
  n_cores = n_cores,
  base_seed = 2703,
  progress_prefix = "Permutation"
)

# Record end time and report
perm_end_time <- Sys.time()
perm_duration <- difftime(perm_end_time, perm_start_time, units = "mins")
cat(sprintf("\nPermutation testing completed in %.2f minutes\n", perm_duration))

# Apply FDR correction
cat("Applying FDR (Benjamini-Hochberg) correction for multiple testing\n")
perm_results <- apply_fdr_correction(perm_results)

# Summarize significant predictors
cat("\n=== Significant Predictor Summary ===\n")
sig_summary <- summarize_significant_predictors(perm_results, fdr_threshold = 0.05)

cat(sprintf("Number of significant predictor-gene associations (FDR < 0.05): %d\n",
            nrow(sig_summary$significant_predictors)))
cat(sprintf("Number of genes with at least one significant predictor: %d\n",
            nrow(sig_summary$gene_summary)))

cat("\nBreakdown of TESTED predictors by type:\n")
print(sig_summary$tested_summary)

cat("\nBreakdown of SIGNIFICANT predictors by type:\n")
print(sig_summary$predictor_type_summary)


################################################################################
# PART 3: VISUALIZATION AND OUTPUT
################################################################################

cat("\n=== PART 3: Generating Visualizations and Output Files ===\n")

# 1. Summary plots
cat("Creating summary plots...\n")
summary_plots <- plot_permutation_summary(perm_results, output_dir, species_code)

# 2. Per-gene permutation plots (for top genes)
cat("Creating per-gene permutation test plots...\n")
top_genes_for_plots <- sig_summary$gene_summary %>%
  arrange(desc(N_significant)) %>%
  head(20) %>%
  pull(Gene)

for (gene in top_genes_for_plots) {
  plot_permutation_results(gene, perm_results, top_n = 20, output_dir, species_code)
}

# 3. Save full permutation results
cat("Saving permutation test results...\n")

write.csv(perm_results, 
          file.path(output_dir, paste0(species_code, "_permutation_results_full.csv")),
          row.names = FALSE)

write.csv(sig_summary$significant_predictors,
          file.path(output_dir, paste0(species_code, "_significant_predictors.csv")),
          row.names = FALSE)

write.csv(sig_summary$gene_summary,
          file.path(output_dir, paste0(species_code, "_gene_summary.csv")),
          row.names = FALSE)

write.csv(sig_summary$predictor_type_summary,
          file.path(output_dir, paste0(species_code, "_predictor_type_summary.csv")),
          row.names = FALSE)

# 4. Create enhanced top predictors table with significance
cat("Enhancing top predictors table with significance information...\n")

get_predictor_type <- function(predictor_name) {
  if (startsWith(predictor_name, "Cluster")) return("miRNA")
  if (startsWith(predictor_name, "lncRNA")) return("lncRNA")
  if (startsWith(predictor_name, "CpG")) return("CpG")
  return("Other")
}

color_palette <- c(miRNA = "#E69F00", lncRNA = "#0072B2", CpG = "#009E73")

all_MP_highperf <- merge(all_MP_highperf, MP_summary_highperf[, c("Feature", "Mean_R2", "SD_R2")], by = "Feature")
all_features_highR2 <- all_MP_highperf %>%
  filter(Mean_R2 > 0.5) %>%
  pull(Feature) %>%
  unique()

get_feature_importance_for_feature <- function(model) {
  coefs <- as.matrix(coef(model, s = "lambda.min"))[-1, , drop = FALSE]
  data.frame(Predictor = rownames(coefs), Importance = abs(as.numeric(coefs)))
}

top_predictors_with_sig <- list()

for (target_feature in all_features_highR2) {
  
  gene_perm <- perm_results %>%
    filter(Gene == target_feature) %>%
    select(Predictor, P_value, P_adj)
  
  importance_all_reps <- lapply(TM_list, function(rep_entry) {
    model <- rep_entry[[target_feature]]
    if (!is.null(model)) {
      get_feature_importance_for_feature(model)
    } else {
      NULL
    }
  }) %>% bind_rows()
  
  if (nrow(importance_all_reps) == 0) next
  
  mean_importance <- importance_all_reps %>%
    group_by(Predictor) %>%
    summarise(MeanImportance = mean(Importance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(MeanImportance)) %>%
    slice_head(n = 20) %>%
    mutate(Type = sapply(Predictor, get_predictor_type))
  
  mean_importance <- mean_importance %>%
    left_join(gene_perm, by = "Predictor") %>%
    mutate(
      Significant = P_adj < 0.05,
      Significance_Label = case_when(
        P_adj < 0.001 ~ "***",
        P_adj < 0.01 ~ "**",
        P_adj < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  
  p <- ggplot(mean_importance, aes(x = reorder(Predictor, MeanImportance), 
                                   y = MeanImportance, fill = Type)) +
    geom_bar(stat = "identity", aes(alpha = Significant)) +
    geom_text(aes(label = Significance_Label), hjust = -0.5, size = 4) +
    scale_fill_manual(values = color_palette) +
    scale_alpha_manual(values = c("FALSE" = 0.4, "TRUE" = 1)) +
    coord_flip() +
    theme_minimal() +
    labs(
      title = sprintf("Top Predictors for %s", target_feature),
      subtitle = "Significance: * p<0.05, ** p<0.01, *** p<0.001 (FDR-corrected)",
      x = "Predictor",
      y = "Mean Importance (across replicates)",
      alpha = "Significant\n(FDR < 0.05)"
    )
  
  outfile <- file.path(output_dir, paste0(species_code, "_gene_plots/"),
                       paste0(target_feature, "_predictors_with_sig.png"))
  ggsave(outfile, plot = p, width = 8, height = 8)
  
  mean_importance$Feature <- target_feature
  top_predictors_with_sig[[target_feature]] <- mean_importance
}

top_predictors_df <- bind_rows(top_predictors_with_sig)
write.csv(top_predictors_df, 
          file.path(output_dir, paste0(species_code, "_top_predictors_with_significance.csv")),
          row.names = FALSE)


################################################################################
# FINAL SUMMARY
################################################################################

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("End time:", as.character(Sys.time()), "\n\n")

cat("Timing summary:\n")
cat(sprintf("  - Bootstrapping: %.2f minutes\n", as.numeric(bootstrap_duration)))
cat(sprintf("  - Permutation testing: %.2f minutes\n", as.numeric(perm_duration)))
cat(sprintf("  - Total runtime: %.2f minutes\n", 
            as.numeric(bootstrap_duration) + as.numeric(perm_duration)))

cat("\nOutput files generated:\n")
cat(sprintf("  1. %s_permutation_results_full.csv - All predictor-gene p-values\n", species_code))
cat(sprintf("  2. %s_significant_predictors.csv - FDR < 0.05 associations only\n", species_code))
cat(sprintf("  3. %s_gene_summary.csv - Summary per gene\n", species_code))
cat(sprintf("  4. %s_predictor_type_summary.csv - Summary by predictor type\n", species_code))
cat(sprintf("  5. %s_top_predictors_with_significance.csv - Top predictors with p-values\n", species_code))
cat(sprintf("  6. Various plots in %s/\n", output_dir))

cat("\nKey statistics:\n")
cat(sprintf("  - Total genes analyzed: %d\n", ncol(vsd_genes_t)))
cat(sprintf("  - Well-predicted genes (R2 > %.2f): %d\n", r2_threshold, ncol(vsd_high_perf_t)))
cat(sprintf("  - Significant predictor-gene associations: %d\n", nrow(sig_summary$significant_predictors)))
cat(sprintf("  - Bootstrap replicates: %d (round 1), %d (round 2)\n", bootstrap1_reps, bootstrap2_reps))
cat(sprintf("  - Permutations per gene: %d\n", n_permutations))
cat(sprintf("  - Cores used: %d\n", n_cores))

cat("\nStatistical notes:\n")
cat("  - P-values were calculated using permutation testing\n")
cat("  - Multiple testing correction: Benjamini-Hochberg FDR (within each gene)\n")
cat("  - Significance threshold: FDR < 0.05\n")

cat(sprintf("\nLog file saved to: %s\n", log_file_path))

cat("\n=== END OF SCRIPT ===\n")

################################################################################
# CLOSE LOGGING
################################################################################
close_logging()

# Restore original cat and print functions
cat <- original_cat
print <- original_print

cat(sprintf("Analysis complete. Log saved to: %s\n", log_file_path))