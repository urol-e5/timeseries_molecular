#!/usr/bin/env Rscript

################################################################################
# Elastic Net Regression with Stability Selection for Predictor Significance
################################################################################
# ============================================================================
# 26.7 = a modified copy of 26.6 with the following changes (search "[26.7]"):
#   (1) optional variance pre-filter: reduce the lncRNA & WGBS predictor layers
#       to their top-N most-variable features (arg 15, var_filter_n; 0 = off).
#   (2) BOOTSTRAPPING REPLACED BY OOF. The 80/20 train/test bootstrap (Round 1 +
#       Round 2) is removed and replaced by colony-grouped leave-one-colony-out
#       OUT-OF-FOLD (OOF) R^2  (random 80/20 splits leak a colony's timepoints 
#       across train/test). Per-gene OOF R^2 and per-sample OOF predicted expression 
#       are saved to CSV.
#   (3) PERFORMANCE GATING REMOVED. 26.6 used a first bootstrap round to drop
#       lowly-predicted genes before stability selection; 26.7 runs stability
#       selection on ALL input genes. => pass a FOCUSED genes_file (e.g. a gene
#       set)
#   (4) plotting is OPTIONAL, controlled by arg 13 `plot_enabled` (default FALSE).
#       When off, every table needed to regenerate plots post-run is written to
#       CSV. When on (pass 1), all 26.6-style PNGs are generated in-line.
#   REMOVED vs 26.6: the args r2_threshold, bootstrap1_reps, bootstrap2_reps (and the
#   26.7-interim run_oof flag), and the bootstrap helper functions (train_models,
#   train_models_split, run_parallel_bootstrap, etc.) -- all orphaned by the OOF switch.
# ============================================================================
################################################################################

cat("=== Elastic Net with Stability Selection (Parallel) ===\n")
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
suppressMessages(suppressWarnings(library(parallel)))
suppressMessages(suppressWarnings(library(doParallel)))

# Check for optional progress bar package
has_pbapply <- requireNamespace("pbapply", quietly = TRUE)
if (has_pbapply) library(pbapply)

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
################################################################################

# Store original cat and print
original_cat <- base::cat
original_print <- base::print

# Log file connection (will be opened after output_dir is known)
log_file_con <- NULL
log_file_path <- NULL

init_logging <- function(output_dir, prefix = "ElasticNet_stabsel") {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  log_file_path <<- file.path(output_dir, paste0(prefix, "_log_", timestamp, ".txt"))
  log_file_con <<- file(log_file_path, open = "wt")
  original_cat(sprintf("Logging to: %s\n", log_file_path))
  writeLines(paste0("=== Log started: ", Sys.time(), " ===\n"), log_file_con)
  flush(log_file_con)
  return(log_file_path)
}

close_logging <- function() {
  if (!is.null(log_file_con)) {
    writeLines(paste0("\n=== Log closed: ", Sys.time(), " ==="), log_file_con)
    flush(log_file_con)
    close(log_file_con)
    log_file_con <<- NULL
  }
}

tee_cat <- function(..., file = "", sep = " ", fill = FALSE, labels = NULL, append = FALSE) {
  msg <- paste0(...)
  original_cat(msg)
  if (!is.null(log_file_con) && isOpen(log_file_con)) {
    writeLines(msg, log_file_con, sep = "")
    flush(log_file_con)
  }
}

tee_print <- function(x, ...) {
  original_print(x, ...)
  if (!is.null(log_file_con) && isOpen(log_file_con)) {
    capture.output(original_print(x, ...), file = log_file_con, append = TRUE)
    flush(log_file_con)
  }
}

print_progress <- function(current, total, prefix = "Progress") {
  pct <- ifelse(total > 0, current / total * 100, 0)
  original_cat(sprintf("\r%s: %d/%d (%.1f%%)", prefix, current, total, pct))
  if (current >= total) original_cat("\n")
  
  if (!is.null(log_file_con) && isOpen(log_file_con) &&
      (current == 0 || current == total || current %% max(1, total %/% 10) == 0)) {
    writeLines(sprintf("%s: %d/%d (%.1f%%)", prefix, current, total, pct),
               log_file_con)
    flush(log_file_con)
  }
}

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
# 9  n_subsamples (optional; default 100) -- number of stability-selection subsample iterations
# 10 n_cores (optional; default 1)
# 11 pi_thr (optional; default 0.6) -- selection probability threshold
# 12 var_filter_n (optional; default 0) -- [26.7] keep top-N most-variable lncRNA & WGBS predictors (0 = off)
# 13 plot_enabled  (optional; default 0) -- [26.7] generate PNG plots in-line (1 = on, 0 = off)


min_expected <- 6
if (length(args) < min_expected) {
  stop(sprintf(paste0("Usage: Rscript script.R <genes_file> <miRNA_file> ",
                      "<lncRNA_file> <WGBS_file> <metadata_file> <output_dir> ",
                      "[excluded_samples_raw] [alpha] [n_subsamples] [n_cores] [pi_thr] [var_filter_n] [plot_enabled]\n",
                      "You provided %d args."), length(args)))
}

genes_file           <- args[1]
miRNA_file           <- args[2]
lncRNA_file          <- args[3]
WGBS_file            <- args[4]
metadata_file        <- args[5]
output_dir           <- args[6]
excluded_samples_raw <- ifelse(length(args) >= 7, args[7], NA)
alpha_value          <- ifelse(length(args) >= 8, as.numeric(args[8]), 0.5)
n_subsamples         <- ifelse(length(args) >= 9, as.integer(args[9]), 100)           # [26.7] was arg 12
n_cores              <- ifelse(length(args) >= 10, as.integer(args[10]), 1)           # [26.7] was arg 13
pi_thr               <- ifelse(length(args) >= 11, as.numeric(args[11]), 0.6)         # [26.7] was arg 14
var_filter_n         <- ifelse(length(args) >= 12, as.integer(args[12]), 0L)          # [26.7] was arg 15
if (is.na(var_filter_n) || var_filter_n < 0) var_filter_n <- 0L                       # [26.7]
plot_enabled         <- ifelse(length(args) >= 13, as.logical(as.integer(args[13])), FALSE)  # [26.7]
if (is.na(plot_enabled)) plot_enabled <- FALSE                                        # [26.7]

# Validate pi_thr
if (pi_thr <= 0.5 || pi_thr > 1.0) {
  stop("pi_thr (selection threshold) must be in (0.5, 1.0]. ",
       "Values <= 0.5 violate the Meinshausen-Bühlmann PFER bound assumptions.")
}

# Detect available cores
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
cat("\nRunning EN model with stability selection:\n")
cat(" Response genes:     ", genes_file, "\n")
cat(" Predictor miRNA:    ", miRNA_file, "\n")
cat(" Predictor lncRNA:   ", lncRNA_file, "\n")
cat(" Predictor WGBS:     ", WGBS_file, "\n")
cat(" Metadata:           ", metadata_file, "\n")
cat(" Output dir:         ", output_dir, "\n")
cat(" Excluded samples:   ", paste(excluded_samples, collapse = ", "), "\n")
cat(" Alpha:              ", alpha_value, "\n")
cat(" N subsamples:       ", n_subsamples, "\n")
cat(" N cores:            ", n_cores, "\n")
cat(" Max available cores:", max_cores, "\n")
cat(" Selection threshold:", pi_thr, "\n")
cat(" Variance filter N:  ", var_filter_n, " (0 = off)\n")   # [26.7]
cat(" Plotting enabled:  ", plot_enabled, "\n")               # [26.7]

# Create output directories
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

################################################################################
# INITIALIZE LOGGING
################################################################################
log_file_path <- init_logging(output_dir, prefix = "ElasticNet_stabsel")

# Override cat and print to use tee versions
cat <- tee_cat
print <- tee_print

cat("\nSession info:\n")
print(sessionInfo())
cat("Start time:", as.character(Sys.time()), "\n\n")

################################################################################
# SHARED PREDICTOR TYPE UTILITIES
################################################################################
# Centralized functions for classifying predictor types and assigning colors.
# Used throughout the script for consistency across all plots and summaries.
################################################################################

#' Classify predictor type from feature name (tidyverse-style, for use in mutate)
classify_type <- function(name) {
  case_when(
    grepl("^Cluster", name) ~ "miRNA",
    grepl("^lncRNA", name) ~ "lncRNA",
    grepl("^FUN", name)    ~ "Gene Body Meth",
    grepl("^CpG", name)    ~ "CpG",
    TRUE                   ~ "Other"
  )
}

#' Classify predictor type from a single feature name (base R, for use in sapply)
get_predictor_type <- function(predictor_name) {
  if (startsWith(predictor_name, "Cluster")) return("miRNA")
  if (startsWith(predictor_name, "lncRNA")) return("lncRNA")
  if (startsWith(predictor_name, "FUN")) return("Gene Body Meth")
  if (startsWith(predictor_name, "CpG")) return("CpG")
  return("Other")
}

#' Color palette for predictor types (consistent across all plots)
color_palette <- c(
  miRNA            = "#E69F00",
  lncRNA           = "#0072B2",
  CpG              = "#009E73",
  `Gene Body Meth` = "#CC79A7",
  Other            = "#999999"
)

################################################################################
# FUNCTION DEFINITIONS -- OOF R^2 & STABILITY SELECTION HELPERS
################################################################################
# [26.7] The 26.6 train/test bootstrap helpers (train_models,
# get_feature_importance, evaluate_model_performance, train_models_split,
# run_parallel_bootstrap) have been removed -- bootstrapping is replaced by
# colony-grouped OOF R^2 (calculate_oof_r2, below).
################################################################################

cat("Defining model functions\n")

#' [26.7] Colony-grouped leave-one-colony-out OUT-OF-FOLD (OOF) R^2
#'
#' For each gene: hold out one whole colony at a time (all of its timepoints),
#' fit cv.glmnet on the remaining colonies, and predict the held-out colony.
#' Repeated measures of a colony therefore never span train/test
#' The pooled out-of-fold predictions give:
#'   OOF_R2     = 1 - SS_res/SS_tot  (coefficient of determination; can be < 0,
#'                meaning worse than predicting the mean).
#'   OOF_R2_cor = squared Pearson correlation of (obs, OOF pred) -- reported for
#'                comparability with this script's bootstrap R2 (which is cor^2).
#' Also returns the per-sample OOF predicted-expression matrix (samples x genes).
#' Units: both R2 metrics are dimensionless; predictions are on the z-scaled
#' (VST) response scale.
calculate_oof_r2 <- function(response_features, predictor_features, colony,
                             alpha = alpha_value, nfolds_inner = 5,
                             n_cores = 1, base_seed = 5703) {
  X_all <- as.matrix(predictor_features)
  folds <- unique(colony)
  genes <- colnames(response_features)

  oof_one_gene <- function(g) {
    y <- response_features[[g]]
    pred <- rep(NA_real_, length(y))
    for (cf in folds) {
      te <- which(colony == cf); tr <- which(colony != cf)
      if (length(unique(y[tr])) < 3) next          # need variation to fit
      m <- tryCatch(cv.glmnet(X_all[tr, , drop = FALSE], y[tr],
                              alpha = alpha, nfolds = nfolds_inner),
                    error = function(e) NULL)
      if (is.null(m)) next
      pred[te] <- as.numeric(predict(m, X_all[te, , drop = FALSE], s = "lambda.min"))
    }
    ok <- !is.na(pred)
    if (sum(ok) < 5) return(list(R2 = NA_real_, R2cor = NA_real_, pred = pred))
    r2  <- 1 - sum((y[ok] - pred[ok])^2) / sum((y[ok] - mean(y[ok]))^2)
    rc  <- suppressWarnings(cor(y[ok], pred[ok]))
    list(R2 = r2, R2cor = if (is.na(rc)) NA_real_ else rc^2, pred = pred)
  }

  set.seed(base_seed)
  res <- mclapply(genes, oof_one_gene, mc.cores = n_cores)
  names(res) <- genes
  ok <- !vapply(res, function(z) is.null(z) || inherits(z, "try-error"), logical(1))
  res <- res[ok]

  oof_r2 <- data.frame(
    Feature    = names(res),
    OOF_R2     = vapply(res, function(z) z$R2,    numeric(1)),
    OOF_R2_cor = vapply(res, function(z) z$R2cor, numeric(1)),
    row.names  = NULL)
  oof_pred <- do.call(cbind, lapply(res, function(z) z$pred))   # samples x genes
  rownames(oof_pred) <- rownames(predictor_features)
  list(oof_r2 = oof_r2, oof_pred = oof_pred)
}


################################################################################
# STABILITY SELECTION FUNCTIONS
################################################################################
#
# These functions implement stability selection (Meinshausen & Bühlmann, 2010)
# with complementary pairs (Shah & Samworth, 2013) for each gene independently.
#
################################################################################

#' Calculate Stability Selection Probabilities for One Gene
#'
#' For a single gene, repeatedly subsample ~50% of observations, fit EN with
#' cv.glmnet (letting each subsample choose its own lambda via CV), and record
#' which predictors have non-zero coefficients. Returns selection probabilities.
#'
#' Uses complementary pairs (Shah & Samworth, 2013): each random half-split
#' produces two subsamples (the selected half and its complement), so
#' n_subsamples iterations produce 2 * n_subsamples model fits.
#'
#' @param y Response vector (gene expression)
#' @param X Predictor matrix
#' @param n_subsamples Number of subsample iterations (each produces 2 fits)
#' @param alpha Elastic net mixing parameter
#' @param n_cores Number of cores for parallel subsampling
#' @param seed Random seed for reproducibility
#' @param subsample_frac Fraction of observations per subsample (default 0.5)
#' @return Data frame with columns: Predictor, Selection_Prob, Mean_Coef,
#'         Mean_Abs_Coef, N_selected, N_total_fits
#'
calculate_stability_selection <- function(y, X, n_subsamples = 100, alpha = 0.5,
                                          n_cores = 1, seed = NULL,
                                          subsample_frac = 0.8) {
  
  if (!is.null(seed)) set.seed(seed)
  
  X <- as.matrix(X)
  n_obs <- nrow(X)
  n_predictors <- ncol(X)
  predictor_names <- colnames(X)
  subsample_size <- floor(n_obs * subsample_frac)
  
  # Adjsted complementary pare logic for subsample_frac =/ 0.5
  # only fit one model:
  # Pre-generate all subsample indices for reproducibility
  subsample_indices <- vector("list", n_subsamples)
  for (i in seq_len(n_subsamples)) {
    idx <- sample(n_obs, subsample_size)
    subsample_indices[[i]] <- list(half_a = idx)
  }
  
  # Function to fit EN on a subsample and return selection indicator + coefficients
  fit_subsample <- function(sample_idx) {
    y_sub <- y[sample_idx]
    X_sub <- X[sample_idx, , drop = FALSE]
    
    # Skip if response is constant
    if (length(unique(y_sub)) <= 1) {
      return(list(selected = rep(FALSE, n_predictors),
                  coefs = rep(0, n_predictors)))
    }
    
    model <- tryCatch({
      cv.glmnet(X_sub, y_sub, alpha = alpha, nfolds = 3)
    }, error = function(e) {
      return(NULL)
    })
    
    if (is.null(model)) {
      return(list(selected = rep(FALSE, n_predictors),
                  coefs = rep(0, n_predictors)))
    }
    
    # Extract coefficients at the model's own lambda.min
    coefs <- as.numeric(coef(model, s = "lambda.min")[-1])
    selected <- abs(coefs) > 1e-6
    
    return(list(selected = selected, coefs = coefs))
  }
  
  # Run single subsample (no complementary pair)
  run_one_pair <- function(pair_idx) {
    if (!is.null(seed)) set.seed(seed + pair_idx)
    pair <- subsample_indices[[pair_idx]]
    result_a <- fit_subsample(pair$half_a)
    return(list(result_a))
  }
  

  
  
  # Execute in parallel or sequentially
  if (n_cores > 1 && n_subsamples >= n_cores) {
    pair_results <- mclapply(seq_len(n_subsamples), run_one_pair, mc.cores = n_cores)
  } else {
    pair_results <- lapply(seq_len(n_subsamples), run_one_pair)
  }
  
  # Flatten: each pair contributes 2 fits
  all_results <- unlist(pair_results, recursive = FALSE)
  n_total_fits <- length(all_results)
  
  # Aggregate selection indicators and coefficients
  selection_matrix <- do.call(rbind, lapply(all_results, function(r) r$selected))
  coef_matrix <- do.call(rbind, lapply(all_results, function(r) r$coefs))
  
  # Calculate selection probabilities and mean coefficients
  selection_prob <- colMeans(selection_matrix, na.rm = TRUE)
  mean_coef <- colMeans(coef_matrix, na.rm = TRUE)
  mean_abs_coef <- colMeans(abs(coef_matrix), na.rm = TRUE)
  n_selected <- colSums(selection_matrix, na.rm = TRUE)
  
  results <- data.frame(
    Predictor = predictor_names,
    Selection_Prob = selection_prob,
    Mean_Coef = mean_coef,
    Mean_Abs_Coef = mean_abs_coef,
    N_Selected = as.integer(n_selected),
    N_Total_Fits = as.integer(n_total_fits),
    stringsAsFactors = FALSE
  )
  
  return(results)
}


#' Wrapper to Calculate Stability Selection for All Genes (Parallelized)
#'
#' Distributes genes across cores for parallel processing.
#' Each gene's stability selection runs independently.
#'
#' @param response_features Data frame with genes in columns
#' @param predictor_features Data frame with predictors in columns
#' @param n_subsamples Number of subsample iterations per gene
#' @param alpha Elastic net alpha parameter
#' @param n_cores Number of cores for parallel processing
#' @param base_seed Base seed for reproducibility
#' @param progress_prefix Prefix for progress display
#' @return Data frame with Gene, Predictor, Selection_Prob, etc.
#'
calculate_all_stability_selection <- function(response_features, predictor_features,
                                              n_subsamples = 100, alpha = 0.5,
                                              n_cores = 1, base_seed = 2703,
                                              progress_prefix = "StabSel") {
  
  gene_names <- colnames(response_features)
  n_genes <- length(gene_names)
  
  cat(sprintf("\nCalculating stability selection for %d genes with %d subsample iterations each\n",
              n_genes, n_subsamples))
  cat(sprintf("Total predictors per gene: %d\n", ncol(predictor_features)))
  cat(sprintf("Using %d core(s) for gene-level parallelization\n", n_cores))
  cat(sprintf("Complementary pairs: each iteration produces 2 fits (%d total fits per gene)\n",
              2 * n_subsamples))
  
  X <- as.matrix(predictor_features)
  
  process_gene <- function(gene_idx) {
    tryCatch({
      gene_name <- gene_names[gene_idx]
      y <- response_features[[gene_name]]
      gene_seed <- base_seed + gene_idx * 1000
      
      # Use n_cores=1 here -- we parallelize across genes, not within
      result <- calculate_stability_selection(
        y, X,
        n_subsamples = n_subsamples,
        alpha = alpha,
        n_cores = 1,
        seed = gene_seed
      )
      result$Gene <- gene_name
      return(result)
    }, error = function(e) {
      warning(sprintf("Stability selection failed for gene %d (%s): %s",
                      gene_idx, gene_names[gene_idx], e$message))
      return(NULL)
    })
  }
  
  # Process genes in parallel or sequentially
  if (n_cores > 1) {
    cat(sprintf("%s: Running in parallel mode across %d genes...\n",
                progress_prefix, n_genes))
    log_progress_milestone(0, n_genes, progress_prefix)
    
    # Wrap parallel call in tryCatch with sequential fallback
    results_list <- tryCatch({
      if (has_pbapply) {
        pboptions(type = "txt", style = 3, char = "=")
        pblapply(1:n_genes, process_gene, cl = n_cores)
      } else {
        mclapply(1:n_genes, process_gene, mc.cores = n_cores)
      }
    }, error = function(e) {
      cat(sprintf("\nWARNING: Parallel stability selection failed: %s\n", e$message))
      cat("Falling back to sequential execution...\n")
      NULL
    })
    
    # If parallel failed entirely, run sequentially
    if (is.null(results_list)) {
      results_list <- vector("list", n_genes)
      for (i in seq_along(gene_names)) {
        print_progress(i - 1, n_genes, prefix = paste0(progress_prefix, " (sequential fallback)"))
        results_list[[i]] <- process_gene(i)
      }
      print_progress(n_genes, n_genes, prefix = paste0(progress_prefix, " (sequential fallback)"))
    }
    log_progress_milestone(n_genes, n_genes, progress_prefix)
  } else {
    cat(sprintf("%s: Running in sequential mode...\n", progress_prefix))
    results_list <- vector("list", n_genes)
    for (i in seq_along(gene_names)) {
      print_progress(i - 1, n_genes, prefix = progress_prefix)
      results_list[[i]] <- process_gene(i)
    }
    print_progress(n_genes, n_genes, prefix = progress_prefix)
  }
  
  # Validate results: mclapply returns try-error objects on worker failure
  is_valid_result <- vapply(results_list, function(x) is.data.frame(x), logical(1))
  n_failed <- sum(!is_valid_result)
  if (n_failed > 0) {
    cat(sprintf("WARNING: %d/%d genes failed stability selection and were excluded.\n",
                n_failed, n_genes))
    failed_names <- gene_names[!is_valid_result]
    cat(sprintf("  Failed genes: %s\n",
                paste(head(failed_names, 10), collapse = ", ")))
    if (n_failed > 10) cat(sprintf("  ... and %d more\n", n_failed - 10))
  }
  results_list <- results_list[is_valid_result]
  
  if (length(results_list) == 0) {
    stop("All genes failed stability selection. Check data and model parameters.")
  }
  
  all_results <- bind_rows(results_list)
  
  # Print diagnostic summary
  n_total <- nrow(all_results)
  n_ever_selected <- sum(all_results$N_Selected > 0, na.rm = TRUE)
  avg_q <- all_results %>%
    group_by(Gene) %>%
    summarize(q = sum(Selection_Prob > 0), .groups = "drop") %>%
    pull(q) %>%
    mean()
  
  cat(sprintf("\n=== Stability Selection Summary ===\n"))
  cat(sprintf("Total predictor-gene combinations: %d\n", n_total))
  cat(sprintf("Predictors ever selected (in at least 1 fit): %d (%.2f%%)\n",
              n_ever_selected, 100 * n_ever_selected / n_total))
  cat(sprintf("Average predictors with any selection per gene: %.1f\n", avg_q))
  cat(sprintf("Subsample iterations per gene: %d (producing %d complementary-pair fits)\n",
              n_subsamples, 2 * n_subsamples))
  
  # Reorder columns
  all_results <- all_results %>%
    select(Gene, Predictor, Selection_Prob, Mean_Coef, Mean_Abs_Coef,
           N_Selected, N_Total_Fits)
  
  return(all_results)
}


#' Apply PFER Bound and Classify Stable Predictors
#'
#' Uses the Meinshausen & Bühlmann (2010) per-family error rate (PFER) bound
#' to provide error control. Predictors with selection probability above pi_thr 
#' are declared "stable."
#'
#' The PFER bound guarantees:
#'   E(V) <= q^2 / ((2 * pi_thr - 1) * p)
#' where V is the number of falsely selected variables.
#'
#' This does NOT require separate multiple testing correction (like BH-FDR),
#' because the PFER bound already controls the expected number of false
#' positives across all predictors simultaneously.
#'
#' @param stabsel_results Data frame from calculate_all_stability_selection
#' @param pi_thr Selection probability threshold (must be > 0.5)
#' @return Data frame with additional Stable column and PFER diagnostics
#'
apply_pfer_bound <- function(stabsel_results, pi_thr = 0.6) {
  
  # Calculate PFER bound per gene
  pfer_by_gene <- stabsel_results %>%
    group_by(Gene) %>%
    summarize(
      p = n(),                                         # total predictors
      q = sum(Selection_Prob > 0),                     # ever-selected predictors
      mean_q = mean(N_Selected) / mean(N_Total_Fits),  # avg selected per fit (as fraction)
      avg_selected_per_fit = mean_q * p,               # avg count selected per fit
      PFER_bound = avg_selected_per_fit^2 / ((2 * pi_thr - 1) * p),
      N_stable = sum(Selection_Prob >= pi_thr),
      .groups = "drop"
    )
  
  cat(sprintf("\nPFER bound applied with selection threshold pi = %.2f\n", pi_thr))
  cat(sprintf("Theoretical guarantee: E(false selections) <= PFER_bound per gene\n\n"))
  
  cat("Per-gene PFER diagnostics:\n")
  for (i in seq_len(nrow(pfer_by_gene))) {
    row <- pfer_by_gene[i, ]
    cat(sprintf("  %s: PFER <= %.3f | Stable predictors: %d | Avg selected/fit: %.1f | p = %d\n",
                row$Gene, row$PFER_bound, row$N_stable, row$avg_selected_per_fit, row$p))
  }
  
  # Classify predictors as stable
  stabsel_results <- stabsel_results %>%
    mutate(Stable = Selection_Prob >= pi_thr)
  
  # Attach PFER info per gene
  stabsel_results <- stabsel_results %>%
    left_join(pfer_by_gene %>% select(Gene, PFER_bound), by = "Gene")
  
  return(stabsel_results)
}


#' Summarize Stable Predictors
#'
#' @param stabsel_results Data frame with stability selection results
#' @param pi_thr Selection probability threshold
#' @return List with summary statistics and stable predictor table
#'
summarize_stable_predictors <- function(stabsel_results, pi_thr = 0.6) {
  
  stable_predictors <- stabsel_results %>%
    filter(Stable) %>%
    arrange(Gene, desc(Selection_Prob))
  
  gene_summary <- stable_predictors %>%
    group_by(Gene) %>%
    summarize(
      N_stable = n(),
      Top_predictor = Predictor[which.max(Selection_Prob)],
      Max_Selection_Prob = max(Selection_Prob),
      Mean_Selection_Prob = mean(Selection_Prob),
      .groups = "drop"
    ) %>%
    arrange(desc(N_stable))
  
  tested_summary <- stabsel_results %>%
    filter(N_Selected > 0) %>%  # ever selected in at least one fit
    mutate(Predictor_Type = classify_type(Predictor)) %>%
    group_by(Predictor_Type) %>%
    summarize(
      N_ever_selected = n(),
      N_stable = sum(Stable),
      Pct_stable = 100 * N_stable / N_ever_selected,
      .groups = "drop"
    )
  
  predictor_type_summary <- stable_predictors %>%
    mutate(Predictor_Type = classify_type(Predictor)) %>%
    group_by(Predictor_Type) %>%
    summarize(
      N_stable = n(),
      N_unique_predictors = n_distinct(Predictor),
      N_genes_affected = n_distinct(Gene),
      Mean_Selection_Prob = mean(Selection_Prob),
      .groups = "drop"
    )
  
  return(list(
    stable_predictors = stable_predictors,
    gene_summary = gene_summary,
    tested_summary = tested_summary,
    predictor_type_summary = predictor_type_summary
  ))
}


################################################################################
# STABILITY SELECTION VISUALIZATION FUNCTIONS
################################################################################

#' Plot Stability Selection Results for a Single Gene
plot_stability_results <- function(gene_name, stabsel_results, top_n = 20,
                                   output_dir, species_code, pi_thr = 0.6) {
  
  gene_data <- stabsel_results %>%
    filter(Gene == gene_name) %>%
    arrange(desc(Selection_Prob)) %>%
    head(top_n)
  
  if (nrow(gene_data) == 0) return(NULL)
  
  gene_data <- gene_data %>%
    mutate(
      Stable = Selection_Prob >= pi_thr,
      Label = sprintf("%s\n(%.0f%%)", Predictor, Selection_Prob * 100)
    )
  
  p <- ggplot(gene_data, aes(x = reorder(Predictor, Selection_Prob),
                             y = Selection_Prob)) +
    geom_point(aes(color = Stable), size = 3) +
    geom_hline(yintercept = pi_thr, linetype = "dashed", color = "red",
               linewidth = 0.8) +
    annotate("text", x = 1, y = pi_thr + 0.03,
             label = sprintf("Threshold = %.0f%%", pi_thr * 100),
             hjust = 0, color = "red", size = 3) +
    scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "steelblue")) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
    coord_flip() +
    labs(
      title = sprintf("Stability Selection: %s", gene_name),
      subtitle = sprintf("Top %d predictors by selection probability", top_n),
      x = "Predictor",
      y = "Selection Probability",
      color = sprintf("Stable\n(>= %.0f%%)", pi_thr * 100)
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  outfile <- file.path(output_dir, paste0(species_code, "_gene_plots/"),
                       paste0(gene_name, "_stability_selection.png"))
  ggsave(outfile, plot = p, width = 8, height = 8)
  
  return(p)
}


#' Plot Summary of Stability Selection Results Across All Genes
plot_stability_summary <- function(stabsel_results, output_dir, species_code,
                                   pi_thr = 0.6) {
  
  # Only include predictors that were ever selected (non-zero in at least 1 fit)
  ever_selected <- stabsel_results %>% filter(N_Selected > 0)
  
  # 1. Distribution of selection probabilities
  p1 <- ggplot(ever_selected, aes(x = Selection_Prob)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "white") +
    geom_vline(xintercept = pi_thr, linetype = "dashed", color = "red") +
    annotate("text", x = pi_thr + 0.02, y = Inf,
             label = sprintf("Threshold = %.0f%%", pi_thr * 100),
             hjust = 0, vjust = 2, color = "red", size = 3.5) +
    scale_x_continuous(labels = scales::percent_format()) +
    labs(
      title = "Distribution of Selection Probabilities",
      subtitle = "Among predictors selected at least once",
      x = "Selection Probability",
      y = "Count"
    ) +
    theme_minimal()
  
  outfile <- file.path(output_dir, paste0(species_code, "_selection_prob_distribution.png"))
  ggsave(outfile, plot = p1, width = 7, height = 5)
  
  # 2. Selection probability vs. mean coefficient magnitude
  p2 <- ggplot(ever_selected, aes(x = Mean_Abs_Coef, y = Selection_Prob)) +
    geom_point(aes(color = Stable), alpha = 0.3, size = 0.5) +
    geom_hline(yintercept = pi_thr, linetype = "dashed", color = "red") +
    scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "steelblue")) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(
      title = "Selection Probability vs. Mean Coefficient Magnitude",
      x = "Mean |Coefficient| Across Subsamples",
      y = "Selection Probability",
      color = sprintf("Stable\n(>= %.0f%%)", pi_thr * 100)
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  outfile <- file.path(output_dir, paste0(species_code, "_stabsel_scatter.png"))
  ggsave(outfile, plot = p2, width = 8, height = 6)
  
  # 3. Summary by predictor type
  predictor_summary <- stabsel_results %>%
    mutate(
      Predictor_Type = classify_type(Predictor)
    ) %>%
    group_by(Predictor_Type) %>%
    summarize(
      N_total = n(),
      N_ever_selected = sum(N_Selected > 0),
      N_stable = sum(Stable, na.rm = TRUE),
      Pct_stable = N_stable / N_total * 100,
      .groups = "drop"
    )
  
  p3 <- ggplot(predictor_summary, aes(x = Predictor_Type, y = Pct_stable,
                                      fill = Predictor_Type)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)",
                                  Pct_stable, N_stable, N_total)),
              vjust = -0.5, size = 3) +
    scale_fill_manual(values = color_palette) +
    labs(
      title = sprintf("Percentage of Stable Associations by Predictor Type (>= %.0f%%)",
                      pi_thr * 100),
      x = "Predictor Type",
      y = sprintf("%% Stable (Selection Prob >= %.0f%%)", pi_thr * 100)
    ) +
    theme_minimal() +
    theme(legend.position = "none") +
    ylim(0, max(predictor_summary$Pct_stable, na.rm = TRUE) * 1.3 + 1)
  
  outfile <- file.path(output_dir, paste0(species_code, "_predictor_type_summary.png"))
  ggsave(outfile, plot = p3, width = 6, height = 5)
  
  return(list(p1 = p1, p2 = p2, p3 = p3, predictor_summary = predictor_summary))
}


################################################################################
# ADDITIONAL VISUALIZATION FUNCTIONS
################################################################################
# These functions generate additional plots:
#   1. Stacked bar chart of predictor type composition per gene (3 orderings)
#   2. Heatmap of predictor-type summary x gene
#   [26.7] (the bootstrap R^2 plot was removed along with bootstrapping)
#
# All plots are generated in two versions:
#   - _sig: incorporating stability selection significance information
#   - _nosig: showing all non-zero predictors without significance filtering
################################################################################


#' Build Predictor Composition Data for Stacked Bar Charts
#'
#' Aggregates summed |coefficient| by predictor type for each gene.
#' Can filter to stable predictors only or include all non-zero predictors.
#'
#' @param stabsel_results Full stability selection results
#' @param MP_summary_highperf Performance summary with Mean_R2 per gene
#' @param stable_only Logical: if TRUE, only include stable predictors
#' @param pi_thr Selection probability threshold for defining "stable"
#' @return Data frame with Gene, Predictor_Type, Summed_Abs_Coef, Mean_R2
#'
build_composition_data <- function(stabsel_results, MP_summary_highperf,
                                   stable_only = FALSE, pi_thr = 0.6) {
  
  if (stable_only) {
    dat <- stabsel_results %>% filter(Stable == TRUE)
  } else {
    dat <- stabsel_results %>% filter(N_Selected > 0)
  }
  
  comp <- dat %>%
    mutate(Predictor_Type = classify_type(Predictor)) %>%
    group_by(Gene, Predictor_Type) %>%
    summarize(Summed_Abs_Coef = sum(Mean_Abs_Coef, na.rm = TRUE), .groups = "drop")
  
  # Attach R2
  comp <- comp %>%
    left_join(MP_summary_highperf %>% select(Feature, Mean_R2),
              by = c("Gene" = "Feature"))
  
  # Determine predominant predictor type per gene
  predominant <- comp %>%
    group_by(Gene) %>%
    slice_max(Summed_Abs_Coef, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(Gene, Predominant_Type = Predictor_Type)
  
  comp <- comp %>% left_join(predominant, by = "Gene")
  
  # Ensure factor ordering for predictor types
  comp$Predictor_Type <- factor(comp$Predictor_Type,
                                levels = names(color_palette))
  
  return(comp)
}


#' Plot Stacked Bar Chart of Predictor Composition
#'
#' Genes on x-axis, stacked bars show predictor type contribution.
#' Sizing adapts to gene count: bar widths and text shrink with more genes.
#'
#' @param comp Composition data frame from build_composition_data
#' @param ordering One of "alphabetical", "predominant", "r2"
#' @param output_dir Output directory
#' @param species_code Species code for filename
#' @param suffix Filename suffix (e.g., "_sig" or "_nosig")
#' @param subtitle_extra Extra text for subtitle
#'
plot_stacked_bar <- function(comp, ordering, output_dir, species_code,
                             suffix = "", subtitle_extra = "") {
  
  # Determine gene order
  gene_order <- switch(ordering,
                       "alphabetical" = {
                         sort(unique(comp$Gene))
                       },
                       "predominant" = {
                         # Order by predominant type, then alphabetically within type
                         gene_meta <- comp %>%
                           select(Gene, Predominant_Type) %>%
                           distinct() %>%
                           arrange(Predominant_Type, Gene)
                         gene_meta$Gene
                       },
                       "r2" = {
                         gene_meta <- comp %>%
                           select(Gene, Mean_R2) %>%
                           distinct() %>%
                           arrange(Mean_R2)
                         gene_meta$Gene
                       }
  )
  
  comp$Gene <- factor(comp$Gene, levels = gene_order)
  
  n_genes <- length(unique(comp$Gene))
  
  # Adaptive sizing: width scales with gene count, within bounds
  plot_width  <- max(8, min(n_genes * 0.18, 50))
  plot_height <- 6
  
  # Adaptive text sizing
  text_size <- ifelse(n_genes > 200, 3, ifelse(n_genes > 100, 4, 6))
  # Adaptive bar width: default is 0.9, shrink for very dense plots
  bar_width <- ifelse(n_genes > 200, 0.95, 0.9)
  # Adaptive R2 diamond size
  diamond_size <- ifelse(n_genes > 200, 1, ifelse(n_genes > 100, 1.5, 2.5))
  
  p <- ggplot(comp, aes(x = Gene, y = Summed_Abs_Coef, fill = Predictor_Type)) +
    geom_bar(stat = "identity", width = bar_width) +
    scale_fill_manual(values = color_palette, drop = FALSE) +
    labs(
      title = sprintf("Predictor Type Composition per Gene (%s)", species_code),
      subtitle = paste0("Ordered by: ", ordering, subtitle_extra),
      x = "Gene",
      y = "Summed |Coefficient|",
      fill = "Predictor Type"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                 size = text_size),
      legend.position = "bottom"
    )
  
  # For R2 ordering, overlay R2 on secondary axis
  if (ordering == "r2") {
    r2_data <- comp %>%
      select(Gene, Mean_R2) %>%
      distinct()
    
    # Scale R2 to the same range as Summed_Abs_Coef for dual axis
    max_coef <- max(comp$Summed_Abs_Coef, na.rm = TRUE)
    if (max_coef == 0) max_coef <- 1
    
    p <- p +
      geom_point(data = r2_data,
                 aes(x = Gene, y = Mean_R2 * max_coef, fill = NULL),
                 color = "black", shape = 18, size = diamond_size,
                 inherit.aes = FALSE) +
      scale_y_continuous(
        sec.axis = sec_axis(~ . / max_coef,
                            name = expression(R^2),
                            labels = scales::number_format(accuracy = 0.01))
      ) +
      labs(subtitle = paste0("Ordered by R² (diamonds = R²)", subtitle_extra))
  }
  
  outfile <- file.path(output_dir,
                       paste0(species_code, "_stacked_bar_", ordering, suffix, ".png"))
  ggsave(outfile, plot = p, width = plot_width, height = plot_height,
         limitsize = FALSE)
  cat(sprintf("  Saved: %s\n", outfile))
  
  return(p)
}


#' Plot Heatmap of Predictor Type Summary per Gene
#'
#' Horizontal layout: columns = genes (clustered), rows = predictor types.
#' Cell color encodes summed |coefficient|. Hierarchical clustering on columns
#' groups genes with similar regulatory profiles.
#' Width adapts to gene count for readable cell sizes.
#'
#' @param comp Composition data from build_composition_data
#' @param output_dir Output directory
#' @param species_code Species code for filename
#' @param suffix Filename suffix (e.g., "_sig" or "_nosig")
#' @param subtitle_extra Extra text for the title
#'
plot_predictor_heatmap <- function(comp, output_dir, species_code,
                                   suffix = "", subtitle_extra = "") {
  
  # Pivot to wide format: genes x predictor types
  heat_mat <- comp %>%
    select(Gene, Predictor_Type, Summed_Abs_Coef) %>%
    pivot_wider(names_from = Predictor_Type, values_from = Summed_Abs_Coef,
                values_fill = 0) %>%
    column_to_rownames("Gene") %>%
    as.matrix()
  
  # Ensure all predictor types are present as columns
  all_types <- names(color_palette)
  for (tp in all_types) {
    if (!tp %in% colnames(heat_mat)) {
      heat_mat <- cbind(heat_mat, setNames(data.frame(rep(0, nrow(heat_mat))), tp))
    }
  }
  heat_mat <- heat_mat[, intersect(all_types, colnames(heat_mat)), drop = FALSE]
  
  if (nrow(heat_mat) < 2) {
    cat("  Skipping heatmap: fewer than 2 genes\n")
    return(NULL)
  }
  
  n_genes <- nrow(heat_mat)
  
  # Transpose: predictor types as rows, genes as columns
  heat_mat_t <- t(heat_mat)
  
  # Adaptive sizing
  plot_width  <- max(8, min(n_genes * 0.18, 50))
  plot_height <- max(3, nrow(heat_mat_t) * 0.6 + 2)  # rows are few (predictor types)
  
  # Adaptive font sizes
  gene_fontsize <- ifelse(n_genes > 200, 3, ifelse(n_genes > 100, 4, 6))
  type_fontsize <- 10
  
  outfile <- file.path(output_dir,
                       paste0(species_code, "_predictor_heatmap", suffix, ".png"))
  png(outfile, width = plot_width, height = plot_height, units = "in", res = 200)
  pheatmap(heat_mat_t,
           scale = "none",
           cluster_rows = FALSE,          # predictor types: keep fixed order
           cluster_cols = TRUE,           # genes: cluster by regulatory profile
           color = colorRampPalette(c("white", "#2166AC"))(100),
           main = paste0("Predictor Type Profile per Gene (", species_code, ")\n",
                         subtitle_extra),
           fontsize_row = type_fontsize,
           fontsize_col = gene_fontsize,
           angle_col = 90,
           border_color = NA)
  dev.off()
  cat(sprintf("  Saved: %s\n", outfile))
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
rownames(WGBS) <- WGBS$gene_id
WGBS <- WGBS %>% select(-gene_id)
colnames(WGBS) <- gsub("\\.", "-", colnames(WGBS))

### Load and format metadata ###
cat("Loading metadata table\n")

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

dir.create(paste0(output_dir, "/", species_code, "_gene_plots"), showWarnings = FALSE)

### Rename columns
#colnames(genes) <- metadata$Sample[match(colnames(genes), metadata$AzentaSampleName)]
#colnames(miRNA) <- sub("_.*", "", colnames(miRNA))
#colnames(miRNA) <- metadata$Sample[match(colnames(miRNA), metadata$AzentaSampleName)]
#colnames(lncRNA) <- sub("...data.", "", colnames(lncRNA))
#colnames(lncRNA) <- sub(".sorted.bam", "", colnames(lncRNA))
#colnames(lncRNA) <- metadata$Sample[match(colnames(lncRNA), metadata$AzentaSampleName)]

### WGBS column name handling
cat("WGBS column name handling", "\n")
WGBS_colnames <- colnames(WGBS)
WGBS_colnames <- gsub("^X", "", WGBS_colnames)
WGBS_colnames <- gsub("\\.", "-", WGBS_colnames)
colnames(WGBS) <- WGBS_colnames

### Filter data sets
cat("Filtering datasets", "\n")
genes <- genes %>% select(-any_of(excluded_samples))
miRNA <- miRNA %>% select(-any_of(excluded_samples))
lncRNA <- lncRNA %>% select(-any_of(excluded_samples))
WGBS <- WGBS %>% select(-any_of(excluded_samples))
metadata <- metadata[!metadata$Sample %in% excluded_samples, ]

### Filter low-count features
cat("filtering low count features", "\n")
genes_filt <- genes[rowSums(genes > 0) >= ncol(genes) * 0.5, ]
miRNA_filt <- miRNA[rowSums(miRNA > 0) >= ncol(miRNA) * 0.5, ]
lncRNA_filt <- lncRNA[rowSums(lncRNA > 0) >= ncol(lncRNA) * 0.5, ]
WGBS_filt <- WGBS[rowSums(WGBS > 0) >= ncol(WGBS) * 0.9, ]

# Ensure integer counts
cat("Ensuring integer counts", "\n")
genes_filt <- genes_filt %>% mutate(across(where(is.numeric), round))
miRNA_filt <- miRNA_filt %>% mutate(across(where(is.numeric), round))

cat("Filtered dimensions:\n")
cat(sprintf("  Genes: %d x %d\n", nrow(genes_filt), ncol(genes_filt)))
cat(sprintf("  miRNA: %d x %d\n", nrow(miRNA_filt), ncol(miRNA_filt)))
cat(sprintf("  lncRNA: %d x %d\n", nrow(lncRNA_filt), ncol(lncRNA_filt)))
cat(sprintf("  WGBS: %d x %d\n", nrow(WGBS_filt), ncol(WGBS_filt)))

### Order columns
cat("Ensuring all dfs have identical column names and orders", "\n")

genes_filt <- genes_filt[, metadata$Sample]
miRNA_filt <- miRNA_filt[, metadata$Sample]
lncRNA_filt <- lncRNA_filt[, metadata$Sample]
WGBS_filt <- WGBS_filt[, metadata$Sample]

identical(colnames(genes_filt), colnames(miRNA_filt))
identical(colnames(genes_filt), colnames(lncRNA_filt))
identical(colnames(genes_filt), colnames(WGBS_filt))

### Variance stabilization
cat("Applying variance stabilization\n")

dds_genes <- suppressMessages(DESeqDataSetFromMatrix(countData = genes_filt,
                                                     colData = metadata,
                                                     design = ~Timepoint+ColonyID))
vsd_genes <- suppressMessages(varianceStabilizingTransformation(dds_genes, blind = TRUE))
vsd_genes <- assay(vsd_genes)

dds_miRNA <- suppressMessages(DESeqDataSetFromMatrix(countData = miRNA_filt,
                                                     colData = metadata,
                                                     design = ~Timepoint+ColonyID))
vsd_miRNA <- suppressMessages(varianceStabilizingTransformation(dds_miRNA, blind = TRUE))
vsd_miRNA <- assay(vsd_miRNA)

lncRNA_filt <- lncRNA_filt %>% mutate(across(where(is.numeric), round))
dds_lncRNA <- suppressMessages(DESeqDataSetFromMatrix(countData = lncRNA_filt,
                                                      colData = metadata,
                                                      design = ~Timepoint+ColonyID))
vsd_lncRNA <- assay(vst(dds_lncRNA, blind = TRUE))

WGBS_prop <- WGBS_filt / 100
WGBS_prop[WGBS_prop == 0] <- 1e-6
WGBS_prop[WGBS_prop == 1] <- 1 - 1e-6
vsd_WGBS <- log2(WGBS_prop / (1 - WGBS_prop))

### [26.7] Optional variance pre-filter of the high-dimensional predictor layers.
# Reduce lncRNA and WGBS to their top-N most-variable features for more tractable
# p >> n elastic net / stability selection. Variance is computed on the
# VST/logit values before z-scaling (after scaling every feature has variance 1).
# This is an unsupervised filter (it never sees the response), so it cannot leak
# information into the cross-validated / OOF performance estimates
# (cf. Hastie, Tibshirani & Friedman 2009, "Elements of Statistical Learning", 18.1).
# miRNA is always kept whole; genes (the response) are never filtered.
if (var_filter_n > 0) {
  top_var_rows <- function(mat, n) {
    if (nrow(mat) <= n) return(mat)
    v <- apply(mat, 1, var)
    mat[order(v, decreasing = TRUE)[seq_len(n)], , drop = FALSE]
  }
  cat(sprintf("[26.7] Variance pre-filter: keeping top %d most-variable lncRNA and WGBS features\n",
              var_filter_n))
  cat(sprintf("  Before: lncRNA %d, WGBS %d\n", nrow(vsd_lncRNA), nrow(vsd_WGBS)))
  vsd_lncRNA <- top_var_rows(vsd_lncRNA, var_filter_n)
  vsd_WGBS   <- top_var_rows(vsd_WGBS,   var_filter_n)
  cat(sprintf("  After:  lncRNA %d, WGBS %d, miRNA %d (kept whole)\n",
              nrow(vsd_lncRNA), nrow(vsd_WGBS), nrow(vsd_miRNA)))
}

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

vsd_genes_t <- t(vsd_genes_scaled)
cat("Gene set dimensions:", dim(vsd_genes_t), "\n")

pred_counts <- as.data.frame(pred_counts)
vsd_genes_t <- as.data.frame(vsd_genes_t)

common_samples <- intersect(rownames(vsd_genes_t), rownames(pred_counts))
pred_counts <- pred_counts[common_samples, ]
vsd_genes_t <- vsd_genes_t[common_samples, ]


################################################################################
# DATA VALIDATION AND CLEANING
################################################################################
# After variance stabilization and z-scaling, check for problems that will
# cause cv.glmnet to fail. Common issues:
#   - scale() on a zero-variance feature produces NaN
#   - Inf/-Inf from log-transforms of boundary values
#   - Constant features after VST
#   - Too few complete samples
#
# These are caught here with diagnostics
################################################################################

cat("\n=== Data Validation and Cleaning ===\n")

# --- Predictor matrix: remove NaN, Inf, and constant columns ---
pred_mat <- as.matrix(pred_counts)

# Identify problematic predictor columns
pred_has_nan    <- apply(pred_mat, 2, function(x) any(is.nan(x)))
pred_has_inf    <- apply(pred_mat, 2, function(x) any(is.infinite(x)))
pred_has_na     <- apply(pred_mat, 2, function(x) any(is.na(x)))
pred_is_const   <- apply(pred_mat, 2, function(x) sd(x, na.rm = TRUE) == 0 | is.na(sd(x, na.rm = TRUE)))

pred_problems <- pred_has_nan | pred_has_inf | pred_has_na | pred_is_const

n_nan   <- sum(pred_has_nan)
n_inf   <- sum(pred_has_inf)
n_na    <- sum(pred_has_na)
n_const <- sum(pred_is_const & !pred_has_nan & !pred_has_na)
n_total_bad <- sum(pred_problems)

cat(sprintf("Predictor matrix: %d samples x %d features\n",
            nrow(pred_mat), ncol(pred_mat)))
cat(sprintf("  Columns with NaN:           %d\n", n_nan))
cat(sprintf("  Columns with Inf:           %d\n", n_inf))
cat(sprintf("  Columns with NA:            %d\n", n_na))
cat(sprintf("  Constant columns (sd = 0):  %d\n", n_const))
cat(sprintf("  Total problematic columns:  %d\n", n_total_bad))

if (n_total_bad > 0) {
  # Report which predictor types are affected
  bad_names <- colnames(pred_mat)[pred_problems]
  bad_types <- sapply(bad_names, get_predictor_type)
  bad_type_table <- table(bad_types)
  cat("  Problematic columns by type:\n")
  for (tp in names(bad_type_table)) {
    cat(sprintf("    %s: %d\n", tp, bad_type_table[tp]))
  }
  
  # Remove problematic columns
  pred_counts <- pred_counts[, !pred_problems, drop = FALSE]
  cat(sprintf("  Removed %d problematic predictors. Remaining: %d\n",
              n_total_bad, ncol(pred_counts)))
}

# --- Response matrix: remove NaN, Inf, and constant genes ---
gene_mat <- as.matrix(vsd_genes_t)

gene_has_nan  <- apply(gene_mat, 2, function(x) any(is.nan(x)))
gene_has_inf  <- apply(gene_mat, 2, function(x) any(is.infinite(x)))
gene_has_na   <- apply(gene_mat, 2, function(x) any(is.na(x)))
gene_is_const <- apply(gene_mat, 2, function(x) sd(x, na.rm = TRUE) == 0 | is.na(sd(x, na.rm = TRUE)))

gene_problems <- gene_has_nan | gene_has_inf | gene_has_na | gene_is_const

n_gene_bad <- sum(gene_problems)

cat(sprintf("\nResponse matrix: %d samples x %d genes\n",
            nrow(gene_mat), ncol(gene_mat)))
cat(sprintf("  Genes with NaN/Inf/NA/constant values: %d\n", n_gene_bad))

if (n_gene_bad > 0) {
  removed_genes <- colnames(gene_mat)[gene_problems]
  cat(sprintf("  Removed genes: %s\n",
              paste(head(removed_genes, 20), collapse = ", ")))
  if (length(removed_genes) > 20) {
    cat(sprintf("  ... and %d more\n", length(removed_genes) - 20))
  }
  vsd_genes_t <- vsd_genes_t[, !gene_problems, drop = FALSE]
}

# --- Check for sufficient samples ---
n_samples <- nrow(pred_counts)
cat(sprintf("\nSamples available: %d\n", n_samples))
if (n_samples < 10) {
  warning("Very few samples (< 10). Model fitting may be unreliable.")
}

# --- Verify no remaining NaN/Inf ---
remaining_bad_pred <- sum(!is.finite(as.matrix(pred_counts)))
remaining_bad_gene <- sum(!is.finite(as.matrix(vsd_genes_t)))
if (remaining_bad_pred > 0 || remaining_bad_gene > 0) {
  stop(sprintf("Data still contains %d non-finite predictor values and %d non-finite gene values after cleaning. Investigate data upstream.",
               remaining_bad_pred, remaining_bad_gene))
}

cat(sprintf("\nCleaned data dimensions:\n"))
cat(sprintf("  Predictors: %d samples x %d features\n",
            nrow(pred_counts), ncol(pred_counts)))
cat(sprintf("  Genes:      %d samples x %d genes\n",
            nrow(vsd_genes_t), ncol(vsd_genes_t)))

# --- Pre-flight test: run ONE sequential model fit to catch errors early ---
cat("\nPre-flight test: fitting one gene sequentially to verify data compatibility...\n")
preflight_gene <- colnames(vsd_genes_t)[1]
preflight_result <- tryCatch({
  y_test <- vsd_genes_t[[preflight_gene]]
  X_test <- as.matrix(pred_counts)
  model_test <- cv.glmnet(X_test, y_test, alpha = alpha_value, nfolds = 3)
  preds_test <- predict(model_test, X_test, s = "lambda.min")
  r2_test <- cor(y_test, preds_test)^2
  cat(sprintf("  Pre-flight passed: gene '%s', R2 = %.4f, %d non-zero coefficients\n",
              preflight_gene, r2_test,
              sum(abs(as.numeric(coef(model_test, s = "lambda.min")[-1])) > 1e-10)))
  TRUE
}, error = function(e) {
  cat(sprintf("  Pre-flight FAILED for gene '%s': %s\n", preflight_gene, e$message))
  cat("  This error would cause all parallel workers to fail.\n")
  cat("  Investigating further...\n")
  
  # Diagnostic: check the predictor matrix directly
  X_check <- as.matrix(pred_counts)
  cat(sprintf("    Predictor matrix class: %s\n", class(X_check)))
  cat(sprintf("    Any NA: %s, Any NaN: %s, Any Inf: %s\n",
              any(is.na(X_check)), any(is.nan(X_check)), any(is.infinite(X_check))))
  cat(sprintf("    Range: [%g, %g]\n", min(X_check, na.rm = TRUE), max(X_check, na.rm = TRUE)))
  cat(sprintf("    Dimensions: %d x %d\n", nrow(X_check), ncol(X_check)))
  FALSE
})

if (!preflight_result) {
  stop("Pre-flight test failed. Fix data issues above before proceeding.")
}

cat("=== Data validation complete ===\n\n")


################################################################################
# PART 1: COLONY-GROUPED OUT-OF-FOLD (OOF) R^2   [26.7 -- REPLACES BOOTSTRAPPING]
################################################################################
# [26.7] MAJOR CHANGE vs 26.6:
#   * The 80/20 train/test bootstrapping (Round 1 + Round 2) is removed and
#     replaced by leave-one-colony-out out-of-fold (OOF) R^2. Random 80/20 splits
#     can scatter a colony's repeated timepoints across train and test (leakage) and
#     over-estimate performance; OOF holds out one whole colony at a time, giving
#     repeated-measures-safe estimate.
#   * The mid-way performance gate is removed.

cat("\n=== PART 1: Colony-grouped leave-one-colony-out OOF R^2 (all genes) ===\n")

oof_start_time <- Sys.time()

# Colony label per sample, aligned to the response/predictor rows (= common_samples).
colony_vec <- as.character(metadata$ColonyID[match(rownames(vsd_genes_t), metadata$Sample)])
n_oof_colonies <- length(unique(colony_vec[!is.na(colony_vec)]))
if (any(is.na(colony_vec))) {
  cat(sprintf("  WARNING: %d samples could not be matched to a ColonyID and are dropped from OOF.\n",
              sum(is.na(colony_vec))))
}
cat(sprintf("  %d genes, %d colonies (leave-one-colony-out, no performance gating)\n",
            ncol(vsd_genes_t), n_oof_colonies))

oof_out <- calculate_oof_r2(
  response_features  = vsd_genes_t,
  predictor_features = pred_counts,
  colony             = colony_vec,
  alpha              = alpha_value,
  nfolds_inner       = 5,
  n_cores            = n_cores,
  base_seed          = 5703
)

# Per-gene OOF R^2 table + per-sample OOF predicted-expression matrix
write.csv(oof_out$oof_r2,
          file.path(output_dir, paste0(species_code, "_oof_r2.csv")), row.names = FALSE)
write.csv(data.frame(Sample = rownames(oof_out$oof_pred), oof_out$oof_pred, check.names = FALSE),
          file.path(output_dir, paste0(species_code, "_oof_predicted_expression.csv")), row.names = FALSE)

oof_duration <- difftime(Sys.time(), oof_start_time, units = "mins")
cat(sprintf("\nOOF complete in %.2f min. median OOF_R2 = %.3f (coefficient of determination; can be < 0)\n",
            as.numeric(oof_duration), median(oof_out$oof_r2$OOF_R2, na.rm = TRUE)))

# [26.7] Per-gene R^2 summary used downstream (e.g. composition plots). This now
# comes from OOF instead of the removed Round-2 bootstrap (MP_summary_highperf).
# The column name `Mean_R2` is kept for drop-in compatibility with build_composition_data().
oof_summary <- data.frame(Feature = oof_out$oof_r2$Feature,
                          Mean_R2 = oof_out$oof_r2$OOF_R2,
                          SD_R2   = NA_real_,
                          stringsAsFactors = FALSE)

################################################################################
# PART 2: STABILITY SELECTION (REPLACES PERMUTATION TESTING)
################################################################################

cat("\n=== PART 2: Stability Selection for Predictor Reliability ===\n")
cat("This step assesses how reliably each predictor is selected across random\n")
cat("subsamples of the data, providing a measure of predictor importance.\n")
cat("Unlike permutation testing, stability selection avoids the degenerate null\n")
cat("distribution problem that arises with sparse penalized regression.\n\n")

genes_for_stabsel <- colnames(vsd_genes_t)
cat(sprintf("Running stability selection for %d genes (all input genes; no gating)\n",
            length(genes_for_stabsel)))

stabsel_start_time <- Sys.time()

# Calculate stability selection probabilities (PARALLEL)
stabsel_results <- calculate_all_stability_selection(
  response_features = vsd_genes_t,
  predictor_features = pred_counts,
  n_subsamples = n_subsamples,
  alpha = alpha_value,
  n_cores = n_cores,
  base_seed = 2703,
  progress_prefix = "StabSel"
)

stabsel_end_time <- Sys.time()
stabsel_duration <- difftime(stabsel_end_time, stabsel_start_time, units = "mins")
cat(sprintf("\nStability selection completed in %.2f minutes\n", stabsel_duration))

# Apply PFER bound and classify stable predictors
cat(sprintf("Applying PFER bound with selection threshold pi = %.2f\n", pi_thr))
stabsel_results <- apply_pfer_bound(stabsel_results, pi_thr = pi_thr)

# Summarize stable predictors
cat(sprintf("\n=== Stable Predictor Summary (Selection Prob >= %.0f%%) ===\n",
            pi_thr * 100))
stab_summary <- summarize_stable_predictors(stabsel_results, pi_thr = pi_thr)

cat(sprintf("Number of stable predictor-gene associations: %d\n",
            nrow(stab_summary$stable_predictors)))
cat(sprintf("Number of genes with at least one stable predictor: %d\n",
            nrow(stab_summary$gene_summary)))

cat("\nBreakdown of EVER-SELECTED predictors by type:\n")
print(stab_summary$tested_summary)

cat("\nBreakdown of STABLE predictors by type:\n")
print(stab_summary$predictor_type_summary)


################################################################################
# PART 3: VISUALIZATION AND OUTPUT
################################################################################

cat("\n=== PART 3: Generating Visualizations and Output Files ===\n")

if (plot_enabled) {  # [26.7] plotting optional via arg 13; data always in *_stabsel_results_full.csv / *_gene_summary.csv
# 1. Summary plots
cat("Creating summary plots...\n")
summary_plots <- plot_stability_summary(stabsel_results, output_dir, species_code,
                                        pi_thr = pi_thr)

# 2. Per-gene stability selection plots (for top genes)
cat("Creating per-gene stability selection plots...\n")
top_genes_for_plots <- stab_summary$gene_summary %>%
  arrange(desc(N_stable)) %>%
  head(20) %>%
  pull(Gene)

for (gene in top_genes_for_plots) {
  plot_stability_results(gene, stabsel_results, top_n = 20,
                         output_dir, species_code, pi_thr = pi_thr)
}
}  # [26.7] end optional plotting

# 3. Save full stability selection results
cat("Saving stability selection results...\n")

write.csv(stabsel_results,
          file.path(output_dir, paste0(species_code, "_stabsel_results_full.csv")),
          row.names = FALSE)

write.csv(stab_summary$stable_predictors,
          file.path(output_dir, paste0(species_code, "_stable_predictors.csv")),
          row.names = FALSE)

write.csv(stab_summary$gene_summary,
          file.path(output_dir, paste0(species_code, "_gene_summary.csv")),
          row.names = FALSE)

write.csv(stab_summary$predictor_type_summary,
          file.path(output_dir, paste0(species_code, "_predictor_type_summary.csv")),
          row.names = FALSE)

# 4. [26.7] Top predictors per gene -- now derived from STABILITY SELECTION
#    (Mean_Abs_Coef across subsamples + Selection_Prob), NOT from the removed
#    bootstrap models (TM_list). Each gene is annotated with its OOF R^2. The
#    per-gene plot is removed; the table is saved for post-run plotting.
cat("Building top-predictors-with-stability table from stability-selection results...\n")

oof_r2_lookup <- setNames(oof_out$oof_r2$OOF_R2, oof_out$oof_r2$Feature)

top_predictors_df <- stabsel_results %>%
  filter(N_Selected > 0) %>%                       # ever-selected predictors only
  group_by(Gene) %>%
  arrange(desc(Mean_Abs_Coef), .by_group = TRUE) %>%
  slice_head(n = 20) %>%                           # top 20 predictors per gene
  ungroup() %>%
  mutate(
    Type            = sapply(Predictor, get_predictor_type),
    OOF_R2          = unname(oof_r2_lookup[Gene]),
    Stability_Label = case_when(
      Selection_Prob >= 0.90  ~ "***",
      Selection_Prob >= 0.75  ~ "**",
      Selection_Prob >= pi_thr ~ "*",
      TRUE ~ ""
    )
  ) %>%
  rename(Feature = Gene)

write.csv(top_predictors_df,
          file.path(output_dir, paste0(species_code, "_top_predictors_with_stability.csv")),
          row.names = FALSE)


################################################################################
# PART 3.5: STACKED BAR CHARTS AND HEATMAPS
################################################################################
# Generate predictor composition plots in two versions each:
#   _sig  = only stable predictors (selection prob >= pi_thr)
#   _nosig = all non-zero predictors (no significance filtering)
################################################################################

cat("\n=== PART 3.5: Stacked Bar Charts and Heatmaps ===\n")

# --- Build composition data for both versions ---
cat("Building composition data (with significance filtering)...\n")
comp_sig <- build_composition_data(stabsel_results, oof_summary,   # [26.7] OOF R^2, was MP_summary_highperf
                                   stable_only = TRUE, pi_thr = pi_thr)

cat("Building composition data (without significance filtering)...\n")
comp_nosig <- build_composition_data(stabsel_results, oof_summary, # [26.7] OOF R^2, was MP_summary_highperf
                                     stable_only = FALSE, pi_thr = pi_thr)

# [26.7] Save composition tables (inputs the disabled stacked-bar/heatmap plots need)
write.csv(comp_sig,   file.path(output_dir, paste0(species_code, "_composition_sig.csv")),   row.names = FALSE)
write.csv(comp_nosig, file.path(output_dir, paste0(species_code, "_composition_nosig.csv")), row.names = FALSE)
cat("[26.7] Saved composition tables (_composition_sig.csv, _composition_nosig.csv)\n")

if (plot_enabled) {  # [26.7] plotting optional via arg 13; data in *_composition_{sig,nosig}.csv
# --- Stacked bar charts: WITH significance ---
cat("\nGenerating stacked bar charts (stable predictors only)...\n")
if (nrow(comp_sig) > 0) {
  sig_subtitle <- sprintf(" | Stable predictors only (>= %.0f%% selection prob)", pi_thr * 100)
  plot_stacked_bar(comp_sig, "alphabetical", output_dir, species_code,
                   suffix = "_sig", subtitle_extra = sig_subtitle)
  plot_stacked_bar(comp_sig, "predominant", output_dir, species_code,
                   suffix = "_sig", subtitle_extra = sig_subtitle)
  plot_stacked_bar(comp_sig, "r2", output_dir, species_code,
                   suffix = "_sig", subtitle_extra = sig_subtitle)
} else {
  cat("  No stable predictors found -- skipping significance-filtered stacked bars\n")
}

# --- Stacked bar charts: WITHOUT significance ---
cat("Generating stacked bar charts (all non-zero predictors)...\n")
if (nrow(comp_nosig) > 0) {
  nosig_subtitle <- " | All non-zero predictors"
  plot_stacked_bar(comp_nosig, "alphabetical", output_dir, species_code,
                   suffix = "_nosig", subtitle_extra = nosig_subtitle)
  plot_stacked_bar(comp_nosig, "predominant", output_dir, species_code,
                   suffix = "_nosig", subtitle_extra = nosig_subtitle)
  plot_stacked_bar(comp_nosig, "r2", output_dir, species_code,
                   suffix = "_nosig", subtitle_extra = nosig_subtitle)
} else {
  cat("  No non-zero predictors found -- skipping unfiltered stacked bars\n")
}

# --- Heatmaps: WITH significance ---
cat("\nGenerating predictor-type heatmaps...\n")
if (nrow(comp_sig) > 0) {
  plot_predictor_heatmap(comp_sig, output_dir, species_code,
                         suffix = "_sig",
                         subtitle_extra = sprintf("Stable predictors only (>= %.0f%%)",
                                                  pi_thr * 100))
} else {
  cat("  Skipping significance-filtered heatmap (no stable predictors)\n")
}

# --- Heatmaps: WITHOUT significance ---
if (nrow(comp_nosig) > 0) {
  plot_predictor_heatmap(comp_nosig, output_dir, species_code,
                         suffix = "_nosig",
                         subtitle_extra = "All non-zero predictors")
} else {
  cat("  Skipping unfiltered heatmap (no non-zero predictors)\n")
}
}  # [26.7] end optional plotting (Part 3.5)


################################################################################
# FINAL SUMMARY
################################################################################

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("End time:", as.character(Sys.time()), "\n\n")

cat("Timing summary:\n")                                          # [26.7]
cat(sprintf("  - OOF R^2: %.2f minutes\n", as.numeric(oof_duration)))
cat(sprintf("  - Stability selection: %.2f minutes\n", as.numeric(stabsel_duration)))
cat(sprintf("  - Total runtime: %.2f minutes\n",
            as.numeric(oof_duration) + as.numeric(stabsel_duration)))

if (plot_enabled) {  # [26.7]
cat("\nOutput files generated (CSV + PNG plots):\n")                                    # [26.7]
} else {
cat("\nOutput files generated (CSV only; plotting disabled -- pass 1 as arg 13 to enable):\n")  # [26.7]
}
cat(sprintf("  1. %s_oof_r2.csv - per-gene colony-grouped OOF R^2 (OOF_R2 = coef. of determination; OOF_R2_cor = cor^2)\n", species_code))
cat(sprintf("  2. %s_oof_predicted_expression.csv - per-sample OOF predicted expression\n", species_code))
cat(sprintf("  3. %s_stabsel_results_full.csv - all predictor-gene selection probabilities\n", species_code))
cat(sprintf("  4. %s_stable_predictors.csv - stable associations only (>= %.0f%%)\n", species_code, pi_thr * 100))
cat(sprintf("  5. %s_gene_summary.csv - summary per gene\n", species_code))
cat(sprintf("  6. %s_predictor_type_summary.csv - summary by predictor type\n", species_code))
cat(sprintf("  7. %s_top_predictors_with_stability.csv - top predictors/gene (from stability selection) + OOF R^2\n", species_code))
cat(sprintf("  8. %s_composition_sig.csv / _composition_nosig.csv - predictor-type composition (for stacked-bar/heatmap plots)\n", species_code))

cat("\nKey statistics:\n")                                        # [26.7]
cat(sprintf("  - Total genes analyzed (no performance gating): %d\n", ncol(vsd_genes_t)))
cat(sprintf("  - Median OOF R^2: %.3f\n", median(oof_out$oof_r2$OOF_R2, na.rm = TRUE)))
cat(sprintf("  - Stable predictor-gene associations: %d\n",
            nrow(stab_summary$stable_predictors)))
cat(sprintf("  - Subsample iterations: %d (producing %d complementary-pair fits)\n",
            n_subsamples, 2 * n_subsamples))
cat(sprintf("  - Variance pre-filter N: %d (0 = off)\n", var_filter_n))
cat(sprintf("  - Selection threshold: %.0f%%\n", pi_thr * 100))
cat(sprintf("  - Cores used: %d\n", n_cores))
cat(sprintf("  - Plotting enabled: %s\n", plot_enabled))                               # [26.7]

cat("\nStatistical notes:\n")
cat("  - Predictor reliability assessed via stability selection\n")
cat("    (Meinshausen & Bühlmann, 2010; Shah & Samworth, 2013)\n")
cat("  - Error control: PFER bound (expected false selections per gene)\n")
cat(sprintf("  - Selection threshold: pi = %.2f\n", pi_thr))
cat("  - No separate multiple testing correction needed (PFER bound is simultaneous)\n")

cat(sprintf("\nLog file saved to: %s\n", log_file_path))

cat("\n=== END OF SCRIPT ===\n")

################################################################################
# CLOSE LOGGING
################################################################################
close_logging()

cat <- original_cat
print <- original_print

cat(sprintf("Analysis complete. Log saved to: %s\n", log_file_path))