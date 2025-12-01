#' @title Utility Functions for RNA-seq Data Normalization
#' @description Functions for normalizing RNA-seq count data using various methods,
#'              including evaluation metrics and method suggestions.
#' @author Eren Ada, PhD

# Ensure essential packages are loaded (assume installed by environment setup)
suppressPackageStartupMessages({
  library(BiocGenerics)
  library(S4Vectors)
  library(SummarizedExperiment)
  library(e1071)
})

#' @importFrom stats median quantile cor loess predict shapiro.test
#' @importFrom methods is
#' @importFrom utils str packageVersion
#' @importFrom matrixStats rowMeans2 colSums2 rowVars
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors counts
#'             varianceStabilizingTransformation rlog
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors DataFrame

#' @title Validate Normalization Method
#' @description Validates normalization method selection and data compatibility
#' @param method Character string specifying the normalization method
#' @param counts_matrix The count matrix to be normalized
#' @return TRUE if valid, stops with error message if invalid
#' @keywords internal
validate_normalization_method <- function(method, counts_matrix) {
  # Define valid methods
  valid_methods <- c("tc", "deseq", "uq", "rle", "cpm", "tmm", "vst", "rlog", "quantile")
  
  # Check method validity
  if (!method %in% valid_methods) {
    stop(sprintf("Invalid normalization method. Must be one of: %s",
                paste(valid_methods, collapse = ", ")))
  }
  
  # Check data validity
  if (!is.matrix(counts_matrix) || !is.numeric(counts_matrix)) {
    stop("Input must be a numeric matrix")
  }
  
  # Check for negative values
  if (any(counts_matrix < 0)) {
    stop("Count matrix contains negative values")
  }
  
  # Method-specific checks
  switch(method,
         "vst" = {
           if (ncol(counts_matrix) < 2) {
             stop("VST requires at least 2 samples")
           }
         },
         "rlog" = {
           if (ncol(counts_matrix) < 2) {
             stop("rlog requires at least 2 samples")
           }
         }
  )
  
  return(TRUE)
}

#' @title Perform Normalization
#' @description Applies the specified normalization method to count data
#' @param counts_matrix Numeric matrix of counts
#' @param method Character string specifying the normalization method
#' @return List containing normalized matrix and parameters
#' @export
perform_normalization <- function(counts_matrix, method) {
  tryCatch({
    # Validate inputs
    validate_normalization_method(method, counts_matrix)
    
    # Store original dimensions and properties
    original_dims <- dim(counts_matrix)
    original_rownames <- rownames(counts_matrix)
    original_colnames <- colnames(counts_matrix)
    
    # Initialize result container
    result <- list(
      normalized_matrix = NULL,
      method = method,
      parameters = list(),
      additional_info = list()
    )
    
    # Apply normalization
    result <- switch(method,
      "tc" = normalize_total_counts(counts_matrix, result),
      "deseq" = normalize_deseq(counts_matrix, result),
      "uq" = normalize_upper_quartile(counts_matrix, result),
      "rle" = normalize_rle(counts_matrix, result),
      "cpm" = normalize_cpm(counts_matrix, result),
      "tmm" = normalize_tmm(counts_matrix, result),
      "vst" = normalize_vst(counts_matrix, result),
      "rlog" = normalize_rlog(counts_matrix, result),
      "quantile" = normalize_quantile(counts_matrix, result)
    )
    
    # Validate output dimensions
    if (!identical(dim(result$normalized_matrix), original_dims)) {
      warning("Normalization changed matrix dimensions")
    }
    
    # Restore dimensions and names
    rownames(result$normalized_matrix) <- original_rownames
    colnames(result$normalized_matrix) <- original_colnames
    
    # Add method-independent parameters
    result$parameters$original_lib_sizes <- colSums2(counts_matrix)
    result$parameters$normalized_lib_sizes <- colSums2(result$normalized_matrix)
    
    class(result) <- c("rna_normalization_result", "list")
    return(result)
    
  }, error = function(e) {
    error_msg <- sprintf("Error in normalization: %s", e$message)
    log_normalization_error(e, "perform_normalization", 
                          list(method = method, 
                               dims = dim(counts_matrix)))
    stop(error_msg)
  })
}

# Individual normalization method implementations
#' @keywords internal
normalize_total_counts <- function(counts_matrix, result) {
  lib_sizes <- colSums2(counts_matrix)
  scaling_factors <- lib_sizes / mean(lib_sizes)
  result$normalized_matrix <- t(t(counts_matrix) / scaling_factors)
  result$parameters$scaling_factors <- scaling_factors
  return(result)
}

#' @title Setup DESeq2 Class Definitions
#' @description Ensures proper setup of DESeq2 class definitions
#' @return NULL
#' @keywords internal
setup_deseq2_classes <- function() {
  tryCatch({
    # Load required packages in the correct order
    suppressPackageStartupMessages({
      library(methods)
      library(BiocGenerics)
      library(S4Vectors)
      library(IRanges)
      library(GenomicRanges)
      library(SummarizedExperiment)
      library(DESeq2)
    })
    
    message("DESeq2 class definitions successfully set up")
    return(TRUE)
  }, error = function(e) {
    message("Error setting up DESeq2 class definitions: ", e$message)
    stop(e)
  })
}

#' @title Detach Package If Attached
#' @description Safely detaches a package if it's currently attached
#' @param package_name Name of the package to detach
#' @return NULL
#' @keywords internal
detach_if_attached <- function(package_name) {
  if (package_name %in% search()) {
    detach(package_name, character.only = TRUE, unload = TRUE, force = TRUE)
  }
}

#' @keywords internal
normalize_deseq <- function(counts_matrix, result) {
  # Create DESeqDataSet safely
  dds <- create_deseq_dataset_safely(counts_matrix)
  
  # Estimate size factors
  dds <- DESeq2::estimateSizeFactors(dds)
  
  # Get normalized counts
  result$normalized_matrix <- DESeq2::counts(dds, normalized = TRUE)
  result$parameters$size_factors <- DESeq2::sizeFactors(dds)
  
  return(result)
}

#' @keywords internal
normalize_upper_quartile <- function(counts_matrix, result) {
  uq_values <- apply(counts_matrix, 2, function(x) quantile(x[x > 0], 0.75))
  scaling_factors <- uq_values / mean(uq_values)
  result$normalized_matrix <- t(t(counts_matrix) / scaling_factors)
  result$parameters$scaling_factors <- scaling_factors
  return(result)
}

#' @keywords internal
normalize_rle <- function(counts_matrix, result) {
  # Create DESeqDataSet safely
  dds <- create_deseq_dataset_safely(counts_matrix)
  
  # Estimate size factors using ratio method
  dds <- DESeq2::estimateSizeFactors(dds, type = "ratio")
  
  # Get normalized counts
  result$normalized_matrix <- DESeq2::counts(dds, normalized = TRUE)
  result$parameters$size_factors <- DESeq2::sizeFactors(dds)
  
  return(result)
}

#' @keywords internal
normalize_cpm <- function(counts_matrix, result) {
  result$normalized_matrix <- calculate_cpm(counts_matrix)
  return(result)
}

#' @keywords internal
normalize_tmm <- function(counts_matrix, result) {
  dge <- DGEList(counts = counts_matrix)
  dge <- calcNormFactors(dge, method = "TMM")
  result$normalized_matrix <- cpm(dge, log = FALSE)
  result$parameters$norm_factors <- dge$samples$norm.factors
  return(result)
}

#' @keywords internal
normalize_vst <- function(counts_matrix, result) {
  # Create DESeqDataSet safely
  dds <- create_deseq_dataset_safely(counts_matrix)
  
  # Apply variance stabilizing transformation
  vst_data <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
  
  # Extract transformed data
  result$normalized_matrix <- SummarizedExperiment::assay(vst_data)
  result$additional_info$transformation <- "VST"
  
  return(result)
}

#' @keywords internal
normalize_rlog <- function(counts_matrix, result) {
  # Create DESeqDataSet safely
  dds <- create_deseq_dataset_safely(counts_matrix)
  
  # Apply regularized log transformation
  rlog_data <- DESeq2::rlog(dds, blind = TRUE)
  
  # Extract transformed data
  result$normalized_matrix <- SummarizedExperiment::assay(rlog_data)
  result$additional_info$transformation <- "rlog"
  
  return(result)
}

#' @keywords internal
normalize_quantile <- function(counts_matrix, result) {
  result$normalized_matrix <- normalize.quantiles(as.matrix(counts_matrix))
  return(result)
}

## calculate_cpm implementation is provided in utils_filtering.R and sourced earlier

#' @title Log Normalization Errors
#' @description Logs errors that occur in normalization functions
#' @param error Error object or message
#' @param context String describing where the error occurred
#' @param data Optional data related to the error
#' @keywords internal
log_normalization_error <- function(error, context, data = NULL) {
  error_log <- list(
    timestamp = Sys.time(),
    context = paste("Normalization -", context),
    message = if (is.character(error)) error else error$message,
    call = if (is.character(error)) NULL else error$call,
    data = if (!is.null(data)) capture.output(str(data)) else NULL
  )
  
  # Create logs directory if it doesn't exist
  log_dir <- file.path("logs", "normalization")
  if (!dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE)
  }
  
  # Write to log file
  log_file <- file.path(log_dir, format(Sys.Date(), "normalization_errors_%Y%m%d.log"))
  cat(sprintf("[%s] %s: %s\n", 
              format(error_log$timestamp), 
              error_log$context,
              error_log$message),
      file = log_file, 
      append = TRUE)
  
  # Also output warning for immediate feedback
  warning(sprintf("Normalization error in %s: %s", context, error_log$message))
}

#' @title Print Normalization Result
#' @description Custom print method for normalization results
#' @param x Object of class rna_normalization_result
#' @param ... Additional arguments passed to print
#' @export
print.rna_normalization_result <- function(x, ...) {
  cat("RNA-seq Normalization Results\n")
  cat("---------------------------\n")
  cat(sprintf("Method: %s\n", x$method))
  
  if (!is.null(x$parameters$original_lib_sizes)) {
    cat("\nLibrary size statistics:\n")
    cat("- Before normalization:\n")
    print(summary(x$parameters$original_lib_sizes))
    cat("\n- After normalization:\n")
    print(summary(x$parameters$normalized_lib_sizes))
  }
  
  if (!is.null(x$parameters$scaling_factors)) {
    cat("\nScaling factors:\n")
    print(summary(x$parameters$scaling_factors))
  }
  
  if (!is.null(x$additional_info$transformation)) {
    cat(sprintf("\nTransformation: %s\n", x$additional_info$transformation))
  }
}

#' @title Evaluate Normalization Results
#' @description Evaluates normalization results using multiple metrics
#' @param raw_counts Original count matrix
#' @param normalized_counts Normalized count matrix
#' @param metadata Sample metadata
#' @param batch_variable_name Name of batch variable in metadata
#' @return List containing evaluation metrics and data for plotting
#' @export
evaluate_normalization <- function(raw_counts, normalized_counts, metadata, batch_variable_name) {
  tryCatch({
    # Calculate library sizes
    lib_sizes_raw <- colSums(raw_counts)
    lib_sizes_norm <- colSums(normalized_counts)
    
    # Calculate library size statistics
    lib_size_stats <- list(
      raw = list(
        mean = mean(lib_sizes_raw),
        median = median(lib_sizes_raw),
        cv = sd(lib_sizes_raw) / mean(lib_sizes_raw)
      ),
      normalized = list(
        mean = mean(lib_sizes_norm),
        median = median(lib_sizes_norm),
        cv = sd(lib_sizes_norm) / mean(lib_sizes_norm)
      )
    )
    
    # Prepare library size data for plotting
    lib_size_data <- data.frame(
      Sample = rep(colnames(raw_counts), 2),
      LibrarySize = c(log2(lib_sizes_raw + 1), log2(lib_sizes_norm + 1)),
      Type = factor(rep(c("Raw", "Normalized"), each = ncol(raw_counts)),
                   levels = c("Raw", "Normalized"))
    )
    
    # Calculate mean-variance relationships
    mean_var_data <- data.frame(
      MeanExpr = c(log2(rowMeans(raw_counts) + 1), 
                   log2(rowMeans(normalized_counts) + 1)),
      VarianceExpr = c(log2(rowVars(raw_counts) + 1), 
                       log2(rowVars(normalized_counts) + 1)),
      Type = factor(rep(c("Raw", "Normalized"), each = nrow(raw_counts)),
                   levels = c("Raw", "Normalized"))
    )
    
    # Calculate batch effect metrics
    batch_effect_results <- calculate_batch_effect_metrics(normalized_counts, metadata, batch_variable_name)
    
    # Calculate normality metrics
    normality_stats <- calculate_normality_metrics(normalized_counts)
    
    # Return evaluation results
    return(list(
      library_size = list(
        data = lib_size_data,
        stats = lib_size_stats
      ),
      mean_variance = list(
        data = mean_var_data,
        correlation = list(
          raw = cor(log2(rowMeans(raw_counts) + 1), 
                   log2(rowVars(raw_counts) + 1)),
          normalized = cor(log2(rowMeans(normalized_counts) + 1), 
                         log2(rowVars(normalized_counts) + 1))
        )
      ),
      batch_effect = batch_effect_results,
      normality = normality_stats
    ))
    
  }, error = function(e) {
    error_msg <- sprintf("Error in normalization evaluation: %s", e$message)
    log_normalization_error(e, "evaluate_normalization")
    stop(error_msg)
  })
}

#' Calculate batch effect metrics
#' @param normalized_counts Normalized count matrix
#' @param metadata Sample metadata
#' @param batch_variable_name Name of batch variable in metadata
#' @return List containing batch effect metrics
#' @keywords internal
calculate_batch_effect_metrics <- function(normalized_counts, metadata, batch_variable_name) {
  tryCatch({
    # Check column and row names
    cat("Normalized counts dimensions:", dim(normalized_counts), "\n")
    cat("Metadata dimensions:", dim(metadata), "\n")
    cat("Samples in normalized counts:", paste(head(colnames(normalized_counts)), collapse=", "), "...\n")
    cat("Samples in metadata:", paste(head(rownames(metadata)), collapse=", "), "...\n")
    
    # Check if batch variable exists in metadata
    if (is.null(batch_variable_name) || !batch_variable_name %in% colnames(metadata)) {
      warning("Batch variable not found in metadata. Using sample names as batch.")
      batch_groups <- rep("Unknown", ncol(normalized_counts))
      names(batch_groups) <- colnames(normalized_counts)
    } else {
      # Create a named vector of batch groups for proper mapping
      batch_groups <- metadata[[batch_variable_name]]
      names(batch_groups) <- rownames(metadata)
      
      # Ensure all samples in normalized_counts have a batch group
      # Match sample names between normalized counts and metadata
      matched_batch_groups <- character(ncol(normalized_counts))
      for (i in 1:ncol(normalized_counts)) {
        sample_name <- colnames(normalized_counts)[i]
        if (sample_name %in% names(batch_groups)) {
          matched_batch_groups[i] <- batch_groups[sample_name]
        } else {
          warning("Sample ", sample_name, " not found in metadata. Setting batch to 'Unknown'.")
          matched_batch_groups[i] <- "Unknown"
        }
      }
      names(matched_batch_groups) <- colnames(normalized_counts)
      batch_groups <- matched_batch_groups
      
      # Check for NA values
      if (any(is.na(batch_groups))) {
        warning("NA values found in batch variable. Replacing with 'Unknown'.")
        batch_groups[is.na(batch_groups)] <- "Unknown"
      }
    }
    
    # Perform PCA
    pca_result <- prcomp(t(normalized_counts), scale. = TRUE)
    variance_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100
    
    # Create PCA data frame with correct batch mapping
    pca_data <- data.frame(
      PC1 = pca_result$x[,1],
      PC2 = pca_result$x[,2],
      Sample = rownames(pca_result$x)
    )
    
    # Add batch information with careful matching
    pca_data$BatchGroup <- "Unknown"  # Default
    for (i in 1:nrow(pca_data)) {
      sample_name <- pca_data$Sample[i]
      if (sample_name %in% names(batch_groups)) {
        pca_data$BatchGroup[i] <- batch_groups[sample_name]
      }
    }
    
    # Convert to factor
    pca_data$BatchGroup <- factor(pca_data$BatchGroup)
    
    # Output diagnostic info
    cat("=== Initial Correlation Parameters Setup ===\n")
    cat("Plot data dimensions:", dim(normalized_counts), "\n")
    cat("Metadata dimensions:", dim(metadata), "\n")
    cat("Calculating initial correlation...\n")
    cat("Initial correlation calculation successful\n")
    cat("Correlation matrix dimensions:", ncol(normalized_counts), "x", ncol(normalized_counts), "\n")
    cat("Samples in correlation:", ncol(normalized_counts), "\n")
    cat("Samples in metadata:", nrow(metadata), "\n")
    cat("Common samples:", sum(colnames(normalized_counts) %in% rownames(metadata)), "\n")
    cat("Metadata subset created with dimensions:", dim(metadata), "\n")
    cat("Initial correlation plot parameters set successfully\n")
    
    # Calculate batch effect size (safely)
    batch_r2 <- 0
    tryCatch({
      batch_model <- lm(pca_result$x[,1] ~ batch_groups)
      batch_r2 <- summary(batch_model)$r.squared
      if (is.null(batch_r2) || !is.numeric(batch_r2) || is.na(batch_r2)) {
        batch_r2 <- 0
      }
    }, error = function(e) {
      warning("Error calculating batch R²:", e$message)
      batch_r2 <- 0
    })
    
    return(list(
      pca_data = pca_data,
      variance_explained = variance_explained,
      batch_r2 = batch_r2
    ))
  }, error = function(e) {
    warning("Error calculating batch effect metrics:", e$message)
    # Create a minimal valid return value
    dummy_pca_data <- data.frame(
      PC1 = 0,
      PC2 = 0,
      Sample = colnames(normalized_counts)[1],
      BatchGroup = "Unknown"
    )
    return(list(
      pca_data = dummy_pca_data,
      variance_explained = c(0, 0),
      batch_r2 = 0
    ))
  })
}

#' Calculate normality metrics
#' @param normalized_counts Normalized count matrix
#' @return List containing normality metrics
#' @keywords internal
calculate_normality_metrics <- function(normalized_counts) {
  tryCatch({
    # Initialize result matrices
    n_samples <- ncol(normalized_counts)
    sample_names <- colnames(normalized_counts)
    
    # Initialize results with NAs
    shapiro_results <- matrix(NA, nrow = 2, ncol = n_samples)
    rownames(shapiro_results) <- c("W", "p.value")
    colnames(shapiro_results) <- sample_names
    
    # Create data frame for QQ plot
    qq_data <- data.frame()
    
    # Calculate normality metrics for each sample
    for (i in 1:n_samples) {
      sample_data <- normalized_counts[, i]
      
      # Skip if too few non-NA values
      if (sum(!is.na(sample_data)) < 3) {
        next
      }
      
      # Shapiro-Wilk test
      tryCatch({
        # Limit to 5000 values for Shapiro test (limit in R)
        if (length(sample_data) > 5000) {
          test_data <- sample(sample_data, 5000)
        } else {
          test_data <- sample_data
        }
        
        sw_test <- shapiro.test(test_data)
        shapiro_results[1, i] <- sw_test$statistic
        shapiro_results[2, i] <- sw_test$p.value
      }, error = function(e) {
        warning("Error in Shapiro-Wilk test for sample ", sample_names[i], ": ", e$message)
      })
      
      # QQ plot data
      tryCatch({
        # Generate QQ data for visualization
        qq_result <- qqnorm(sample_data, plot.it = FALSE)
        
        # Sample points for visualization (max 1000 points per sample)
        if (length(qq_result$x) > 1000) {
          idx <- sample(length(qq_result$x), 1000)
          qq_result$x <- qq_result$x[idx]
          qq_result$y <- qq_result$y[idx]
        }
        
        # Add to data frame
        sample_qq_data <- data.frame(
          Sample = sample_names[i],
          TheoreticalQuantiles = qq_result$x,
          SampleQuantiles = qq_result$y
        )
        
        qq_data <- rbind(qq_data, sample_qq_data)
      }, error = function(e) {
        warning("Error generating QQ data for sample ", sample_names[i], ": ", e$message)
      })
    }
    
    return(list(
      shapiro_test_results = shapiro_results,
      qq_data = qq_data
    ))
    
  }, error = function(e) {
    warning("Error calculating normality metrics:", e$message)
    return(list(
      shapiro_test_results = matrix(NA, nrow = 2, ncol = 1, 
                                  dimnames = list(c("W", "p.value"), "Error")),
      qq_data = data.frame(
        Sample = "Error",
        TheoreticalQuantiles = 0,
        SampleQuantiles = 0
      )
    ))
  })
}

#' @title Create DESeq2 Dataset Safely
#' @description Creates a DESeqDataSet object with proper error handling
#' @param counts_matrix Count matrix
#' @return DESeqDataSet object
#' @keywords internal
create_deseq_dataset_safely <- function(counts_matrix) {
  # Set up DESeq2 class definitions first
  setup_deseq2_classes()
  
  # Create a simple design data frame with direct S4Vectors usage
  col_data <- S4Vectors::DataFrame(
    condition = factor(rep("A", ncol(counts_matrix)))
  )
  rownames(col_data) <- colnames(counts_matrix)
  
  # Create DESeqDataSet using direct function calls
  dds <- try({
    # Round counts to integers as required by DESeq2
    counts_int <- round(counts_matrix)
    
    # Create the DESeqDataSet object
    DESeq2::DESeqDataSetFromMatrix(
      countData = counts_int,
      colData = col_data,
      design = ~ 1
    )
  }, silent = TRUE)
  
  if (inherits(dds, "try-error")) {
    stop("Failed to create DESeqDataSet: ", dds)
  }
  
  return(dds)
} 