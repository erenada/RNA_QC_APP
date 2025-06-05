#' @title Utility Functions for Gene Filtering
#' @description Functions for filtering low-expressed genes in RNA-seq data,
#'              with support for group-aware filtering and comprehensive statistics.
#' @author Eren Ada, PhD
#' @importFrom stats median quantile
#' @importFrom methods is
#' @importFrom utils str
#' @importFrom matrixStats rowMeans2 colSums2

#' @title Validate Filtering Inputs
#' @description Validates filtering parameters before processing
#' @param min_expression Numeric threshold for expression
#' @param min_samples Numeric minimum sample count
#' @param expression_unit Character string ("raw" or "cpm")
#' @return TRUE if valid, stops with error message if invalid
#' @keywords internal
validate_filtering_inputs <- function(min_expression, min_samples, expression_unit) {
  # Check min_expression
  if (!is.numeric(min_expression) || min_expression < 0) {
    stop("Minimum expression threshold must be a non-negative number")
  }

  # Check min_samples
  if (!is.numeric(min_samples) || min_samples < 1 || !is.finite(min_samples)) {
    stop("Minimum samples must be a positive integer")
  }

  # Check expression_unit
  if (!expression_unit %in% c("raw", "cpm")) {
    stop("Expression unit must be either 'raw' or 'cpm'")
  }

  return(TRUE)
}

#' @title Filter Low Expressed Genes
#' @description Filters genes based on minimum expression in a minimum number of samples,
#'              optionally considering experimental groups defined in metadata.
#' @param counts_matrix Numeric matrix (genes x samples) of raw counts
#' @param metadata_df Data frame (samples x attributes) for sample metadata
#' @param group_column_name Character string, name of the column in metadata_df to use for grouping
#' @param min_expression Numeric, the minimum expression threshold
#' @param expression_unit Character, "raw" for raw counts or "cpm" for Counts Per Million
#' @param min_samples_in_group_or_total Numeric, the minimum number of samples that must meet the threshold
#' @return A list containing filtering results and statistics
#' @export
filter_low_expressed_genes_grouped <- function(counts_matrix, 
                                             metadata_df, 
                                             group_column_name = NULL, 
                                             min_expression, 
                                             expression_unit = "raw", 
                                             min_samples_in_group_or_total) {
  tryCatch({
    # Input Validation
    validate_filtering_inputs(min_expression, min_samples_in_group_or_total, expression_unit)
    validate_matrix_metadata_compatibility(counts_matrix, metadata_df)
    
    # Debug: Print sample names from counts and metadata
    message("Counts matrix sample names: ", paste(colnames(counts_matrix), collapse=", "))
    message("Metadata sample names: ", paste(rownames(metadata_df), collapse=", "))
    
    if (!is.null(group_column_name) && group_column_name != "None") {
      validate_group_column(group_column_name, metadata_df)
      
      # Debug: Print group information
      group_counts <- table(metadata_df[[group_column_name]])
      message("Group counts: ")
      print(group_counts)
    }

    # Convert counts_matrix to matrix if it's a data.frame
    if (is.data.frame(counts_matrix)) {
      counts_matrix <- as.matrix(counts_matrix)
    }

    # Store initial state
    genes_before <- nrow(counts_matrix)
    initial_stats <- calculate_expression_statistics(counts_matrix)

    # Prepare matrix for filtering
    matrix_to_filter <- if (expression_unit == "cpm") {
      calculate_cpm(counts_matrix)
    } else {
      counts_matrix
    }

    # Initialize tracking
    keep_gene_indices <- logical(nrow(matrix_to_filter))
    filtering_details <- list()

    if (!is.null(group_column_name) && group_column_name != "None") {
      # Group-aware filtering
      unique_groups <- unique(metadata_df[[group_column_name]])
      group_stats <- list()
      
      # Debug: Print group processing information
      message("Processing groups: ", paste(unique_groups, collapse=", "))

      for (group_level in unique_groups) {
        samples_in_current_group <- rownames(metadata_df[metadata_df[[group_column_name]] == group_level, ])
        message(sprintf("Group %s has %d samples: %s", 
                      group_level, 
                      length(samples_in_current_group),
                      paste(samples_in_current_group, collapse=", ")))
        
        counts_for_group <- matrix_to_filter[, samples_in_current_group, drop = FALSE]
        
        # Check group size
        if (ncol(counts_for_group) < min_samples_in_group_or_total) {
          warning(sprintf("Group %s has fewer samples (%d) than minimum required (%d)",
                        group_level, ncol(counts_for_group), min_samples_in_group_or_total))
          group_stats[[group_level]] <- list(
            samples = length(samples_in_current_group),
            warning = "Insufficient samples"
          )
          next
        }

        # Calculate passing genes for this group
        samples_meeting_threshold <- rowSums(counts_for_group >= min_expression)
        genes_passing_in_group <- samples_meeting_threshold >= min_samples_in_group_or_total
        
        # Update tracking
        keep_gene_indices <- keep_gene_indices | genes_passing_in_group
        group_stats[[group_level]] <- list(
          samples = length(samples_in_current_group),
          genes_passing = sum(genes_passing_in_group),
          expression_distribution = summary(rowMeans2(counts_for_group))
        )
      }
      
      filtering_details$group_stats <- group_stats
    } else {
      # Global filtering
      samples_meeting_threshold <- rowSums(matrix_to_filter >= min_expression)
      keep_gene_indices <- samples_meeting_threshold >= min_samples_in_group_or_total
      
      filtering_details$global_stats <- list(
        samples_total = ncol(matrix_to_filter),
        expression_distribution = summary(rowMeans2(matrix_to_filter))
      )
    }

    # Apply filtering to original counts
    filtered_matrix <- counts_matrix[keep_gene_indices, , drop = FALSE]
    genes_after <- nrow(filtered_matrix)
    genes_removed <- genes_before - genes_after

    # Calculate post-filtering statistics
    final_stats <- calculate_expression_statistics(filtered_matrix)

    # Prepare detailed summary
    filtering_summary <- list(
      filtered_matrix = filtered_matrix,
      genes_before = genes_before,
      genes_after = genes_after,
      genes_removed = genes_removed,
      filtering_criteria = list(
        min_expression = min_expression,
        expression_unit = expression_unit,
        min_samples = min_samples_in_group_or_total,
        group_aware = !is.null(group_column_name) && group_column_name != "None"
      ),
      statistics = list(
        initial = initial_stats,
        final = final_stats
      ),
      details = filtering_details
    )

    class(filtering_summary) <- c("rna_filtering_result", "list")
    return(filtering_summary)

  }, error = function(e) {
    error_msg <- sprintf("Error in gene filtering: %s", e$message)
    log_filtering_error(e, "filter_low_expressed_genes", 
                       list(min_expr = min_expression, 
                            min_samples = min_samples_in_group_or_total,
                            expression_unit = expression_unit))
    stop(error_msg)
  })
}

#' @title Calculate Expression Statistics
#' @description Calculates basic statistics for a count matrix
#' @param count_matrix Numeric matrix of counts
#' @return List of statistics
#' @keywords internal
calculate_expression_statistics <- function(count_matrix) {
  # Convert to matrix if it's a data.frame
  if (is.data.frame(count_matrix)) {
    count_matrix <- as.matrix(count_matrix)
  }
  
  mean_counts <- rowMeans2(count_matrix)
  return(list(
    total_counts = colSums2(count_matrix),
    mean_per_gene = mean(mean_counts),
    median_per_gene = median(mean_counts),
    zero_proportion = mean(count_matrix == 0),
    quantiles = quantile(mean_counts, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
  ))
}

#' @title Validate Matrix Metadata Compatibility
#' @description Validates that count matrix and metadata match
#' @param counts_matrix Count matrix to validate
#' @param metadata_df Metadata data frame to validate
#' @return TRUE if valid, stops with error message if invalid
#' @keywords internal
validate_matrix_metadata_compatibility <- function(counts_matrix, metadata_df) {
  if (!all(colnames(counts_matrix) %in% rownames(metadata_df))) {
    missing_samples <- setdiff(colnames(counts_matrix), rownames(metadata_df))
    stop(sprintf("Not all samples in count matrix are present in metadata. Missing: %s",
                paste(missing_samples, collapse = ", ")))
  }
  if (!all(rownames(metadata_df) %in% colnames(counts_matrix))) {
    extra_samples <- setdiff(rownames(metadata_df), colnames(counts_matrix))
    warning(sprintf("Some metadata samples are not present in count matrix: %s",
                   paste(extra_samples, collapse = ", ")))
  }
  return(TRUE)
}

#' @title Validate Group Column
#' @description Validates group column in metadata
#' @param group_column_name Name of the group column
#' @param metadata_df Metadata data frame
#' @return TRUE if valid, stops with error message if invalid
#' @keywords internal
validate_group_column <- function(group_column_name, metadata_df) {
  if (!group_column_name %in% colnames(metadata_df)) {
    stop(sprintf("Group column '%s' not found in metadata", group_column_name))
  }
  unique_groups <- unique(metadata_df[[group_column_name]])
  if (length(unique_groups) < 2) {
    stop(sprintf("Group column '%s' must have at least 2 unique values, found: %s",
                group_column_name, paste(unique_groups, collapse = ", ")))
  }
  return(TRUE)
}

#' @title Calculate CPM
#' @description Calculates Counts Per Million
#' @param raw_counts_matrix Raw counts matrix
#' @return CPM matrix
#' @keywords internal
calculate_cpm <- function(raw_counts_matrix) {
  # Convert to matrix if it's a data.frame
  if (is.data.frame(raw_counts_matrix)) {
    raw_counts_matrix <- as.matrix(raw_counts_matrix)
  }
  
  library_sizes <- colSums2(raw_counts_matrix)
  if (any(library_sizes == 0)) {
    warning("Some samples have zero total counts. CPM calculation may produce NaN/Inf values.")
  }
  cpm_matrix <- t(t(raw_counts_matrix) / library_sizes * 1e6)
  return(cpm_matrix)
}

#' @title Log Filtering Errors
#' @description Logs errors that occur in filtering functions
#' @param error Error object or message
#' @param context String describing where the error occurred
#' @param data Optional data related to the error
#' @keywords internal
log_filtering_error <- function(error, context, data = NULL) {
  error_log <- list(
    timestamp = Sys.time(),
    context = paste("Filtering -", context),
    message = if (is.character(error)) error else error$message,
    call = if (is.character(error)) NULL else error$call,
    data = if (!is.null(data)) capture.output(str(data)) else NULL
  )
  
  # Create logs directory if it doesn't exist
  log_dir <- file.path("logs", "filtering")
  if (!dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE)
  }
  
  # Write to log file
  log_file <- file.path(log_dir, format(Sys.Date(), "filtering_errors_%Y%m%d.log"))
  cat(sprintf("[%s] %s: %s\n", 
              format(error_log$timestamp), 
              error_log$context,
              error_log$message),
      file = log_file, 
      append = TRUE)
  
  # Also output warning for immediate feedback
  warning(sprintf("Filtering error in %s: %s", context, error_log$message))
}

#' @title Print Filtering Result
#' @description Custom print method for filtering results
#' @param x Object of class rna_filtering_result
#' @param ... Additional arguments passed to print
#' @export
print.rna_filtering_result <- function(x, ...) {
  # Helper function to format numbers without scientific notation
  format_number <- function(num) {
    if (num >= 1000) {
      sprintf("%.2f", num)
    } else if (num >= 10) {
      sprintf("%.3f", num)
    } else if (num >= 1) {
      sprintf("%.6f", num)
    } else {
      sprintf("%.6f", num)
    }
  }
  
  cat("RNA-seq Filtering Results\n")
  cat("------------------------\n")
  cat(sprintf("Initial genes: %d\n", x$genes_before))
  cat(sprintf("Genes after filtering: %d\n", x$genes_after))
  cat(sprintf("Genes removed: %d (%.1f%%)\n", 
              x$genes_removed, 
              x$genes_removed/x$genes_before * 100))
  cat("\nFiltering criteria:\n")
  cat(sprintf("- Minimum expression: %s (%s)\n", 
              x$filtering_criteria$min_expression,
              x$filtering_criteria$expression_unit))
  cat(sprintf("- Minimum samples: %d\n", x$filtering_criteria$min_samples))
  cat(sprintf("- Group-aware: %s\n", x$filtering_criteria$group_aware))
  
  cat("\nExpression statistics:\n")
  cat("- Before filtering:\n")
  
  # Format quantiles without scientific notation
  before_quantiles <- x$statistics$initial$quantiles
  formatted_before <- sapply(before_quantiles, format_number)
  names(formatted_before) <- names(before_quantiles)
  print(formatted_before, quote = FALSE)
  
  cat("\n- After filtering:\n")
  
  # Format quantiles without scientific notation
  after_quantiles <- x$statistics$final$quantiles
  formatted_after <- sapply(after_quantiles, format_number)
  names(formatted_after) <- names(after_quantiles)
  print(formatted_after, quote = FALSE)
} 