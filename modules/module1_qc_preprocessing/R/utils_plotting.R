#' @title Utility Functions for Plotting in Filtering & Normalization Tab
#' @description A collection of plotting functions for RNA-seq data visualization,
#'              specifically for the filtering and normalization steps.
#' @author Eren Ada, PhD
#' @importFrom ggplot2 theme_bw theme element_text element_blank labs aes geom_bar geom_text
#'             scale_fill_brewer scale_color_brewer geom_density geom_point geom_smooth
#'             stat_ellipse annotate theme_void element_rect margin scale_x_continuous
#'             facet_wrap geom_abline
#' @importFrom plotly ggplotly
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats quantile qqnorm
#' @importFrom utils packageVersion
#' @importFrom grDevices colorRampPalette
#' @importFrom dplyr %>%
#' @importFrom scales comma_format
#' @importFrom matrixStats rowVars
#' @importFrom tidyr gather

#' @title Get RNA Processing Theme
#' @description Returns a consistent theme for all plots in the RNA processing module
#' @return A ggplot2 theme object
#' @export
get_theme_rna_processing <- function() {
  theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )
}

#' @title Get Color Palette
#' @description Returns a consistent color palette for plots
#' @param n Number of colors needed
#' @return A vector of color codes
#' @export
get_color_palette <- function(n) {
  if (n <= 8) {
    return(brewer.pal(max(3, n), "Set2")[1:n])  # Ensures minimum palette size
  } else {
    return(colorRampPalette(brewer.pal(8, "Set2"))(n))
  }
}

#' @title Plot Filtering Impact
#' @description Creates a visual summary of filtering impact on gene counts
#' @param filtering_summary List containing filtering results and statistics
#' @return A ggplot object
#' @export
plot_filtering_impact <- function(filtering_summary) {
  tryCatch({
    # Prepare data
    impact_data <- data.frame(
      stage = factor(c("Pre-filtering", "Post-filtering"),
                    levels = c("Pre-filtering", "Post-filtering")),
      genes = c(filtering_summary$genes_before, filtering_summary$genes_after),
      zero_prop = c(
        filtering_summary$statistics$initial$zero_proportion,
        filtering_summary$statistics$final$zero_proportion
      )
    )
    
    # Create plot
    p <- ggplot() +
      geom_bar(data = impact_data,
               aes(x = stage, y = genes, fill = stage),
               stat = "identity", width = 0.7) +
      geom_text(data = impact_data,
                aes(x = stage, y = genes,
                    label = sprintf("%s\n(%s%% zeros)",
                                  format(genes, big.mark = ",", scientific = FALSE),
                                  round(zero_prop * 100, 1))),
                vjust = -0.5) +
      scale_fill_brewer(palette = "Set2") +
      scale_y_continuous(labels = comma_format()) +
      labs(title = "Impact of Filtering",
           y = "Number of Genes",
           x = "") +
      get_theme_rna_processing() +
      theme(legend.position = "none")
    
    return(p)
    
  }, error = function(e) {
    log_plotting_error(e, "filtering_impact")
    return(create_error_plot("Could not create filtering impact plot"))
  })
}

#' @title Plot Library Size Comparison
#' @description Creates library size comparison plot between conditions
#' @param lib_size_data data.frame with library sizes pre/post normalization
#' @return A plotly object
#' @export
plot_lib_size_comparison <- function(lib_size_data) {
  tryCatch({
    p <- ggplot(lib_size_data, 
                aes(x = LibrarySize, fill = Type, color = Type)) +
      geom_density(alpha = 0.5) +
      scale_fill_brewer(palette = "Set2") +
      scale_color_brewer(palette = "Set2") +
      scale_x_continuous(labels = comma_format()) +
      labs(title = "Library Size Distribution",
           x = "Log2(Library Size + 1)",
           y = "Density") +
      get_theme_rna_processing()
    
    return(ggplotly(p, tooltip = c("x", "y", "Type")))
    
  }, error = function(e) {
    log_plotting_error(e, "library_size_comparison")
    return(create_error_plot("Could not create library size comparison plot"))
  })
}

#' @title Plot Mean-Variance Relationship
#' @description Creates mean-variance relationship plot
#' @param mean_variance_data data.frame with mean and variance values
#' @return A plotly object
#' @export
plot_mean_variance <- function(mean_variance_data) {
  tryCatch({
    p <- ggplot(mean_variance_data, 
                aes(x = MeanExpr, y = VarianceExpr, color = Type)) +
      geom_point(alpha = 0.3, size = 1) +
      geom_smooth(method = "loess", se = FALSE) +
      scale_color_brewer(palette = "Set2") +
      scale_x_continuous(labels = comma_format()) +
      scale_y_continuous(labels = comma_format()) +
      labs(title = "Mean-Variance Relationship",
           x = "Log2(Mean Expression)",
           y = "Log2(Variance)") +
      get_theme_rna_processing()
    
    return(ggplotly(p, tooltip = c("x", "y", "Type")))
    
  }, error = function(e) {
    log_plotting_error(e, "mean_variance")
    return(create_error_plot("Could not create mean-variance plot"))
  })
}

#' @title Plot PCA by Batch
#' @description Creates PCA plot colored by batch
#' @param pca_data data.frame with PC coordinates and batch information
#' @param variance_explained vector of variance explained by PCs
#' @return A plotly object
#' @export
plot_pca_batch <- function(pca_data, variance_explained) {
  tryCatch({
    p <- ggplot(pca_data, 
                aes(x = PC1, y = PC2, 
                    color = BatchGroup, 
                    text = Sample)) +
      geom_point(size = 3) +
      stat_ellipse(aes(group = BatchGroup), 
                   type = "norm", 
                   level = 0.95, 
                   show.legend = FALSE) +
      scale_color_brewer(palette = "Set2") +
      labs(title = "PCA of Normalized Data",
           x = sprintf("PC1 (%s%%)", round(variance_explained[1], 1)),
           y = sprintf("PC2 (%s%%)", round(variance_explained[2], 1)),
           color = "Batch Group") +
      get_theme_rna_processing()
    
    return(ggplotly(p, tooltip = c("text", "x", "y", "color")))
    
  }, error = function(e) {
    log_plotting_error(e, "pca_batch")
    return(create_error_plot("Could not create PCA batch plot"))
  })
}

#' @title Plot Normality Q-Q
#' @description Creates Q-Q plots for normalized data
#' @param normalized_data Normalized count matrix
#' @param sample_subset Optional vector of sample names to plot
#' @return A plotly object
#' @export
plot_normality_qq <- function(normalized_data, sample_subset = NULL) {
  tryCatch({
    # Select samples to plot
    if (is.null(sample_subset)) {
      sample_subset <- colnames(normalized_data)[1:min(ncol(normalized_data), 9)]
    }
    
    # Prepare data for plotting
    qq_data <- data.frame()
    for (sample in sample_subset) {
      qq_result <- qqnorm(normalized_data[, sample], plot.it = FALSE)
      qq_data <- rbind(qq_data, data.frame(
        Sample = sample,
        TheoreticalQuantiles = qq_result$x,
        SampleQuantiles = qq_result$y
      ))
    }
    
    # Create plot
    p <- ggplot(qq_data, 
                aes(x = TheoreticalQuantiles, 
                    y = SampleQuantiles, 
                    color = Sample)) +
      geom_point(alpha = 0.6) +
      geom_abline(intercept = 0, slope = 1, 
                  linetype = "dashed", 
                  color = "gray50") +
      facet_wrap(~Sample) +
      scale_color_brewer(palette = "Set2") +
      labs(title = "Normal Q-Q Plots",
           x = "Theoretical Quantiles",
           y = "Sample Quantiles") +
      get_theme_rna_processing() +
      theme(legend.position = "none")
    
    return(ggplotly(p, tooltip = c("x", "y", "Sample")))
    
  }, error = function(e) {
    log_plotting_error(e, "normality_qq", list(samples = sample_subset))
    return(create_error_plot("Could not create normality Q-Q plot"))
  })
}

#' @title Create Error Plot
#' @description Creates a simple error message plot
#' @param error_message Error message to display
#' @return A ggplot object
#' @export
create_error_plot <- function(error_message = "Error creating plot") {
  # Create a simple error plot
  p <- ggplot() + 
    annotate("text", x = 0.5, y = 0.5, 
             label = error_message,
             size = 5, color = "red") +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "gray95"),
      plot.margin = margin(20, 20, 20, 20)
    )
  
  return(p)
}

#' @title Log Plotting Errors
#' @description Logs errors that occur in plotting functions
#' @param error Error object or message
#' @param context String describing where the error occurred
#' @param data Optional data related to the error
#' @keywords internal
log_plotting_error <- function(error, context, data = NULL) {
  error_log <- list(
    timestamp = Sys.time(),
    context = paste("Plotting -", context),
    message = if (is.character(error)) error else error$message,
    call = if (is.character(error)) NULL else error$call,
    data = if (!is.null(data)) capture.output(str(data)) else NULL
  )
  
  # Create logs directory if it doesn't exist
  log_dir <- file.path("logs", "plotting")
  if (!dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE)
  }
  
  # Write to log file
  log_file <- file.path(log_dir, format(Sys.Date(), "plotting_errors_%Y%m%d.log"))
  cat(sprintf("[%s] %s: %s\n", 
              format(error_log$timestamp), 
              error_log$context,
              error_log$message),
      file = log_file, 
      append = TRUE)
  
  # Also output warning for immediate feedback
  warning(sprintf("Plotting error in %s: %s", context, error_log$message))
}

#' @title Plot Normalization Library Size Evaluation
#' @description Creates boxplots comparing library sizes before and after normalization
#' @param eval_results List containing evaluation results
#' @return A plotly object
#' @export
plot_normalization_boxplots <- function(eval_results) {
  tryCatch({
    # Extract library size data
    lib_size_data <- eval_results$library_size$data
    
    p <- ggplot(lib_size_data, 
                aes(x = Sample, y = LibrarySize, fill = Type)) +
      geom_boxplot() +
      scale_fill_brewer(palette = "Set2") +
      labs(title = "Library Size Distribution Before/After Normalization",
           x = "Sample",
           y = "Log2(Library Size + 1)") +
      get_theme_rna_processing() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    return(ggplotly(p))
  }, error = function(e) {
    log_plotting_error(e, "normalization_boxplots")
    return(create_error_plot("Could not create library size boxplots"))
  })
}

#' @title Plot Normalization Density
#' @description Creates density plots for normalized data
#' @param eval_results List containing normalization evaluation results
#' @return A ggplot object
#' @export
plot_normalization_density <- function(eval_results) {
  tryCatch({
    # Add debug output
    cat("Attempting to create normalization density plot\n")
    cat("Structure of eval_results:\n")
    str(eval_results, max.level = 2)
    
    # Check if mean_variance data exists
    if (is.null(eval_results) || is.null(eval_results$mean_variance) || 
        is.null(eval_results$mean_variance$data)) {
      warning("Mean-variance data not found in evaluation results")
      return(create_error_plot("Mean-variance data not available"))
    }
    
    # Extract data
    plot_data <- eval_results$mean_variance$data
    
    # Create plot
    p <- ggplot(plot_data, 
                aes(x = MeanExpr, y = VarianceExpr, color = Type)) +
      geom_point(alpha = 0.3, size = 1) +
      geom_smooth(method = "loess", se = FALSE) +
      scale_color_brewer(palette = "Set2") +
      scale_x_continuous(labels = comma_format()) +
      scale_y_continuous(labels = comma_format()) +
      labs(title = "Mean-Variance Relationship",
           x = "Log2(Mean Expression)",
           y = "Log2(Variance)") +
      get_theme_rna_processing()
    
    return(p)
    
  }, error = function(e) {
    cat("Error in normalization density plot:", e$message, "\n")
    return(create_error_plot(sprintf("Error creating normalization density plot: %s", e$message)))
  })
}

#' @title Plot Normalization PCA
#' @description Creates PCA plot for normalized data
#' @param eval_results List containing normalization evaluation results
#' @return A ggplot object
#' @export
plot_normalization_pca <- function(eval_results) {
  tryCatch({
    # Add debug output
    cat("Attempting to create normalization PCA plot\n")
    cat("Structure of eval_results:\n")
    str(eval_results, max.level = 2)
    
    # Check if batch effect data exists
    if (is.null(eval_results) || is.null(eval_results$batch_effect) || 
        is.null(eval_results$batch_effect$pca_data)) {
      warning("Batch effect data not found in evaluation results")
      return(create_error_plot("Batch effect data not available"))
    }
    
    # Extract data
    pca_data <- eval_results$batch_effect$pca_data
    variance_explained <- eval_results$batch_effect$variance_explained
    
    # Print detailed info for debugging
    cat("PCA data dimensions:", dim(pca_data), "\n")
    cat("Batch groups:", levels(pca_data$BatchGroup), "\n")
    cat("Number of samples per batch group:\n")
    print(table(pca_data$BatchGroup))
    
    # Ensure BatchGroup is a factor with proper levels
    pca_data$BatchGroup <- factor(pca_data$BatchGroup)
    
    # Check if we have enough data points for ellipses
    batch_counts <- table(pca_data$BatchGroup)
    can_draw_ellipses <- any(batch_counts >= 3)
    
    # Get a better color palette for many groups
    n_groups <- length(unique(pca_data$BatchGroup))
    if (n_groups > 8) {
      # Use a color palette that can handle many groups
      color_palette <- colorRampPalette(brewer.pal(8, "Dark2"))(n_groups)
    } else {
      color_palette <- brewer.pal(max(3, n_groups), "Set2")[1:n_groups]
    }
    
    # Create plot
    p <- ggplot(pca_data, 
                aes(x = PC1, y = PC2, 
                    color = BatchGroup, 
                    text = Sample)) +
      geom_point(size = 3)
    
    # Only add ellipses if there are enough points
    if (can_draw_ellipses) {
      # Identify groups with enough samples for ellipses
      valid_groups <- names(batch_counts[batch_counts >= 3])
      
      # Add ellipses only for groups with enough samples
      if (length(valid_groups) > 0) {
        p <- p + stat_ellipse(data = subset(pca_data, BatchGroup %in% valid_groups),
                            aes(group = BatchGroup), 
                            type = "norm", 
                            level = 0.95, 
                            show.legend = FALSE)
      }
      cat("Added ellipses for groups with sufficient samples\n")
    } else {
      cat("Not enough points per batch to calculate ellipses\n")
    }
    
    p <- p + 
      scale_color_manual(values = color_palette) +
      labs(title = "PCA of Normalized Data",
           x = sprintf("PC1 (%s%%)", round(variance_explained[1], 1)),
           y = sprintf("PC2 (%s%%)", round(variance_explained[2], 1)),
           color = "Batch Group") +
      get_theme_rna_processing() +
      theme(
        legend.position = "right",
        legend.key.size = unit(0.5, "cm")
      )
    
    return(p)
    
  }, error = function(e) {
    cat("Error in normalization PCA plot:", e$message, "\n")
    return(create_error_plot(sprintf("Error creating PCA plot: %s", e$message)))
  })
}

#' @title Print Normalization Statistics
#' @description Prints various statistics about the normalization results
#' @param eval_results List containing evaluation results
#' @return Printed output (invisible)
#' @export
print_normalization_stats <- function(eval_results) {
  tryCatch({
    # Check if eval_results has the required structure
    if (is.null(eval_results) || !is.list(eval_results)) {
      cat("No normalization evaluation results available\n")
      return(invisible(NULL))
    }
    
    cat("=== Normalization Statistics ===\n\n")
    
    # Library Size Statistics
    if (!is.null(eval_results$library_size) && 
        !is.null(eval_results$library_size$stats) &&
        !is.null(eval_results$library_size$stats$cv_raw) &&
        !is.null(eval_results$library_size$stats$cv_norm) &&
        is.numeric(eval_results$library_size$stats$cv_raw) &&
        is.numeric(eval_results$library_size$stats$cv_norm)) {
      
      cat("Library Size Statistics:\n")
      cat("- CV before normalization:", 
          round(eval_results$library_size$stats$cv_raw, 3), "\n")
      cat("- CV after normalization:", 
          round(eval_results$library_size$stats$cv_norm, 3), "\n\n")
    } else {
      cat("Library Size Statistics: Not available\n\n")
    }
    
    # Mean-Variance Relationship
    if (!is.null(eval_results$mean_variance) && 
        !is.null(eval_results$mean_variance$correlation) &&
        !is.null(eval_results$mean_variance$correlation$raw) &&
        !is.null(eval_results$mean_variance$correlation$normalized) &&
        is.numeric(eval_results$mean_variance$correlation$raw) &&
        is.numeric(eval_results$mean_variance$correlation$normalized)) {
      
      cat("Mean-Variance Correlation:\n")
      cat("- Before normalization:", 
          round(eval_results$mean_variance$correlation$raw, 3), "\n")
      cat("- After normalization:", 
          round(eval_results$mean_variance$correlation$normalized, 3), "\n\n")
    } else {
      cat("Mean-Variance Correlation: Not available\n\n")
    }
    
    # Batch Effect
    if (!is.null(eval_results$batch_effect) && 
        !is.null(eval_results$batch_effect$batch_r2) &&
        is.numeric(eval_results$batch_effect$batch_r2)) {
      
      cat("Batch Effect Size (R-squared):", 
          round(eval_results$batch_effect$batch_r2, 3), "\n\n")
    } else {
      cat("Batch Effect Size: Not available\n\n")
    }
    
    # Normality Statistics (showing first few samples)
    if (!is.null(eval_results$normality) && 
        !is.null(eval_results$normality$shapiro_test_results) &&
        is.matrix(eval_results$normality$shapiro_test_results) &&
        ncol(eval_results$normality$shapiro_test_results) > 0 &&
        nrow(eval_results$normality$shapiro_test_results) > 0 &&
        !all(is.na(eval_results$normality$shapiro_test_results))) {
      
      cat("Normality Statistics (first few samples):\n")
      n_cols <- min(5, ncol(eval_results$normality$shapiro_test_results))
      norm_stats <- eval_results$normality$shapiro_test_results[, 1:n_cols, drop = FALSE]
      
      # Check if we have valid numeric values before printing
      if (all(sapply(norm_stats, is.numeric))) {
        print(round(norm_stats, 3))
      } else {
        cat("Normality test results contain non-numeric values\n")
        # Print what we can
        for (col in colnames(norm_stats)) {
          cat(sprintf("%s: ", col))
          if (all(is.na(norm_stats[, col]))) {
            cat("All NA values\n")
          } else {
            cat("Mixed or invalid values\n")
          }
        }
      }
    } else {
      cat("Normality Statistics: Not available or contains only NA values\n")
    }
    
    invisible(NULL)
  }, error = function(e) {
    cat("Error printing normalization statistics:", e$message, "\n")
  })
}

#' @title Get Library Size Statistics
#' @description Returns formatted library size statistics
#' @param eval_results List containing evaluation results
#' @return Character string with statistics
#' @export
get_lib_size_stats <- function(eval_results) {
  tryCatch({
    stats <- eval_results$library_size$stats
    sprintf(
      paste0("Library Size Statistics:\n\n",
      "Coefficient of Variation (CV):\n",
      "- Before normalization: %.3f\n",
      "- After normalization: %.3f\n\n",
      "Improvement: %.1f%%"),
      stats$cv_raw,
      stats$cv_norm,
      100 * (stats$cv_raw - stats$cv_norm) / stats$cv_raw
    )
  }, error = function(e) {
    return("Error calculating library size statistics")
  })
}

#' @title Get Effects Statistics
#' @description Returns formatted normalization effects statistics
#' @param eval_results List containing evaluation results
#' @return Character string with statistics
#' @export
get_effects_stats <- function(eval_results) {
  tryCatch({
    # Check if eval_results has the required structure
    if (is.null(eval_results) || !is.list(eval_results) || 
        is.null(eval_results$mean_variance) || !is.list(eval_results$mean_variance) ||
        is.null(eval_results$mean_variance$correlation)) {
      return("Mean-variance relationship data not available")
    }
    
    corr <- eval_results$mean_variance$correlation
    
    # Check if correlation values are numeric
    if (!is.numeric(corr$raw) || !is.numeric(corr$normalized)) {
      return("Mean-variance correlation values are not numeric")
    }
    
    # Calculate percent change safely
    percent_change <- if (corr$raw != 0) {
      100 * abs(corr$raw - corr$normalized) / abs(corr$raw)
    } else {
      NA_real_
    }
    
    # Format the output with NA handling
    result <- sprintf(
      paste0("Mean-Variance Relationship:\n\n",
      "Correlation coefficient:\n",
      "- Before normalization: %.3f\n",
      "- After normalization: %.3f\n"),
      corr$raw,
      corr$normalized
    )
    
    # Add percent change if available
    if (!is.na(percent_change)) {
      result <- paste0(result, sprintf("\nChange: %.1f%%", percent_change))
    }
    
    return(result)
  }, error = function(e) {
    return(paste("Error calculating normalization effects statistics:", e$message))
  })
}

#' @title Get Batch Effect Statistics
#' @description Returns formatted batch effect statistics
#' @param eval_results List containing evaluation results
#' @return Character string with statistics
#' @export
get_batch_stats <- function(eval_results) {
  tryCatch({
    # Check if eval_results has the required structure
    if (is.null(eval_results) || !is.list(eval_results) || 
        is.null(eval_results$batch_effect) || !is.list(eval_results$batch_effect) ||
        is.null(eval_results$batch_effect$batch_r2)) {
      return("Batch effect data not available")
    }
    
    r2 <- eval_results$batch_effect$batch_r2
    
    # Check if R² is numeric
    if (!is.numeric(r2)) {
      return("Batch effect R² value is not numeric")
    }
    
    sprintf(
      paste0("Batch Effect Analysis:\n\n",
      "R-squared value: %.3f\n\n",
      "Interpretation:\n",
      "- R² < 0.1: Low batch effect\n",
      "- 0.1 ≤ R² < 0.3: Moderate batch effect\n",
      "- R² ≥ 0.3: Strong batch effect"),
      r2
    )
  }, error = function(e) {
    return(paste("Error calculating batch effect statistics:", e$message))
  })
}

# Library Size Plotting Functions

#' Plot library size distribution
#' @param lib_size_data data.frame with library sizes pre/post normalization
#' @return A ggplot object
#' @export
plot_library_size_distribution <- function(lib_size_data) {
  tryCatch({
    # Ensure data is properly formatted
    if (!all(c("Sample", "LibrarySize", "Type") %in% colnames(lib_size_data))) {
      stop("Required columns not found in lib_size_data: Sample, LibrarySize, Type")
    }
    
    # Calculate summary statistics for each sample and type
    lib_stats <- lib_size_data %>%
      group_by(Sample, Type) %>%
      summarise(
        Mean = mean(LibrarySize),
        SD = sd(LibrarySize),
        .groups = 'drop'
      )
    
    # Create the plot
    p <- ggplot(lib_stats, 
                aes(x = reorder(Sample, -Mean), 
                    y = Mean, 
                    fill = Type)) +
      geom_bar(stat = "identity", 
               position = position_dodge(width = 0.8),
               width = 0.7,
               alpha = 0.7) +
      geom_errorbar(aes(ymin = Mean - SD, 
                       ymax = Mean + SD),
                   position = position_dodge(width = 0.8),
                   width = 0.25) +
      scale_fill_brewer(palette = "Set2") +
      labs(title = "Library Size Distribution",
           x = "Sample",
           y = "Library Size (log2)") +
      get_theme_rna_processing() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_line(color = "gray95"),
        plot.title = element_text(hjust = 0.5)
      )
    
    return(p)
    
  }, error = function(e) {
    warning("Error in library size distribution plot:", e$message)
    return(NULL)
  })
}

#' Format library size statistics
#' @param stats List containing library size statistics
#' @return Formatted string with statistics
#' @export
format_library_size_stats <- function(stats) {
  tryCatch({
    # Check if stats has the required structure
    if (is.null(stats) || !is.list(stats) || 
        !all(c("raw", "normalized") %in% names(stats))) {
      return("Library size statistics not available in the expected format")
    }
    
    # Check if raw and normalized stats have the required fields
    required_fields <- c("mean", "median", "cv")
    if (!all(required_fields %in% names(stats$raw)) || 
        !all(required_fields %in% names(stats$normalized))) {
      return("Library size statistics missing required fields (mean, median, cv)")
    }
    
    # Check if values are numeric
    for (field in required_fields) {
      if (!is.numeric(stats$raw[[field]]) || !is.numeric(stats$normalized[[field]])) {
        return(paste("Non-numeric", field, "value found in library size statistics"))
      }
    }
    
    # Format numbers with appropriate scaling
    format_number <- function(x) {
      if (is.null(x) || !is.numeric(x) || !is.finite(x)) return("NA")
      if (x >= 1e6) {
        sprintf("%.2fM", x/1e6)
      } else if (x >= 1e3) {
        sprintf("%.2fK", x/1e3)
      } else {
        sprintf("%.2f", x)
      }
    }
    
    # Calculate impact metrics safely
    cv_reduction <- if (stats$raw$cv != 0) {
      100 * (stats$raw$cv - stats$normalized$cv) / stats$raw$cv
    } else {
      NA_real_
    }
    
    var_stabilization <- if (stats$raw$cv != 0) {
      100 * (1 - stats$normalized$cv/stats$raw$cv)
    } else {
      NA_real_
    }
    
    # Create the formatted string
    stats_text <- sprintf(
      paste0("Library Size Statistics:\n\n%s\n%s\n"),
      
      # Raw data section
      sprintf("Raw Data:\n- Mean: %s\n- Median: %s\n- CV: %.3f",
              format_number(stats$raw$mean),
              format_number(stats$raw$median),
              stats$raw$cv),
      
      # Normalized data section
      sprintf("Normalized Data:\n- Mean: %s\n- Median: %s\n- CV: %.3f",
              format_number(stats$normalized$mean),
              format_number(stats$normalized$median),
              stats$normalized$cv)
    )
    
    # Add impact section if metrics are available
    if (!is.na(cv_reduction) && !is.na(var_stabilization)) {
      stats_text <- paste0(stats_text, 
                          sprintf("\nNormalization Impact:\n- CV Reduction: %.1f%%\n- Variance Stabilization: %.1f%%",
                                 cv_reduction, var_stabilization))
    }
    
    return(stats_text)
    
  }, error = function(e) {
    return(paste("Error formatting library size statistics:", e$message))
  })
}

#' @title Plot Filtering Gene Counts
#' @description Creates a bar plot comparing gene counts before and after filtering
#' @param initial_genes Number of initial genes
#' @param filtered_genes Number of genes after filtering
#' @return A ggplot object
#' @export
plot_filtering_gene_counts <- function(initial_genes, filtered_genes) {
  tryCatch({
    # Input validation
    if (is.null(initial_genes) || is.null(filtered_genes)) {
      stop("Both initial_genes and filtered_genes must be provided")
    }
    
    # Calculate removed genes and percentages
    removed_genes <- initial_genes - filtered_genes
    kept_percent <- round(filtered_genes/initial_genes * 100, 1)
    removed_percent <- round(removed_genes/initial_genes * 100, 1)
    
    # Create data frame for plotting
    plot_data <- data.frame(
      Stage = factor(c("Initial", "After Filtering", "Removed"),
                    levels = c("Initial", "After Filtering", "Removed")),
      Count = c(initial_genes, filtered_genes, removed_genes),
      Percentage = c(100, kept_percent, removed_percent)
    )
    
    # Create custom color palette
    stage_colors <- c("#66C2A5", "#FC8D62", "#8DA0CB")
    
    # Create plot
    p <- ggplot(plot_data, aes(x = Stage, y = Count, fill = Stage)) +
      geom_bar(stat = "identity", width = 0.7) +
      geom_text(aes(label = sprintf("%s\n(%.1f%%)", 
                                   format(Count, big.mark = ",", scientific = FALSE),
                                   Percentage)),
                position = position_stack(vjust = 0.5),
                size = 4) +
      scale_fill_manual(values = stage_colors) +
      scale_y_continuous(labels = scales::comma,
                        expand = expansion(mult = c(0, 0.2))) +
      labs(title = "Gene Counts Before and After Filtering",
           subtitle = sprintf("Removed %s genes (%.1f%%)", 
                            format(removed_genes, big.mark = ","), 
                            removed_percent),
           y = "Number of Genes",
           x = "") +
      get_theme_rna_processing() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.subtitle = element_text(size = 10, color = "gray30")
      )
    
    return(p)
    
  }, error = function(e) {
    warning("Error in filtering gene counts plot:", e$message)
    return(NULL)
  })
}

#' @title Plot Expression Distribution
#' @description Creates line plots comparing expression distribution before and after filtering
#' @param before_stats Named numeric vector of expression statistics before filtering
#' @param after_stats Named numeric vector of expression statistics after filtering
#' @return A ggplot object
#' @export
plot_expression_distribution <- function(before_stats, after_stats) {
  tryCatch({
    # Input validation
    if (is.null(before_stats) || is.null(after_stats)) {
      stop("Both before_stats and after_stats must be provided")
    }
    
    # Handle named vectors and ensure we have numeric values
    if (is.list(before_stats) && "quantiles" %in% names(before_stats)) {
      before_stats <- before_stats$quantiles
    }
    
    if (is.list(after_stats) && "quantiles" %in% names(after_stats)) {
      after_stats <- after_stats$quantiles
    }
    
    # Extract values and ensure they are numeric
    before_values <- as.numeric(before_stats)
    after_values <- as.numeric(after_stats)
    
    # Get percentile names, either from the named vector or create default ones
    if (!is.null(names(before_stats)) && length(names(before_stats)) == length(before_values)) {
      percentiles <- names(before_stats)
    } else {
      percentiles <- c("10%", "25%", "50%", "75%", "90%")
    }
    
    # Create data frame for plotting
    plot_data <- data.frame(
      Percentile = factor(percentiles, levels = percentiles),
      Before = before_values,
      After = after_values
    )
    
    # Calculate fold changes
    plot_data$FoldChange <- after_values / before_values
    
    # Reshape data for plotting
    plot_data_long <- tidyr::gather(plot_data, 
                                  key = "Stage", 
                                  value = "Expression", 
                                  Before:After)
    
    # Custom colors
    stage_colors <- c("#66C2A5", "#FC8D62")
    
    # Create plot
    p <- ggplot(plot_data_long, 
                aes(x = Percentile, y = Expression, 
                    color = Stage, group = Stage)) +
      geom_line(size = 1) +
      geom_point(size = 3) +
      geom_text(data = plot_data,
                aes(x = Percentile, 
                    y = pmax(Before, After) * 1.1,
                    label = sprintf("%.1fx", FoldChange)),
                color = "gray30",
                size = 3) +
      scale_color_manual(values = stage_colors,
                        labels = c("Before Filtering", "After Filtering")) +
      scale_y_log10(labels = scales::comma,
                    expand = expansion(mult = c(0.1, 0.2))) +
      labs(title = "Expression Distribution Before vs After Filtering",
           subtitle = "Fold change shown above each percentile",
           y = "Expression Level (log10 scale)",
           x = "Percentile",
           color = "Stage") +
      get_theme_rna_processing() +
      theme(
        legend.position = "top",
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.subtitle = element_text(size = 10, color = "gray30")
      )
    
    return(p)
    
  }, error = function(e) {
    warning("Error in expression distribution plot:", e$message)
    return(create_error_plot(sprintf("Error: %s", e$message)))
  })
}

#' @title Plot Expression Density
#' @description Creates density plots comparing expression distribution before and after filtering
#' @param before_data Numeric vector or matrix of expression values before filtering
#' @param after_data Numeric vector or matrix of expression values after filtering
#' @return A ggplot object
#' @export
plot_filtering_density <- function(before_data, after_data) {
  tryCatch({
    # Input validation
    if (is.null(before_data) || is.null(after_data)) {
      stop("Both before_data and after_data must be provided")
    }
    
    # Convert to vectors if matrices or data frames
    if (is.matrix(before_data) || is.data.frame(before_data)) {
      before_data <- as.vector(as.matrix(before_data))
    }
    
    if (is.matrix(after_data) || is.data.frame(after_data)) {
      after_data <- as.vector(as.matrix(after_data))
    }
    
    # Apply log2 transformation (adding small value to avoid log(0))
    before_log2 <- log2(before_data + 1)
    after_log2 <- log2(after_data + 1)
    
    # Create data frame for plotting
    plot_data <- data.frame(
      Expression = c(before_log2, after_log2),
      Type = factor(c(rep("Before Filtering", length(before_log2)),
                     rep("After Filtering", length(after_log2))),
                   levels = c("Before Filtering", "After Filtering"))
    )
    
    # Create plot
    p <- ggplot(plot_data, aes(x = Expression, fill = Type, color = Type)) +
      geom_density(alpha = 0.5) +
      scale_fill_manual(values = c("Before Filtering" = "#66C2A5", "After Filtering" = "#FC8D62")) +
      scale_color_manual(values = c("Before Filtering" = "#66C2A5", "After Filtering" = "#FC8D62")) +
      labs(title = "Expression Distribution Before vs After Filtering",
           x = "Log2(Expression + 1)",
           y = "Density") +
      get_theme_rna_processing() +
      theme(
        legend.position = "top",
        plot.title = element_text(hjust = 0.5)
      )
    
    return(p)
    
  }, error = function(e) {
    warning("Error in filtering density plot:", e$message)
    return(create_error_plot(sprintf("Error: %s", e$message)))
  })
}

#' @title Plot QQ
#' @description Creates a Q-Q plot for normality assessment
#' @param qq_data Data frame with theoretical and sample quantiles
#' @return A ggplot object
#' @export
plot_qq <- function(qq_data) {
  tryCatch({
    # Add debug output
    cat("Attempting to create QQ plot\n")
    
    # Check if qq_data exists and has required columns
    if (is.null(qq_data) || !is.data.frame(qq_data) ||
        !all(c("Sample", "TheoreticalQuantiles", "SampleQuantiles") %in% colnames(qq_data))) {
      warning("QQ data not available or missing required columns")
      return(create_error_plot("QQ plot data not available in the expected format"))
    }
    
    # Check if we have enough samples to plot
    unique_samples <- unique(qq_data$Sample)
    n_samples <- length(unique_samples)
    
    if (n_samples == 0) {
      return(create_error_plot("No samples available for QQ plot"))
    }
    
    # Limit to first 9 samples if there are more
    if (n_samples > 9) {
      qq_data <- qq_data[qq_data$Sample %in% unique_samples[1:9], ]
      cat("Limiting QQ plot to first 9 samples\n")
    }
    
    # Create plot
    p <- ggplot(qq_data, 
                aes(x = TheoreticalQuantiles, y = SampleQuantiles, color = Sample)) +
      geom_point(size = 1, alpha = 0.7) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
      facet_wrap(~ Sample, scales = "free") +
      labs(title = "Normal Q-Q Plot",
           x = "Theoretical Quantiles",
           y = "Sample Quantiles") +
      get_theme_rna_processing() +
      theme(
        legend.position = "none",
        strip.background = element_rect(fill = "gray95"),
        strip.text = element_text(face = "bold")
      )
    
    return(p)
    
  }, error = function(e) {
    cat("Error in QQ plot:", e$message, "\n")
    return(create_error_plot(sprintf("Error creating QQ plot: %s", e$message)))
  })
} 