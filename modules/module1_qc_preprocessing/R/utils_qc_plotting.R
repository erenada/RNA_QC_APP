# utils_qc_plotting.R
# QC & Pre-processing Tool
# Utility Functions for QC Plotting and Analysis
# Author: Eren Ada, PhD
# Date: 05/13/2024

#' @import shiny
#' @import ggplot2
#' @import plotly
#' @import pheatmap
#' @import DESeq2
#' @import dplyr
#' @import moments
#' @import SummarizedExperiment
#' @importFrom ggrepel geom_text_repel
#' @importFrom stats cor cor.test prcomp p.adjust shapiro.test qqnorm qqline density
#' @importFrom matrixStats rowVars
#' @importFrom RColorBrewer brewer.pal
#' @import grid
#' @importFrom grDevices colorRampPalette
NULL

# --- 1. General Metrics & Plots Helper Functions ---

#' Generate Library Size Distribution Plot
#' @param raw_counts_data Matrix of raw counts (genes x samples)
#' @param current_view_data Matrix of currently processed counts
#' @param plot_type Character string, "boxplot" or "barplot"
#' @param use_log_scale Logical, whether to use log scale for y-axis
#' @return A plotly object
#' @export
plot_library_sizes_func <- function(raw_counts_data, current_view_data, plot_type = "boxplot", use_log_scale = FALSE) {
  # Calculate library sizes
  raw_lib_sizes <- colSums(raw_counts_data)
  current_lib_sizes <- colSums(current_view_data)
  
  # Create data frame for plotting
  plot_df <- data.frame(
    Sample = rep(names(raw_lib_sizes), 2),
    LibrarySize = c(raw_lib_sizes, current_lib_sizes),
    DataType = rep(c("Raw", "Processed"), each = length(raw_lib_sizes))
  )
  
  # Create base plot
  p <- ggplot(plot_df, aes(x = Sample, y = LibrarySize, fill = DataType))
  
  # Add appropriate geometry based on plot type
  if (plot_type == "boxplot") {
    p <- p + geom_boxplot()
  } else {
    p <- p + geom_bar(stat = "identity", position = "dodge")
  }
  
  # Add styling
  p <- p + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Library Size Distribution",
         y = "Library Size",
         x = "Sample")
  
  # Apply log scale if requested
  if (use_log_scale) {
    p <- p + scale_y_log10()
  }
  
  # Convert to plotly
  ggplotly(p)
}

#' Generate Gene Detection Rates Plot
#' @param counts_data Matrix of counts
#' @param detection_threshold Numeric threshold
#' @return A plotly object
#' @export
plot_gene_detection_func <- function(counts_data, detection_threshold = 1) {
  # Calculate detection rates
  detected_genes <- colSums(counts_data > detection_threshold)
  total_genes <- nrow(counts_data)
  detection_rates <- detected_genes / total_genes * 100
  
  # Create data frame for plotting
  plot_df <- data.frame(
    Sample = names(detected_genes),
    DetectedGenes = detected_genes,
    DetectionRate = detection_rates
  )
  
  # Create plot
  p <- ggplot(plot_df, aes(x = Sample, y = DetectionRate)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0("Gene Detection Rates (threshold = ", detection_threshold, ")"),
         y = "Detection Rate (%)",
         x = "Sample")
  
  # Convert to plotly
  ggplotly(p)
}

#' Calculate QC Summary Statistics
#' @param raw_counts_data Raw count matrix
#' @param current_view_data Processed count matrix
#' @param detection_threshold Detection threshold
#' @return List of summary statistics
#' @export
calculate_qc_summary_stats_func <- function(raw_counts_data, current_view_data, detection_threshold = 1) {
  # Calculate library size statistics
  raw_lib_sizes <- colSums(raw_counts_data)
  current_lib_sizes <- colSums(current_view_data)
  
  # Calculate detection statistics
  detected_genes <- colSums(current_view_data > detection_threshold)
  detection_rates <- detected_genes / nrow(current_view_data) * 100
  
  # Compile statistics
  list(
    raw_library_sizes = list(
      min = min(raw_lib_sizes),
      max = max(raw_lib_sizes),
      mean = mean(raw_lib_sizes),
      median = median(raw_lib_sizes)
    ),
    current_library_sizes = list(
      min = min(current_lib_sizes),
      max = max(current_lib_sizes),
      mean = mean(current_lib_sizes),
      median = median(current_lib_sizes)
    ),
    gene_detection = list(
      min_rate = min(detection_rates),
      max_rate = max(detection_rates),
      mean_rate = mean(detection_rates),
      median_rate = median(detection_rates),
      threshold_used = detection_threshold
    )
  )
}

# --- 2. Sample Similarity Plots Helper Functions ---

#' Perform PCA on Count Data
#' @param counts_data Count matrix
#' @param center Center the data
#' @param scale. Scale the data
#' @param n_top_genes Number of top variable genes
#' @return List with PCA results
#' @export
perform_pca_func <- function(counts_data, center = TRUE, scale. = TRUE, n_top_genes = NULL) {
  # Select top variable genes if specified
  if (!is.null(n_top_genes)) {
    gene_vars <- rowVars(as.matrix(counts_data))
    top_genes <- order(gene_vars, decreasing = TRUE)[1:min(n_top_genes, length(gene_vars))]
    counts_data <- counts_data[top_genes, ]
  }
  
  # Check if any rows have zero variance and remove them to avoid PCA errors
  row_vars <- apply(counts_data, 1, var)
  if(any(row_vars == 0)) {
    message("Removing ", sum(row_vars == 0), " rows with zero variance before PCA")
    counts_data <- counts_data[row_vars > 0, ]
  }
  
  # Perform PCA
  pca_res <- prcomp(t(counts_data), center = center, scale. = scale.)
  
  # Calculate variance explained
  var_explained <- (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100
  
  # Return results
  list(
    pca_object = pca_res,
    var_explained = var_explained
  )
}

#' Generate PCA Plot
#' @param pca_results PCA results from perform_pca_func
#' @param pc_x_choice X-axis PC
#' @param pc_y_choice Y-axis PC
#' @param pc_z_choice Z-axis PC (for 3D plot)
#' @param metadata_df Sample metadata
#' @param color_var Color variable
#' @param show_labels Show sample labels
#' @param plot_type Either "2d" or "3d"
#' @return A plotly object
#' @export
plot_pca_func <- function(pca_results, pc_x_choice, pc_y_choice, pc_z_choice = NULL,
                         metadata_df = NULL, color_var = NULL,
                         show_labels = FALSE, plot_type = "2d") {
  # Validate inputs
  if (is.null(pca_results) || is.null(pca_results$pca_object) || is.null(pca_results$var_explained)) {
    stop("Invalid PCA results provided")
  }
  
  # Extract PC scores
  pc_data <- as.data.frame(pca_results$pca_object$x)
  
  # Get PC numbers
  pc_x_num <- as.numeric(gsub("PC", "", pc_x_choice))
  pc_y_num <- as.numeric(gsub("PC", "", pc_y_choice))
  
  # Get variance explained
  var_exp <- pca_results$var_explained
  
  # Validate PC numbers
  if (pc_x_num > ncol(pc_data) || pc_y_num > ncol(pc_data)) {
    stop("Requested PC number exceeds available components")
  }
  
  # Create base plot data
  plot_data <- data.frame(
    PC1 = pc_data[, pc_x_num],
    PC2 = pc_data[, pc_y_num],
    Sample = rownames(pc_data)
  )
  
  # Initialize group column
  plot_data$group <- "All Samples"
  
  # Add group information if available and valid
  if (!is.null(color_var) && !is.null(metadata_df)) {
    if (color_var == "Sample") {
      plot_data$group <- plot_data$Sample
    } else if (color_var != "none" && color_var %in% colnames(metadata_df)) {
      # Ensure sample names are character type
      plot_data$Sample <- as.character(plot_data$Sample)
      # Get the color variable values, ensuring proper matching of sample names
      matched_samples <- intersect(plot_data$Sample, rownames(metadata_df))
      if (length(matched_samples) > 0) {
        plot_data$group[plot_data$Sample %in% matched_samples] <- 
          metadata_df[matched_samples, color_var]
      }
    }
  }
  
  # Handle NA values in group
  plot_data$group[is.na(plot_data$group)] <- "NA"
  
  if (plot_type == "2d") {
    # Create the basic 2D plot with ggplot2
    p <- ggplot(plot_data, aes(x = PC1, y = PC2, color = group, text = Sample)) +
      geom_point(size = 3) +
      labs(title = "PCA Plot",
           x = sprintf("PC%d (%.1f%%)", pc_x_num, var_exp[pc_x_num]),
           y = sprintf("PC%d (%.1f%%)", pc_y_num, var_exp[pc_y_num])) +
      theme_minimal() +
      theme(legend.title = element_blank())
    
    # Add labels if requested
    if (show_labels) {
      p <- p + geom_text_repel(aes(label = Sample), size = 3, box.padding = 0.5, show.legend = FALSE)
    }
    
    # Convert to plotly with proper hover text
    p <- ggplotly(p, tooltip = c("text", "x", "y"))
    
  } else {
    # Validate PC3 for 3D plot
    if (is.null(pc_z_choice)) {
      stop("Z-axis PC choice required for 3D plot")
    }
    pc_z_num <- as.numeric(gsub("PC", "", pc_z_choice))
    if (pc_z_num > ncol(pc_data)) {
      stop("Requested PC number for Z-axis exceeds available components")
    }
    
    # Add PC3 to plot data
    plot_data$PC3 <- pc_data[, pc_z_num]
    
    # Create 3D plot directly with plotly
    p <- plot_ly(plot_data,
                x = ~PC1, y = ~PC2, z = ~PC3,
                color = ~group,
                text = ~Sample,
                type = "scatter3d",
                mode = "markers",
                marker = list(size = 5)) %>%
      layout(scene = list(
        xaxis = list(title = sprintf("PC%d (%.1f%%)", pc_x_num, var_exp[pc_x_num])),
        yaxis = list(title = sprintf("PC%d (%.1f%%)", pc_y_num, var_exp[pc_y_num])),
        zaxis = list(title = sprintf("PC%d (%.1f%%)", pc_z_num, var_exp[pc_z_num]))
      )) %>%
      layout(title = "3D PCA Plot",
             showlegend = TRUE)
  }
  
  return(p)
}

#' Calculate Sample Correlation Matrix
#' @param data_matrix Count matrix
#' @param method Correlation method
#' @return List with correlation results
#' @export
calculate_correlation_func <- function(data_matrix, method = "pearson") {
  # Calculate correlation matrix
  cor_matrix <- cor(data_matrix, method = method)
  
  # Calculate p-values with appropriate handling for ties
  n_samples <- ncol(data_matrix)
  p_matrix <- matrix(NA, n_samples, n_samples)
  
  for(i in 1:n_samples) {
    for(j in 1:n_samples) {
      if(i != j) {
        # Use try-catch to handle potential warnings and errors
        test_result <- tryCatch({
          if (method == "spearman") {
            # For Spearman, use asymptotic approximation which handles ties better
            stats::cor.test(data_matrix[,i], data_matrix[,j], 
                          method = method, 
                          exact = FALSE,
                          continuity = TRUE)
          } else {
            # For Pearson, no need for exact test
            stats::cor.test(data_matrix[,i], data_matrix[,j], 
                          method = method)
          }
        }, warning = function(w) {
          # Return the test anyway if there's just a warning
          if (method == "spearman") {
            stats::cor.test(data_matrix[,i], data_matrix[,j], 
                          method = method, 
                          exact = FALSE,
                          continuity = TRUE)
          } else {
            stats::cor.test(data_matrix[,i], data_matrix[,j], 
                          method = method)
          }
        }, error = function(e) {
          # Return NA if there's an actual error
          list(p.value = NA)
        })
        
        p_matrix[i,j] <- test_result$p.value
      }
    }
  }
  
  # Adjust p-values for multiple testing
  p_adjusted <- matrix(p.adjust(p_matrix, method = "BH"), n_samples, n_samples)
  
  # Set row and column names
  rownames(cor_matrix) <- colnames(cor_matrix) <- colnames(data_matrix)
  rownames(p_adjusted) <- colnames(p_adjusted) <- colnames(data_matrix)
  
  return(list(
    correlation = cor_matrix,
    p_values = p_matrix,
    p_adjusted = p_adjusted
  ))
}

#' Generate Sample Correlation Heatmap
#' @param cor_results Correlation results
#' @param metadata_df Sample metadata
#' @param annotation_cols Annotation columns
#' @param color_scheme Color scheme
#' @param show_values Show correlation values
#' @param show_significance Show significance
#' @param auto_scale_colors Auto-scale colors for high correlations
#' @return A pheatmap plot
#' @export
plot_correlation_heatmap_func <- function(cor_results, metadata_df = NULL, 
                                        annotation_cols = NULL, color_scheme = "default",
                                        show_values = FALSE, show_significance = TRUE,
                                        auto_scale_colors = TRUE) {
  
  # Debug logging
  message("Debug - plot_correlation_heatmap_func:")
  message("Input correlation matrix dimensions: ", paste(dim(cor_results$correlation), collapse=" x "))
  if (!is.null(metadata_df)) {
    message("Input metadata dimensions: ", paste(dim(metadata_df), collapse=" x "))
    message("Annotation columns: ", paste(annotation_cols, collapse=", "))
  } else {
    message("No metadata provided for annotations")
  }
  message("Color scheme: ", color_scheme)
  message("Show values: ", show_values)
  message("Show significance: ", show_significance)
  message("Auto-scale colors: ", auto_scale_colors)
  
  # Validate correlation results
  if (is.null(cor_results) || is.null(cor_results$correlation)) {
    stop("Invalid correlation results: NULL or missing correlation matrix")
  }
  
  cor_matrix <- cor_results$correlation
  
  # Verify matrix format
  if (!is.matrix(cor_matrix)) {
    message("Warning: cor_matrix is not a matrix, converting...")
    cor_matrix <- as.matrix(cor_matrix)
  }
  if (!is.numeric(cor_matrix)) {
    stop("Error: cor_matrix must be numeric")
  }
  
  # Verify matrix dimensions
  if (nrow(cor_matrix) < 2 || ncol(cor_matrix) < 2) {
    stop("Error: correlation matrix must have at least 2x2 dimensions")
  }
  
  # Verify that rownames and colnames match
  if (!identical(rownames(cor_matrix), colnames(cor_matrix))) {
    message("Warning: row and column names in correlation matrix don't match, fixing...")
    # Use row names for both to ensure consistency
    colnames(cor_matrix) <- rownames(cor_matrix)
  }
  
  # Prepare color scheme with better handling of color palettes
  colors <- switch(color_scheme,
    "default" = colorRampPalette(c("#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0", 
                                  "#ffffff", "#fddbc7", "#f4a582", "#d6604d", 
                                  "#b2182b", "#67001f"))(200),
    "blue_red" = colorRampPalette(c("navy", "white", "firebrick3"))(200),
    "purple_orange" = colorRampPalette(c("#4A148C", "white", "#E65100"))(200),
    "ocean" = colorRampPalette(c("#006064", "white", "#1B5E20"))(200),
    # Default fallback
    colorRampPalette(c("#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0", 
                      "#ffffff", "#fddbc7", "#f4a582", "#d6604d", 
                      "#b2182b", "#67001f"))(200)
  )
  
  # Prepare annotations if provided
  annotation_data <- NULL
  annotation_colors <- NULL
  
  if (!is.null(metadata_df) && !is.null(annotation_cols) && length(annotation_cols) > 0) {
    tryCatch({
      # Check if all annotation columns exist in the metadata
      valid_cols <- annotation_cols[annotation_cols %in% colnames(metadata_df)]
      message("Valid annotation columns: ", paste(valid_cols, collapse=", "))
      
      if (length(valid_cols) > 0) {
        # Get samples that exist in both correlation matrix and metadata
        cor_samples <- rownames(cor_matrix)
        meta_samples <- rownames(metadata_df)
        
        if (length(meta_samples) == 0) {
          message("Warning: metadata has no rownames, cannot match samples")
          annotation_data <- NULL
        } else {
          common_samples <- intersect(cor_samples, meta_samples)
          message("Samples in correlation: ", length(cor_samples))
          message("Samples in metadata: ", length(meta_samples))
          message("Common samples: ", length(common_samples))
          
          if (length(common_samples) > 0) {
            # Create metadata subset with only common samples and valid columns
            annotation_data <- metadata_df[common_samples, valid_cols, drop = FALSE]
            
            # Ensure rownames are set correctly and match correlation matrix
            if (is.null(rownames(annotation_data))) {
              rownames(annotation_data) <- common_samples
            }
            
            # Create annotation colors
            annotation_colors <- list()
            for (col in valid_cols) {
              if (is.factor(annotation_data[[col]]) || is.character(annotation_data[[col]])) {
                unique_values <- unique(annotation_data[[col]])
                n_values <- length(unique_values)
                
                # Choose appropriate color palette based on number of values
                if (n_values <= 9) {
                  # Use qualitative palette for categorical data
                  pal <- brewer.pal(max(3, n_values), "Set1")[1:n_values]
                  names(pal) <- unique_values
                  annotation_colors[[col]] <- pal
                }
              }
            }
            
            message("Prepared annotation data dimensions: ", paste(dim(annotation_data), collapse=" x "))
          } else {
            message("No matching samples between correlation matrix and metadata")
            annotation_data <- NULL
          }
        }
      } else {
        message("No valid annotation columns found in metadata")
        annotation_data <- NULL
      }
    }, error = function(e) {
      message("Error preparing annotations: ", e$message)
      annotation_data <- NULL
      annotation_colors <- NULL
    })
  }
  
  # Prepare display numbers with proper handling of missing values
  display_numbers <- FALSE
  if (show_values) {
    if (show_significance && !is.null(cor_results$p_adjusted)) {
      display_numbers <- matrix(
        sprintf("%.2f%s",
                cor_matrix,
                ifelse(!is.na(cor_results$p_adjusted) & cor_results$p_adjusted < 0.05, "*", "")),
        nrow = nrow(cor_matrix)
      )
      rownames(display_numbers) <- rownames(cor_matrix)
      colnames(display_numbers) <- colnames(cor_matrix)
    } else {
      display_numbers <- round(cor_matrix, 2)
    }
  } else if (show_significance && !is.null(cor_results$p_adjusted)) {
    display_numbers <- matrix(
      ifelse(!is.na(cor_results$p_adjusted) & cor_results$p_adjusted < 0.05, "*", ""),
      nrow = nrow(cor_matrix)
    )
    rownames(display_numbers) <- rownames(cor_matrix)
    colnames(display_numbers) <- colnames(cor_matrix)
  }
  
  # Create heatmap using corrplot (simplified - no fallbacks needed)
  result <- tryCatch({
    message("Using corrplot for correlation heatmap rendering...")
    
    # Determine appropriate color range based on correlation values
    cor_range <- range(cor_matrix[upper.tri(cor_matrix)], na.rm = TRUE)
    message("Correlation range: ", sprintf("%.3f to %.3f", cor_range[1], cor_range[2]))
    message("Range span: ", sprintf("%.3f", cor_range[2] - cor_range[1]))
    message("Auto-scale colors setting: ", auto_scale_colors)
    
    # Create a display matrix for better visualization of high correlations
    display_matrix <- cor_matrix
    
    # Choose color scheme and scaling based on correlation values
    if ((cor_range[1] > 0.8 || (cor_range[1] > 0.7 && (cor_range[2] - cor_range[1]) < 0.3)) && auto_scale_colors) {
      # High correlations with small range - rescale for better visualization
      message("Using high-correlation rescaling for better visualization")
      message("Condition met: min_cor > 0.8 OR (min_cor > 0.7 AND range < 0.3)")
      
      # Calculate an appropriate rescaling range based on the actual data
      range_span <- cor_range[2] - cor_range[1]
      
      # If the range is very small (< 0.1), expand it more aggressively
      if (range_span < 0.1) {
        # Expand the range to show more color variation while keeping it realistic
        expansion_factor <- max(0.2, 0.1 / range_span)  # At least 0.2 expansion
        expanded_range <- range_span * expansion_factor
        
        # Center the expanded range around the actual range
        center_point <- (cor_range[1] + cor_range[2]) / 2
        cor_min_rescaled <- max(center_point - expanded_range/2, 0.3)  # Don't go below 0.3
        cor_max_rescaled <- min(center_point + expanded_range/2, 1.0)  # Don't exceed 1.0
      } else {
        # For larger ranges, use a smaller expansion
        buffer <- range_span * 0.2  # Add 20% buffer on each side
        cor_min_rescaled <- max(cor_range[1] - buffer, 0.3)
        cor_max_rescaled <- min(cor_range[2] + buffer, 1.0)
      }
      
      message("Calculated rescaling range: ", sprintf("%.3f to %.3f", cor_min_rescaled, cor_max_rescaled))
      
      # Linear rescaling to map actual range to the rescaled range
      display_matrix <- cor_min_rescaled + (cor_max_rescaled - cor_min_rescaled) * 
        (cor_matrix - cor_range[1]) / (cor_range[2] - cor_range[1])
      
      # Keep diagonal at max value for proper correlation interpretation
      diag(display_matrix) <- cor_max_rescaled
      
      # Use the selected color scheme
      col_scheme <- colors
      
      message("Original correlation range: ", sprintf("%.4f to %.4f", cor_range[1], cor_range[2]))
      message("Rescaled display range: ", sprintf("%.3f to %.3f", cor_min_rescaled, cor_max_rescaled))
      message("Using color scheme: ", color_scheme)
      
      # Set up corrplot with rescaled data
      corrplot::corrplot(
        display_matrix,
        method = "color",
        type = "full",
        order = "hclust",
        tl.cex = 0.8,
        tl.col = "black",
        tl.srt = 45,
        col = col_scheme,
        title = "Sample Correlation Heatmap",
        mar = c(0,0,2,0),
        addCoef.col = if(show_values) "black" else NULL,
        number.cex = 0.7,
        is.corr = FALSE  # Set to FALSE since we're using rescaled values
      )
      
      # Add a note about the actual correlation values using base graphics
      mtext("Note: Colors enhanced for high correlation visualization", 
            side = 1, line = 3, cex = 0.8, col = "gray50")
      mtext(sprintf("Actual range: %.4f to %.4f | Display range: %.3f to %.3f", 
            cor_range[1], cor_range[2], cor_min_rescaled, cor_max_rescaled), 
            side = 1, line = 4, cex = 0.8, col = "gray50")
      
    } else {
      # Standard range - use original matrix
      message("Using standard correlation visualization")
      message("Condition NOT met for high correlation scaling")
      display_matrix <- cor_matrix
      col_scheme <- colors
      message("Using color scheme: ", color_scheme)
      message("Using standard color limits: -1 to 1")
      
      # Set up corrplot parameters with standard scaling
      corrplot::corrplot(
        display_matrix,
        method = "color",
        type = "full",
        order = "hclust",
        tl.cex = 0.8,
        tl.col = "black",
        tl.srt = 45,
        col = col_scheme,
        title = "Sample Correlation Heatmap",
        mar = c(0,0,2,0),
        addCoef.col = if(show_values) "black" else NULL,
        number.cex = 0.7,
        is.corr = TRUE  # Set to TRUE for standard correlation scaling
      )
    }
    
    message("Corrplot created successfully")
    return("corrplot_success")
    
  }, error = function(e) {
    message("Error in corrplot rendering: ", e$message)
    message("Error class: ", class(e))
    message("Error call: ", deparse(e$call))
    
    # Simple error plot
    plot(NULL, xlim=c(0,1), ylim=c(0,1), ann=FALSE, axes=FALSE)
    text(0.5, 0.5, paste("Error rendering correlation heatmap:", e$message), 
         cex=1.2, col="red")
    text(0.5, 0.3, "Please check your data and try again", 
         cex=1, col="black")
    return(NULL)
  })
  
  # Return the result
  return(result)
}

# --- 3. Normality Assessment Helper Functions ---

#' Assess Normality of Samples
#' @param counts_data Count matrix
#' @param samples_to_test Samples to test (defaults to all samples)
#' @param test_method Normality test method ("shapiro", "ks", or "ad")
#' @param significance_level Significance level for tests (default 0.05)
#' @return List with normality assessment
#' @export
assess_normality_func <- function(counts_data, samples_to_test = NULL, 
                                test_method = "shapiro", significance_level = 0.05) {
  # Initialize results list
  sample_stats <- list()
  
  # Use all samples if none specified
  if (is.null(samples_to_test)) {
    samples_to_test <- colnames(counts_data)
  }
  
  # Assess each selected sample
  for (sample in samples_to_test) {
    if (!(sample %in% colnames(counts_data))) {
      next  # Skip samples not in the data
    }
    
    sample_data <- counts_data[, sample]
    
    # Calculate basic statistics
    skew <- moments::skewness(sample_data)
    kurt <- moments::kurtosis(sample_data)
    
    # Perform normality test based on method
    p_value <- NA
    test_name <- "Unknown"
    
    if (test_method == "shapiro") {
      test_name <- "Shapiro-Wilk"
      # Shapiro-Wilk has sample size limitations
      if (length(sample_data) < 5000) {
        p_value <- shapiro.test(sample_data)$p.value
      }
    } else if (test_method == "ks") {
      test_name <- "Kolmogorov-Smirnov"
      # KS test against normal distribution
      p_value <- ks.test(scale(sample_data), "pnorm")$p.value
    } else if (test_method == "ad") {
      test_name <- "Anderson-Darling"
      # AD test requires nortest package
      if (requireNamespace("nortest", quietly = TRUE)) {
        p_value <- nortest::ad.test(sample_data)$p.value
      }
    }
    
    # Determine deviation severity based on thresholds
    skew_deviation <- "Normal range"
    if (abs(skew) > 3) {
      skew_deviation <- "Severe deviation"
    } else if (abs(skew) > 2) {
      skew_deviation <- "Moderate deviation"
    }
    
    kurt_deviation <- "Normal range"
    if (abs(kurt - 3) > 10) {
      kurt_deviation <- "Severe deviation"
    } else if (abs(kurt - 3) > 7) {
      kurt_deviation <- "Moderate deviation"
    }
    
    p_value_deviation <- "Normal range"
    if (!is.na(p_value)) {
      if (p_value < 0.001) {
        p_value_deviation <- "Severe deviation"
      } else if (p_value < 0.01) {
        p_value_deviation <- "Moderate deviation"
      }
    }
    
    # Determine overall normality status
    has_severe_deviation <- (abs(skew) > 3) || (abs(kurt - 3) > 10) || (!is.na(p_value) && p_value < 0.001)
    has_moderate_deviation <- (abs(skew) > 2) || (abs(kurt - 3) > 7) || (!is.na(p_value) && p_value < 0.01)
    
    # Store results
    sample_stats[[sample]] <- list(
      skewness = skew,
      skewness_status = skew_deviation,
      kurtosis = kurt,
      kurtosis_status = kurt_deviation,
      test_name = test_name,
      p_value = p_value,
      p_value_status = p_value_deviation,
      is_normal = !is.na(p_value) && p_value > significance_level,
      has_severe_deviation = has_severe_deviation,
      has_moderate_deviation = has_moderate_deviation
    )
  }
  
  # Count samples with different deviation levels
  total_samples <- length(sample_stats)
  normal_count <- sum(sapply(sample_stats, function(x) !x$has_moderate_deviation && !x$has_severe_deviation))
  moderate_count <- sum(sapply(sample_stats, function(x) x$has_moderate_deviation && !x$has_severe_deviation))
  severe_count <- sum(sapply(sample_stats, function(x) x$has_severe_deviation))
  
  # Calculate percentages
  normal_pct <- round(100 * normal_count / total_samples, 1)
  moderate_pct <- round(100 * moderate_count / total_samples, 1)
  severe_pct <- round(100 * severe_count / total_samples, 1)
  
  # Generate detailed recommendation
  recommendation <- NULL
  if (severe_pct > 20) {
    recommendation <- list(
      summary = sprintf("STRONGLY RECOMMEND using Spearman correlation"),
      reasoning = sprintf("%.1f%% of samples (%d of %d) show severe deviations from normality", 
                         severe_pct, severe_count, total_samples),
      correlation_method = "spearman"
    )
  } else if (moderate_pct + severe_pct > 30) {
    recommendation <- list(
      summary = sprintf("RECOMMEND using Spearman correlation"),
      reasoning = sprintf("%.1f%% of samples (%d of %d) show moderate or severe deviations from normality", 
                         moderate_pct + severe_pct, moderate_count + severe_count, total_samples),
      correlation_method = "spearman"
    )
  } else {
    recommendation <- list(
      summary = sprintf("RECOMMEND using Pearson correlation"),
      reasoning = sprintf("%.1f%% of samples (%d of %d) appear approximately normally distributed", 
                         normal_pct, normal_count, total_samples),
      correlation_method = "pearson"
    )
  }
  
  # Create summary statistics
  summary_stats <- list(
    total_samples = total_samples,
    normal_count = normal_count,
    normal_pct = normal_pct,
    moderate_count = moderate_count,
    moderate_pct = moderate_pct,
    severe_count = severe_count,
    severe_pct = severe_pct,
    test_method = test_method,
    significance_level = significance_level
  )
  
  list(
    sample_stats = sample_stats,
    summary_stats = summary_stats,
    recommendation = recommendation
  )
}

#' Generate Q-Q Plots for Normality Assessment
#' @param counts_data Count matrix
#' @param samples_to_plot Samples to plot
#' @return A ggplot object
#' @export
plot_normality_qq_func <- function(counts_data, samples_to_plot = NULL) {
  # Select samples to plot
  if (is.null(samples_to_plot)) {
    samples_to_plot <- colnames(counts_data)[1:min(9, ncol(counts_data))]
  }
  
  # Limit to a reasonable number of samples for visualization
  if (length(samples_to_plot) > 9) {
    samples_to_plot <- samples_to_plot[1:9]
  }
  
  # Create plot data
  plot_data <- lapply(samples_to_plot, function(sample) {
    theoretical_q <- qqnorm(counts_data[, sample], plot.it = FALSE)
    data.frame(
      Sample = sample,
      Theoretical = theoretical_q$x,
      Observed = theoretical_q$y
    )
  })
  plot_df <- do.call(rbind, plot_data)
  
  # Create plot with improved styling
  ggplot(plot_df, aes(x = Theoretical, y = Observed)) +
    geom_point(alpha = 0.5, color = "steelblue") +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    facet_wrap(~Sample, scales = "free") +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "lightblue", color = "black"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    labs(title = "Normal Q-Q Plots",
         subtitle = "Points should follow the red line if normally distributed",
         x = "Theoretical Quantiles",
         y = "Sample Quantiles")
}

#' Generate Density Plots for Normality Assessment
#' @param counts_data Count matrix
#' @param samples_to_plot Samples to plot
#' @return A ggplot object
#' @export
plot_normality_density_func <- function(counts_data, samples_to_plot = NULL) {
  # Select samples to plot
  if (is.null(samples_to_plot)) {
    samples_to_plot <- colnames(counts_data)[1:min(9, ncol(counts_data))]
  }
  
  # Limit to a reasonable number of samples for visualization
  if (length(samples_to_plot) > 9) {
    samples_to_plot <- samples_to_plot[1:9]
  }
  
  # Create plot data with normal distribution overlay
  plot_data <- list()
  normal_curves <- list()
  
  for (sample in samples_to_plot) {
    sample_data <- counts_data[, sample]
    d <- density(sample_data)
    
    # Add density data
    plot_data[[length(plot_data) + 1]] <- data.frame(
      Sample = sample,
      x = d$x,
      y = d$y,
      Type = "Actual"
    )
    
    # Calculate normal distribution parameters
    mu <- mean(sample_data)
    sigma <- sd(sample_data)
    x_range <- seq(min(d$x), max(d$x), length.out = 100)
    y_normal <- dnorm(x_range, mean = mu, sd = sigma)
    
    # Add normal curve data
    normal_curves[[length(normal_curves) + 1]] <- data.frame(
      Sample = sample,
      x = x_range,
      y = y_normal,
      Type = "Normal"
    )
  }
  
  # Combine data
  plot_df <- do.call(rbind, c(plot_data, normal_curves))
  
  # Create plot with both density and normal distribution
  ggplot(plot_df, aes(x = x, y = y, color = Type, linetype = Type)) +
    geom_line(size = 1) +
    scale_color_manual(values = c("Actual" = "steelblue", "Normal" = "red")) +
    scale_linetype_manual(values = c("Actual" = "solid", "Normal" = "dashed")) +
    facet_wrap(~Sample, scales = "free") +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "lightblue", color = "black"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    ) +
    labs(title = "Density Plots with Normal Distribution Overlay",
         subtitle = "Blue line: actual distribution, Red dashed: theoretical normal",
         x = "Expression Value",
         y = "Density")
}

#' Generate Combined Histogram and Density Plots for Normality Assessment
#' @param counts_data Count matrix
#' @param samples_to_plot Samples to plot
#' @return A ggplot object
#' @export
plot_normality_hist_density_func <- function(counts_data, samples_to_plot = NULL) {
  # Select samples to plot
  if (is.null(samples_to_plot)) {
    samples_to_plot <- colnames(counts_data)[1:min(6, ncol(counts_data))]
  }
  
  # Limit to a reasonable number of samples for visualization
  if (length(samples_to_plot) > 6) {
    samples_to_plot <- samples_to_plot[1:6]
  }
  
  # Create plot data with additional info for normal curves
  all_data <- data.frame()
  normal_curves <- data.frame()
  
  for (sample in samples_to_plot) {
    sample_data <- counts_data[, sample]
    
    # Add histogram data
    all_data <- rbind(all_data, data.frame(
      Sample = sample,
      Value = sample_data
    ))
    
    # Calculate normal distribution parameters for this sample
    mu <- mean(sample_data)
    sigma <- sd(sample_data)
    x_range <- seq(min(sample_data), max(sample_data), length.out = 100)
    y_normal <- dnorm(x_range, mean = mu, sd = sigma)
    
    # Add normal curve data
    normal_curves <- rbind(normal_curves, data.frame(
      Sample = sample,
      x = x_range,
      y = y_normal
    ))
  }
  
  # Create plot with histogram and density overlay
  p <- ggplot(all_data, aes(x = Value)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "darkblue", alpha = 0.7) +
    geom_density(color = "red", size = 1) +
    geom_line(data = normal_curves, aes(x = x, y = y), color = "darkgreen", size = 1, linetype = "dashed") +
    facet_wrap(~Sample, scales = "free") +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "lightblue", color = "black"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    ) +
    labs(title = "Histogram with Density Overlay",
         subtitle = "Red line: actual density, Green dashed: theoretical normal",
         x = "Expression Value",
         y = "Density")
  
  return(p)
}

#' Save Heatmap Plot to File
#' @param cor_results Correlation results
#' @param file Output file
#' @param metadata_df Metadata
#' @param annotation_cols Annotation columns
#' @param color_scheme Color scheme
#' @export
save_heatmap_plot_func <- function(cor_results, file, metadata_df, annotation_cols, color_scheme) {
  pdf(file, width = 10, height = 8)
  plot_correlation_heatmap_func(cor_results, metadata_df, annotation_cols, color_scheme)
  dev.off()
}

#' Generate QC Report
#' @param raw_counts Raw counts
#' @param processed_data Processed data
#' @param metadata Metadata
#' @param normality_results Normality results
#' @param pca_results PCA results
#' @param output_file Output file
#' @export
generate_qc_report_func <- function(raw_counts, processed_data, metadata, 
                                  normality_results, pca_results, output_file) {
  # This function would generate an HTML report using rmarkdown
  # Implementation would depend on specific reporting requirements
  # For now, we'll create a basic report template
  rmarkdown::render(
    system.file("rmd", "qc_report_template.Rmd", package = "RNAProcessing"),
    output_file = output_file,
    params = list(
      raw_counts = raw_counts,
      processed_data = processed_data,
      metadata = metadata,
      normality_results = normality_results,
      pca_results = pca_results
    )
  )
} 