# server_tab2_qc_plots.R
# QC & Pre-processing Tool
# Tab 2: Quality Control Plots & Summaries - Server Logic
# Author: Eren Ada, PhD
# Date: 05/13/2024

#' @import shiny
#' @import ggplot2
#' @import plotly
#' @import DESeq2
#' @import edgeR
#' @import dplyr
#' @import grid
#' @importFrom stats cor prcomp
#' @importFrom S4Vectors DataFrame
NULL

#' Server Module for Quality Control Plots & Summaries tab
#' @param id Module ID
#' @param shared_data A reactive list containing processed data and metadata
#' @return A shiny server module
#' @export
mod_qc_plots_server <- function(id, shared_data) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive values for plot data management
    plot_data <- reactiveVal(NULL)
    current_data_description <- reactiveVal("No data processed yet.")
    raw_counts_for_tab <- reactiveVal(NULL)
    
    # Update plot data when shared data changes
    observe({
      if (isTRUE(getOption("app.debug"))) {
        message("\n=== Data Transformation Debug ===")
        message("Checking shared_data availability:")
        message("- data_processed flag: ", !is.null(shared_data$data_processed))
        message("- processed_counts available: ", !is.null(shared_data$processed_counts))
        message("- processed_metadata available: ", !is.null(shared_data$processed_metadata))
      }
      
      req(shared_data$data_processed)
      req(shared_data$processed_counts)
      req(shared_data$processed_metadata)
      
      counts_matrix <- shared_data$processed_counts
      if (isTRUE(getOption("app.debug"))) {
        message("\nCounts Matrix Info:")
        message("- Dimensions: ", paste(dim(counts_matrix), collapse=" x "))
        message("- Class: ", class(counts_matrix))
        message("- Sample names: ", paste(head(colnames(counts_matrix)), collapse=", "))
      }
      
      raw_counts_for_tab(counts_matrix)
      
      if (!is.null(shared_data$normalized_counts) && 
          !is.null(shared_data$normalization_status_flag) && 
          shared_data$normalization_status_flag) {
        message("\nUsing normalized counts")
        plot_data(shared_data$normalized_counts)
        current_data_description(
          paste("Displaying plots for:", 
                shared_data$normalization_method_used, 
                "normalized data")
        )
      } else {
        message("\nPerforming log2(CPM + 1) transformation")
        # Default: Use log2(CPM + 1) transformation
        if (!is.null(counts_matrix) && all(dim(counts_matrix) > 0)) {
          message("- Converting to CPM")
          cpm_counts <- edgeR::cpm(counts_matrix, log = FALSE)
          message("- Applying log2 transformation")
          log2_cpm_plus_1 <- log2(cpm_counts + 1)
          plot_data(log2_cpm_plus_1)
          message("- Transformation complete")
          message("- Transformed data dimensions: ", paste(dim(log2_cpm_plus_1), collapse=" x "))
          current_data_description(
            "Displaying plots for: log2(CPM+1) transformed data (pre-normalization)"
          )
        } else {
          plot_data(NULL)
          current_data_description("Waiting for valid count data...")
        }
      }
    })
    
    # Data Status Text
    output$data_status_text <- renderUI({
      tags$p(
        tags$em(current_data_description()),
        style = "text-align: center; margin-bottom: 15px;"
      )
    })
    
    # Library Size Distribution Plot
    output$library_size_plot <- renderPlotly({
      req(raw_counts_for_tab(), 
          plot_data(), 
          input$lib_size_plot_type, 
          !is.null(input$lib_size_log_scale))
      
      plot_library_sizes_func(
        raw_counts_data = raw_counts_for_tab(),
        current_view_data = plot_data(),
        plot_type = input$lib_size_plot_type,
        use_log_scale = input$lib_size_log_scale
      )
    })
    
    # Gene Detection Rates Plot
    output$gene_detection_plot <- renderPlotly({
      req(plot_data(), !is.null(input$gene_detection_threshold))
      
      plot_gene_detection_func(
        counts_data = plot_data(),
        detection_threshold = input$gene_detection_threshold
      )
    })
    
    # Summary Statistics
    output$summary_stats_table <- DT::renderDataTable({
      req(raw_counts_for_tab(), plot_data())
      
      # Calculate statistics
      stats_list <- calculate_qc_summary_stats_func(
        raw_counts_data = raw_counts_for_tab(),
        current_view_data = plot_data(),
        detection_threshold = input$gene_detection_threshold
      )
      
      # Create a data frame for display
      data.frame(
        Metric = c(
          "Raw Library Size (Min)", "Raw Library Size (Max)",
          "Raw Library Size (Mean)", "Raw Library Size (Median)",
          "Processed Library Size (Min)", "Processed Library Size (Max)",
          "Processed Library Size (Mean)", "Processed Library Size (Median)",
          "Gene Detection Rate (Min %)", "Gene Detection Rate (Max %)",
          "Gene Detection Rate (Mean %)", "Gene Detection Rate (Median %)"
        ),
        Value = c(
          format(stats_list$raw_library_sizes$min, scientific = FALSE, big.mark = ","),
          format(stats_list$raw_library_sizes$max, scientific = FALSE, big.mark = ","),
          format(stats_list$raw_library_sizes$mean, scientific = FALSE, big.mark = ","),
          format(stats_list$raw_library_sizes$median, scientific = FALSE, big.mark = ","),
          format(stats_list$current_library_sizes$min, scientific = FALSE, big.mark = ","),
          format(stats_list$current_library_sizes$max, scientific = FALSE, big.mark = ","),
          format(stats_list$current_library_sizes$mean, scientific = FALSE, big.mark = ","),
          format(stats_list$current_library_sizes$median, scientific = FALSE, big.mark = ","),
          sprintf("%.2f", stats_list$gene_detection$min_rate),
          sprintf("%.2f", stats_list$gene_detection$max_rate),
          sprintf("%.2f", stats_list$gene_detection$mean_rate),
          sprintf("%.2f", stats_list$gene_detection$median_rate)
        )
      )
    }, options = list(
      pageLength = 12,
      dom = 't',
      ordering = FALSE
    ))
    
    # Initialize PCA submodule using current plot_data
    pca_module <- mod_pca_server("pca", shared_data, plot_data)
    
    # PCA results storage provided by submodule
    
    # Correlation annotations UI moved to correlation submodule
    
    # Reactive value to store parameters for the correlation plot
    correlation_plot_params <- reactiveVal(NULL)
    
    # Observer to set initial correlation plot parameters
    observe({
      message("\n=== Initial Correlation Parameters Setup ===")
      
      # Require necessary data
      req(plot_data(), shared_data$processed_metadata, input$correlation_method)
      
      current_plot_data <- plot_data()
      current_metadata <- shared_data$processed_metadata
      
      message("Plot data dimensions: ", paste(dim(current_plot_data), collapse=" x "))
      message("Metadata dimensions: ", paste(dim(current_metadata), collapse=" x "))
      
      # Validate data
      if (is.null(current_plot_data) || !is.matrix(current_plot_data) || ncol(current_plot_data) < 2) {
        message("ERROR: Invalid plot data")
        return()
      }
      
      if (is.null(current_metadata) || !is.data.frame(current_metadata)) {
        message("ERROR: Invalid metadata")
        return()
      }
      
      # Calculate correlation
      message("Calculating initial correlation...")
      cor_res <- tryCatch({
        calculate_correlation_func(
          data_matrix = current_plot_data,
          method = input$correlation_method
        )
      }, error = function(e) {
        message("Error in initial correlation calculation: ", e$message)
        return(NULL)
      })
      
      if (is.null(cor_res) || is.null(cor_res$correlation)) {
        message("ERROR: Correlation calculation failed")
        return()
      }
      
      message("Initial correlation calculation successful")
      message("Correlation matrix dimensions: ", paste(dim(cor_res$correlation), collapse=" x "))
      
      # Match samples between correlation and metadata
      samples_in_cor <- rownames(cor_res$correlation)
      samples_in_meta <- rownames(current_metadata)
      
      message("Samples in correlation: ", length(samples_in_cor))
      message("Samples in metadata: ", length(samples_in_meta))
      
      # Find common samples
      common_samples <- intersect(samples_in_cor, samples_in_meta)
      message("Common samples: ", length(common_samples))
      
      # Create metadata subset for annotations
      metadata_for_annot <- NULL
      if (length(common_samples) > 0) {
        metadata_for_annot <- current_metadata[common_samples, , drop = FALSE]
        message("Metadata subset created with dimensions: ", paste(dim(metadata_for_annot), collapse=" x "))
      } else {
        message("WARNING: No matching samples found in metadata for annotation")
      }
      
      # Set initial correlation plot parameters
      correlation_plot_params(list(
        cor_results = cor_res,
        metadata_for_annot = metadata_for_annot
      ))
      
      # Set initial correlation results for statistics
      correlation_results(cor_res)
      
      message("Initial correlation plot parameters set successfully")
    })
    
    # Dynamic UI for metadata-based filters
    output$filter_criteria <- renderUI({
      message("\n=== Rendering Filter Criteria UI ===")
      req(shared_data$processed_metadata)
      
      metadata_cols <- colnames(shared_data$processed_metadata)
      message("Available filter columns: ", paste(metadata_cols, collapse=", "))
      
      checkboxGroupInput(
        session$ns("filter_criteria"),
        "Select Filters to Apply:",
        choices = metadata_cols,
        selected = NULL
      )
    })
    
    # Dynamic filter UI generator
    output$dynamic_filters <- renderUI({
      message("\n=== Rendering Dynamic Filters UI ===")
      req(shared_data$processed_metadata, input$filter_criteria)
      
      message("Selected filter criteria: ", paste(input$filter_criteria, collapse=", "))
      
      filter_list <- lapply(input$filter_criteria, function(col) {
        unique_values <- unique(shared_data$processed_metadata[[col]])
        message("Column ", col, " has ", length(unique_values), " unique values")
        
        selectInput(
          session$ns(paste0("filter_", col)),
          paste("Select", col, ":"),
          choices = unique_values,
          selected = unique_values,  # Select all values by default
          multiple = TRUE
        )
      })
      
      do.call(tagList, filter_list)
    })
    
    # Filtered sample selection with improved validation
    selected_samples <- reactive({
      message("\n=== Computing Selected Samples ===")
      req(shared_data$processed_metadata)
      
      # Start with all samples
      all_samples <- rownames(shared_data$processed_metadata)
      message("Total available samples: ", length(all_samples))
      
      # If no filters selected, return all samples
      if (is.null(input$filter_criteria) || length(input$filter_criteria) == 0) {
        message("No filters selected, returning all samples")
        return(all_samples)
      }
      
      # Apply each filter
      selected <- all_samples
      for (col in input$filter_criteria) {
        filter_input_name <- paste0("filter_", col)
        
        # Check if the filter input exists and has values
        if (!is.null(input[[filter_input_name]]) && length(input[[filter_input_name]]) > 0) {
          message("Applying filter for ", col, ": ", paste(input[[filter_input_name]], collapse=", "))
          
          # Get samples that match the filter
          matching_samples <- rownames(shared_data$processed_metadata)[shared_data$processed_metadata[[col]] %in% input[[filter_input_name]]]
          message("Samples matching this filter: ", length(matching_samples))
          
          # Update selected samples
          selected <- intersect(selected, matching_samples)
          message("Samples remaining after this filter: ", length(selected))
        } else {
          message("Filter ", col, " has no values selected, skipping")
        }
      }
      
      # Ensure we have valid sample names
      if (length(selected) == 0) {
        message("WARNING: No samples match all filters, returning all samples")
        return(all_samples)
      }
      
      message("Final selected samples: ", length(selected))
      return(selected)
    })
    
    # Correlation Statistics Summary
    output$correlation_stats_text <- renderUI({
      req(correlation_results(), selected_samples()) # This uses the correlation_results reactiveVal
      cor_res <- correlation_results()
      data_desc <- correlation_data_description()

      # Add a check for valid correlation results before trying to access them
      if (is.null(cor_res) || is.null(cor_res$correlation) || 
          nrow(cor_res$correlation) == 0 || ncol(cor_res$correlation) == 0 ||
          (nrow(cor_res$correlation) == 1 && ncol(cor_res$correlation) == 1 && is.na(cor_res$correlation[1,1])) # check for single NA matrix
         ) {
        return(tags$div(
          style = "margin-top: 15px; padding: 10px; background-color: #f8f9fa; border-radius: 5px;",
          tags$h4("Correlation Summary Statistics:", style = "margin-bottom: 15px;"),
          tags$p("Not enough samples selected or data available to calculate correlation statistics.")
        ))
      }
      
      # Calculate summary statistics
      cor_values <- cor_res$correlation[upper.tri(cor_res$correlation)]
      # Ensure cor_values is not empty before calculating stats
      if(length(cor_values) == 0) {
        mean_cor <- NA
        median_cor <- NA
        min_cor <- NA
        max_cor <- NA
      } else {
        mean_cor <- mean(cor_values, na.rm = TRUE)
        median_cor <- median(cor_values, na.rm = TRUE)
        min_cor <- min(cor_values, na.rm = TRUE)
        max_cor <- max(cor_values, na.rm = TRUE)
      }
      
      # Calculate significant correlations
      # Ensure p_adjusted is not NULL and has same dimensions
      if (!is.null(cor_res$p_adjusted) && all(dim(cor_res$p_adjusted) == dim(cor_res$correlation))) {
        sig_cors <- sum(cor_res$p_adjusted < 0.05 & upper.tri(cor_res$correlation), na.rm = TRUE)
      } else {
        sig_cors <- NA
      }
      total_cors <- sum(upper.tri(cor_res$correlation))
      
      # Get filter information
      filter_info <- if (length(input$filter_criteria) > 0) {
        paste(sapply(input$filter_criteria, function(col) {
          selected <- input[[paste0("filter_", col)]]
          sprintf("%s: %s", col, paste(selected, collapse = ", "))
        }), collapse = "<br>")
      } else {
        "No filters applied"
      }
      
      tagList(
        tags$div(
          style = "margin-top: 15px; padding: 10px; background-color: #f8f9fa; border-radius: 5px;",
          tags$h4("Correlation Summary Statistics:", 
                 style = "margin-bottom: 15px;"),
          tags$div(
            style = "margin-bottom: 10px; padding: 8px; background-color: #e3f2fd; border-radius: 4px;",
            tags$strong("Data Source: "), tags$span(data_desc, style = "font-weight: normal;")
          ),
          tags$div(
            style = "margin-bottom: 15px;",
            tags$strong("Applied Filters:"),
            tags$p(HTML(filter_info))
          ),
          fluidRow(
            column(width = 6,
              tags$ul(
                style = "list-style-type: none; padding-left: 0;",
                tags$li(sprintf("Mean correlation: %.3f", mean_cor)),
                tags$li(sprintf("Median correlation: %.3f", median_cor)),
                tags$li(sprintf("Range: %.3f to %.3f", min_cor, max_cor))
              )
            ),
            column(width = 6,
              tags$ul(
                style = "list-style-type: none; padding-left: 0;",
                tags$li(sprintf("Significant correlations: %d", sig_cors)),
                tags$li(sprintf("Total correlations: %d", total_cors)),
                tags$li(sprintf("Significance rate: %.1f%%", 100 * sig_cors/total_cors))
              )
            )
          ),
          tags$div(
            style = "margin-top: 15px;",
            tags$p(sprintf("Number of samples included: %d", length(selected_samples())))
          )
        )
      )
    })
    
    # Removed refresh handler; use explicit update actions per section
    
    # Removed data source guidance; only processed/log2(CPM+1) or normalized data are used
    
    # Observer to handle correlation plot updates
    observeEvent(input$update_correlation, {
      message("\n=== Update Correlation Button Clicked ===")
      
      # Validate required inputs
      req(shared_data$processed_metadata)
      
      # Determine data source: use normalized counts if present; otherwise use current plot data (log2(CPM+1))
      correlation_data <- NULL
      data_description <- ""
      if (!is.null(shared_data$normalized_counts) && 
          !is.null(shared_data$normalization_status_flag) && 
          shared_data$normalization_status_flag) {
        correlation_data <- shared_data$normalized_counts
        data_description <- paste(shared_data$normalization_method_used, "normalized data")
        if (isTRUE(getOption("app.debug"))) message("Using normalized data for correlation analysis")
      } else {
        req(plot_data())
        correlation_data <- plot_data()
        data_description <- "log2(CPM+1) transformed data"
        if (isTRUE(getOption("app.debug"))) message("Using log2(CPM+1) transformed data for correlation analysis")
      }
      
      # Validate correlation data
      if (is.null(correlation_data) || all(dim(correlation_data) == 0)) {
        showNotification("No valid data available for correlation analysis", 
                        type = "error", duration = 5)
        return()
      }
      
      # Get current metadata
      current_metadata <- shared_data$processed_metadata
      
      message("Update button - Correlation data dimensions: ", paste(dim(correlation_data), collapse=" x "))
      message("Update button - Data description: ", data_description)
      message("Update button - Metadata dimensions: ", paste(dim(current_metadata), collapse=" x "))
      
      message("Starting correlation update process...")
      
      # Apply filters to get the correlation data
      cor_data_filtered <- tryCatch({
        # Identify active filter criteria
        active_filters <- input$filter_criteria
        message("=== FILTER APPLICATION DEBUG ===")
        message("Active filter criteria: ", paste(active_filters, collapse=", "))
        
        # Apply each selected filter
        filtered_data <- correlation_data
        filtered_metadata <- current_metadata
        
        message("Starting data dimensions:")
        message("- Data: ", paste(dim(filtered_data), collapse=" x "))
        message("- Metadata: ", paste(dim(filtered_metadata), collapse=" x "))
        
        # Apply filters based on selected criteria
        for (filter_type in active_filters) {
          if (filter_type == "Sample") {
            selected_samples <- input[[paste0("filter_", filter_type)]]
            if (!is.null(selected_samples) && length(selected_samples) > 0) {
              message("Applying Sample filter: ", paste(selected_samples, collapse=", "))
              # Filter data to selected samples
              common_samples <- intersect(colnames(filtered_data), selected_samples)
              message("Found ", length(common_samples), " common samples out of ", length(selected_samples), " selected")
              if (length(common_samples) > 0) {
                filtered_data <- filtered_data[, common_samples, drop = FALSE]
                # Also filter metadata to match
                if (any(rownames(filtered_metadata) %in% common_samples)) {
                  filtered_metadata <- filtered_metadata[rownames(filtered_metadata) %in% common_samples, , drop = FALSE]
                }
                message("After Sample filter - Data: ", paste(dim(filtered_data), collapse=" x "))
              }
            }
          } else {
            # Handle other filter types (colonization, time, cell_type, group)
            selected_values <- input[[paste0("filter_", filter_type)]]
            if (!is.null(selected_values) && length(selected_values) > 0) {
              message("Applying ", filter_type, " filter: ", paste(selected_values, collapse=", "))
              # Find samples that match the filter criteria
              if (filter_type %in% colnames(filtered_metadata)) {
                matching_samples <- rownames(filtered_metadata)[filtered_metadata[[filter_type]] %in% selected_values]
                message("Found ", length(matching_samples), " samples matching ", filter_type, " criteria")
                common_samples <- intersect(colnames(filtered_data), matching_samples)
                message("Found ", length(common_samples), " common samples in data")
                if (length(common_samples) > 0) {
                  filtered_data <- filtered_data[, common_samples, drop = FALSE]
                  filtered_metadata <- filtered_metadata[rownames(filtered_metadata) %in% common_samples, , drop = FALSE]
                  message("After ", filter_type, " filter - Data: ", paste(dim(filtered_data), collapse=" x "))
                }
              } else {
                message("Warning: Filter column '", filter_type, "' not found in metadata")
              }
            }
          }
        }
        
        message("=== FINAL FILTERING RESULTS ===")
        message("Final filtered data dimensions: ", paste(dim(filtered_data), collapse=" x "))
        message("Final metadata dimensions: ", paste(dim(filtered_metadata), collapse=" x "))
        message("Final sample names: ", paste(colnames(filtered_data), collapse=", "))
        
        # Ensure we have enough samples for correlation
        if (ncol(filtered_data) < 2) {
          message("ERROR: Not enough samples for correlation (need at least 2)")
          showNotification("Not enough samples selected for correlation analysis (minimum 2 required)", 
                          type = "error", duration = 5)
          return(NULL)
        }
        
        # Calculate correlation
        message("=== CORRELATION CALCULATION ===")
        message("Calculating correlation using method: ", input$correlation_method)
        message("Input data for correlation: ", paste(dim(filtered_data), collapse=" x "))
        
        cor_results <- calculate_correlation_func(filtered_data, input$correlation_method, compute_p = isTRUE(input$compute_p_values))
        
        if (!is.null(cor_results) && !is.null(cor_results$correlation)) {
          message("Correlation calculation successful")
          message("Correlation matrix dimensions: ", paste(dim(cor_results$correlation), collapse=" x "))
          
          # Calculate and report correlation range for scaling debug
          cor_values <- cor_results$correlation[upper.tri(cor_results$correlation)]
          if (length(cor_values) > 0) {
            cor_range <- range(cor_values, na.rm = TRUE)
            message("=== CORRELATION RANGE FOR SCALING ===")
            message("Correlation range: ", sprintf("%.4f to %.4f", cor_range[1], cor_range[2]))
            message("Range span: ", sprintf("%.4f", cor_range[2] - cor_range[1]))
            message("Min correlation > 0.8: ", cor_range[1] > 0.8)
            message("Min correlation > 0.7 AND range < 0.3: ", 
                   cor_range[1] > 0.7 && (cor_range[2] - cor_range[1]) < 0.3)
            message("Will trigger high-correlation scaling: ", 
                   (cor_range[1] > 0.8 || (cor_range[1] > 0.7 && (cor_range[2] - cor_range[1]) < 0.3)))
          }
        } else {
          message("ERROR: Correlation calculation failed")
          return(NULL)
        }
        
        cor_results
      }, error = function(e) {
        message("Error in correlation calculation: ", e$message)
        message("Error traceback: ", paste(traceback(), collapse="\n"))
        showNotification(paste("Error calculating correlation:", e$message), 
                        type = "error", duration = 5)
        return(NULL)
      })
      
      # Validate correlation results
      if (is.null(cor_data_filtered)) {
        message("Correlation calculation failed - stopping update")
        return()
      }
      
      # Prepare metadata for annotations (matching correlation matrix samples)
      final_metadata_for_plot <- NULL
      metadata_for_plot <- shared_data$processed_metadata
      
      if (!is.null(metadata_for_plot) && nrow(metadata_for_plot) > 0) {
        # Match samples between correlation and metadata
        cor_samples <- rownames(cor_data_filtered$correlation)
        meta_samples <- rownames(metadata_for_plot)
        common_samples <- intersect(cor_samples, meta_samples)
        
        if (length(common_samples) > 0) {
          # Reorder metadata to match correlation matrix
          final_metadata_for_plot <- metadata_for_plot[common_samples, , drop = FALSE]
          message("Final metadata for annotations: ", nrow(final_metadata_for_plot), " samples")
        } else {
          message("No matching samples between correlation and metadata")
        }
      }
      
      # Update the reactive correlation parameters
      message("Updating correlation_plot_params reactive value...")
      correlation_plot_params(list(
        cor_results = cor_data_filtered,
        metadata_for_annot = final_metadata_for_plot,
        filter_timestamp = Sys.time(),  # Add timestamp to ensure reactivity
        filter_info = list(
          active_filters = active_filters,
          sample_count = ncol(correlation_data),
          filtered_sample_count = if(!is.null(cor_data_filtered)) ncol(cor_data_filtered$correlation) else 0
        )
      ))
      message("correlation_plot_params updated successfully")
      
      # Update correlation results for stats display
      correlation_results(cor_data_filtered)
      correlation_data_description(data_description)
      
      # Show success notification
      showNotification("Correlation plot updated successfully!", 
                      type = "message", duration = 3)
      
      message("Update correlation process completed")
    })
    
    # Correlation Heatmap Plot - Single rendering function
    output$correlation_heatmap <- renderPlot({
      # Add explicit dependency on update button to ensure re-rendering
      input$update_correlation
      
      message("\n=== Correlation Heatmap Rendering ===")
      message("Update button counter: ", input$update_correlation)
      
      # Isolate the parameter retrieval to avoid unwanted reactivity
      plot_params <- isolate({
        tryCatch({
          params <- correlation_plot_params()
          message("- plot_params retrieved successfully")
          params
        }, error = function(e) {
          message("Error getting plot_params: ", e$message)
          NULL
        })
      })
      
      message("Checking plot parameters:")
      message("- plot_params available: ", !is.null(plot_params))
      if (!is.null(plot_params)) {
        message("- cor_results available: ", !is.null(plot_params$cor_results))
        if (!is.null(plot_params$cor_results)) {
          message("- correlation matrix dimensions: ", paste(dim(plot_params$cor_results$correlation), collapse=" x "))
          
          # Show current filter information
          if (!is.null(plot_params$filter_info)) {
            message("- filter info available: TRUE")
            message("- active filters: ", paste(plot_params$filter_info$active_filters, collapse=", "))
            message("- original sample count: ", plot_params$filter_info$sample_count)
            message("- filtered sample count: ", plot_params$filter_info$filtered_sample_count)
            message("- filter timestamp: ", plot_params$filter_timestamp)
          }
          
          # Show correlation range for current filtered data
          cor_values <- plot_params$cor_results$correlation[upper.tri(plot_params$cor_results$correlation)]
          if (length(cor_values) > 0) {
            cor_range <- range(cor_values, na.rm = TRUE)
            message("=== CURRENT CORRELATION RANGE ===")
            message("- correlation range: ", sprintf("%.4f to %.4f", cor_range[1], cor_range[2]))
            message("- range span: ", sprintf("%.4f", cor_range[2] - cor_range[1]))
            high_corr_scaling <- (cor_range[1] > 0.8 || (cor_range[1] > 0.7 && (cor_range[2] - cor_range[1]) < 0.3))
            message("- will use high-correlation scaling: ", high_corr_scaling)
            if (high_corr_scaling) {
              message("- expected legend range: approximately -0.3 to 1.0")
            } else {
              message("- expected legend range: -1.0 to 1.0")
            }
          }
        }
        message("- metadata_for_annot available: ", !is.null(plot_params$metadata_for_annot))
        if (!is.null(plot_params$metadata_for_annot)) {
          message("- metadata dimensions: ", paste(dim(plot_params$metadata_for_annot), collapse=" x "))
        }
      }
      message("- color_scheme available: ", !is.null(input$cor_color_scheme))
      message("- show_cor_values: ", !is.null(input$show_cor_values))
      message("- show_significance: ", !is.null(input$show_significance))
      
      # Check for required inputs step by step
      if (input$update_correlation == 0) {
        message("Initial state - no correlation calculated yet")
        plot(NULL, xlim=c(0,1), ylim=c(0,1), ann=FALSE, axes=FALSE)
        text(0.5, 0.5, "Click 'Update Correlation Plot' to generate heatmap", cex=1.2)
        return()
      }
      
      if (is.null(plot_params)) {
        message("ERROR: plot_params is NULL")
        plot(NULL, xlim=c(0,1), ylim=c(0,1), ann=FALSE, axes=FALSE)
        text(0.5, 0.5, "Correlation data not available. Please update the correlation plot.", cex=1.2)
        return()
      }
      
      if (is.null(plot_params$cor_results)) {
        message("ERROR: cor_results is NULL")
        plot(NULL, xlim=c(0,1), ylim=c(0,1), ann=FALSE, axes=FALSE)
        text(0.5, 0.5, "Correlation results not available. Please recalculate.", cex=1.2)
        return()
      }
      
      if (is.null(plot_params$cor_results$correlation)) {
        message("ERROR: correlation matrix is NULL")
        plot(NULL, xlim=c(0,1), ylim=c(0,1), ann=FALSE, axes=FALSE)
        text(0.5, 0.5, "Correlation matrix not available. Please recalculate.", cex=1.2)
        return()
      }
      
      if (is.null(input$cor_color_scheme)) {
        message("ERROR: color scheme is NULL")
        plot(NULL, xlim=c(0,1), ylim=c(0,1), ann=FALSE, axes=FALSE)
        text(0.5, 0.5, "Color scheme not selected. Please select a color scheme.", cex=1.2)
        return()
      }
      
      message("\nPlot rendering started")
      
      # Validate correlation matrix dimensions
      if (is.null(plot_params$cor_results$correlation) || 
          nrow(plot_params$cor_results$correlation) < 2 || 
          ncol(plot_params$cor_results$correlation) < 2) {
        message("ERROR: Invalid correlation matrix dimensions")
        plot(NULL, xlim=c(0,1), ylim=c(0,1), ann=FALSE, axes=FALSE)
        text(0.5, 0.5, "Not enough data to display correlation heatmap.", cex=1.2)
        return()
      }
      
      message("\nRendering correlation heatmap:")
      message("- Correlation matrix dimensions: ", 
              paste(dim(plot_params$cor_results$correlation), collapse=" x "))
      if (!is.null(plot_params$metadata_for_annot)) {
        message("- Metadata dimensions: ", 
              paste(dim(plot_params$metadata_for_annot), collapse=" x "))
      } else {
        message("- No metadata available for annotations")
      }
      
      # Render the plot with comprehensive error handling
      tryCatch({
        message("About to call plot_correlation_heatmap_func...")
        
        # Safely get input values with defaults
        color_scheme <- if (!is.null(input$cor_color_scheme)) input$cor_color_scheme else "pheatmap_default"
        show_values <- if (!is.null(input$show_cor_values)) input$show_cor_values else FALSE
        show_significance <- if (!is.null(input$show_significance)) input$show_significance else TRUE
        auto_scale_colors <- if (!is.null(input$auto_scale_colors)) input$auto_scale_colors else TRUE
        
        message("Using parameters:")
        message("- color_scheme: ", color_scheme)
        message("- show_values: ", show_values)
        message("- show_significance: ", show_significance)
        message("- auto_scale_colors: ", auto_scale_colors)
        
        plot_correlation_heatmap_func(
          cor_results = plot_params$cor_results,
          metadata_df = plot_params$metadata_for_annot,
          annotation_cols = input$cor_annotations,
          color_scheme = color_scheme,
          show_significance = show_significance,
          show_values = show_values,
          auto_scale_colors = auto_scale_colors
        )
        message("plot_correlation_heatmap_func completed successfully")
      }, error = function(e) {
        message("Error in plot_correlation_heatmap_func: ", e$message)
        message("Error class: ", class(e))
        message("Error call: ", deparse(e$call))
        
        # Create error plot
        plot(NULL, xlim=c(0,1), ylim=c(0,1), ann=FALSE, axes=FALSE)
        text(0.5, 0.5, paste("Error rendering correlation heatmap:", e$message), 
             cex=1.2, col="red")
        text(0.5, 0.3, "Check console for detailed error information", 
             cex=1, col="black")
      })
      message("Heatmap rendering completed")
    }, height = 550)
    
    # Initialize Correlation submodule using current plot_data
    corr_module <- mod_correlation_server("corr", shared_data, plot_data)
    
    # Download Handlers (correlation only; PCA downloads moved into PCA submodule)
    output$download_heatmap <- downloadHandler(
      filename = function() {
        paste0("correlation_heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
      },
      content = function(file) {
        # Use current correlation results from submodule
        cor_results <- corr_module$correlation_results()
        if (is.null(cor_results)) cor_results <- calculate_correlation_func(plot_data(), input$correlation_method)
        
        save_heatmap_plot_func(
          cor_results,
          file,
          metadata_df = shared_data$processed_metadata,
          annotation_cols = NULL,
          color_scheme = input$cor_color_scheme
        )
      }
    )
    
    # PCA statistics download is handled inside the PCA submodule
    
    # Add reactive value to track QC completion status
    qc_completion_status <- reactiveVal(FALSE)
    
    # Observer to update QC completion status based on required plots and analyses
    observe({
      # Check if all necessary QC components are available
      req(plot_data(),
          pca_module$pca_results(),
          corr_module$correlation_results())
      
      # Update completion status
      qc_completion_status(TRUE)
    })
    
    # QC completion status output for the proceed section
    output$qc_completion_message <- renderUI({
      source("modules/module1_qc_preprocessing/R/shared_ui_components.R", local = TRUE)
      
      if (qc_completion_status()) {
        create_status_message(
          "QC analysis complete. You can continue to filtering and normalization.",
          status = "success"
        )
      } else {
        create_status_message(
          "Please complete all QC analyses before proceeding.",
          status = "warning"
        )
      }
    })
    
    # Expose QC completion status and data to parent module
    return(list(
      qc_completion_status = qc_completion_status,
      processed_data = plot_data,
      qc_results = list(
        pca = pca_module$pca_results,
        correlation = corr_module$correlation_results
      )
    ))
  })
} 