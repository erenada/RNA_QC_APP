#' @title Filtering & Normalization Server Module
#' @description Server logic for Tab 3: Filtering & Normalization
#' @author Eren Ada, PhD
#' @importFrom shiny moduleServer reactive reactiveVal observeEvent req withProgress
#'             incProgress showNotification updateSelectInput renderPrint
#' @importFrom stats sd
#' @importFrom utils str
#' @importFrom plotly renderPlotly ggplotly

#' @export
mod_filtering_normalization_server <- function(id, shared_data) {
  moduleServer(id, function(input, output, session) {
    
    # Source utility functions first
    source("modules/module1_qc_preprocessing/R/utils_normalization.R")
    
    # Set up DESeq2 class definitions
    tryCatch({
      setup_deseq2_classes()
    }, error = function(e) {
      showNotification(
        paste("Error setting up DESeq2 class definitions:", e$message),
        type = "error",
        duration = 10  # Auto-dismiss after 10 seconds
      )
    })
    
    # --- Reactive Values ---
    current_counts_for_tab <- reactiveVal(NULL)    # Holds counts passed to this tab
    filtered_counts_data <- reactiveVal(NULL)      # Stores the result after filtering
    normalization_evaluation_results <- reactiveVal(NULL) # Stores evaluation results
    data_state <- reactiveVal("raw")              # Tracks data state: "raw", "filtered", "normalized"
    current_operation <- reactiveVal(NULL)         # Tracks ongoing operations for progress
    
    # --- Integration with Previous Tabs ---
    observe({
      req(shared_data$data_processed)
      req(shared_data$processed_counts)
      req(shared_data$processed_metadata)
      
      # Reset states when input data changes
      current_counts_for_tab(shared_data$processed_counts)
      filtered_counts_data(NULL)
      data_state("raw")
      shared_data$normalized_counts <- NULL
      shared_data$normalization_status_flag <- FALSE
      shared_data$normalization_method_used <- NULL
      normalization_evaluation_results(NULL)
      
      # Get unique groups and update minimum samples suggestion
      if (!is.null(shared_data$processed_metadata)) {
        metadata_cols <- colnames(shared_data$processed_metadata)
        updateSelectInput(session, "filter_group_column", 
                         choices = c("None", metadata_cols), 
                         selected = "None")
        
        # If Treatment column exists, suggest appropriate minimum samples
        if ("Treatment" %in% metadata_cols) {
          group_sizes <- table(shared_data$processed_metadata$Treatment)
          min_group_size <- min(group_sizes)
          updateNumericInput(session, "filter_min_samples",
                           value = min(2, min_group_size))
        }
      }
      
      showNotification("Data loaded from previous tab", type = "message", duration = 5)
    })
    
    # --- Dynamic UI Population with Error Handling ---
    observe({
      tryCatch({
        req(shared_data$processed_metadata)
        metadata_cols <- colnames(shared_data$processed_metadata)
        updateSelectInput(session, "filter_group_column", 
                         choices = c("None", metadata_cols), 
                         selected = "None")
        updateSelectInput(session, "norm_eval_batch_variable_select", 
                         choices = metadata_cols, 
                         selected = if (length(metadata_cols) > 0) metadata_cols[1] else NULL)
      }, error = function(e) {
        showNotification(paste("Error updating UI:", e$message), type = "error", duration = 10)
      })
    })
    
    # --- Gene Filtering Logic with Progress ---
    observeEvent(input$apply_filters_btn, {
      withProgress(message = "Applying filters...", value = 0, {
        tryCatch({
          req(current_counts_for_tab())
          req(shared_data$processed_metadata)
          
          # Update progress
          incProgress(0.2, detail = "Preparing data...")
          
          counts_to_filter <- current_counts_for_tab()
          metadata <- shared_data$processed_metadata
          group_col <- if (input$filter_group_column == "None") NULL else input$filter_group_column
          
          # Validate inputs
          validate_filtering_inputs(input$filter_min_expr_threshold, 
                                  input$filter_min_samples,
                                  input$filter_expr_unit)
          
          incProgress(0.4, detail = "Filtering genes...")
          
          result_filtered_counts <- filter_low_expressed_genes_grouped(
            counts_matrix = counts_to_filter,
            metadata_df = metadata,
            group_column_name = group_col,
            min_expression = input$filter_min_expr_threshold,
            expression_unit = input$filter_expr_unit,
            min_samples_in_group_or_total = input$filter_min_samples
          )
          
          # Update states
          filtered_counts_data(result_filtered_counts$filtered_matrix)
          data_state("filtered")
          reset_normalization_state()
          
          incProgress(0.8, detail = "Updating displays...")
          
          # Update UI elements
          update_filtering_summary(result_filtered_counts)
          update_data_preview()
          
          showNotification("Filtering completed successfully", type = "message", duration = 5)
        }, error = function(e) {
          showNotification(paste("Error during filtering:", e$message), type = "error", duration = 10)
          log_error(e, "filtering")
        })
      })
    })
    
    # --- Reset Filtering Logic ---
    observeEvent(input$reset_filters_btn, {
      tryCatch({
        # Reset UI input fields to their defaults
        updateSelectInput(session, "filter_group_column", selected = "None")
        updateNumericInput(session, "filter_min_expr_threshold", value = 10)
        updateRadioButtons(session, "filter_expr_unit", selected = "raw")
        updateNumericInput(session, "filter_min_samples", value = 2)
        
        # Clear any previously filtered data
        filtered_counts_data(NULL)
        
        # Revert data state to raw (or the state of current_counts_for_tab)
        data_state("raw") # Assumes current_counts_for_tab() holds the data before filtering in this tab
        
        # Clear normalization state as well, as the input data for it might change
        reset_normalization_state()
        
        # Update UI elements
        output$filtering_summary_text <- renderPrint({ cat("Filtering settings have been reset. No filters applied from this panel.") })
        output$filtering_impact_plot <- renderPlotly({ NULL }) # Clear plot
        output$filtered_data_preview_dt <- renderDT({ NULL }) # Clear DT
        update_data_preview() # Update the main data preview at the bottom
        
        showNotification("Filtering settings and any applied filters from this panel have been reset.", 
                         type = "message", duration = 5)
      }, error = function(e) {
        showNotification(paste("Error resetting filters:", e$message), type = "error", duration = 10)
        log_error(e, "reset_filters")
      })
    })
    
    # --- Normalization Application & Evaluation ---
    observeEvent(input$apply_normalization_btn, {
      withProgress(message = "Applying normalization...", value = 0, {
        tryCatch({
          counts_for_norm <- get_current_data()
          req(counts_for_norm)
          selected_method <- input$norm_method_select
          
          # Validate method selection
          validate_normalization_method(selected_method, counts_for_norm)
          
          incProgress(0.3, detail = "Normalizing data...")
          
          norm_result <- perform_normalization(
            counts_for_norm, 
            method = selected_method
          )
          
          incProgress(0.6, detail = "Evaluating results...")
          
          # Update shared data
          shared_data$normalized_counts <- norm_result$normalized_matrix
          shared_data$normalization_status_flag <- TRUE
          shared_data$normalization_method_used <- selected_method
          data_state("normalized")
          
          # Evaluate results
          eval_results <- evaluate_normalization(
            raw_counts = counts_for_norm,
            normalized_counts = norm_result$normalized_matrix,
            metadata = shared_data$processed_metadata,
            batch_variable_name = input$norm_eval_batch_variable_select
          )
          
          normalization_evaluation_results(eval_results)
          
          incProgress(0.9, detail = "Updating visualizations...")
          
          # Update all plots and stats
          update_all_evaluation_plots()
          trigger_qc_plots_update()  # Updates plots in QC tab
          
          showNotification(paste(selected_method, "normalization completed"), type = "message", duration = 5)
        }, error = function(e) {
          showNotification(paste("Error during normalization:", e$message), type = "error", duration = 10)
          log_error(e, "normalization")
        })
      })
    })
    
    # --- Reset Normalization Logic ---
    observeEvent(input$reset_normalization_btn, {
      tryCatch({
        # Reset normalization method selection to default
        updateSelectInput(session, "norm_method_select", selected = "deseq")
        
        # Clear all normalization-related data and state
        reset_normalization_state()
        
        # Revert data state to the last stable state (filtered if available, otherwise raw)
        if (!is.null(filtered_counts_data())) {
          data_state("filtered")
        } else {
          data_state("raw")
        }
        
        # Clear all evaluation plots and statistics
        output$norm_eval_lib_size_plot <- renderPlot({
          create_error_plot("No normalization evaluation data available. Apply normalization first.")
        })
        
        output$norm_eval_lib_size_stats_text <- renderPrint({
          cat("No normalization statistics available. Apply normalization to see evaluation results.")
        })
        
        output$norm_eval_effects_plot <- renderPlot({
          create_error_plot("No normalization evaluation data available. Apply normalization first.")
        })
        
        output$norm_eval_effects_stats_text <- renderPrint({
          cat("No normalization effects statistics available. Apply normalization to see evaluation results.")
        })
        
        output$norm_eval_batch_plot <- renderPlot({
          create_error_plot("No normalization evaluation data available. Apply normalization first.")
        })
        
        output$norm_eval_batch_stats_text <- renderPrint({
          cat("No batch effect statistics available. Apply normalization to see evaluation results.")
        })
        
        output$post_norm_quick_normality_summary <- renderPrint({
          cat("No normality summary available. Apply normalization to see post-normalization normality assessment.")
        })
        
        output$normality_qq_plot <- renderPlot({
          create_error_plot("No normalization evaluation data available. Apply normalization first.")
        })
        
        # Update main data preview to reflect current state
        update_data_preview()
        
        showNotification(
          "Normalization settings and results have been reset. Data reverted to pre-normalization state.",
          type = "message", 
          duration = 5
        )
        
        # Trigger QC plots update to reflect the reset state
        trigger_qc_plots_update()
        
      }, error = function(e) {
        showNotification(paste("Error resetting normalization:", e$message), type = "error", duration = 10)
        log_error(e, "reset_normalization")
      })
    })
    
    # --- Helper Functions ---
    get_current_data <- reactive({
      switch(data_state(),
             "raw" = current_counts_for_tab(),
             "filtered" = filtered_counts_data(),
             "normalized" = shared_data$normalized_counts
      )
    })
    
    reset_normalization_state <- function() {
      shared_data$normalized_counts <- NULL
      shared_data$normalization_status_flag <- FALSE
      shared_data$normalization_method_used <- NULL
      normalization_evaluation_results(NULL)
    }
    
    log_error <- function(error, context) {
      error_log <- list(
        timestamp = Sys.time(),
        context = context,
        message = error$message,
        call = error$call,
        data_state = data_state()
      )
      # Write to error log file
      log_dir <- file.path("logs", "filtering_normalization")
      if (!dir.exists(log_dir)) {
        dir.create(log_dir, recursive = TRUE)
      }
      log_file <- file.path(log_dir, format(Sys.Date(), "errors_%Y%m%d.log"))
      cat(sprintf("[%s] %s: %s\n", 
                  format(error_log$timestamp), 
                  error_log$context,
                  error_log$message),
          file = log_file, 
          append = TRUE)
    }
    
    trigger_qc_plots_update <- function() {
      # Trigger updates in QC plots tab if the function exists
      if (exists("mod_qc_plots_update")) {
        mod_qc_plots_update()
      }
    }
    
    # --- UI Update Functions ---
    update_filtering_summary <- function(result) {
      # Update text summary with error handling
      output$filtering_summary_text <- renderPrint({
        tryCatch({
          print(result)  # Uses the custom print method for filtering results
        }, error = function(e) {
          cat("Error displaying filtering summary:", e$message)
        })
      })
      
      # Update gene counts plot with error handling
      output$filtering_gene_counts_plot <- renderPlot({
        tryCatch({
          validate(need(!is.null(result$genes_before) && !is.null(result$genes_after),
                       "Required gene count data not available"))
          
          p <- plot_filtering_gene_counts(
            initial_genes = result$genes_before,
            filtered_genes = result$genes_after
          )
          
          if (!is.null(p)) {
            p
          } else {
            create_error_plot("Error creating gene counts plot")
          }
        }, error = function(e) {
          create_error_plot(sprintf("Error: %s", e$message))
        })
      })
      
      # Update expression distribution plot with error handling
      output$filtering_expression_dist_plot <- renderPlot({
        tryCatch({
          validate(need(!is.null(result$filtered_matrix), 
                       "Filtered data not available"))
          
          # Get the original counts data
          original_counts <- current_counts_for_tab()
          
          # Create density plot comparing before and after filtering
          p <- plot_filtering_density(
            before_data = original_counts,
            after_data = result$filtered_matrix
          )
          
          if (!is.null(p)) {
            p
          } else {
            create_error_plot("Error creating expression distribution plot")
          }
        }, error = function(e) {
          create_error_plot(sprintf("Error: %s", e$message))
        })
      })
      
      # Update data preview with error handling
      output$filtered_data_preview_dt <- renderDT({
        tryCatch({
          validate(need(!is.null(result$filtered_matrix), "Filtered data not available"))
          
          datatable(
            head(result$filtered_matrix, n = 10),
            options = list(
              scrollX = TRUE,
              dom = 't',
              pageLength = 10
            )
          )
        }, error = function(e) {
          NULL
        })
      })
    }
    
    update_data_preview <- function() {
      output$main_data_preview_dt <- renderDT({
        current_data <- get_current_data()
        if (is.null(current_data)) {
          return(NULL)
        }
        
        preview_data <- current_data[1:min(100, nrow(current_data)), 1:min(20, ncol(current_data)), drop = FALSE]
        
        datatable(
          preview_data,
          options = list(
            scrollX = TRUE,
            pageLength = 10,
            lengthMenu = list(c(10, 25, 50), c('10', '25', '50')),
            dom = 'lrtip'
          ),
          rownames = TRUE,
          class = 'cell-border stripe'
        ) %>%
          formatRound(columns = 1:ncol(preview_data), digits = 2)
      })
    }
    
    update_all_evaluation_plots <- function() {
      eval_results <- normalization_evaluation_results()
      if (is.null(eval_results)) return()
      
      # Add debug output
      cat("Updating evaluation plots with available data\n")
      cat("Structure of eval_results:\n")
      str(eval_results, max.level = 1)
      
      # --- Library Size Evaluation Tab ---
      output$norm_eval_lib_size_plot <- renderPlot({
        cat("Rendering library size plot\n")
        plot_library_size_distribution(eval_results$library_size$data)
      })
      
      output$norm_eval_lib_size_stats_text <- renderPrint({
        cat("Rendering library size stats\n")
        stats_text <- format_library_size_stats(eval_results$library_size$stats)
        cat(stats_text)
      })
      
      # --- Normalization Effects Tab ---
      output$norm_eval_effects_plot <- renderPlot({
        cat("Rendering normalization effects plot\n")
        plot_normalization_density(eval_results)
      })
      
      output$norm_eval_effects_stats_text <- renderPrint({
        cat("Rendering normalization effects stats\n")
        cat(get_effects_stats(eval_results))
      })
      
      # --- Batch Effect Check Tab ---
      output$norm_eval_batch_plot <- renderPlot({
        cat("Rendering batch effect plot\n")
        plot_normalization_pca(eval_results)
      })
      
      output$norm_eval_batch_stats_text <- renderPrint({
        cat("Rendering batch effect stats\n")
        cat(get_batch_stats(eval_results))
      })
      
      # --- Normality Check Tab ---
      output$post_norm_quick_normality_summary <- renderPrint({
        cat("Rendering normality summary\n")
        cat("Normality Test Results (Shapiro-Wilk):\n\n")
        
        if (!is.null(eval_results$normality) && 
            !is.null(eval_results$normality$shapiro_test_results)) {
          # Print first few samples
          print(head(eval_results$normality$shapiro_test_results, n = 5))
        } else {
          cat("Normality test results not available\n")
        }
      })
      
      output$normality_qq_plot <- renderPlot({
        cat("Rendering normality QQ plot\n")
        
        if (!is.null(eval_results$normality) && 
            !is.null(eval_results$normality$qq_data)) {
          plot_qq(eval_results$normality$qq_data)
        } else {
          create_error_plot("Normality QQ data not available")
        }
      })
    }
    
    # --- Normality Assessment Logic ---
    # Normality Assessment
    normality_assessment_results_val <- reactiveVal(NULL)
    
    # Sample selection UI for normality assessment
    output$normality_sample_selection <- renderUI({
      current_data <- get_current_data()
      req(current_data)
      
      # Get sample names from the data
      sample_names <- colnames(current_data)
      
      # Update the selectizeInput choices
      updateSelectizeInput(session, "normality_samples_to_test",
        choices = sample_names,
        selected = sample_names[1:min(5, length(sample_names))]
      )
    })
    
    # Update normality assessment when button is clicked
    observeEvent(input$update_normality, {
      current_data <- get_current_data()
      req(current_data, input$normality_test_method, input$normality_significance)
      
      # Get selected samples (or use all if none selected)
      selected_samples <- input$normality_samples_to_test
      if (is.null(selected_samples) || length(selected_samples) == 0) {
        selected_samples <- colnames(current_data)
      }
      
      # Run normality assessment with selected parameters
      results <- assess_normality_func(
        counts_data = current_data,
        samples_to_test = selected_samples,
        test_method = input$normality_test_method,
        significance_level = as.numeric(input$normality_significance)
      )
      
      normality_assessment_results_val(results)
    })
    
    # Initialize normality assessment with default data
    observe({
      current_data <- get_current_data()
      req(current_data)
      results <- assess_normality_func(current_data)
      normality_assessment_results_val(results)
      
      # Update sample choices
      sample_names <- colnames(current_data)
      updateSelectizeInput(session, "normality_samples_to_test",
        choices = sample_names,
        selected = sample_names[1:min(5, length(sample_names))]
      )
    })
    
    # Normality Summary Text - Concise version
    output$normality_summary_text <- renderPrint({
      req(normality_assessment_results_val())
      
      results <- normality_assessment_results_val()
      summary_stats <- results$summary_stats
      
      cat("=== Quick Assessment Summary ===\n\n")
      cat(sprintf("Total samples analyzed: %d\n", summary_stats$total_samples))
      cat(sprintf("Normal distribution: %d samples (%.1f%%)\n", 
                 summary_stats$normal_count, summary_stats$normal_pct))
      cat(sprintf("Moderate deviation: %d samples (%.1f%%)\n", 
                 summary_stats$moderate_count, summary_stats$moderate_pct))
      cat(sprintf("Severe deviation: %d samples (%.1f%%)\n", 
                 summary_stats$severe_count, summary_stats$severe_pct))
      
      cat("\n=== Test Details ===\n")
      cat(sprintf("Method: %s test\n", summary_stats$test_method))
      cat(sprintf("Significance level: α = %s\n", summary_stats$significance_level))
      
      cat("\n=== Recommendation ===\n")
      cat(sprintf("Suggested correlation method: %s\n", 
                 results$recommendation$correlation_method))
      
      if (summary_stats$total_samples <= 10) {
        cat("\n=== Sample-by-Sample Results ===\n")
        sample_names <- names(results$sample_stats)
        for (sample in sample_names) {
          stats <- results$sample_stats[[sample]]
          overall_status <- if (stats$has_severe_deviation) {
            "Severe deviation"
          } else if (stats$has_moderate_deviation) {
            "Moderate deviation"
          } else {
            "Normal"
          }
          cat(sprintf("%s: %s\n", sample, overall_status))
        }
      } else {
        cat(sprintf("\nFor detailed per-sample results (%d samples),\n", summary_stats$total_samples))
        cat("see the 'Detailed Results' tab.\n")
      }
    })
    
    # Normality Recommendation Text
    output$normality_recommendation_text <- renderUI({
      req(normality_assessment_results_val())
      
      results <- normality_assessment_results_val()
      
      # Create a styled recommendation box
      recommendation_color <- switch(
        results$recommendation$correlation_method,
        "pearson" = "#28a745",  # Green for Pearson
        "spearman" = "#dc3545", # Red for Spearman (when strongly recommended)
        "#ffc107"               # Yellow/amber as default
      )
      
      # Create the recommendation box with appropriate styling
      div(
        style = paste0(
          "margin: 15px 0; padding: 15px; border-radius: 5px; ",
          "background-color: ", recommendation_color, "20; ", # 20 for transparency
          "border-left: 5px solid ", recommendation_color, "; ",
          "color: #333;"
        ),
        h4(style = "margin-top: 0; color: #333;", "Correlation Method Recommendation"),
        p(style = "font-weight: bold; font-size: 16px;", results$recommendation$summary),
        p(results$recommendation$reasoning),
        hr(style = "margin: 10px 0; border-color: #ddd;"),
        p(style = "font-style: italic; margin-bottom: 0;", 
          "Based on normality assessment using ", strong(results$summary_stats$test_method), 
          " test with significance level α = ", results$summary_stats$significance_level)
      )
    })
    
    # Normality Q-Q Plots
    output$normality_qq_plots <- renderPlot({
      current_data <- get_current_data()
      req(current_data)
      
      # Use selected samples if available
      samples_to_plot <- input$normality_samples_to_test
      if (is.null(samples_to_plot) || length(samples_to_plot) == 0) {
        # Default to first few samples if none selected
        samples_to_plot <- colnames(current_data)[1:min(6, ncol(current_data))]
      }
      
      plot_normality_qq_func(counts_data = current_data, samples_to_plot = samples_to_plot)
    }, height = 600)
    
    # Normality Density Plots
    output$normality_density_plots <- renderPlot({
      current_data <- get_current_data()
      req(current_data)
      
      # Use selected samples if available
      samples_to_plot <- input$normality_samples_to_test
      if (is.null(samples_to_plot) || length(samples_to_plot) == 0) {
        # Default to first few samples if none selected
        samples_to_plot <- colnames(current_data)[1:min(6, ncol(current_data))]
      }
      
      plot_normality_density_func(counts_data = current_data, samples_to_plot = samples_to_plot)
    }, height = 600)
    
    # Normality Test Results
    output$normality_test_results <- renderPrint({
      req(normality_assessment_results_val())
      
      results <- normality_assessment_results_val()
      
      # Print overall assessment
      cat("Normality Assessment Summary:\n")
      cat("==========================\n\n")
      
      # Print summary statistics
      summary_stats <- results$summary_stats
      cat("Overall Assessment:\n")
      cat("------------------\n")
      cat(sprintf("Total samples analyzed: %d\n", summary_stats$total_samples))
      cat(sprintf("Samples with normal distribution: %d (%.1f%%)\n", 
                 summary_stats$normal_count, summary_stats$normal_pct))
      cat(sprintf("Samples with moderate deviation: %d (%.1f%%)\n", 
                 summary_stats$moderate_count, summary_stats$moderate_pct))
      cat(sprintf("Samples with severe deviation: %d (%.1f%%)\n", 
                 summary_stats$severe_count, summary_stats$severe_pct))
      
      # Print recommendation
      cat("\nRecommendation:\n")
      cat("--------------\n")
      cat(results$recommendation$summary, "\n")
      cat("Reasoning: ", results$recommendation$reasoning, "\n")
      
      # Print assessment criteria
      cat("\nAssessment Criteria:\n")
      cat("------------------\n")
      cat("Moderate deviation if ANY of:\n")
      cat("- |Skewness| > 2\n")
      cat("- |Kurtosis - 3| > 7\n")
      cat("- ", results$summary_stats$test_method, " p < 0.01\n\n")
      cat("Severe deviation if ANY of:\n")
      cat("- |Skewness| > 3\n")
      cat("- |Kurtosis - 3| > 10\n")
      cat("- ", results$summary_stats$test_method, " p < 0.001\n")
      
      # Print detailed results for first few samples
      cat("\nDetailed Results (First 5 Samples):\n")
      cat("--------------------------------\n")
      
      sample_names <- names(results$sample_stats)
      for (sample in sample_names[1:min(5, length(sample_names))]) {
        stats <- results$sample_stats[[sample]]
        cat(sprintf("\n%s:\n", sample))
        cat(sprintf("Skewness: %.2f (%s)\n", stats$skewness, stats$skewness_status))
        cat(sprintf("Kurtosis: %.2f (%s)\n", stats$kurtosis, stats$kurtosis_status))
        
        if (!is.na(stats$p_value)) {
          cat(sprintf("%s p-value: %.3e (%s)\n", 
                     stats$test_name, stats$p_value, stats$p_value_status))
        } else {
          cat(sprintf("%s p-value: NA (not calculated)\n", stats$test_name))
        }
        
        # Overall assessment for this sample
        if (stats$has_severe_deviation) {
          cat("Overall: Severe deviation from normality\n")
        } else if (stats$has_moderate_deviation) {
          cat("Overall: Moderate deviation from normality\n")
        } else {
          cat("Overall: Approximately normal\n")
        }
      }
    })
    
    # Add histogram with density plots
    output$normality_hist_plots <- renderPlot({
      current_data <- get_current_data()
      req(current_data)
      
      # Use selected samples if available
      samples_to_plot <- input$normality_samples_to_test
      if (is.null(samples_to_plot) || length(samples_to_plot) == 0) {
        # Default to first few samples if none selected
        samples_to_plot <- colnames(current_data)[1:min(6, ncol(current_data))]
      }
      
      plot_normality_hist_density_func(counts_data = current_data, samples_to_plot = samples_to_plot)
    }, height = 600)
    
    # --- Download Handlers ---
    output$download_current_data <- downloadHandler(
      filename = function() {
        paste0("filtered_normalized_data_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
      },
      content = function(file) {
        # Get current data state
        current_data <- get_current_data()
        req(current_data)
        
        # Create output data frame with gene IDs as first column
        output_data <- data.frame(
          Gene_ID = rownames(current_data),
          current_data,
          check.names = FALSE
        )
        
        # Write to CSV
        write.csv(output_data, file, row.names = FALSE)
      }
    )
    
    output$download_session_report <- downloadHandler(
      filename = function() {
        paste0("analysis_report_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".html")
      },
      content = function(file) {
        # Get current data and states
        current_data <- get_current_data()
        req(current_data)
        
        # Create report content
        report_content <- list(
          title = "RNA-Seq Data Analysis Report",
          author = "Eren Ada, PhD",
          date = format(Sys.time(), "%Y-%m-%d"),
          data_state = data_state(),
          filtering_params = if(!is.null(input$filter_min_samples)) {
            list(
              min_samples = input$filter_min_samples,
              min_counts = input$filter_min_counts,
              group_column = input$filter_group_column
            )
          } else NULL,
          normalization_params = if(!is.null(input$norm_method_select)) {
            list(
              method = input$norm_method_select,
              evaluation_results = normalization_evaluation_results()
            )
          } else NULL,
          data_summary = list(
            dimensions = dim(current_data),
            total_counts = sum(current_data),
            mean_counts = mean(current_data),
            median_counts = median(current_data)
          )
        )
        
        # Generate HTML report
        rmarkdown::render(
          input = "modules/module1_qc_preprocessing/R/templates/analysis_report_template.Rmd",
          output_file = file,
          params = list(
            report_content = report_content
          ),
          quiet = TRUE
        )
      }
    )
    
    # Filtering & Normalization completion message for proceed section
    output$filtering_normalization_completion_message <- renderUI({
      source("modules/module1_qc_preprocessing/R/shared_ui_components.R", local = TRUE)
      
      if (shared_data$normalization_status_flag) {
        create_status_message(
          paste("Data filtering and normalization complete using", shared_data$normalization_method_used, "method. Ready for export and downstream analysis."),
          status = "success"
        )
      } else if (!is.null(filtered_counts_data())) {
        create_status_message(
          "Data filtering complete. You can proceed to export your filtered data, or optionally apply normalization first.",
          status = "success"
        )
      } else if (!is.null(current_counts_for_tab())) {
        create_status_message(
          "Data loaded and ready for export. Filtering and normalization are optional processing steps.",
          status = "success"
        )
      } else {
        create_status_message(
          "Please complete the previous steps before accessing this tab.",
          status = "error"
        )
      }
    })
    
    # Handle proceed to next step button
    observeEvent(input$proceed_to_next, {
      if (shared_data$normalization_status_flag) {
        # Navigate to next tab (when implemented) or show download options
        showNotification(
          "Filtering and normalization complete! Use the download buttons to export your processed data.",
          type = "message",
          duration = 5
        )
        
        # For now, since there's no next tab, just highlight the download section
        # In the future, this would navigate to a results/export tab
        
      } else if (!is.null(filtered_counts_data())) {
        # Allow proceeding with filtering only (normalization is optional)
        showNotification(
          "Filtering complete! Use the download buttons to export your filtered data. Note: normalization was not applied.",
          type = "message",
          duration = 5
        )
        
        # For now, since there's no next tab, just highlight the download section
        # In the future, this would navigate to a results/export tab
        
      } else if (!is.null(current_counts_for_tab())) {
        # Allow proceeding with just the loaded data (no filtering or normalization applied)
        showNotification(
          "Data ready for export! Use the download buttons to export your data. Note: no additional filtering or normalization was applied in this step.",
          type = "message",
          duration = 5
        )
        
        # For now, since there's no next tab, just highlight the download section
        # In the future, this would navigate to a results/export tab
        
      } else {
        showNotification(
          "No data available. Please complete the previous steps first.",
          type = "error",
          duration = 5
        )
      }
    })
    
    # Return completion status for use in main app
    return(list(
      filtering_complete = reactive({ !is.null(filtered_counts_data()) }),
      normalization_complete = reactive({ shared_data$normalization_status_flag }),
      current_data = get_current_data,
      data_state = data_state
    ))
  })
}

# --- Helper Functions ---
create_error_plot <- function(error_message) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, 
             label = error_message,
             color = "red",
             size = 5) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(20, 20, 20, 20)
    )
} 