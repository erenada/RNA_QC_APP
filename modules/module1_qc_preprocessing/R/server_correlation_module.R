# server_correlation_module.R
# Correlation Submodule Server

#' @title Server for Correlation Submodule
#' @param id Module ID
#' @param shared_data ReactiveValues list (metadata, normalized data flags)
#' @param plot_data_reactive Reactive returning current data matrix for correlation
#' @return A list with reactive accessors: correlation_results(), data_description()
#' @export
mod_correlation_server <- function(id, shared_data, plot_data_reactive) {
  moduleServer(id, function(input, output, session) {
    # Utilities expected: calculate_correlation_func, plot_correlation_heatmap_func


    # Filter criteria UI
    output$filter_criteria <- renderUI({
      req(shared_data$processed_metadata)
      metadata_cols <- colnames(shared_data$processed_metadata)
      checkboxGroupInput(session$ns("filter_criteria"), "Select Filters to Apply:", choices = metadata_cols, selected = NULL)
    })

    # Dynamic filters UI
    output$dynamic_filters <- renderUI({
      req(shared_data$processed_metadata, input$filter_criteria)
      filter_list <- lapply(input$filter_criteria, function(col) {
        unique_values <- unique(shared_data$processed_metadata[[col]])
        selectInput(session$ns(paste0("filter_", col)), paste("Select", col, ":"), choices = unique_values, selected = unique_values, multiple = TRUE)
      })
      do.call(tagList, filter_list)
    })

    # Selected samples reactive
    selected_samples <- reactive({
      req(shared_data$processed_metadata)
      all_samples <- rownames(shared_data$processed_metadata)
      if (is.null(input$filter_criteria) || length(input$filter_criteria) == 0) return(all_samples)
      selected <- all_samples
      for (col in input$filter_criteria) {
        vals <- input[[paste0("filter_", col)]]
        if (!is.null(vals) && length(vals) > 0 && col %in% colnames(shared_data$processed_metadata)) {
          matching <- rownames(shared_data$processed_metadata)[shared_data$processed_metadata[[col]] %in% vals]
          selected <- intersect(selected, matching)
        }
      }
      if (length(selected) == 0) all_samples else selected
    })

    # Correlation parameters and results
    correlation_plot_params <- reactiveVal(NULL)
    correlation_results <- reactiveVal(NULL)
    data_description <- reactiveVal("No data selected")

    observeEvent(input$update_correlation, {
      req(shared_data$processed_metadata)

      # Choose data: normalized if available and flagged, else plot_data_reactive()
      if (!is.null(shared_data$normalized_counts) && isTRUE(shared_data$normalization_status_flag)) {
        correlation_data <- shared_data$normalized_counts
        desc <- paste(shared_data$normalization_method_used, "normalized data")
      } else {
        req(plot_data_reactive())
        correlation_data <- plot_data_reactive()
        desc <- "log2(CPM+1) transformed data"
      }

      # Apply filters to columns (samples)
      meta <- shared_data$processed_metadata
      keep_samples <- selected_samples()
      common <- intersect(colnames(correlation_data), keep_samples)
      if (length(common) < 2) {
        showNotification("Not enough samples selected for correlation analysis (minimum 2 required)", type = "error", duration = 5)
        return()
      }
      correlation_data <- correlation_data[, common, drop = FALSE]
      meta <- meta[rownames(meta) %in% common, , drop = FALSE]

      # Calculate correlation
      cor_res <- tryCatch({
        calculate_correlation_func(correlation_data, input$correlation_method, compute_p = isTRUE(input$compute_p_values))
      }, error = function(e) {
        showNotification(paste("Error calculating correlation:", e$message), type = "error", duration = 5)
        NULL
      })
      req(cor_res)

      correlation_plot_params(list(cor_results = cor_res, metadata_for_annot = meta))
      correlation_results(cor_res)
      data_description(desc)
      showNotification("Correlation plot updated successfully!", type = "message", duration = 3)
    })

    # Heatmap render
    output$correlation_heatmap <- renderPlot({
      input$update_correlation
      params <- isolate(correlation_plot_params())
      if (is.null(params) || is.null(params$cor_results) || is.null(params$cor_results$correlation)) {
        plot(NULL, xlim=c(0,1), ylim=c(0,1), ann=FALSE, axes=FALSE)
        text(0.5, 0.5, "Click 'Update Correlation Plot' to generate heatmap", cex=1.2)
        return()
      }
      color_scheme <- if (!is.null(input$cor_color_scheme)) input$cor_color_scheme else "default"
      show_values <- isTRUE(input$show_cor_values)
      show_significance <- isTRUE(input$show_significance)
      auto_scale_colors <- isTRUE(input$auto_scale_colors)
      plot_correlation_heatmap_func(
        cor_results = params$cor_results,
        metadata_df = params$metadata_for_annot,
        annotation_cols = NULL,
        color_scheme = color_scheme,
        show_significance = show_significance,
        show_values = show_values,
        auto_scale_colors = auto_scale_colors
      )
    }, height = 550)

    # Stats summary
    output$correlation_stats_text <- renderUI({
      req(correlation_results())
      cor_res <- correlation_results()
      # Guard invalid correlation matrix
      if (is.null(cor_res$correlation) || nrow(cor_res$correlation) < 2 || ncol(cor_res$correlation) < 2) {
        return(tags$div(style = "margin-top: 15px; padding: 10px; background-color: #f8f9fa; border-radius: 5px;",
                        tags$h4("Correlation Summary Statistics:", style = "margin-bottom: 15px;"),
                        tags$p("Not enough samples selected or data available to calculate correlation statistics.")))
      }
      cor_values <- cor_res$correlation[upper.tri(cor_res$correlation)]
      mean_cor <- if (length(cor_values) > 0) mean(cor_values, na.rm = TRUE) else NA
      median_cor <- if (length(cor_values) > 0) median(cor_values, na.rm = TRUE) else NA
      min_cor <- if (length(cor_values) > 0) min(cor_values, na.rm = TRUE) else NA
      max_cor <- if (length(cor_values) > 0) max(cor_values, na.rm = TRUE) else NA
      if (!is.null(cor_res$p_adjusted) && all(dim(cor_res$p_adjusted) == dim(cor_res$correlation))) {
        sig_cors <- sum(cor_res$p_adjusted < 0.05 & upper.tri(cor_res$correlation), na.rm = TRUE)
      } else {
        sig_cors <- NA
      }
      total_cors <- sum(upper.tri(cor_res$correlation))
      tagList(
        tags$div(style = "margin-top: 15px; padding: 10px; background-color: #f8f9fa; border-radius: 5px;",
          tags$h4("Correlation Summary Statistics:", style = "margin-bottom: 15px;"),
          tags$div(style = "margin-bottom: 10px; padding: 8px; background-color: #e3f2fd; border-radius: 4px;",
                   tags$strong("Data Source: "), tags$span(isolate(data_description()), style = "font-weight: normal;")),
          fluidRow(
            column(width = 6,
              tags$ul(style = "list-style-type: none; padding-left: 0;",
                tags$li(sprintf("Mean correlation: %.3f", mean_cor)),
                tags$li(sprintf("Median correlation: %.3f", median_cor)),
                tags$li(sprintf("Range: %.3f to %.3f", min_cor, max_cor))
              )
            ),
            column(width = 6,
              tags$ul(style = "list-style-type: none; padding-left: 0;",
                tags$li(sprintf("Significant correlations: %s", ifelse(is.na(sig_cors), "NA", sig_cors))),
                tags$li(sprintf("Total correlations: %d", total_cors)),
                tags$li(sprintf("Significance rate: %s", ifelse(is.na(sig_cors), "NA", sprintf("%.1f%%", 100*sig_cors/total_cors))))
              )
            )
          )
        )
      )
    })

    # Download handler
    output$download_heatmap <- downloadHandler(
      filename = function() paste0("correlation_heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) {
        cor_results <- correlation_results()
        if (is.null(cor_results)) {
          cor_results <- calculate_correlation_func(plot_data_reactive(), input$correlation_method)
        }
        save_heatmap_plot_func(cor_results, file, metadata_df = shared_data$processed_metadata, annotation_cols = NULL, color_scheme = input$cor_color_scheme)
      }
    )

    # Expose public API
    return(list(
      correlation_results = reactive({ correlation_results() }),
      data_description = reactive({ data_description() })
    ))
  })
}
