# server_pca_module.R
# PCA Submodule Server

#' @title Server for PCA Submodule
#' @param id Module ID
#' @param shared_data ReactiveValues list with processed_metadata
#' @param plot_data_reactive A reactive returning the current matrix used for PCA
#' @return A list with reactive accessors, including pca_results()
#' @export
mod_pca_server <- function(id, shared_data, plot_data_reactive) {
  moduleServer(id, function(input, output, session) {
    # Import utilities from utils_qc_plotting
    # Requires perform_pca_func and plot_pca_func to be available in scope

    # Expose color-by UI based on metadata
    output$pca_color_by_ui <- renderUI({
      req(shared_data$processed_metadata)
      metadata_cols <- colnames(shared_data$processed_metadata)
      default_selection <- if (length(metadata_cols) > 0) metadata_cols[1] else "none"
      selectInput(
        session$ns("pca_color_by"),
        "Color by:",
        choices = c("None" = "none", "Sample" = "Sample", metadata_cols),
        selected = default_selection
      )
    })

    # PCA results storage
    pca_results_data <- reactiveVal(NULL)

    # 2D PCA
    output$pca_plot_2d <- renderPlotly({
      req(plot_data_reactive(), input$pca_pc_x, input$pca_pc_y)
      pca_res <- tryCatch({
        perform_pca_func(counts_data = plot_data_reactive(), center = TRUE, scale. = TRUE)
      }, error = function(e) {
        message("Error in PCA calculation: ", e$message)
        NULL
      })
      req(pca_res)
      pca_results_data(pca_res)

      color_var <- NULL
      if (!is.null(input$pca_color_by) && input$pca_color_by != "none") color_var <- input$pca_color_by

      plot_pca_func(
        pca_results = pca_res,
        pc_x_choice = input$pca_pc_x,
        pc_y_choice = input$pca_pc_y,
        metadata_df = shared_data$processed_metadata,
        color_var = color_var,
        show_labels = input$pca_show_labels,
        plot_type = "2d"
      )
    })

    # 3D PCA
    output$pca_plot_3d <- renderPlotly({
      req(plot_data_reactive(), input$pca_pc_x, input$pca_pc_y)
      pca_res <- pca_results_data()
      if (is.null(pca_res)) {
        pca_res <- tryCatch({
          perform_pca_func(counts_data = plot_data_reactive(), center = TRUE, scale. = TRUE)
        }, error = function(e) {
          message("Error in PCA calculation: ", e$message)
          NULL
        })
        req(pca_res)
        pca_results_data(pca_res)
      }

      pc_z_choice <- paste0("PC", as.numeric(gsub("PC", "", input$pca_pc_y)) + 1)
      if (pc_z_choice == input$pca_pc_x || pc_z_choice == input$pca_pc_y) pc_z_choice <- "PC3"

      color_var <- NULL
      if (!is.null(input$pca_color_by) && input$pca_color_by != "none") color_var <- input$pca_color_by

      plot_pca_func(
        pca_results = pca_res,
        pc_x_choice = input$pca_pc_x,
        pc_y_choice = input$pca_pc_y,
        pc_z_choice = pc_z_choice,
        metadata_df = shared_data$processed_metadata,
        color_var = color_var,
        show_labels = FALSE,
        plot_type = "3d"
      )
    })

    # Variance explained UI
    output$pca_variance_explained_text <- renderUI({
      req(pca_results_data())
      pca_res_obj <- pca_results_data()
      var_explained <- pca_res_obj$var_explained[1:min(10, length(pca_res_obj$var_explained))]
      cum_var_explained <- cumsum(var_explained)
      tagList(
        h4("PCA Statistics", class = "section-title"),
        fluidRow(
          column(width = 6,
            tags$div(style = "margin-top: 15px; padding: 10px; background-color: #f8f9fa; border-radius: 5px; height: 100%;",
              tags$p(strong("Individual PC Contributions:"), style = "margin-bottom: 10px;"),
              tags$ul(style = "list-style-type: none; padding-left: 0;",
                lapply(1:length(var_explained), function(i) {
                  tags$li(sprintf("PC%d: %.1f%% variance", i, var_explained[i]), style = "margin-bottom: 5px;")
                })
              )
            )
          ),
          column(width = 6,
            tags$div(style = "margin-top: 15px; padding: 10px; background-color: #f8f9fa; border-radius: 5px; height: 100%;",
              tags$p(strong("Cumulative Variance Explained:"), style = "margin-bottom: 10px;"),
              tags$ul(style = "list-style-type: none; padding-left: 0;",
                lapply(1:length(cum_var_explained), function(i) {
                  tags$li(sprintf("PC1 to PC%d: %.1f%%", i, cum_var_explained[i]), style = "margin-bottom: 5px;")
                })
              )
            )
          )
        )
      )
    })

    # Downloads
    output$download_pca_plot <- downloadHandler(
      filename = function() paste0("pca_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
      content = function(file) {
        req(pca_results_data())
        pca_res <- pca_results_data()
        plot_data <- data.frame(
          PC1 = pca_res$pca_object$x[, as.numeric(gsub("PC", "", input$pca_pc_x))],
          PC2 = pca_res$pca_object$x[, as.numeric(gsub("PC", "", input$pca_pc_y))],
          Sample = rownames(pca_res$pca_object$x)
        )
        if (!is.null(input$pca_color_by) && input$pca_color_by != "none") {
          if (input$pca_color_by == "Sample") {
            plot_data$group <- plot_data$Sample
          } else if (input$pca_color_by %in% colnames(shared_data$processed_metadata)) {
            plot_data$group <- shared_data$processed_metadata[as.character(plot_data$Sample), input$pca_color_by]
          } else {
            plot_data$group <- "All Samples"
          }
        } else {
          plot_data$group <- "All Samples"
        }
        plot_data$group[is.na(plot_data$group)] <- "NA"
        p <- ggplot(plot_data, aes(x = PC1, y = PC2, color = group)) +
          geom_point(size = 3) +
          labs(title = "PCA Plot",
               x = sprintf("PC%s (%.1f%%)", gsub("PC", "", input$pca_pc_x), pca_res$var_explained[as.numeric(gsub("PC", "", input$pca_pc_x))]),
               y = sprintf("PC%s (%.1f%%)", gsub("PC", "", input$pca_pc_y), pca_res$var_explained[as.numeric(gsub("PC", "", input$pca_pc_y))])) +
          theme_minimal() +
          theme(legend.title = element_blank())
        if (input$pca_show_labels) {
          p <- p + ggrepel::geom_text_repel(aes(label = Sample), size = 3, box.padding = 0.5)
        }
        ggsave(file, plot = p, width = 8, height = 6)
      }
    )

    output$download_pca_stats <- downloadHandler(
      filename = function() paste0("pca_statistics_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip"),
      content = function(file) {
        req(pca_results_data())
        pca_res <- pca_results_data()
        pc_x_num <- as.numeric(gsub("PC", "", input$pca_pc_x))
        pc_y_num <- as.numeric(gsub("PC", "", input$pca_pc_y))
        pc_z_num <- as.numeric(gsub("PC", "", input$pca_pc_y)) + 1
        if (pc_z_num == pc_x_num || pc_z_num == pc_y_num) pc_z_num <- 3
        var_explained <- pca_res$var_explained[1:min(20, length(pca_res$var_explained))]
        cum_var_explained <- cumsum(var_explained)
        pca_stats <- data.frame(
          PC = paste0("PC", 1:length(var_explained)),
          Individual_Variance_Pct = round(var_explained, 2),
          Cumulative_Variance_Pct = round(cum_var_explained, 2),
          row.names = NULL
        )
        loadings <- pca_res$pca_object$rotation[,1:length(var_explained)]
        contributions <- sweep(loadings^2, 2, pca_res$pca_object$sdev[1:length(var_explained)]^2, "*")
        top_genes_list <- list()
        for (i in 1:min(10, ncol(contributions))) {
          pc_contributions <- contributions[, i]
          sorted_idx <- order(pc_contributions, decreasing = TRUE)
          top_idx <- sorted_idx[1:min(50, length(sorted_idx))]
          top_genes_df <- data.frame(
            Gene = rownames(loadings)[top_idx],
            Contribution_Pct = round(pc_contributions[top_idx]/sum(pc_contributions)*100, 2)
          )
          top_genes_list[[paste0("PC", i, "_Top_Genes")]] <- top_genes_df
        }
        pc_display_info <- data.frame(
          View = c("2D_X_Axis", "2D_Y_Axis", "3D_Z_Axis"),
          PC = c(paste0("PC", pc_x_num), paste0("PC", pc_y_num), paste0("PC", pc_z_num)),
          Variance_Explained_Pct = c(round(var_explained[pc_x_num], 2), round(var_explained[pc_y_num], 2), round(var_explained[pc_z_num], 2))
        )
        output_list <- list(PCA_Display_Info = pc_display_info, PCA_Variance_Statistics = pca_stats)
        output_list <- c(output_list, top_genes_list)
        
        # Create a dedicated folder to avoid deep nested paths in ZIP
        temp_dir <- tempdir()
        stats_dir_name <- "pca_statistics"
        stats_dir <- file.path(temp_dir, stats_dir_name)
        if (dir.exists(stats_dir)) {
          unlink(list.files(stats_dir, full.names = TRUE), recursive = TRUE, force = TRUE)
        } else {
          dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)
        }
        
        # Write all CSVs into the dedicated folder
        for (name in names(output_list)) {
          csv_file <- file.path(stats_dir, paste0(name, ".csv"))
          write.csv(output_list[[name]], file = csv_file, row.names = FALSE)
        }
        
        # Zip the folder so the ZIP has a single top-level directory
        owd <- setwd(temp_dir)
        on.exit(setwd(owd), add = TRUE)
        zip_file <- file.path(temp_dir, paste0(stats_dir_name, ".zip"))
        # Zip the whole folder to avoid extra nested paths
        utils::zip(zip_file, files = stats_dir_name)
        file.copy(zip_file, file, overwrite = TRUE)
      }
    )

    # Expose public API
    return(list(
      pca_results = reactive({ pca_results_data() })
    ))
  })
}
