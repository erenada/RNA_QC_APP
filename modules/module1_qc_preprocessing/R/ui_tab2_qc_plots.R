# ui_tab2_qc_plots.R
# QC & Pre-processing Tool
# Tab 2: Quality Control Plots & Summaries
# Author: Eren Ada, PhD
# Date: 05/13/2024

#' @import shiny
#' @import plotly
#' @import DT
#' @import ggrepel
#' @importFrom stats cor
NULL

# Source shared UI components
source("modules/module1_qc_preprocessing/R/shared_ui_components.R")

#' @title UI Module for Quality Control Plots & Summaries tab
#' @description Creates the UI elements for the QC plots tab
#' @param id Module ID
#' @return A shiny UI module
#' @export
mod_qc_plots_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    # Data Status Indicator
    fluidRow(
      column(width = 12,
        div(class = "data-status",
          uiOutput(ns("data_status_text"))
        )
      )
    ),
    
    # Main Layout with Tabs
    tabsetPanel(
      # Tab 1: Basic QC Metrics
      tabPanel("Basic QC Metrics",
        fluidRow(
          # Controls
          column(width = 3,
            wellPanel(
              h4("Plot Controls", class = "section-title"),
              
              # Library Size Distribution Controls
              h5("Library Size Distribution"),
              selectInput(ns("lib_size_plot_type"), 
                         "Plot Type:", 
                         choices = c("Boxplot" = "boxplot", 
                                   "Barplot" = "barplot"),
                         selected = "barplot"),
              checkboxInput(ns("lib_size_log_scale"), 
                           "Log Scale Y-axis", 
                           value = TRUE),
              
              # Gene Detection Controls
              tags$hr(),
              h5("Gene Detection Rates"),
              numericInput(ns("gene_detection_threshold"),
                          "Detection Threshold (Min Count):",
                          value = 1,
                          min = 0)
            )
          ),
          
          # Plots
          column(width = 9,
            h4("General Data Overview", class = "section-title"),
            fluidRow(
              column(width = 6,
                h5("Library Size Distribution"),
                plotlyOutput(ns("library_size_plot"),
                            height = "400px")
              ),
              column(width = 6,
                h5("Gene Detection Rates"),
                plotlyOutput(ns("gene_detection_plot"),
                            height = "400px")
              )
            ),
            
            # Summary Statistics Table
            tags$hr(),
            h4("Summary Statistics", class = "section-title"),
            fluidRow(
              column(width = 12,
                DT::dataTableOutput(ns("summary_stats_table"))
              )
            )
          )
        )
      ),
      
      # Tab 2: Sample Similarity Analysis
      tabPanel("Sample Similarity Analysis",
        fluidRow(
          # Controls
          column(width = 3,
            wellPanel(
              h4("PCA Settings", class = "section-title"),
              selectInput(ns("pca_pc_x"),
                         "X-axis PC:",
                         choices = paste0("PC", 1:10),
                         selected = "PC1"),
              selectInput(ns("pca_pc_y"),
                         "Y-axis PC:",
                         choices = paste0("PC", 1:10),
                         selected = "PC2"),
              uiOutput(ns("pca_color_by_ui")),
              checkboxInput(ns("pca_show_labels"),
                           "Show Sample Labels",
                           value = FALSE)
            )
          ),
          
          # PCA Plots
          column(width = 9,
            h4("Principal Component Analysis", class = "section-title"),
            fluidRow(
              column(width = 6,
                h5("2D PCA Plot"),
                plotlyOutput(ns("pca_plot_2d"),
                            height = "500px")
              ),
              column(width = 6,
                h5("3D PCA Plot"),
                plotlyOutput(ns("pca_plot_3d"),
                            height = "500px")
              )
            ),
            fluidRow(
              column(width = 12,
                div(style = "display: flex; justify-content: space-between; align-items: center;",
                  uiOutput(ns("pca_variance_explained_text")),
                  downloadButton(ns("download_pca_stats"), "Download PCA Statistics (2D & 3D)", 
                               class = "btn-info btn-sm")
                )
              )
            ),
            # PCA Downloads Section
            tags$hr(),
            fluidRow(
              column(width = 12,
                div(class = "text-right",
                  downloadButton(ns("download_pca_plot"),
                               "Download PCA Plot",
                               class = "btn-info btn-sm")
                )
              )
            )
          )
        )
      ),
      
      # Tab 3: Sample Correlation Analysis
      tabPanel("Sample Correlation Analysis",
        fluidRow(
          # Controls
          column(width = 3,
            wellPanel(
              h4("Correlation Settings", class = "section-title"),
              # Data source selection
              selectInput(ns("correlation_data_source"),
                         "Data Source:",
                         choices = c(
                           "log2(CPM+1) transformed data" = "log2_cpm",
                           "Raw count matrix" = "raw_counts"
                         ),
                         selected = "log2_cpm"),
              tags$div(
                style = "font-size: 12px; color: #d32f2f; margin-top: 5px; margin-bottom: 10px; border: 1px solid #d32f2f; background-color: #ffeaa7; padding: 8px; border-radius: 4px;",
                tags$strong("Warning:"), 
                tags$br(),
                tags$span("If using raw count matrix, ensure you select the appropriate correlation method (typically Spearman). Consider applying normalization in the Filtering & Normalization tab first.")
              ),
              selectInput(ns("correlation_method"),
                         "Correlation Method:",
                         choices = c("Pearson" = "pearson",
                                   "Spearman" = "spearman"),
                         selected = "pearson"),
              tags$div(
                style = "font-size: 12px; color: #666; margin-top: 5px; margin-bottom: 10px;",
                tags$strong("Guidance:"), 
                tags$br(),
                tags$span("• Pearson: Recommended for log-transformed data"),
                tags$br(),
                tags$span("• Spearman: Use if you suspect outliers or prefer rank-based correlation")
              ),
              selectInput(ns("cor_color_scheme"),
                         "Color Scheme:",
                         choices = c(
                           "Default" = "default",
                           "Blue-Red" = "blue_red",
                           "Purple-Orange" = "purple_orange",
                           "Ocean" = "ocean"
                         ),
                         selected = "default"),
              tags$hr(),
              h4("Sample Filtering", class = "section-title"),
              # Dynamic metadata filtering section
              uiOutput(ns("filter_criteria")),
              uiOutput(ns("dynamic_filters")),
              tags$hr(),
              h4("Display Options", class = "section-title"),
              checkboxInput(ns("show_cor_values"), 
                          "Show Correlation Values", 
                          value = FALSE),
              checkboxInput(ns("show_significance"), 
                          "Show Significance Stars", 
                          value = TRUE),
              checkboxInput(ns("auto_scale_colors"), 
                          "Auto-scale Colors for High Correlations", 
                          value = TRUE),
              actionButton(ns("update_correlation"), 
                         "Update Correlation Plot",
                         class = "btn-primary")
            )
          ),
          
          # Correlation Plot and Stats
          column(width = 9,
            h4("Sample-to-Sample Correlation", class = "section-title"),
            plotOutput(ns("correlation_heatmap"),
                      height = "600px"),
            uiOutput(ns("correlation_stats_text")),
            # Correlation Downloads Section
            tags$hr(),
            fluidRow(
              column(width = 12,
                div(class = "text-right",
                  downloadButton(ns("download_heatmap"),
                               "Download Heatmap",
                               class = "btn-info btn-sm")
                )
              )
            )
          )
        )
      ),
      
      # Tab 4: Normality Assessment has been moved to Filtering & Normalization tab
    ),
    
    # Add proceed section at the bottom
    create_proceed_section(
      ns = ns,
      status_output_id = "qc_completion_message",
      button_id = "proceed_to_filtering",
      button_text = "Continue to Filtering & Normalization",
      button_class = "btn-success",
      section_title = "Proceed to Next Step"
    )
  )
}

# Example server call (commented out as it belongs in a different file)
# mod_qc_plots_server("qc_plots_module", shared_data_reactive_values) 