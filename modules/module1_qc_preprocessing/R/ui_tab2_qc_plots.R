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
        mod_pca_ui(ns("pca"))
      ),
      
      # Tab 3: Sample Correlation Analysis
      tabPanel("Sample Correlation Analysis",
        mod_correlation_ui(ns("corr"))
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