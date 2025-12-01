# ui_tab1.R
# Module 1: QC & Pre-processing Tool
# Tab 1: Data Input & Validation
# Author: Eren Ada, PhD
# Date: 05/13/2024

#' @import shiny
#' @import DT
#' @import shinydashboard
NULL

# Source shared UI components
source("modules/module1_qc_preprocessing/R/shared_ui_components.R")

library(shiny)
library(DT)

#' UI Module for Data Input & Validation tab
#' @param id Module ID
#' @importFrom shiny NS tagList fluidRow column wellPanel h4 div tags icon fileInput actionButton conditionalPanel selectInput helpText downloadButton DTOutput uiOutput
#' @importFrom DT renderDT datatable
mod_input_validation_ui <- function(id) {
  # Create a namespace function specific to this instance
  ns <- shiny::NS(id)
  
  tagList(
    fluidRow(
      # Column 1: File Inputs and Initial Analysis
      column(width = 4,
        wellPanel(
          # Section 1: File Uploads
          h4("1. Upload Data Files", class = "section-title"),
          
          # Info about supported formats
          div(class = "format-info",
            tags$p(
              icon("info-circle"), 
              "Supported formats: CSV (.csv) and TSV (.tsv, .txt)",
              tags$br(),
              "File format and separator will be automatically detected."
            )
          ),
          
          # Extended information about file requirements
          div(class = "file-requirements",
            h5("Count Matrix Requirements:", class = "mt-3"),
            tags$ul(
              tags$li("First column: Gene IDs/names"),
              tags$li("Column headers: Sample names"),
              tags$li("Values: Raw counts (non-negative integers)"),
              tags$li("No missing values allowed")
            ),
            
            h5("Metadata Requirements:", class = "mt-3"),
            tags$ul(
              tags$li("First column: Sample names (matching count matrix)"),
              tags$li("Additional columns: Sample attributes (e.g., treatment, condition)"),
              tags$li("Each sample must have complete information"),
              tags$li("Sample names must be unique")
            ),
            
            # Expandable example format section
            tags$div(class = "example-format-section",
              tags$button(
                "Show Example Format",
                onclick = "$(this).next().slideToggle(); $(this).text($(this).text() == 'Show Example Format' ? 'Hide Example Format' : 'Show Example Format');",
                class = "btn btn-sm btn-info mt-3"
              ),
              tags$div(class = "file-example mt-2", style = "display: none;",
                h5("Example Format:"),
                tags$p("Count Matrix:"),
                tags$pre(
                  "Gene_ID,Sample1,Sample2,Sample3\n",
                  "GENE1,123,456,789\n",
                  "GENE2,234,567,890\n",
                  "..."
                ),
                tags$p("Metadata:"),
                tags$pre(
                  "Sample,Treatment,Condition\n",
                  "Sample1,Control,WT\n",
                  "Sample2,Treatment,KO\n",
                  "Sample3,Control,KO\n",
                  "..."
                )
              )
            )
          ),
          
          fileInput(ns("count_matrix_file"),
                   "Upload Count Matrix",
                   multiple = FALSE,
                   accept = c(".csv", ".tsv", ".txt")),
          
          fileInput(ns("metadata_file"),
                   "Upload Sample Metadata",
                   multiple = FALSE,
                   accept = c(".csv", ".tsv", ".txt")),
          
          # Initial Analysis Button
          div(style = "margin-top: 20px;",
              actionButton(ns("analyze_files_btn"),
                         "Validate Files",
                         icon = icon("check-circle"),
                         class = "btn btn-primary btn-lg btn-block"))
        ),
        
        # Section 2: Data Analysis Results (Initially Hidden)
        conditionalPanel(
          condition = sprintf("output['%s'] == true", ns("files_analyzed")),
          wellPanel(
            h4("2. Data Analysis Results", class = "section-title"),
            
            # Display detected settings
            uiOutput(ns("detected_settings")),
            
            # Add download original data button
            downloadButton(ns("download_original_data"), 
                         "Download Original Count Matrix",
                         class = "btn-info btn-sm"),
            
            # Display data issues (focusing on duplicates)
            conditionalPanel(
              condition = sprintf("output['%s'] == true", ns("has_duplicates")),
              div(class = "alert alert-warning mt-3",
                h4(icon("exclamation-triangle"), " Duplicate Genes Detected", 
                   class = "alert-heading"),
                p("Found duplicate gene IDs that need to be resolved before proceeding."),
                
                # Duplicate handling options in a clean card
                div(class = "card mt-3",
                  div(class = "card-body",
                    h5(class = "card-title", "Duplicate Gene Handling"),
                    
                    # Strategy selector
                    selectInput(ns("duplicate_strategy"),
                              "Select handling method:",
                              choices = c(
                                "Sum Counts for Duplicates" = "sum",
                                "Keep Entry with Highest Average Expression" = "keep_highest_avg",
                                "Rename with Suffix (e.g., GENE_1)" = "rename",
                                "Remove All Duplicate Entries" = "remove_all"
                              ),
                              selected = "sum"),
                    
                    # Preview of what the selected strategy will do
                    div(class = "alert alert-info",
                      h5(icon("info-circle"), "Preview of selected strategy:"),
                      p(class = "strategy-preview"),
                      uiOutput(ns("duplicate_preview")),
                      p(class = "text-muted", "This change will be applied when you click 'Process Data'")
                    ),
                    
                    # Duplicate entries preview
                    div(class = "duplicate-preview mt-3",
                      h5("Preview of Duplicated Entries:"),
                      DTOutput(ns("duplicate_genes_table")),
                      downloadButton(ns("download_duplicate_stats"), 
                                  "Download Complete Statistics",
                                  class = "btn-sm btn-info mt-2")
                    )
                  )
                )
              )
            ),
            
            # Process Data Button
            div(style = "margin-top: 20px;",
                actionButton(ns("process_data_btn"),
                           "Process Data",
                           icon = icon("play-circle"),
                           class = "btn btn-success btn-lg btn-block"))
          )
        )
      ),
      
      # Column 2: Data Preview and Validation Summary
      column(width = 8,
        wellPanel(
          h4("Data Preview & Validation Summary", class = "section-title"),
          
          # Initial file analysis results
          conditionalPanel(
            condition = sprintf("output['%s'] == true", ns("files_analyzed")),
            
            # Data Summary Section with better space utilization
            wellPanel(
              h4("Data Summary"),
              fluidRow(class = "equal-height-row",
                column(width = 8, class = "equal-height-col",
                  uiOutput(ns("data_summary"))
                ),
                column(width = 4, class = "equal-height-col",
                  conditionalPanel(
                    condition = sprintf("output['%s'] == true", ns("data_processed")),
                    h5("Processed Data Summary"),
                    uiOutput(ns("processed_data_summary"))
                  )
                )
              )
            ),
            
            # Library Size Plot
            wellPanel(
              h4("Library Size Plot"),
              plotlyOutput(ns("library_size_plot"), height = "400px")
            ),
            
            # Data Previews back to back (vertical stacking)
            h4("Count Matrix Preview"),
            DTOutput(ns("count_matrix_preview")),
            
            h4("Metadata Preview"),
            DTOutput(ns("metadata_preview"))
          ),
          
          # Final processing results
          conditionalPanel(
            condition = sprintf("output['%s'] == true", ns("data_processed")),
            
            # Validation Summary
            h4("Processing Summary"),
            div(
              style = "margin: 15px 0;",
              downloadButton(ns("download_processed_matrix"), "Download Processed Count Matrix",
                           class = "btn-info")
            ),
            uiOutput(ns("validation_summary")),
            
            # Next Step Button
            create_proceed_section(
              ns = ns,
              status_output_id = "processing_completion_message",
              button_id = "next_step_btn", 
              button_text = "Continue to QC Analysis",
              button_class = "btn-success",
              section_title = "Proceed to Next Step"
            )
          ),
          
          # Message when no data is loaded
          conditionalPanel(
            condition = sprintf("!output['%s']", ns("files_analyzed")),
            div(
              class = "alert alert-info",
              icon("info-circle"),
              "Please upload your files and click 'Analyze Files' to begin."
            )
          )
        )
      )
    )
  )
} 