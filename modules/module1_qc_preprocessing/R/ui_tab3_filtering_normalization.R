#' @title Filtering & Normalization UI Module
#' @description UI components for Tab 3: Filtering & Normalization
#' @author Eren Ada, PhD
#' @importFrom shiny NS tagList fluidRow column h3 h4 uiOutput wellPanel div selectInput
#'             numericInput radioButtons helpText actionButton icon verbatimTextOutput
#'             plotlyOutput downloadButton DTOutput tabsetPanel tabPanel tags
#' @importFrom DT DTOutput
#' @export

# Source shared UI components
source("modules/module1_qc_preprocessing/R/shared_ui_components.R")

mod_filtering_normalization_ui <- function(id) {
  ns <- NS(id)  # Namespace function
  
  tagList(
    # --- Header Section ---
    fluidRow(
      column(width = 12,
        h3("Filtering & Normalization", 
           class = "tab-title"),
        uiOutput(ns("data_status_text")),  # Shows current data state
        tags$hr()
      )
    ),
    
    # --- Gene Filtering Section ---
    fluidRow(
      # Controls for Filtering
      column(width = 4,
        wellPanel(
          h4("Gene Filtering Settings"),
          
          # Filtering Guidance Text
          div(class = "help-text",
            h5("About Gene Filtering:"),
            p("Pre-filtering removes genes with very low counts across samples. 
              This can improve the power of downstream differential expression analysis 
              by reducing the number of tests performed (multiple testing burden). 
              Genes with consistently low counts are unlikely to be statistically significant 
              and can introduce noise."),
            h5("Common Strategies:"),
            tags$ul(
              tags$li(strong("Minimum Expression Threshold:"), " Set a minimum count/CPM/TPM 
                       that a gene must reach in a minimum number of samples (or samples within a group). 
                       For example, keep genes with > 10 counts in at least 'n' samples, where 'n' 
                       could be the size of your smallest experimental group."),
              tags$li(strong("Group-based Filtering (using Metadata):"), " If you have experimental groups 
                       (e.g., 'Treatment', 'Control'), consider filtering based on expression within 
                       at least one group. This helps retain genes that might be specifically active 
                       in certain conditions."),
              tags$li(strong("DESeq2's Independent Filtering:"), " Note that DESeq2 automatically 
                       performs independent filtering based on mean normalized counts when you call 
                       the results() function. This is generally recommended.")
            ),
            h5("Considerations:"),
            tags$ul(
              tags$li("Avoid overly aggressive filtering, which might remove biologically relevant genes."),
              tags$li("The ideal filtering criteria can depend on sequencing depth, number of samples, 
                       and overall data quality."),
              tags$li("Explore your data (e.g., using PCA plots) before and after filtering to assess the impact.")
            ),
            p("Refer to the DESeq2 vignette section on 'Pre-filtering' and 'Independent filtering' 
              for more details.")
          ),
          hr(), # Add a horizontal line for separation
          div(class = "filtering-controls",
            selectInput(ns("filter_group_column"),
              "Metadata Grouping Column:",
              choices = c("None"),  # Will be dynamically updated
              selected = "None"
            ),
            numericInput(ns("filter_min_expr_threshold"),
              "Minimum Expression Threshold:",
              value = 10,
              min = 0
            ),
            radioButtons(ns("filter_expr_unit"),
              "Expression Unit for Threshold:",
              choices = c(
                "Raw Counts" = "raw",
                "CPM" = "cpm"
              ),
              selected = "raw"
            ),
            numericInput(ns("filter_min_samples"),
              "Minimum Samples Meeting Threshold:",
              value = 2,
              min = 1
            ),
            helpText(
              "If grouping, applies to min samples *within at least one group*.",
              "If no grouping, applies to total samples."
            ),
            div(class = "action-buttons",
              actionButton(ns("apply_filters_btn"),
                "Apply Filters",
                class = "btn btn-primary btn-lg",
                icon = icon("filter")
              ),
              actionButton(ns("reset_filters_btn"),
                "Reset",
                class = "btn btn-default btn-lg",
                icon = icon("undo")
              )
            )
          )
        )
      ),
      
      # Feedback for Filtering
      column(width = 8,
        wellPanel(
          h4("Filtering Summary"),
          verbatimTextOutput(ns("filtering_summary_text")),
          
          # Plots Section with Better Alignment
          div(class = "filtering-plots-container",
            style = "margin-top: 15px;",
            fluidRow(
              column(width = 6,
                div(class = "plot-container",
                  style = "padding-right: 10px;",
                  h5("Gene Counts Before and After Filtering", 
                     style = "text-align: center; margin-bottom: 10px; font-weight: bold;"),
                  plotOutput(ns("filtering_gene_counts_plot"),
                            height = "350px")
                )
              ),
              column(width = 6,
                div(class = "plot-container",
                  style = "padding-left: 10px;",
                  h5("Expression Distribution Before vs After Filtering", 
                     style = "text-align: center; margin-bottom: 10px; font-weight: bold;"),
                  plotOutput(ns("filtering_expression_dist_plot"),
                            height = "350px")
                )
              )
            )
          ),
          
          # Data Preview Section
          div(style = "margin-top: 20px;",
            h5("Filtered Data Preview"),
            DTOutput(ns("filtered_data_preview_dt"))
          )
        )
      )
    ),
    
    # --- Normalization Section ---
    fluidRow(
      # Controls for Normalization
      column(width = 4,
        wellPanel(
          h4("Normalization Settings"),
          div(class = "normalization-controls",
            selectInput(ns("norm_method_select"),
              "Select Normalization Method:",
              choices = list(
                "Library Size Based" = c(
                  "Total Count (TC)" = "tc",
                  "Median of Ratios (DESeq2)" = "deseq",
                  "Upper Quartile (UQ)" = "uq",
                  "Relative Log Expression (RLE)" = "rle"
                ),
                "Distribution Transformations" = c(
                  "CPM (Counts Per Million)" = "cpm",
                  "TMM (Trimmed Mean of M-values)" = "tmm",
                  "VST (DESeq2)" = "vst",
                  "rlog (DESeq2)" = "rlog",
                  "Quantile Normalization" = "quantile"
                )
              ),
              selected = "deseq"
            ),
            helpText("DESeq2 is selected as default - it's robust for most RNA-seq datasets. You can change the method and compare results using the evaluation plots."),
            div(class = "action-buttons",
              actionButton(ns("apply_normalization_btn"),
                "Apply Normalization & Evaluate",
                class = "btn btn-success btn-lg",
                icon = icon("play")
              ),
              actionButton(ns("reset_normalization_btn"),
                "Reset",
                class = "btn btn-default btn-lg",
                icon = icon("undo")
              )
            )
          )
        )
      ),
      
      # Feedback and Evaluation for Normalization
      column(width = 8,
        wellPanel(
          h4("Normalization Evaluation"),
          tabsetPanel(id = ns("normalization_evaluation_tabs"),
            # Library Size Evaluation Tab
            tabPanel("Library Size Evaluation",
              plotOutput(ns("norm_eval_lib_size_plot"),
                        height = "400px"),
              verbatimTextOutput(ns("norm_eval_lib_size_stats_text"))
            ),
            
            # Normalization Effects Tab
            tabPanel("Normalization Effects",
              plotOutput(ns("norm_eval_effects_plot"),
                        height = "400px"),
              verbatimTextOutput(ns("norm_eval_effects_stats_text"))
            ),
            
            # Batch Effect Check Tab
            tabPanel("Batch Effect Check",
              div(class = "batch-effect-controls",
                selectInput(ns("norm_eval_batch_variable_select"),
                  "Select Batch Variable:",
                  choices = NULL  # Will be dynamically populated
                )
              ),
              plotOutput(ns("norm_eval_batch_plot"),
                        height = "400px"),
              verbatimTextOutput(ns("norm_eval_batch_stats_text"))
            ),
            
            # Normality Check Tab
            tabPanel("Normality Check",
              verbatimTextOutput(ns("post_norm_quick_normality_summary")),
              plotOutput(ns("normality_qq_plot"),
                        height = "300px")
            ),
            
            # Data Distribution Analysis Tab
            tabPanel("Data Distribution Analysis",
              fluidRow(
                column(width = 3,
                  h5("Distribution Analysis Settings"),
                  selectizeInput(ns("normality_samples_to_test"),
                    "Select Samples to Test:",
                    choices = NULL,
                    multiple = TRUE,
                    options = list(
                      placeholder = "Select samples...",
                      plugins = list("remove_button"),
                      maxItems = 10
                    )
                  ),
                  selectInput(ns("normality_test_method"),
                    "Normality Test Method:",
                    choices = c(
                      "Shapiro-Wilk" = "shapiro",
                      "Kolmogorov-Smirnov" = "ks",
                      "Anderson-Darling" = "ad"
                    ),
                    selected = "shapiro"
                  ),
                  selectInput(ns("normality_significance"),
                    "Significance Level:",
                    choices = c(
                      "0.05" = 0.05,
                      "0.01" = 0.01,
                      "0.001" = 0.001
                    ),
                    selected = 0.05
                  ),
                  actionButton(ns("update_normality"),
                    "Update Analysis",
                    class = "btn-primary",
                    icon = icon("refresh")
                  )
                ),
                column(width = 9,
                  h5("Normality Assessment Summary"),
                  uiOutput(ns("normality_recommendation_text")),
                  div(class = "normality-summary-container",
                    verbatimTextOutput(ns("normality_summary_text"))
                  )
                )
              ),
              
              # Distribution Plots
              fluidRow(
                column(width = 12,
                  tabsetPanel(
                    tabPanel("Q-Q Plots",
                      plotOutput(ns("normality_qq_plots"), height = "600px")
                    ),
                    tabPanel("Density Plots", 
                      plotOutput(ns("normality_density_plots"), height = "600px")
                    ),
                    tabPanel("Histograms",
                      plotOutput(ns("normality_hist_plots"), height = "600px")
                    ),
                    tabPanel("Detailed Results",
                      div(class = "detailed-results-container",
                        verbatimTextOutput(ns("normality_test_results"))
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    ),
    
    # --- Data Preview Section ---
    fluidRow(
      column(width = 12,
        wellPanel(
          div(class = "data-preview-header",
            h4("Current Data Preview"),
            span(class = "data-state-indicator",
                 textOutput(ns("current_data_state")))
          ),
          helpText("Displaying data after applied filtering and/or normalization steps."),
          div(class = "data-preview-controls",
            downloadButton(ns("download_current_data"),
                         "Download Current Data"),
            NULL
          ),
          DTOutput(ns("main_data_preview_dt"))
        )
      )
    ),
    
    # --- Proceed Section ---
    create_proceed_section(
      ns = ns,
      status_output_id = "filtering_normalization_completion_message",
      button_id = "proceed_to_next",
      button_text = "Continue to Download & Export",
      button_class = "btn-success",
      section_title = "Proceed to Next Step"
    ),
    
    # --- Custom CSS ---
    tags$head(
      tags$style(HTML("
        /* Tab Styling */
        .tab-title {
          color: #2c3e50;
          margin-bottom: 20px;
        }
        
        /* Control Panel Styling */
        .filtering-controls, .normalization-controls {
          padding: 10px 0;
        }
        
        .action-buttons {
          margin-top: 15px;
          display: flex;
          gap: 10px;
        }
        
        .action-buttons .btn {
          flex: 1;
        }
        
        /* Filtering Plots Container Styling */
        .filtering-plots-container {
          border: 1px solid #dee2e6;
          border-radius: 8px;
          padding: 15px;
          background-color: #fafbfc;
        }
        
        .plot-container {
          background: white;
          border-radius: 6px;
          border: 1px solid #e9ecef;
          padding: 10px;
          box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        }
        
        .plot-container h5 {
          margin-top: 0;
          margin-bottom: 15px;
          color: #495057;
          font-size: 0.95em;
          border-bottom: 1px solid #e9ecef;
          padding-bottom: 8px;
        }
        
        /* Data Preview Styling */
        .data-preview-header {
          display: flex;
          justify-content: space-between;
          align-items: center;
          margin-bottom: 15px;
        }
        
        .data-state-indicator {
          font-style: italic;
          color: #666;
          background: #f8f9fa;
          padding: 5px 10px;
          border-radius: 4px;
        }
        
        .data-preview-controls {
          margin: 10px 0;
          display: flex;
          gap: 10px;
        }
        
        /* Well Panel Enhancements */
        .well {
          border-radius: 8px;
          box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        
        /* Tab Panel Styling */
        .nav-tabs {
          border-bottom: 2px solid #dee2e6;
        }
        
        .nav-tabs > li.active > a {
          border-bottom-color: #2c3e50;
        }
        
        /* Plot Container Styling */
        .plotly {
          border: 1px solid #dee2e6;
          border-radius: 4px;
          padding: 10px;
        }
        
        /* Help Text Styling */
        .help-block {
          color: #6c757d;
          font-size: 0.9em;
          margin-top: 5px;
        }
        
        /* Normality Summary Styling */
        .normality-summary-container {
          max-height: 300px;
          overflow-y: auto;
          background-color: #f8f9fa;
          border: 1px solid #dee2e6;
          border-radius: 4px;
          padding: 10px;
          margin: 10px 0;
        }
        
        .normality-summary-container pre {
          background: transparent;
          border: none;
          margin: 0;
          padding: 0;
          font-size: 0.9em;
          line-height: 1.4;
        }
        
        /* Detailed Results Tab Styling */
        .detailed-results-container {
          max-height: 500px;
          overflow-y: auto;
          background-color: #f8f9fa;
          border: 1px solid #dee2e6;
          border-radius: 4px;
          padding: 15px;
          margin: 10px 0;
        }
        
        .detailed-results-container pre {
          background: transparent;
          border: none;
          margin: 0;
          padding: 0;
          font-size: 0.85em;
          line-height: 1.3;
        }
        
        /* Button States */
        .btn-primary:hover, .btn-success:hover, .btn-info:hover {
          opacity: 0.9;
        }
        
        /* Responsive Adjustments */
        @media (max-width: 768px) {
          .action-buttons {
            flex-direction: column;
          }
          
          .data-preview-controls {
            flex-direction: column;
          }
          
          .filtering-plots-container .col-sm-6 {
            margin-bottom: 15px;
          }
          
          .plot-container {
            margin-bottom: 10px;
          }
        }
      "))
    )
  )
} 