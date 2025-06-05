# app.R
# QC & Pre-processing Tool
# Author: Eren Ada, PhD
# Date: 05/29/2024

# Load essential packages first
suppressPackageStartupMessages({
  # Core Bioconductor infrastructure
  library(methods)
  library(BiocGenerics)
  library(S4Vectors)
  library(IRanges)
  library(GenomeInfoDb)
  library(GenomicRanges)
  library(SummarizedExperiment)
  
  # Fix for DESeq2 on R 4.5.0 - "superclass 'ExpData' not defined" error
  setOldClass("ExpData")
  
  library(DESeq2)
})

# Function to check and install Bioconductor packages
check_and_install_bioc_packages <- function(packages) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) {
    message("Installing missing Bioconductor packages: ", paste(new_packages, collapse = ", "))
    BiocManager::install(new_packages)
  }
}

# Function to check and install CRAN packages
check_and_install_packages <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) {
    message("Installing missing CRAN packages: ", paste(new_packages, collapse = ", "))
    install.packages(new_packages)
  }
}

# Required CRAN packages
required_packages <- c(
  # Shiny and UI packages
  "shiny",
  "shinythemes",
  "shinyWidgets",
  "DT",
  "htmltools",
  
  # Data manipulation packages
  "dplyr",
  "tidyr",
  "readr",
  
  # Visualization packages
  "ggplot2",
  "plotly",
  "scales",
  "pheatmap",
  "ggrepel",
  "viridis",
  "RColorBrewer",
  "gridExtra",  # Added for plot arrangements
  "cowplot",    # Added for plot theme modifications
  "grid",       # Added for pheatmap rendering in Shiny
  "corrplot",   # Added for more reliable correlation plotting in Shiny
  
  # Statistical packages
  "moments",
  "e1071",  # Added for skewness and kurtosis calculations
  "matrixStats",  # Added for row/column statistics
  
  # Additional packages for filtering & normalization
  "preprocessCore"
)

# Required Bioconductor packages
required_bioc_packages <- c(
  # Core Bioconductor infrastructure
  "BiocGenerics",
  "S4Vectors",
  "IRanges",
  "GenomeInfoDb",
  "GenomicRanges",
  "SummarizedExperiment",
  
  # RNA-seq specific packages
  "DESeq2",
  "edgeR",
  "limma"  # Added for normalization methods
)

# Check and install required packages
check_and_install_packages(required_packages)
check_and_install_bioc_packages(required_bioc_packages)

# Set up DESeq2 class definitions
message("\nSetting up DESeq2 class definitions...")
tryCatch({
  if (!isClass("ExpData")) {
    setClass("ExpData",
             slots = c(
               assays = "SimpleList",
               colData = "DataFrame",
               metadata = "list"
             ),
             contains = "VIRTUAL"
    )
  }
  
  # Load remaining required packages
  suppressPackageStartupMessages({
    # Load CRAN packages in groups
    # UI packages
    library(shiny)
    library(shinythemes)
    library(shinyWidgets)
    library(DT)
    library(htmltools)
    
    # Data manipulation
    library(dplyr)
    library(tidyr)
    library(readr)
    
    # Visualization
    library(ggplot2)
    library(plotly)
    library(scales)
    library(pheatmap)
    library(ggrepel)
    library(viridis)
    library(RColorBrewer)
    library(grid)  # Added for pheatmap rendering in Shiny
    library(corrplot)  # Added for more reliable correlation plotting in Shiny
    
    # Statistical
    library(moments)
    library(e1071)
    
    # Additional
    library(preprocessCore)
    library(methods)  # Added to ensure method dispatch works properly
    
    # Load remaining Bioconductor packages
    library(edgeR)
    library(matrixStats)
  })
  
  message("Package and class setup completed successfully")
}, error = function(e) {
  message("Error in setup: ", e$message)
  stop(e)
})

# Configure Shiny options
options(shiny.maxRequestSize = 1000 * 1024^2)  # Set max upload size to 1000MB (1GB)

# Source module files
source("modules/module1_qc_preprocessing/R/ui_tab1.R")
source("modules/module1_qc_preprocessing/R/server_tab1.R")
source("modules/module1_qc_preprocessing/R/ui_tab2_qc_plots.R")
source("modules/module1_qc_preprocessing/R/server_tab2_qc_plots.R")
source("modules/module1_qc_preprocessing/R/utils_qc_plotting.R")

# Source new filtering and normalization modules
source("modules/module1_qc_preprocessing/R/ui_tab3_filtering_normalization.R")
source("modules/module1_qc_preprocessing/R/server_tab3_filtering_normalization.R")
source("modules/module1_qc_preprocessing/R/utils_filtering.R")
source("modules/module1_qc_preprocessing/R/utils_normalization.R")
source("modules/module1_qc_preprocessing/R/utils_plotting.R")

# Define UI
ui <- fluidPage(
  theme = shinytheme("flatly"),
  
  # Include custom CSS
  includeCSS("modules/module1_qc_preprocessing/www/styles.css"),
  
  # App title
  titlePanel("QC & Pre-processing Tool"),
  
  # Main navigation structure
  navbarPage(
    title = NULL,
    id = "mainNav",
    
    # Data Input & Validation Tab
    tabPanel(
      "Data Input & Validation",
      value = "data_input",
      mod_input_validation_ui("qc_input")
    ),
    
    # QC Plots Tab
    tabPanel(
      "QC Plots & Summaries",
      value = "qc_plots",
      mod_qc_plots_ui("qc_plots_module")
    ),
    
    # Filtering & Normalization Tab
    tabPanel(
      "Filtering & Normalization",
      value = "filtering_normalization",
      mod_filtering_normalization_ui("filtering_normalization_module")
    ),
    
    # About tab
    tabPanel("About",
             value = "about",
             fluidRow(
               column(8,
                      h3("Modular Bulk RNA-seq Analysis Tool Suite"),
                      h4("QC & Pre-processing Module", style = "color: #666; margin-top: -10px;"),
                      
                      div(style = "background-color: #e8f4fd; padding: 15px; border-radius: 8px; border-left: 4px solid #1f77b4; margin: 15px 0;",
                          p(strong("Version:"), "2.0.0 |", 
                            strong("Author:"), "Eren Ada, PhD |",
                            strong("Institution:"), "Harvard Medical School, Department of Immunology |",
                            strong("Date:"), "6/5/2025"),
                          p(strong("License:"), "CC BY-NC-SA 4.0 (Non-commercial research use) |",
                            strong("Repository:"), a("GitHub", href = "https://github.com/hms-immunology/RNA_QC_APP", target = "_blank"))
                      ),
                      
                      h4("About This Tool"),
                      p("This application is part of a comprehensive, modular toolkit designed for bulk RNA-sequencing analysis.",
                        "The QC & Pre-processing module provides rigorous quality control assessment and data preparation",
                        "capabilities, ensuring your RNA-seq count data meets the standards required for robust downstream analysis."),
                      
                      p("Built specifically for researchers conducting bulk RNA-seq experiments, this tool bridges the gap between",
                        "raw count matrices and analysis-ready datasets. It implements best practices from the bioinformatics",
                        "community and provides intuitive visualizations to guide data-driven decisions throughout the pre-processing workflow."),
                      
                      div(style = "background-color: #fff3cd; padding: 12px; border-radius: 5px; border-left: 3px solid #ffc107; margin: 15px 0;",
                                                     p(strong("Perfect for:"), 
                            "Biologists, bioinformaticians, and researchers working with bulk RNA-seq data who need",
                            "comprehensive QC assessment and standardized pre-processing before differential expression analysis.",
                            style = "margin: 0;")),
                      
                      hr(),
                      
                      h4("Quick Start"),
                      p("New to this tool?", 
                        a("Follow our 15-minute Quick Start Guide", 
                          href = "docs/user_manual/quick_start.html", target = "_blank"),
                        "to get up and running with example data."),
                      
                                             h4("Core Capabilities"),
                      div(
                        style = "display: flex; flex-wrap: wrap; gap: 20px;",
                        div(
                          style = "flex: 1; min-width: 220px; background-color: #f8f9fa; padding: 15px; border-radius: 8px;",
                          h5("Data Import & Validation", style = "color: #495057; margin-top: 0;"),
                          tags$ul(
                            tags$li("Flexible CSV import with auto-detection"),
                            tags$li("Comprehensive data validation & QC checks"),
                            tags$li("Duplicate gene detection & resolution"),
                            tags$li("Sample-metadata consistency verification"),
                            tags$li("Missing data identification & handling")
                          )
                        ),
                        div(
                          style = "flex: 1; min-width: 220px; background-color: #f8f9fa; padding: 15px; border-radius: 8px;",
                          h5("Quality Control Analysis", style = "color: #495057; margin-top: 0;"),
                          tags$ul(
                            tags$li("Library size distribution analysis"),
                            tags$li("Gene expression distribution profiling"),
                            tags$li("Sample-to-sample correlation heatmaps"),
                            tags$li("Interactive 2D & 3D PCA visualization"),
                            tags$li("Outlier detection & flagging"),
                            tags$li("Comprehensive QC summary reports")
                          )
                        ),
                        div(
                          style = "flex: 1; min-width: 220px; background-color: #f8f9fa; padding: 15px; border-radius: 8px;",
                          h5("Pre-processing & Export", style = "color: #495057; margin-top: 0;"),
                          tags$ul(
                            tags$li("Multiple gene filtering strategies"),
                            tags$li("Various normalization methods (CPM, TMM, VST, etc.)"),
                            tags$li("Before/after comparison visualizations"),
                            tags$li("Publication-ready plot exports"),
                            tags$li("Processed data export (CSV, RDS)"),
                            tags$li("Comprehensive analysis reports")
                          )
                        )
                      ),
                      
                      hr(),
                      
                                             h4("Analysis Workflow"),
                      div(
                        style = "background-color: #e8f5e8; padding: 20px; border-radius: 8px; margin: 15px 0; border-left: 4px solid #28a745;",
                        div(style = "text-align: center; margin-bottom: 15px;",
                            h5("Step-by-Step Process", style = "margin: 0; color: #155724;")),
                                                 p(style = "margin: 0; font-family: 'Courier New', monospace; text-align: center; font-size: 14px; color: #155724;",
                           "Data Input → Validation → QC Analysis → Filtering → Normalization → Export"),
                        hr(style = "margin: 15px 0; border-color: #c3e6cb;"),
                        div(style = "font-size: 13px; color: #155724;",
                            p(strong("Typical Session:"), "15-30 minutes | ",
                              strong("Data Size:"), "Up to 50K genes × 200 samples | ",
                              strong("Output:"), "Analysis-ready datasets + QC reports", 
                              style = "margin: 0;"))
                      ),
                      
                                             h4("Use Cases & Applications"),
                      div(style = "background-color: #fff; border: 1px solid #dee2e6; border-radius: 8px; padding: 15px; margin: 15px 0;",
                          tags$ul(style = "margin-bottom: 0;",
                            tags$li(strong("Experimental Design Validation:"), "Assess sample groupings and identify batch effects"),
                            tags$li(strong("Data Quality Assessment:"), "Evaluate library preparation quality and sequencing depth"),
                            tags$li(strong("Outlier Detection:"), "Identify problematic samples before downstream analysis"),
                            tags$li(strong("Normalization Strategy:"), "Choose appropriate normalization methods based on data characteristics"),
                            tags$li(strong("Publication Preparation:"), "Generate high-quality QC plots for methods sections")
                          )),
                      
                      hr(),
                      
                                             h4("Documentation & Resources"),
                      div(style = "display: flex; flex-wrap: wrap; gap: 15px; margin: 15px 0;",
                          div(style = "flex: 1; min-width: 200px;",
                              h5("User Guides", style = "color: #495057; margin-bottom: 10px;"),
                              tags$ul(style = "padding-left: 20px;",
                                tags$li(a("User Manual Overview", 
                                         href = "docs/user_manual/README.html", target = "_blank")),
                                tags$li(a("Quick Start Guide (15 min)", 
                                         href = "docs/user_manual/quick_start.html", target = "_blank")),
                                tags$li(a("Data Input & Validation", 
                                         href = "docs/user_manual/data_input.html", target = "_blank")),
                                tags$li(a("QC Plots Interpretation", 
                                         href = "docs/user_manual/qc_plots_interpretation_guide.html", target = "_blank"))
                              )),
                          div(style = "flex: 1; min-width: 200px;",
                              h5("Advanced Topics", style = "color: #495057; margin-bottom: 10px;"),
                              tags$ul(style = "padding-left: 20px;",
                                tags$li(a("Filtering & Normalization", 
                                         href = "docs/user_manual/filtering_normalization_guide.html", target = "_blank")),
                                tags$li(a("Export & Download Options", 
                                         href = "docs/user_manual/export_download_guide.html", target = "_blank")),
                                tags$li(a("Technical Requirements", 
                                         href = "docs/user_manual/technical_requirements.html", target = "_blank")),
                                tags$li(a("Troubleshooting Guide", 
                                         href = "docs/user_manual/troubleshooting.html", target = "_blank"))
                              ))
                      ),
                      
                      div(style = "background-color: #d1ecf1; padding: 15px; border-radius: 8px; border-left: 4px solid #17a2b8; margin: 15px 0;",
                          h5("Getting Started", style = "margin-top: 0; color: #0c5460;"),
                          p("New users should start with the", 
                            a("Quick Start Guide", href = "docs/user_manual/quick_start.html", target = "_blank", style = "font-weight: bold;"), 
                            "which walks through a complete analysis using provided example data. For specific questions,", 
                            "consult the detailed user manual sections above.",
                            style = "margin-bottom: 0; color: #0c5460;"))
                      
               ),
               column(4,
                      wellPanel(
                        h4("System Information"),
                        verbatimTextOutput("system_info"),
                        hr(),
                        h4("Session Status"),
                        verbatimTextOutput("session_status"),
                        hr(),
                        div(style = "text-align: center; margin-top: 15px;",
                            h5("Quick Links", style = "margin-bottom: 10px;"),
                            p(
                              a("Documentation", href = "docs/user_manual/README.html", target = "_blank", 
                                class = "btn btn-info btn-sm", style = "margin: 2px;"), br(),
                              a("Report Issue", href = "https://github.com/hms-immunology/RNA_QC_APP/issues", target = "_blank", 
                                class = "btn btn-warning btn-sm", style = "margin: 2px;"), br(),
                              a("Contact Author", href = "mailto:erenada@gmail.com", 
                                class = "btn btn-secondary btn-sm", style = "margin: 2px;"),
                              style = "margin: 0;"
                            ))
                      )
               )
             ))
  )
)

# Define server
server <- function(input, output, session) {
  # Initialize shared reactive values
  shared_data <- reactiveValues(
    processed_counts = NULL,
    processed_metadata = NULL,
    data_processed = FALSE,
    normalized_counts = NULL,
    normalization_status_flag = FALSE,
    normalization_method_used = NULL,
    qc_completed = FALSE,  # New flag for QC completion
    qc_results = NULL      # Store QC results
  )
  
  # Initialize QC module
  qc_results <- mod_input_validation_server("qc_input", session)
  
  # Observe QC results for use in other tabs
  observe({
    req(qc_results$data_processed())
    
    # Access processed data
    counts <- qc_results$processed_counts()
    metadata <- qc_results$processed_metadata()
    
    # Store in shared values for use across tabs
    shared_data$processed_counts <- counts
    shared_data$processed_metadata <- metadata
    shared_data$data_processed <- TRUE
    
    # Print summary to console for debugging
    message("Data processed successfully")
    message("Counts dimension: ", paste(dim(counts), collapse = " x "))
    message("Metadata dimension: ", paste(dim(metadata), collapse = " x "))
  })
  
  # Initialize QC Plots module when data is ready
  qc_plots_results <- mod_qc_plots_server("qc_plots_module", shared_data)
  
  # Observe QC completion status
  observe({
    req(qc_plots_results)
    shared_data$qc_completed <- qc_plots_results$qc_completion_status()
    if (!is.null(qc_plots_results$qc_results)) {
      shared_data$qc_results <- qc_plots_results$qc_results
    }
  })
  
  # Initialize Filtering & Normalization module when data is ready
  observe({
    req(shared_data$data_processed)
    mod_filtering_normalization_server("filtering_normalization_module", shared_data)
  })
  
  # Observer for the proceed button in QC Plots tab
  observeEvent(input[["qc_plots_module-proceed_to_filtering"]], {
    req(shared_data$qc_completed)
    updateTabsetPanel(session, "mainNav", selected = "filtering_normalization")
  })
  
  # Observer for tab switching with validation
  observeEvent(input$mainNav, {
    if(input$mainNav == "qc_plots") {
      req(shared_data$data_processed)
      message("QC Plots tab activated with processed data")
      message("Counts dimension: ", paste(dim(shared_data$processed_counts), collapse = " x "))
      message("Metadata dimension: ", paste(dim(shared_data$processed_metadata), collapse = " x "))
    }
    else if(input$mainNav == "filtering_normalization") {
      # Validate that QC is completed before allowing access
      if (!shared_data$qc_completed) {
        showNotification(
          "Please complete the QC analysis before proceeding to filtering and normalization.",
          type = "warning",
          duration = 5
        )
        updateTabsetPanel(session, "mainNav", selected = "qc_plots")
      } else {
        message("Filtering & Normalization tab activated with processed data")
        message("Counts dimension: ", paste(dim(shared_data$processed_counts), collapse = " x "))
        message("Metadata dimension: ", paste(dim(shared_data$processed_metadata), collapse = " x "))
        if (!is.null(shared_data$qc_results)) {
          message("QC results available for use in filtering & normalization")
        }
      }
    }
  })
  
  # System information output for About tab
  output$system_info <- renderText({
    paste(
      paste("R Version:", R.version.string),
      paste("Platform:", Sys.info()["sysname"]),
      paste("Memory:", format(object.size(ls(envir = .GlobalEnv)), units = "MB")),
      paste("Packages Loaded:", length(.packages())),
      sep = "\n"
    )
  })
  
  # Session status output for About tab
  output$session_status <- renderText({
    status_lines <- c()
    
    if (shared_data$data_processed) {
      status_lines <- c(status_lines, "✓ Data loaded and validated")
      if (!is.null(shared_data$processed_counts)) {
        dims <- dim(shared_data$processed_counts)
        status_lines <- c(status_lines, paste("  ", dims[1], "genes,", dims[2], "samples"))
      }
    } else {
      status_lines <- c(status_lines, "○ No data loaded")
    }
    
    if (shared_data$qc_completed) {
      status_lines <- c(status_lines, "✓ QC analysis completed")
    } else {
      status_lines <- c(status_lines, "○ QC analysis pending")
    }
    
    if (shared_data$normalization_status_flag) {
      status_lines <- c(status_lines, paste("✓ Normalized with", shared_data$normalization_method_used))
    } else {
      status_lines <- c(status_lines, "○ Normalization pending")
    }
    
    paste(status_lines, collapse = "\n")
  })
}

# Run the application
shinyApp(ui = ui, server = server) 