# ui_about.R

#' Create the About tab UI
#' @return A shiny tabPanel containing the About content
#' @export
create_about_tab <- function() {
  tabPanel("About",
           value = "about",
           fluidRow(
             column(8,
                    h3("Modular Bulk RNA-seq Analysis Tool Suite"),
                    h4("QC & Pre-processing Module", style = "color: #666; margin-top: -10px;"),
                    
                    div(class = "about-box info",
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
                    
                    div(class = "about-box warning",
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
                      class = "about-box workflow",
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
                    
                    div(class = "about-box info",
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
}



