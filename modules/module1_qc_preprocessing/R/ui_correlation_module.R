# ui_correlation_module.R
# Correlation Submodule UI

#' @title UI for Correlation Submodule
#' @param id Module ID
#' @return Shiny UI elements for correlation analysis
#' @export
mod_correlation_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(width = 3,
        wellPanel(
          h4("Correlation Settings", class = "section-title"),
          checkboxInput(ns("compute_p_values"),
                        "Compute p-values (slower)",
                        value = TRUE),
          helpText("Computes pairwise correlations; enabling p-values can be slow for many samples."),
          selectInput(ns("correlation_method"),
                      "Correlation Method:",
                      choices = c("Pearson" = "pearson", "Spearman" = "spearman"),
                      selected = "pearson"),
          tags$hr(),
          h4("Sample Filtering", class = "section-title"),
          uiOutput(ns("filter_criteria")),
          uiOutput(ns("dynamic_filters")),
          tags$hr(),
          h4("Display Options", class = "section-title"),
          selectInput(ns("cor_color_scheme"),
                      "Color Scheme:",
                      choices = c("Default" = "default","Blue-Red" = "blue_red","Purple-Orange" = "purple_orange","Ocean" = "ocean"),
                      selected = "default"),
          checkboxInput(ns("show_cor_values"), "Show Correlation Values", value = FALSE),
          checkboxInput(ns("show_significance"), "Show Significance Stars", value = TRUE),
          checkboxInput(ns("auto_scale_colors"), "Auto-scale Colors for High Correlations", value = TRUE),
          actionButton(ns("update_correlation"), "Update Correlation Plot", class = "btn btn-primary btn-lg")
        )
      ),
      column(width = 9,
        h4("Sample-to-Sample Correlation", class = "section-title"),
        plotOutput(ns("correlation_heatmap"), height = "600px"),
        uiOutput(ns("correlation_stats_text")),
        tags$hr(),
        fluidRow(
          column(width = 12,
            div(class = "text-right",
              downloadButton(ns("download_heatmap"), "Download Heatmap", class = "btn-info btn-sm")
            )
          )
        )
      )
    )
  )
}
