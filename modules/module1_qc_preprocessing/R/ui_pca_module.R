# ui_pca_module.R
# PCA Submodule UI

#' @title UI for PCA Submodule
#' @param id Module ID
#' @return Shiny UI elements for PCA
#' @export
mod_pca_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      # Controls
      column(width = 3,
        wellPanel(
          h4("PCA Settings", class = "section-title"),
          helpText("PCA is computed on transformed data (centered and scaled). Large datasets may take time."),
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
      # Plots and variance stats
      column(width = 9,
        h4("Principal Component Analysis", class = "section-title"),
        fluidRow(
          column(width = 6,
            h5("2D PCA Plot"),
            plotlyOutput(ns("pca_plot_2d"), height = "500px")
          ),
          column(width = 6,
            h5("3D PCA Plot"),
            plotlyOutput(ns("pca_plot_3d"), height = "500px")
          )
        ),
        fluidRow(
          column(width = 12,
            div(style = "display: flex; justify-content: space-between; align-items: center;",
              uiOutput(ns("pca_variance_explained_text")),
              downloadButton(ns("download_pca_stats"), "Download PCA Statistics (2D & 3D)", class = "btn-info btn-sm")
            )
          )
        ),
        tags$hr(),
        fluidRow(
          column(width = 12,
            div(class = "text-right",
              downloadButton(ns("download_pca_plot"), "Download PCA Plot", class = "btn-info btn-sm")
            )
          )
        )
      )
    )
  )
}
