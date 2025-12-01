# app_ui.R

#' Build the application UI
#' @return Shiny UI object
#' @export
app_ui <- function() {
  fluidPage(
    theme = shinytheme("flatly"),
    includeCSS("modules/module1_qc_preprocessing/www/styles.css"),
    titlePanel("QC & Pre-processing Tool"),
    navbarPage(
      title = NULL,
      id = "mainNav",
      tabPanel("Data Input & Validation", value = "data_input", mod_input_validation_ui("qc_input")),
      tabPanel("QC Plots & Summaries", value = "qc_plots", mod_qc_plots_ui("qc_plots_module")),
      tabPanel("Filtering & Normalization", value = "filtering_normalization", mod_filtering_normalization_ui("filtering_normalization_module")),
      create_about_tab()
    )
  )
}



