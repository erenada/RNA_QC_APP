# app_server.R

#' Application server logic
#' @param input Shiny input
#' @param output Shiny output
#' @param session Shiny session
#' @export
app_server <- function(input, output, session) {
  # Initialize shared reactive values
  shared_data <- reactiveValues(
    processed_counts = NULL,
    processed_metadata = NULL,
    data_processed = FALSE,
    normalized_counts = NULL,
    normalization_status_flag = FALSE,
    normalization_method_used = NULL,
    qc_completed = FALSE,
    qc_results = NULL
  )

  # Initialize QC module
  qc_results <- mod_input_validation_server("qc_input", session)

  # Observe QC results for use in other tabs
  observe({
    req(qc_results$data_processed())
    counts <- qc_results$processed_counts()
    metadata <- qc_results$processed_metadata()
    shared_data$processed_counts <- counts
    shared_data$processed_metadata <- metadata
    shared_data$data_processed <- TRUE
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

  # Proceed button
  observeEvent(input[["qc_plots_module-proceed_to_filtering"]], {
    req(shared_data$qc_completed)
    updateNavbarPage(session, "mainNav", selected = "filtering_normalization")
  })

  # Tab guards
  observeEvent(input$mainNav, {
    if (input$mainNav == "qc_plots") {
      req(shared_data$data_processed)
      message("QC Plots tab activated with processed data")
      message("Counts dimension: ", paste(dim(shared_data$processed_counts), collapse = " x "))
      message("Metadata dimension: ", paste(dim(shared_data$processed_metadata), collapse = " x "))
    } else if (input$mainNav == "filtering_normalization") {
      if (!shared_data$qc_completed) {
        updateNavbarPage(session, "mainNav", selected = "qc_plots")
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

  # About outputs
  output$system_info <- renderText({
    paste(
      paste("R Version:", R.version.string),
      paste("Platform:", Sys.info()["sysname"]),
      paste("Memory:", format(object.size(ls(envir = .GlobalEnv)), units = "MB")),
      paste("Packages Loaded:", length(.packages())),
      sep = "\n"
    )
  })

  output$session_status <- renderText({
    status_lines <- c()
    if (shared_data$data_processed) {
      status_lines <- c(status_lines, "Data loaded and validated")
      if (!is.null(shared_data$processed_counts)) {
        dims <- dim(shared_data$processed_counts)
        status_lines <- c(status_lines, paste("  ", dims[1], "genes,", dims[2], "samples"))
      }
    } else {
      status_lines <- c(status_lines, "No data loaded")
    }
    if (shared_data$qc_completed) {
      status_lines <- c(status_lines, "QC analysis completed")
    } else {
      status_lines <- c(status_lines, "QC analysis pending")
    }
    if (shared_data$normalization_status_flag) {
      status_lines <- c(status_lines, paste("Normalized with", shared_data$normalization_method_used))
    } else {
      status_lines <- c(status_lines, "Normalization pending")
    }
    paste(status_lines, collapse = "\n")
  })
}



