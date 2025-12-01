#' @title Shared UI Components for QC Module
#' @description Reusable UI components for consistent styling across tabs
#' @author Eren Ada, PhD
#' @date 05/29/2025

#' @import shiny
#' @import shinydashboard
NULL

#' @title Create Standardized Proceed Section
#' @description Creates a consistent proceed section with status message and button
#' @param ns Namespace function for the module
#' @param status_output_id The output ID for the status message
#' @param button_id The button ID for the proceed button
#' @param button_text The text to display on the button
#' @param button_icon The icon to use (default: "arrow-right")
#' @param button_class The CSS class for the button (default: "btn-success")
#' @param section_title The title for the section (default: "Proceed to Next Step")
#' @return A div element with standardized proceed section
#' @export
create_proceed_section <- function(ns, 
                                  status_output_id, 
                                  button_id, 
                                  button_text,
                                  button_icon = "arrow-right",
                                  button_class = "btn-success",
                                  section_title = "Proceed to Next Step") {
  
  div(
    class = "proceed-section",
    style = "margin-top: 30px; padding: 20px; background-color: #f8f9fa; border-radius: 8px; border-left: 4px solid #28a745;",
    
    # Section header
    h4(section_title, 
       class = "section-title",
       style = "margin-bottom: 15px; color: #2c3e50;"),
    
    # Main content area
    div(
      class = "proceed-content",
      style = "display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 15px;",
      
      # Status message area
      div(
        class = "proceed-status",
        style = "flex-grow: 1; min-width: 300px;",
        uiOutput(ns(status_output_id))
      ),
      
      # Button area
      div(
        class = "proceed-button",
        style = "flex-shrink: 0;",
        actionButton(
          ns(button_id),
          button_text,
          icon = icon(button_icon),
          class = paste("btn-lg", button_class),
          style = "min-width: 250px; font-weight: 600;"
        )
      )
    )
  )
}

#' @title Create Status Message with Icon
#' @description Creates a standardized status message with icon and styling
#' @param message The message text
#' @param status The status type: "success", "warning", "error", "info"
#' @param icon_name The icon name (will auto-select if NULL)
#' @return A div with styled status message
#' @export
create_status_message <- function(message, status = "success", icon_name = NULL) {
  
  # Auto-select icons based on status
  if (is.null(icon_name)) {
    icon_name <- switch(status,
                       "success" = "check-circle",
                       "warning" = "exclamation-triangle", 
                       "error" = "times-circle",
                       "info" = "info-circle",
                       "check-circle")
  }
  
  # Color mapping
  color_map <- c(
    "success" = "#28a745",
    "warning" = "#ffc107", 
    "error" = "#dc3545",
    "info" = "#17a2b8"
  )
  
  color <- if (!is.null(color_map[[status]])) color_map[[status]] else color_map[["info"]]
  
  div(
    class = paste("status-message", paste0("status-", status)),
    style = paste0("color: ", color, "; font-size: 16px; font-weight: 500; display: flex; align-items: center; gap: 8px;"),
    icon(icon_name, style = paste0("color: ", color, ";")),
    span(message)
  )
}

#' @title Create Completion Checker
#' @description Helper function to create completion status reactive
#' @param completion_conditions List of reactive expressions that must be TRUE
#' @return Reactive expression returning completion status
#' @export
create_completion_checker <- function(completion_conditions) {
  reactive({
    all(sapply(completion_conditions, function(condition) {
      tryCatch({
        isTRUE(condition())
      }, error = function(e) FALSE)
    }))
  })
}

# CSS for proceed sections (to be included in main CSS file)