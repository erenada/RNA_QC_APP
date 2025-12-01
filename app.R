# app.R
# QC & Pre-processing Tool
# Author: Eren Ada, PhD
# Date: 05/29/2024

# Centralized configuration and packages
source("R/options.R")
source("R/packages.R")

source("R/modules_source.R")

ui <- app_ui()

server <- app_server

# Run the application
shinyApp(ui = ui, server = server)
