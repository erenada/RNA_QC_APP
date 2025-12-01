# options.R
# Application-wide options

# Set max upload size to 1000MB (1GB)
options(shiny.maxRequestSize = 1000 * 1024^2)

# Debug flag (set APP_DEBUG=true in environment to enable verbose logging)
options(app.debug = isTRUE(as.logical(Sys.getenv("APP_DEBUG", "FALSE"))))



