# Developer Documentation

**Author:** Eren Ada, PhD  
**Date:** 5/30/2025  
**Version:** 2.0.0

## Project Overview

The RNA-seq QC & Pre-processing Tool is a modular Shiny application designed for comprehensive quality control and pre-processing of bulk RNA-sequencing data. This documentation covers the technical architecture, development practices, and guidelines for contributors.

## Table of Contents

1. [Architecture Overview](#architecture-overview)
2. [Module Structure](#module-structure)
3. [Code Organization](#code-organization)
4. [Development Guidelines](#development-guidelines)
5. [Testing Framework](#testing-framework)
6. [Deployment](#deployment)
7. [Contributing](#contributing)

## Architecture Overview

### Design Principles
- **Modularity**: Each major functionality is encapsulated in separate modules
- **Reactivity**: Efficient data flow using Shiny's reactive programming
- **Scalability**: Designed to handle large datasets with memory optimization
- **Maintainability**: Clear separation of UI, server logic, and utility functions

### Technology Stack
- **Framework**: R Shiny 1.7+
- **Data Processing**: Bioconductor packages (DESeq2, edgeR, limma)
- **Visualization**: ggplot2, plotly, pheatmap
- **UI Components**: shinythemes, shinyWidgets, DT

### Application Flow
```
app.R (main)
├── UI Definition
│   ├── Data Input & Validation Tab
│   ├── QC Plots & Summaries Tab
│   ├── Filtering & Normalization Tab
│   └── About Tab
└── Server Logic
    ├── Shared Reactive Values
    ├── Module Servers
    └── Tab Navigation Logic
```

## Module Structure

### Core Modules

#### 1. Data Input & Validation (`mod_input_validation`)
**Location**: `modules/module1_qc_preprocessing/R/`
- `ui_tab1.R`: File upload and configuration UI
- `server_tab1.R`: Data processing and validation logic

**Key Functions**:
- File upload handling
- Data format validation
- Duplicate gene detection and handling
- Sample consistency checking

#### 2. QC Plots & Summaries (`mod_qc_plots`)
**Location**: `modules/module1_qc_preprocessing/R/`
- `ui_tab2_qc_plots.R`: QC visualization UI
- `server_tab2_qc_plots.R`: Plot generation logic
- `utils_qc_plotting.R`: Plotting utility functions

**Key Functions**:
- Library size analysis
- Expression distribution plots
- Sample correlation heatmaps
- PCA analysis (2D and 3D)

#### 3. Filtering & Normalization (`mod_filtering_normalization`)
**Location**: `modules/module1_qc_preprocessing/R/`
- `ui_tab3_filtering_normalization.R`: Processing controls UI
- `server_tab3_filtering_normalization.R`: Processing logic
- `utils_filtering.R`: Gene filtering functions
- `utils_normalization.R`: Normalization methods
- `utils_plotting.R`: Before/after comparison plots

**Key Functions**:
- Multiple filtering strategies
- Various normalization methods
- Impact visualization
- Data export functionality

### Module Communication

Modules communicate through shared reactive values:

```r
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
```

## Code Organization

### Directory Structure
```
APP/
├── app.R                           # Main application file
├── modules/
│   └── module1_qc_preprocessing/
│       ├── R/                      # Module R files
│       ├── tests/                  # Unit tests
│       ├── templates/              # Code templates
│       └── www/                    # Static assets
├── docs/                          # Documentation
├── example_data/                  # Sample datasets
├── logs/                          # Application logs
└── inst/                          # Package assets
```

### Naming Conventions

#### Files
- **UI files**: `ui_[module]_[description].R`
- **Server files**: `server_[module]_[description].R`
- **Utility files**: `utils_[functionality].R`
- **Test files**: `test_[module]_[function].R`

#### Functions
- **Module UI**: `mod_[module_name]_ui()`
- **Module Server**: `mod_[module_name]_server()`
- **Utility functions**: `[verb]_[noun]()`
- **Validation functions**: `validate_[data_type]()`

#### Variables
- **Reactive values**: `snake_case`
- **Input IDs**: `snake_case`
- **Constants**: `UPPER_CASE`

### Coding Standards

#### R Code Style
```r
# Function definition
process_counts <- function(counts_matrix, 
                          metadata, 
                          min_count = 10) {
  # Validate inputs
  stopifnot(is.matrix(counts_matrix))
  stopifnot(is.data.frame(metadata))
  
  # Process data
  filtered_counts <- filter_low_counts(
    counts = counts_matrix,
    threshold = min_count
  )
  
  return(filtered_counts)
}
```

#### Shiny Module Pattern
```r
# UI Module
mod_example_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Module Title"),
    fileInput(ns("file"), "Upload File"),
    verbatimTextOutput(ns("status"))
  )
}

# Server Module
mod_example_server <- function(id, shared_data) {
  moduleServer(id, function(input, output, session) {
    # Reactive logic here
    output$status <- renderText({
      if (shared_data$data_processed) {
        "Data ready"
      } else {
        "No data loaded"
      }
    })
    
    # Return reactive values
    return(reactive({
      list(
        status = "complete"
      )
    }))
  })
}
```

## Development Guidelines

### Getting Started

1. **Clone Repository**
   ```bash
   git clone [repository-url]
   cd APP
   ```

2. **Install Dependencies**
   ```r
   # The app will automatically install missing packages
   # For development, install additional packages:
   install.packages(c("testthat", "devtools", "roxygen2"))
   ```

3. **Run Application**
   ```r
   shiny::runApp("app.R")
   ```

### Development Workflow

1. **Feature Development**
   - Create feature branch from main
   - Implement changes in appropriate module
   - Add tests for new functionality
   - Update documentation

2. **Testing**
   ```r
   # Run tests
   testthat::test_dir("modules/module1_qc_preprocessing/tests/")
   
   # Test with example data
   # Use files in example_data/ directory
   ```

3. **Documentation**
   - Update function documentation
   - Add examples to user manual
   - Update version numbers

### Error Handling

#### User-Facing Errors
```r
# Use showNotification for user feedback
if (!file.exists(input$file$datapath)) {
  showNotification(
    "File upload failed. Please try again.",
    type = "error",
    duration = 5
  )
  return(NULL)
}
```

#### Validation Functions
```r
validate_count_matrix <- function(counts) {
  errors <- c()
  
  if (!is.matrix(counts) && !is.data.frame(counts)) {
    errors <- c(errors, "Count data must be a matrix or data frame")
  }
  
  if (any(counts < 0, na.rm = TRUE)) {
    errors <- c(errors, "Count values must be non-negative")
  }
  
  if (length(errors) > 0) {
    stop(paste(errors, collapse = "; "))
  }
  
  return(TRUE)
}
```

### Performance Optimization

#### Memory Management
```r
# Use gc() after large operations
process_large_dataset <- function(data) {
  result <- heavy_computation(data)
  gc()  # Force garbage collection
  return(result)
}

# Avoid copying large objects unnecessarily
# Use in-place operations when possible
```

#### Reactive Optimization
```r
# Use req() to prevent unnecessary computation
filtered_data <- reactive({
  req(input$file)
  req(shared_data$data_processed)
  
  apply_filters(shared_data$processed_counts)
})

# Use isolate() to prevent reactive dependencies
observe({
  req(input$update_button)
  isolate({
    # Non-reactive code here
  })
})
```

## Testing Framework

### Test Structure
```
modules/module1_qc_preprocessing/tests/
├── testthat/
│   ├── test_data_validation.R
│   ├── test_filtering.R
│   ├── test_normalization.R
│   └── test_plotting.R
└── test_data/
    ├── example_counts.csv
    └── example_metadata.csv
```

### Test Examples
```r
# Test file: test_data_validation.R
test_that("validate_count_matrix handles valid input", {
  # Create test data
  counts <- matrix(
    c(10, 20, 30, 5, 15, 25),
    nrow = 2,
    dimnames = list(c("Gene1", "Gene2"), c("Sample1", "Sample2", "Sample3"))
  )
  
  # Test validation
  expect_true(validate_count_matrix(counts))
})

test_that("validate_count_matrix rejects negative values", {
  counts <- matrix(c(-1, 5, 10, 20), nrow = 2)
  expect_error(validate_count_matrix(counts))
})
```

### Running Tests
```r
# Run all tests
testthat::test_dir("modules/module1_qc_preprocessing/tests/")

# Run specific test file
testthat::test_file("modules/module1_qc_preprocessing/tests/testthat/test_data_validation.R")

# Run with coverage
covr::package_coverage(".")
```

## Deployment

### Local Deployment
```r
# Standard local deployment
shiny::runApp("app.R")

# With specific host/port
shiny::runApp("app.R", host = "0.0.0.0", port = 3838)
```

### Server Deployment

#### Shiny Server
1. Copy application to `/srv/shiny-server/`
2. Ensure all dependencies are installed
3. Configure `shiny-server.conf`

#### Docker Deployment
```dockerfile
FROM rocker/shiny:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev

# Install R packages
COPY install_packages.R /tmp/
RUN Rscript /tmp/install_packages.R

# Copy application
COPY . /srv/shiny-server/rna-qc/
COPY shiny-server.conf /etc/shiny-server/

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
```

## Contributing

### Pull Request Process

1. **Fork & Branch**
   ```bash
   git checkout -b feature/new-functionality
   ```

2. **Development**
   - Follow coding standards
   - Add appropriate tests
   - Update documentation

3. **Testing**
   ```r
   # Run full test suite
   testthat::test_dir("modules/module1_qc_preprocessing/tests/")
   
   # Manual testing with example data
   ```

4. **Documentation**
   - Update function documentation
   - Add user manual sections if needed
   - Update version information

5. **Submission**
   - Create pull request with detailed description
   - Include test results
   - Reference any related issues

### Code Review Guidelines

#### For Reviewers
- Check code follows established patterns
- Verify tests cover new functionality
- Ensure documentation is updated
- Test with example data

#### For Contributors
- Provide clear commit messages
- Include rationale for design decisions
- Respond promptly to review feedback
- Keep pull requests focused and manageable

### Issue Reporting

When reporting bugs or requesting features:

1. **Use Issue Templates**
2. **Provide System Information**
   - R version
   - Operating system
   - Browser (for UI issues)

3. **Include Reproducible Example**
   - Minimal code to reproduce issue
   - Sample data if relevant
   - Expected vs actual behavior

4. **Categorize Appropriately**
   - Bug report
   - Feature request
   - Documentation improvement
   - Performance issue

---

For questions about development, consult this documentation first, then reach out to the development team. 