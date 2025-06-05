# Module 1: QC & Pre-processing

## Overview
This module provides quality control and pre-processing functionality for RNA-seq data analysis. It handles data input validation, duplicate gene management, and basic QC metrics calculation.

## Directory Structure
```
module1_qc_preprocessing/
├── DESCRIPTION          # Package dependencies and metadata
├── LICENSE             # MIT License
├── README.md           # This file
├── app.R              # Main application file
├── R/                  # R source files
│   ├── ui_tab1.R      # UI components for data input tab
│   └── server_tab1.R  # Server logic for data input tab
├── tests/              # Test files
│   └── testthat/      # Unit tests
├── www/               # Web assets
│   └── styles.css     # Global CSS styles
```

## Installation & Setup

1. First, install R (>= 4.1.0) and RStudio (recommended).

2. Install required packages:
```R
# Install required CRAN packages
install.packages(c(
  "shiny",
  "shinythemes",
  "shinyWidgets",
  "DT",
  "dplyr",
  "tidyr",
  "readr",
  "matrixStats",
  "ggplot2",
  "plotly",
  "htmltools"
))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("SummarizedExperiment", "edgeR"))
```

## Running the App

### Method 1: Using RStudio
1. Open RStudio
2. Navigate to the `module1_qc_preprocessing` directory
3. Open `app.R`
4. Click the "Run App" button in RStudio, or run:
```R
shiny::runApp()
```

### Method 2: From R Console
1. Set your working directory to the module folder:
```R
setwd("/path/to/APP/modules/module1_qc_preprocessing")
shiny::runApp()
```

### Method 3: From Terminal
```bash
cd /path/to/APP/modules/module1_qc_preprocessing
R -e "shiny::runApp()"
```

## Input Data Format

### Count Matrix File
- CSV, TSV, or TXT format
- Rows: genes
- Columns: samples
- Optional: gene IDs in a separate column

Example:
```
GeneID,Sample1,Sample2,Sample3
GENE1,100,150,200
GENE2,50,75,100
```

### Metadata File
- CSV, TSV, or TXT format
- Rows: samples
- Columns: sample attributes
- Must include a column matching sample names in count matrix

Example:
```
SampleID,Condition,Treatment
Sample1,Control,A
Sample2,Treatment,B
Sample3,Treatment,A
```

## Features
- File upload for count matrix and metadata
- Flexible parsing options
- Duplicate gene handling strategies
- Data validation and QC metrics
- Interactive data previews
- Comprehensive validation summary

## Development
- Author: Eren Ada, PhD
- Version: 0.1.0
- License: MIT

## Testing
Run tests using:
```R
testthat::test_dir("tests/testthat")
```

## Troubleshooting

### Common Issues
1. **Package Loading Errors**
   - Ensure all packages are installed correctly
   - Check package versions match requirements in DESCRIPTION file

2. **File Upload Issues**
   - Verify file format (CSV, TSV, TXT)
   - Check file encoding (UTF-8 recommended)
   - Ensure file size is reasonable

3. **Data Format Problems**
   - Confirm count matrix contains numeric values
   - Verify sample names match between count matrix and metadata
   - Check for missing values or special characters

### Getting Help
- Check the console for error messages
- Review the validation summary in the app
- File an issue on the project repository 