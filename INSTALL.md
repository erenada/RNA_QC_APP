# Installation Guide

This guide provides step-by-step instructions for installing and setting up the RNA Processing APP.

## System Requirements

### Minimum Requirements
- **R Version:** 4.0.0 or higher
- **RAM:** 4GB minimum
- **Storage:** 2GB free space
- **Operating System:** Windows 10+, macOS 10.14+, or Linux (Ubuntu 18.04+)
- **Browser:** Chrome, Firefox, Safari, or Edge (latest versions)

### Recommended Requirements
- **R Version:** 4.3.0 or higher
- **RAM:** 8GB or more
- **Storage:** 5GB free space
- **Internet Connection:** Required for initial package installation

## Installation Methods

### Method 1: Quick Setup (Recommended)

1. **Clone the repository:**
   ```bash
   git clone https://github.com/hms-immunology/RNA_QC_APP.git
   cd RNA_QC_APP
   ```

2. **Launch the application:**
   ```r
   # Open R/RStudio and navigate to the project directory
   setwd("path/to/RNA_QC_APP")
   
   # Launch the app (packages install automatically)
   shiny::runApp("app.R")
   ```

The application will automatically detect and install missing packages on first launch.

### Method 2: Manual Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/hms-immunology/RNA_QC_APP.git
   cd RNA_QC_APP
   ```

2. **Install R dependencies:**
   ```r
   # Install package managers
   if (!require("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   
   # Set Bioconductor version
   BiocManager::install(version = "3.17")
   
   # Install CRAN packages
   cran_packages <- c(
     "shiny", "shinythemes", "shinyWidgets", "DT", "htmltools",
     "dplyr", "tidyr", "readr", "ggplot2", "plotly", "scales",
     "pheatmap", "ggrepel", "viridis", "RColorBrewer", "gridExtra",
     "cowplot", "moments", "e1071", "matrixStats", "preprocessCore"
   )
   install.packages(cran_packages)
   
   # Install Bioconductor packages
   bioc_packages <- c(
     "BiocGenerics", "S4Vectors", "IRanges", "GenomeInfoDb",
     "GenomicRanges", "SummarizedExperiment", "DESeq2", "edgeR", "limma"
   )
   BiocManager::install(bioc_packages)
   ```

3. **Launch the application:**
   ```r
   shiny::runApp("app.R")
   ```

## Platform-Specific Instructions

### Windows

1. **Install R:**
   - Download R from [CRAN](https://cran.r-project.org/bin/windows/base/)
   - Install RStudio (optional but recommended) from [RStudio](https://www.rstudio.com/products/rstudio/download/)

2. **Install Rtools (if needed):**
   - Download from [CRAN Rtools](https://cran.r-project.org/bin/windows/Rtools/)
   - Required for compiling some packages

3. **Follow installation methods above**

### macOS

1. **Install R:**
   - Download R from [CRAN](https://cran.r-project.org/bin/macosx/)
   - Install RStudio (optional) from [RStudio](https://www.rstudio.com/products/rstudio/download/)

2. **Install Xcode Command Line Tools:**
   ```bash
   xcode-select --install
   ```

3. **Install dependencies using Homebrew (optional):**
   ```bash
   brew install openssl libxml2 libgit2
   ```

4. **Follow installation methods above**

### Linux (Ubuntu/Debian)

1. **Install R:**
   ```bash
   sudo apt update
   sudo apt install r-base r-base-dev
   ```

2. **Install system dependencies:**
   ```bash
   sudo apt install libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev
   ```

3. **Install RStudio (optional):**
   ```bash
   wget https://download1.rstudio.org/desktop/bionic/amd64/rstudio-2023.06.1-524-amd64.deb
   sudo dpkg -i rstudio-2023.06.1-524-amd64.deb
   ```

4. **Follow installation methods above**

## Troubleshooting Installation

### Common Issues

#### Package Installation Fails
```r
# Clear package cache and retry
remove.packages("problematic_package")
install.packages("problematic_package", dependencies = TRUE)
```

#### BiocManager Issues
```r
# Update BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
```

#### Memory Issues
```r
# Increase memory limit (Windows)
memory.limit(size = 8000)

# Check memory usage
gc()
```

#### Compilation Errors (Windows)
- Ensure Rtools is properly installed
- Check PATH environment variable includes Rtools

#### SSL/TLS Issues
```r
# For HTTPS/SSL issues
options(download.file.method = "libcurl")
```

### Platform-Specific Troubleshooting

#### Windows
- Install Rtools if compilation errors occur
- Use RStudio for better package management
- Ensure antivirus is not blocking R

#### macOS
- Install Xcode command line tools
- Use Homebrew for system dependencies
- Check for ARM vs Intel compatibility (M1/M2 Macs)

#### Linux
- Install development packages (`-dev` packages)
- Ensure proper permissions for R library directory
- Use system package manager when possible

## Verification

After installation, verify the setup:

```r
# Check R version
R.version.string

# Check if key packages are available
required_packages <- c("shiny", "DESeq2", "edgeR", "ggplot2", "DT")
available_packages <- sapply(required_packages, requireNamespace, quietly = TRUE)
print(available_packages)

# Test app loading
source("app.R", echo = FALSE)
print("Installation successful!")
```

## Development Setup

For development work:

```r
# Install additional development packages
install.packages(c("testthat", "devtools", "roxygen2", "styler"))

# Install shinytest2 for UI testing
install.packages("shinytest2")
```

## Getting Help

If you encounter issues:

1. Check the [Troubleshooting Guide](docs/user_manual/troubleshooting.md)
2. Verify [Technical Requirements](docs/user_manual/technical_requirements.md)
3. Try with a fresh R session
4. Create an issue on GitHub with system information

## Next Steps

After successful installation:

1. Read the [Quick Start Guide](docs/user_manual/quick_start.md)
2. Try the example data in `example_data/`
3. Review the [User Manual](docs/user_manual/README.md)

---

**Author:** Eren Ada, PhD  
**Last Updated:** 6/3/2025