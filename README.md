# Modular Bulk RNA-seq Analysis RShiny Tool Suite

[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-sa/4.0/)
[![R Version](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)](https://www.r-project.org/)

**Author:** Eren Ada, PhD  
**Date:** 12/01/2025  
**Current Version:** 2.0.0

## Project Overview

This project develops a comprehensive and modular bulk RNA-sequencing (RNA-seq) analysis tool using RShiny. The suite consists of interconnected RShiny applications, each dedicated to a specific stage of the RNA-seq analysis workflow.

### Current Status: QC & Pre-processing Tool (Complete)

Production-ready QC and pre-processing for bulk RNA‑seq count data. See Screenshots below for a brief tour.

## Module 1: QC & Pre-processing Tool (Complete)

**Status:** Complete and Production-Ready  
**Location:** `app.R` (main application)

### Features
- Data Input & Validation
- Quality Control: library size, expression distribution, correlation heatmaps, 2D/3D PCA
- Filtering & Normalization: multiple strategies and evaluation views
- Export: processed data and publication-ready plots

### Quick Start
```r
# Install dependencies (once)
Rscript scripts/install_dependencies.R

# Launch from repository root
shiny::runApp("app.R")
```

For detailed usage instructions, see the [Quick Start Guide](docs/user_manual/quick_start.md).

<!-- concise overview retained above -->

## Documentation

### For Users
- **[Complete User Manual](docs/user_manual/README.md)** - Comprehensive guide to all features
- **[Quick Start Guide](docs/user_manual/quick_start.md)** - 15-minute walkthrough
- **[Technical Requirements](docs/user_manual/technical_requirements.md)** - System requirements and setup
- **[Troubleshooting Guide](docs/user_manual/troubleshooting.md)** - Common issues and solutions

### For Developers
- **[Developer Documentation](docs/developer/README.md)** - Architecture and development guidelines
- **[Testing Framework](docs/developer/README.md#testing-framework)** - Testing guidelines and examples
- **[Project Status](docs/project_status.md)** - Current development status

## System Requirements

- **R Version**: 4.0.0 or higher (4.3.0+ recommended)
- **Operating System**: Windows 10+, macOS 10.14+, or Linux (Ubuntu 18.04+)
- **RAM**: 4GB minimum (8GB recommended)
- **Browser**: Chrome, Firefox, Safari, or Edge (latest versions)

## Installation & Setup

### Install Dependencies (Recommended)
Install all required CRAN and Bioconductor packages:

```r
Rscript scripts/install_dependencies.R

# Then launch the application
shiny::runApp("app.R")
```

### Manual Setup
For development or troubleshooting:

```r
# Install package manager
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install core dependencies
BiocManager::install(c("DESeq2", "edgeR"))
install.packages(c("shiny", "DT", "plotly", "ggplot2"))

# Note: The app uses additional packages (ggrepel, corrplot, shinythemes, shinyWidgets, 
# dplyr, tidyr, readr, viridis, RColorBrewer, scales, moments, e1071, matrixStats, 
# and Bioconductor preprocessCore, SummarizedExperiment and friends). 
# To avoid omissions, prefer running: Rscript scripts/install_dependencies.R
```

## Screenshots & Interface Overview
Key views:
- Data Input & Validation
- QC: PCA (2D/3D), correlation heatmaps
- Filtering & Normalization with evaluation

Screenshots:
- Validation: `screenshots/validation.png`
- PCA: `screenshots/pca.png`
- Correlation: `screenshots/sample_correlation.png`
- Filtering/Normalization: `screenshots/filtering_normalization.png`
- About: `screenshots/about.png`

## Example Usage

### Basic Workflow
1. **Data Upload**: Load count matrix and metadata (CSV format)
2. **Validation**: Review data validation results and handle any issues
3. **Quality Control**: Examine QC plots and assess data quality
4. **Processing**: Apply filtering and normalization based on QC results
5. **Export**: Download processed data and reports

### Example Data
Sample datasets are provided in the `example_data/` directory for testing and learning:

- **`example_counts.csv`** - Anonymized count matrix (5,001 genes × 24 samples)
- **`example_metadata.csv`** - Corresponding sample metadata with experimental design
 

 

 

## Contributing

Contributions are welcome. See the [Developer Documentation](docs/developer/README.md) for setup, standards, testing, and PR guidance.

### Quick Contribution Steps
1. Fork the repository
2. Create a feature branch
3. Follow coding standards in [Developer Guide](docs/developer/README.md)
4. Add tests for new functionality
5. Update documentation
6. Submit pull request

## Support & Community

### Getting Help
1. Check the [Troubleshooting Guide](docs/user_manual/troubleshooting.md)
2. Review [Technical Requirements](docs/user_manual/technical_requirements.md)
3. Try with example data to isolate issues
4. Create an issue with detailed information

 

## License

This project is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this tool in your research, please cite:

```bibtex
@software{ada2025rna,
  author = {Ada, Eren},
  title = {Modular Bulk RNA-seq Analysis RShiny Tool Suite},
  year = {2025},
  version = {2.0.0},
  url = {https://github.com/hms-immunology/RNA_QC_APP}
}
```

## Acknowledgments

This project builds upon the excellent work of the Bioconductor community, particularly the DESeq2, edgeR, and limma packages for RNA-seq analysis.

 