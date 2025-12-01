#!/usr/bin/env Rscript

# Install all packages required to run the QC & Pre-processing Shiny app
# Usage:
# - From repo root: Rscript scripts/install_dependencies.R
# - Or within R: source("scripts/install_dependencies.R")

message("Starting dependency installation...")

# Ensure BiocManager is available
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  message("Installing BiocManager from CRAN...")
  install.packages("BiocManager", repos = "https://cran.r-project.org")
}

# CRAN packages required
cran_packages <- c(
  "shiny",
  "shinythemes",
  "shinyWidgets",
  "DT",
  "htmltools",
  "dplyr",
  "tidyr",
  "readr",
  "ggplot2",
  "plotly",
  "scales",
  "ggrepel",
  "viridis",
  "RColorBrewer",
  "corrplot",
  "moments",
  "e1071",
  "matrixStats",
  "nortest" # used via nortest::ad.test
)

# Bioconductor packages required
bioc_packages <- c(
  "BiocGenerics",
  "S4Vectors",
  "IRanges",
  "GenomeInfoDb",
  "GenomicRanges",
  "SummarizedExperiment",
  "DESeq2",
  "edgeR",
  "preprocessCore"
)

installed <- rownames(installed.packages())

# Install missing CRAN packages
missing_cran <- setdiff(cran_packages, installed)
if (length(missing_cran) > 0) {
  message("Installing CRAN packages: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran, repos = "https://cran.r-project.org", dependencies = TRUE)
} else {
  message("All required CRAN packages are already installed.")
}

# Install missing Bioconductor packages
missing_bioc <- setdiff(bioc_packages, installed)
if (length(missing_bioc) > 0) {
  message("Installing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
} else {
  message("All required Bioconductor packages are already installed.")
}

message("Dependency installation completed.")


