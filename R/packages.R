# packages.R
# Centralized package loading for the Shiny app

suppressPackageStartupMessages({
  # Core Bioconductor infrastructure
  library(methods)
  library(BiocGenerics)
  library(S4Vectors)
  library(IRanges)
  library(GenomeInfoDb)
  library(GenomicRanges)
  library(SummarizedExperiment)
  
  library(DESeq2)
})

# Load remaining required packages
tryCatch({
  suppressPackageStartupMessages({
    # UI packages
    library(shiny)
    library(shinythemes)
    library(shinyWidgets)
    library(DT)
    library(htmltools)
    
    # Data manipulation
    library(dplyr)
    library(tidyr)
    library(readr)
    
    # Visualization
    library(ggplot2)
    library(plotly)
    library(scales)
    # pheatmap not used; corrplot is used for heatmaps
    library(ggrepel)
    library(viridis)
    library(RColorBrewer)
    library(grid)
    library(corrplot)
    
    # Statistical
    library(moments)
    library(e1071)
    
    # Additional
    library(preprocessCore)
    library(methods)  # ensure method dispatch
    
    # Remaining Bioconductor packages
    library(edgeR)
    library(matrixStats)
  })
  message("Package setup completed successfully")
}, error = function(e) {
  message("Error in setup: ", e$message)
  stop(e)
})



