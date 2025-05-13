# Enhanced PCA Analysis with 3D visualization
# This script creates detailed PCA plots including 3D visualization

# Load required libraries
library(DESeq2)
library(tidyverse)
library(plotly)
library(htmlwidgets)

# Read QC objects
hp_qc <- readRDS("data/processed/hippocampus_qc.rds")
m_qc <- readRDS("data/processed/meninges_qc.rds")

# Function to create group labels
create_group_labels <- function(metadata) {
    paste0(metadata$Condition, "_", metadata$Age, "m_", metadata$Sex)
}

# Function to perform PCA with filtering options
perform_enhanced_pca <- function(dds, tissue, min_count = 10, min_samples = 5) {
    # Original PCA without filtering
    vsd_no_filter <- vst(dds, blind=FALSE)
    pca_no_filter <- prcomp(t(assay(vsd_no_filter)))
    
    # Filter low count genes
    keep <- rowSums(counts(dds) >= min_count) >= min_samples
    dds_filtered <- dds[keep,]
    
    # PCA with filtering
    vsd_filtered <- vst(dds_filtered, blind=FALSE)
    pca_filtered <- prcomp(t(assay(vsd_filtered)))
    
    # Create group labels
    groups <- create_group_labels(colData(dds))
    
    # Calculate variance explained
    var_explained_no_filter <- round(100 * pca_no_filter$sdev^2 / sum(pca_no_filter$sdev^2), 1)
    var_explained_filtered <- round(100 * pca_filtered$sdev^2 / sum(pca_filtered$sdev^2), 1)
    
    # Create 3D PCA plots
    # Unfiltered data
    p1 <- plot_ly(x = pca_no_filter$x[,1], 
                  y = pca_no_filter$x[,2], 
                  z = pca_no_filter$x[,3],
                  type = "scatter3d",
                  mode = "markers",
                  color = groups,
                  text = colnames(dds),
                  marker = list(size = 8)) %>%
          layout(title = paste0(tissue, ": 3D PCA (Unfiltered)"),
                 scene = list(
                     xaxis = list(title = paste0("PC1 (", var_explained_no_filter[1], "%)")),
                     yaxis = list(title = paste0("PC2 (", var_explained_no_filter[2], "%)")),
                     zaxis = list(title = paste0("PC3 (", var_explained_no_filter[3], "%)"))
                 ))
    
    # Filtered data
    p2 <- plot_ly(x = pca_filtered$x[,1], 
                  y = pca_filtered$x[,2], 
                  z = pca_filtered$x[,3],
                  type = "scatter3d",
                  mode = "markers",
                  color = groups,
                  text = colnames(dds),
                  marker = list(size = 8)) %>%
          layout(title = paste0(tissue, ": 3D PCA (Filtered)"),
                 scene = list(
                     xaxis = list(title = paste0("PC1 (", var_explained_filtered[1], "%)")),
                     yaxis = list(title = paste0("PC2 (", var_explained_filtered[2], "%)")),
                     zaxis = list(title = paste0("PC3 (", var_explained_filtered[3], "%)"))
                 ))
    
    # Save plots
    saveWidget(p1, file = sprintf("results/figures/qc/figures/%s_3D_PCA_unfiltered.html", tissue))
    saveWidget(p2, file = sprintf("results/figures/qc/figures/%s_3D_PCA_filtered.html", tissue))
    
    # Return filtering statistics
    stats <- data.frame(
        total_genes = nrow(dds),
        genes_after_filtering = nrow(dds_filtered),
        percent_kept = round(100 * nrow(dds_filtered) / nrow(dds), 1)
    )
    
    return(list(
        pca_no_filter = pca_no_filter,
        pca_filtered = pca_filtered,
        var_explained_no_filter = var_explained_no_filter,
        var_explained_filtered = var_explained_filtered,
        filtering_stats = stats
    ))
}

# Perform enhanced PCA for both tissues
hp_pca <- perform_enhanced_pca(hp_qc$dds, "Hippocampus", min_count = 10, min_samples = 5)
m_pca <- perform_enhanced_pca(m_qc$dds, "Meninges", min_count = 10, min_samples = 5)

# Save filtering statistics
write.csv(rbind(
    cbind(Tissue = "Hippocampus", hp_pca$filtering_stats),
    cbind(Tissue = "Meninges", m_pca$filtering_stats)
), file = "results/tables/qc/basic_stats/filtering_statistics.csv", row.names = FALSE)

cat("\nEnhanced PCA analysis completed.\n")
cat("3D PCA plots saved as interactive HTML files in results/figures/qc/figures/\n")
cat("Filtering statistics saved in results/tables/qc/basic_stats/filtering_statistics.csv\n") 