# Quality Control Analysis
# This script performs QC analysis on the RNA-seq data for both tissues

# Load required libraries
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)

# Create output directories
dir.create("results/figures/qc/basic_stats", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/qc/distributions", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# Read DESeq2 objects
dds_hp <- readRDS("data/processed/dds_hippocampus.rds")
dds_m <- readRDS("data/processed/dds_meninges.rds")

# Function to perform QC analysis
perform_qc <- function(dds, tissue) {
    # 1. Basic count statistics
    count_stats <- data.frame(
        total_counts = colSums(counts(dds)),
        detected_genes = colSums(counts(dds) > 0)
    )
    
    # Save count statistics
    write.csv(count_stats, 
              file = sprintf("results/tables/%s_count_statistics.csv", tissue))
    
    # Plot count distributions
    pdf(sprintf("results/figures/qc/basic_stats/%s_count_distributions.pdf", tissue))
    par(mfrow=c(1,2))
    
    # Total counts per sample
    barplot(count_stats$total_counts/1e6, 
            main = sprintf("%s: Total Counts per Sample", tissue),
            ylab = "Millions of Reads",
            las = 2)
    
    # Detected genes per sample
    barplot(count_stats$detected_genes,
            main = sprintf("%s: Detected Genes per Sample", tissue),
            ylab = "Number of Genes",
            las = 2)
    dev.off()
    
    # 2. Variance stabilizing transformation
    vsd <- vst(dds, blind=FALSE)
    
    # 3. PCA analysis
    pcaData <- plotPCA(vsd, intgroup = c("Condition", "Age", "Sex"), returnData = TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    
    # Create PCA plots with different color schemes
    pca_plots <- list()
    
    # PCA by Condition
    pca_plots[["Condition"]] <- ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Sex)) +
        geom_point(size=3) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        ggtitle(sprintf("%s: PCA by Condition", tissue)) +
        theme_minimal()
    
    # PCA by Age
    pca_plots[["Age"]] <- ggplot(pcaData, aes(PC1, PC2, color=Age, shape=Sex)) +
        geom_point(size=3) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        ggtitle(sprintf("%s: PCA by Age", tissue)) +
        theme_minimal()
    
    # Save PCA plots
    pdf(sprintf("results/figures/qc/basic_stats/%s_pca_plots.pdf", tissue))
    for(plot in pca_plots) {
        print(plot)
    }
    dev.off()
    
    # 4. Sample-to-sample distances
    sampleDists <- dist(t(assay(vsd)))
    sampleDistMatrix <- as.matrix(sampleDists)
    
    # Create annotation columns for the heatmap
    annotation_col <- data.frame(
        Condition = dds$Condition,
        Age = dds$Age,
        Sex = dds$Sex
    )
    rownames(annotation_col) <- colnames(dds)
    
    # Create color schemes
    colors <- list(
        Condition = c(Control="#1f77b4", AD="#ff7f0e"),
        Age = c("5"="#2ca02c", "9"="#d62728"),
        Sex = c(M="#9467bd", F="#8c564b")
    )
    
    # Sample distance heatmap
    pdf(sprintf("results/figures/qc/basic_stats/%s_sample_distances.pdf", tissue))
    pheatmap(sampleDistMatrix,
             clustering_distance_rows = sampleDists,
             clustering_distance_cols = sampleDists,
             annotation_col = annotation_col,
             annotation_colors = colors,
             main = sprintf("%s: Sample-to-Sample Distances", tissue))
    dev.off()
    
    # 5. Gene expression distributions
    # First run DESeq
    dds <- estimateSizeFactors(dds)
    dds <- estimateDispersions(dds)
    
    pdf(sprintf("results/figures/qc/distributions/%s_expression_distributions.pdf", tissue))
    plotDispEsts(dds, main=sprintf("%s: Dispersion Estimates", tissue))
    dev.off()
    
    return(list(vsd=vsd, pcaData=pcaData, dds=dds))
}

# Perform QC for both tissues
hp_qc <- perform_qc(dds_hp, "Hippocampus")
m_qc <- perform_qc(dds_m, "Meninges")

# Save QC objects for future use
saveRDS(hp_qc, "data/processed/hippocampus_qc.rds")
saveRDS(m_qc, "data/processed/meninges_qc.rds")

cat("\nQC analysis completed. Results saved in results/figures/qc/ and results/tables/\n") 