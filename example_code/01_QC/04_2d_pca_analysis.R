# 2D PCA Analysis with Multiple Factor Visualizations
# This script creates detailed 2D PCA plots with different factor combinations

# Load required libraries
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

# Create output directories
dir.create("results/figures/qc/pca", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# Read QC objects
hp_qc <- readRDS("data/processed/hippocampus_qc.rds")
m_qc <- readRDS("data/processed/meninges_qc.rds")

# Function to create 2D PCA plots
create_2d_pca_plots <- function(dds, tissue, min_count = 10, min_samples = 5) {
    # Original PCA without filtering
    vsd_no_filter <- vst(dds, blind=FALSE)
    pca_no_filter <- prcomp(t(assay(vsd_no_filter)))
    
    # Filter low count genes
    keep <- rowSums(counts(dds) >= min_count) >= min_samples
    dds_filtered <- dds[keep,]
    
    # PCA with filtering
    vsd_filtered <- vst(dds_filtered, blind=FALSE)
    pca_filtered <- prcomp(t(assay(vsd_filtered)))
    
    # Create data frames for plotting
    create_plot_df <- function(pca_obj, metadata) {
        var_explained <- round(100 * pca_obj$sdev^2 / sum(pca_obj$sdev^2), 1)
        data.frame(
            PC1 = pca_obj$x[,1],
            PC2 = pca_obj$x[,2],
            PC3 = pca_obj$x[,3],
            Condition = metadata$Condition,
            Age = metadata$Age,
            Sex = metadata$Sex,
            Group = paste0(metadata$Condition, "_", metadata$Age, "m_", metadata$Sex)
        )
    }
    
    df_no_filter <- create_plot_df(pca_no_filter, colData(dds))
    df_filtered <- create_plot_df(pca_filtered, colData(dds))
    
    # Calculate variance explained
    var_explained_no_filter <- round(100 * pca_no_filter$sdev^2 / sum(pca_no_filter$sdev^2), 1)
    var_explained_filtered <- round(100 * pca_filtered$sdev^2 / sum(pca_filtered$sdev^2), 1)
    
    # Create different plot versions
    plot_versions <- list()
    
    # Function to create base plot
    create_base_plot <- function(data, var_explained, title_suffix) {
        list(
            # 1. Condition only
            ggplot(data, aes(x=PC1, y=PC2, color=Condition)) +
                geom_point(size=3) +
                theme_minimal() +
                labs(
                    x = paste0("PC1 (", var_explained[1], "%)"),
                    y = paste0("PC2 (", var_explained[2], "%)"),
                    title = paste0(tissue, ": PCA by Condition ", title_suffix)
                ),
            
            # 2. Condition + Sex shapes
            ggplot(data, aes(x=PC1, y=PC2, color=Condition, shape=Sex)) +
                geom_point(size=3) +
                theme_minimal() +
                labs(
                    x = paste0("PC1 (", var_explained[1], "%)"),
                    y = paste0("PC2 (", var_explained[2], "%)"),
                    title = paste0(tissue, ": PCA by Condition and Sex ", title_suffix)
                ),
            
            # 3. Age + Sex
            ggplot(data, aes(x=PC1, y=PC2, color=Age, shape=Sex)) +
                geom_point(size=3) +
                theme_minimal() +
                labs(
                    x = paste0("PC1 (", var_explained[1], "%)"),
                    y = paste0("PC2 (", var_explained[2], "%)"),
                    title = paste0(tissue, ": PCA by Age and Sex ", title_suffix)
                ),
            
            # 4. All groups
            ggplot(data, aes(x=PC1, y=PC2, color=Group)) +
                geom_point(size=3) +
                theme_minimal() +
                labs(
                    x = paste0("PC1 (", var_explained[1], "%)"),
                    y = paste0("PC2 (", var_explained[2], "%)"),
                    title = paste0(tissue, ": PCA by Groups ", title_suffix)
                ) +
                theme(legend.position = "right")
        )
    }
    
    # Create plots for both filtered and unfiltered data
    plots_no_filter <- create_base_plot(df_no_filter, var_explained_no_filter, "(Unfiltered)")
    plots_filtered <- create_base_plot(df_filtered, var_explained_filtered, "(Filtered)")
    
    # Save plots
    pdf(sprintf("results/figures/qc/pca/%s_2D_PCA_plots.pdf", tissue), width=12, height=10)
    # Unfiltered plots
    grid.arrange(grobs=plots_no_filter, ncol=2, 
                top=grid::textGrob(paste0(tissue, " - Unfiltered Data"), gp=grid::gpar(fontsize=14)))
    # Filtered plots
    grid.arrange(grobs=plots_filtered, ncol=2,
                top=grid::textGrob(paste0(tissue, " - Filtered Data"), gp=grid::gpar(fontsize=14)))
    dev.off()
    
    # Return variance explained
    return(list(
        var_explained_no_filter = var_explained_no_filter,
        var_explained_filtered = var_explained_filtered
    ))
}

# Create plots for both tissues
hp_pca_stats <- create_2d_pca_plots(hp_qc$dds, "Hippocampus", min_count = 10, min_samples = 5)
m_pca_stats <- create_2d_pca_plots(m_qc$dds, "Meninges", min_count = 10, min_samples = 5)

# Save PCA statistics
pca_stats <- data.frame(
    Tissue = rep(c("Hippocampus", "Meninges"), each=2),
    Data_Type = rep(c("Unfiltered", "Filtered"), 2),
    PC1_var = c(hp_pca_stats$var_explained_no_filter[1], 
                hp_pca_stats$var_explained_filtered[1],
                m_pca_stats$var_explained_no_filter[1], 
                m_pca_stats$var_explained_filtered[1]),
    PC2_var = c(hp_pca_stats$var_explained_no_filter[2], 
                hp_pca_stats$var_explained_filtered[2],
                m_pca_stats$var_explained_no_filter[2], 
                m_pca_stats$var_explained_filtered[2])
)

write.csv(pca_stats, "results/tables/pca_variance_explained.csv", row.names=FALSE)

cat("\n2D PCA analysis completed.\n")
cat("Plots saved in results/figures/qc/pca/\n")
cat("PCA statistics saved in results/tables/pca_variance_explained.csv\n") 