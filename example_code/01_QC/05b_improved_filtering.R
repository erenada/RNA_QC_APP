#!/usr/bin/env Rscript

# Improved Filtering Script with Condition-Specific Criteria
# Author: Eren Ada, PhD
# Date: April 2025
# Purpose: Implement condition-specific filtering to better capture down-regulated genes
#         while maintaining data quality and providing comprehensive QC metrics

# Load required libraries
suppressPackageStartupMessages({
    library(tidyverse)
    library(edgeR)
    library(ggplot2)
    library(gridExtra)
    library(pheatmap)
    library(reshape2)
    library(yaml)
    library(VennDiagram)
})

# Create output directories
output_dir <- "results/qc/improved_filtering"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Initialize logging
log_file <- file.path(output_dir, "filtering_log.txt")
write_log <- function(message) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    entry <- paste0("[", timestamp, "] ", message, "\n")
    cat(entry, file = log_file, append = TRUE)
    cat(message, "\n")
}

write_log("=== Starting Improved Filtering Analysis ===")

# Read the raw count matrix
write_log("Reading count matrix...")
counts <- read.csv("data/processed/counts_matrix_summed.csv")
write_log(sprintf("Initial number of genes: %d", nrow(counts)))

# Read metadata
write_log("Reading metadata...")
metadata <- readRDS("data/processed/metadata_processed.rds")

# Define experimental factors
tissues <- c("CE", "DRG", "MNG", "Sciat", "SCord")  # Using tissue codes from metadata
conditions <- c("Mock", "Inf")
genotypes <- c("B6", "IL10KO")

# Function to calculate CPM
calculate_cpm <- function(count_matrix) {
    lib_sizes <- colSums(count_matrix)
    t(t(count_matrix) / lib_sizes * 1e6)
}

# Function to get sample indices
get_sample_indices <- function(metadata, tissue, condition, genotype) {
    metadata %>%
        filter(
            Tissue == tissue,
            Treatment == condition,
            Genotype == genotype
        ) %>%
        pull(Samples)  # Updated to use correct column name
}

# Function to calculate coefficient of variation
calc_cv <- function(x) sd(x) / mean(x)

# Function to perform condition-specific filtering
filter_condition_specific <- function(counts, metadata, tissue) {
    write_log(sprintf("\nProcessing tissue: %s", tissue))
    
    # Initialize results storage
    passing_genes <- matrix(FALSE, nrow = nrow(counts), ncol = 4)
    colnames(passing_genes) <- c("B6_Mock", "B6_Inf", "IL10KO_Mock", "IL10KO_Inf")
    
    # Calculate CPM
    count_matrix <- as.matrix(counts[, -1])  # Remove gene name column
    cpm_matrix <- calculate_cpm(count_matrix)
    
    # Process each condition combination
    for(genotype in genotypes) {
        for(condition in conditions) {
            # Get relevant samples
            samples <- get_sample_indices(metadata, tissue, condition, genotype)
            
            if(length(samples) < 2) {
                write_log(sprintf("Warning: Insufficient samples for %s_%s", genotype, condition))
                next
            }
            
            # Get CPM for these samples
            condition_cpm <- cpm_matrix[, samples]
            
            # Apply less stringent criteria (CPM > 0.5 in 2/3 replicates)
            pass_threshold <- rowSums(condition_cpm > 0.5) >= 2
            
            # Store results
            col_name <- paste0(genotype, "_", condition)
            passing_genes[, col_name] <- pass_threshold
        }
    }
    
    # A gene passes if it's well-expressed in ANY condition
    genes_to_keep <- rowSums(passing_genes) > 0
    
    # Calculate additional metrics for passing genes
    cv_stats <- data.frame(
        GeneName = counts$GeneName,
        row.names = counts$GeneName
    )
    
    for(genotype in genotypes) {
        for(condition in conditions) {
            samples <- get_sample_indices(metadata, tissue, condition, genotype)
            condition_cpm <- cpm_matrix[, samples]
            cv_stats[[paste0(genotype, "_", condition, "_CV")]] <- 
                apply(condition_cpm, 1, calc_cv)
        }
    }
    
    list(
        passing_genes = genes_to_keep,
        cv_stats = cv_stats,
        pass_matrix = passing_genes
    )
}

# Process each tissue
tissue_results <- list()
for(tissue in tissues) {
    tissue_results[[tissue]] <- filter_condition_specific(counts, metadata, tissue)
}

# Combine results across tissues
all_passing_genes <- unique(unlist(lapply(tissue_results, function(x) 
    which(x$passing_genes))))

# Filter the count matrix
counts_filtered_new <- counts[all_passing_genes, ]

# Create comparison plots
write_log("\nGenerating comparison plots...")

# 1. Venn diagram of genes passing in different conditions
pdf(file.path(output_dir, "condition_overlap.pdf"))
for(tissue in tissues) {
    pass_matrix <- tissue_results[[tissue]]$pass_matrix
    
    # Create Venn diagram
    venn.plot <- venn.diagram(
        x = list(
            B6_Mock = which(pass_matrix[, "B6_Mock"]),
            B6_Inf = which(pass_matrix[, "B6_Inf"]),
            IL10KO_Mock = which(pass_matrix[, "IL10KO_Mock"]),
            IL10KO_Inf = which(pass_matrix[, "IL10KO_Inf"])
        ),
        filename = NULL,
        main = paste("Gene Overlap -", tissue),
        fill = c("red", "blue", "green", "yellow"),
        alpha = 0.5,
        cex = 1,
        cat.cex = 0.8,
        margin = 0.1
    )
    
    # Draw the plot
    grid.draw(venn.plot)
}
dev.off()

# 2. Expression distribution plots
plot_expression_distribution <- function(counts_old, counts_new, tissue) {
    old_cpm <- calculate_cpm(as.matrix(counts_old[, -1]))
    new_cpm <- calculate_cpm(as.matrix(counts_new[, -1]))
    
    old_data <- data.frame(CPM = log2(as.vector(old_cpm) + 1),
                          Type = "Original Filtering")
    new_data <- data.frame(CPM = log2(as.vector(new_cpm) + 1),
                          Type = "New Filtering")
    
    combined_data <- rbind(old_data, new_data)
    
    ggplot(combined_data, aes(x = CPM, fill = Type)) +
        geom_density(alpha = 0.5) +
        theme_minimal() +
        labs(title = paste("Expression Distribution -", tissue),
             x = "log2(CPM + 1)",
             y = "Density")
}

# 3. CV distribution plots
plot_cv_distribution <- function(tissue_result, tissue) {
    cv_data <- reshape2::melt(tissue_result$cv_stats[, grep("_CV$", colnames(tissue_result$cv_stats))])
    
    ggplot(cv_data, aes(x = value, fill = variable)) +
        geom_density(alpha = 0.5) +
        theme_minimal() +
        labs(title = paste("CV Distribution -", tissue),
             x = "Coefficient of Variation",
             y = "Density")
}

# Generate comparison metrics
write_log("\nCalculating comparison metrics...")

# Read original filtered data
counts_filtered_old <- read.csv("data/processed/counts_final_filtered.csv")

comparison_stats <- data.frame(
    Metric = c(
        "Total Genes - Original",
        "Total Genes - New",
        "Unique to Original",
        "Unique to New",
        "Common to Both"
    ),
    Value = c(
        nrow(counts_filtered_old),
        nrow(counts_filtered_new),
        sum(!counts_filtered_old$GeneName %in% counts_filtered_new$GeneName),
        sum(!counts_filtered_new$GeneName %in% counts_filtered_old$GeneName),
        sum(counts_filtered_new$GeneName %in% counts_filtered_old$GeneName)
    )
)

# Save results
write_log("\nSaving results...")

# Save new filtered counts
write.csv(counts_filtered_new,
          "data/processed/counts_final_filtered_improved.csv",
          row.names = FALSE)

# Save comparison metrics
write.csv(comparison_stats,
          file.path(output_dir, "filtering_comparison_metrics.csv"),
          row.names = FALSE)

# Generate comprehensive report
sink(file.path(output_dir, "filtering_report.md"))

cat("# Improved Filtering Analysis Report\n")
cat("## Author: Eren Ada, PhD\n")
cat(paste("## Date:", format(Sys.Date(), "%B %d, %Y"), "\n\n"))

cat("## 1. Filtering Strategy\n")
cat("### Modified Criteria:\n")
cat("- CPM > 0.5 in at least 2/3 replicates in ANY condition\n")
cat("- Condition-specific filtering to capture condition-specific expression\n")
cat("- No minimum threshold for low-expression condition\n\n")

cat("## 2. Comparison Metrics\n")
print(comparison_stats)

cat("\n## 3. Quality Control\n")
cat("- Expression distribution plots generated\n")
cat("- CV distribution plots generated\n")
cat("- Condition overlap analysis performed\n\n")

cat("## 4. Files Generated\n")
cat("1. Filtered count matrix: data/processed/counts_final_filtered_improved.csv\n")
cat("2. QC plots: results/qc/improved_filtering/\n")
cat("3. Comparison metrics: results/qc/improved_filtering/filtering_comparison_metrics.csv\n")
cat("4. This report\n")

sink()

write_log("\nAnalysis completed successfully!") 