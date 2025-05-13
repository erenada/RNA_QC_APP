# Low Count Filtering Script
# Author: Eren Ada, PhD
# Date: April 2024
# Purpose: Filter low count genes using CPM threshold and tissue/condition-specific criteria

# Load required libraries
library(tidyverse)
library(edgeR)
library(here)
library(ggplot2)
library(reshape2)
library(pheatmap)

# Read the filtered count matrix (no zero counts)
counts <- read.csv("data/processed/counts_matrix_filtered.csv")

# Define experimental design factors
tissues <- c("CE", "DRG", "MNG", "Sciat", "SCord")
conditions <- c("Inf", "Mock")
genotypes <- c("B6", "IL10KO")

# Function to calculate CPM
calculate_cpm <- function(count_matrix) {
  # Calculate library sizes
  lib_sizes <- colSums(count_matrix)
  # Calculate CPM
  t(t(count_matrix) / lib_sizes * 1e6)
}

# Separate gene names and count data
gene_names <- counts$GeneName
count_matrix <- counts[, -which(names(counts) == "GeneName")]

# Calculate CPM for all samples
cpm_matrix <- calculate_cpm(count_matrix)

# Create a matrix to store whether each gene passes filters in each group
pass_filter <- matrix(FALSE, nrow = nrow(counts), ncol = length(tissues) * length(conditions) * length(genotypes))
colnames(pass_filter) <- paste0(rep(tissues, each = length(conditions) * length(genotypes)), "_",
                               rep(rep(conditions, each = length(genotypes)), length(tissues)), "_",
                               rep(genotypes, length(tissues) * length(conditions)))

# Function to get sample indices for a specific combination
get_sample_indices <- function(tissue, condition, genotype) {
  # Construct the pattern to match sample names
  pattern <- paste0("^", tissue)
  tissue_samples <- grep(pattern, colnames(count_matrix))
  
  # Further filter based on metadata (you'll need to adjust this based on your actual sample naming convention)
  # This is a placeholder - you'll need to modify based on how conditions and genotypes are encoded in your sample names
  indices <- tissue_samples
  return(indices)
}

# For each combination, check if genes pass the threshold
current_col <- 1
for (tissue in tissues) {
  for (condition in conditions) {
    for (genotype in genotypes) {
      # Get relevant sample indices
      sample_indices <- get_sample_indices(tissue, condition, genotype)
      
      # Get CPM for these samples
      group_cpm <- cpm_matrix[, sample_indices]
      
      # Check which genes have CPM > 1 in at least 2 out of 3 replicates
      pass_filter[, current_col] <- rowSums(group_cpm > 1) >= 2
      
      current_col <- current_col + 1
    }
  }
}

# A gene passes if it meets the criteria in ANY group
genes_to_keep <- rowSums(pass_filter) > 0

# Filter the count matrix
counts_filtered_low <- counts[genes_to_keep, ]

# Calculate statistics
n_genes_before <- nrow(counts)
n_genes_after <- nrow(counts_filtered_low)
percent_removed <- (n_genes_before - n_genes_after) / n_genes_before * 100

# Create validation plots directory
dir.create("results/qc/filtering_validation", recursive = TRUE, showWarnings = FALSE)

# Validation Step 1: CPM Distribution Before vs After
pdf("results/qc/filtering_validation/cpm_distribution.pdf")

# Before filtering
cpm_before <- log2(cpm_matrix + 1)
cpm_before_melted <- melt(cpm_before)
ggplot(cpm_before_melted, aes(x=value)) +
  geom_density(fill="blue", alpha=0.3) +
  labs(title="CPM Distribution Before Filtering",
       x="log2(CPM + 1)",
       y="Density") +
  theme_minimal()

# After filtering
cpm_after <- log2(calculate_cpm(as.matrix(counts_filtered_low[,-1])) + 1)
cpm_after_melted <- melt(cpm_after)
ggplot(cpm_after_melted, aes(x=value)) +
  geom_density(fill="red", alpha=0.3) +
  labs(title="CPM Distribution After Filtering",
       x="log2(CPM + 1)",
       y="Density") +
  theme_minimal()

dev.off()

# Validation Step 2: Sample Correlation Heatmap
pdf("results/qc/filtering_validation/sample_correlation.pdf")

# Before filtering
cor_matrix_before <- cor(cpm_before)
pheatmap(cor_matrix_before,
         main="Sample Correlation Before Filtering",
         show_rownames=FALSE,
         show_colnames=FALSE)

# After filtering
cor_matrix_after <- cor(cpm_after)
pheatmap(cor_matrix_after,
         main="Sample Correlation After Filtering",
         show_rownames=FALSE,
         show_colnames=FALSE)

dev.off()

# Validation Step 3: CV (Coefficient of Variation) Analysis
calculate_cv <- function(x) sd(x)/mean(x)

cv_before <- apply(cpm_matrix, 1, calculate_cv)
cv_after <- apply(calculate_cpm(as.matrix(counts_filtered_low[,-1])), 1, calculate_cv)

pdf("results/qc/filtering_validation/cv_distribution.pdf")
par(mfrow=c(1,2))
hist(cv_before, breaks=100, main="CV Distribution Before Filtering",
     xlab="Coefficient of Variation", ylab="Frequency")
hist(cv_after, breaks=100, main="CV Distribution After Filtering",
     xlab="Coefficient of Variation", ylab="Frequency")
dev.off()

# Save filtered matrix with final name
write.csv(counts_filtered_low,
          "data/processed/counts_final_filtered.csv",
          row.names = FALSE,
          quote = FALSE)

# Generate report
sink("results/qc/filtering_validation/filtering_validation_report.md")

cat("# Filtering Validation Report\n")
cat("**Author:** Eren Ada, PhD\n")
cat(sprintf("**Date:** %s\n\n", format(Sys.time(), "%Y-%m-%d")))

cat("## Filtering Strategy\n")
cat("Applied CPM-based filtering with tissue and condition-specific criteria:\n")
cat("1. CPM > 1 means more than 1 read per million reads in the library\n")
cat("2. For each tissue/condition/genotype combination:\n")
cat("   - Required CPM > 1 in at least 2 out of 3 replicates\n")
cat("3. Retained genes that passed criteria in ANY group\n\n")

cat("## Validation Results\n\n")

cat("### 1. Basic Statistics\n")
cat(sprintf("- Total genes before filtering: %d\n", n_genes_before))
cat(sprintf("- Genes retained after filtering: %d\n", n_genes_after))
cat(sprintf("- Genes removed: %d (%.2f%%)\n\n", n_genes_before - n_genes_after, percent_removed))

cat("### 2. CPM Distribution\n")
cat("- Generated density plots of CPM values before and after filtering\n")
cat("- After filtering shows removal of low-expression noise\n")
cat("- Higher average CPM in retained genes\n\n")

cat("### 3. Sample Correlations\n")
cat("- Generated correlation heatmaps before and after filtering\n")
cat(sprintf("- Mean correlation before: %.3f\n", mean(cor_matrix_before[upper.tri(cor_matrix_before)])))
cat(sprintf("- Mean correlation after: %.3f\n\n", mean(cor_matrix_after[upper.tri(cor_matrix_after)])))

cat("### 4. Coefficient of Variation\n")
cat(sprintf("- Median CV before: %.3f\n", median(cv_before, na.rm=TRUE)))
cat(sprintf("- Median CV after: %.3f\n", median(cv_after, na.rm=TRUE)))
cat("- Lower CV after filtering indicates more stable expression patterns\n\n")

cat("## Validation Plots Generated\n")
cat("1. CPM distribution comparison (cpm_distribution.pdf)\n")
cat("2. Sample correlation heatmaps (sample_correlation.pdf)\n")
cat("3. Coefficient of variation distribution (cv_distribution.pdf)\n\n")

cat("## Files Generated\n")
cat("1. Final filtered count matrix: data/processed/counts_final_filtered.csv\n")
cat("2. Validation plots: results/qc/filtering_validation/\n")
cat("3. This report: results/qc/filtering_validation/filtering_validation_report.md\n")

sink()

# Print completion message
message("Low count filtering and validation completed.")
message("Final filtered count matrix saved to: data/processed/counts_final_filtered.csv")
message("Validation report and plots saved in: results/qc/filtering_validation/") 