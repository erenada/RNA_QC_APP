# Filter Zero Counts Script
# Author: Eren Ada, PhD
# Date: April 2024
# Purpose: Remove genes with zero counts across all samples and analyze count distributions

# Load required libraries
library(tidyverse)
library(here)
library(ggplot2)
library(scales)

# Read the summed count matrix
counts <- read.csv("data/processed/counts_matrix_summed.csv")

# Separate gene names and count data
gene_names <- counts$GeneName
count_matrix <- counts[, -which(names(counts) == "GeneName")]

# Find genes with zero counts across all samples
zero_counts <- rowSums(count_matrix) == 0
n_zero_genes <- sum(zero_counts)

# Filter out zero-count genes
counts_filtered <- counts[!zero_counts, ]

# Calculate some statistics for the filtered dataset
n_genes_before <- nrow(counts)
n_genes_after <- nrow(counts_filtered)
percent_removed <- (n_genes_before - n_genes_after) / n_genes_before * 100

# Analyze count distribution
count_stats <- data.frame(
  mean_counts = rowMeans(count_matrix[!zero_counts, ]),
  max_counts = apply(count_matrix[!zero_counts, ], 1, max),
  non_zero_samples = rowSums(count_matrix[!zero_counts, ] > 0)
)

# Save filtered count matrix
write.csv(counts_filtered, 
          "data/processed/counts_matrix_filtered.csv",
          row.names = FALSE,
          quote = FALSE)

# Create visualizations directory
dir.create("results/qc/count_distribution", recursive = TRUE, showWarnings = FALSE)

# Create count distribution plots
pdf("results/qc/count_distribution/count_distribution_plots.pdf")

# 1. Distribution of mean counts
ggplot(count_stats, aes(x = mean_counts)) +
  geom_histogram(bins = 100) +
  scale_x_log10(labels = comma) +
  theme_minimal() +
  labs(title = "Distribution of Mean Counts per Gene",
       x = "Mean Counts (log10 scale)",
       y = "Number of Genes")

# 2. Distribution of non-zero samples
ggplot(count_stats, aes(x = non_zero_samples)) +
  geom_histogram(bins = 60) +
  theme_minimal() +
  labs(title = "Number of Samples with Non-zero Counts per Gene",
       x = "Number of Samples",
       y = "Number of Genes")

# 3. Scatter plot of mean vs max counts
ggplot(count_stats, aes(x = mean_counts, y = max_counts)) +
  geom_point(alpha = 0.1) +
  scale_x_log10(labels = comma) +
  scale_y_log10(labels = comma) +
  theme_minimal() +
  labs(title = "Mean vs Maximum Counts per Gene",
       x = "Mean Counts (log10 scale)",
       y = "Maximum Counts (log10 scale)")

dev.off()

# Generate report
sink("results/qc/count_distribution/filtering_report.md")

cat("# Count Filtering Analysis Report\n")
cat("**Author:** Eren Ada, PhD\n")
cat(sprintf("**Date:** %s\n\n", format(Sys.time(), "%Y-%m-%d")))

cat("## Summary Statistics\n")
cat(sprintf("- Total genes before filtering: %d\n", n_genes_before))
cat(sprintf("- Genes with zero counts across all samples: %d\n", n_zero_genes))
cat(sprintf("- Genes retained after filtering: %d\n", n_genes_after))
cat(sprintf("- Percentage of genes removed: %.2f%%\n\n", percent_removed))

cat("## Count Distribution Analysis\n")
cat("### Overall Statistics\n")
cat(sprintf("- Minimum non-zero mean count: %.2f\n", min(count_stats$mean_counts)))
cat(sprintf("- Median mean count: %.2f\n", median(count_stats$mean_counts)))
cat(sprintf("- Mean of mean counts: %.2f\n", mean(count_stats$mean_counts)))
cat(sprintf("- Maximum mean count: %.2f\n\n", max(count_stats$mean_counts)))

cat("### Sample Detection Rates\n")
detection_summary <- table(count_stats$non_zero_samples)
cat("Number of genes detected in X samples:\n")
for(i in sort(as.numeric(names(detection_summary)))) {
  cat(sprintf("- %d samples: %d genes\n", i, detection_summary[as.character(i)]))
}

cat("\n## Suggestions for Low Count Filtering\n")
cat("Several approaches can be considered for filtering low count genes:\n\n")

cat("1. **Minimum Count Threshold**\n")
cat("   - Filter genes with a minimum mean count across all samples\n")
cat("   - Common thresholds: 1, 5, or 10 counts per sample\n")
cat("   - Pros: Simple to implement and understand\n")
cat("   - Cons: Doesn't consider sample-specific variations\n\n")

cat("2. **Minimum Samples Threshold**\n")
cat("   - Keep genes with counts above X in at least Y samples\n")
cat("   - Example: >10 counts in at least 3 samples per condition\n")
cat("   - Pros: Considers sample-specific detection\n")
cat("   - Cons: May need different thresholds for different conditions\n\n")

cat("3. **CPM-based Filtering**\n")
cat("   - Filter based on Counts Per Million (CPM)\n")
cat("   - Common threshold: CPM > 1 in at least X samples\n")
cat("   - Pros: Accounts for library size differences\n")
cat("   - Cons: May still retain some noise in high-count samples\n\n")

cat("4. **Variance-based Filtering**\n")
cat("   - Remove genes with low variance across samples\n")
cat("   - Can be combined with mean count thresholds\n")
cat("   - Pros: Focuses on genes with biological variation\n")
cat("   - Cons: May remove consistently expressed housekeeping genes\n\n")

cat("## Files Generated\n")
cat("1. Filtered count matrix: data/processed/counts_matrix_filtered.csv\n")
cat("2. Count distribution plots: results/qc/count_distribution/count_distribution_plots.pdf\n")
cat("3. This report: results/qc/count_distribution/filtering_report.md\n\n")

cat("## Next Steps\n")
cat("1. Review the count distribution plots\n")
cat("2. Choose appropriate low count filtering strategy based on:\n")
cat("   - Experimental design\n")
cat("   - Number of replicates per condition\n")
cat("   - Expected expression patterns\n")
cat("3. Implement chosen filtering strategy\n")

sink()

# Print completion message
message("Zero count filtering completed.")
message("Filtered count matrix saved to: data/processed/counts_matrix_filtered.csv")
message("Report and plots saved in: results/qc/count_distribution/") 