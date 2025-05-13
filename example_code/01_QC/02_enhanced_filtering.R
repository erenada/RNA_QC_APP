# scripts/02_enhanced_filtering.R
#
# Author: Eren Ada, PhD
# Date: 05/08/2025
#
# Description: 
# This script performs group-aware pre-filtering on a DESeqDataSet object.
# The goal is to retain genes that show sufficient expression in at least
# one experimental group, improving the power and reliability of downstream
# differential expression analysis.
#
# Input files:
# - DESeqDataSet object (results/dea/deseq2_objects/dds_unprocessed.rds)
#
# Output files:
# - Filtered DESeqDataSet object (results/dea/deseq2_objects/dds_filtered.rds)
# - Filtered count matrix as CSV (results/dea/counts/filtered_counts_matrix.csv)
# - QC report for filtering (logs/02_filtering_qc_report.csv)
# - Filtering log (logs/02_filtering_log.txt)

# ========================= #
# 0. Setup and Preparation #
# ========================= #

# Load required packages
suppressPackageStartupMessages({
  # Core analysis packages
  library(DESeq2)
  # Fix for DESeq2 on R 4.5.0 - "superclass 'ExpData' not defined" error
  setOldClass("ExpData")
  
  # Data manipulation packages
  library(tidyverse)
  library(dplyr)
  
  # Visualization packages
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(vsn)  # For variance stabilization plots
  
  # Project organization
  library(here)
})

# Initialize here
here::i_am("scripts/02_enhanced_filtering.R")

# Set up file paths
paths <- list(
  dds_unprocessed = here("results", "dea", "deseq2_objects", "dds_unprocessed.rds"),
  dds_filtered_output = here("results", "dea", "deseq2_objects", "dds_filtered.rds"),
  counts_filtered_output = here("results", "dea", "counts", "filtered_counts_matrix.csv"),
  log_file = here("logs", "02_filtering_log.txt"),
  qc_report_file = here("logs", "02_filtering_qc_report.csv")
)

# Start logging
dir.create(dirname(paths$log_file), recursive = TRUE, showWarnings = FALSE)
sink(paths$log_file, append = FALSE, split = TRUE)
cat("=== Gene Filtering Log ===\n")
cat("Script started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# Record Git commit hash if available
git_hash <- try(system("git rev-parse HEAD", intern = TRUE), silent = TRUE)
if (!inherits(git_hash, "try-error")) {
  cat("Git commit:", git_hash, "\n")
}

# Record R version
cat("\nR version:", R.version.string, "\n")

# ==================================== #
# 1. Load Unprocessed DESeqDataSet     #
# ==================================== #

cat("\nLoading unprocessed DESeqDataSet object...\n")
if (!file.exists(paths$dds_unprocessed)) {
  stop("Input DESeqDataSet file not found: ", paths$dds_unprocessed)
}

dds <- readRDS(paths$dds_unprocessed)

cat("Initial number of genes:", nrow(dds), "\n")
cat("Initial number of samples:", ncol(dds), "\n")
cat("\nSample groups:\n")
print(table(dds$group))

# ================================================= #
# 2. Define Filtering Strategy & Parameters          #
# ================================================= #

# Determine the smallest experimental group size
smallest_group_size <- min(table(dds$group))
cat("\nSmallest experimental group size (N):", smallest_group_size, "\n")

# Set filtering parameters
min_count_per_sample <- 10  # Minimum count threshold per gene in a sample
num_samples_threshold <- smallest_group_size  # Number of samples that must meet the threshold

cat("Minimum count per sample (X):", min_count_per_sample, "\n")
cat("Minimum samples meeting threshold (Y):", num_samples_threshold, "\n")

# ==================================== #
# 3. Apply Filter                      #
# ==================================== #

cat("\nApplying gene filter...\n")

# Get count matrix
counts_matrix <- counts(dds)

# Identify genes to keep (genes with >= min_count_per_sample in >= num_samples_threshold samples)
genes_to_keep <- rowSums(counts_matrix >= min_count_per_sample) >= num_samples_threshold

# Calculate filtering statistics
num_genes_initial <- nrow(dds)
num_genes_kept <- sum(genes_to_keep)
num_genes_removed <- num_genes_initial - num_genes_kept
percent_kept <- (num_genes_kept / num_genes_initial) * 100

# Print filtering statistics
cat("\nFiltering Statistics:\n")
cat("Number of genes before filtering:", num_genes_initial, "\n")
cat("Number of genes passing filter:", num_genes_kept, "\n")
cat("Number of genes removed:", num_genes_removed, "\n")
cat("Percentage of genes retained:", sprintf("%.2f%%", percent_kept), "\n")

# Subset the DESeqDataSet
dds_filtered <- dds[genes_to_keep, ]

# ============================================ #
# 4. Save Filtered DESeqDataSet & QC Report   #
# ============================================ #

# Save filtered DESeqDataSet
cat("\nSaving filtered DESeqDataSet object...\n")
dir.create(dirname(paths$dds_filtered_output), recursive = TRUE, showWarnings = FALSE)
saveRDS(dds_filtered, file = paths$dds_filtered_output)
cat("Filtered DESeqDataSet saved to:", paths$dds_filtered_output, "\n")

# Save filtered count matrix as CSV
cat("\nSaving filtered count matrix as CSV...\n")
dir.create(dirname(paths$counts_filtered_output), recursive = TRUE, showWarnings = FALSE)
filtered_counts <- counts(dds_filtered)
# Add gene names as a column instead of rownames for better compatibility
filtered_counts_df <- as.data.frame(filtered_counts) %>%
  rownames_to_column("gene_id")
write.csv(filtered_counts_df, file = paths$counts_filtered_output, row.names = FALSE)
cat("Filtered count matrix saved to:", paths$counts_filtered_output, "\n")

# Generate QC report
qc_filtering_summary <- data.frame(
  Metric = c(
    "Initial number of genes",
    "Smallest group size (N)",
    "Min count threshold per sample (X)",
    "Min samples meeting threshold (Y)",
    "Genes passing filter",
    "Genes removed",
    "Percentage of genes retained"
  ),
  Value = c(
    num_genes_initial,
    smallest_group_size,
    min_count_per_sample,
    num_samples_threshold,
    num_genes_kept,
    num_genes_removed,
    sprintf("%.2f%%", percent_kept)
  )
)

# Save QC report
cat("\nGenerating QC report for filtering...\n")
dir.create(dirname(paths$qc_report_file), recursive = TRUE, showWarnings = FALSE)
write.csv(qc_filtering_summary, paths$qc_report_file, row.names = FALSE)
cat("Filtering QC report saved to:", paths$qc_report_file, "\n")

# Optional: Generate and save distribution plots
# You can add plotting code here if desired

# ======================== #
# 5. Finalize Logging     #
# ======================== #

cat("\nSession Info:\n")
print(sessionInfo())
cat("\nScript finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
sink() # Stop logging

cat("Filtering script completed successfully.\n")
cat("- Filtered DDS saved to:", paths$dds_filtered_output, "\n")
cat("- Filtered count matrix saved to:", paths$counts_filtered_output, "\n")
cat("- QC report saved to:", paths$qc_report_file, "\n")
cat("- Log file saved to:", paths$log_file, "\n") 