# Duplicate Gene Analysis Script
# Author: Eren Ada, PhD
# Date: April 2024
# Purpose: Analyze duplicate genes in detail and create deduplicated matrix

# Load required libraries
library(tidyverse)
library(here)
library(DESeq2)
library(pheatmap)

# Set working directory to the project root
setwd("/Users/eren/Desktop/HMS/EricLoo/eric_rnaseq")

# Read the original data
counts_raw <- read.csv("data/raw/all_rna_counts.csv")

# Find duplicate genes
gene_names <- counts_raw$GeneName
duplicate_genes <- unique(gene_names[duplicated(gene_names)])
n_duplicates <- length(duplicate_genes)

# Create output directory
dir.create("results/qc/duplicate_analysis", recursive = TRUE, showWarnings = FALSE)

# Function to analyze a duplicate gene
analyze_duplicate <- function(gene_name, counts_df) {
  # Get rows for this gene
  gene_rows <- which(counts_df$GeneName == gene_name)
  gene_data <- counts_df[gene_rows, ]
  
  # Calculate summary statistics for each duplicate
  summaries <- list()
  for (i in seq_along(gene_rows)) {
    counts <- as.numeric(gene_data[i, -1])  # exclude gene name column
    summaries[[i]] <- c(
      mean = mean(counts),
      sd = sd(counts),
      median = median(counts),
      zeros = sum(counts == 0),
      total = sum(counts)
    )
  }
  
  return(list(
    rows = gene_rows,
    data = gene_data,
    summaries = summaries
  ))
}

# Function to handle duplicate genes and create unique identifiers
create_unique_gene_names <- function(counts_df) {
  # Create a copy of the input data
  counts_new <- counts_df
  
  # Get gene names
  gene_names <- counts_df$GeneName
  
  # Find duplicates
  duplicated_indices <- duplicated(gene_names)
  duplicate_genes <- unique(gene_names[duplicated(gene_names)])
  
  # For each duplicate gene
  for (gene in duplicate_genes) {
    # Find all occurrences
    indices <- which(gene_names == gene)
    
    # Add suffix to create unique names
    for (i in seq_along(indices)) {
      gene_names[indices[i]] <- sprintf("%s_%d", gene, i)
    }
  }
  
  # Update gene names in the dataframe
  counts_new$GeneName <- gene_names
  
  # Ensure count columns are integers
  count_cols <- 2:ncol(counts_new)  # Assuming first column is GeneName
  counts_new[,count_cols] <- lapply(counts_new[,count_cols], function(x) as.integer(round(as.numeric(x))))
  
  return(counts_new)
}

# Create deduplicated count matrix
counts_unique <- create_unique_gene_names(counts_raw)

# Save deduplicated matrix with integer counts
options(scipen = 999)  # Prevent scientific notation
write.csv(counts_unique, "data/processed/counts_matrix_unique_genes.csv", row.names = FALSE)

# Open the main report file
sink("results/qc/duplicate_analysis/duplicate_genes_detailed.md")

cat("# Detailed Analysis of Duplicate Genes\n")
cat(sprintf("**Date:** %s\n\n", format(Sys.time(), "%Y-%m-%d")))

cat(sprintf("Found %d genes with duplicate entries.\n\n", n_duplicates))

# Analyze each duplicate gene
all_summaries <- list()
for (gene in duplicate_genes) {
  cat(sprintf("\n## %s\n\n", gene))
  
  analysis <- analyze_duplicate(gene, counts_raw)
  all_summaries[[gene]] <- analysis$summaries
  
  cat("### Expression Summary\n\n")
  cat("| Version | Mean | SD | Median | Zero Counts | Total Counts |\n")
  cat("|---------|------|----|---------| ------------|---------------|\n")
  
  for (i in seq_along(analysis$summaries)) {
    stats <- analysis$summaries[[i]]
    cat(sprintf("| %s_%d | %.2f | %.2f | %.2f | %d | %.0f |\n",
                gene, i, stats["mean"], stats["sd"], 
                stats["median"], stats["zeros"], stats["total"]))
  }
  
  # Calculate correlation between duplicates if there are exactly 2 versions
  if (length(analysis$summaries) == 2) {
    counts1 <- as.numeric(analysis$data[1, -1])
    counts2 <- as.numeric(analysis$data[2, -1])
    cor_val <- cor(counts1, counts2)
    cat(sprintf("\nCorrelation between duplicates: %.3f\n", cor_val))
  }
  
  cat("\n### Possible Reasons for Duplication\n")
  # Check for common duplication patterns
  if (all(sapply(analysis$summaries, function(x) x["total"]) > 0)) {
    cat("- Both/all versions show expression, suggesting possible gene duplication or alternative transcripts\n")
  }
  if (any(sapply(analysis$summaries, function(x) x["total"]) == 0)) {
    cat("- Some versions show no expression, suggesting possible pseudogenes or annotation errors\n")
  }
  if (length(grep("Rik|Gm|[0-9]$", gene)) > 0) {
    cat("- Gene name suggests predicted gene or RIKEN sequence, which often have uncertain annotations\n")
  }
}

# Overall summary
cat("\n## Overall Summary\n\n")
cat(sprintf("- Total number of duplicate genes: %d\n", n_duplicates))

# Categorize duplicates
zero_exp_pairs <- 0
both_exp_pairs <- 0
mixed_exp_pairs <- 0

for (gene in names(all_summaries)) {
  totals <- sapply(all_summaries[[gene]], function(x) x["total"])
  if (all(totals == 0)) zero_exp_pairs <- zero_exp_pairs + 1
  else if (all(totals > 0)) both_exp_pairs <- both_exp_pairs + 1
  else mixed_exp_pairs <- mixed_exp_pairs + 1
}

cat("\n### Expression Patterns\n")
cat(sprintf("- Genes where all duplicates show expression: %d\n", both_exp_pairs))
cat(sprintf("- Genes where no duplicates show expression: %d\n", zero_exp_pairs))
cat(sprintf("- Genes with mixed expression patterns: %d\n", mixed_exp_pairs))

cat("\n### Common Patterns\n")
cat(sprintf("- RIKEN genes (Rik): %d\n", sum(grepl("Rik", duplicate_genes))))
cat(sprintf("- Predicted genes (Gm): %d\n", sum(grepl("Gm", duplicate_genes))))
cat(sprintf("- Pseudogenes (ps): %d\n", sum(grepl("ps", duplicate_genes))))

cat("\n## Handling Strategy\n")
cat("Duplicate genes were handled by:\n")
cat("1. Identifying all duplicate gene entries\n")
cat("2. Appending numeric suffixes to create unique identifiers (e.g., Gene_1, Gene_2)\n")
cat("3. Maintaining all versions separately for downstream analysis\n")
cat("4. Recording duplicate information for reference during biological interpretation\n")
cat(sprintf("\nDeduplicated count matrix saved as: data/processed/counts_matrix_unique_genes.csv\n"))

sink()

# Save the detailed analysis results
saveRDS(all_summaries, "results/qc/duplicate_analysis/duplicate_gene_summaries.rds")

# Create a log entry
sink("logs/duplicate_analysis.log")
cat(sprintf("Duplicate gene analysis completed: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Found %d duplicate genes\n", n_duplicates))
cat(sprintf("Results saved in results/qc/duplicate_analysis/\n"))
cat(sprintf("Deduplicated count matrix saved as: data/processed/counts_matrix_unique_genes.csv\n"))
sink() 