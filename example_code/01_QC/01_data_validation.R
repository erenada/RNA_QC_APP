# Data Validation Script for RNA-seq Analysis
# Author: Eren Ada, PhD
# Date: April 2024
# Purpose: Validate raw count data and metadata for DEG analysis

# Load required libraries
library(tidyverse)
library(here)
library(DESeq2)
library(edgeR)  # for handling duplicate genes
library(skimr)  # for data summaries
library(logger) # for logging

# Set working directory to the project root
setwd("/Users/eren/Desktop/HMS/EricLoo/eric_rnaseq")

# Set up logging
dir.create("logs", showWarnings = FALSE)
log_appender(appender_file("logs/data_validation.log"))
log_info("Starting data validation")

# Create necessary directories
dir.create("results/qc", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("docs/methods", recursive = TRUE, showWarnings = FALSE)

#' Function to check for and handle duplicate gene names
#' @param counts_df DataFrame with gene names as rownames
#' @return List containing processed dataframe and duplicate genes info
handle_duplicates <- function(counts_df) {
  # Get the gene names
  gene_names <- counts_df$GeneName
  
  # Find duplicates
  duplicate_genes <- gene_names[duplicated(gene_names)]
  
  if (length(duplicate_genes) > 0) {
    log_warn(sprintf("Found %d duplicate gene names", length(duplicate_genes)))
    
    # Create a data frame with gene names and row numbers
    gene_info <- data.frame(
      gene_name = gene_names,
      row_num = 1:length(gene_names)
    )
    
    # For each duplicate, append a suffix
    for (gene in unique(duplicate_genes)) {
      dups <- gene_info$row_num[gene_info$gene_name == gene]
      for (i in seq_along(dups)) {
        gene_info$gene_name[dups[i]] <- sprintf("%s_%d", gene, i)
      }
    }
    
    # Update the gene names in the counts dataframe
    counts_df$GeneName <- gene_info$gene_name
    
    log_info("Duplicate genes have been renamed with numeric suffixes")
    
    # Write out the duplicate gene information
    sink("results/qc/duplicate_genes.txt")
    cat("Duplicate Genes Report\n")
    cat("=====================\n\n")
    cat("The following genes had duplicate entries:\n\n")
    for (gene in unique(duplicate_genes)) {
      cat(sprintf("- %s\n", gene))
    }
    sink()
  }
  
  # Now we can safely set the row names
  counts_df_final <- counts_df %>%
    column_to_rownames("GeneName")
  
  return(list(
    counts = counts_df_final,
    duplicates = duplicate_genes
  ))
}

#' Function to validate sample matching
#' @param counts_df Count matrix
#' @param metadata_df Metadata dataframe
#' @return Boolean indicating if samples match
validate_samples <- function(counts_df, metadata_df) {
  count_samples <- colnames(counts_df)
  meta_samples <- metadata_df$Samples  # Updated to match your metadata column name
  
  # Check if all samples in counts exist in metadata
  counts_in_meta <- all(count_samples %in% meta_samples)
  # Check if all samples in metadata exist in counts
  meta_in_counts <- all(meta_samples %in% count_samples)
  
  if (!counts_in_meta || !meta_in_counts) {
    log_error("Sample mismatch between count matrix and metadata!")
    return(FALSE)
  }
  
  return(TRUE)
}

# Main execution
tryCatch({
  # Read data
  log_info("Reading count matrix and metadata")
  counts <- read.csv("data/raw/all_rna_counts.csv")
  metadata <- read.csv("data/raw/all_rna_metadata.csv")
  
  # Basic data structure checks
  log_info("Checking data structure")
  log_info(sprintf("Count matrix dimensions: %d genes x %d samples", 
                  nrow(counts)-1, ncol(counts)-1))  # -1 for the gene name column
  log_info(sprintf("Metadata dimensions: %d samples x %d columns", 
                  nrow(metadata), ncol(metadata)))
  
  # Handle duplicate genes
  processed_data <- handle_duplicates(counts)
  counts <- processed_data$counts
  duplicate_genes <- processed_data$duplicates
  
  # Validate sample matching
  if (!validate_samples(counts, metadata)) {
    stop("Sample validation failed")
  }
  
  # Basic count matrix QC
  log_info("Performing basic QC checks")
  
  # Check for zero-count genes
  zero_counts <- rowSums(counts == 0) == ncol(counts)
  log_info(sprintf("Found %d genes with zero counts across all samples", 
                  sum(zero_counts)))
  
  # Library sizes
  lib_sizes <- colSums(as.matrix(counts))
  log_info(sprintf("Library size range: %.0f to %.0f", 
                  min(lib_sizes), max(lib_sizes)))
  
  # Generate QC summary
  sink("results/qc/data_validation_summary.txt")
  cat("RNA-seq Data Validation Summary\n")
  cat("==============================\n\n")
  
  cat("1. Data Dimensions\n")
  cat(sprintf("   - Genes: %d\n", nrow(counts)))
  cat(sprintf("   - Samples: %d\n", ncol(counts)))
  cat(sprintf("   - Zero-count genes: %d\n", sum(zero_counts)))
  
  cat("\n2. Library Sizes\n")
  cat(sprintf("   - Minimum: %.0f\n", min(lib_sizes)))
  cat(sprintf("   - Maximum: %.0f\n", max(lib_sizes)))
  cat(sprintf("   - Median: %.0f\n", median(lib_sizes)))
  
  cat("\n3. Sample Groups\n")
  cat("By Tissue:\n")
  print(table(metadata$Tissue))
  cat("\nBy Treatment:\n")
  print(table(metadata$Treatment))
  cat("\nBy Genotype:\n")
  print(table(metadata$Genotype))
  
  # Additional QC metrics
  cat("\n4. Count Distribution\n")
  count_summary <- summary(as.numeric(as.matrix(counts)))
  print(count_summary)
  
  cat("\n5. Sample Information\n")
  cat("Samples per group combination:\n")
  print(table(metadata$Tissue, metadata$Treatment, metadata$Genotype))
  
  sink()
  
  # Save processed data
  log_info("Saving processed data")
  saveRDS(counts, "data/processed/counts_processed.rds")
  saveRDS(metadata, "data/processed/metadata_processed.rds")
  
  # Create a methods documentation entry
  sink("docs/methods/data_processing.md", append = TRUE)
  cat("\n## Data Validation and Processing\n")
  cat(sprintf("Date: %s\n\n", format(Sys.time(), "%Y-%m-%d")))
  cat("### Data Processing Steps:\n")
  cat("1. Loaded raw count matrix and metadata\n")
  cat(sprintf("2. Processed %d genes across %d samples\n", nrow(counts), ncol(counts)))
  if (length(duplicate_genes) > 0) {
    cat(sprintf("3. Renamed %d duplicate gene entries with numeric suffixes\n", 
                length(duplicate_genes)))
  }
  cat(sprintf("4. Identified %d genes with zero counts across all samples\n", 
              sum(zero_counts)))
  cat("\n### Sample Groups:\n")
  cat("- Tissues:", paste(unique(metadata$Tissue), collapse = ", "), "\n")
  cat("- Treatments:", paste(unique(metadata$Treatment), collapse = ", "), "\n")
  cat("- Genotypes:", paste(unique(metadata$Genotype), collapse = ", "), "\n")
  sink()
  
  log_info("Data validation completed successfully")
  
}, error = function(e) {
  log_error(sprintf("Error in data validation: %s", e$message))
  stop(e)
}) 