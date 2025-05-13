# Data Validation and Quality Control Analysis
# This script validates sample metadata against count matrices and performs initial QC

# Load required libraries
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)

# Set seed for reproducibility
set.seed(42)

# Create output directories if they don't exist
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

# Read sample metadata
hp_metadata <- read.csv("data/Samples_HP_simplified.csv")
m_metadata <- read.csv("data/Samples_M_simplified.csv")

# Read count matrices
hp_counts <- read.csv("data/Counts_HP_simplified.csv")
m_counts <- read.csv("data/Counts_M_simplified.csv")

# Data Validation Functions
validate_samples <- function(counts, metadata, tissue) {
    # Get sample names
    count_samples <- colnames(counts)[-1]  # Assuming first column is gene names/ids
    metadata_samples <- metadata$Sample
    
    # Check for missing samples
    missing_in_counts <- setdiff(metadata_samples, count_samples)
    missing_in_metadata <- setdiff(count_samples, metadata_samples)
    
    # Print validation results
    cat(sprintf("\nValidation Results for %s:\n", tissue))
    cat("Samples in count matrix:", length(count_samples), "\n")
    cat("Samples in metadata:", length(metadata_samples), "\n")
    
    if(length(missing_in_counts) > 0) {
        cat("\nWARNING: Samples in metadata but missing in counts:\n")
        print(missing_in_counts)
    }
    
    if(length(missing_in_metadata) > 0) {
        cat("\nWARNING: Samples in counts but missing in metadata:\n")
        print(missing_in_metadata)
    }
    
    if(length(missing_in_counts) == 0 && length(missing_in_metadata) == 0) {
        cat("\nAll samples match between metadata and counts.\n")
    }
    
    # Return TRUE if validation passes
    return(length(missing_in_counts) == 0 && length(missing_in_metadata) == 0)
}

# Validate Hippocampus samples
hp_valid <- validate_samples(hp_counts, hp_metadata, "Hippocampus")

# Validate Meninges samples
m_valid <- validate_samples(m_counts, m_metadata, "Meninges")

# If validation passes, proceed with QC
if(hp_valid && m_valid) {
    cat("\nValidation passed for both tissues. Proceeding with QC analysis...\n")
    
    # Function to prepare DESeq2 object
    prepare_deseq <- function(counts, metadata, tissue) {
        # Ensure samples are in the same order
        metadata <- metadata[match(colnames(counts)[-1], metadata$Sample),]
        
        # Create DESeq2 object
        dds <- DESeqDataSetFromMatrix(
            countData = counts[,-1],  # Remove first column (gene names)
            colData = metadata,
            design = ~ Sex + Age + Condition
        )
        
        # Add gene names as row names
        rownames(dds) <- counts[,1]
        
        return(dds)
    }
    
    # Create DESeq2 objects
    dds_hp <- prepare_deseq(hp_counts, hp_metadata, "Hippocampus")
    dds_m <- prepare_deseq(m_counts, m_metadata, "Meninges")
    
    # Save validation results
    sink("logs/data_validation_results.txt")
    cat("Data Validation Results\n")
    cat("======================\n\n")
    cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
    validate_samples(hp_counts, hp_metadata, "Hippocampus")
    validate_samples(m_counts, m_metadata, "Meninges")
    sink()
    
    # Save DESeq2 objects for future use
    saveRDS(dds_hp, "data/processed/dds_hippocampus.rds")
    saveRDS(dds_m, "data/processed/dds_meninges.rds")
    
    cat("\nValidation results saved to logs/data_validation_results.txt")
    cat("\nDESeq2 objects saved to data/processed/")
} else {
    stop("Validation failed. Please check the sample names in count matrices and metadata.")
} 