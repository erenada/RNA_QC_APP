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
dir.create("results/figures/qc", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables/qc", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

# Read sample metadata
metadata <- read.csv("data/metadata/larissa_samples.csv")

# Print column names in metadata to check
cat("Columns in metadata:", paste(colnames(metadata), collapse=", "), "\n")

# Read count matrix
counts <- read.csv("data/raw/larissa_counts.csv")

# Print basic information about the data
cat("Count matrix dimensions:", dim(counts), "\n")
cat("Number of samples:", ncol(counts) - 1, "\n")  # Subtract 1 for gene column
cat("Number of genes:", nrow(counts), "\n")
cat("Number of samples in metadata:", nrow(metadata), "\n")

# Extract sample names from count matrix
count_samples <- colnames(counts)[-1]  # Exclude the first column (gene names)

# Extract sample names from metadata
metadata_samples <- metadata$Sample

# Check for missing samples
missing_in_counts <- setdiff(metadata$ID, count_samples)
missing_in_metadata <- setdiff(count_samples, metadata$ID)

# Print validation results
cat("\nValidation Results:\n")
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

# Print unique tissues in metadata
cat("\nUnique tissues in metadata:", paste(unique(metadata$Tissue), collapse=", "), "\n")

# Function to prepare data for a specific tissue
prepare_tissue_data <- function(tissue_name) {
    cat("\nProcessing tissue:", tissue_name, "\n")
    
    # Filter metadata for the specific tissue
    tissue_metadata <- metadata[metadata$Tissue == tissue_name, ]
    cat("Number of samples for this tissue:", nrow(tissue_metadata), "\n")
    
    # Get sample IDs for this tissue
    tissue_samples <- tissue_metadata$ID
    cat("Sample IDs for this tissue:", paste(head(tissue_samples, 3), collapse=", "), "...\n")
    
    # Check if all tissue samples exist in the count matrix
    missing_samples <- setdiff(tissue_samples, colnames(counts)[-1])
    if(length(missing_samples) > 0) {
        cat("WARNING: Some tissue samples are missing from the count matrix:\n")
        print(missing_samples)
        tissue_samples <- intersect(tissue_samples, colnames(counts)[-1])
        cat("Proceeding with", length(tissue_samples), "samples that are present in both metadata and counts.\n")
    }
    
    # Create a subset of the count matrix for this tissue
    tissue_counts <- counts[, c("Gene", tissue_samples)]
    cat("Dimensions of tissue count matrix:", dim(tissue_counts), "\n")
    
    # Print information about the tissue data
    cat("Tissue:", tissue_name, "\n")
    cat("Number of samples:", nrow(tissue_metadata), "\n")
    
    # Create a DESeq2 object for this tissue
    # First, set up the count matrix with gene names as row names
    count_matrix <- as.matrix(tissue_counts[, -1])
    
    # Check if the counts are integers
    if(!all(count_matrix == round(count_matrix))) {
        cat("WARNING: Count matrix contains non-integer values (likely from RSEM).\n")
        cat("Rounding counts to integers for DESeq2 compatibility...\n")
        
        # Round the counts to integers
        count_matrix <- round(count_matrix)
    }
    
    rownames(count_matrix) <- tissue_counts$Gene
    colnames(count_matrix) <- tissue_samples
    
    # Ensure tissue_metadata rows match count_matrix columns
    tissue_metadata <- tissue_metadata[match(colnames(count_matrix), tissue_metadata$ID), ]
    
    # Check if the sample IDs match between count matrix and metadata
    if(!all(colnames(count_matrix) == tissue_metadata$ID)) {
        cat("ERROR: Sample IDs don't match between count matrix and metadata!\n")
        cat("Count matrix column names:", paste(head(colnames(count_matrix), 3), collapse=", "), "...\n")
        cat("Metadata IDs:", paste(head(tissue_metadata$ID, 3), collapse=", "), "...\n")
        stop("Sample ID mismatch")
    }
    
    # Check if required columns exist in metadata
    required_cols <- c("Sex", "Day", "Treatment")
    missing_cols <- setdiff(required_cols, colnames(tissue_metadata))
    if(length(missing_cols) > 0) {
        cat("ERROR: Missing required columns in metadata:", paste(missing_cols, collapse=", "), "\n")
        cat("Available columns:", paste(colnames(tissue_metadata), collapse=", "), "\n")
        stop("Missing required columns in metadata")
    }
    
    # Ensure factors are properly set
    tissue_metadata$Sex <- as.factor(tissue_metadata$Sex)
    tissue_metadata$Day <- as.factor(tissue_metadata$Day)
    tissue_metadata$Treatment <- as.factor(tissue_metadata$Treatment)
    
    # Print factor levels
    cat("Sex levels:", paste(levels(tissue_metadata$Sex), collapse=", "), "\n")
    cat("Day levels:", paste(levels(tissue_metadata$Day), collapse=", "), "\n")
    cat("Treatment levels:", paste(levels(tissue_metadata$Treatment), collapse=", "), "\n")
    
    # Create the DESeq2 object
    cat("Creating DESeq2 object...\n")
    tryCatch({
        dds <- DESeqDataSetFromMatrix(
            countData = count_matrix,
            colData = tissue_metadata,
            design = ~ Sex + Day + Treatment
        )
        cat("DESeq2 object created successfully.\n")
    }, error = function(e) {
        cat("ERROR creating DESeq2 object:", conditionMessage(e), "\n")
        print(head(tissue_metadata))
        stop("Failed to create DESeq2 object")
    })
    
    # Pre-filter low count genes
    cat("Filtering low count genes...\n")
    keep <- rowSums(counts(dds) >= 10) >= 5
    dds <- dds[keep, ]
    
    # Print filtering results
    cat("Genes after filtering (≥10 counts in ≥5 samples):", nrow(dds), "\n")
    cat("Percentage of genes retained:", round(nrow(dds) / nrow(counts) * 100, 1), "%\n")
    
    # Save the DESeq2 object
    file_path <- paste0("data/processed/dds_", tolower(tissue_name), ".rds")
    saveRDS(dds, file_path)
    cat("DESeq2 object saved to", file_path, "\n")
    
    return(dds)
}

# Process each tissue type
tissues <- unique(metadata$Tissue)
dds_list <- list()

cat("\nProcessing", length(tissues), "tissue types:", paste(tissues, collapse=", "), "\n")

for(tissue in tissues) {
    cat("\n==== Processing tissue:", tissue, "====\n")
    tryCatch({
        dds_list[[tissue]] <- prepare_tissue_data(tissue)
        cat("Successfully processed tissue:", tissue, "\n")
    }, error = function(e) {
        cat("ERROR processing tissue", tissue, ":", conditionMessage(e), "\n")
    })
}

# Save validation results to a log file
sink("logs/data_validation_results.txt")
cat("Data Validation Results\n")
cat("======================\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Count matrix dimensions:", dim(counts), "\n")
cat("Number of samples:", ncol(counts) - 1, "\n")
cat("Number of genes:", nrow(counts), "\n")
cat("Number of samples in metadata:", nrow(metadata), "\n\n")

cat("Validation Results:\n")
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

for(tissue in tissues) {
    if(tissue %in% names(dds_list)) {
        cat("\nTissue:", tissue, "\n")
        cat("Number of samples:", sum(metadata$Tissue == tissue), "\n")
        cat("Genes after filtering (≥10 counts in ≥5 samples):", nrow(dds_list[[tissue]]), "\n")
        cat("Percentage of genes retained:", round(nrow(dds_list[[tissue]]) / nrow(counts) * 100, 1), "%\n")
    } else {
        cat("\nTissue:", tissue, "- Processing failed\n")
    }
}
sink()

cat("\nData validation complete. Results saved to logs/data_validation_results.txt\n")
cat("DESeq2 objects saved to data/processed/\n") 