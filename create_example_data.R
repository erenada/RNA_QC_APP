#!/usr/bin/env Rscript
# Script to create anonymized example data for RNA Processing App
# Author: Eren Ada, PhD
# Date: 6/5/2025

library(dplyr)
library(readr)

# Set seed for reproducibility
set.seed(42)

cat("Creating anonymized example data...\n")

# Read original data
samples_df <- read_csv("example_data/test_samples.csv")
cat("Original samples:", nrow(samples_df), "\n")

# Create anonymized metadata
anonymized_samples <- samples_df %>%
  mutate(
    # Anonymize strain names
    Strain = case_when(
      Strain == "B6" ~ "StrainA",
      Strain == "BALBc" ~ "StrainB",
      TRUE ~ Strain
    ),
    # Anonymize treatment names
    Treatment = case_when(
      Treatment == "CTRL" ~ "Control",
      Treatment == "FA" ~ "TreatmentX",
      TRUE ~ Treatment
    ),
    # Anonymize tissue names
    Tissue = case_when(
      Tissue == "NG" ~ "Tissue1",
      Tissue == "T-DRG" ~ "Tissue2", 
      Tissue == "L-DRG" ~ "Tissue3",
      TRUE ~ Tissue
    ),
    # Create new anonymized sample names
    Sample_New = paste0("Sample_", sprintf("%02d", row_number()))
  )

# Create mapping of old to new sample names
sample_mapping <- setNames(anonymized_samples$Sample_New, anonymized_samples$Sample)

# Select a subset of samples for the example (12 samples total - 2 per condition)
selected_samples <- anonymized_samples %>%
  group_by(Strain, Treatment, Tissue) %>%
  slice_head(n = 2) %>%  # Take first 2 samples from each group
  ungroup()

cat("Selected samples for example:", nrow(selected_samples), "\n")

# Save anonymized metadata
final_metadata <- selected_samples %>%
  select(Sample = Sample_New, Strain, Treatment, Tissue)

write_csv(final_metadata, "example_data/example_metadata.csv")
cat("Saved anonymized metadata to example_data/example_metadata.csv\n")

# Now process the count matrix
cat("Processing count matrix...\n")

# Read header to get column names
header_line <- readLines("example_data/test_counts.csv", n = 1)
column_names <- strsplit(header_line, ",")[[1]]

# Get the sample columns (all except first "Gene" column)
sample_columns <- column_names[-1]
selected_sample_columns <- intersect(sample_columns, selected_samples$Sample)

cat("Found", length(selected_sample_columns), "matching sample columns\n")

# Read the full count matrix (this might take a moment)
cat("Reading count matrix (this may take a few moments)...\n")
counts_df <- read_csv("example_data/test_counts.csv")

# Select only the samples we want to keep
selected_columns <- c("Gene", selected_sample_columns)
counts_subset <- counts_df %>%
  select(all_of(selected_columns))

# Rename sample columns to anonymized names
for (old_name in selected_sample_columns) {
  if (old_name %in% names(sample_mapping)) {
    new_name <- sample_mapping[old_name]
    counts_subset <- counts_subset %>%
      rename(!!new_name := !!old_name)
  }
}

# Subsample genes to make the dataset smaller (keep ~5000 genes)
# Filter out low-count genes first, then sample
cat("Subsampling genes...\n")

# Calculate row sums and filter genes with reasonable expression
gene_sums <- rowSums(counts_subset[,-1])
high_expr_genes <- which(gene_sums >= quantile(gene_sums, 0.3))  # Keep top 70% by expression

# Sample 5000 genes from the high-expression genes
n_genes_to_keep <- min(5000, length(high_expr_genes))
sampled_gene_indices <- sample(high_expr_genes, n_genes_to_keep)

# Create final count matrix
final_counts <- counts_subset[c(1, sampled_gene_indices + 1), ]  # +1 because Gene column is first
final_counts <- final_counts[order(final_counts$Gene), ]  # Sort by gene name

cat("Final dataset: ", nrow(final_counts), "genes x", ncol(final_counts)-1, "samples\n")

# Save the example count matrix
write_csv(final_counts, "example_data/example_counts.csv")
cat("Saved anonymized count matrix to example_data/example_counts.csv\n")

# Print summary
cat("\n=== SUMMARY ===\n")
cat("Created anonymized example dataset:\n")
cat("- Metadata: example_data/example_metadata.csv (", nrow(final_metadata), "samples )\n")
cat("- Counts: example_data/example_counts.csv (", nrow(final_counts), "genes x", ncol(final_counts)-1, "samples )\n")
cat("\nAnonymization mapping:\n")
cat("- Strains: B6 -> StrainA, BALBc -> StrainB\n")
cat("- Treatments: CTRL -> Control, FA -> TreatmentX\n")
cat("- Tissues: NG -> Tissue1, T-DRG -> Tissue2, L-DRG -> Tissue3\n")
cat("- Sample names: anonymized to Sample_01, Sample_02, etc.\n")

cat("\nExample data creation completed!\n") 