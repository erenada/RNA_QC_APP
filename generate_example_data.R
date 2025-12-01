# Generate Example RNA-seq Count Data
# Creates realistic count data for demonstration purposes

set.seed(12345)  # For reproducibility

# Read metadata
metadata <- read.csv("example_data/example_metadata.csv", stringsAsFactors = FALSE)
samples <- metadata$Sample

# Parameters
n_genes <- 5000
n_samples <- length(samples)

# Generate gene names
gene_names <- paste0("GENE_", sprintf("%05d", 1:n_genes))

# Create realistic count data
# Different expression patterns based on experimental design
counts_matrix <- matrix(0, nrow = n_genes, ncol = n_samples)
rownames(counts_matrix) <- gene_names
colnames(counts_matrix) <- samples

# Generate realistic RNA-seq counts
for (i in 1:n_genes) {
  # Base expression level
  base_expr <- rpois(1, lambda = 100)
  
  for (j in 1:n_samples) {
    strain <- metadata$Strain[j]
    treatment <- metadata$Treatment[j]
    tissue <- metadata$Tissue[j]
    
    # Modify expression based on conditions
    lambda <- base_expr
    
    # Strain effect (some genes)
    if (i <= 1000 && strain == "StrainB") {
      lambda <- lambda * runif(1, 0.5, 2.0)
    }
    
    # Treatment effect (some genes)
    if (i <= 1500 && i > 500 && treatment == "TreatmentX") {
      lambda <- lambda * runif(1, 0.3, 3.0)
    }
    
    # Tissue effect (some genes)
    if (i <= 2000 && i > 1000) {
      if (tissue == "Tissue2") lambda <- lambda * runif(1, 0.4, 2.5)
      if (tissue == "Tissue3") lambda <- lambda * runif(1, 0.6, 1.8)
    }
    
    # Add biological variation
    lambda <- lambda * runif(1, 0.8, 1.2)
    
    # Generate count with overdispersion
    counts_matrix[i, j] <- rnbinom(1, size = 10, mu = max(1, lambda))
  }
}

# Add some very lowly expressed genes
low_expr_genes <- sample(1:n_genes, 1000)
for (i in low_expr_genes) {
  counts_matrix[i, ] <- rpois(n_samples, lambda = runif(1, 0, 5))
}

# Convert to data frame and add gene column
counts_df <- data.frame(Gene = gene_names, counts_matrix, stringsAsFactors = FALSE)

# Write to file
write.csv(counts_df, "example_data/example_counts.csv", row.names = FALSE)

cat("Example data generated successfully!\n")
cat("Files created:\n")
cat("- example_data/example_metadata.csv (24 samples)\n")
cat("- example_data/example_counts.csv (5,000 genes x 24 samples)\n") 