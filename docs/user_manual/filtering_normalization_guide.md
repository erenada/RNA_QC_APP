# Filtering & Normalization Guide

**Author:** Eren Ada, PhD  
**Date:** 6/4/2025  
**Version:** 2.0  

## Overview

The **Filtering & Normalization** tab provides comprehensive tools for preparing RNA-seq data for downstream analysis. This tab performs two critical steps: filtering low-expressed genes to reduce multiple testing burden and noise, and applying normalization methods to account for technical differences between samples. The tab operates on processed data from the Data Input & Validation stage and automatically updates QC plots to reflect the effects of filtering and normalization.

## Tab Organization

The Filtering & Normalization workflow is organized into several main sections:

### 1. Gene Filtering Settings & Results
### 2. Normalization Methods & Evaluation
### 3. Comprehensive Normality Assessment
### 4. Data Preview & Export Options

---

## Data Status Indicator

At the top of the tab, a status indicator shows the current data state:

- **"Raw processed data"** - Using data from validation step
- **"Filtered data"** - After gene filtering is applied
- **"Normalized data"** - After normalization is applied

---

## Gene Filtering Section

### Purpose of Gene Filtering

Gene filtering removes genes with very low expression across samples, which:

- **Improves statistical power** by reducing multiple testing burden
- **Reduces noise** from unreliably detected genes
- **Focuses analysis** on biologically relevant features
- **Accelerates processing** by reducing data size

### Filtering Strategy Guidelines

**Minimum Expression Threshold:** Set a count/CPM value that genes must reach to be considered "expressed"

- **Conservative:** 10+ counts (removes noise while retaining most genes)
- **Moderate:** 5+ counts (retains more lowly expressed genes)
- **Liberal:** 1+ counts (minimal filtering)

**Group-based vs. Global Filtering:**

- **Group-based:** Requires expression in at least one experimental group (recommended)
- **Global:** Requires expression across all samples regardless of condition

### Filtering Controls

**Metadata Grouping Column:**

- **"None":** Apply filtering globally across all samples
- **Select metadata variable:** Apply group-aware filtering (e.g., Treatment, Condition)
- **Recommendation:** Use treatment/condition groups when available

**Minimum Expression Threshold:**

- Numeric value for the expression cutoff
- **Default:** 10 (suitable for most datasets)
- **Range:** 0+ (higher values = more stringent filtering)

**Expression Unit:**

- **Raw Counts:** Use actual count values (affected by library size differences)
- **CPM (Counts Per Million):** Normalized for library size (more consistent across samples)
- **Recommendation:** Use CPM for more consistent filtering

**Minimum Samples Meeting Threshold:**

- Number of samples that must meet the expression threshold
- **If grouping:** Samples within at least one experimental group
- **If no grouping:** Total samples across the entire dataset
- **Recommendation:** Set to smallest group size or 2-3 samples minimum

### Filtering Results

**Filtering Summary:**

- Genes before and after filtering
- Percentage of genes removed
- Group-specific statistics (if applicable)
- Expression distribution changes

**Visualization:**

- **Gene Count Plot:** Bar chart showing genes before/after filtering
- **Expression Distribution Plot:** Density plots comparing expression distributions
- **Filtered Data Preview:** Interactive table of the filtered dataset

### Best Practices for Filtering

**Conservative Approach:**

- Start with default settings (≥10 counts in ≥2 samples)
- Examine filtering impact plots
- Adjust if too many/few genes are removed

**Group-aware Filtering:**

- Use experimental groups (Treatment, Condition) when available
- Ensures genes active in specific conditions are retained
- More biologically relevant than global filtering

**Validation:**

- Check that filtering doesn't remove too many genes (typically retain 60-80%)
- Ensure known important genes are not filtered out
- Review expression distribution changes

---

## Normalization Section

### Purpose of Normalization

Normalization corrects for technical differences between samples:

- **Library size differences** from sequencing depth variation
- **Composition bias** from highly expressed genes
- **Technical artifacts** from sample preparation or sequencing

### Available Normalization Methods

#### Library Size Based Methods

**Total Count (TC):**

- Scales samples to equal total counts
- **Best for:** Datasets with minimal composition bias
- **Limitations:** Sensitive to highly expressed genes

**Median of Ratios (DESeq2):**

- Uses geometric mean of ratios for size factor calculation
- **Best for:** Most RNA-seq datasets (robust default choice)
- **Advantages:** Handles composition bias well, widely validated

**Upper Quartile (UQ):**

- Normalizes using 75th percentile of non-zero counts
- **Best for:** Datasets with many zeros or low counts
- **Advantages:** Less sensitive to highly expressed genes than total count

**Relative Log Expression (RLE):**

- Similar to DESeq2 but uses median ratios
- **Best for:** Datasets with stable reference genes
- **Advantages:** Robust to outlier genes

#### Distribution Transformation Methods

**CPM (Counts Per Million):**

- Simple library size normalization
- **Best for:** Initial exploration, correlation analysis
- **Formula:** (counts / total counts) × 1,000,000

**TMM (Trimmed Mean of M-values):**

- edgeR's normalization method using trimmed means
- **Best for:** Differential expression analysis with edgeR
- **Advantages:** Robust to extreme values

**VST (Variance Stabilizing Transformation):**

- DESeq2's variance stabilization for count data
- **Best for:** Exploratory analysis, PCA, heatmaps
- **Advantages:** Removes mean-variance dependence

**rlog (Regularized Log Transformation):**

- DESeq2's log transformation with shrinkage
- **Best for:** Small datasets (<30 samples), visualization
- **Advantages:** Better handling of low counts than log2

**Quantile Normalization:**

- Forces all samples to have identical distributions
- **Best for:** Microarray-style normalization (use with caution)
- **Limitations:** Strong distributional assumptions

### Normalization Evaluation

The tab provides comprehensive evaluation of normalization effectiveness through five sub-tabs:

#### Library Size Evaluation

**Purpose:** Assess whether normalization successfully reduces library size variation

**Metrics:**

- **Coefficient of Variation (CV):** Lower CV after normalization indicates success
- **Mean and median library sizes:** Before vs. after comparison
- **Visualization:** Density plots of log2 library sizes

**Interpretation:**

- **Good normalization:** CV reduction of 50%+ 
- **Excellent normalization:** CV reduction of 70%+
- **Poor normalization:** Little to no CV reduction

#### Normalization Effects

**Purpose:** Evaluate the impact on mean-variance relationships

**Metrics:**

- **Mean-variance correlation:** Should decrease after normalization
- **Scatter plots:** Show mean vs. variance before/after normalization
- **LOESS smoothing:** Highlights trends in mean-variance relationships

**Interpretation:**

- **Good normalization:** Reduced mean-variance correlation
- **VST/rlog:** Should show flat variance across expression levels
- **Count-based methods:** Some mean-variance relationship may remain

#### Batch Effect Check

**Purpose:** Assess whether normalization reduces technical batch effects

**Controls:**

- **Batch Variable Selection:** Choose metadata column representing technical batches
- **PCA Visualization:** Samples colored by batch variable

**Metrics:**

- **R² from batch effects:** Lower values indicate reduced batch effects
- **Visual separation:** Less clustering by batch in PCA

**Interpretation:**

- **Good normalization:** Biological groups separate more than technical batches
- **Poor normalization:** Strong clustering by technical factors
- **May require:** Additional batch correction methods

#### Normality Check

**Purpose:** Quick assessment of post-normalization data normality

**Metrics:**

- **Summary statistics:** Skewness and kurtosis across samples
- **Q-Q plots:** Visual assessment of normality
- **Shapiro-Wilk tests:** Statistical normality testing

**Applications:**

- **Parametric methods:** Require approximately normal data
- **Correlation choice:** Influences Pearson vs. Spearman selection

#### Data Distribution Analysis

**Purpose:** Comprehensive normality assessment with multiple visualization and testing options

**Controls:**

- **Sample Selection:** Choose specific samples for detailed analysis
- **Test Methods:** Shapiro-Wilk, Kolmogorov-Smirnov, Anderson-Darling
- **Significance Levels:** 0.05, 0.01, 0.001

**Visualizations:**

- **Q-Q Plots:** Compare sample distributions to normal distribution
- **Density Plots:** Actual vs. theoretical normal distributions
- **Histograms:** Combined histogram and density overlays
- **Detailed Results:** Comprehensive statistical test results

**Automated Recommendations:**

The system provides guidance for:

- **Correlation method choice:** Pearson vs. Spearman based on normality
- **Downstream analysis approaches:** Parametric vs. non-parametric methods
- **Additional transformation needs:** If normality is poor

---

## Normalization Method Selection Guide

### For Most RNA-seq Datasets

**Recommended:** DESeq2 (Median of Ratios)

- Well-validated and robust
- Handles composition bias effectively
- Suitable for differential expression analysis

### For Exploratory Analysis

**Recommended:** VST or rlog (for visualization), CPM (for correlation)

- VST: Better for larger datasets (>30 samples)
- rlog: Better for smaller datasets (<30 samples)
- CPM: Simple and interpretable

### For Specific Pipelines

**edgeR workflow:** TMM normalization

**Correlation analysis:** CPM with log2 transformation

**Machine learning:** Quantile normalization (with caution)

### Special Considerations

**Small datasets (<10 samples):** Consider rlog or conservative approaches

**Large datasets (>100 samples):** VST typically preferred over rlog

**Compositionally biased data:** Avoid Total Count, prefer DESeq2 or TMM

**Batch effects present:** May need additional batch correction after normalization

---

## Quality Assessment Guidelines

### Excellent Normalization Indicators

- Library size CV reduction >70%
- Mean-variance correlation significantly reduced
- Biological groups separate clearly in PCA
- Batch effects minimized
- Improved data normality (if using VST/rlog)

### Good Normalization Indicators

- Library size CV reduction 50-70%
- Some improvement in mean-variance relationship
- Biological signal preserved or enhanced
- Batch effects reduced but may still be present

### Concerning Normalization Indicators

- Little to no CV reduction (<30%)
- Mean-variance relationship unchanged or worsened
- Loss of biological signal
- Increased batch effects
- Severe distributional artifacts

---

## Troubleshooting Common Issues

### Poor Normalization Performance

**Possible Causes:**

- Inappropriate method for data type
- Severe composition bias
- Strong batch effects
- Highly heterogeneous sample groups

**Solutions:**

- Try alternative normalization methods
- Consider more stringent pre-filtering
- Investigate batch correction methods
- Examine sample preparation protocols

### High Library Size Variation After Normalization

**Possible Causes:**

- Method not designed for library size equalization
- Severe technical artifacts
- Inappropriate data transformation

**Solutions:**

- Use library size-focused methods (TC, DESeq2, UQ)
- Investigate underlying technical issues
- Consider additional QC steps

### Loss of Biological Signal

**Possible Causes:**

- Over-normalization
- Inappropriate method for experimental design
- Batch effects confounded with biology

**Solutions:**

- Try less aggressive normalization
- Use method appropriate for experimental design
- Careful batch effect analysis

---

## Integration with Workflow

### From QC Plots & Summaries

- Normality assessment guides normalization method choice
- Sample correlation patterns inform filtering stringency
- PCA results help evaluate batch effects

### To Downstream Analysis

- Normalized data automatically updates QC plots
- Filtered gene set reduces multiple testing burden
- Normalization method choice affects downstream statistical approaches

### Data Export Options

- **Current Data Download:** Filtered and/or normalized data in CSV format
- **Analysis Report:** Comprehensive summary of filtering and normalization steps
- **Interactive Preview:** Real-time view of current data state

---

## Best Practices Summary

### Before Starting

1. **Review QC results** to understand data characteristics
2. **Identify experimental groups** for group-aware filtering
3. **Check for batch effects** that may influence normalization choice

### During Filtering

1. **Start conservatively** with default parameters
2. **Use group-aware filtering** when experimental groups are available
3. **Monitor filtering impact** using provided visualizations
4. **Ensure adequate gene retention** (typically 60-80% of genes)

### During Normalization

1. **Choose appropriate method** based on data characteristics and intended use
2. **Evaluate normalization effectiveness** using all provided metrics
3. **Check for preserved biological signal** in PCA plots
4. **Assess batch effect reduction** if batch variables are available

### After Processing

1. **Review updated QC plots** to confirm improvements
2. **Document method choices** and rationale
3. **Export processed data** with appropriate naming
4. **Proceed to downstream analysis** with confidence

---

## Technical Notes

### Data Processing Pipeline

- **Input:** Processed counts from Data Input & Validation
- **Filtering:** Applied to raw counts, then normalization applied to filtered data
- **State tracking:** System maintains data provenance throughout processing

### Computational Considerations

- **Memory usage:** Scales with dataset size; large datasets may require more time
- **Method complexity:** VST/rlog computationally intensive for large datasets
- **Parallel processing:** Some methods utilize multiple cores when available

### Statistical Methods

- **Size factor calculation:** Uses geometric means and ratios (DESeq2, RLE)
- **Normality testing:** Multiple statistical tests with appropriate corrections
- **Variance stabilization:** Empirical Bayes approaches (VST, rlog)

---

## Support and Resources

For additional help:

- Consult the **[QC Plots & Summaries Interpretation Guide](qc_plots_interpretation_guide.html)** for understanding evaluation metrics
- Review **[Data Input & Validation Guide](data_input.html)** for upstream data quality issues
- Check the **DESeq2** and **edgeR** package documentation for method-specific details
- Refer to the troubleshooting section for common problems and solutions

*This guide provides comprehensive coverage of filtering and normalization options. Method selection should be based on your specific experimental design, data characteristics, and downstream analysis plans. When in doubt, DESeq2 normalization provides a robust starting point for most RNA-seq datasets.* 