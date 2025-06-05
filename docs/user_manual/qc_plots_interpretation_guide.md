# QC Plots & Summaries Interpretation Guide

**Author:** Eren Ada, PhD  
**Date:** 6/4/2025  
**Version:** 2.0  

## Overview

The **QC Plots & Summaries** tab provides comprehensive quality control visualizations and statistical analyses to assess RNA-seq data quality. This tab operates on your processed data from the Data Input & Validation stage and automatically updates when normalization is applied. The analyses help identify potential issues, batch effects, outliers, and guide downstream analysis decisions.

## Tab Organization

The QC analysis is organized into four main sub-tabs, each focusing on different aspects of data quality:

### 1. Basic QC Metrics
### 2. Sample Similarity Analysis  
### 3. Sample Correlation Analysis
### 4. Normality Assessment (available in Filtering & Normalization tab)

---

## Data Status Indicator

At the top of each tab, you'll see a status indicator showing which data is currently being analyzed:

- **"Displaying plots for: log2(CPM+1) transformed data (pre-normalization)"** - Initial analysis on processed raw data
- **"Displaying plots for: [Method] normalized data"** - Analysis after normalization is applied

---

## Basic QC Metrics Tab

### Library Size Distribution

**Purpose:** Assesses the sequencing depth consistency across samples and compares raw vs. processed data.

**Controls:**

- **Plot Type:** Choose between boxplot or barplot visualization
- **Log Scale Y-axis:** Enable logarithmic scaling for better visualization of large dynamic ranges

**Interpretation:**

- **Consistent libraries:** Similar library sizes across samples indicate good sequencing consistency
- **Large variations:** May indicate technical issues or require normalization
- **Raw vs. Processed comparison:** Shows the effect of data processing steps

**Typical Issues:**

- Libraries differing by >10-fold may need normalization
- Extremely small libraries (<100K reads) may be low-quality samples

### Gene Detection Rates

**Purpose:** Evaluates how many genes are detectable in each sample at a given threshold.

**Controls:**

- **Detection Threshold:** Minimum count value to consider a gene "detected" (default: 1)

**Interpretation:**

- **High detection rates (>60%):** Good sample quality with broad gene coverage
- **Low detection rates (<40%):** May indicate degraded RNA, low sequencing depth, or contamination
- **Consistent rates across samples:** Good experimental consistency

**Typical Ranges:**

- **Excellent:** 70-85% detection rate
- **Good:** 60-70% detection rate  
- **Concerning:** <50% detection rate

### Summary Statistics Table

**Content:** Comprehensive statistics comparing raw and processed data:

- Library size metrics (min, max, mean, median)
- Gene detection statistics across all samples
- Allows quick identification of outlier samples

---

## Sample Similarity Analysis Tab

### Principal Component Analysis (PCA)

**Purpose:** Reduces dimensionality to visualize sample relationships and identify patterns, outliers, or batch effects.

**Controls:**

- **X-axis PC / Y-axis PC:** Choose which principal components to display (PC1-PC10)
- **Color by:** Color points by metadata variables or sample names
- **Show Sample Labels:** Toggle sample name display on 2D plots

**Features:**

- **2D PCA Plot:** Standard scatter plot for easy interpretation
- **3D PCA Plot:** Interactive 3D visualization for exploring additional dimensions
- **Variance Explained:** Percentage of variance captured by each PC

**Interpretation:**

- **PC1 and PC2:** Usually capture the most biological variation
- **Tight clustering:** Samples from same condition should cluster together
- **Outliers:** Points far from others may indicate poor quality or batch effects
- **Separation patterns:** Should align with experimental design, not technical factors

**Typical Results:**

- **Good:** Clear separation between biological conditions, tight within-group clustering
- **Concerning:** No clear pattern, excessive scatter, or clustering by technical factors

**Download Options:**

- PCA plot (high-resolution image)
- PCA statistics and variance explained data

---

## Sample Correlation Analysis Tab

### Correlation Heatmap

**Purpose:** Visualizes pairwise correlations between all samples to identify outliers and assess overall data quality.

**Correlation Settings:**

- **Data Source:**
  - **log2(CPM+1) transformed data:** Default option, recommended for most analyses. Uses counts per million normalization with log2 transformation.
  - **Raw count matrix:** Uses untransformed count data. Requires careful consideration of correlation method choice.
- **Correlation Method:**
  - **Pearson:** Use for normally distributed, log-transformed data (recommended with log2(CPM+1) data)
  - **Spearman:** Use for non-normal data, raw counts, or when outliers are present (recommended with raw count matrix)
- **Color Scheme:** Multiple palettes available for visualization preferences

**Important Warning:** When using raw count matrix, ensure you select the appropriate correlation method (typically Spearman) and consider applying normalization in the Filtering & Normalization tab first. Raw count data often violates assumptions of Pearson correlation due to non-normal distributions and heteroscedasticity.

**Data Source Recommendations:**

- **Use log2(CPM+1) transformed data when:**
  - Performing standard correlation analysis
  - Data appears reasonably normal after transformation
  - You want to use Pearson correlation
- **Use raw count matrix when:**
  - Specifically required for your analysis pipeline
  - Comparing with results from other tools that use raw counts
  - Using non-parametric methods throughout your analysis

**Advanced Filtering:**

- **Sample Filtering:** Filter samples based on metadata variables
- **Dynamic Filters:** Create custom sample subsets for focused analysis

**Display Options:**

- **Show Correlation Values:** Display numerical correlation coefficients
- **Show Significance Stars:** Mark statistically significant correlations
- **Auto-scale Colors:** Optimize color range for high correlation datasets

**Interpretation:**

- **High correlations (>0.85):** Expected for good quality samples
- **Low correlations (<0.7):** May indicate poor quality or batch effects
- **Block patterns:** Should reflect experimental design
- **Outlier samples:** Consistently low correlations with all other samples

**Correlation Summary Statistics:**

- Data source used for correlation analysis (log2(CPM+1) transformed data, raw count matrix, or normalized data)
- Mean and median correlation values
- Range of correlations observed
- Number and percentage of significant correlations
- Sample filtering information

**Download Options:**

- High-resolution heatmap (PDF format)
- Correlation matrix data

---

## Normality Assessment

**Note:** Full normality assessment is available in the Filtering & Normalization tab, with comprehensive tools for:

### Distribution Analysis

- **Q-Q Plots:** Compare sample distributions to normal distribution
- **Density Plots:** Visualize actual vs. theoretical normal distributions  
- **Histogram Overlays:** Combined histogram and density analysis

### Statistical Testing

- **Multiple Test Methods:**
  - Shapiro-Wilk test (default, most sensitive)
  - Kolmogorov-Smirnov test
  - Anderson-Darling test
- **Customizable significance levels:** 0.05, 0.01, 0.001

### Normality Recommendations

The system provides automated recommendations for:

- **Correlation method choice:** Pearson vs. Spearman based on normality
- **Downstream analysis guidance:** Transformation recommendations
- **Sample-specific assessments:** Individual sample quality evaluation

---

## Data Quality Guidelines

### Excellent Quality Indicators

- Library sizes within 2-fold range
- Gene detection rates >70%
- Sample correlations >0.85
- Clear PCA clustering by biological groups
- Normal or near-normal data distributions

### Good Quality Indicators  

- Library sizes within 5-fold range
- Gene detection rates 60-70%
- Sample correlations >0.75
- Identifiable patterns in PCA
- Moderate deviations from normality

### Concerning Quality Indicators

- Library sizes >10-fold differences
- Gene detection rates <50%
- Sample correlations <0.7
- No clear PCA patterns or excessive scatter
- Severe non-normality

---

## Troubleshooting Common Issues

### Low Sample Correlations

**Possible Causes:**

- Batch effects
- Sample degradation
- Contamination
- Incorrect sample labeling

**Solutions:**

- Check metadata for batch variables
- Consider batch correction methods
- Review sample preparation protocols
- Verify sample identity

### Poor PCA Clustering

**Possible Causes:**

- Weak biological signal
- Batch effects dominating
- Heterogeneous sample groups
- Technical noise

**Solutions:**

- Filter low-expression genes more stringently
- Apply normalization methods
- Check for batch effects in metadata
- Consider additional QC filtering

### Non-Normal Distributions

**Implications:**

- May affect downstream statistical tests
- Could indicate data transformation needs
- Might suggest correlation method choice

**Actions:**

- Use Spearman correlation instead of Pearson
- Consider alternative transformations (VST, rlog)
- Check for outliers or contamination

---

## Best Practices

### Before Normalization

1. **Review all basic metrics** to identify obvious quality issues
2. **Examine PCA plots** for biological vs. technical patterns  
3. **Check sample correlations** to identify outliers
4. **Assess normality** to guide method choices

### After Normalization

1. **Compare before/after plots** to evaluate normalization effectiveness
2. **Re-examine correlations** to ensure improvement
3. **Check if biological patterns** become clearer in PCA
4. **Verify normality improvements** if using parametric methods

### Decision Making

- **Remove samples** with consistently poor metrics across multiple measures
- **Apply normalization** if library size differences are large
- **Use robust methods** (Spearman, rank-based) for non-normal data
- **Document decisions** and rationale for reproducibility

---

## Integration with Other Tabs

### From Data Input & Validation

- Processed count matrix and metadata are automatically loaded
- Validation issues should be resolved before QC analysis
- Library size plots provide continuity between tabs

### To Filtering & Normalization  

- QC results inform filtering and normalization parameter choices
- Normality assessment guides transformation method selection
- Poor quality samples may be excluded before normalization

### Throughout the Workflow

- QC plots automatically update when normalization is applied
- Comparative analysis shows the effect of processing steps
- Quality metrics guide downstream analysis decisions

---

## Technical Notes

### Data Transformations

- **Initial analysis:** Uses log2(CPM+1) transformation of processed counts
- **Post-normalization:** Uses the specific normalization method applied
- **PCA calculation:** Performed on transposed data matrix (samples as observations)

### Statistical Methods

- **Correlation calculations:** Support both Pearson and Spearman methods
- **Significance testing:** Multiple testing correction applied where appropriate
- **PCA:** Uses prcomp() with centering and optional scaling

### Performance Considerations

- **Large datasets:** Some plots may take time to render with many samples
- **Memory usage:** PCA and correlation calculations scale with dataset size
- **Interactive features:** 3D PCA and plotly graphics require modern web browsers

---

## Support and Resources

For additional help:

- Consult the **Filtering & Normalization Guide** for normalization-specific QC
- Review **Data Input & Validation Guide** for upstream quality issues
- Check the application logs for detailed processing information
- Refer to the methodology documentation for statistical details

*This guide provides comprehensive coverage of QC interpretation. For specific research contexts, consider consulting with a bioinformatics specialist to ensure appropriate quality thresholds and analysis strategies.* 