# Quick Start Guide

**Author:** Eren Ada, PhD  
**Date:** 6/5/2025  
**Version:** 2.0.0

**Estimated time:** 15-20 minutes

This guide walks you through a complete analysis using example data to familiarize you with the tool's workflow.

## Before You Begin

### System Requirements
- R version 4.0 or higher (4.3.0+ recommended)
- At least 4GB RAM (8GB recommended for larger datasets)
- Modern web browser (Chrome, Firefox, Safari, or Edge)

### Data Format Requirements
- **Count Matrix**: CSV file with genes as rows, samples as columns
- **Metadata**: CSV file with samples as rows, experimental factors as columns
- Sample names must match between count matrix and metadata exactly

## Step 1: Launch the Application

1. Open R/RStudio
2. Navigate to the APP directory
3. Run: `shiny::runApp("app.R")`
4. The application will open in your default web browser

## Step 2: Data Input & Validation

### Upload Your Files
1. Navigate to the **"Data Input & Validation"** tab
2. **Upload Count Matrix**: Click button and select your CSV file
3. **Upload Sample Metadata**: Click button and select your metadata file

### Automatic File Processing
The application will automatically:
- **Detect file formats** (CSV, TSV) and separators
- **Validate data types** (ensure non-negative integer counts)
- **Round decimal values** to integers with notification
- **Check sample consistency** between count matrix and metadata

### Review Validation Results
After upload, examine:
- **✅ Success messages** for completed validation steps
- **⚠️ Warnings** for automatically resolved issues
- **❌ Errors** requiring user intervention
- **📊 Data summary** with dataset characteristics

### Handle Duplicate Genes (if detected)
If duplicate gene IDs are found, choose from four strategies:
1. **Sum Counts for Duplicates** (Recommended) - Combines expression values
2. **Keep Highest Expression** - Retains duplicate with highest average
3. **Rename with Suffix** - Makes unique names (GENE1_1, GENE1_2)
4. **Remove All Duplicates** - Eliminates all problematic genes

**Proceed only when validation shows "Data successfully processed"**

### Available Downloads
- **Original Count Matrix**: Unprocessed backup data
- **Processed Count Matrix**: After duplicate handling and validation
- **Duplicate Statistics**: Detailed duplicate gene analysis (if applicable)

## Step 3: Quality Control Analysis

Navigate to **"QC Plots & Summaries"** tab for comprehensive data quality assessment.

### Basic QC Metrics
**Library Size Distribution:**
- **Bar plots/Box plots**: Total counts per sample
- **Log scale option**: Better visualization of large ranges
- **Interpretation**: Libraries within 2-5x range are good; >10x differences concerning

**Gene Detection Rates:**
- **Detection threshold**: Adjustable minimum count (default: 1)
- **Good quality**: >60% detection rate
- **Excellent quality**: 70-85% detection rate
- **Concerning**: <50% detection rate

### Sample Similarity Analysis (PCA)
**Principal Component Analysis:**
- **2D PCA Plot**: Standard scatter plot for interpretation
- **3D PCA Plot**: Interactive visualization for additional dimensions
- **Axis Selection**: Choose PC combinations (PC1-PC10)
- **Color Grouping**: Color by metadata variables or sample names
- **Sample Labels**: Toggle display on/off

**Expected Results:**
- **Good**: Clear separation between biological conditions
- **Concerning**: No clear patterns or clustering by technical factors

**Available Downloads:**
- **PCA Plot**: High-resolution PDF format
- **PCA Statistics**: Variance explained and loadings data (CSV)

### Sample Correlation Analysis
**Correlation Heatmap:**
- **Data Source Options**:
  - log2(CPM+1) transformed data (recommended)
  - Raw count matrix (use with Spearman correlation)
- **Correlation Methods**:
  - Pearson (for log-transformed data)
  - Spearman (for raw counts or non-normal data)
- **Multiple color schemes** available

**Quality Indicators:**
- **High correlations (>0.85)**: Good quality samples
- **Low correlations (<0.7)**: May indicate quality issues or batch effects
- **Block patterns**: Should reflect experimental design

**Available Downloads:**
- **Correlation Heatmap**: High-resolution PDF format
- **Correlation Matrix**: Statistical data (CSV)

**Click "Continue to Filtering & Normalization" when satisfied with QC**

## Step 4: Filtering & Normalization

### Gene Filtering
**Purpose**: Remove low-expressed genes to reduce noise and multiple testing burden

**Filtering Parameters:**
1. **Grouping Column**: Select metadata variable for group-aware filtering
2. **Expression Threshold**: Minimum count or CPM value (default: 10)
3. **Expression Unit**: Raw counts vs. CPM (recommended: CPM)
4. **Minimum Samples**: Number of samples that must meet threshold

**Filtering Strategies:**
- **Conservative**: ≥10 counts in ≥2 samples (retains 60-80% of genes)
- **Moderate**: ≥5 counts in group-aware setting
- **Liberal**: ≥1 count for minimal filtering

**Review Results:**
- **Gene count plots**: Before/after filtering comparison
- **Expression distributions**: Impact visualization
- **Filtered data preview**: Interactive table

### Normalization
**Purpose**: Correct for technical differences (library size, composition bias)

**Available Methods:**

**Library Size-Based:**
- **DESeq2 (Median of Ratios)**: Recommended for most datasets
- **Total Count (TC)**: Simple library size scaling
- **Upper Quartile (UQ)**: Robust to highly expressed genes
- **Relative Log Expression (RLE)**: Median ratio method

**Distribution Transformations:**
- **CPM**: Counts per million (simple, for exploration)
- **TMM**: Trimmed mean of M-values (for edgeR workflows)
- **VST**: Variance stabilizing transformation (for visualization, >30 samples)
- **rlog**: Regularized log transformation (for small datasets, <30 samples)

**Normalization Evaluation:**
The application provides comprehensive assessment through:
- **Library Size Evaluation**: CV reduction analysis
- **Mean-Variance Relationship**: Before/after scatter plots
- **Batch Effect Assessment**: PCA colored by batch variables
- **Distribution Analysis**: Q-Q plots, density plots, histograms
- **Statistical Testing**: Multiple normality tests with recommendations

**Quality Indicators:**
- **Excellent**: CV reduction >70%, clear biological separation in PCA
- **Good**: CV reduction 50-70%, improved mean-variance relationship
- **Concerning**: <30% CV reduction, loss of biological signal

## Step 5: Export Results

### Data Downloads
- **Current Data**: Filtered/normalized counts with timestamps
  - Format: CSV with gene IDs as first column
  - Filename: `filtered_normalized_data_YYYYMMDD_HHMMSS.csv`

### Report Generation
- **Analysis Report**: Comprehensive HTML report including:
  - Dataset overview and processing steps
  - Filtering and normalization parameters
  - Session information for reproducibility
  - Filename: `analysis_report_YYYYMMDD_HHMMSS.html`

### Visualization Downloads
All plots available as publication-quality PDFs with timestamps:
- **PCA plots**: Vector graphics, scalable
- **Correlation heatmaps**: High-resolution format
- **Statistical data**: CSV format for further analysis

## Example Workflow Timeline

| Step | Time | Key Actions |
|------|------|-------------|
| Data Upload | 2-3 min | Upload files, review validation |
| Duplicate Handling | 1-2 min | Choose strategy if needed |
| QC Analysis | 5-8 min | Examine all QC metrics and plots |
| Filtering | 3-5 min | Set parameters, review impact |
| Normalization | 3-5 min | Select method, evaluate results |
| Export | 1-2 min | Download processed data and reports |

## Common First-Time Issues & Solutions

### Data Upload Problems
- **File format**: Ensure CSV/TSV with standard delimiters
- **Sample names**: Must match exactly (case-sensitive, no extra spaces)
- **Missing values**: Replace NA/NULL with 0 or remove before upload
- **File size**: Large files may take time; ensure stable connection

### QC Interpretation Guidance
- **Outlier samples**: Low correlation with all others, extreme library sizes
- **Batch effects**: Unwanted clustering in PCA by technical factors
- **Poor clustering**: May indicate weak biological signal or technical noise
- **Low detection rates**: Could indicate RNA degradation or low sequencing depth

### Filtering & Normalization Decisions
- **Too stringent filtering**: Losing >50% of genes, may remove biological signal
- **Too lenient filtering**: Keeping mostly zero-count genes, increases noise
- **Group-aware vs. global**: Use experimental groups when available
- **Normalization method**: DESeq2 for most cases, VST/rlog for visualization

### Performance Optimization
- **Large datasets**: Close other applications, allow extra processing time
- **Browser issues**: Use Chrome/Firefox, disable extensions if needed
- **Memory concerns**: Filter aggressively for very large datasets

## Quality Assessment Guidelines

### Excellent Quality Indicators
- Library sizes within 2-fold range
- Gene detection rates >70%
- Sample correlations >0.85
- Clear PCA clustering by biological groups
- Successful normalization (CV reduction >70%)

### Concerning Quality Indicators
- Library sizes >10-fold differences
- Gene detection rates <50%
- Sample correlations <0.7
- No clear PCA patterns or technical clustering
- Poor normalization effectiveness

## Next Steps After Quick Start

### For Comprehensive Understanding
1. **Detailed Guides**: Review specific documentation for each step:
   - [Data Input & Validation Guide](data_input.html) - Advanced options and troubleshooting
   - [QC Plots Interpretation Guide](qc_plots_interpretation_guide.html) - Detailed quality assessment
   - [Filtering & Normalization Guide](filtering_normalization_guide.html) - Method selection and evaluation
   - [Export & Download Guide](export_download_guide.html) - Data management and sharing

2. **Advanced Features**: Explore additional functionality:
   - Custom sample filtering in correlation analysis
   - Detailed normality assessment tools
   - Multiple PCA plot combinations
   - Comprehensive batch effect evaluation

### For Downstream Analysis
- **DESeq2 workflows**: Use filtered, non-normalized counts
- **edgeR workflows**: Use TMM-normalized data
- **Visualization/clustering**: Use VST or rlog transformed data
- **Custom analysis**: Export session data for reproducibility

## Getting Help

### Built-in Resources
- **Validation messages**: Address all warnings and errors
- **Interactive help**: Hover over UI elements for tooltips
- **Example data**: Test with provided datasets in `example_data/`

### Documentation
- **[Troubleshooting Guide](troubleshooting.html)**: Common issues and solutions
- **[Technical Requirements](technical_requirements.html)**: System setup and optimization
- **[User Manual Overview](README.html)**: Feature overview and navigation guide

### Support
- **Contact**: Eren Ada, PhD - erenada@gmail.com
- **GitHub**: https://github.com/erenada
- **Issue Reporting**: Include system info, error messages, and dataset characteristics

---

**Congratulations!** You've completed the quick start workflow. Your data is now quality-controlled, filtered, normalized, and ready for downstream differential expression analysis or other RNA-seq workflows. 