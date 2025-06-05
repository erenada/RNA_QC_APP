# Export & Download Guide

**Author:** Eren Ada, PhD  
**Date:** 6/5/2025  
**Version:** 2.0.0

## Overview

The Export & Download functionality is integrated throughout the RNA-seq QC & Pre-processing Tool, allowing users to download processed data, visualizations, and reports at multiple stages of the analysis workflow. This guide covers all available export options and best practices for data preservation and reproducibility.

## Table of Contents

1. [Export Locations](#export-locations)
2. [Data Export Options](#data-export-options)
3. [Visualization Downloads](#visualization-downloads)
4. [Report Generation](#report-generation)
5. [Session Data](#session-data)
6. [File Formats & Compatibility](#file-formats--compatibility)
7. [Best Practices](#best-practices)

## Export Locations

Export functionality is available in multiple tabs:

### Data Input & Validation Tab
- Original and processed count matrices
- Duplicate gene statistics
- Validation reports

### QC Plots & Summaries Tab
- All quality control visualizations
- PCA data and statistics
- Correlation matrices
- Summary statistics

### Filtering & Normalization Tab
- Filtered count matrices
- Normalized data in multiple formats
- Before/after comparison data
- Normalization diagnostic plots

## Data Export Options

### Count Matrices

#### **Original Count Matrix**
- **Content**: Unprocessed data as originally uploaded
- **Format**: CSV with gene IDs as first column
- **Use Case**: Backup, alternative analysis pipelines
- **Location**: Data Input & Validation tab

#### **Processed Count Matrix**
- **Content**: Data after duplicate handling and validation
- **Format**: CSV ready for downstream analysis
- **Use Case**: Input for other RNA-seq tools
- **Location**: Data Input & Validation tab

#### **Filtered Count Matrix**
- **Content**: Data after gene filtering applied
- **Format**: CSV with filtered gene set
- **Use Case**: Reduced dataset for focused analysis
- **Location**: Filtering & Normalization tab

#### **Filtered/Normalized Count Matrix**
- **Content**: Filtered and normalized expression data (current data state)
- **Format**: CSV format with timestamps
- **Filename**: `filtered_normalized_data_YYYYMMDD_HHMMSS.csv`
- **Use Case**: Ready for differential expression analysis
- **Location**: Filtering & Normalization tab

### Metadata Export

#### **Original Metadata**
- **Content**: Sample information as uploaded
- **Format**: CSV format
- **Use Case**: Documentation, alternative analyses

#### **Processed Metadata**
- **Content**: Validated and potentially modified metadata
- **Format**: CSV with consistent formatting
- **Use Case**: Downstream analysis with processed data

## Visualization Downloads

### Individual Plot Downloads

All visualizations can be downloaded in multiple formats:

#### **Supported Formats**
- **PDF**: Vector graphics for publication (primary format)
- **CSV**: For statistical data (PCA statistics)
- **Standard Naming**: Downloads include timestamps for version control

#### **Available Downloads**

**From QC Plots & Summaries Tab:**
- **PCA Plot**: 2D PCA visualization as PDF
- **Correlation Heatmap**: Sample correlation heatmap as PDF  
- **PCA Statistics**: Variance explained and loadings as CSV

**Automatic Features:**
- All plots generated at publication quality
- Timestamps included in filenames
- Vector format for scalability

### Plot Customization Options

#### **PCA Plot Customization**
- **Axis Selection**: Choose any PC combination (PC1-PC10)
- **Color Grouping**: Color by metadata variables or sample names
- **Label Display**: Toggle sample name labels on/off
- **Automated Settings**: Publication-ready dimensions (8x6 inches)

#### **Correlation Heatmap Options**
- **Data Source**: Raw counts vs. log2(CPM+1) transformed data
- **Correlation Method**: Pearson vs. Spearman correlation
- **Color Schemes**: Multiple visualization palettes available
- **Sample Filtering**: Filter by metadata for focused analysis

## Report Generation

### Analysis Report

**Content Includes:**
- Dataset overview and statistics
- Current data state and processing steps
- Filtering parameters (if applied)
- Normalization method and parameters (if applied)
- Session information for reproducibility

**Format:**
- **HTML**: Interactive report with R session information
- **Filename**: `analysis_report_YYYYMMDD_HHMMSS.html`
- **Location**: Filtering & Normalization tab

### Analysis Parameter Log

**Content Includes:**
- All parameter settings used
- Processing steps applied
- Timestamp information
- Software version details

**Use Case:** Reproducibility and method documentation

## Session Data & Reproducibility

### Analysis Report Components

**HTML Report includes:**
- Complete session information (sessionInfo())
- R version and package versions
- Analysis timestamp
- Parameter settings used for filtering and normalization
- Data processing steps applied

### Manual Session Saving

For advanced users who need complete reproducibility:
```r
# Save your current workspace manually
save.image("analysis_session.RData")

# Or save specific objects
saveRDS(processed_data, "processed_counts.rds")
saveRDS(metadata, "sample_metadata.rds")
```

## File Formats & Compatibility

### Count Data Format

#### **CSV Format (Standard)**
All data exports use CSV format with the following structure:
```
Gene_ID,Sample1,Sample2,Sample3
GENE1,123,456,789
GENE2,234,567,890
```

**Key Features:**
- **Universal Compatibility**: Readable by all analysis software
- **Gene IDs**: Always included as first column
- **Timestamps**: Filenames include date/time for version tracking
- **No Row Names**: Clean format for downstream import

### Downstream Tool Compatibility

#### **DESeq2 Compatible**
- Count matrices with integer values
- Metadata with sample information
- Recommended: Use filtered, non-normalized counts

#### **edgeR Compatible**
- DGEList objects (via .rds export)
- Count matrices with library size information
- Recommended: Use TMM-normalized data

#### **limma Compatible**
- Log-transformed normalized data
- Recommended: Use VST or rlog normalized data

## Best Practices

### Data Organization

#### **Folder Structure Recommendation**
```
Project_Name/
├── 01_raw_data/
│   ├── count_matrix.csv
│   └── metadata.csv
├── 02_processed_data/
│   ├── filtered_counts.csv
│   ├── normalized_counts.csv
│   └── processing_log.txt
├── 03_qc_results/
│   ├── plots/
│   └── qc_report.html
└── 04_session_data/
    └── analysis_session.rds
```

#### **File Naming Convention**
- Use descriptive names with dates
- Include processing stage in filename
- Example: `counts_filtered_normalized_20250605.csv`

### Version Control

#### **Track Analysis Versions**
- Save parameter settings for each analysis
- Document any manual modifications
- Keep original data separate from processed data

#### **Reproducibility Checklist**
- ✓ Save all parameter settings
- ✓ Export session information
- ✓ Document any manual steps
- ✓ Keep analysis log/notes

### Quality Assurance

#### **Before Download**
- Review all QC plots for expected patterns
- Verify filtering parameters are appropriate
- Check normalization effectiveness
- Ensure sample metadata is complete

#### **After Download**
- Verify file integrity (check file sizes)
- Test loading data in target analysis software
- Document any format conversions needed

### Storage Considerations

#### **File Sizes**
- Raw count matrices: 10-500MB typically
- Normalized data: Similar to raw data
- QC plots: 1-10MB per plot
- Session data: Can be very large (>1GB)

#### **Backup Strategy**
- Keep multiple copies of processed data
- Store on different media/locations
- Include analysis parameters with backups

## Troubleshooting Export Issues

### Common Problems

#### **Download Fails**
- **Cause**: File too large, browser timeout
- **Solution**: Try smaller subsets, use direct file save

#### **Incorrect File Format**
- **Cause**: Format not suitable for downstream tool
- **Solution**: Check target software requirements

#### **Missing Data in Export**
- **Cause**: Processing step not completed
- **Solution**: Ensure all analysis steps are finished

### Browser-Specific Issues

#### **Chrome/Edge**
- Generally most reliable for downloads
- May block very large files

#### **Firefox**
- Good compatibility
- Check download settings for large files

#### **Safari**
- May have issues with certain file types
- Verify file associations

## Integration with Downstream Analysis

### Preparing for Differential Expression

#### **For DESeq2**
1. Export filtered, non-normalized counts
2. Export complete metadata
3. Use integer count values
4. Include all sample information

#### **For edgeR**
1. Export TMM-normalized data or raw counts
2. Include library size information
3. Export as DGEList object if possible

#### **For limma**
1. Export log-transformed normalized data
2. Use VST or rlog transformation
3. Include sample metadata with batch information

### Documentation for Collaborators

#### **Include with Data Exports**
- Analysis parameter summary
- QC report highlights
- Known issues or limitations
- Recommended next steps

---

**Related Documentation:**
- [Data Input & Validation](data_input.html) for initial data processing
- [QC Plots & Summaries](qc_plots_interpretation_guide.html) for quality assessment
- [Filtering & Normalization](filtering_normalization_guide.html) for data processing options 