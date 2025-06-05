# Data Input & Validation Guide

**Author:** Eren Ada, PhD  
**Date:** 6/4/2025  
**Version:** 2.0.0

## Overview

The Data Input & Validation tab is the first step in the RNA-seq QC & Pre-processing workflow. This tab handles the upload, validation, and initial processing of your count matrix and sample metadata files. It ensures data integrity and consistency before proceeding to quality control analysis.

## Table of Contents

1. [File Requirements](#file-requirements)
2. [Upload Process](#upload-process)
3. [Data Validation](#data-validation)
4. [Duplicate Gene Handling](#duplicate-gene-handling)
5. [Common Issues & Solutions](#common-issues--solutions)
6. [Data Processing Options](#data-processing-options)
7. [Export & Download](#export--download)
8. [Best Practices](#best-practices)

## File Requirements

### Count Matrix Requirements

Your count matrix file must meet the following specifications:

#### **File Format**
- **Supported formats**: CSV (.csv) and TSV (.tsv, .txt)
- **Separator detection**: Automatic (comma, tab, or semicolon)
- **Encoding**: UTF-8 (recommended)

#### **Structure Requirements**
- **First column**: Gene IDs or gene names
  - Must be unique (duplicates will be detected and handled)
  - Can be any identifier (Ensembl, Gene Symbol, etc.)
  - No missing values allowed
- **Column headers**: Sample names as column headers
  - Must match sample names in metadata file exactly
  - No spaces or special characters recommended
- **Data values**: Raw count data
  - Must be non-negative integers
  - Decimal values will be automatically rounded
  - No missing values (NA, NULL, empty cells) allowed

#### **Example Count Matrix Format**
```
Gene_ID,Sample1,Sample2,Sample3,Sample4
GENE1,123,456,789,234
GENE2,234,567,890,345
GENE3,345,678,901,456
GENE4,0,12,45,23
```

### Metadata Requirements

Your sample metadata file must include:

#### **File Format**
- **Supported formats**: CSV (.csv) and TSV (.tsv, .txt)
- **Separator detection**: Automatic
- **Encoding**: UTF-8 (recommended)

#### **Structure Requirements**
- **First column**: Sample names (must match count matrix headers exactly)
- **Additional columns**: Sample attributes and experimental factors
  - Treatment groups, conditions, batch information, etc.
  - Each sample must have complete information
  - No missing values in critical columns

#### **Example Metadata Format**
```
Sample,Treatment,Condition,Batch,Age
Sample1,Control,WT,1,8
Sample2,Treatment,KO,1,8
Sample3,Control,WT,2,12
Sample4,Treatment,KO,2,12
```

## Upload Process

### Step 1: File Upload

1. **Upload Count Matrix**
   - Click "Upload Count Matrix" button
   - Select your count matrix file (.csv, .tsv, or .txt)
   - File format will be automatically detected

2. **Upload Sample Metadata**
   - Click "Upload Sample Metadata" button
   - Select your metadata file
   - Ensure sample names match count matrix exactly

### Step 2: Initial Validation

1. **Click "Validate Files"**
   - Initiates comprehensive validation process
   - Automatic format detection and file parsing
   - Data integrity checks and validation reports

2. **Review Detection Results**
   - File formats and separators detected
   - Any automatic corrections applied (e.g., rounding decimals)
   - Initial data summary statistics

## Data Validation

### Automatic Validation Checks

The application performs comprehensive validation including:

#### **1. File Format Validation**
- **Separator Detection**: Automatically detects comma, tab, or semicolon separators
- **Format Consistency**: Ensures consistent formatting throughout files
- **Character Encoding**: Validates proper file encoding

#### **2. Data Type Validation**
- **Integer Counts**: Ensures count data contains non-negative integers
- **Decimal Handling**: Automatically rounds decimal values with notification
- **Missing Values**: Checks for and reports any missing data

#### **3. Sample Consistency Check**
- **Name Matching**: Verifies sample names match between count matrix and metadata
- **Sample Count**: Ensures same number of samples in both files
- **Detailed Reporting**: Provides specific information about mismatches

#### **4. Data Integrity Checks**
- **Duplicate Detection**: Identifies duplicate gene IDs
- **Zero-Count Genes**: Reports genes with zero counts across all samples
- **Matrix Dimensions**: Validates data structure and dimensions

### Validation Results Display

After validation, you'll see:

- **✅ Success Messages**: Confirmation of successful validation steps
- **⚠️ Warnings**: Non-critical issues that were automatically resolved
- **❌ Errors**: Critical issues requiring user intervention
- **📊 Data Summary**: Overview of dataset characteristics

## Duplicate Gene Handling

### Detection and Identification

When duplicate gene IDs are detected, the application provides:

- **Duplicate Count**: Number of genes with duplicate entries
- **Preview Table**: Shows all duplicate entries with expression levels
- **Impact Assessment**: Estimates effect on downstream analysis

### Resolution Strategies

Choose from four duplicate handling methods:

#### **1. Sum Counts for Duplicates (Recommended)**
- **Description**: Combines counts from all duplicate entries
- **Use Case**: When duplicates represent the same biological entity
- **Result**: Single entry per gene with summed expression values
- **Example**: 
  ```
  GENE1_entry1: [100, 200, 150]
  GENE1_entry2: [50, 100, 75]
  Result: GENE1: [150, 300, 225]
  ```

#### **2. Keep Entry with Highest Average Expression**
- **Description**: Retains only the duplicate with highest mean expression
- **Use Case**: When one entry represents the "best" measurement
- **Result**: Single entry per gene (highest expressing variant)
- **Advantage**: Preserves original count distributions

#### **3. Rename with Suffix**
- **Description**: Makes gene names unique by adding suffixes
- **Use Case**: When you want to keep all entries as separate entities
- **Result**: GENE1 becomes GENE1_1, GENE1_2, etc.
- **Note**: Increases total gene count

#### **4. Remove All Duplicate Entries**
- **Description**: Eliminates all entries for genes with duplicates
- **Use Case**: Conservative approach when unsure about duplicate source
- **Result**: Complete removal of problematic genes
- **Caution**: May result in significant data loss

### Making the Choice

**Recommended approach**:
1. **Examine duplicate preview table** to understand the nature of duplicates
2. **Consider biological context** (e.g., alternative splicing, gene families)
3. **Check expression levels** - are they similar or very different?
4. **Default to "Sum Counts"** for most RNA-seq datasets

## Common Issues & Solutions

### File Upload Problems

#### **Issue**: File not recognized or fails to upload
**Solutions**:
- Ensure file has correct extension (.csv, .tsv, .txt)
- Check file size (very large files may take time)
- Verify file is not corrupted
- Try saving file in a different format

#### **Issue**: Separator not detected correctly
**Solutions**:
- Manually verify your file uses standard separators (comma, tab)
- Avoid mixed separators within the same file
- Check for hidden characters or unusual encoding

### Data Format Problems

#### **Issue**: Sample name mismatches
**Solutions**:
- Compare sample names in both files character by character
- Check for extra spaces, case differences, or special characters
- Use consistent naming conventions
- Verify no duplicate sample names exist

#### **Issue**: Non-integer count values
**Solutions**:
- The application automatically rounds decimal values
- Review the notification to understand what was changed
- Consider whether decimal values indicate a problem with data generation

#### **Issue**: Missing values (NA/NULL)
**Solutions**:
- Identify source of missing values in original data
- Replace missing values with 0 for counts (if appropriate)
- Remove samples or genes with excessive missing data before upload

### Validation Failures

#### **Issue**: No overlapping sample names
**Solutions**:
- Check sample naming conventions in both files
- Verify the correct columns are being used for sample identification
- Look for systematic differences (e.g., prefixes, suffixes)

#### **Issue**: Excessive duplicate genes
**Solutions**:
- Investigate source of duplicates (gene annotation issues?)
- Consider if this is expected for your data type
- Choose appropriate duplicate handling strategy

## Data Processing Options

### Zero-Count Gene Removal

**Option**: Remove genes with zero counts across all samples
- **Default**: Not selected (genes retained)
- **When to use**: To reduce dataset size and focus on expressed genes
- **Impact**: Reduces total gene count, may improve computational efficiency

### Additional Processing Features

- **Original Data Download**: Download unprocessed count matrix
- **Automatic Rounding**: Non-integer values rounded to integers
- **Format Standardization**: Consistent data structure for downstream analysis

## Export & Download

### Available Downloads

#### **1. Original Count Matrix**
- **Content**: Unprocessed data as uploaded
- **Format**: CSV with gene IDs as first column
- **Use**: Backup or alternative analysis pipelines

#### **2. Processed Count Matrix**
- **Content**: Data after duplicate handling and processing
- **Format**: CSV ready for downstream analysis
- **Use**: Input for other RNA-seq analysis tools

#### **3. Duplicate Statistics**
- **Content**: Detailed information about duplicate genes
- **Format**: CSV with duplicate gene analysis
- **Use**: Documentation and quality assessment

### Data Previews

The interface provides:
- **Count Matrix Preview**: First 10 rows and all columns
- **Metadata Preview**: All metadata with search and filter options
- **Library Size Plot**: Interactive visualization of total counts per sample
- **Processing Summary**: Detailed validation and processing results

## Best Practices

### File Preparation

1. **Consistent Naming**: Use consistent, descriptive sample names
2. **Clean Data**: Remove or handle missing values before upload
3. **Standard Formats**: Use standard CSV/TSV formats with UTF-8 encoding
4. **Backup Originals**: Keep copies of original files
5. **Document Changes**: Note any preprocessing steps applied

### Quality Checks

1. **Review All Validation Messages**: Address warnings and errors
2. **Examine Duplicate Genes**: Understand the source and choose appropriate handling
3. **Check Sample Consistency**: Verify all samples have matching metadata
4. **Validate Library Sizes**: Look for extreme outliers in the library size plot

### Data Processing Decisions

1. **Duplicate Handling**: Start with "Sum Counts" unless you have specific reasons for other approaches
2. **Zero-Count Genes**: Consider biological relevance before removal
3. **Sample Exclusions**: Note any samples that might need exclusion based on validation results

### Documentation

1. **Save Processing Summary**: Download and keep validation reports
2. **Record Decisions**: Document duplicate handling and other processing choices
3. **Version Control**: Keep track of different versions of your processed data

## Troubleshooting

### Performance Issues

- **Large Files**: Allow extra time for upload and processing of large datasets
- **Memory Usage**: Close other applications if experiencing slow performance
- **Browser Issues**: Try refreshing or using a different browser

### Data Quality Concerns

- **Low Sample Numbers**: Ensure you have adequate replicates for your experimental design
- **High Duplicate Rates**: Investigate potential issues with gene annotation or data generation
- **Extreme Library Size Variation**: Consider whether batch effects or technical issues are present

### Getting Help

If you encounter persistent issues:
1. Check the [Troubleshooting Guide](troubleshooting.html)
2. Review the [Technical Requirements](technical_requirements.html)
3. Try with the provided example data to isolate the issue
4. Document specific error messages and dataset characteristics

---

**Next Step**: After successful data validation and processing, proceed to the [QC Plots & Summaries](qc_plots_interpretation_guide.html) tab to assess data quality through comprehensive visualizations. 