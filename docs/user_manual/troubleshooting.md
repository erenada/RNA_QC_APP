# Troubleshooting Guide

**Date:** 6/5/2025  
**Version:** 2.0.0

This guide addresses common issues and their solutions. For technical requirements, see [Technical Requirements](technical_requirements.html).

## Quick Diagnostic Checklist

Before diving into specific issues:
1. ✓ Check that your system meets [minimum requirements](technical_requirements.html)
2. ✓ Ensure you have a stable internet connection for initial setup
3. ✓ Verify your data files are in proper CSV format
4. ✓ Close unnecessary applications to free up memory
5. ✓ Try with example data to isolate the issue

## Application Launch Issues

### Problem: App Won't Start
**Symptoms:** Error messages when running `shiny::runApp("app.R")`

**Solutions:**
1. **Check R Version:**
   ```r
   R.version.string  # Should be 4.0.0 or higher
   ```

2. **Update R if needed:**
   - Download latest R from [r-project.org](https://www.r-project.org/)
   - Restart RStudio after installation

3. **Clear Package Cache:**
   ```r
   # Remove problematic packages
   remove.packages(c("shiny", "DESeq2"))
   # Restart R session, then run the app again
   ```

### Problem: Package Installation Fails
**Symptoms:** Error messages about missing packages or installation failures

**Solutions:**
1. **Check Internet Connection:**
   ```r
   # Test CRAN connection
   available.packages()[1:5, 1:2]
   ```

2. **Manual Package Installation:**
   ```r
   # Install Bioconductor manager first
   if (!require("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   
   # Install core packages manually
   BiocManager::install(c("DESeq2", "edgeR"))
   install.packages(c("shiny", "DT", "plotly"))
   ```

3. **Clear Library and Reinstall:**
   ```r
   # Check library location
   .libPaths()
   # Remove and reinstall problematic packages
   ```

### Problem: Browser Doesn't Open
**Symptoms:** App starts but browser doesn't launch automatically

**Solutions:**
1. **Manual Browser Launch:**
   - Look for the message "Listening on http://127.0.0.1:XXXX"
   - Copy this URL into your browser manually

2. **Specify Browser:**
   ```r
   options(browser = "chrome")  # or "firefox", "safari"
   shiny::runApp("app.R")
   ```

## Data Upload Issues

### Problem: File Upload Fails
**Symptoms:** "File upload failed" or empty file display

**Solutions:**
1. **Check File Size:**
   - Maximum file size is 1GB by default
   - For larger files, modify `options(shiny.maxRequestSize = ...)` in `app.R`

2. **Verify File Format:**
   - Must be CSV format
   - Check for proper delimiters (comma, semicolon, tab)
   - Ensure consistent encoding (UTF-8 recommended)

3. **Test File Locally:**
   ```r
   # Test if R can read your file
   test_data <- read.csv("your_file.csv", row.names = 1)
   head(test_data)
   ```

### Problem: Data Validation Errors
**Symptoms:** Red error messages in validation output

**Common Issues and Fixes:**

#### Sample Name Mismatches
```
Error: Sample names don't match between count matrix and metadata
```
**Solution:**
- Ensure exact spelling match (case-sensitive)
- Check for extra spaces or special characters
- Use consistent naming convention

#### Missing or Invalid Gene IDs
```
Error: Duplicate gene names found
```
**Solution:**
- Choose "Make unique" option for duplicate handling
- Or pre-process your data to handle duplicates manually

#### Non-numeric Count Data
```
Error: Non-numeric values detected in count matrix
```
**Solution:**
- Check for text values in count columns
- Ensure all count values are non-negative integers
- Remove any NA or missing values

### Problem: Metadata Issues
**Symptoms:** Errors during metadata processing

**Solutions:**
1. **Check Column Names:**
   - No special characters in column names
   - Use underscores instead of spaces
   - Avoid starting with numbers

2. **Verify Factor Levels:**
   - Ensure categorical variables have consistent spelling
   - Check for missing values (NA)
   - Use simple naming (avoid special characters)

## Quality Control Issues

### Problem: QC Plots Don't Generate
**Symptoms:** Blank plots or error messages in QC tab

**Solutions:**
1. **Memory Issues:**
   ```r
   # Check available memory
   memory.size()  # Windows
   # Close other applications and try again
   ```

2. **Data Size Issues:**
   - Try with a subset of your data first
   - Reduce number of genes if >50,000

3. **Graphics Device Issues:**
   ```r
   # Reset graphics device
   dev.off()
   # Try generating plots again
   ```

### Problem: PCA Plots Show Unexpected Patterns
**Symptoms:** Samples don't cluster as expected

**This is often biological, not technical. Consider:**
1. **Batch Effects:** Look for clustering by batch/processing date
2. **Outlier Samples:** Identify samples with low correlation
3. **Experimental Design:** Verify metadata matches experimental setup

### Problem: Interactive Plots Don't Work
**Symptoms:** 3D PCA plots or plotly visualizations don't respond

**Solutions:**
1. **Browser Compatibility:**
   - Use Chrome or Firefox for best performance
   - Enable JavaScript in browser settings
   - Disable browser extensions temporarily

2. **Memory Issues:**
   - Reduce dataset size for interactive plotting
   - Close other browser tabs

## Filtering and Normalization Issues

### Problem: Too Many Genes Filtered Out
**Symptoms:** >80% of genes removed during filtering

**Solutions:**
1. **Adjust Filter Parameters:**
   - Lower minimum count threshold
   - Reduce minimum sample requirement
   - Use group-aware filtering for complex designs

2. **Check Data Quality:**
   - Very low counts might indicate sequencing issues
   - Review library sizes in QC plots

### Problem: Normalization Fails
**Symptoms:** Error messages during normalization step

**Solutions:**
1. **Check for Zero-Count Samples:**
   ```r
   # In R console, check your data
   colSums(your_count_data) == 0  # Should all be FALSE
   ```

2. **Memory Issues:**
   - Try simpler normalization methods first (CPM)
   - Process smaller subsets if needed

3. **Package Issues:**
   ```r
   # Reinstall normalization packages
   BiocManager::install(c("DESeq2", "edgeR", "limma"))
   ```

## Performance Issues

### Problem: App Runs Very Slowly
**Symptoms:** Long processing times, unresponsive interface

**Solutions:**
1. **Reduce Data Size:**
   - Pre-filter genes with zero counts across all samples
   - Use subset of samples for initial exploration

2. **Optimize System:**
   - Close unnecessary applications
   - Increase virtual memory (Windows)
   - Use latest R version

3. **Browser Optimization:**
   - Use Chrome for best performance
   - Clear browser cache
   - Disable unnecessary browser extensions

### Problem: Memory Errors
**Symptoms:** "Cannot allocate vector of size..." errors

**Solutions:**
1. **Increase Memory Limit (Windows):**
   ```r
   memory.limit(size = 8000)  # Set to 8GB
   ```

2. **Optimize R Session:**
   ```r
   # Clean environment
   rm(list = ls())
   gc()  # Garbage collection
   ```

3. **Use Data Subsets:**
   - Analyze smaller portions of data
   - Filter data before upload

## Export and Download Issues

### Problem: Downloads Don't Start
**Symptoms:** Click download buttons but nothing happens

**Solutions:**
1. **Browser Settings:**
   - Check if downloads are blocked
   - Verify download folder permissions
   - Try different browser

2. **Large File Issues:**
   - Large files may take time to prepare
   - Check browser's download manager
   - Ensure sufficient disk space

### Problem: Exported Files Are Corrupted
**Symptoms:** Cannot open downloaded files

**Solutions:**
1. **Wait for Complete Download:**
   - Ensure download fully completes before opening
   - Check file size matches expected size

2. **File Format Issues:**
   - Try opening with appropriate software
   - Check file extension matches content

## Getting Additional Help

### Information to Provide When Seeking Help
1. **System Information:**
   - Operating system and version
   - R version (`R.version.string`)
   - Browser and version

2. **Error Details:**
   - Exact error message
   - When the error occurs
   - Steps to reproduce

3. **Data Information:**
   - Dataset size (genes × samples)
   - File formats used
   - Any special characteristics of your data

### Self-Diagnosis Steps
1. **Try with Example Data:**
   - Use provided example files to isolate issues
   - If examples work, the issue is likely data-specific

2. **Check R Console:**
   - Look for additional error messages in R console
   - These often provide more detailed information

3. **Simplify the Analysis:**
   - Start with smaller datasets
   - Use basic settings first, then add complexity

### When to Start Over
Sometimes it's faster to restart:
1. Restart R session completely
2. Clear browser cache and cookies
3. Try with fresh data files
4. Use minimal settings to start

Remember: Many "errors" are actually data quality issues that the app is correctly identifying. Review the validation messages carefully before assuming there's a technical problem. 