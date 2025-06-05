# Technical Requirements

**Date:** 6/5/2025  
**Version:** 2.0.0

## System Requirements

### Minimum Requirements
- **Operating System**: Windows 10, macOS 10.14, or Linux (Ubuntu 18.04+)
- **R Version**: 4.0.0 or higher
- **RAM**: 4GB (8GB recommended for datasets >10,000 genes)
- **Storage**: 2GB free space (additional space needed for data and outputs)
- **Browser**: Chrome 90+, Firefox 88+, Safari 14+, or Edge 90+

### Recommended Requirements
- **R Version**: 4.3.0 or higher
- **RAM**: 8GB or more
- **CPU**: Multi-core processor (2+ cores)
- **Storage**: 5GB+ free space for large datasets

## Software Dependencies

### R Packages (Automatically Installed)

#### CRAN Packages
```r
# Shiny and UI
"shiny", "shinythemes", "shinyWidgets", "DT", "htmltools"

# Data manipulation
"dplyr", "tidyr", "readr"

# Visualization
"ggplot2", "plotly", "scales", "pheatmap", "ggrepel", 
"viridis", "RColorBrewer", "gridExtra", "cowplot"

# Statistical analysis
"moments", "e1071", "matrixStats", "preprocessCore"
```

#### Bioconductor Packages
```r
# Core infrastructure
"BiocGenerics", "S4Vectors", "IRanges", "GenomeInfoDb", 
"GenomicRanges", "SummarizedExperiment"

# RNA-seq analysis
"DESeq2", "edgeR", "limma"
```

### Installation Notes
- The application automatically checks for and installs missing packages
- First launch may take 10-15 minutes for package installation
- Internet connection required for initial setup

## Data Requirements

### Count Matrix Format
- **File Type**: CSV format
- **Structure**: Genes as rows, samples as columns
- **Gene IDs**: First column or row names (Ensembl IDs, Gene Symbols, etc.)
- **Values**: Non-negative integers (raw counts)
- **Size Limits**: Up to 1GB file size (configurable)

#### Example Structure:
```
Gene_ID,Sample1,Sample2,Sample3,Sample4
ENSG00000000003,1043,2156,987,1234
ENSG00000000005,0,1,0,2
ENSG00000000419,2145,3456,2098,2876
```

### Metadata Format
- **File Type**: CSV format
- **Structure**: Samples as rows, experimental factors as columns
- **Sample IDs**: Must match count matrix column names exactly
- **Factors**: Categorical variables (Treatment, Genotype, Time, etc.)

#### Example Structure:
```
Sample_ID,Treatment,Genotype,Batch
Sample1,Control,WT,1
Sample2,Treated,WT,1
Sample3,Control,KO,2
Sample4,Treated,KO,2
```

## Performance Considerations

### Dataset Size Guidelines
| Genes | Samples | RAM Needed | Processing Time |
|-------|---------|------------|-----------------|
| <10K  | <20     | 4GB        | <5 minutes      |
| 10K-20K | 20-50   | 8GB        | 5-15 minutes    |
| 20K+  | 50+     | 16GB+      | 15+ minutes     |

### Memory Usage
- **Data Loading**: ~3x file size in RAM
- **QC Analysis**: Additional 2-4GB for plots and calculations
- **Normalization**: 2-3x processed data size

### Browser Performance
- **Chrome/Edge**: Best performance for large datasets
- **Firefox**: Good performance, may be slower with 3D plots
- **Safari**: Limited performance with very large datasets

## Network Requirements

### Initial Setup
- Internet connection required for package installation
- ~500MB download for all dependencies

### Runtime
- No internet connection required after setup
- Local file access only

## Security Considerations

### Data Privacy
- All data processing occurs locally
- No data transmitted to external servers
- Files remain on your local system

### File Access
- Application reads files from selected directories only
- No system file access beyond user-selected folders
- Temporary files created in R session directory

## Troubleshooting System Issues

### Common Problems

#### Package Installation Failures
```r
# Clear package cache and reinstall
remove.packages("problematic_package")
install.packages("problematic_package")
```

#### Memory Errors
- Reduce dataset size
- Close other applications
- Increase virtual memory (Windows)
- Use data filtering before analysis

#### Browser Issues
- Clear browser cache
- Try different browser
- Disable browser extensions
- Check browser JavaScript settings

### Performance Optimization

#### For Large Datasets
1. **Pre-filter data**: Remove zero-count genes before upload
2. **Use subsets**: Analyze smaller groups separately
3. **Increase memory**: Close unnecessary applications
4. **Use batch mode**: Process multiple smaller files

#### R Configuration
```r
# Increase memory limit (Windows)
memory.limit(size = 8000)

# Optimize garbage collection
gc()

# Check memory usage
memory.size()
```

## Compatibility Notes

### R Version Compatibility
- **R 4.0.x**: Minimum supported version
- **R 4.1.x**: Fully supported
- **R 4.2.x**: Recommended
- **R 4.3.x+**: Latest features available

### Operating System Notes
- **Windows**: Full compatibility, use R 4.0+
- **macOS**: Full compatibility, may need Xcode tools
- **Linux**: Full compatibility, install system dependencies

### Known Issues
- **R 4.5.0**: DESeq2 compatibility issues resolved in app
- **macOS Big Sur+**: Potential graphics rendering issues (use X11)
- **Windows 7**: Not supported due to R version requirements

## Support and Updates

### Getting Help
1. Check this technical guide first
2. Review [Troubleshooting Guide](troubleshooting.html)
3. Ensure system meets minimum requirements
4. Try with example data to isolate issues

### Version Updates
- New versions may require package updates
- Check release notes for compatibility changes
- Backup important analysis results before updating 