# RNA-seq QC & Pre-processing Tool - User Manual

**Author:** Eren Ada, PhD  
**Date:** 6/5/2025  
**Version:** 2.0.0

## Table of Contents

1. [Quick Start Guide](quick_start.html) - 15-minute walkthrough with example data
2. [Data Input & Validation](data_input.html) - File upload, validation, and duplicate handling
3. [Quality Control Plots](qc_plots_interpretation_guide.html) - Comprehensive QC assessment and interpretation
4. [Filtering & Normalization](filtering_normalization_guide.html) - Gene filtering and normalization methods
5. [Export & Download](export_download_guide.html) - Data export and report generation
6. [Troubleshooting](troubleshooting.html) - Common issues and solutions
7. [Technical Requirements](technical_requirements.html) - System requirements and setup

## Overview

The RNA-seq QC & Pre-processing Tool is a comprehensive Shiny application designed for quality control and pre-processing of bulk RNA-sequencing count data. This tool prepares clean, well-characterized data for downstream analysis in specialized differential expression and pathway analysis tools.

**Key Capabilities:**
- Comprehensive data validation and quality assessment
- Advanced filtering and normalization with evaluation metrics
- Interactive visualizations with publication-quality exports
- Detailed reporting for reproducibility and documentation

## Comprehensive Feature Set

### Data Handling & Validation
- **Flexible Import**: CSV, TSV file support with automatic format detection
- **Comprehensive Validation**: Sample consistency, data type, and integrity checks
- **Intelligent Duplicate Handling**: Four strategies including sum, keep highest, rename, or remove
- **Automatic Processing**: Decimal rounding, missing value detection, format standardization
- **Quality Reporting**: Detailed validation summaries with actionable recommendations

### Advanced Quality Control
- **Library Size Analysis**: Distribution plots with log scaling and outlier detection
- **Gene Detection Rates**: Configurable thresholds with quality benchmarks
- **Interactive PCA**: 2D/3D plots with PC selection (PC1-PC10) and metadata coloring
- **Correlation Analysis**: Multiple methods (Pearson/Spearman) with data source options
- **Sample Filtering**: Dynamic metadata-based filtering for focused analysis
- **Statistical Assessment**: Comprehensive normality testing with multiple methods

### Sophisticated Filtering & Normalization
- **Group-Aware Filtering**: Experimental design-based gene filtering strategies
- **Multiple Expression Units**: Raw counts vs. CPM-based thresholds
- **Comprehensive Normalization**: Eight methods including DESeq2, TMM, VST, rlog
- **Evaluation Framework**: CV reduction, mean-variance relationships, batch effect assessment
- **Distribution Analysis**: Q-Q plots, density analysis, and statistical testing
- **Method Recommendations**: Data-driven guidance for optimal normalization choice

### Professional Export & Reporting
- **Timestamped Data Export**: Filtered/normalized counts with processing provenance
- **Comprehensive Reports**: HTML analysis reports with session information
- **Publication-Quality Plots**: PDF vector graphics for all visualizations
- **Statistical Data**: CSV exports of PCA statistics and correlation matrices
- **Reproducibility Support**: Complete parameter tracking and session documentation

## Detailed Workflow

### Phase 1: Data Input & Validation (5-10 minutes)
1. **File Upload**: Drag-and-drop or browse for count matrix and metadata files
2. **Automatic Validation**: Format detection, data type validation, sample consistency
3. **Duplicate Resolution**: Interactive handling of duplicate gene IDs with preview
4. **Quality Assessment**: Initial library size and detection rate evaluation
5. **Data Export**: Original and processed data download options

### Phase 2: Quality Control Analysis (10-15 minutes)
1. **Basic Metrics**: Library size distributions and gene detection rates
2. **Sample Similarity**: PCA analysis with interactive 2D/3D visualizations
3. **Correlation Assessment**: Heatmap generation with multiple data sources
4. **Outlier Identification**: Statistical and visual outlier detection
5. **Export Options**: High-resolution plots and statistical summaries

### Phase 3: Filtering & Normalization (10-20 minutes)
1. **Strategic Filtering**: Group-aware gene filtering with impact visualization
2. **Method Selection**: Choose from eight normalization approaches
3. **Comprehensive Evaluation**: Multi-metric assessment of normalization effectiveness
4. **Distribution Analysis**: Detailed normality testing and recommendations
5. **Final Export**: Processed data with complete analysis documentation

## Quality Assessment Framework

### Data Quality Indicators

**Excellent Quality:**
- Library sizes within 2-fold range
- Gene detection rates >70%
- Sample correlations >0.85
- Clear biological clustering in PCA
- Successful normalization (CV reduction >70%)

**Good Quality:**
- Library sizes within 5-fold range
- Gene detection rates 60-70%
- Sample correlations >0.75
- Identifiable biological patterns
- Moderate normalization improvement (50-70% CV reduction)

**Concerning Quality:**
- Library sizes >10-fold differences
- Gene detection rates <50%
- Sample correlations <0.7
- No clear biological patterns or technical clustering
- Poor normalization effectiveness (<30% CV reduction)

### Decision Support Framework

**Filtering Recommendations:**
- **Conservative (Recommended)**: ≥10 counts in ≥2 samples
- **Moderate**: ≥5 counts with group-aware filtering
- **Liberal**: ≥1 count for minimal filtering
- **Target Retention**: 60-80% of original genes

**Normalization Guidance:**
- **Most RNA-seq datasets**: DESeq2 (Median of Ratios)
- **edgeR workflows**: TMM normalization
- **Visualization/clustering**: VST (>30 samples) or rlog (<30 samples)
- **Exploratory analysis**: CPM with log transformation

## Advanced Features

### Interactive Analysis
- **Real-time Updates**: QC plots automatically refresh after normalization
- **Dynamic Filtering**: Live preview of filtering impact on gene counts
- **Custom Coloring**: PCA and correlation plots colored by any metadata variable
- **Zoom & Selection**: Interactive plot exploration with plotly integration

### Statistical Rigor
- **Multiple Testing**: Comprehensive normality assessment with three test methods
- **Significance Levels**: Configurable α levels (0.05, 0.01, 0.001)
- **Correlation Significance**: Statistical testing with visual significance indicators
- **Batch Effect Analysis**: Quantitative assessment of technical factor influence

### Reproducibility Features
- **Parameter Tracking**: Complete log of all analysis settings
- **Session Documentation**: R version, package versions, and system information
- **Timestamped Outputs**: All exports include processing date/time
- **Analysis Reports**: HTML summaries with embedded statistics and parameters

## Integration & Compatibility

### Downstream Tool Preparation
- **DESeq2**: Filtered, non-normalized integer counts with metadata
- **edgeR**: TMM-normalized data or raw counts with library size information
- **limma**: Log-transformed data (VST/rlog) with batch information
- **Custom Workflows**: Flexible CSV exports compatible with any analysis pipeline

### File Format Standards
- **Universal CSV**: Gene IDs as first column, clean headers, no row names
- **Metadata Preservation**: Complete sample annotation export
- **Statistical Summaries**: Structured data for meta-analysis or reporting

## Performance & Scalability

### System Optimization
- **Memory Management**: Efficient handling of large datasets (>20K genes, >50 samples)
- **Processing Speed**: Optimized algorithms with progress indicators
- **Browser Compatibility**: Full functionality across modern browsers
- **Interactive Performance**: Responsive UI even with complex datasets

### Dataset Guidelines
| Genes | Samples | RAM Needed | Processing Time | Recommendation |
|-------|---------|------------|-----------------|----------------|
| <10K  | <20     | 4GB        | <5 minutes      | Standard settings |
| 10K-20K | 20-50   | 8GB        | 5-15 minutes    | Monitor progress |
| 20K+  | 50+     | 16GB+      | 15+ minutes     | Consider filtering |

## Getting Started

### New Users
1. **Start Here**: [Quick Start Guide](quick_start.html) - Complete 15-minute tutorial
2. **System Setup**: [Technical Requirements](technical_requirements.html) - Installation and configuration
3. **Example Data**: Use provided datasets in `example_data/` directory

### Experienced Users
1. **Advanced Features**: [Filtering & Normalization Guide](filtering_normalization_guide.html)
2. **Quality Assessment**: [QC Plots Interpretation Guide](qc_plots_interpretation_guide.html)
3. **Data Management**: [Export & Download Guide](export_download_guide.html)

## Support & Resources

### Documentation
- **Comprehensive Guides**: Detailed documentation for each analysis phase
- **Best Practices**: Evidence-based recommendations for method selection
- **Troubleshooting**: Common issues with step-by-step solutions
- **Method References**: Links to original publications and validation studies

### Community Support
- **Contact**: Eren Ada, PhD - erenada@gmail.com
- **GitHub Repository**: https://github.com/erenada
- **Issue Reporting**: Bug reports and feature requests welcome
- **Collaboration**: Open to community contributions and feedback

### Quality Assurance
- **Validation Testing**: Extensive testing with diverse datasets
- **Method Verification**: Comparison with established workflows
- **Continuous Improvement**: Regular updates based on user feedback
- **Documentation Maintenance**: Synchronized with software updates

---

**Production Ready**: This tool has been thoroughly tested and is ready for research use. All methods are based on established best practices in the RNA-seq analysis community.

*This documentation corresponds to RNA-seq QC & Pre-processing Tool v2.0.0* 