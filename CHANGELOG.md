# Changelog

All notable changes to the RNA Processing APP project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] - 2025-06-03

### Added
- Complete QC & Pre-processing Tool (Module 1)
- Comprehensive data input and validation system
- Multiple filtering strategies (count-based, CPM-based, group-aware)
- Multiple normalization methods (CPM, TMM, RLE, VST, RLOG)
- Interactive quality control visualizations
- Real-time filtering and normalization impact visualization
- Export capabilities for processed data and plots
- Professional documentation suite
- User manual with quick start guide
- Developer documentation
- Troubleshooting guide
- Technical requirements documentation

### Features
- **Data Input & Validation**
  - Flexible CSV upload with configurable delimiters
  - Automatic duplicate gene detection and handling
  - Sample consistency validation
  - Support for various gene ID formats

- **Quality Control**
  - Library size distribution analysis
  - Gene expression distribution plots
  - Sample-to-sample correlation heatmaps
  - Interactive 2D and 3D PCA plots
  - Comprehensive QC summary reports

- **Filtering & Normalization**
  - Count-based filtering with customizable thresholds
  - CPM-based filtering with group awareness
  - Multiple normalization algorithms
  - Before/after comparison visualizations
  - Real-time impact assessment

- **Export & Integration**
  - Multiple export formats (CSV, RDS, PDF)
  - Publication-ready plots
  - Comprehensive analysis reports
  - Downstream analysis compatibility

### Technical Improvements
- Modular architecture for scalability
- Comprehensive error handling
- Automatic package installation
- Cross-platform compatibility
- Performance optimizations for large datasets

## [1.0.0] - 2024-05-29

### Added
- Initial release of QC & Pre-processing Tool
- Basic data input functionality
- Core quality control visualizations
- Simple filtering and normalization options
- Basic export capabilities

### Infrastructure
- Project structure establishment
- Initial documentation framework
- Basic testing structure

## [Unreleased]

### Planned
- Module 2: Differential Gene Expression Analysis Tool
- Module 3: Pathway & Functional Enrichment Analysis Tool
- Module 4: Advanced Visualization Tool
- Enhanced performance optimizations
- Batch processing capabilities
- Integration between modules

### Fixed
- **Sample-to-Sample Correlation Heatmap Display Issue** (06/03/2025 - Eren Ada, PhD)
  - Resolved duplicate `output$correlation_heatmap` definition that was preventing heatmaps from displaying
  - Fixed conflicting renderPlot handlers in server_tab2_qc_plots.R
  - Correlation heatmaps now display correctly with proper error handling and validation
  - Enhanced pheatmap rendering with explicit grid.draw() calls for proper Shiny compatibility
  - Added fallback base R heatmap for cases where pheatmap fails to render
  - Improved debugging and error reporting for correlation plot issues

- **Download Buttons UI Issue** (06/03/2025 - Eren Ada, PhD)
  - Fixed download buttons appearing on all sub-tabs of QC Plots & Summaries
  - Moved download buttons to their respective tabs where they are relevant
  - "Download PCA Plot" and "Download PCA Statistics" buttons now only appear on Sample Similarity Analysis tab
  - "Download Heatmap" button now only appears on Sample Correlation Analysis tab

### Removed
- **Download QC Report Functionality** (06/03/2025 - Eren Ada, PhD)
  - Removed `download_qc_report` button from QC Plots & Summaries tab
  - Removed corresponding server handler for QC report generation
  - Simplified download section to include only PCA plot and correlation heatmap downloads

### Technical Details
- Fixed duplicate output definition in `modules/module1_qc_preprocessing/R/server_tab2_qc_plots.R` lines 722 and 794
- Removed QC report UI button from `modules/module1_qc_preprocessing/R/ui_tab2_qc_plots.R`
- Maintained PCA statistics download functionality
- Enhanced correlation plot parameter validation and error messaging
- Improved UI organization by placing download buttons within relevant tabs only

---

**Author:** Eren Ada, PhD  
**Project:** Modular Bulk RNA-seq Analysis RShiny Tool Suite 