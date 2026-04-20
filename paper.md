---
title: 'RNA-seq Quality Control and Pre-processing Tool: An Interactive Shiny Application for Bulk RNA-seq Data Validation and Normalization'
tags:
  - R
  - RNA-seq
  - quality control
  - transcriptomics
  - bioinformatics
  - Shiny
authors:
  - name: Eren Ada
    orcid: 0000-0001-9141-8282
    affiliation: 1
affiliations:
  - name: Department of Immunology, Harvard Medical School
    index: 1
date: 2026-04-20
bibliography: paper.bib
---

# Summary

The `RNA_QC_APP` is an interactive R Shiny application that enables researchers to perform rigorous quality control (QC) and pre-processing of bulk RNA sequencing (RNA-seq) count data without writing code. Through a browser-based graphical interface, users upload a count matrix and a sample metadata file, and the application guides them through a structured workflow: data validation, exploratory QC visualization (library size distributions, expression density plots, sample-to-sample correlation heatmaps, and two- and three-dimensional principal component analysis), low-count gene filtering, normalization with several established algorithms, and export of processed matrices and publication-ready figures. The application wraps validated Bioconductor methods [@huber2015orchestrating] such as `DESeq2` [@love2014deseq2] and `edgeR` [@robinson2010edger] behind a consistent interface, so that investigators can inspect the effect of each step in real time and produce reproducible pre-processing outputs that are ready for downstream differential expression or pathway analyses.

# Statement of Need

Bulk RNA-seq is now a standard tool in biology and medicine, yet rigorous pre-processing remains a nontrivial barrier. Count validation, outlier detection, filtering of low-expression genes, and normalization strongly influence the statistical power and interpretability of downstream analyses, but they typically require proficiency in R or Python, familiarity with Bioconductor workflows, and enough statistical background to choose among methods such as `DESeq2`'s median-of-ratios, `edgeR`'s trimmed mean of M-values (TMM) [@robinson2010scaling], and variance-stabilizing transformations [@love2014deseq2]. Biologists, immunologists, and clinicians who generate or receive RNA-seq data often lack the time or training to assemble such a workflow from scripts.

Existing tools address parts of this problem but leave a gap. `DESeq2` [@love2014deseq2], `edgeR` [@robinson2010edger], and `limma` [@ritchie2015limma] are authoritative but are command-line R packages. `MultiQC` [@ewels2016multiqc] aggregates upstream alignment and read-level QC reports, but does not offer interactive count-matrix QC, filtering, or normalization. `RSeQC` [@wang2012rseqc] similarly targets read-level QC on BAM files. `RNA_QC_APP` fills this gap by providing an integrated, no-code environment that takes the user from a raw count matrix through to a validated, filtered, and normalized dataset, with interactive diagnostics at every stage. The target audience is wet-lab biologists and clinician-scientists who need defensible, reproducible pre-processing without adopting a scripting workflow, as well as educators teaching RNA-seq analysis. This application was developed in direct response to recurring requests from wet-lab collaborators at Harvard Medical School who generate bulk RNA-seq data but do not maintain scripting workflows, and it has since been used internally to support multiple ongoing immunology projects.

# State of the Field

Several Shiny and web-based tools interactively explore RNA-seq data. `iSEE` [@rue2018isee] provides a highly configurable interface for exploring `SummarizedExperiment` objects after analysis. `pcaExplorer` [@marini2019pcaexplorer] focuses on interactive principal component analysis and functional interpretation of RNA-seq results, and `DEBrowser` [@kucukural2019debrowser] offers a Shiny interface for differential expression analysis and visualization. These tools excel at post-hoc exploration but generally assume that pre-processing (count validation, filtering, normalization choice) has already been performed by a bioinformatician.

`RNA_QC_APP` is distinct in two respects. First, it is explicitly dedicated to the pre-processing stage, bringing validation, QC visualization, filtering, and normalization into a single linear workflow with consistent state management. Second, it is designed as the first component of a modular Shiny suite; the processed output is structured to be consumed directly by a companion downstream application (`RNA_DEG_APP`) for differential expression analysis, so that an end-to-end no-code pipeline is possible without leaving the browser.

# Software Design

`RNA_QC_APP` is implemented in R using the Shiny framework [@chang2023shiny] and follows a modular architecture. The top-level `app.R` wires together three core files in `R/` (`app_ui.R`, `app_server.R`, `modules_source.R`) that assemble the interface from a set of Shiny modules under `modules/module1_qc_preprocessing/R/`. Each analytical tab (input and validation, QC plots, filtering and normalization, about) has paired `ui_*` and `server_*` files, and shared logic is factored into utility files (`utils_filtering.R`, `utils_normalization.R`, `utils_plotting.R`, `utils_qc_plotting.R`). This separation keeps individual components under a manageable size, allows isolated testing of each module with `testthat`, and makes the application extensible to additional stages of analysis.

State is coordinated through a single `reactiveValues` container (`shared_data`) in the server that tracks the processed count matrix, the processed metadata, normalization status, and QC completion flags. Reactive programming propagates changes across tabs: uploading new data invalidates downstream views, and navigation is guarded so that users cannot reach filtering and normalization until validation and QC are complete. This reactive pipeline mirrors the scientific workflow and prevents common user errors.

The application accepts comma- or tab-separated count matrices (genes in rows, samples in columns) and a matching sample metadata table. Validation logic detects duplicate gene identifiers, non-numeric entries, and mismatches between sample columns and metadata rows, and offers configurable strategies to resolve each. QC visualizations are rendered with `ggplot2` [@wickham2016ggplot2] and made interactive with `plotly` [@sievert2020plotly]. Filtering supports count- and CPM-based thresholds with optional group awareness, and normalization wraps library-size scaling, upper-quartile, RLE, CPM, TMM (via `edgeR`), median-of-ratios and variance-stabilizing transformations (via `DESeq2`), rlog (via `DESeq2`), and quantile normalization (via `preprocessCore`). UI polish is provided by `shinythemes`, `shinyWidgets`, and `DT`. Outputs include normalized matrices as CSV or RDS, and publication-ready figures as PDF.

# Research Impact Statement

`RNA_QC_APP` was developed in the Department of Immunology at Harvard Medical School in response to recurring requests from wet-lab collaborators who generate bulk RNA-seq data but do not maintain a scripting workflow. By making validated Bioconductor methods accessible through a graphical interface, the application lowers the technical barrier to reproducible pre-processing and allows investigators to take ownership of their own QC decisions rather than outsourcing them. Packaged example data (a 5,001-gene by 24-sample count matrix with matching metadata) shortens the onboarding time for new users and for teaching settings. Although developed locally, the tool is open-source and designed for global use in immunology, cancer biology, clinical transcriptomics, and other fields that rely on bulk RNA-seq. It is the upstream component of a modular suite; the downstream companion application, `RNA_DEG_APP`, consumes its output for differential expression analysis, enabling an end-to-end no-code pipeline from raw counts to differentially expressed genes.

# AI Usage Disclosure

AI-assisted coding tools, including GitHub Copilot, OpenAI ChatGPT, Anthropic Claude, and the Cursor editor, were used during development of `RNA_QC_APP` and during the drafting of this manuscript. Their role was limited to code completion, refactoring suggestions, documentation scaffolding, and prose editing. All scientific decisions, algorithmic choices, and interpretations of results remain the responsibility of the author. Every code change and every section of this paper was reviewed, tested, and validated by the author prior to inclusion, and the final work reflects the author's judgment.

# Acknowledgements

The author thanks the Bioconductor community [@huber2015orchestrating], and in particular the authors and maintainers of `DESeq2` [@love2014deseq2], `edgeR` [@robinson2010edger], and `limma` [@ritchie2015limma], whose methods form the statistical backbone of this application. The author also thanks the Shiny [@chang2023shiny] and `plotly` [@sievert2020plotly] developer communities, and members of the Department of Immunology at Harvard Medical School for feedback during early testing.

# References
