
This repository contains the bioinformatics pipelines and tools developed and applied to study the role of the small nuclear RNA in metastatic castration-resistant prostate cancer (mCRPC). The study explores minor intron retention, alternative splicing, and its transcriptomic landscape, to assess its potential as a therapeutic target.

## Overview

Metastatic castration-resistant prostate cancer (mCRPC) is a significant challenge in oncology due to its resistance to androgen deprivation therapy (ADT) and the severe side effects associated with newer therapies. This project investigates U6atac, an snRNA involved in minor intron splicing, as a potential selective therapeutic target. The study uses two primary pipelines:

1. **Intron Retention & Alternative Splicing Pipeline**  
   A pipeline for analyzing minor intron retention and alternative splicing, optimized for RNA sequencing data. (Private)
   
2. **Spatial Omics Pipeline**  
   A workflow for analyzing transcriptomic data from patient-derived tissue samples using GeoMx Digital Spatial Profiling (DSP) technology.

---

## How It Works

### 1. **Intron Retention & Alternative Splicing Pipeline**
This pipeline detects and quantifies minor intron retention and alternative splicing events. It uses RNA sequencing data and consists of:

- **Scripts**: Includes BASH scripts and a Conda environment for integration with computing clusters.
- **Filtering Criteria**: Retained introns must pass thresholds for exon-intron boundary reads, splice site alignment, and intron coverage.
- **Output**: Mis-splicing index (MSI) and categorized splicing events (CAT1–CAT9) for each minor intron.

#### Data Sources
- RNA-seq data from patient-derived glioblastoma neurospheres and prostate cancer cell lines (C4-2, MCF-10A MYC-ER, etc.).

#### Future Enhancements
The pipeline can be adapted to analyze major introns with improved parallelization.

---

### 2. **Spatial Omics Pipeline**
This pipeline analyzes spatial transcriptomic data obtained using GeoMx DSP technology. It identifies differentially expressed genes and deconvolutes cell type composition in high vs. low U6atac expression regions.

#### Key Features
- **Data Preparation**: Filters low-quality samples and applies `DESeq2` and `Limma` for differential expression analysis.
- **Deconvolution**: Uses the `SpatialDecon` R package with prostate-specific profile matrices.

#### Outputs
- Gene expression profiles for regions of interest (ROIs) in the tumour and tumour microenvironment (TME).
- Retention and alternative splicing events count. 
---

### 3. **Minor Intron Database Scraper**
A Python tool automates querying the Minor Intron Database for efficient classification of genes of interest based on intron type.

---

## Installation (Basic)

1. Clone this repository:
   ```bash
   git clone https://github.com/XiaoyueLenax/U6atac_Landscape.git
   cd U6atac_Landscape


2. Set up the Conda environment for the intron retention pipeline:
    ```bash
    conda env create -f environment.yml
    conda activate environment
