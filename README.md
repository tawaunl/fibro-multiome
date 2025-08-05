# Fibrotic Memory of Fibroblasts: Single-Cell Multiome Analysis

## Overview

This repository contains scripts for analyzing single-cell ATAC-seq and RNA-seq multiome data from human lung fibroblasts to investigate fibrotic memory mechanisms. The analysis compares normal and IPF (Idiopathic Pulmonary Fibrosis) fibroblasts under different treatment conditions (No Treatment, OSM, TGF-β, IL-13) to understand epigenetic priming and transcriptional responses.

## Project Structure

```         
AnalysisScripts/
├── Analysis.Rmd              # Main multiome analysis workflow
├── BatchCorrection.R         # Harmony batch correction and cell scoring
├── CellLabels.R              # Cell type annotation and module scoring
├── DA_Analysis.R             # Differential accessibility analysis
├── DE_byTreatment.R          # Differential expression analysis
├── DifferentialAbundance.R   # Cell abundance analysis
├── Fig3FinalAnalysis.R       # Final analysis for manuscript figures
└── README.md                 # This file
```

## Environment Setup

### Required R Packages

#### Core Analysis

``` r
# Single-cell analysis
library(ArchR)           # Main multiome analysis framework
library(BSgenome.Hsapiens.UCSC.hg38)  # Human genome reference

# Data manipulation and visualization
library(grid)
library(gridExtra)
library(dsassembly)
library(ggplot2)
library(cowplot)
library(ComplexHeatmap)
library(pheatmap)

# Statistical analysis
library(edgeR)
library(VennDiagram)
```

#### Specialized Packages

``` r
# Motif analysis
library(chromVARmotifs)  # Install: devtools::install_github("GreenleafLab/chromVARmotifs")

# 4-way comparison plots
library(gg4way)

# Gene annotation
library(org.Hs.eg.db)

# Color schemes
library(RColorBrewer)
library(viridis)
library(wesanderson)

# Data processing
library(tidyverse)
library(dplyr)
library(ggrepel)
library(ggrastr)
```

#### Alternative Analysis (if using Seurat/Signac)

``` r
library(Seurat)
library(Signac)
```

### External Dependencies

1.  **MACS2** for peak calling:

    ``` bash
    conda create -n MACS python=3.8
    conda activate MACS
    conda install -c bioconda macs2
    ```

2.  **ArchR Genome Annotations**:

    ``` r
    addArchRGenome("hg38")  # Downloads automatically
    ```

## Input Data

The analysis expects **Cell Ranger ARC** output containing:

1.  Filtered feature barcode matrices

2.  ATAC-seq fragments file

**Sample Structure**:

Normal fibroblasts: 4 samples (No Treatment, OSM, TGF-β, IL-13)

IPF fibroblasts: 4 samples (No Treatment, OSM, TGF-β, IL-13)

## Analysis Workflow

### 1. Initial Setup and QC (`Analysis.Rmd`)

**Purpose**: Main analysis pipeline for multiome data processing

**Key Steps**: - Load and create ArchR project - Dimensionality reduction (LSI) for ATAC and RNA - Combined dimensionality reduction - UMAP embedding generation - Cluster identification - Peak-to-gene linkage analysis - Marker gene/peak identification

**Outputs**: - ArchR project saved to `Save-ProjMulti2/` - UMAP plots for different modalities - Marker gene heatmaps - Peak-to-gene linkage results

### 2. Batch Correction (`BatchCorrection.R`)

**Purpose**: Remove batch effects and perform cell type scoring

**Key Steps**: - Harmony batch correction by sample - Cell type scoring using predefined gene modules - Cell type annotation based on marker expression - Generate corrected UMAP embeddings

**Outputs**: - Harmony-corrected embeddings - Cell type annotations - Module score plots

### 3. Cell Type Annotation (`CellLabels.R`)

**Purpose**: Detailed fibroblast subtype identification

**Gene Modules**: - **Pan-Fibroblast**: DCN, DPT, FBLN1, FN1, PDGFRA, etc. - **Adventitial**: CD34, PI16, CADM3, CLEC3B, etc. - **Alveolar**: CD82, CES1, FGFR4, GPC3, etc. - **Myofibroblast**: CDH11, ASPN, CTHRC1, POSTN, etc. - **Smooth Muscle**: ACTA2, MYH11, TAGLN, etc. - **Pericyte**: COX4I2, CSPG4, RGS5, etc.

**Outputs**: - Cell type assignments - Module score heatmaps - UMAP plots colored by cell types

### 4. Differential Accessibility Analysis (`DA_Analysis.R`)

**Purpose**: Identify treatment-specific chromatin accessibility changes

**Comparisons**: - IL-13 vs No Treatment (IPF vs Normal) - OSM vs No Treatment (IPF vs Normal)\
- TGF-β vs No Treatment (IPF vs Normal) - IPF vs Normal (baseline)

**Outputs**: - Volcano plots for each comparison - 4-way comparison plots - Venn diagrams of overlapping peaks - TFBS enrichment analysis - Browser tracks for key regions

### 5. Differential Expression Analysis (`DE_byTreatment.R`)

**Purpose**: Identify treatment-induced gene expression changes

**Key Features**: - Pairwise comparisons for each treatment - 4-way analysis plots using `gg4way` - Identification of IPF-specific responses - Venn diagram analysis of overlapping genes

**Outputs**: - Expression volcano plots - Treatment-specific gene lists - Overlap analysis between treatments

### 6. Cell Abundance Analysis (`DifferentialAbundance.R`)

**Purpose**: Analyze changes in cell type proportions

**Methods**: - Cluster resolution optimization - Stacked bar plots by condition - EdgeR-based abundance testing - Cellularity percentage calculations

**Outputs**: - Cell proportion plots - Statistical tests for abundance changes - Resolution comparison plots

### 7. Manuscript Figure Generation (`Fig3FinalAnalysis.R`)

**Purpose**: Generate publication-ready figures

**Key Analyses**:

#### Panel 1: TF Motif Enrichment

-   Motif enrichment in IPF vs Normal
-   Heatmap of top enriched motifs
-   Focus on priming-related transcription factors

#### Panel 2: Priming Signature Analysis

#### **Gene Modules**: 

**Transcription Factors**: FOXC2, FOXP1, CEBPB, STAT3, etc. -

**Profibrotic Effectors**: TGFB1, IL6, PDGFRA, FN1, etc. -

**Contractility/Survival**: ACTA2, ITGA5, MMP2, etc. -

**Senescence/Stress**: CDKN1A, GDF15, SERPINE1, etc. -

**EMT/Matrix**: ZEB1, VIM, LOX, LOXL2, etc. -

**Cytokines**: CXCL12, CCL2, IL1B, etc. -

**Chromatin Remodeling**: HMGA1, EZH2, BRD4, etc.

#### Panel 3: IPF Treatment Response

-   Treatment-induced TF activation in IPF
-   chromVAR deviation scores
-   Expression vs motif activity correlation

#### Panel 4: Normal Fibroblast Response

-   Differential responses in normal fibroblasts
-   Comparison with IPF responses
-   TF priming index calculation

## Citation

If you use this analysis pipeline, please cite: - ArchR: Greenleaf et al. (Nature Genetics, 2021) - Additional methods as appropriate for your publication

## Contact

For questions about this analysis pipeline, please contact:

Data, Analyis, interpretation: Rakshin Kharwadkar [kharwadr\@gene.com](mailto:kharwadr@gene.com){.email} (First author)

Code: Tawaun Lucas (Bioinformatics support)
