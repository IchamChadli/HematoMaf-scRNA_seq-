# HematoMaf-scRNA_seq

## Single-Cell RNA-seq Analysis of Hematopoietic Stem Cells in Mafb Knockout Mice


<p align="center">
  <img width="297" height="341" alt="image" src="https://github.com/user-attachments/assets/c003fee4-223f-4479-8bb4-ed1c65f1ee94" />
</p>

## Project Overview

This repository contains the single-cell RNA-seq (scRNA-seq) analysis pipeline developed
for the study of hematopoietic stem cell (HSC) biology in the context of Mafb knockout (KO)
mice. The primary objective is to characterize the transcriptional heterogeneity of
hematopoietic progenitor populations and to understand how Mafb deficiency and
M-CSF stimulation alter HSC commitment dynamics and myeloid differentiation.

### Biological Context

The analysis focuses on four experimental conditions:
- WT_NT — Wild-type, non-treated
- WT_MCSF — Wild-type, M-CSF stimulated (20µg)
- KO_NT — Mafb knockout, non-treated
- KO_MCSF — Mafb knockout, M-CSF stimulated (20µg)

Cells were sequenced using the 10x Genomics Chromium platform and aligned with
Cell Ranger. The dataset covers the full hematopoietic hierarchy, including HSCs,
MPPs, and committed progenitors (myeloid, erythroid, lymphoid lineages).

---

## Pipeline Overview

### 1. Data Loading and Quality Control
- Loading of 10x Genomics count matrices via Seurat Read10X
- Merging of multiple conditions into a single Seurat object
- QC filtering:
  - Removal of cells with 0% mitochondrial gene expression (likely empty droplets)
  - Removal of cells with low RNA counts and low detected features
  - Doublet detection via scDblFinder

### 2. Normalization and Dimensionality Reduction
- Log-normalization with NormalizeData
- Highly variable gene selection with FindVariableFeatures
- Scaling and PCA with ScaleData and RunPCA
- Cell cycle regression
- UMAP embedding with RunUMAP

### 3. Clustering
- Graph-based clustering with FindNeighbors and FindClusters
- Cluster annotation using marker genes and reference-based tools

### 4. Cell Type Annotation — HemaScribe and HemaScape

A key component of this pipeline is the use of HemaScribe and HemaScape
from the Rabadan Lab (https://github.com/RabadanLab/HemaScribe), two tools
specifically designed for hematopoietic single-cell data.

#### HemaScribe
HemaScribe performs automated cell type annotation of hematopoietic cells by
projecting cells onto a reference hematopoietic atlas. It assigns three levels
of annotation:
- broad.annot — broad cell type (HSC, myeloid, lymphoid, erythroid...)
- fine.annot — fine-grained HSPC subtypes (LT-HSC, ST-HSC, MPP2, MPP3, MPP4...)
- combined.annot — combined annotation merging broad and fine resolution

#### HemaScape
HemaScape computes a pseudotime trajectory and an Elastic Embedding (EE)
dimensionality reduction specifically designed for hematopoietic differentiation.
It uses Symphony-based batch correction and provides:
- pseudotime_pred — pseudotime score from 0 (stem) to 1 (mature)
- branch_pred — hematopoietic branch assignment
- density_cluster_pred — density-based cluster prediction
- branch_segment_clusters_pred — refined branch segment clusters
- EE reduction — Elastic Embedding coordinates for visualization

The EE reduction provides a biologically meaningful layout where the stem cell
compartment (HSC/ST-HSC) is positioned at one end and committed progenitors
radiate outward along their respective lineage branches.

### 5. Differential Expression Analysis
- FindAllMarkers for identification of cluster-specific marker genes
- FindMarkers per cluster and per condition for differential gene expression
  between conditions (e.g., WT_NT vs KO_NT) within each annotated cell type
- Results filtered by adjusted p-value and log2 fold change thresholds

### 6. Functional Enrichment Analysis
- GO Biological Process and Molecular Function enrichment via clusterProfiler
- Semantic similarity reduction via rrvgo with scatter plots and treemaps

### 7. Module Scoring
- AddModuleScore applied with bulk RNA-seq signatures to score single cells:
  - Apoptosis, Inflammation, Aging, Quiescence gene sets
  - Bulk-derived DEG signatures (up/down in KO_NT, up/down in MCSF)
- Visualization via density ridge plots split by condition

---

## Dependencies

All required R packages and their versions are listed in 00_generalDeps.R.

---

## Usage

Set the working directory and launch the full report by sourcing:
launch_reports_compilation.R

The pipeline generates an interactive HTML report containing all QC metrics,
clustering results, cell type annotations, differential expression tables,
and functional enrichment analyses.

---

## Citation

If you use HemaScribe or HemaScape in your work, please cite:
Rabadan Lab. HemaScribe: automated annotation of hematopoietic single-cell data.
https://github.com/RabadanLab/HemaScribe

---

## Author

Icham Chadli
Master's graduate in Health Biology, specialization in Biomarkers and AI
CIML — Dr. Toby Lawrence's Lab, Marseille, France
