# MPhil Project

This repository contains analysis scripts, processed data outputs, and supporting files for the MPhil thesis:  

**_Single-cell transcriptomic- and network-based identification of stage-specific modulators of brown adipocyte differentiation for anti-obesity drug repurposing_**.

## Overview
The project integrates **single-nucleus RNA-sequencing** analysis of human perivascular adipose tissue with **network-based drug prioritisation** to identify existing compounds that could promote brown adipocyte differentiation, a potential therapeutic strategy for obesity.

## Repository structure
- **`scripts/`** – R and Python scripts for snRNAseq data processing, pathway enrichment analysis, and network analysis.  
- **`results/`** – Figures, tables, and processed outputs generated during the analysis.  

## Data availability
Raw snRNA-seq data are publicly available from NCBI GEO:  
[GSE164528 – Defining the lineage of thermogenic perivascular adipose tissue](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164528)  

## Reproducibility
- All analyses were run in **R (v4.4.1)** and **Python (v3.9.12)**; see **`requirements.txt`** for required Python packages and versions.
- To reproduce results:
  1. Download raw data from GEO (link above).
  2. Follow the workflow in `scripts/` as described in thesis *Methods*.
  3. Output files will be saved in `results/`.
