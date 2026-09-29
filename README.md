# Longitudinal Transcriptomic Analysis of FOLFIRI Response in CRC PDX Models

R code for the longitudinal RNA-seq analysis in my MSc dissertation, *Characterising Transcriptomic Determinants of FOLFIRI/5-FU Response in Colorectal Cancer using Patient-Derived Xenograft Models* (MSc Bioinformatics and Computational Genomics, Queen's University Belfast, 2026).

## Overview

PDX models were sampled at three timepoints: **Placebo**, **24 hours** and **6 weeks** of FOLFIRI treatment. Response (PR vs PD) was classified from the 3–6 week average tumour volume change.

Three analyses are included:

1. **Treatment effect**: paired comparisons, matched by PDX model, ignoring response (24h vs Placebo, 6w vs Placebo, 6w vs 24h)
2. **PR vs PD within each arm**: unpaired comparisons at Placebo, 24h and 6w
3. **Treatment × Response interaction**: Placebo vs 6w (`~ Response + Treatment + Response:Treatment`)

Each analysis includes:
- Differential expression (DESeq2)
- Pathway enrichment: GSEA, ssGSEA and ORA with Hallmark, KEGG and GO gene sets
- Transcription factor activity (DoRothEA + decoupleR ULM)
- Visualisation: volcano plots, GSEA bar plots, TF activity plots, waterfall plots

## Repository structure

```
├── 01_treatment_effect/     # Paired timepoint comparisons
├── 02_PR_vs_PD/             # Placebo, 24h and 6w response comparisons
├── 03_interaction/          # Treatment × Response model
└── README.md
```

## Requirements

- **R:** 4.6.1
- **Differential expression:** DESeq2, glmGamPoi
- **Enrichment:** clusterProfiler, msigdbr, GSVA
- **TF activity:** dorothea, decoupleR
- **Plotting:** ggplot2, ggrepel, patchwork, ComplexHeatmap
- **Data handling:** tidyverse

## Usage

1. Update the paths at the top of each script (`clinicals`, `pipelines`, `output`).
2. Run the scripts in the folder order above.
3. Results tables (CSV) and figures (PNG) are written to the `output` directory.

## Data availability

Raw count data and clinical metadata are not included. The PDX models come from the XENTURION biobank (Candiolo Cancer Institute, Italy; Leto et al., 2024, *Nature Communications*), and data access is subject to the original data providers.

## Author

Janmita Kaverimane Umesh, MSc Bioinformatics and Computational Genomics, Queen's University Belfast
Supervisor: Prof. Simon McDade
