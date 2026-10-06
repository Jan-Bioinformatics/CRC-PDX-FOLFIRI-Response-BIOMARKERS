# Transcriptomic Determinants of FOLFIRI / 5-FU Response in Colorectal Cancer PDX Models

R analysis code for my MSc dissertation (Bioinformatics and Computational Genomics, Queen's University Belfast, 2026). The project uses bulk RNA-seq data from colorectal cancer (CRC) liver-metastasis **patient-derived xenograft (PDX)** models to characterise transcriptional differences between responders and non-responders to FOLFIRI / 5-FU.

- **Responders** are **Partial Response (PR)** models.
- **Non-responders** are **Progressive Disease (PD)** models.

The analysis has two arms:

| Arm | Question | Design |
|---|---|---|
| **Baseline (cross-sectional)** | Do untreated tumours that will respond differ from those that won't? | Unpaired PR vs PD on baseline samples; response averaged over weeks 3–6 |
| **Longitudinal** | How does treatment change the transcriptome over time, and does that change differ between PR and PD? | The same PDX models sampled at Placebo, 24 h and 6 weeks of treatment |

---

## Repository contents

| Script | Analysis | Design formula |
|---|---|---|
| `Baseline_PR_vs_PD.R` | Baseline PR vs PD (with SD as context) | `~ Response` |
| `Long_overall_timepoint.R` | Treatment effect over time: 24h vs Placebo, 6wks vs Placebo, 6wks vs 24h | `~ Case_ID + Treatment` (paired) |
| `Long_Response_PR_vs_PD.R` | PR vs PD within each arm: Placebo, 24h, 6wks | `~ Response_3_6wk` (unpaired) |
| `Long_Interaction_Placebo_vs_6wks.R` | Does the Placebo → 6wks change differ between PR and PD? | `~ Response + Treatment + Response:Treatment` |

### Suggested run order

1. `Baseline_PR_vs_PD.R`. This runs on its own and has no dependency on the other scripts.
2. `Long_overall_timepoint.R`
3. `Long_Response_PR_vs_PD.R`
4. `Long_Interaction_Placebo_vs_6wks.R` (run last, as planned in the study design)

Each longitudinal script reads the same counts and metadata files, so scripts 2–4 do not depend on each other's outputs.

---

## Workflow at a glance

Every script follows the same layout:

```
Header           purpose, workflow, inputs, outputs
0. Setup         libraries, project paths
1. Config        every threshold / parameter in one place
2. Load data     counts + metadata, sample alignment
3. Resources     MSigDB gene sets, DoRothEA / CollecTRI regulons
4. Helpers       reusable, documented functions
5+ Steps         analysis in numbered steps
Last             sessionInfo()
```

### Methods used

| Step | Method | Package |
|---|---|---|
| Batch correction (baseline) | ComBat-seq on raw counts | `sva` |
| Differential expression | Negative-binomial GLM, Wald test, BH-adjusted | `DESeq2` (`glmGamPoi` for paired models) |
| LFC shrinkage (baseline) | apeglm | `apeglm` |
| Sample QC | VST + PCA, Cook's distance, library-size sensitivity re-fits | `DESeq2` |
| Global structure (baseline) | PERMANOVA, PERMDISP, NMDS | `vegan` |
| Pre-ranked GSEA | Wald statistic ranking; Hallmark, KEGG (legacy), GO (C5), C6 | `clusterProfiler`, `msigdbr` |
| Over-representation | Hypergeometric test of DEGs | `clusterProfiler::enricher` |
| Single-sample scores | ssGSEA + Wilcoxon / Kruskal-Wallis | `GSVA` |
| Dual GSEA consensus | Pre-ranked NES vs ssGSEA delta / difference-in-differences | custom |
| TF activity | Univariate linear model on DoRothEA (A/B/C) and CollecTRI regulons | `decoupleR`, `dorothea` |
| Gene-level validation | edgeR TMM log2CPM boxplots | `edgeR` |
| Subtype context | CMS / CRIS association, Fisher's exact test (descriptive only) | base R |

---

## Conventions

- **Direction:** a positive log2FC, NES, ssGSEA delta or TF score means **higher in the first-named group** (PR in response analyses; the later timepoint in treatment analyses). PD and Placebo are the reference levels.
- **Main thresholds:**
  - DEGs: `padj < 0.05` and `|log2FC| > 0.585` (1.5-fold); the interaction model uses `|log2FC| > 1`.
  - Pre-filter: genes with ≥ 10 total reads.
  - Pathways and TFs: BH-adjusted `p < 0.05` unless stated in the script header.
- **Seed:** `set.seed(123)` is used throughout for reproducibility.
- **IDs:**
  - `Case_ID` (e.g. `CRC0081`) identifies the patient / PDX model line.
  - The full sample barcode (e.g. `CRC0081LMX0B02202TUMFF0100`) identifies the individual specimen.

---

## Expected project structure

The scripts use relative paths from a project root (set with `setwd()` at the top of each script):

```
PDX_PROJECT/
├── 0_data/
│   ├── PDX_Baseline/
│   │   ├── counts/raw/        GEX_raw_counts_Human.csv, selected_metadata.csv
│   │   ├── clinical_data/     clinical_data.xlsx
│   │   └── annotations/       CMS.csv, CRIS.csv, YA_DTP_signatures.xlsx
│   └── PDX_Longitudinal/
│       └── counts/            CRC_graft_raw_counts_matrix_12022024.csv,
│                              Complete_GraftCRC_colData_19022024.csv
├── 2_pipelines/               fitted DESeq2 objects (.rds), full DE tables
└── 3_output/                  figures and result tables (created by the scripts)
```

To run on another machine, edit the `setwd()` line and the path constants in **Section 0** of each script. Nothing else needs changing.

---

## Requirements

- **R:** ≥ 4.3
- **Bioconductor:** `DESeq2`, `glmGamPoi`, `apeglm`, `sva`, `edgeR`, `clusterProfiler`, `enrichplot`, `GSVA`, `decoupleR`, `dorothea`, `ComplexHeatmap`
- **CRAN:** `tidyverse`, `msigdbr`, `ggrepel`, `ggpubr`, `patchwork`, `circlize`, `matrixStats`, `vegan`, `readxl`, `scales`, `tidygraph`, `ggraph`

The longitudinal scripts contain an `INSTALL_PACKAGES` switch for a first-time setup. Exact package versions are printed by `sessionInfo()` at the end of every script.

---

## Data availability

Raw counts and clinical metadata are **not included** in this repository. The PDX models and sequencing data come from the Candiolo mCRC biobank (Turin). Access is subject to the data owners' approval.

---

## Notes and limitations

- **Small group sizes:** PR samples are few in some comparisons, notably 24h. Results there should be read as exploratory.
- **Repeated samples:** at baseline, a few models contribute more than one sample, and these are treated as independent.
- **Outlier checks:** sensitivity re-fits are included in `Long_overall_timepoint.R` for one low-depth outlier sample and for the lowest-depth libraries.
- **HBA2:** this gene appears in the interaction model, and the script includes a check because it may reflect blood contamination.
- **Subtypes:** CMS / CRIS subtype associations are descriptive only.

---

## Author

**Dr Janmita Kaverimane Umesh**, MSc Bioinformatics and Computational Genomics, Queen's University Belfast

- GitHub: [Jan-Bioinformatics](https://github.com/Jan-Bioinformatics)
- LinkedIn: [dr-k-u-janmita](www.linkedin.com/in/dr-janmita-kaverimane-umesh-3393a431a)
