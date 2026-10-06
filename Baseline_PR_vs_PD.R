################################################################################
#
#  Script      : Baseline_PR_vs_PD.R
#  Project     : Transcriptomic profiling of FOLFIRI / 5-FU response in
#                colorectal cancer (CRC) patient-derived xenograft (PDX) models
#  Analysis    : Baseline (untreated LMX) transcriptome — Partial Response (PR)
#                vs Progressive Disease (PD), response averaged over weeks 3-6
#  Author      : Janmita Kaverimane Umesh (MSc Bioinformatics, QUB)
#
#  Workflow
#  ---------------------------------------------------------------------------
#    1.  Load counts + metadata; ComBat-seq batch correction (all samples)
#    2.  Build cohorts: PR/PD (main) and PR/SD/PD (context)
#    3.  SD placement — 3-group PCA, PERMANOVA, centroid distances, PERMDISP
#    4.  Cohort overview — sample counts, waterfalls, batch and CMS/CRIS tables
#    5.  DESeq2 (~ Response, PD = reference)
#    6.  DEG summaries — threshold counts, MA plot, volcano, DEG tables
#    7.  Heatmaps of DEGs (z-scored) + DEG signature score by CMS
#    8.  Global structure — PERMANOVA (top variable genes), NMDS, PERMDISP
#    9.  Gene-level validation — edgeR TMM log2CPM boxplots
#    10. Pre-ranked GSEA — Hallmark, KEGG, C6 oncogenic, C5 GO
#    11. TF activity — DoRothEA (A/B/C) and CollecTRI, ULM, TF-DEG networks
#    12. Dual GSEA consensus — pre-ranked GSEA NES vs ssGSEA delta
#    13. Drug-tolerant persister (DTP) signatures — ssGSEA UP / DOWN
#    14. TEAD1 / TEAD4 sample-level activity
#
#  Direction convention: positive log2FC / NES / delta / TF score = higher in PR.
#
#  Inputs
#    - 0_data/PDX_Baseline/counts/raw/GEX_raw_counts_Human.csv : raw counts
#    - 0_data/PDX_Baseline/counts/raw/selected_metadata.csv    : sample metadata
#    - 0_data/PDX_Baseline/annotations/CMS.csv, CRIS.csv       : subtype calls
#    - 0_data/PDX_Baseline/annotations/YA_DTP_signatures.xlsx  : DTP signatures
#    - 0_data/PDX_Baseline/clinical_data/clinical_data.xlsx    : tumour response
#    - CollecTRI network exported to CSV (CT_NET_FILE)
#
#  Outputs
#    - 0_data/.../counts/raw/GEX_corrected.csv : ComBat-seq corrected counts
#    - 2_pipelines/PDX_Baseline/DESeq2/        : DE table, normalised counts
#    - 3_output/PDX_Baseline/DESeq2/3-6av_PR_V_PD/Final/ : figures and tables
#
#  Notes
#    - All thresholds and parameters are unchanged from the original script.
#
################################################################################


#===============================================================================
# 0. SESSION SETUP
#===============================================================================

rm(list = ls(all.names = TRUE))     # start from a clean environment
gc()                                # release memory from previous sessions

# ---- 0.1 Libraries ----
# Core DE and batch correction
library(DESeq2)          # differential expression (apeglm must be installed for lfcShrink)
library(sva)             # ComBat_seq batch correction for counts
library(edgeR)           # TMM-normalised CPM for gene-level boxplots
# Enrichment and regulatory analysis
library(clusterProfiler) # pre-ranked GSEA
library(enrichplot)      # gseaplot2 running-score curves
library(msigdbr)         # MSigDB gene-set collections
library(GSVA)            # single-sample GSEA (ssGSEA)
library(dorothea)        # TF -> target regulons
library(decoupleR)       # TF activity inference (ULM)
# Multivariate statistics
library(vegan)           # adonis2 (PERMANOVA), betadisper (PERMDISP), metaMDS (NMDS)
library(matrixStats)     # fast row-wise SD / variance
# Data handling
library(tidyverse)       # dplyr, tidyr, tibble, stringr, ggplot2
library(readxl)
# Plotting
library(ggrepel)         # non-overlapping labels
library(ggpubr)          # stat_compare_means on boxplots
library(scales)          # axis label formatting
library(ComplexHeatmap)  # heatmaps (also provides a pheatmap()-compatible wrapper)
library(circlize)        # colorRamp2 for heatmap colour scales
library(tidygraph)       # graph objects for TF-target networks
library(ggraph)          # network plotting

# ---- 0.2 Project paths ----
setwd("C:/Users/Admin/OneDrive/Desktop/R studio Files/CRC project/PDX_PROJECT")

COUNTS_DIR    <- "./0_data/PDX_Baseline/counts/raw/"
CLINICAL_DIR  <- "./0_data/PDX_Baseline/clinical_data/"
ANNOT_DIR     <- "./0_data/PDX_Baseline/annotations/"
PIPELINES_DIR <- "./2_pipelines/PDX_Baseline/DESeq2/"
OUTPUT_DIR    <- "./3_output/PDX_Baseline/DESeq2/3-6av_PR_V_PD/Final/"

GEX_FILE       <- paste0(COUNTS_DIR, "GEX_raw_counts_Human.csv")       # 613 samples
METADATA_FILE  <- paste0(COUNTS_DIR, "selected_metadata.csv")
CORRECTED_FILE <- paste0(COUNTS_DIR, "GEX_corrected.csv")
CLINICAL_FILE  <- paste0(CLINICAL_DIR, "clinical_data.xlsx")
DTP_FILE       <- paste0(ANNOT_DIR, "YA_DTP_signatures.xlsx")
CT_NET_FILE    <- "C:/Users/Admin/Downloads/ct_net.csv"                # CollecTRI network (exported)

dir.create(PIPELINES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_DIR,    recursive = TRUE, showWarnings = FALSE)


#===============================================================================
# 1. ANALYSIS CONFIGURATION
#===============================================================================

PREFIX          <- "3-6w_PR.PD_batch"   # prefix for most output files
WEEK            <- "3-6"
SEED            <- 123

# ---- 1.1 Cohort definition ----
RESP_THRESHOLD  <- 35                   # % volume change reference line (waterfalls)
EXCLUDED_MODELS <- "CRC0464"            # excluded from the response table
BASELINE_TYPES  <- c("LMX_BASALE", "LMX_BASALE.1")   # untreated liver-metastasis xenografts
RUN_COMBAT_SEQ  <- TRUE                 # FALSE = reuse the saved GEX_corrected.csv

# ---- 1.2 Differential expression ----
MIN_TOTAL_CNT   <- 10                   # pre-filter: keep genes with >= 10 reads in total
PADJ_CUTOFF     <- 0.05                 # BH-adjusted significance threshold
LFC_CUTOFF      <- 0.585                # |log2FC| > 0.585  ==  fold change > 1.5
TOP_N_HEATMAP   <- 50                   # DEGs in the top-gene heatmaps

# ---- 1.3 Multivariate structure ----
N_TOP_VAR_GENES <- 2000                 # most variable genes for PERMANOVA / NMDS
N_PERM          <- 999                  # permutations (PERMANOVA, PERMDISP)

# ---- 1.4 Enrichment ----
CONSENSUS_PADJ  <- 0.1                  # dual-GSEA consensus threshold (both methods)
TEAD_TARGETS    <- c("FRMD6", "GLIS3")  # TEAD target genes shown in CPM boxplots

# ---- 1.5 Colours ----
RESP_COLS       <- c("PD" = "#F8766D", "PR" = "#4DBBD5")                     # waterfalls, boxplots
RESP_COLS_SD    <- c("PD" = "#F8766D", "SD" = "#FDFD96", "PR" = "#4DBBD5")
RESP_POINT_COLS <- c("PD" = "#c9463a", "PR" = "#2a8fa0")                     # jittered points
HEAT_COL_FUN    <- colorRamp2(c(-2, 0, 2), c("#BA55D3", "#FFFFF0", "#2E8B57"))  # z-score scale


#===============================================================================
# 2. LOAD DATA AND BATCH CORRECTION
#===============================================================================

# ---- 2.1 Raw counts (genes x samples) ----
GEX <- read.csv(GEX_FILE, row.names = 1, check.names = FALSE)
cat("Successfully loaded GEX with", nrow(GEX), "genes.\n")

# ---- 2.2 Sample metadata (sequencing batch, sample type) ----
colData <- read.csv(METADATA_FILE, row.names = 1)
colData$batch    <- as.factor(colData$batch)
colData$batch_id <- as.factor(colData$batch)

# ---- 2.3 Pre-computed CMS and CRIS subtype calls ----
CMS  <- read.csv(paste0(ANNOT_DIR, "CMS.csv"),  row.names = 1)
CRIS <- read.csv(paste0(ANNOT_DIR, "CRIS.csv"), row.names = 1)
colData$CMS  <- CMS[rownames(colData), 1]
colData$CRIS <- CRIS[rownames(colData), 1]

# ---- 2.4 ComBat-seq across ALL samples (negative-binomial; output stays counts) ----
# group = NULL: no biological covariate is protected, so correction is
# response-blind (response labels exist for only a subset of samples).
if (RUN_COMBAT_SEQ) {
  full_samples  <- intersect(colnames(GEX), rownames(colData))
  GEX_mat       <- as.matrix(GEX[, full_samples])
  batches       <- as.factor(colData[full_samples, "batch"])
  GEX_corrected <- ComBat_seq(counts = GEX_mat, batch = batches, group = NULL)
  write.csv(GEX_corrected, file = CORRECTED_FILE, row.names = TRUE)
}

# Re-read from disk so both branches use identical input
GEX_corrected <- read.csv(CORRECTED_FILE, row.names = 1)


#===============================================================================
# 3. RESPONSE DATA AND COHORTS
#===============================================================================

# Read the response table, drop missing values and excluded models, and keep
# the requested response classes (PD is always the first / reference level).
read_response <- function(levels) {
  read_excel(CLINICAL_FILE, sheet = 1) %>%
    dplyr::select(IDs = 1, perc = 13, Response = 14) %>%     # model ID, % volume change, class
    dplyr::filter(
      !is.na(IDs) & IDs != "NA" & IDs != "N/A" & !(IDs %in% EXCLUDED_MODELS),
      !is.na(perc) & perc != "NA" & perc != "N/A",
      !is.na(Response) & Response != "NA" & Response != "N/A"
    ) %>%
    dplyr::filter(Response %in% levels) %>%
    dplyr::mutate(perc     = as.numeric(perc),
                  Response = factor(Response, levels = levels))
}

# Match response labels to baseline RNA-seq samples. The first 7 characters
# of a sample name are the PDX model ID (e.g. "CRC0081"). Only untreated LMX
# baseline samples are kept; models can contribute >1 baseline sample.
build_cohort <- function(resp_tbl, counts, coldata, label) {

  model_ids     <- unique(substr(colnames(counts), 1, 7))
  selected_resp <- resp_tbl %>% dplyr::filter(IDs %in% model_ids)
  cat("\n[", label, "] models with response:", nrow(resp_tbl),
      "| with counts:", nrow(selected_resp),
      "| without counts:", nrow(resp_tbl) - nrow(selected_resp), "\n")

  # ---- All RNA-seq columns for those models, with sample type ----
  matching_cols  <- colnames(counts)[substr(colnames(counts), 1, 7) %in% selected_resp$IDs]
  sample_mapping <- data.frame(LongName = matching_cols,
                               IDs      = substr(matching_cols, 1, 7),
                               Type     = coldata[matching_cols, "type"],
                               stringsAsFactors = FALSE)

  # ---- Keep baseline (untreated) samples only ----
  keep_cols <- sample_mapping %>% dplyr::filter(Type %in% BASELINE_TYPES) %>% dplyr::pull(LongName)

  # ---- Subset counts and metadata; attach response ----
  cohort_counts <- counts[, keep_cols]
  cohort_meta   <- coldata %>% dplyr::filter(rownames(coldata) %in% keep_cols)
  cohort_meta   <- cohort_meta %>% dplyr::mutate(IDs = substr(rownames(cohort_meta), 1, 7))
  sample_ids    <- rownames(cohort_meta)
  cohort_meta   <- cohort_meta %>% dplyr::left_join(selected_resp, by = "IDs")
  rownames(cohort_meta) <- sample_ids                          # left_join drops row names

  # ---- Align metadata rows to count columns (required by DESeq2) ----
  if (!setequal(rownames(cohort_meta), colnames(cohort_counts))) {
    stop("Mismatch between colData and GEX column names!")
  }
  cohort_meta <- cohort_meta[colnames(cohort_counts), ]
  if (!all(rownames(cohort_meta) == colnames(cohort_counts))) stop("Reordering failed! Investigate further.")
  cat("Reordering successful. Data aligned.\n")
  cat("Filtered GEX dimensions:", dim(cohort_counts), "\n")
  cat("Filtered colData dimensions:", dim(cohort_meta), "\n")

  # ---- Models with more than one baseline sample ----
  dup_ids <- unique(cohort_meta$IDs[duplicated(cohort_meta$IDs)])
  cat("Duplicate Patient IDs found:\n")
  print(table(cohort_meta$IDs[cohort_meta$IDs %in% dup_ids]))

  list(counts = cohort_counts, meta = cohort_meta)
}

# ---- 3.1 Main cohort: PR vs PD ----
resp         <- read_response(c("PD", "PR"))
cohort_prpd  <- build_cohort(resp, GEX_corrected, colData, "PR/PD")
filtered_GEX   <- cohort_prpd$counts
merged_colData <- cohort_prpd$meta
merged_colData$Response <- factor(merged_colData$Response, levels = c("PD", "PR"))   # PD = reference

# ---- 3.2 Context cohort: PR vs SD vs PD ----
resp2        <- read_response(c("PD", "SD", "PR"))
cohort_sd    <- build_cohort(resp2, GEX_corrected, colData, "PR/SD/PD")
filtered_GEX2   <- cohort_sd$counts
merged_colData2 <- cohort_sd$meta
merged_colData2$Response <- factor(merged_colData2$Response, levels = c("PD", "SD", "PR"))

# ---- 3.3 Group sizes (used in plot subtitles) ----
N_PR <- sum(merged_colData$Response == "PR")
N_PD <- sum(merged_colData$Response == "PD")
N_SD <- sum(merged_colData2$Response == "SD")


#===============================================================================
# 4. REFERENCE RESOURCES
#===============================================================================

# ---- 4.1 MSigDB gene sets (TERM2GENE format: gs_name, gene_symbol) ----
h_t2g    <- msigdbr(species = "Homo sapiens", collection = "H")  %>% dplyr::select(gs_name, gene_symbol)
kegg_t2g <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY") %>%
  dplyr::select(gs_name, gene_symbol)
c5_t2g   <- msigdbr(species = "Homo sapiens", collection = "C5") %>% dplyr::select(gs_name, gene_symbol)
c6_t2g   <- msigdbr(species = "Homo sapiens", collection = "C6") %>% dplyr::select(gs_name, gene_symbol)

# ---- 4.2 DoRothEA regulons (A/B/C confidence) ----
data("dorothea_hs", package = "dorothea")
dorothea_net <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B", "C")) %>%
  dplyr::rename(source = tf) %>%                   # decoupleR expects 'source'
  dplyr::mutate(weight = as.numeric(mor))          # mode of regulation: +1 / -1
cat("DoRothEA rows (A/B/C):", nrow(dorothea_net), "\n")

# ---- 4.3 CollecTRI regulons (pre-exported CSV: source, target, mor) ----
ct_net <- read.csv(CT_NET_FILE, stringsAsFactors = FALSE)
stopifnot(all(c("source", "target", "mor") %in% colnames(ct_net)))
ct_net$mor <- as.numeric(ct_net$mor)
cat("CollecTRI rows:", nrow(ct_net), "\n")

# ---- 4.4 Drug-tolerant persister (DTP) signatures (HCT116 & SW, day 14 vs day 0) ----
dtp_up_genes   <- read_excel(DTP_FILE, sheet = "DTP_Up_HCT_and_SW")$name
dtp_down_genes <- read_excel(DTP_FILE, sheet = "DTP_down_HCT_and_SW")$name
dtp_up_genes   <- dtp_up_genes[!is.na(dtp_up_genes) & dtp_up_genes != ""]
dtp_down_genes <- dtp_down_genes[!is.na(dtp_down_genes) & dtp_down_genes != ""]

dtp_t2g <- data.frame(
  gs_name = c(rep("DTP_Up_Signature",   length(dtp_up_genes)),
              rep("DTP_Down_Signature", length(dtp_down_genes))),
  gene    = c(dtp_up_genes, dtp_down_genes)
)
print(table(dtp_t2g$gs_name))


#===============================================================================
# 5. HELPER FUNCTIONS
#===============================================================================

# -----------------------------------------------------------------------------
# 5.1 Cohort plots
# -----------------------------------------------------------------------------

# "PR(n = 15) / PD(n = 39)" style count string, in the order given
count_string <- function(meta, levels) {
  paste(sprintf("%s(n = %d)", levels, sapply(levels, function(l) sum(meta$Response == l))),
        collapse = " / ")
}

# One bar per sample, ordered by % tumour volume change
plot_waterfall <- function(meta, palette, title, xlab, count_levels, file) {
  plot_data <- meta %>%
    dplyr::mutate(SampleID = rownames(.)) %>%
    dplyr::arrange(dplyr::desc(perc)) %>%
    dplyr::mutate(SampleID = factor(SampleID, levels = SampleID))

  p <- ggplot(plot_data, aes(x = SampleID, y = perc, fill = Response)) +
    geom_col(width = 1, color = "white", linewidth = 0.05) +
    geom_hline(yintercept = RESP_THRESHOLD, linetype = "dashed", color = "grey", linewidth = 0.5) +
    scale_fill_manual(values = palette) +
    labs(title    = title,
         subtitle = paste("Threshold:", RESP_THRESHOLD, "% | n =", nrow(plot_data), "|",
                          count_string(meta, count_levels)),
         x = xlab, y = "% Tumour Volume Change", fill = "Status") +
    theme_classic() +
    theme(axis.text.x     = element_blank(),
          axis.ticks.x    = element_blank(),
          legend.position = "top",
          plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
          plot.subtitle   = element_text(size = 10, hjust = 0.5),
          axis.title      = element_text(size = 10))

  print(p)
  ggsave(filename = file, plot = p, width = 8, height = 5, dpi = 300)
  invisible(p)
}

# Bar chart of samples per response group
plot_sample_counts <- function(dds_obj, recode_map, fill_cols, title, file) {
  counts_df <- as.data.frame(colData(dds_obj)) %>%
    dplyr::count(Response) %>%
    dplyr::mutate(Response = dplyr::recode(Response, !!!recode_map))

  p <- ggplot(counts_df, aes(x = Response, y = n, fill = Response)) +
    geom_col(width = 0.45, colour = "white") +
    geom_text(aes(label = n), vjust = -0.5, size = 2.8) +
    scale_fill_manual(values = fill_cols) +
    labs(title = title, x = NULL, y = "Number of PDX samples") +
    theme_classic() +
    theme(legend.position = "none", plot.title = element_text(size = 14, hjust = 0.5))
  ggsave(file, p, width = 6, height = 3.2, dpi = 600)
  invisible(p)
}


# -----------------------------------------------------------------------------
# 5.2 Heatmap utilities
# -----------------------------------------------------------------------------

# Row-wise z-score (mean 0, SD 1 per gene); zero-variance genes set to 0
row_zscore <- function(mat) {
  mat   <- as.matrix(mat)
  means <- rowMeans(mat, na.rm = TRUE)
  sds   <- matrixStats::rowSds(mat, na.rm = TRUE)
  sds[sds == 0] <- NA
  z <- (mat - means) / sds
  z[is.na(z)] <- 0
  z
}

# Harmonise missing / unclassified subtype labels
clean_unclassified <- function(x) {
  x <- as.character(x)
  ifelse(is.na(x) | x == "" | x == "UNC", "Unclassified", x)
}


# -----------------------------------------------------------------------------
# 5.3 Gene-level boxplots (edgeR log2CPM)
# -----------------------------------------------------------------------------

# Faceted PR vs PD boxplots for a gene list, annotated with DESeq2 padj stars
# and log2FC. Facets keep the order of 'genes'.
plot_cpm_boxplots <- function(genes, log_cpm, sig_df, meta, title, subtitle, file) {

  plot_df <- as.data.frame(t(log_cpm[genes, , drop = FALSE])) %>%
    rownames_to_column("SampleID") %>%
    dplyr::left_join(meta %>% rownames_to_column("SampleID") %>% dplyr::select(SampleID, Response),
                     by = "SampleID") %>%
    pivot_longer(cols = -c(SampleID, Response), names_to = "Gene", values_to = "log2CPM") %>%
    dplyr::mutate(Gene     = factor(Gene, levels = genes),
                  Response = factor(Response, levels = c("PD", "PR")))

  # DESeq2 statistics for facet labels
  padj_labels <- sig_df %>%
    dplyr::filter(rownames(.) %in% genes) %>%
    dplyr::mutate(Gene       = rownames(.),
                  padj_label = dplyr::case_when(padj < 0.001 ~ "***",
                                                padj < 0.01  ~ "**",
                                                padj < 0.05  ~ "*",
                                                TRUE         ~ "ns"),
                  padj_text  = paste0("padj = ", formatC(padj, format = "e", digits = 2)),
                  fc_label   = paste0("log2FC = ", round(log2FoldChange, 2))) %>%
    dplyr::select(Gene, padj_label, padj_text, fc_label)

  plot_df <- plot_df %>%
    dplyr::left_join(padj_labels, by = "Gene") %>%
    dplyr::mutate(Gene_label = paste0(Gene, "\n(", padj_text, ")"))

  label_order <- plot_df %>% dplyr::distinct(Gene, Gene_label) %>%
    dplyr::arrange(match(Gene, genes)) %>% dplyr::pull(Gene_label)
  plot_df$Gene_label <- factor(plot_df$Gene_label, levels = label_order)

  p <- ggplot(plot_df, aes(x = Response, y = log2CPM, fill = Response)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.75, colour = "grey30", linewidth = 0.4) +
    geom_jitter(aes(colour = Response), width = 0.15, size = 1.5, alpha = 0.8, show.legend = FALSE) +
    geom_text(data = plot_df %>% dplyr::distinct(Gene_label, padj_label),
              aes(x = 1.5, y = Inf, label = padj_label),
              inherit.aes = FALSE, vjust = 1.5, size = 5, fontface = "bold") +
    geom_text(data = plot_df %>% dplyr::distinct(Gene_label, fc_label),
              aes(x = 1.5, y = Inf, label = fc_label),
              inherit.aes = FALSE, vjust = 3.5, size = 3, colour = "black") +
    scale_fill_manual(values = RESP_COLS,
                      labels = c("PD" = "Progressive Disease", "PR" = "Partial Response")) +
    scale_colour_manual(values = RESP_POINT_COLS) +
    facet_wrap(~ Gene_label, scales = "free_y", ncol = 3) +
    labs(title = title, subtitle = subtitle, x = NULL, y = "Log2 CPM", fill = "Response") +
    theme_classic(base_size = 11) +
    theme(plot.title       = element_text(face = "bold", hjust = 0.5, size = 13),
          plot.subtitle    = element_text(hjust = 0.5, colour = "grey40", size = 9),
          strip.text       = element_text(face = "bold", size = 9),
          strip.background = element_rect(fill = "grey95", colour = NA),
          legend.position  = "bottom",
          axis.text.x      = element_text(size = 9),
          panel.spacing    = unit(0.8, "lines"))

  print(p)
  ggsave(filename = file, plot = p, width = 10, height = 10, dpi = 500)
  invisible(p)
}


# -----------------------------------------------------------------------------
# 5.4 Pre-ranked GSEA (Hallmark / KEGG / C6)
# -----------------------------------------------------------------------------

# Shared theme for the single-colour GSEA bar plots
gsea_bar_theme <- function(title_size, x_title_size = 14, y_text_size = 8) {
  theme(axis.line        = element_line(colour = "black"),
        plot.title       = element_text(hjust = 0.5, size = title_size, face = "bold"),
        axis.title.y     = element_text(size = 14, colour = "black", face = "bold"),
        axis.text.y      = element_text(size = y_text_size, colour = "black"),
        axis.title.x     = element_text(size = x_title_size, colour = "black", face = "bold"),
        axis.text.x      = element_text(size = 12, colour = "black"),
        panel.background = element_blank(), panel.grid.minor = element_blank(),
        panel.border     = element_blank(), strip.background = element_blank())
}

# Run GSEA and produce, for one collection:
#   *_GSEA_RESULTS.csv            all tested sets
#   *_GSEA_barplot.png            all sets, significant ones highlighted
#   *_GSEA_Pathways_Sorted.csv    significant sets, NES high -> low
#   significant-only bar plot
#   running-score curve for the top-NES significant set
run_gsea_collection <- function(gene_list, t2g, cfg) {

  gsea_obj <- GSEA(gene_list, exponent = 1, minGSSize = 10, maxGSSize = 10000,
                   pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = t2g,
                   verbose = TRUE, seed = TRUE)
  df <- as.data.frame(gsea_obj)
  write.csv(df, file = paste0(OUTPUT_DIR, PREFIX, "_", cfg$tag, "_GSEA_RESULTS.csv"))

  # ---- All pathways, significant (p.adjust < 0.05) in green ----
  png(file = paste0(OUTPUT_DIR, PREFIX, "_", cfg$tag, "_GSEA_barplot.png"),
      width = cfg$all_size[1], height = cfg$all_size[2], units = "in", res = 500)
  print(
    ggplot(df, aes(reorder(Description, NES), NES)) +
      geom_col(aes(fill = p.adjust < PADJ_CUTOFF)) +
      coord_flip() +
      scale_fill_manual(values = c("#EBF1EC", "#29BF4E")) +
      labs(x = "Pathway", y = "Normalized Enrichment Score", title = cfg$all_title) +
      gsea_bar_theme(title_size = 16, x_title_size = cfg$all_x_title_size, y_text_size = 6) +
      theme(axis.line = element_line())
  )
  dev.off()

  # ---- Significant pathways sorted by NES ----
  sig_df    <- df[df$p.adjust < PADJ_CUTOFF, ]
  sorted_df <- sig_df[order(-sig_df$NES), ]
  write.csv(sorted_df, file = paste0(OUTPUT_DIR, PREFIX, "_", cfg$tag, "_GSEA_Pathways_Sorted.csv"))

  if (nrow(sig_df) > 0) {
    png(file = paste0(OUTPUT_DIR, PREFIX, "_", cfg$tag, cfg$sig_suffix),
        width = 11, height = cfg$sig_height, units = "in", res = 500)
    print(
      ggplot(sig_df, aes(reorder(Description, NES), NES)) +
        geom_col(fill = "#29BF4E", width = 0.7) +
        coord_flip() +
        labs(x = "Pathway", y = "Normalized Enrichment Score", title = cfg$sig_title) +
        gsea_bar_theme(title_size = cfg$sig_title_size)
    )
    dev.off()

    # ---- Running-score curve for the top (highest NES) significant set ----
    top_id <- sorted_df$ID[1]
    png(file = paste0(OUTPUT_DIR, PREFIX, "_", cfg$curve_tag, "_", gsub("[[:punct:]]", "_", top_id), ".png"),
        width = 10, height = 8, units = "in", res = 500)
    print(gseaplot2(gsea_obj, top_id, title = paste0(cfg$curve_title, top_id),
                    color = "green", base_size = 15, subplots = 1:3))
    dev.off()
  } else {
    message("No significant ", cfg$tag, " pathways found (p.adjust < 0.05).")
  }

  list(obj = gsea_obj, all = df, sig = sorted_df)
}


# -----------------------------------------------------------------------------
# 5.5 TF activity pipeline (DoRothEA or CollecTRI)
# -----------------------------------------------------------------------------

# For one regulon network:
#   ULM on the DESeq2 Wald statistic -> BH padj -> bar plot (top 10 each way)
#   targets of significant TFs present in the data
#   per-TF grouped target lists, TF-DEG overlap tables, top-10 TF network
# A 'confidence' column (DoRothEA only) is carried through where present.
run_tf_pipeline <- function(net, res, sig_df, files, bar_subtitle, net_subtitle) {

  # ---- Gene-level statistic as a 1-column matrix ----
  stat_vec <- res$stat
  names(stat_vec) <- rownames(res)
  stat_vec <- stat_vec[!is.na(stat_vec)]
  mat <- matrix(stat_vec, ncol = 1, dimnames = list(names(stat_vec), PREFIX))
  cat("Genes in matrix:", nrow(mat), "\n")

  # ---- ULM + BH correction across all TFs tested ----
  acts <- decoupleR::run_ulm(mat = mat, net = net, .source = "source",
                             .target = "target", .mor = "mor", minsize = 5) %>%
    dplyr::mutate(padj = p.adjust(p_value, method = "BH"))
  write.csv(acts, paste0(OUTPUT_DIR, PREFIX, files$acts), row.names = FALSE)
  cat("TFs tested:", nrow(acts), "| padj < 0.05:", sum(acts$padj < PADJ_CUTOFF), "\n")

  # ---- Bar plot: top 10 most positive and top 10 most negative ----
  tf_sig <- acts %>% dplyr::filter(padj < PADJ_CUTOFF) %>% dplyr::arrange(dplyr::desc(score))
  n_show <- min(10, nrow(tf_sig))
  tf_plot <- dplyr::bind_rows(head(tf_sig, n_show), tail(tf_sig, n_show)) %>%
    dplyr::distinct(source, .keep_all = TRUE) %>%
    dplyr::mutate(direction = ifelse(score > 0, "Active in PR", "Active in PD"),
                  source    = factor(source, levels = source[order(score)]))

  p_bar <- ggplot(tf_plot, aes(x = source, y = score, fill = direction)) +
    geom_col(width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c("Active in PR" = "#4DBBD5", "Active in PD" = "#F8766D")) +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey30") +
    labs(title = paste("TF activity (ULM) —", PREFIX), subtitle = bar_subtitle,
         x = "Transcription factor", y = "Activity score (ULM)", fill = NULL) +
    theme_classic(base_size = 11) +
    theme(plot.title      = element_text(face = "bold", hjust = 0.5, size = 16),
          plot.subtitle   = element_text(hjust = 0.5, colour = "grey40", size = 14),
          legend.position = "top",
          axis.text.x     = element_text(size = 12, color = "black"),
          axis.text.y     = element_text(size = 12, color = "black", face = "bold"),
          axis.title.x    = element_text(face = "bold", size = 14),
          axis.title.y    = element_text(face = "bold", size = 14),
          legend.text     = element_text(size = 12))
  print(p_bar)
  ggsave(paste0(OUTPUT_DIR, PREFIX, files$bar), p_bar, width = 11, height = 10, dpi = 500)

  # ---- Targets of significant TFs that are present in the data ----
  sig_targets <- net %>%
    dplyr::filter(source %in% tf_sig$source, target %in% rownames(mat)) %>%
    dplyr::select(TranscriptionFactor = source, TargetGene = target, ModeOfRegulation = mor,
                  any_of(c(Confidence = "confidence"))) %>%
    dplyr::arrange(TranscriptionFactor, dplyr::desc(ModeOfRegulation))
  cat("Network interactions:", nrow(net), "| targets of significant TFs in data:", nrow(sig_targets), "\n")
  write.csv(sig_targets, paste0(OUTPUT_DIR, PREFIX, files$targets), row.names = FALSE)

  # ---- Grouped target lists per TF and mode of regulation ----
  tf_stats <- acts %>%
    dplyr::select(TranscriptionFactor = source, TF_Activity_Score = score,
                  TF_p_value = p_value, TF_padj = padj)
  group_cols <- c("TranscriptionFactor", "TF_Activity_Score", "TF_p_value", "TF_padj", "ModeOfRegulation")

  grouped <- sig_targets %>%
    dplyr::left_join(tf_stats, by = "TranscriptionFactor") %>%
    dplyr::group_by(dplyr::across(all_of(group_cols))) %>%
    dplyr::summarise(dplyr::across(any_of("Confidence"), ~ paste(unique(.x), collapse = ",")),
                     TargetGenes = paste(TargetGene, collapse = ", "),
                     GeneCount   = dplyr::n(),
                     .groups     = "drop") %>%
    dplyr::arrange(TF_padj, TranscriptionFactor, dplyr::desc(ModeOfRegulation))
  write.csv(grouped, paste0(OUTPUT_DIR, PREFIX, files$grouped), row.names = FALSE)

  # ---- How many DEGs are covered by the network at all? ----
  deg_genes <- rownames(sig_df)
  cat("DEGs with TF connections:", length(intersect(deg_genes, unique(net$target))),
      "| without:", length(setdiff(deg_genes, unique(net$target))), "\n")

  # ---- Network: top 10 TFs by padj, DEG targets highlighted ----
  top_tfs <- tf_sig %>% dplyr::arrange(padj) %>% head(10) %>% dplyr::pull(source)
  edges   <- sig_targets %>%
    dplyr::filter(TranscriptionFactor %in% top_tfs) %>%
    dplyr::select(from = TranscriptionFactor, to = TargetGene, ModeOfRegulation)
  node_df <- dplyr::bind_rows(data.frame(name = unique(edges$from), type = "TF"),
                              data.frame(name = unique(edges$to),   type = "Target Gene")) %>%
    dplyr::mutate(Highlight_Group = dplyr::case_when(
      type == "TF"                                ~ "Transcription Factor Hub",
      type == "Target Gene" & name %in% deg_genes ~ "Significant DEG Target",
      TRUE                                        ~ "Non-DEG Target Gene"))

  graph <- tbl_graph(nodes = node_df, edges = edges, directed = TRUE) %>%
    activate("nodes") %>%
    dplyr::filter(!node_is_isolated())

  p_net <- ggraph(graph, layout = "stress") +
    geom_edge_diagonal(aes(color = factor(ModeOfRegulation)), alpha = 0.3, width = 0.4) +
    geom_node_point(aes(size = type, fill = Highlight_Group), shape = 21, color = "black", stroke = 0.3) +
    geom_node_text(aes(label = name, filter = (type == "TF")),
                   repel = TRUE, fontface = "bold", size = 6.0, color = "black", bg.color = "white", bg.r = 0.15) +
    geom_node_text(aes(label = name, filter = (Highlight_Group == "Significant DEG Target")),
                   repel = TRUE, fontface = "italic", size = 5, color = "red", bg.color = "white", bg.r = 0.08) +
    scale_edge_color_manual(values = c("-1" = "#E06666", "1" = "#6FA8DC"),
                            labels = c("-1" = "Repression (-1)", "1" = "Activation (+1)")) +
    scale_fill_manual(values = c("Transcription Factor Hub" = "#FFD966",
                                 "Significant DEG Target"   = "red",
                                 "Non-DEG Target Gene"      = "#E0E0E0")) +
    scale_size_manual(values = c("TF" = 7.0, "Target Gene" = 2.5), guide = "none") +
    labs(title = paste("Transcriptional Regulatory Landscape —", PREFIX), subtitle = net_subtitle,
         fill = "Biological Feature", edge_color = "Regulatory Impact") +
    theme_void() +
    theme(plot.title      = element_text(face = "bold", hjust = 0.5, size = 18, margin = margin(t = 25, b = 10)),
          plot.subtitle   = element_text(hjust = 0.5, colour = "grey30", size = 14, margin = margin(b = 15)),
          legend.position = "right",
          legend.title    = element_text(face = "bold", size = 14),
          legend.text     = element_text(size = 13))

  png(filename = paste0(OUTPUT_DIR, PREFIX, files$network), width = 15, height = 14, units = "in", res = 500)
  print(p_net)
  dev.off()

  # ---- TF-DEG overlap: grouped ----
  overlap <- sig_targets %>%
    dplyr::filter(TargetGene %in% deg_genes) %>%
    dplyr::left_join(tf_stats, by = "TranscriptionFactor") %>%
    dplyr::group_by(dplyr::across(all_of(group_cols))) %>%
    dplyr::summarise(dplyr::across(any_of("Confidence"), ~ paste(unique(.x), collapse = ", ")),
                     Overlap_DEG_Genes = paste(TargetGene, collapse = ", "),
                     DEG_Count         = dplyr::n(),
                     .groups           = "drop") %>%
    dplyr::arrange(TF_padj, TranscriptionFactor, dplyr::desc(ModeOfRegulation))
  print(head(overlap))
  write.csv(overlap, paste0(OUTPUT_DIR, PREFIX, files$overlap), row.names = FALSE)

  # ---- TF-DEG overlap: one row per TF-gene pair, with gene direction ----
  individual <- sig_targets %>%
    dplyr::filter(TargetGene %in% deg_genes) %>%
    dplyr::left_join(tf_stats, by = "TranscriptionFactor") %>%
    dplyr::mutate(Gene_Log2FC    = sig_df[TargetGene, "log2FoldChange"],
                  Gene_Higher_In = ifelse(Gene_Log2FC > 0, "PR (Responders)", "PD (Non-Responders)")) %>%
    dplyr::select(TranscriptionFactor, TF_Activity_Score, TF_p_value, TF_padj, TargetGene,
                  any_of("Confidence"), ModeOfRegulation, Gene_Log2FC, Gene_Higher_In) %>%
    dplyr::arrange(TF_padj, TranscriptionFactor, dplyr::desc(ModeOfRegulation))
  print(head(individual))
  write.csv(individual, paste0(OUTPUT_DIR, PREFIX, files$individual), row.names = FALSE)

  invisible(list(acts = acts, sig = tf_sig, targets = sig_targets))
}


# -----------------------------------------------------------------------------
# 5.6 Dual GSEA consensus (pre-ranked GSEA + ssGSEA)
# -----------------------------------------------------------------------------

# Pre-ranked GSEA used in the consensus step. A tiny random offset
# (sd = 1e-9) breaks tied statistics; eps = 0 gives exact small p-values.
run_pairwise_gsea <- function(res, t2g, label, min_size = 10, max_size = 500) {
  if (requireNamespace("BiocParallel", quietly = TRUE)) {
    BiocParallel::register(BiocParallel::SerialParam())   # single-threaded, avoids worker locks
  }

  gene_list <- res$stat
  names(gene_list) <- rownames(res)
  gene_list <- gene_list[!is.na(gene_list) & !is.na(names(gene_list))]
  gene_list <- gene_list[!duplicated(names(gene_list))]
  set.seed(SEED)                                          # reproducible tie-breaking
  gene_list <- gene_list + rnorm(length(gene_list), mean = 0, sd = 1e-9)
  gene_list <- sort(gene_list, decreasing = TRUE)

  gsea_obj <- GSEA(gene_list, exponent = 1, minGSSize = min_size, maxGSSize = max_size,
                   pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = t2g,
                   eps = 0, verbose = TRUE, seed = TRUE)
  list(result_obj = gsea_obj, df = as.data.frame(gsea_obj) %>% dplyr::mutate(Label = label))
}

# ssGSEA on log2(normalised counts + 1), then Wilcoxon PR vs PD per set.
# delta = mean(PR) - mean(PD).
run_ssgsea_stats <- function(norm_counts, meta, t2g, label, gene_col = "gene_symbol", min_size = 10) {
  gene_sets <- split(t2g[[gene_col]], t2g$gs_name)

  set.seed(SEED)
  param  <- ssgseaParam(exprData = log2(norm_counts + 1), geneSets = gene_sets,
                        normalize = TRUE, minSize = min_size, maxSize = 10000)
  scores <- gsva(param, verbose = FALSE)

  stats <- as.data.frame(scores) %>%
    rownames_to_column("Pathway") %>%
    pivot_longer(-Pathway, names_to = "SampleID", values_to = "ssGSEA_score") %>%
    dplyr::left_join(meta %>% rownames_to_column("SampleID") %>% dplyr::select(SampleID, Response),
                     by = "SampleID") %>%
    dplyr::group_by(Pathway) %>%
    dplyr::summarise(mean_PR = mean(ssGSEA_score[Response == "PR"], na.rm = TRUE),
                     mean_PD = mean(ssGSEA_score[Response == "PD"], na.rm = TRUE),
                     delta   = mean_PR - mean_PD,
                     p_value = wilcox.test(ssGSEA_score[Response == "PR"],
                                           ssGSEA_score[Response == "PD"])$p.value,
                     .groups = "drop") %>%
    dplyr::mutate(padj = p.adjust(p_value, method = "BH"), Label = label) %>%
    dplyr::arrange(padj)

  list(scores = scores, stats = stats)
}

# Join both methods and classify each pathway:
#   Consensus Activated / Suppressed : padj < 0.1 in both, same direction
#   Single-Method Significant        : padj < 0.1 in one method only (or disagreeing)
#   Not Significant                  : neither
# style = "standard" (MSigDB collections) or "custom" (DTP signatures plot style)
build_dual_consensus <- function(pairwise_df, ssgsea_stats_df, name_prefix, strip_pattern,
                                 out_dir, out_prefix, style = "standard") {

  dual_df <- dplyr::inner_join(
    pairwise_df %>% dplyr::select(Pathway = ID, Pairwise_NES = NES, Pairwise_padj = p.adjust),
    ssgsea_stats_df %>% dplyr::select(Pathway, ssGSEA_Delta = delta, ssGSEA_padj = padj),
    by = "Pathway"
  ) %>%
    dplyr::mutate(
      Direction_Agrees      = (Pairwise_NES * ssGSEA_Delta) > 0,
      Dual_Enrichment_Score = sign(Pairwise_NES) * sqrt(abs(Pairwise_NES * ssGSEA_Delta)),
      Consensus_Status = dplyr::case_when(
        Pairwise_padj < CONSENSUS_PADJ & ssGSEA_padj < CONSENSUS_PADJ & Direction_Agrees & Pairwise_NES > 0 ~ "Consensus Activated (PR)",
        Pairwise_padj < CONSENSUS_PADJ & ssGSEA_padj < CONSENSUS_PADJ & Direction_Agrees & Pairwise_NES < 0 ~ "Consensus Suppressed (PD)",
        Pairwise_padj < CONSENSUS_PADJ | ssGSEA_padj < CONSENSUS_PADJ ~ "Single-Method Significant",
        TRUE ~ "Not Significant"),
      Pathway_Clean = Pathway %>% gsub(strip_pattern, "", .) %>% gsub("_", " ", .)
    ) %>%
    # Draw order: grey first, consensus last (on top)
    dplyr::arrange(factor(Consensus_Status, levels = c("Not Significant", "Single-Method Significant",
                                                       "Consensus Activated (PR)", "Consensus Suppressed (PD)")))

  write.csv(dual_df, paste0(out_dir, out_prefix, "_", name_prefix, "_Dual_GSEA_Consensus_Matrix.csv"),
            row.names = FALSE)
  cat("\n[", name_prefix, "] Consensus status counts:\n")
  print(table(dual_df$Consensus_Status))

  is_std <- style == "standard"
  p <- ggplot(dual_df, aes(x = Pairwise_NES, y = ssGSEA_Delta)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(aes(color = Consensus_Status), size = if (is_std) 2.5 else 4, alpha = 0.8) +
    scale_color_manual(values = c("Consensus Activated (PR)"  = "#2b8cbe",
                                  "Consensus Suppressed (PD)" = "#de2d26",
                                  "Single-Method Significant" = "#FDAA48",
                                  "Not Significant"           = "grey85")) +
    theme_minimal(base_size = 13) +
    labs(title    = if (is_std) paste0("Dual GSEA Alignment Matrix — ", name_prefix, " | ", out_prefix)
                    else        paste0("Dual GSEA Alignment Matrix — ", name_prefix),
         subtitle = if (is_std) "Pairwise GSEA (NES) vs ssGSEA (Delta) — Consensus dots brought to the front layer"
                    else        NULL,
         x = if (is_std) "Pairwise GSEA — Normalized Enrichment Score"
             else        "Pairwise GSEA — Normalized Enrichment Score (NES)",
         y = "ssGSEA — Delta Enrichment Score (Mean PR - Mean PD)",
         color = "Consensus Status") +
    theme(plot.title      = element_text(face = "bold", hjust = 0.5, size = 14),
          plot.subtitle   = element_text(hjust = 0.5, color = "grey40", size = 10),
          legend.position = "bottom")

  print(p)
  ggsave(paste0(out_dir, out_prefix, "_", name_prefix, "_Dual_GSEA_Scatter.png"), p,
         width = if (is_std) 9 else 7, height = if (is_std) 7 else 6, dpi = if (is_std) 500 else 300)
  dual_df
}


#===============================================================================
# 6. STEP 1 — WHERE DOES STABLE DISEASE (SD) SIT? (PR / SD / PD context)
#===============================================================================

# ---- 6.1 DESeq2 object for the 3-group cohort (used for VST only) ----
dds2 <- DESeqDataSetFromMatrix(countData = filtered_GEX2, colData = merged_colData2, design = ~ Response)
dds2 <- dds2[rowSums(counts(dds2)) >= MIN_TOTAL_CNT, ]
vsd2 <- varianceStabilizingTransformation(dds2, blind = FALSE)

# ---- 6.2 PCA ----
pca_data   <- plotPCA(vsd2, intgroup = "Response", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p_pca_3groups <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Response)) +
  geom_point(size = 3.5, alpha = 0.8) +
  scale_color_manual(values = c("PD" = "#F8766D", "SD" = "#8EE081", "PR" = "#4DBBD5")) +
  labs(title    = "PCA Alignment: PR vs. SD vs. PD",
       subtitle = "Evaluating biological placement of Stable Disease (SD)",
       x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic() +
  theme(plot.title      = element_text(face = "bold", hjust = 0.5, size = 14),
        plot.subtitle   = element_text(hjust = 0.5, color = "grey40", size = 10),
        legend.position = "right")
print(p_pca_3groups)
ggsave(filename = paste0(OUTPUT_DIR, PREFIX, "_PCA_with_SD.png"), plot = p_pca_3groups,
       width = 7, height = 5, dpi = 500)

# ---- 6.3 Global PERMANOVA: do the 3 groups differ in overall expression? ----
expr_mat2    <- t(assay(vsd2))                         # samples x genes
dist_matrix3 <- dist(expr_mat2, method = "euclidean")
meta3        <- as.data.frame(colData(vsd2))

set.seed(SEED)
global_permanova <- adonis2(dist_matrix3 ~ Response, data = meta3, permutations = N_PERM)
print("--- GLOBAL PERMANOVA RESULTS ---")
print(global_permanova)

# ---- 6.4 Is SD closer to PR or to PD? (centroid distances in PC1-PC2) ----
centroids_3g <- pca_data %>%
  dplyr::group_by(Response) %>%
  dplyr::summarise(mean_PC1 = mean(PC1), mean_PC2 = mean(PC2))
print("--- GROUP CENTROIDS ---")
print(centroids_3g)

calc_dist <- function(g1, g2) {
  c1 <- centroids_3g %>% dplyr::filter(Response == g1)
  c2 <- centroids_3g %>% dplyr::filter(Response == g2)
  sqrt((c1$mean_PC1 - c2$mean_PC1)^2 + (c1$mean_PC2 - c2$mean_PC2)^2)
}
cat("\nDistance from SD to PD:", calc_dist("SD", "PD"), "\n")
cat("Distance from SD to PR:", calc_dist("SD", "PR"), "\n")
cat("Distance from PR to PD:", calc_dist("PR", "PD"), "\n")

# ---- 6.5 PERMDISP: does any group have more variable expression? ----
dispersion_test3 <- betadisper(dist_matrix3, group = meta3$Response, type = "centroid")
dispersion_p     <- permutest(dispersion_test3, permutations = N_PERM)
print("--- MULTIVARIATE DISPERSION RESULTS ---")
print(dispersion_p)


#===============================================================================
# 7. STEP 2 — COHORT OVERVIEW
#===============================================================================

# ---- 7.1 DESeq2 object for the main PR vs PD cohort ----
dds <- DESeqDataSetFromMatrix(countData = filtered_GEX, colData = merged_colData, design = ~ Response)
dds <- dds[rowSums(counts(dds)) >= MIN_TOTAL_CNT, ]

# ---- 7.2 Sample counts ----
plot_sample_counts(dds,
                   recode_map = c("PD" = "Progressive-Disease", "PR" = "Partial-Response"),
                   fill_cols  = c("Progressive-Disease" = "#B39DDB", "Partial-Response" = "#A5D6A7"),
                   title      = paste("Sample counts — PR Vs. PD at: Average", PREFIX),
                   file       = paste0(OUTPUT_DIR, PREFIX, "_Sample_Counts.png"))

plot_sample_counts(dds2,
                   recode_map = c("PD" = "Progressive-Disease", "SD" = "Stable-Disease", "PR" = "Partial-Response"),
                   fill_cols  = c("Progressive-Disease" = "#B39DDB", "Stable-Disease" = "#FDFD96",
                                  "Partial-Response" = "#A5D6A7"),
                   title      = "Sample counts - PR Vs. SD Vs. PD at 3-6week average",
                   file       = paste0(OUTPUT_DIR, PREFIX, "_Sample_Counts_SD.png"))

# ---- 7.3 Waterfall plots ----
plot_waterfall(merged_colData, RESP_COLS,
               title        = "Waterfall Plot: Irinotecan Response at 3-6 weeks average(PR Vs. PD)",
               xlab         = "PDX Models (Ordered by Response PR Vs. PD)",
               count_levels = c("PR", "PD"),
               file         = paste0(OUTPUT_DIR, PREFIX, "_waterfall.png"))

plot_waterfall(merged_colData2, RESP_COLS_SD,
               title        = "Waterfall Plot: Irinotecan Response at 3-6 weeks average (PR Vs. SD Vs. PD)",
               xlab         = "PDX Models (Ordered by Response)",
               count_levels = c("PR", "SD", "PD"),
               file         = paste0(OUTPUT_DIR, PREFIX, "_waterfall_SD.png"))

# ---- 7.4 Batch x response balance (confounding check) ----
print(table(merged_colData$batch_id, merged_colData$Response))

# ---- 7.5 CMS / CRIS subtype composition by response (displayed, not saved) ----
p_subtype <- merged_colData %>%
  dplyr::mutate(CMS  = ifelse(CMS  == "" | is.na(CMS)  | CMS  == "Unclassified", "UNC", CMS),
                CRIS = ifelse(CRIS == "" | is.na(CRIS) | CRIS == "Unclassified", "UNC", CRIS)) %>%
  dplyr::select(Response, CMS, CRIS) %>%
  pivot_longer(cols = c(CMS, CRIS), names_to = "Classification_System", values_to = "Subtype") %>%
  dplyr::filter(!is.na(Response) & Subtype != "NA") %>%
  ggplot(aes(x = Response, fill = Subtype)) +
  geom_bar(position = "fill", color = "black", linewidth = 0.3, width = 0.5) +
  facet_wrap(~ Classification_System, scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = c("CMS1" = "#66C2A5", "CMS2" = "#FC8D62", "CMS3" = "#8DA0CB", "CMS4" = "#E78AC3",
                               "CRIS-A" = "#A6D854", "CRIS-B" = "#FFD92F", "CRIS-C" = "#E5C494",
                               "CRIS-D" = "#B3B3B3", "CRIS-E" = "#999999",
                               "UNC" = "#4D4D4D")) +                      # unclassified: dark grey
  labs(title = "Subtype Distribution by Response Group",
       subtitle = sprintf("PR (n=%d) vs. PD (n=%d)", N_PR, N_PD),
       x = "Clinical Response", y = "Percentage (%)", fill = "Subtype Cluster") +
  theme_classic(base_size = 14) +
  theme(strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 14),
        plot.title       = element_text(face = "bold", hjust = 0.5),
        plot.subtitle    = element_text(hjust = 0.5, color = "gray40"),
        panel.spacing    = unit(2, "lines"))
print(p_subtype)

# ---- 7.6 Subtype vs response association (Fisher, Monte-Carlo p, B = 10,000) ----
cms_fisher  <- fisher.test(table(merged_colData$CMS,  merged_colData$Response), simulate.p.value = TRUE, B = 10000)
cat("CMS Subtype vs Response Exact p-value:", cms_fisher$p.value, "\n")
cris_fisher <- fisher.test(table(merged_colData$CRIS, merged_colData$Response), simulate.p.value = TRUE, B = 10000)
cat("CRIS Subtype vs Response Exact p-value:", cris_fisher$p.value, "\n")


#===============================================================================
# 8. STEP 3 — DIFFERENTIAL EXPRESSION (DESeq2, PR vs PD)
#===============================================================================

# ---- 8.1 Fit and extract (log2FC = PR / PD) ----
set.seed(SEED)
dds <- DESeq(dds)
res <- results(dds, contrast = c("Response", "PR", "PD"))
print(resultsNames(dds))
write.csv(as.data.frame(res), paste0(PIPELINES_DIR, PREFIX, "_DESeq2_results.csv"))

# ---- 8.2 Variance-stabilised expression (PERMANOVA / NMDS / TEAD steps) ----
vsd <- vst(dds, blind = FALSE)

# ---- 8.3 apeglm-shrunken LFCs — quick-look volcano (displayed, not saved) ----
# Shrinkage pulls noisy, low-count fold changes towards zero.
res_lfc <- lfcShrink(dds, coef = "Response_PR_vs_PD", type = "apeglm")

lfc_df <- as.data.frame(res_lfc) %>% dplyr::mutate(Gene = rownames(.))
lfc_top10 <- lfc_df %>% dplyr::filter(!is.na(pvalue)) %>% dplyr::arrange(pvalue) %>% head(10) %>% dplyr::pull(Gene)
lfc_df <- lfc_df %>% dplyr::mutate(Label = ifelse(Gene %in% lfc_top10, Gene, ""))

print(
  ggplot(lfc_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(aes(color = (padj < 0.05 & abs(log2FoldChange) > 1)), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("darkgray", "red")) +
    geom_text_repel(aes(label = Label), size = 4, max.overlaps = Inf, box.padding = 0.5) +
    theme_minimal() +
    labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10(p-value)") +
    theme(legend.position = "none")
)


#===============================================================================
# 9. STEP 4 — DEG SUMMARIES
#===============================================================================

# ---- 9.1 DEG counts across a grid of thresholds ----
resOrdered  <- res[order(res$padj), ]
padj_levels <- c(0.2, 0.1, 0.05)
fc_levels   <- c(0, 0.585, 1)
fc_labels   <- c("no FC filter", "FC > 1.5", "FC > 2")

threshold_summary <- expand.grid(padj_cutoff = padj_levels, fc_cutoff = fc_levels) %>%
  dplyr::mutate(
    Count_Total = mapply(function(p, fc) sum(resOrdered$padj < p & abs(resOrdered$log2FoldChange) > fc, na.rm = TRUE), padj_cutoff, fc_cutoff),
    Count_UP    = mapply(function(p, fc) sum(resOrdered$padj < p & resOrdered$log2FoldChange >  fc, na.rm = TRUE), padj_cutoff, fc_cutoff),
    Count_DOWN  = mapply(function(p, fc) sum(resOrdered$padj < p & resOrdered$log2FoldChange < -fc, na.rm = TRUE), padj_cutoff, fc_cutoff),
    fc_label    = fc_labels[match(fc_cutoff, fc_levels)],
    Threshold   = factor(paste0("padj<", padj_cutoff,
                                ifelse(fc_cutoff == 0, "", paste0(" & ", fc_labels[match(fc_cutoff, fc_levels)])))),
    Group       = paste0("padj < ", padj_cutoff)
  )

# Unadjusted p-value rows for reference
pvalue_rows <- data.frame(
  padj_cutoff = NA, fc_cutoff = 0,
  Count_Total = c(sum(resOrdered$pvalue < 0.05, na.rm = TRUE), sum(resOrdered$pvalue < 0.01, na.rm = TRUE)),
  Count_UP    = c(sum(resOrdered$pvalue < 0.05 & resOrdered$log2FoldChange > 0, na.rm = TRUE),
                  sum(resOrdered$pvalue < 0.01 & resOrdered$log2FoldChange > 0, na.rm = TRUE)),
  Count_DOWN  = c(sum(resOrdered$pvalue < 0.05 & resOrdered$log2FoldChange < 0, na.rm = TRUE),
                  sum(resOrdered$pvalue < 0.01 & resOrdered$log2FoldChange < 0, na.rm = TRUE)),
  fc_label = "no FC filter", Threshold = factor(c("pvalue<0.05", "pvalue<0.01")), Group = "pvalue"
)
threshold_summary <- dplyr::bind_rows(pvalue_rows, threshold_summary)

threshold_cols <- c("pvalue" = "#CC79A7", "padj < 0.2" = "#5DC863FF",
                    "padj < 0.1" = "#0072B2", "padj < 0.05" = "#F46A25")

p_total <- ggplot(threshold_summary, aes(x = Threshold, y = Count_Total, fill = Group)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_text(aes(label = Count_Total), vjust = -0.4, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = threshold_cols) +
  labs(title = paste("DEG counts by threshold —", PREFIX), x = "Threshold", y = "Number of genes", fill = "Group") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(OUTPUT_DIR, PREFIX, "_Gene_count_by_Threshold.png"), p_total, width = 9, height = 5, dpi = 300)

threshold_long <- threshold_summary %>%
  pivot_longer(c(Count_UP, Count_DOWN), names_to = "Direction", values_to = "Count") %>%
  dplyr::mutate(Direction  = dplyr::recode(Direction, "Count_UP" = "UP (higher in PR)", "Count_DOWN" = "DOWN (lower in PR)"),
                Count_plot = ifelse(grepl("DOWN", Direction), -Count, Count))

p_updown <- ggplot(threshold_long, aes(x = Threshold, y = Count_plot, fill = Direction)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_text(aes(label = abs(Count), vjust = ifelse(Count_plot >= 0, -0.4, 1.2)), size = 3, fontface = "bold") +
  scale_fill_manual(values = c("UP (higher in PR)" = "#D55E00", "DOWN (lower in PR)" = "#0072B2")) +
  geom_hline(yintercept = 0) + scale_y_continuous(labels = abs) +
  labs(title = paste("DEG direction by threshold —", PREFIX), x = "Threshold", y = "Number of genes") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(OUTPUT_DIR, PREFIX, "_Gene_count_UP_DOWN_by_Threshold.png"), p_updown, width = 9, height = 5.5, dpi = 300)

# ---- 9.2 MA plot ----
ma_data <- as.data.frame(res) %>%
  dplyr::filter(!is.na(padj), !is.na(log2FoldChange), baseMean > 0) %>%
  dplyr::mutate(status = dplyr::case_when(
                  padj < PADJ_CUTOFF & log2FoldChange >  LFC_CUTOFF ~ "UP",
                  padj < PADJ_CUTOFF & log2FoldChange < -LFC_CUTOFF ~ "DOWN",
                  TRUE ~ "NS"),
                gene = rownames(.))

ma_labels <- dplyr::bind_rows(
  ma_data %>% dplyr::filter(status == "UP")   %>% dplyr::arrange(padj) %>% dplyr::slice_head(n = 10),
  ma_data %>% dplyr::filter(status == "DOWN") %>% dplyr::arrange(padj) %>% dplyr::slice_head(n = 10))
n_up   <- sum(ma_data$status == "UP")
n_down <- sum(ma_data$status == "DOWN")
n_ns   <- sum(ma_data$status == "NS")

p_ma <- ggplot(ma_data, aes(x = log10(baseMean), y = log2FoldChange, color = status)) +
  geom_point(data = dplyr::filter(ma_data, status == "NS"), alpha = 0.3, size = 0.8) +
  geom_point(data = dplyr::filter(ma_data, status != "NS"), alpha = 0.8, size = 1.2) +
  scale_color_manual(values = c("NS" = "grey60", "UP" = "#D55E00", "DOWN" = "#0072B2"),
                     labels = c("NS"   = paste0("Not significant (n=", n_ns, ")"),
                                "UP"   = paste0("UP in PR (n=", n_up, ")"),
                                "DOWN" = paste0("DOWN in PR (n=", n_down, ")")), name = NULL) +
  geom_hline(yintercept = 0,           color = "black",   linewidth = 0.5) +
  geom_hline(yintercept =  LFC_CUTOFF, color = "#D55E00", linewidth = 0.4, linetype = "dashed") +
  geom_hline(yintercept = -LFC_CUTOFF, color = "#0072B2", linewidth = 0.4, linetype = "dashed") +
  geom_text_repel(data = ma_labels, aes(label = gene), size = 2.8, max.overlaps = 20,
                  segment.color = "grey50", segment.size = 0.3, box.padding = 0.3, show.legend = FALSE) +
  labs(title    = paste("MA Plot (without batch effect) —", PREFIX, "| PR vs PD"),
       subtitle = paste0("padj < ", PADJ_CUTOFF, "  &  |log2FC| > ", LFC_CUTOFF, " | Top 10 UP/DOWN labelled"),
       x = "Mean expression [log10(baseMean)]", y = "log2 Fold Change (PR / PD)") +
  theme_classic(base_size = 12) +
  theme(plot.title         = element_text(hjust = 0.5, face = "bold", size = 13),
        plot.subtitle      = element_text(hjust = 0.5, color = "grey40", size = 9),
        legend.position    = "bottom",
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3))
ggsave(paste0(OUTPUT_DIR, PREFIX, "_MAPlot.png"), p_ma, width = 10, height = 7, dpi = 500)

# ---- 9.3 Volcano plot (five gene classes, top 20 UP / DOWN labelled) ----
res_filtered <- as.data.frame(res) %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::mutate(
    diffexpressed = dplyr::case_when(
      padj <  PADJ_CUTOFF & log2FoldChange >  LFC_CUTOFF      ~ "UP (Sensitive)",
      padj <  PADJ_CUTOFF & log2FoldChange < -LFC_CUTOFF      ~ "DOWN (Resistant)",
      padj >= PADJ_CUTOFF & abs(log2FoldChange) >  LFC_CUTOFF ~ "High Change, Low Sig",
      padj <  PADJ_CUTOFF & abs(log2FoldChange) <= LFC_CUTOFF ~ "Low Change, High Sig",
      TRUE ~ "NS"),
    delabel = ifelse(diffexpressed %in% c("UP (Sensitive)", "DOWN (Resistant)"), rownames(.), NA)
  )

total_up   <- sum(res_filtered$diffexpressed == "UP (Sensitive)")
total_down <- sum(res_filtered$diffexpressed == "DOWN (Resistant)")
volcano_labels <- rbind(
  res_filtered %>% dplyr::filter(diffexpressed == "UP (Sensitive)")   %>% dplyr::arrange(padj, dplyr::desc(log2FoldChange)) %>% dplyr::slice_head(n = 20),
  res_filtered %>% dplyr::filter(diffexpressed == "DOWN (Resistant)") %>% dplyr::arrange(padj, log2FoldChange)              %>% dplyr::slice_head(n = 20))
ymax <- max(-log10(res_filtered$padj), na.rm = TRUE) + 1

p_volcano <- ggplot(res_filtered, aes(x = log2FoldChange, y = -log10(padj),
                                      color = diffexpressed, size = diffexpressed, alpha = diffexpressed)) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major = element_line(color = "grey95", linewidth = 0.5),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle    = element_text(color = "grey40", size = 10, hjust = 0.5),
        legend.position  = "right",
        legend.title     = element_blank()) +
  geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_point() +
  scale_color_manual(values = c("UP (Sensitive)" = "#C62828", "DOWN (Resistant)" = "#2E7D32",
                                "High Change, Low Sig" = "#FDAA48", "Low Change, High Sig" = "#89CFF0",
                                "NS" = "#CCCCCC")) +
  scale_alpha_manual(values = c("UP (Sensitive)" = 0.85, "DOWN (Resistant)" = 0.85,
                                "High Change, Low Sig" = 0.60, "Low Change, High Sig" = 0.60, "NS" = 0.20)) +
  scale_size_manual(values = c("UP (Sensitive)" = 2.5, "DOWN (Resistant)" = 2.5,
                               "High Change, Low Sig" = 1.8, "Low Change, High Sig" = 1.8, "NS" = 1.0)) +
  geom_label_repel(data = volcano_labels, aes(label = delabel),
                   max.overlaps = 15, size = 3.5, fontface = "italic", fill = "white",
                   label.size = 0.1, label.padding = unit(0.18, "lines"), box.padding = 0.35,
                   segment.color = "grey40", segment.linewidth = 0.3, show.legend = FALSE) +
  annotate("text", x =  Inf, y = ymax, label = paste("Upregulated in PR (Sensitive):", total_up),
           color = "#C62828", fontface = "bold", hjust = 1.1) +
  annotate("text", x = -Inf, y = ymax, label = paste("Upregulated in PD (Resistant):", total_down),
           color = "#2E7D32", fontface = "bold", hjust = -0.1) +
  labs(title = paste("Volcano Plot (without batch effect) —", PREFIX),
       subtitle = sprintf("PR (n=%d) vs PD (n=%d)", N_PR, N_PD),
       x = "log2(Fold Change)", y = "-log10(padj)")
ggsave(paste0(OUTPUT_DIR, PREFIX, "_VolcanoPlot.png"), p_volcano, width = 11, height = 9, dpi = 500)

# ---- 9.4 Normalised counts and DEG table ----
normCount <- counts(dds, normalized = TRUE)                     # size-factor normalised
write.csv(normCount, paste0(PIPELINES_DIR, PREFIX, "_Normalized_counts.csv"), quote = FALSE, row.names = TRUE)

sig <- res_filtered %>% dplyr::filter(diffexpressed %in% c("UP (Sensitive)", "DOWN (Resistant)"))
write.csv(sig, paste0(OUTPUT_DIR, PREFIX, "_DEGs_PR_vs_PD.csv"), quote = FALSE, row.names = TRUE)

deg_up_genes   <- rownames(sig)[sig$diffexpressed == "UP (Sensitive)"]
deg_down_genes <- rownames(sig)[sig$diffexpressed == "DOWN (Resistant)"]

top_sig_genes <- sig %>% dplyr::arrange(padj) %>% head(TOP_N_HEATMAP)
write.csv(top_sig_genes, paste0(OUTPUT_DIR, PREFIX, "_Top_", TOP_N_HEATMAP, "_DEGs_Details.csv"),
          quote = FALSE, row.names = TRUE)


#===============================================================================
# 10. STEP 5 — HEATMAPS OF DEGs
#     Input: log2(normalised counts + 1), z-scored per gene.
#===============================================================================

# ---- 10.1 Z-score matrices (all DEGs and top 50 by padj) ----
allSig_z <- row_zscore(log2(normCount[rownames(normCount) %in% rownames(sig), ] + 1))
allSig_z <- allSig_z[, intersect(colnames(allSig_z), rownames(merged_colData)), drop = FALSE]

topSig_z <- allSig_z[rownames(allSig_z) %in% rownames(top_sig_genes), , drop = FALSE]
topSig_z <- topSig_z[, intersect(colnames(topSig_z), rownames(merged_colData)), drop = FALSE]

# ---- 10.2 pheatmap-style PDFs (ComplexHeatmap::pheatmap wrapper) ----
legacy_heat_cols <- colorRampPalette(c("navy", "white", "firebrick"))(100)
legacy_ann_cols  <- list(Response_Groups = c("PR" = "#009E73", "PD" = "#D55E00"))

# All DEGs: annotate Response (+ CMS, CRIS)
ann_all <- as.data.frame(merged_colData[colnames(allSig_z), c("Response", "CMS", "CRIS"), drop = FALSE])
colnames(ann_all)[1] <- "Response_Groups"
rownames(ann_all)    <- colnames(allSig_z)

pdf(paste0(OUTPUT_DIR, PREFIX, "_Heatmap_Grouped.pdf"), width = 17, height = 14)
ComplexHeatmap::pheatmap(allSig_z, scale = "none",
                         main = "Grouped clustering of differentially expressed genes at 3–6 weeks",
                         fontsize = 16, fontsize_row = 7, color = legacy_heat_cols,
                         show_rownames = FALSE, show_colnames = FALSE,
                         annotation_col = ann_all, annotation_colors = legacy_ann_cols,
                         column_split = as.factor(ann_all$Response_Groups),
                         name = "Z-score", border_color = NA)
dev.off()

pdf(paste0(OUTPUT_DIR, PREFIX, "_Heatmap_Hierarchical.pdf"), width = 17, height = 14)
ComplexHeatmap::pheatmap(allSig_z, scale = "none",
                         main = "Hierarchical clustering of differentially expressed genes at 3–6 weeks",
                         color = legacy_heat_cols, show_rownames = FALSE, show_colnames = FALSE,
                         annotation_col = ann_all, annotation_colors = legacy_ann_cols,
                         name = "Z-score", fontsize_row = 7, fontsize = 16, border_color = NA)
dev.off()

# Top 50 DEGs: annotate Response only
ann_top <- data.frame(Response_Groups = merged_colData[colnames(topSig_z), "Response"],
                      row.names = colnames(topSig_z))

pdf(paste0(OUTPUT_DIR, PREFIX, "_TopGenes_Heatmap_Grouped.pdf"), width = 10, height = 10)
ComplexHeatmap::pheatmap(topSig_z, scale = "none",
                         main = paste0("Grouped Clustering of top", TOP_N_HEATMAP, "DEGs - ", PREFIX),
                         color = legacy_heat_cols, show_rownames = TRUE, show_colnames = FALSE,
                         annotation_col = ann_top, annotation_colors = legacy_ann_cols,
                         column_split = as.factor(ann_top$Response_Groups),
                         name = "Z-score", fontsize_row = 10, border_color = NA)
dev.off()

pdf(paste0(OUTPUT_DIR, PREFIX, "_TopGenes_Heatmap_Hierarchical.pdf"), width = 10, height = 10)
ComplexHeatmap::pheatmap(topSig_z, scale = "none",
                         main = paste0("Hierarchial Clustering of top", TOP_N_HEATMAP, "DEGs - ", PREFIX),
                         color = legacy_heat_cols, show_rownames = TRUE, show_colnames = FALSE,
                         annotation_col = ann_top, annotation_colors = legacy_ann_cols,
                         name = "Z-score", fontsize_row = 10, border_color = NA)
dev.off()

# ---- 10.3 ComplexHeatmap JPGs — annotated by CRIS, Response, CMS ----
annotation_df <- merged_colData[colnames(allSig_z), c("CRIS", "Response", "CMS"), drop = FALSE]
annotation_df$Response <- factor(clean_unclassified(annotation_df$Response), levels = c("PR", "PD", "Unclassified"))
annotation_df$CMS      <- factor(clean_unclassified(annotation_df$CMS),
                                 levels = c("CMS1", "CMS2", "CMS3", "CMS4", "Unclassified"))
annotation_df$CRIS     <- factor(clean_unclassified(annotation_df$CRIS),
                                 levels = c("CRIS-A", "CRIS-B", "CRIS-C", "CRIS-D", "CRIS-E", "Unclassified"))

# Colour-blind-friendly annotation palettes
ann_colors <- list(
  CRIS     = c("CRIS-A" = "#882255", "CRIS-B" = "#44AA99", "CRIS-C" = "#999933",
               "CRIS-D" = "#117733", "CRIS-E" = "#DDCC77", "Unclassified" = "#D3D3D3"),
  Response = c("PR" = "#009E73", "PD" = "#E69F00", "Unclassified" = "#D3D3D3"),
  CMS      = c("CMS1" = "#56B4E9", "CMS2" = "#332288", "CMS3" = "#D55E00",
               "CMS4" = "#CC79A7", "Unclassified" = "#D3D3D3")
)

make_col_annotation <- function(ann) {
  HeatmapAnnotation(CRIS = ann$CRIS, Response = ann$Response, CMS = ann$CMS, col = ann_colors)
}
deg_direction <- function(genes) ifelse(genes %in% deg_up_genes, "UP (Sensitive)", "DOWN (Resistant)")

big_title_gp  <- gpar(fontsize = 24, fontface = "bold")
legend_params <- list(title_gp = gpar(fontsize = 16, fontface = "bold"), labels_gp = gpar(fontsize = 14))

# Heatmap 1 — all DEGs, columns split by Response x CMS, rows split by direction
ht_grouped <- Heatmap(allSig_z, name = "Z-score", col = HEAT_COL_FUN,
                      column_title = "Grouped clustering of DEGs (3-6 weeks)\nannotated by CRIS, Response, CMS",
                      column_title_gp = big_title_gp,
                      column_split = annotation_df[, c("Response", "CMS")],
                      row_split = deg_direction(rownames(allSig_z)), row_title_gp = gpar(fontsize = 16),
                      top_annotation = make_col_annotation(annotation_df),
                      show_row_names = FALSE, show_column_names = FALSE,
                      cluster_row_slices = TRUE, cluster_column_slices = TRUE,
                      heatmap_legend_param = legend_params)
jpeg(paste0(OUTPUT_DIR, PREFIX, "_Heatmap_Grouped_CMS.jpg"), width = 17, height = 17, units = "in", res = 500)
draw(ht_grouped, padding = unit(c(10, 10, 10, 10), "mm"))
dev.off()

# Heatmap 2 — all DEGs, unsupervised clustering
ht_hier <- Heatmap(allSig_z, name = "Z-score", col = HEAT_COL_FUN,
                   column_title = "Hierarchical clustering of DEGs (3-6 weeks)\nannotated by CRIS, Response, CMS",
                   column_title_gp = big_title_gp,
                   top_annotation = make_col_annotation(annotation_df),
                   show_row_names = FALSE, show_column_names = FALSE,
                   heatmap_legend_param = legend_params)
jpeg(paste0(OUTPUT_DIR, PREFIX, "_Heatmap_Hierarchical.jpg"), width = 17, height = 14, units = "in", res = 500)
draw(ht_hier, padding = unit(c(10, 10, 10, 10), "mm"))
dev.off()

# Heatmap 3 — top 50 DEGs, columns split by Response x CRIS
annotation_df_top <- annotation_df[colnames(topSig_z), , drop = FALSE]
jpeg(paste0(OUTPUT_DIR, PREFIX, "_TopGenes_Heatmap_Grouped.jpg"), width = 10, height = 10, units = "in", res = 500)
draw(Heatmap(topSig_z, name = "Z-score", col = HEAT_COL_FUN,
             column_title = paste0("Grouped clustering of top ", TOP_N_HEATMAP, " DEGs - ", PREFIX,
                                   "\nannotated by CRIS, Response, CMS"),
             row_split = deg_direction(rownames(topSig_z)),
             column_split = annotation_df_top[, c("Response", "CRIS")],
             top_annotation = make_col_annotation(annotation_df_top),
             show_row_names = TRUE, show_column_names = FALSE, row_names_gp = gpar(fontsize = 7),
             cluster_row_slices = FALSE, cluster_column_slices = TRUE))
dev.off()

# Heatmap 4 — top 50 DEGs, unsupervised clustering
jpeg(paste0(OUTPUT_DIR, PREFIX, "_TopGenes_Heatmap_Hierarchical.jpg"), width = 10, height = 10, units = "in", res = 500)
draw(Heatmap(topSig_z, name = "Z-score", col = HEAT_COL_FUN,
             column_title = paste0("Hierarchical clustering of top ", TOP_N_HEATMAP, " DEGs - ", PREFIX,
                                   "\nannotated by CRIS, Response, CMS"),
             top_annotation = make_col_annotation(annotation_df_top),
             show_row_names = TRUE, show_column_names = FALSE, row_names_gp = gpar(fontsize = 7)))
dev.off()

# ---- 10.4 UP-DEG signature score (mean z) by CMS subtype ----
sample_scores <- data.frame(
  SampleID   = colnames(allSig_z),
  UP_score   = colMeans(allSig_z[rownames(allSig_z) %in% deg_up_genes,   , drop = FALSE]),
  DOWN_score = colMeans(allSig_z[rownames(allSig_z) %in% deg_down_genes, , drop = FALSE])
) %>%
  dplyr::left_join(merged_colData %>% dplyr::mutate(SampleID = rownames(.)) %>%
                     dplyr::select(SampleID, Response, CMS, CRIS), by = "SampleID")

p_up_cms <- ggplot(sample_scores, aes(x = CMS, y = UP_score, fill = CMS)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.6) +
  stat_compare_means(method = "kruskal.test") +
  theme_classic(base_size = 13) +
  labs(title = "UP (Sensitive) signature score by CMS subtype", y = "Mean Z-score of UP DEGs", x = "CMS")
ggsave(paste0(OUTPUT_DIR, PREFIX, "_UPscore_by_CMS_boxplot.pdf"), p_up_cms, width = 7, height = 5)


#===============================================================================
# 11. STEP 6 — GLOBAL TRANSCRIPTOMIC STRUCTURE (PR vs PD)
#===============================================================================

# ---- 11.1 PERMANOVA on the most variable genes (VST, Euclidean) ----
# Restricting to high-variance genes removes low-count background noise.
rv          <- matrixStats::rowVars(assay(vsd))
top_var_idx <- order(rv, decreasing = TRUE)[1:min(N_TOP_VAR_GENES, length(rv))]
dist_matrix <- dist(t(assay(vsd)[top_var_idx, ]), method = "euclidean")
meta        <- as.data.frame(colData(dds))

set.seed(SEED)
permanova_PRPD <- adonis2(dist_matrix ~ Response, data = meta, permutations = N_PERM)
print(permanova_PRPD)

# Sensitivity of R2 to the number of genes used
for (n_genes in c(500, 1000, 2000, 5000)) {
  top_g     <- order(rv, decreasing = TRUE)[1:min(n_genes, length(rv))]
  test_dist <- dist(t(assay(vsd)[top_g, ]), method = "euclidean")
  res_perm  <- adonis2(test_dist ~ Response, data = meta, permutations = 199)
  cat("Genes:", n_genes, " | R2:", res_perm$R2[1], " | p-val:", res_perm$`Pr(>F)`[1], "\n")
}

# ---- 11.2 NMDS ordination (k = 3 dimensions; first two plotted) ----
set.seed(SEED)                                   # NMDS starts from a random configuration
nmds_results <- metaMDS(dist_matrix, k = 3, trymax = 100, autotransform = FALSE)
print(nmds_results)
cat("NMDS Stress Value is:", nmds_results$stress, "\n")   # < 0.2 is usually considered acceptable

nmds_plot_data <- cbind(as.data.frame(scores(nmds_results, display = "sites")), meta)
write.csv(nmds_plot_data, paste0(OUTPUT_DIR, "NMDS_Coordinates_with_Metadata.csv"), row.names = TRUE)

# Centroid / spider plot
centroids <- nmds_plot_data %>%
  dplyr::group_by(Response) %>%
  dplyr::summarise(Centroid_X = mean(NMDS1, na.rm = TRUE), Centroid_Y = mean(NMDS2, na.rm = TRUE), .groups = "drop")
nmds_centroid_data <- nmds_plot_data %>% dplyr::inner_join(centroids, by = "Response")

centroid_plot <- ggplot(nmds_centroid_data) +
  geom_segment(aes(x = Centroid_X, y = Centroid_Y, xend = NMDS1, yend = NMDS2, color = Response),
               alpha = 0.3, linewidth = 0.6) +
  geom_point(aes(x = NMDS1, y = NMDS2, color = Response), size = 2.5, alpha = 0.7) +
  geom_point(data = centroids, aes(x = Centroid_X, y = Centroid_Y, fill = Response),
             size = 6, shape = 23, color = "black", stroke = 1.5) +
  theme_bw(base_size = 14) +
  labs(title    = "NMDS Centroid & Spider Plot of 5-FU/Irinotecan Response",
       subtitle = paste("Global Partitioning (PERMANOVA p = 0.012) | Stress:", round(nmds_results$stress, 4)),
       x = "NMDS Dimension 1", y = "NMDS Dimension 2",
       color = "Individual Samples", fill = "Group Centroid") +
  theme(plot.title       = element_text(face = "bold", hjust = 0.5),
        plot.subtitle    = element_text(hjust = 0.5, face = "italic"),
        panel.grid.minor = element_blank(),
        legend.position  = "right")
print(centroid_plot)
ggsave(filename = paste0(OUTPUT_DIR, "NMDS_Centroid_Spider_Plot.pdf"), plot = centroid_plot,
       width = 8, height = 6, dpi = 300)

# ---- 11.3 PERMDISP: is a PERMANOVA difference driven by unequal spread? ----
dispersion_test <- betadisper(dist_matrix, group = meta$Response, type = "centroid")
set.seed(SEED)
print(permutest(dispersion_test, permutations = N_PERM))

pdf(file = paste0(OUTPUT_DIR, "NMDS_Dispersion_Boxplot.pdf"), width = 8, height = 6)
boxplot(dispersion_test, main = "Multivariate Dispersion (Distance to Centroid)",
        xlab = "Response Group", ylab = "Distance")
dev.off()


#===============================================================================
# 12. STEP 7 — GENE-LEVEL VALIDATION (edgeR TMM log2CPM)
#     Independent normalisation to confirm top DESeq2 hits are not artefacts.
#===============================================================================

dge         <- calcNormFactors(DGEList(counts = filtered_GEX), method = "TMM")
log_cpm_mat <- cpm(dge, log = TRUE, prior.count = 1)

# ---- 12.1 Top 9 DEGs by padj ----
top9_genes <- sig %>% dplyr::arrange(padj) %>% head(9) %>% rownames()
plot_cpm_boxplots(top9_genes, log_cpm_mat, sig, merged_colData,
                  title    = paste("Top 9 DEGs — Log2 CPM (edgeR TMM) |", PREFIX),
                  subtitle = sprintf("PR (n = %d) vs PD (n = %d) | edgeR TMM-normalised CPM", N_PR, N_PD),
                  file     = paste0(OUTPUT_DIR, PREFIX, "_Top9_DEGs_CPM_Boxplot.png"))

# ---- 12.2 Top 7 DEGs + TEAD targets ----
top7_genes     <- sig %>% dplyr::filter(!rownames(.) %in% TEAD_TARGETS) %>% dplyr::arrange(padj) %>% head(7) %>% rownames()
combined_genes <- c(top7_genes, TEAD_TARGETS)
plot_cpm_boxplots(combined_genes, log_cpm_mat, sig, merged_colData,
                  title    = paste("Top DEGs & TEAD Targets — Log2 CPM (edgeR TMM) |", PREFIX),
                  subtitle = sprintf("PR (n = %d) vs PD (n = %d) | Top 7 DEGs + %s", N_PR, N_PD,
                                     paste(TEAD_TARGETS, collapse = " & ")),
                  file     = paste0(OUTPUT_DIR, PREFIX, "_Top7_and_TEAD_Targets_CPM_Boxplot.png"))


#===============================================================================
# 13. STEP 8 — PRE-RANKED GSEA
#===============================================================================

# ---- 13.1 Ranked list: Wald statistic, genes with a padj (as saved to disk) ----
DESeq_results <- read.csv(paste0(PIPELINES_DIR, PREFIX, "_DESeq2_results.csv"), header = TRUE, row.names = 1)
DESeq_results <- DESeq_results[!is.na(DESeq_results$padj), ]    # drop independent-filtered genes
geneList      <- sort(setNames(DESeq_results$stat, rownames(DESeq_results)), decreasing = TRUE)

# ---- 13.2 Hallmark, KEGG, C6 oncogenic (all sets tested; pvalueCutoff = 1) ----
gsea_h <- run_gsea_collection(geneList, h_t2g, list(
  tag = "Hallmarks", curve_tag = "Hallmarks", curve_title = "Hallmark Top Hit - ",
  all_title = "Enrichment of Hallmark Pathways in Responders Vs. Non-Responders",
  all_size = c(10, 8), all_x_title_size = 14,
  sig_suffix = "_GSEA_barplot_Significant.png", sig_height = 8.5,
  sig_title = "Enrichment of Significant Hallmark Pathways\n in PR Vs. PD", sig_title_size = 14))

gsea_kegg <- run_gsea_collection(geneList, kegg_t2g, list(
  tag = "KEGG", curve_tag = "KEGG", curve_title = "KEGG - ",
  all_title = "Enrichment of KEGG Pathways in PR Vs. PD",
  all_size = c(12, 12), all_x_title_size = 10,
  sig_suffix = "_GSEA_barplot_significant.png", sig_height = 8.5,
  sig_title = "Enrichment of Significant KEGG Pathways\n in PR Vs. PD", sig_title_size = 16))

gsea_c6 <- run_gsea_collection(geneList, c6_t2g, list(
  tag = "C6_Oncogenic", curve_tag = "C6", curve_title = "C6 - ",
  all_title = "Enrichment of Oncogenic Pathways in PR Vs. PD",
  all_size = c(10, 8), all_x_title_size = 14,
  sig_suffix = "_GSEA_barplot_Significant.png", sig_height = 12,
  sig_title = "Enrichment of Significant Oncogenic Pathways\nin PR Vs. PD", sig_title_size = 14))

# ---- 13.3 C5 GO (only significant sets returned: pvalueCutoff = 0.05) ----
c5_gseaResult <- GSEA(geneList, exponent = 1, minGSSize = 10, maxGSSize = 10000,
                      pvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE = c5_t2g,
                      verbose = TRUE, seed = TRUE)
c5_df <- as.data.frame(c5_gseaResult)
write.csv(c5_df, file = paste0(OUTPUT_DIR, PREFIX, "_C5_GO_GSEA_RESULTS.csv"))

c5_df <- c5_df[order(-c5_df$NES), ]
write.csv(c5_df, file = paste0(OUTPUT_DIR, PREFIX, "_C5_GO_GSEA_RESULTS_Sorted by NES.csv"))

if (nrow(c5_df) > 0) {
  # Top 25 by NES (i.e. most positively enriched in PR)
  png(file = paste0(OUTPUT_DIR, PREFIX, "_C5_GO_GSEA_barplot.png"), width = 11, height = 8.5, units = "in", res = 500)
  print(
    ggplot(head(c5_df, 25), aes(reorder(Description, NES), NES)) +
      geom_col(fill = "#29BF4E", width = 0.7) +
      coord_flip() +
      labs(x = "Pathway", y = "Normalized Enrichment Score",
           title = "Enrichment of Top 25 Significant GO Pathways\nin PR Vs. PD") +
      gsea_bar_theme(title_size = 18, x_title_size = 16, y_text_size = 9) +
      theme(axis.title.y = element_text(size = 16, colour = "black", face = "bold"))
  )
  dev.off()

  c5_pathwayID <- c5_df$ID[1]
  png(file = paste0(OUTPUT_DIR, PREFIX, "_C5_GO_TOP_", gsub("[[:punct:]]", "_", c5_pathwayID), ".png"),
      width = 10, height = 8, units = "in", res = 500)
  print(gseaplot2(c5_gseaResult, c5_pathwayID, title = paste0("GO Top Hit - ", c5_pathwayID),
                  color = "green", base_size = 15, subplots = 1:3))
  dev.off()
} else {
  message("No significant C5 pathways found.")
}


#===============================================================================
# 14. STEP 9 — TRANSCRIPTION FACTOR ACTIVITY (ULM on the DESeq2 statistic)
#===============================================================================

# ---- 14.1 DoRothEA (A/B/C) ----
tf_dorothea <- run_tf_pipeline(
  net = dorothea_net, res = res, sig_df = sig,
  files = list(acts       = "_TF_activity_ULM.csv",
               bar        = "_TF_activity_barplot.png",
               targets    = "_Gene_Targets_for_Sig_TFs.csv",
               grouped    = "_Grouped_TF_Targets_With_Significance.csv",
               network    = "_Top_10_TF_DEG_Network.png",
               overlap    = "_TF_Significant_DEG_Overlap_Summary.csv",
               individual = "_TF_Significant_DEG_Individual_With_Direction.csv"),
  bar_subtitle = "DoRothEA A/B/C | positive score = more active in PR | FDR<0.05",
  net_subtitle = "Top 10 TFs by padj, with Highlighted DEGs")
tf_acts <- tf_dorothea$acts          # reused for TEAD labels (Step 12)

# ---- 14.2 CollecTRI ----
tf_collectri <- run_tf_pipeline(
  net = ct_net, res = res, sig_df = sig,
  files = list(acts       = "_CollecTRI_TF_activity_ULM.csv",
               bar        = "_CollecTRI_TF_activity_barplot.png",
               targets    = "_CollecTRI_Gene_Targets_for_Sig_TFs.csv",
               grouped    = "_CollecTRI_Grouped_TF_Targets_With_Significance.csv",
               network    = "_CollecTRI_Top_10_TF_DEG_Network.png",
               overlap    = "_CollecTRI_TF_DEG_Overlap_Summary.csv",
               individual = "_CollecTRI_TF_DEG_Individual_With_Direction.csv"),
  bar_subtitle = "CollecTRI | positive score = more active in PR | FDR<0.05",
  net_subtitle = "CollecTRI | Top 10 TFs by padj, with Highlighted DEGs")


#===============================================================================
# 15. STEP 10 — DUAL GSEA CONSENSUS (pre-ranked GSEA vs ssGSEA)
#     A pathway is a robust hit when the rank-based and the sample-level
#     methods agree in direction and are both significant (padj < 0.1).
#===============================================================================

CONSENSUS_SETS <- list(
  Hallmark     = list(t2g = h_t2g,    strip = "^HALLMARK_"),
  KEGG         = list(t2g = kegg_t2g, strip = "^KEGG_"),
  C5_GO        = list(t2g = c5_t2g,   strip = "^GOBP_|^GOCC_|^GOMF_"),
  C6_Oncogenic = list(t2g = c6_t2g,   strip = "^")            # C6 names share no common prefix
)

dual_results <- list()
for (set_name in names(CONSENSUS_SETS)) {
  cs        <- CONSENSUS_SETS[[set_name]]
  pairwise  <- run_pairwise_gsea(res, cs$t2g, label = set_name)
  ss        <- run_ssgsea_stats(normCount, merged_colData, cs$t2g, label = set_name)
  dual_results[[set_name]] <- build_dual_consensus(pairwise$df, ss$stats, set_name, cs$strip,
                                                   out_dir = OUTPUT_DIR, out_prefix = PREFIX)
}


#===============================================================================
# 16. STEP 11 — DRUG-TOLERANT PERSISTER (DTP) SIGNATURES
#===============================================================================

# ---- 16.1 Quick look: separate ssGSEA scores, distribution and UP-DOWN correlation ----
set.seed(SEED)
sep_ssgsea_scores <- gsva(ssgseaParam(exprData = log2(normCount + 1),
                                      geneSets = list(DTP_UP   = dtp_t2g$gene[dtp_t2g$gs_name == "DTP_Up_Signature"],
                                                      DTP_DOWN = dtp_t2g$gene[dtp_t2g$gs_name == "DTP_Down_Signature"]),
                                      normalize = TRUE, minSize = 5, maxSize = 10000),
                          verbose = FALSE)

sep_ssgsea_long <- sep_ssgsea_scores %>%
  as_tibble(rownames = "Signature") %>%
  pivot_longer(cols = -Signature, names_to = "SampleID", values_to = "ssGSEA_Score") %>%
  dplyr::left_join(merged_colData %>% as_tibble(rownames = "SampleID") %>% dplyr::select(SampleID, Response),
                   by = "SampleID")

print(
  ggplot(sep_ssgsea_long, aes(x = Response, y = ssGSEA_Score, fill = Response)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
    geom_jitter(width = 0.15, size = 1.5, color = "black", alpha = 0.6) +
    facet_wrap(~ Signature, scales = "free_y") +
    scale_fill_manual(values = c("PR" = "#2b8cbe", "PD" = "#de2d26")) +
    scale_y_continuous(labels = label_number(accuracy = 0.01)) +
    theme_bw(base_size = 13) +
    labs(title = "Enrichment Distribution: Up vs Down Signatures",
         subtitle = "Calculated separately per signature across Patient Responses",
         x = "Patient Clinical Response Status", y = "Enrichment Score (ssGSEA)") +
    theme(plot.title = element_text(face = "bold", size = 14),
          strip.text = element_text(face = "bold", size = 12),
          panel.grid.minor = element_blank(), legend.position = "none")
)

# Are the UP and DOWN signatures independent (orthogonal) across samples?
correlation_df <- sep_ssgsea_long %>%
  dplyr::select(SampleID, Signature, ssGSEA_Score, Response) %>%
  pivot_wider(names_from = Signature, values_from = ssGSEA_Score)
cor_test <- cor.test(correlation_df$DTP_UP, correlation_df$DTP_DOWN, method = "pearson")
cor_coef <- round(cor_test$estimate, 3)
cor_p    <- format.pval(cor_test$p.value, digits = 3)
cat("Pearson R (DTP_UP vs DTP_DOWN):", cor_coef, "| p =", cor_p, "\n")

p_cor <- ggplot(correlation_df, aes(x = DTP_UP, y = DTP_DOWN)) +
  geom_point(aes(color = Response), size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE) +
  scale_color_manual(values = c("PR" = "#2b8cbe", "PD" = "#de2d26")) +
  theme_bw(base_size = 13) +
  labs(title = "Signature Orthogonality Check",
       subtitle = paste0("Pearson R = ", cor_coef, " (p = ", cor_p, ")"),
       x = "DTP_UP ssGSEA Score", y = "DTP_DOWN ssGSEA Score")
print(p_cor)
ggsave("./Signature_Correlation_Diagnostic.png", p_cor, width = 6, height = 5, dpi = 300)

# ---- 16.2 Dual consensus for the DTP signatures ----
custom_pairwise <- run_pairwise_gsea(res, dtp_t2g, label = "Custom_DTP", min_size = 5, max_size = 1000)
custom_ssgsea   <- run_ssgsea_stats(normCount, merged_colData, dtp_t2g, label = "Custom_DTP",
                                    gene_col = "gene", min_size = 5)
custom_dual     <- build_dual_consensus(custom_pairwise$df, custom_ssgsea$stats, "Custom_DTP", "",
                                        out_dir = "./", out_prefix = "DTP_Analysis", style = "custom")

# ---- 16.3 Final DTP ssGSEA scores (unique genes) + Wilcoxon PR vs PD ----
set.seed(SEED)
dtp_scores <- gsva(ssgseaParam(exprData = log2(normCount + 1),
                               geneSets = list(DTP_UP = unique(dtp_up_genes), DTP_DOWN = unique(dtp_down_genes)),
                               normalize = TRUE, minSize = 5, maxSize = 10000),
                   verbose = FALSE)

dtp_long <- as.data.frame(dtp_scores) %>%
  rownames_to_column("Signature") %>%
  pivot_longer(-Signature, names_to = "SampleID", values_to = "ssGSEA_score") %>%
  dplyr::left_join(merged_colData %>% rownames_to_column("SampleID") %>% dplyr::select(SampleID, Response),
                   by = "SampleID") %>%
  dplyr::mutate(Response  = factor(Response,  levels = c("PD", "PR")),
                Signature = factor(Signature, levels = c("DTP_UP", "DTP_DOWN")))
write.csv(dtp_long, paste0(OUTPUT_DIR, PREFIX, "_DTP_UP_DOWN_ssGSEA_scores.csv"), row.names = FALSE)

dtp_stats <- dtp_long %>%
  dplyr::group_by(Signature) %>%
  dplyr::summarise(mean_PR = mean(ssGSEA_score[Response == "PR"], na.rm = TRUE),
                   mean_PD = mean(ssGSEA_score[Response == "PD"], na.rm = TRUE),
                   delta   = mean_PR - mean_PD,
                   p_value = wilcox.test(ssGSEA_score[Response == "PR"], ssGSEA_score[Response == "PD"])$p.value,
                   .groups = "drop") %>%
  dplyr::mutate(padj = p.adjust(p_value, method = "BH"))
write.csv(dtp_stats, paste0(OUTPUT_DIR, PREFIX, "_DTP_UP_DOWN_ssGSEA_stats.csv"), row.names = FALSE)
print(dtp_stats)

p_dtp_box <- ggplot(dtp_long, aes(x = Response, y = ssGSEA_score, fill = Response)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.75, colour = "grey30") +
  geom_jitter(aes(colour = Response), width = 0.15, size = 1.5, alpha = 0.8, show.legend = FALSE) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 3.5) +
  facet_wrap(~ Signature, scales = "free_y") +
  scale_fill_manual(values = RESP_COLS, labels = c("PD" = "Progressive Disease", "PR" = "Partial Response")) +
  scale_colour_manual(values = RESP_POINT_COLS) +
  labs(title    = paste("DTP Signature ssGSEA Scores — PR vs PD |", PREFIX),
       subtitle = "Drug-tolerant persister signatures (HCT & SW, day 14 vs day 0), scored independently",
       x = NULL, y = "ssGSEA Enrichment Score", fill = "Response") +
  theme_classic(base_size = 12) +
  theme(plot.title       = element_text(face = "bold", hjust = 0.5, size = 13),
        plot.subtitle    = element_text(hjust = 0.5, colour = "grey40", size = 9),
        strip.text       = element_text(face = "bold", size = 10),
        strip.background = element_rect(fill = "grey95", colour = NA),
        legend.position  = "bottom")
print(p_dtp_box)
ggsave(paste0(OUTPUT_DIR, PREFIX, "_DTP_UP_DOWN_ssGSEA_Boxplot.png"), p_dtp_box, width = 8, height = 5, dpi = 500)

# UP vs DOWN per sample (high UP + low DOWN = persister-like phenotype)
dtp_wide <- as.data.frame(t(dtp_scores)) %>%
  rownames_to_column("SampleID") %>%
  dplyr::left_join(merged_colData %>% rownames_to_column("SampleID") %>% dplyr::select(SampleID, Response),
                   by = "SampleID") %>%
  dplyr::mutate(Response = factor(Response, levels = c("PD", "PR")))

p_dtp_scatter <- ggplot(dtp_wide, aes(x = DTP_UP, y = DTP_DOWN, colour = Response)) +
  geom_point(size = 2.2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey40", linewidth = 0.4, linetype = "dashed") +
  scale_colour_manual(values = RESP_COLS, labels = c("PD" = "Progressive Disease", "PR" = "Partial Response")) +
  labs(title = paste("DTP UP vs DOWN Signature Scores per Sample |", PREFIX),
       x = "DTP_UP ssGSEA score", y = "DTP_DOWN ssGSEA score", colour = "Response") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12), legend.position = "bottom")
print(p_dtp_scatter)
ggsave(paste0(OUTPUT_DIR, PREFIX, "_DTP_UP_vs_DOWN_Scatter.png"), p_dtp_scatter, width = 6.5, height = 5.5, dpi = 500)

# ---- 16.4 Outlier sensitivity: drop 1.5 x IQR outliers per group and re-test ----
flag_boxplot_outliers <- function(df) {
  q1  <- quantile(df$ssGSEA_score, 0.25, na.rm = TRUE)
  q3  <- quantile(df$ssGSEA_score, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  df$Is_Outlier <- df$ssGSEA_score < (q1 - 1.5 * iqr) | df$ssGSEA_score > (q3 + 1.5 * iqr)
  df
}

dtp_flagged <- dtp_long %>%
  dplyr::group_by(Signature, Response) %>%
  dplyr::group_modify(~ flag_boxplot_outliers(.x)) %>%
  dplyr::ungroup()

detected_outliers <- dtp_flagged %>% dplyr::filter(Is_Outlier)
cat("\n=== AUTOMATICALLY DETECTED BOXPLOT OUTLIERS ===\n")
if (nrow(detected_outliers) == 0) {
  cat("No structural outliers detected via the 1.5*IQR rule!\n")
} else {
  print(as.data.frame(detected_outliers %>% dplyr::select(SampleID, Signature, Response, ssGSEA_score)))
}

dtp_stats_clean <- dtp_flagged %>%
  dplyr::filter(!Is_Outlier) %>%
  dplyr::group_by(Signature) %>%
  dplyr::summarise(delta_clean   = mean(ssGSEA_score[Response == "PR"], na.rm = TRUE) -
                                   mean(ssGSEA_score[Response == "PD"], na.rm = TRUE),
                   p_value_clean = wilcox.test(ssGSEA_score[Response == "PR"],
                                               ssGSEA_score[Response == "PD"])$p.value,
                   .groups = "drop")

sensitivity_matrix <- dtp_stats %>%
  dplyr::select(Signature, Original_Delta = delta, Original_p = p_value) %>%
  dplyr::left_join(dtp_stats_clean, by = "Signature")
cat("\n=== DYNAMIC SENSITIVITY MATRIX ===\n")
print(as.data.frame(sensitivity_matrix))


#===============================================================================
# 17. STEP 12 — TEAD1 / TEAD4 SAMPLE-LEVEL ACTIVITY (DoRothEA ULM per sample)
#===============================================================================

# ---- 17.1 On VST expression (displayed, not saved) ----
sample_tf_vst <- decoupleR::run_ulm(mat = assay(vsd), net = dorothea_net,
                                    .source = "source", .target = "target", .mor = "mor")

tead_vst_df <- sample_tf_vst %>%
  dplyr::filter(source %in% c("TEAD1", "TEAD4")) %>%
  dplyr::rename(SampleID = condition, TF = source, Activity = score) %>%
  dplyr::left_join(as.data.frame(colData(vsd)) %>% rownames_to_column("SampleID"), by = "SampleID")

print(
  ggplot(tead_vst_df, aes(x = Response, y = Activity, fill = Response)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, size = 1.5, aes(color = Response)) +
    facet_wrap(~ TF) +
    theme_minimal() +
    scale_fill_manual(values = RESP_COLS) +
    scale_color_manual(values = c("PD" = "#D65A52", "PR" = "#3A8FA3")) +
    labs(title    = "TEAD Transcription Factor Activity Across Patient Groups",
         subtitle = sprintf("Sample-level analysis (PR n=%d vs PD n=%d)", N_PR, N_PD),
         x = "Clinical Response Status", y = "ULM Activity Score (vsd)")
)

# ---- 17.2 On edgeR log2CPM (saved; labelled with contrast-level ULM stats) ----
sample_tf_cpm <- decoupleR::run_ulm(mat = log_cpm_mat, net = dorothea_net,
                                    .source = "source", .target = "target", .mor = "mor")

tead_labels <- tf_acts %>%
  dplyr::filter(source %in% c("TEAD1", "TEAD4")) %>%
  dplyr::rename(TF = source) %>%
  dplyr::mutate(padj_label = dplyr::case_when(padj < 0.001 ~ "***", padj < 0.01 ~ "**",
                                              padj < 0.05  ~ "*",   TRUE ~ "ns"),
                padj_text  = paste0("padj = ", formatC(padj, format = "e", digits = 2)),
                score_text = paste0("Contrast Score = ", round(score, 2))) %>%
  dplyr::select(TF, padj_label, padj_text, score_text)

plot_tead_df <- sample_tf_cpm %>%
  dplyr::filter(source %in% c("TEAD1", "TEAD4")) %>%
  dplyr::rename(SampleID = condition, TF = source, Activity = score) %>%
  dplyr::left_join(merged_colData %>% rownames_to_column("SampleID") %>% dplyr::select(SampleID, Response),
                   by = "SampleID") %>%
  dplyr::mutate(TF       = factor(TF, levels = c("TEAD1", "TEAD4")),
                Response = factor(Response, levels = c("PD", "PR"))) %>%
  dplyr::left_join(tead_labels, by = "TF") %>%
  dplyr::mutate(TF_label = paste0(TF, "\n(", padj_text, ")"))

p_tead_boxplot <- ggplot(plot_tead_df, aes(x = Response, y = Activity, fill = Response)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.75, colour = "grey30", linewidth = 0.4) +
  geom_jitter(aes(colour = Response), width = 0.15, size = 1.5, alpha = 0.8, show.legend = FALSE) +
  geom_text(data = plot_tead_df %>% dplyr::distinct(TF_label, padj_label),
            aes(x = 1.5, y = Inf, label = padj_label), inherit.aes = FALSE, vjust = 1.5, size = 5, fontface = "bold") +
  geom_text(data = plot_tead_df %>% dplyr::distinct(TF_label, score_text),
            aes(x = 1.5, y = Inf, label = score_text), inherit.aes = FALSE, vjust = 3.5, size = 3, colour = "black") +
  scale_fill_manual(values = RESP_COLS, labels = c("PD" = "Progressive Disease", "PR" = "Partial Response")) +
  scale_colour_manual(values = RESP_POINT_COLS) +
  facet_wrap(~ TF_label, scales = "free_y", ncol = 2) +
  labs(title    = paste("TEAD Transcription Factor Activity Summary |", PREFIX),
       subtitle = sprintf("PR (n = %d) vs PD (n = %d) | Sample-level ULM scores using edgeR log2CPM", N_PR, N_PD),
       x = NULL, y = "ULM Activity Score", fill = "Response") +
  theme_classic(base_size = 11) +
  theme(plot.title       = element_text(face = "bold", hjust = 0.5, size = 13),
        plot.subtitle    = element_text(hjust = 0.5, colour = "grey40", size = 9),
        strip.text       = element_text(face = "bold", size = 9),
        strip.background = element_rect(fill = "grey95", colour = NA),
        legend.position  = "bottom",
        axis.text.x      = element_text(size = 9),
        panel.spacing    = unit(0.8, "lines"))
print(p_tead_boxplot)
ggsave(filename = paste0(OUTPUT_DIR, PREFIX, "_TEAD_Activity_edgeR_Boxplot.png"), plot = p_tead_boxplot,
       width = 8, height = 6, dpi = 500)


#===============================================================================
# 18. REPRODUCIBILITY
#===============================================================================

sessionInfo()
