################################################################################
#
#  Script      : Long_overall_timepoint.R
#  Project     : Transcriptomic profiling of FOLFIRI / 5-FU response in
#                colorectal cancer (CRC) patient-derived xenograft (PDX) models
#  Analysis    : Treatment / timepoint effect (block: trt_tf) — paired,
#                within-model comparisons between treatment arms
#  Author      : Janmita Kaverimane Umesh (MSc Bioinformatics, QUB)
#
#  Comparisons (label = "<comparator>_vs_<reference>")
#    - 24h_vs_Placebo
#    - 6wks_vs_Placebo
#    - 6wks_vs_24h
#    DESeq2 contrast = c("Treatment", comparator, reference), so a positive
#    log2FC / NES / TF score = higher in the comparator (first-named) arm.
#
#  Workflow
#  ---------------------------------------------------------------------------
#    1. Load data and QC (sample matching, pair availability, duplicates)
#    2. Matched pairs — PDX models sampled in both arms of a comparison
#    3. Paired DESeq2 (~ Case_ID + Treatment) for each comparison
#    4. Sample-level QC — PCA, PC1 outliers, library size, sensitivity analysis
#    5. MA and volcano plots
#    6. Pre-ranked GSEA (Hallmark / KEGG / GO) + combined figures
#    7. Over-representation analysis (ORA) of DEGs
#    8. ssGSEA per-sample Hallmark scores
#    9. TF activity — DoRothEA (A/B/C) regulons + decoupleR ULM
#
#  Inputs  (0_data/PDX_Longitudinal/counts/)
#    - CRC_graft_raw_counts_matrix_12022024.csv  : raw counts, genes x samples
#    - Complete_GraftCRC_colData_19022024.csv    : sample metadata
#
#  Outputs
#    - 2_pipelines/PDX_Longitudinal/ : fitted DESeq2 objects + full DE tables
#    - 3_output/PDX_Longitudinal/LONG_TF/ : plots, enrichment and TF tables
#
#  Notes
#    - All thresholds and parameters are unchanged from the original script.
#
################################################################################


#===============================================================================
# 0. SESSION SETUP
#===============================================================================

# ---- 0.1 One-off package installation (set TRUE on a new machine) ----
INSTALL_PACKAGES <- FALSE

if (INSTALL_PACKAGES) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(c("DESeq2", "clusterProfiler", "GSVA", "decoupleR", "dorothea",
                         "glmGamPoi"), update = FALSE, ask = FALSE)
  install.packages(c("tidyverse", "msigdbr", "ggrepel", "patchwork"))
}

# ---- 0.2 Libraries ----
# Core DE
library(DESeq2)          # differential expression (NB GLM); glmGamPoi must be installed
# Enrichment
library(clusterProfiler) # GSEA and ORA (enricher)
library(msigdbr)         # MSigDB gene-set collections
library(GSVA)            # single-sample GSEA (ssGSEA)
library(dorothea)        # TF -> target regulons
library(decoupleR)       # TF activity inference (ULM)
# Data handling and plotting
library(tidyverse)       # dplyr, tidyr, tibble, stringr, ggplot2
library(ggrepel)         # non-overlapping gene labels
library(patchwork)       # multi-panel figure assembly

# ---- 0.3 Project paths ----
setwd("C:/Users/40496110/Downloads/PDX_PROJECT")

PIPELINES_DIR <- "./2_pipelines/PDX_Longitudinal/"
OUTPUT_DIR    <- "./3_output/PDX_Longitudinal/LONG_TF/"
COUNTS_DIR    <- "./0_data/PDX_Longitudinal/counts/"

COUNTS_FILE   <- "CRC_graft_raw_counts_matrix_12022024.csv"
COLDATA_FILE  <- "Complete_GraftCRC_colData_19022024.csv"

dir.create(PIPELINES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_DIR,    recursive = TRUE, showWarnings = FALSE)


#===============================================================================
# 1. ANALYSIS CONFIGURATION
#===============================================================================

BLOCK_ID      <- "trt_tf"

# ---- 1.1 Comparisons: label -> c(comparator, reference) ----
TRT_LEVELS    <- c("Placebo", "24h", "6wks")
COMPARISONS   <- list(
  "24h_vs_Placebo"  = c(comparator = "24h",  reference = "Placebo"),
  "6wks_vs_Placebo" = c(comparator = "6wks", reference = "Placebo"),
  "6wks_vs_24h"     = c(comparator = "6wks", reference = "24h")
)

# ---- 1.2 Thresholds and parameters ----
SEED          <- 123
MIN_TOTAL_CNT <- 10          # pre-filter: keep genes with >= 10 reads in total
PADJ_THRESH   <- 0.05        # BH-adjusted significance threshold
FC_THRESH     <- 0.585       # |log2FC| > 0.585  ==  fold change > 1.5
N_PERM_SIMPLE <- 100000      # nPermSimple for every GSEA() call

# ---- 1.3 Plot styling ----
BAR_FILL      <- "#29BF4E"   # single-colour GSEA bar plots
BASE_SIZE     <- 11
TF_UP_COLOR   <- "#00C4CC"   # teal  — active in comparator (first-named arm)
TF_DOWN_COLOR <- "#F4796B"   # coral — active in reference (second-named arm)

# ---- 1.4 Sensitivity analysis ----
# Sample flagged as an extreme PC1 outlier with very low library size in the
# 24h vs Placebo PCA (Section 8). Its whole PDX line is dropped in the
# sensitivity re-fits to check results are not driven by it.
OUTLIER_SAMPLE <- "CRC0081LMX0B02202TUMFF0100"

set.seed(SEED)   # set once, globally; GSEA() additionally uses seed = TRUE


#===============================================================================
# 2. LOAD DATA AND QC
#===============================================================================

# ---- 2.1 Raw counts (genes x samples) and sample metadata ----
GEX <- read.csv(paste0(COUNTS_DIR, COUNTS_FILE),  row.names = 1, check.names = FALSE)
TRT <- read.csv(paste0(COUNTS_DIR, COLDATA_FILE), row.names = 1, check.names = FALSE)

cat("GEX:", nrow(GEX), "genes x", ncol(GEX), "samples\n")
cat("TRT:", nrow(TRT), "rows\n")
cat("Treatment levels:\n")
print(table(TRT$Treatment))

# ---- 2.2 Sample matching between counts and metadata ----
cat("Samples in GEX but not TRT:", length(setdiff(colnames(GEX), rownames(TRT))), "\n")
cat("Samples in TRT but not GEX:", length(setdiff(rownames(TRT), colnames(GEX))), "\n")

# ---- 2.3 Factor coding (Placebo = reference level) ----
TRT$Treatment <- factor(TRT$Treatment, levels = TRT_LEVELS)
TRT$Case_ID   <- as.factor(TRT$Case_ID)

# ---- 2.4 Keep samples with counts; order counts to match metadata ----
TRT <- TRT[rownames(TRT) %in% colnames(GEX), ]
GEX <- GEX[, rownames(TRT)]

cat("TRT post-filtering:", nrow(TRT), "\n")
cat("TRT matches GEX:", setequal(rownames(TRT), colnames(GEX)), "\n")
cat("Treatment counts after filtering:\n")
print(table(TRT$Treatment))

# ---- 2.5 How many PDX models have both arms of each comparison? ----
pairs_present <- TRT %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  dplyr::select(SampleID, Case_ID, Treatment) %>%
  dplyr::group_by(Case_ID) %>%
  dplyr::summarise(has_Placebo = "Placebo" %in% Treatment,
                   has_24h     = "24h"     %in% Treatment,
                   has_6wks    = "6wks"    %in% Treatment,
                   .groups = "drop")

cat("Case_IDs with Placebo & 24h: ",  sum(pairs_present$has_Placebo & pairs_present$has_24h),  "\n")
cat("Case_IDs with Placebo & 6wks: ", sum(pairs_present$has_Placebo & pairs_present$has_6wks), "\n")
cat("Case_IDs with 24h & 6wks: ",     sum(pairs_present$has_24h     & pairs_present$has_6wks), "\n")

# ---- 2.6 Duplicate Case_ID x Treatment combinations (expected: 0) ----
dup_check <- TRT %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  dplyr::count(Case_ID, Treatment) %>%
  dplyr::filter(n > 1)

cat("Case_ID x Treatment combos with >1 sample:", nrow(dup_check), "\n")
print(dup_check)


#===============================================================================
# 3. REFERENCE RESOURCES
#===============================================================================

# ---- 3.1 MSigDB gene sets (TERM2GENE format: gs_name, gene_symbol) ----
h_t2g <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, gene_symbol)

kegg_t2g <- msigdbr(species = "Homo sapiens", collection = "C2",
                    subcollection = "CP:KEGG_LEGACY") %>%
  dplyr::select(gs_name, gene_symbol)

c5_t2g <- msigdbr(species = "Homo sapiens", collection = "C5") %>%
  dplyr::select(gs_name, gene_symbol)   # note: C5 also contains HPO sets

# Per-collection settings
#   gsea_top_n : pathways in the single-comparison GSEA bar plot (NULL = all significant)
#   panel_top_n: pathways in each panel of the combined figure
GENESET_COLLECTIONS <- list(
  Hallmark = list(t2g = h_t2g,    gsea_top_n = NULL, panel_top_n = 12),
  KEGG     = list(t2g = kegg_t2g, gsea_top_n = NULL, panel_top_n = 12),
  GO       = list(t2g = c5_t2g,   gsea_top_n = 25,   panel_top_n = 12)
)

# ---- 3.2 DoRothEA regulons (high / medium confidence: A, B, C) ----
data(dorothea_hs, package = "dorothea")

DOROTHEA_NET <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B", "C")) %>%
  dplyr::rename(source = tf)

cat("DoRothEA columns:", colnames(DOROTHEA_NET), "\n")   # source, confidence, target, mor


#===============================================================================
# 4. HELPER FUNCTIONS
#===============================================================================

# Split "6wks_vs_24h" into c("6wks", "24h") for direction-aware labels
split_label <- function(label) strsplit(label, "_vs_")[[1]]


# -----------------------------------------------------------------------------
# 4.1 Matched pairs and paired differential expression
# -----------------------------------------------------------------------------

# Keep PDX models (Case_ID) present in BOTH arms; one sample per model per arm
# (first if replicated). level_A = reference, level_B = comparator.
build_matched_pairs <- function(trt, gex, level_A, level_B) {

  sub_meta <- trt %>%
    as.data.frame() %>%
    rownames_to_column("SampleID") %>%
    dplyr::filter(Treatment %in% c(level_A, level_B))

  ids_with_both <- sub_meta %>%
    dplyr::group_by(Case_ID) %>%
    dplyr::summarise(has_A = level_A %in% Treatment,
                     has_B = level_B %in% Treatment,
                     .groups = "drop") %>%
    dplyr::filter(has_A & has_B) %>%
    dplyr::pull(Case_ID)

  sub_meta <- sub_meta %>%
    dplyr::filter(Case_ID %in% ids_with_both) %>%
    dplyr::group_by(Case_ID, Treatment) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    column_to_rownames("SampleID")

  sub_meta$Treatment <- factor(sub_meta$Treatment, levels = c(level_A, level_B))
  sub_meta$Case_ID   <- factor(sub_meta$Case_ID)

  list(meta = sub_meta, counts = gex[, rownames(sub_meta)])
}

# Paired design ~ Case_ID + Treatment: Case_ID absorbs model-to-model
# variation, so the Treatment coefficient is the within-model change.
run_deseq2_paired <- function(pair, label, pipelines_dir, min_count = MIN_TOTAL_CNT) {

  dds <- DESeqDataSetFromMatrix(countData = pair$counts,
                                colData   = pair$meta,
                                design    = ~ Case_ID + Treatment)
  dds <- dds[rowSums(counts(dds)) >= min_count, ]
  cat(sprintf("[%s] Genes going into DESeq2: %d | Samples: %d\n", label, nrow(dds), ncol(dds)))

  dds <- DESeq(dds, fitType = "glmGamPoi")   # glmGamPoi: fast, stable dispersion fit

  saveRDS(dds, paste0(pipelines_dir, "dds_", label, "_fitted.rds"))
  cat(sprintf("[%s] Fitted object saved to disk.\n", label))
  dds
}

# Extract comparator vs reference results and save the full table.
extract_results <- function(dds, comparator, reference, pipelines_dir, label) {
  res <- results(dds, contrast = c("Treatment", comparator, reference))

  cat(sprintf("\n[%s] resultsNames:\n", label))
  print(resultsNames(dds))
  summary(res)

  write.csv(as.data.frame(res), paste0(pipelines_dir, label, "_DESeq2_results.csv"))
  cat(sprintf("[%s] Results saved.\n", label))
  res
}


# -----------------------------------------------------------------------------
# 4.2 Sample-level QC
# -----------------------------------------------------------------------------

# VST (design-aware, blind = FALSE) + PCA coloured by arm; saved and returned.
plot_pca_labeled <- function(dds, label, output_dir) {
  vsd        <- vst(dds, blind = FALSE)
  pca_data   <- plotPCA(vsd, intgroup = "Treatment", returnData = TRUE)
  percentVar <- round(100 * attr(pca_data, "percentVar"))

  p <- ggplot(pca_data, aes(PC1, PC2, color = Treatment)) +
    geom_point(size = 2.5, alpha = 0.85) +
    labs(title = paste("PCA:", label),
         x = paste0("PC1: ", percentVar[1], "% variance"),
         y = paste0("PC2: ", percentVar[2], "% variance")) +
    theme_classic(base_size = BASE_SIZE) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(paste0(output_dir, label, "_PCA.png"), p, width = 7, height = 6, dpi = 400)
  print(p)
  list(vsd = vsd, pca_data = pca_data, plot = p)
}

# Report the n samples furthest from the origin on PC1 with their library sizes.
detect_pca_outliers <- function(pca_data, dds, n = 3) {
  pct_var <- round(100 * attr(pca_data, "percentVar"), 1)   # read before subsetting drops the attribute
  top     <- pca_data[order(-abs(pca_data$PC1)), ][seq_len(n), ]

  lib_sizes <- colSums(counts(dds))
  lib_df    <- data.frame(SampleID = names(lib_sizes), TotalReads = lib_sizes)

  cat(sprintf("Variance explained: PC1 = %.1f%% | PC2 = %.1f%%\n", pct_var[1], pct_var[2]))
  cat("\nTop PC1 outliers:\n");                 print(top[, c("name", "PC1", "PC2")])
  cat("\nLibrary size summary:\n");             print(summary(lib_sizes))
  cat("\nLibrary sizes for flagged outliers:\n"); print(lib_df[lib_df$SampleID %in% top$name, ])

  invisible(list(outliers = top, lib_sizes = lib_df, pct_var = pct_var))
}

# Number of genes for which a sample's Cook's distance exceeds the
# F(0.99; p, m - p) cut-off DESeq2 uses to flag count outliers.
count_cooks_outliers <- function(dds, sample_id) {
  p   <- ncol(model.matrix(design(dds), colData(dds)))   # number of model coefficients
  cut <- qf(0.99, p, ncol(dds) - p)
  sum(assays(dds)[["cooks"]][, sample_id] > cut, na.rm = TRUE)
}

# Genes passing padj and |log2FC| thresholds
sig_genes <- function(r) {
  rownames(r)[which(r$padj < PADJ_THRESH & abs(r$log2FoldChange) > FC_THRESH)]
}

# Re-fit after dropping the PDX line(s) that contain bad_samples, then compare
# with the full fit: Spearman correlation of Wald statistics and overlap of
# significant genes. 'suffix' is appended to the label for the saved .rds.
sensitivity_check <- function(pair, bad_samples, comparator, reference, res_full,
                              label, pipelines_dir, suffix) {
  m      <- pair$meta
  m      <- droplevels(m[!(m$Case_ID %in% m[bad_samples, "Case_ID"]), ])   # drop whole line(s)
  pair_s <- list(meta = m, counts = pair$counts[, rownames(m)])
  dds_s  <- run_deseq2_paired(pair_s, paste0(label, suffix), pipelines_dir)
  res_s  <- results(dds_s, contrast = c("Treatment", comparator, reference))

  common <- intersect(rownames(res_full), rownames(res_s))
  rho    <- cor(res_full[common, "stat"], res_s[common, "stat"],
                method = "spearman", use = "complete.obs")

  list(pairs_without = nrow(m) / 2,
       spearman      = rho,
       with          = length(sig_genes(res_full)),
       without       = length(sig_genes(res_s)),
       shared        = length(intersect(sig_genes(res_full), sig_genes(res_s))))
}


# -----------------------------------------------------------------------------
# 4.3 MA and volcano plots
# -----------------------------------------------------------------------------

plot_ma_simple <- function(res_obj, label, output_dir,
                           padj_thresh = PADJ_THRESH, fc_thresh = FC_THRESH) {

  ma_data <- as.data.frame(res_obj) %>%
    dplyr::filter(!is.na(padj), !is.na(log2FoldChange), baseMean > 0) %>%
    dplyr::mutate(status = dplyr::case_when(
      padj < padj_thresh & log2FoldChange >  fc_thresh ~ "UP",
      padj < padj_thresh & log2FoldChange < -fc_thresh ~ "DOWN",
      TRUE ~ "NS"))

  n_up   <- sum(ma_data$status == "UP")
  n_down <- sum(ma_data$status == "DOWN")

  p <- ggplot(ma_data, aes(x = log10(baseMean), y = log2FoldChange, color = status)) +
    geom_point(data = dplyr::filter(ma_data, status == "NS"), alpha = 0.3, size = 0.6) +
    geom_point(data = dplyr::filter(ma_data, status != "NS"), alpha = 0.7, size = 0.8) +
    scale_color_manual(values = c("NS" = "grey70", "UP" = "#D55E00", "DOWN" = "#0072B2"),
                       labels = c("NS"   = "Not sig.",
                                  "UP"   = paste0("UP (n=", n_up, ")"),
                                  "DOWN" = paste0("DOWN (n=", n_down, ")")),
                       name = NULL) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    labs(title    = paste("MA Plot —", label),
         subtitle = paste0("padj < ", padj_thresh, " & |log2FC| > ", fc_thresh),
         x = "Mean expression [log10(baseMean)]", y = "log2 Fold Change") +
    theme_classic(base_size = BASE_SIZE) +
    theme(plot.title      = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle   = element_text(hjust = 0.5, color = "grey40", size = 9),
          legend.position = "bottom")

  ggsave(paste0(output_dir, label, "_MAPlot_simple.png"), p, width = 8, height = 6, dpi = 400)
  print(p)
  p
}

# Direction-aware volcano: corner labels name the arm each side is higher in.
plot_volcano_styled <- function(res_obj, label, output_dir,
                                padj_thresh = PADJ_THRESH, fc_thresh = FC_THRESH,
                                n_label = 5) {

  arms     <- split_label(label)
  up_arm   <- arms[1]   # comparator: positive log2FC
  down_arm <- arms[2]   # reference:  negative log2FC

  res_df <- as.data.frame(res_obj) %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::mutate(diffexpressed = dplyr::case_when(
                    padj < padj_thresh & log2FoldChange >  fc_thresh ~ "UP",
                    padj < padj_thresh & log2FoldChange < -fc_thresh ~ "DOWN",
                    TRUE ~ "NS"),
                  gene = rownames(.))

  n_up   <- sum(res_df$diffexpressed == "UP")
  n_down <- sum(res_df$diffexpressed == "DOWN")

  # Label the n most significant genes on each side
  top_up      <- res_df %>% dplyr::filter(diffexpressed == "UP")   %>% dplyr::arrange(padj) %>% dplyr::slice_head(n = n_label)
  top_down    <- res_df %>% dplyr::filter(diffexpressed == "DOWN") %>% dplyr::arrange(padj) %>% dplyr::slice_head(n = n_label)
  label_genes <- dplyr::bind_rows(top_up, top_down)

  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(data = dplyr::filter(res_df, diffexpressed == "NS"),   color = "grey70",  alpha = 0.3, size = 0.8) +
    geom_point(data = dplyr::filter(res_df, diffexpressed == "DOWN"), color = "#2E7D32", alpha = 0.8, size = 1.4) +
    geom_point(data = dplyr::filter(res_df, diffexpressed == "UP"),   color = "#C62828", alpha = 0.8, size = 1.4) +
    geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed", color = "grey80", linewidth = 0.4) +
    geom_vline(xintercept = c(-fc_thresh, fc_thresh), linetype = "dashed", color = "grey80", linewidth = 0.4) +
    # Counts anchored to the top corners of the panel
    annotate("text", x = -Inf, y = Inf, label = paste0("Higher in ", down_arm, ": ", n_down),
             color = "#2E7D32", fontface = "bold", hjust = -0.05, vjust = 1.5, size = 5.5) +
    annotate("text", x =  Inf, y = Inf, label = paste0("Higher in ", up_arm, ": ", n_up),
             color = "#C62828", fontface = "bold", hjust = 1.05, vjust = 1.5, size = 5.5) +
    geom_label_repel(data = label_genes, aes(label = gene, color = diffexpressed),
                     size = 4.2, fontface = "italic", fill = "white",
                     max.overlaps = 40, segment.size = 0.3, segment.color = "grey40",
                     box.padding = 0.4, label.padding = unit(0.2, "lines"),
                     show.legend = FALSE) +
    scale_color_manual(values = c("UP" = "#C62828", "DOWN" = "#2E7D32")) +
    xlim(-25, 25) +   # fixed symmetric x-range; genes with |log2FC| > 25 are not drawn
    labs(title    = paste("Volcano Plot —", gsub("_", " ", label)),
         subtitle = paste0("padj < ", padj_thresh, " & |log2FC| > ", fc_thresh),
         x = "log2(Fold Change)", y = "-log10(padj)") +
    theme_minimal() +
    theme(plot.title       = element_text(hjust = 0.5, face = "bold", size = 18),
          plot.subtitle    = element_text(hjust = 0.5, color = "grey40", size = 13),
          axis.title.x     = element_text(size = 14, margin = margin(t = 10)),
          axis.title.y     = element_text(size = 14, margin = margin(r = 10)),
          axis.text        = element_text(size = 12),
          axis.line        = element_line(color = "grey80"),
          panel.grid.major = element_line(color = "#f5f5f5"),
          panel.grid.minor = element_blank(),
          legend.position  = "none")

  ggsave(paste0(output_dir, label, "_VolcanoPlot_styled.png"), p, width = 8, height = 8, dpi = 400)
  cat("Saved styled volcano plot for:", label, "\n")
  p
}


# -----------------------------------------------------------------------------
# 4.4 Pre-ranked GSEA
# -----------------------------------------------------------------------------

# Rank genes by the DESeq2 Wald statistic (sign = direction, size = evidence).
build_ranked_list <- function(res_obj) {
  df <- as.data.frame(res_obj)
  df <- df[!is.na(df$stat), ]
  sort(setNames(df$stat, rownames(df)), decreasing = TRUE)
}

# Run GSEA for one comparison x collection; save all / significant tables and
# a single-colour bar plot. pvalueCutoff = 1 returns every tested set.
run_gsea_and_plot <- function(gene_list, t2g, set_name, comparison_label, output_dir,
                              top_n = NULL, padj_thresh = PADJ_THRESH,
                              n_perm_simple = N_PERM_SIMPLE) {

  gsea_res <- GSEA(gene_list,
                   exponent = 1, minGSSize = 10, maxGSSize = 10000,
                   pvalueCutoff = 1, pAdjustMethod = "BH",
                   TERM2GENE = t2g, verbose = TRUE, seed = TRUE,
                   nPermSimple = n_perm_simple)

  df        <- as.data.frame(gsea_res)
  file_stub <- paste0(comparison_label, "_", set_name, "_GSEA")

  # ---- All tested pathways ----
  write.csv(df, paste0(output_dir, file_stub, "_RESULTS.csv"))
  cat(sprintf("[%s | %s] Pathways tested: %d | Significant (padj<%.2f): %d\n",
              comparison_label, set_name, nrow(df), padj_thresh, sum(df$p.adjust < padj_thresh)))

  # ---- Significant pathways, ordered by NES (high -> low) ----
  sig <- df[df$p.adjust < padj_thresh, ]
  sig <- sig[order(-sig$NES), ]
  write.csv(sig, paste0(output_dir, file_stub, "_Sorted.csv"))

  if (nrow(sig) == 0) {
    message(sprintf("[%s | %s] No significant pathways to plot.", comparison_label, set_name))
    return(list(all = df, sig = sig, plot = NULL))
  }

  # ---- Bar plot (top_n by NES, or all significant) ----
  plot_data    <- if (!is.null(top_n)) head(sig, top_n) else sig
  title_suffix <- if (!is.null(top_n)) paste0("Top ", top_n, " ") else "Significant "

  p <- ggplot(plot_data, aes(reorder(Description, NES), NES)) +
    geom_col(fill = BAR_FILL, width = 0.7) +
    coord_flip() +
    labs(x = "Pathway", y = "Normalized Enrichment Score",
         title = paste0(title_suffix, set_name, " Pathways: ", comparison_label)) +
    theme_classic(base_size = BASE_SIZE) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(paste0(output_dir, file_stub, "_barplot.png"), p,
         width = 10, height = max(5, nrow(plot_data) * 0.3), dpi = 400)   # height scales with n
  print(p)

  list(all = df, sig = sig, plot = p)
}

# Display names: drop collection prefix, underscores -> spaces, UPPER CASE, wrap.
clean_pathway_label <- function(x) {
  x %>%
    str_remove("^(GOBP_|GOCC_|GOMF_|HP_|KEGG_|HALLMARK_)") %>%
    str_replace_all("_", " ") %>%
    str_to_upper() %>%
    str_wrap(width = 30)
}

# One panel of the combined figure: top positive and negative NES pathways,
# green = positive (higher in comparator), red = negative (higher in reference).
build_gsea_panel <- function(sig_df, panel_title, top_n = 12) {

  if (nrow(sig_df) == 0) {
    return(ggplot() + theme_void() + labs(title = panel_title, subtitle = "No significant pathways"))
  }

  top_pos <- sig_df %>% dplyr::filter(NES > 0) %>% dplyr::arrange(dplyr::desc(NES)) %>% dplyr::slice_head(n = ceiling(top_n / 2))
  top_neg <- sig_df %>% dplyr::filter(NES < 0) %>% dplyr::arrange(NES)              %>% dplyr::slice_head(n = floor(top_n / 2))

  plot_data <- dplyr::bind_rows(top_pos, top_neg) %>%
    dplyr::mutate(label_clean = clean_pathway_label(Description),
                  label_clean = factor(label_clean,
                                       levels = rev(clean_pathway_label(c(top_pos$Description,
                                                                          rev(top_neg$Description))))),
                  direction   = ifelse(NES > 0, "Positive", "Negative"))

  ggplot(plot_data, aes(label_clean, NES, fill = direction)) +
    geom_col(width = 0.65) +
    coord_flip() +
    scale_fill_manual(values = c("Negative" = "#B22222", "Positive" = "#008000"), guide = "none") +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
    labs(x = NULL, y = "Normalized Enrichment Score", title = panel_title) +
    theme_classic(base_size = 12) +
    theme(plot.title   = element_text(hjust = 0.5, face = "bold", size = 14),
          axis.text.y  = element_text(size = 11, face = "bold", color = "black", lineheight = 0.8),
          axis.title.x = element_text(size = 11, face = "bold"),
          axis.text.x  = element_text(size = 10),
          plot.margin  = margin(5, 10, 5, 5))
}

# One collection, all three comparisons side by side.
plot_gsea_combined <- function(pathway_results, set_name, output_dir, top_n = 12) {

  panels <- lapply(names(pathway_results), function(lbl) {
    build_gsea_panel(pathway_results[[lbl]][[set_name]]$sig, gsub("_", " ", lbl), top_n)
  })

  combined <- wrap_plots(panels, nrow = 1, widths = c(1, 1, 1)) +
    plot_annotation(
      title    = paste0("GSEA — significant ", set_name, " pathways (padj < ", PADJ_THRESH, ")"),
      subtitle = "Ordered by NES across all timepoints",
      theme    = theme(plot.title    = element_text(hjust = 0.5, face = "bold", size = 15),
                       plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey30"))
    )

  ggsave(paste0(output_dir, set_name, "_GSEA_combined_barplot.png"), combined,
         width = 20, height = 8, dpi = 400)
  print(combined)
  combined
}


# -----------------------------------------------------------------------------
# 4.5 Over-representation analysis (ORA)
# -----------------------------------------------------------------------------

# Hypergeometric test of DEGs (padj + |log2FC| thresholds) against the
# universe of genes that received an adjusted p-value.
run_ora <- function(res_obj, t2g, set_name, comparison_label, output_dir,
                    padj_thresh = PADJ_THRESH, fc_thresh = FC_THRESH) {

  deg      <- rownames(res_obj)[which(res_obj$padj < padj_thresh &
                                        abs(res_obj$log2FoldChange) > fc_thresh)]
  universe <- rownames(res_obj)[!is.na(res_obj$padj)]

  ora_res <- enricher(gene = deg, universe = universe, TERM2GENE = t2g,
                      pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH")

  df <- as.data.frame(ora_res)
  write.csv(df, paste0(output_dir, comparison_label, "_", set_name, "_ORA_RESULTS.csv"))

  cat(sprintf("[%s | %s ORA] DEGs tested: %d | Significant (padj<%.2f): %d\n",
              comparison_label, set_name, length(deg), padj_thresh, sum(df$p.adjust < padj_thresh)))
  df
}


# -----------------------------------------------------------------------------
# 4.6 ssGSEA
# -----------------------------------------------------------------------------

# Per-sample pathway scores from VST expression (default ssGSEA settings).
run_ssgsea <- function(dds, t2g, set_name, comparison_label, output_dir) {
  expr_mat  <- assay(vst(dds, blind = FALSE))
  gene_sets <- split(t2g$gene_symbol, t2g$gs_name)

  gsvapar <- ssgseaParam(expr_mat, gene_sets, minSize = 10, maxSize = 10000)
  scores  <- gsva(gsvapar)   # pathways x samples

  write.csv(scores, paste0(output_dir, comparison_label, "_", set_name, "_ssGSEA_scores.csv"))
  cat(sprintf("[%s | %s ssGSEA] Pathways scored: %d | Samples: %d\n",
              comparison_label, set_name, nrow(scores), ncol(scores)))
  scores
}


# -----------------------------------------------------------------------------
# 4.7 Transcription factor activity (DoRothEA + decoupleR ULM)
# -----------------------------------------------------------------------------

# ULM regresses the gene-level Wald statistic on each TF's signed regulon.
# Positive score = TF more active in the comparator arm.
run_tf_activity <- function(res_obj, net, comparison_label, output_dir, minsize = 5) {
  df  <- as.data.frame(res_obj)
  df  <- df[!is.na(df$stat), , drop = FALSE]
  mat <- matrix(df$stat, ncol = 1, dimnames = list(rownames(df), "stat"))

  tf_act <- run_ulm(mat = mat, net = net, .source = "source",
                    .target = "target", .mor = "mor", minsize = minsize)

  write.csv(tf_act, paste0(output_dir, comparison_label, "_TF_activity_DoRothEA.csv"))
  cat(sprintf("[%s] TFs tested: %d\n", comparison_label, nrow(tf_act)))
  tf_act
}

# Single-comparison bar plot: top_n TFs by |score| among raw p < threshold.
plot_tf_activity <- function(tf_act, comparison_label, output_dir,
                             top_n = 20, padj_thresh = PADJ_THRESH) {

  arms <- split_label(comparison_label)
  df   <- tf_act %>% dplyr::filter(p_value < padj_thresh) %>% dplyr::arrange(dplyr::desc(abs(score)))

  cat(sprintf("[%s TF activity] TFs tested: %d | Significant (p<%.2f): %d\n",
              comparison_label, nrow(tf_act), padj_thresh, nrow(df)))

  if (nrow(df) == 0) {
    message(sprintf("[%s] No significant TFs to plot.", comparison_label))
    return(NULL)
  }

  plot_data <- head(df, top_n) %>% dplyr::mutate(direction = ifelse(score > 0, "UP", "DOWN"))
  n_up      <- sum(plot_data$direction == "UP")
  n_down    <- sum(plot_data$direction == "DOWN")

  p <- ggplot(plot_data, aes(reorder(source, score), score, fill = direction)) +
    geom_col(width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c("UP" = "#D55E00", "DOWN" = "#0072B2"),
                      labels = c("UP"   = paste0("Higher in ", arms[1], " (n=", n_up, ")"),
                                 "DOWN" = paste0("Higher in ", arms[2], " (n=", n_down, ")")),
                      name = NULL) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    labs(x = "Transcription factor", y = "Activity score (ULM)",
         title    = paste0("Top ", nrow(plot_data), " TF Activities: ", gsub("_", " ", comparison_label)),
         subtitle = paste0("p < ", padj_thresh, " | positive = higher in ", arms[1])) +
    theme_classic(base_size = BASE_SIZE) +
    theme(plot.title      = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle   = element_text(hjust = 0.5, color = "grey40", size = 9),
          legend.position = "bottom")

  ggsave(paste0(output_dir, comparison_label, "_TF_activity_barplot.png"), p,
         width = 8, height = max(5, nrow(plot_data) * 0.3), dpi = 400)
  print(p)
  p
}

# One panel of the combined TF figure: top_n positive + top_n negative scores.
build_tf_panel <- function(tf_act, comparison_label, top_n = 10, padj_thresh = PADJ_THRESH) {

  arms       <- split_label(comparison_label)
  up_label   <- arms[1]
  down_label <- arms[2]

  sig      <- tf_act %>% dplyr::filter(p_value < padj_thresh)
  top_up   <- sig %>% dplyr::filter(score > 0) %>% dplyr::arrange(dplyr::desc(score)) %>% dplyr::slice_head(n = top_n)
  top_down <- sig %>% dplyr::filter(score < 0) %>% dplyr::arrange(score)              %>% dplyr::slice_head(n = top_n)

  plot_data <- dplyr::bind_rows(top_up, top_down) %>%
    dplyr::mutate(direction = ifelse(score > 0, up_label, down_label),
                  source    = factor(source, levels = rev(c(top_up$source, rev(top_down$source)))))

  ggplot(plot_data, aes(source, score, fill = direction)) +
    geom_col(width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = setNames(c(TF_UP_COLOR, TF_DOWN_COLOR), c(up_label, down_label)),
                      name = NULL) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
    labs(x = NULL, y = "Activity score (ULM)", title = gsub("_", " ", comparison_label)) +
    theme_classic(base_size = BASE_SIZE) +
    theme(plot.title      = element_text(hjust = 0.5, face = "bold", size = 11),
          axis.text.y     = element_text(size = 10, face = "bold"),
          legend.position = "none")
}

# All three comparisons side by side.
plot_tf_combined <- function(tf_results, output_dir, top_n = 10, padj_thresh = PADJ_THRESH) {

  panels <- lapply(names(tf_results), function(lbl) {
    build_tf_panel(tf_results[[lbl]], lbl, top_n, padj_thresh)
  })

  combined <- wrap_plots(panels, nrow = 1) +
    plot_annotation(
      title    = "TF activity (ULM) — DoRothEA A/B/C",
      subtitle = "Positive score = higher activity in the first-named arm",
      theme    = theme(plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
                       plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30"))
    )

  ggsave(paste0(output_dir, "TF_activity_combined_barplot.png"), combined,
         width = 14, height = 6, dpi = 400)
  print(combined)
  combined
}


#===============================================================================
# 5. STEP 1 — MATCHED PAIRS, PAIRED DESeq2, RESULTS (all comparisons)
#===============================================================================

pair_list <- list()
dds_list  <- list()
res_list  <- list()

for (lbl in names(COMPARISONS)) {
  comp <- COMPARISONS[[lbl]]

  # ---- Matched pairs (level_A = reference, level_B = comparator) ----
  pair_list[[lbl]] <- build_matched_pairs(TRT, GEX, level_A = comp[["reference"]],
                                          level_B = comp[["comparator"]])
  cat(lbl, "— samples:", nrow(pair_list[[lbl]]$meta), "\n")
  print(table(pair_list[[lbl]]$meta$Treatment))

  # ---- Fit and extract ----
  dds_list[[lbl]] <- run_deseq2_paired(pair_list[[lbl]], lbl, PIPELINES_DIR)
  res_list[[lbl]] <- extract_results(dds_list[[lbl]], comp[["comparator"]], comp[["reference"]],
                                     PIPELINES_DIR, lbl)
}


#===============================================================================
# 6. STEP 2 — PCA AND OUTLIER SCREEN (all comparisons)
#===============================================================================

pca_results <- lapply(names(dds_list), function(lbl) plot_pca_labeled(dds_list[[lbl]], lbl, OUTPUT_DIR))
names(pca_results) <- names(dds_list)

for (lbl in names(pca_results)) {
  cat("\n=====", lbl, "=====\n")
  detect_pca_outliers(pca_results[[lbl]]$pca_data, dds_list[[lbl]])
}


#===============================================================================
# 7. STEP 3 — OUTLIER SENSITIVITY ANALYSIS
#    Does the flagged low-depth sample (OUTLIER_SAMPLE) drive the results?
#    For each affected comparison: Cook's outlier count, library-size rank,
#    and a re-fit without its PDX line.
#===============================================================================

# ---- 7.1 24h vs Placebo ----
cat("\n===== Sensitivity: 24h_vs_Placebo =====\n")
print(TRT[TRT$Case_ID == TRT[OUTLIER_SAMPLE, "Case_ID"], c("Case_ID", "Treatment")])   # arms for this line
print(sort(colSums(counts(dds_list[["24h_vs_Placebo"]])))[1:5])                       # is it the only very-low library?
cat("Cook's outlier genes for", OUTLIER_SAMPLE, ":",
    count_cooks_outliers(dds_list[["24h_vs_Placebo"]], OUTLIER_SAMPLE), "\n")

print(sensitivity_check(pair_list[["24h_vs_Placebo"]], OUTLIER_SAMPLE, "24h", "Placebo",
                        res_list[["24h_vs_Placebo"]], "24h_vs_Placebo", PIPELINES_DIR,
                        suffix = "_noCRC0081"))

# ---- 7.2 Is the sample also in the other comparisons? ----
cat("In 6wks_vs_Placebo:", OUTLIER_SAMPLE %in% colnames(dds_list[["6wks_vs_Placebo"]]), "\n")
cat("In 6wks_vs_24h:",     OUTLIER_SAMPLE %in% colnames(dds_list[["6wks_vs_24h"]]), "\n")

# ---- 7.3 6wks vs Placebo ----
cat("\n===== Sensitivity: 6wks_vs_Placebo =====\n")
print(sensitivity_check(pair_list[["6wks_vs_Placebo"]], OUTLIER_SAMPLE, "6wks", "Placebo",
                        res_list[["6wks_vs_Placebo"]], "6wks_vs_Placebo", PIPELINES_DIR,
                        suffix = "_noCRC0081"))

dds_6p <- dds_list[["6wks_vs_Placebo"]]
cat("Cook's outlier genes for", OUTLIER_SAMPLE, ":", count_cooks_outliers(dds_6p, OUTLIER_SAMPLE), "\n")
print(sort(colSums(counts(dds_6p)))[1:3])

# ---- 7.4 6wks vs Placebo — drop the three lowest-depth samples + outlier ----
low_depth <- names(sort(colSums(counts(dds_6p)))[1:3])
print(TRT[low_depth, c("Case_ID", "Treatment")])                 # which lines and arms
cat("Median library size:", median(colSums(counts(dds_6p))), "\n")

print(sensitivity_check(pair_list[["6wks_vs_Placebo"]], c(low_depth, OUTLIER_SAMPLE), "6wks", "Placebo",
                        res_list[["6wks_vs_Placebo"]], "6wks_vs_Placebo_lowdepth", PIPELINES_DIR,
                        suffix = "_sens"))

# ---- 7.5 Library size across the whole dataset ----
lib <- colSums(GEX)
print(sort(lib)[1:10])                                             # 10 smallest libraries
print(sapply(c(100e3, 250e3, 500e3, 1e6), function(t) sum(lib < t)))   # samples below each cut-off
print(table(TRT[names(lib)[lib < 500e3], "Treatment"]))            # arm of samples < 500k reads


#===============================================================================
# 8. STEP 4 — MA AND VOLCANO PLOTS (all comparisons)
#===============================================================================

for (lbl in names(res_list)) {
  plot_ma_simple(res_list[[lbl]], lbl, OUTPUT_DIR)
  plot_volcano_styled(res_list[[lbl]], lbl, OUTPUT_DIR)
}


#===============================================================================
# 9. STEP 5 — PRE-RANKED GSEA (all comparisons x collections)
#===============================================================================

pathway_results <- list()

for (lbl in names(res_list)) {
  gene_list <- build_ranked_list(res_list[[lbl]])
  cat(sprintf("\n[%s] Genes in ranked list: %d\n", lbl, length(gene_list)))

  pathway_results[[lbl]] <- lapply(names(GENESET_COLLECTIONS), function(set_name) {
    coll <- GENESET_COLLECTIONS[[set_name]]
    run_gsea_and_plot(gene_list, coll$t2g, set_name, lbl, OUTPUT_DIR, top_n = coll$gsea_top_n)
  })
  names(pathway_results[[lbl]]) <- names(GENESET_COLLECTIONS)
}

# ---- Combined figures: one per collection, three comparisons side by side ----
gsea_combined_plots <- lapply(names(GENESET_COLLECTIONS), function(set_name) {
  plot_gsea_combined(pathway_results, set_name, OUTPUT_DIR,
                     top_n = GENESET_COLLECTIONS[[set_name]]$panel_top_n)
})
names(gsea_combined_plots) <- names(GENESET_COLLECTIONS)


#===============================================================================
# 10. STEP 6 — OVER-REPRESENTATION ANALYSIS (all comparisons x collections)
#===============================================================================

ora_results <- list()

for (lbl in names(res_list)) {
  ora_results[[lbl]] <- lapply(names(GENESET_COLLECTIONS), function(set_name) {
    run_ora(res_list[[lbl]], GENESET_COLLECTIONS[[set_name]]$t2g, set_name, lbl, OUTPUT_DIR)
  })
  names(ora_results[[lbl]]) <- names(GENESET_COLLECTIONS)
}


#===============================================================================
# 11. STEP 7 — ssGSEA HALLMARK SCORES (all comparisons)
#===============================================================================

ssgsea_results <- list()

for (lbl in names(dds_list)) {
  ssgsea_results[[lbl]] <- list(
    Hallmark = run_ssgsea(dds_list[[lbl]], h_t2g, "Hallmark", lbl, OUTPUT_DIR)
  )
}


#===============================================================================
# 12. STEP 8 — TF ACTIVITY (all comparisons)
#===============================================================================

# ---- 12.1 Infer activities ----
tf_results <- list()
for (lbl in names(res_list)) {
  tf_results[[lbl]] <- run_tf_activity(res_list[[lbl]], DOROTHEA_NET, lbl, OUTPUT_DIR)
}

print(tf_results[["24h_vs_Placebo"]] %>% dplyr::arrange(p_value) %>% head(15))   # quick look

# ---- 12.2 Single-comparison bar plots ----
tf_plots <- list()
for (lbl in names(tf_results)) {
  tf_plots[[lbl]] <- plot_tf_activity(tf_results[[lbl]], lbl, OUTPUT_DIR)
}

# ---- 12.3 Combined figure ----
tf_combined_plot <- plot_tf_combined(tf_results, OUTPUT_DIR)


#===============================================================================
# 13. REPRODUCIBILITY
#===============================================================================

sessionInfo()
