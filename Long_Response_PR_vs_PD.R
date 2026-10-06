################################################################################
#
#  Script      : Long_Response_PR_vs_PD.R
#  Project     : Transcriptomic profiling of FOLFIRI / 5-FU response in
#                colorectal cancer (CRC) patient-derived xenograft (PDX) models
#  Analysis    : Longitudinal response — Partial Response (PR) vs
#                Progressive Disease (PD), run separately at three arms /
#                timepoints: Placebo, 24 h and 6 weeks
#  Author      : Janmita Kaverimane Umesh (MSc Bioinformatics, QUB)
#
#  Workflow (identical order for every comparison)
#  ---------------------------------------------------------------------------
#    1. Subset samples to one timepoint and to PR / PD labels only
#    2. Differential expression  — DESeq2 (~ Response, PD = reference)
#    3. Sample-level QC          — Cook's distance, VST + PCA
#    4. Volcano plot             — padj < 0.05 & |log2FC| > 0.585
#    5. Pre-ranked GSEA          — clusterProfiler, MSigDB H / KEGG / GO
#    6. ssGSEA                   — GSVA per-sample scores + Wilcoxon PR vs PD
#    7. TF activity              — DoRothEA (A/B/C) regulons + decoupleR ULM
#
#  Inputs  (0_data/PDX_Longitudinal/counts/)
#    - CRC_graft_raw_counts_matrix_12022024.csv  : raw counts, genes x samples
#    - Complete_GraftCRC_colData_19022024.csv    : sample metadata
#
#  Outputs
#    - 2_pipelines/.../LONG_RESPONSE/           : fitted DESeq2 objects (.rds)
#                                                 + full DE result tables
#    - 3_output/.../LONG_RESPONSE/<comparison>/ : plots and enrichment tables
#
#  Notes
#    - All thresholds and parameters are carried over unchanged from the
#      original per-timepoint scripts. Where those scripts differed, the
#      difference is kept and declared explicitly in Section 2 (COMPARISONS).
#    - Positive log2FC / NES / delta / TF score = higher in PR than in PD.
#
################################################################################


#===============================================================================
# 0. SESSION SETUP
#===============================================================================

rm(list = ls(all.names = TRUE))     # start from a clean environment
gc()                                # release memory from previous sessions

# ---- 0.1 Libraries ----
# Core DE
library(DESeq2)          # differential expression (NB GLM); glmGamPoi must be installed
# Enrichment
library(clusterProfiler) # pre-ranked GSEA
library(msigdbr)         # MSigDB gene-set collections
library(GSVA)            # single-sample GSEA (ssGSEA)
library(dorothea)        # TF -> target regulons
library(decoupleR)       # TF activity inference (ULM)
# Data handling
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
# Plotting
library(ggplot2)
library(ggrepel)         # non-overlapping gene labels on volcano plots
library(patchwork)       # multi-panel figure assembly

# ---- 0.2 Project paths ----
setwd("C:/Users/40496110/Downloads/PDX_PROJECT")

PIPELINES_DIR <- "./2_pipelines/PDX_Longitudinal/LONG_RESPONSE/"
OUTPUT_ROOT   <- "./3_output/PDX_Longitudinal/LONG_RESPONSE/"
COUNTS_DIR    <- "./0_data/PDX_Longitudinal/counts/"
ANNOT_DIR     <- "./0_data/PDX_Longitudinal/annotations/"   # reserved, not used here

COUNTS_FILE   <- "CRC_graft_raw_counts_matrix_12022024.csv"
COLDATA_FILE  <- "Complete_GraftCRC_colData_19022024.csv"

# ---- 0.3 Global analysis constants (shared by all three comparisons) ----
RESPONSE_COL   <- "Response_3_6wk"  # response label measured at 3-6 weeks
MIN_TOTAL_CNT  <- 10                # pre-filter: keep genes with >= 10 reads in total
PADJ_CUTOFF    <- 0.05              # BH-adjusted significance threshold
LFC_CUTOFF     <- 0.585             # |log2FC| > 0.585  ==  fold change > 1.5
N_TOP_PATHWAYS <- 5                 # pathways shown per direction in bar plots
N_TOP_TFS      <- 10                # TFs shown per direction in TF bar plot
SEED           <- 123


#===============================================================================
# 1. COMPARISON CONFIGURATION
#    One entry per timepoint. Everything that differed between the original
#    scripts is declared here so the differences are visible in one place.
#===============================================================================

COMPARISONS <- list(

  # ---------------------------------------------------------------- Placebo ---
  placebo = list(
    treatment   = "Placebo",                       # value in coldata$Treatment
    tag         = "placebo",                       # used in the .rds filename
    prefix      = "Placebo_RES",                   # prefix for every output file
    title_label = "Placebo (PR vs PD)",
    short_label = "Placebo",
    output_dir  = paste0(OUTPUT_ROOT, "long_placebo_pr_pd/"),
    de_csv      = "Placebo_PR_vs_PD_DESeq2_results.csv",

    pca  = list(title = "PCA: Placebo PR vs PD", style = "bold", save = FALSE),

    gsea = list(n_perm   = 100000,                 # nPermSimple
                eps      = 0,                      # exact small p-values
                sig_mode = "recompute_padj",       # adds a BH 'padj' column
                top_mode = "slice",                # top 5 by NES per direction
                individual_plots = FALSE,
                bar_theme = "y13",
                combined_title    = "GSEA — Significant Pathways: Placebo (PR vs PD)",
                combined_subtitle = "Top 5 activated & top 5 suppressed by NES | p.adj < 0.05",
                annot_sizes = c(title = 18, subtitle = 13),
                width = 18, height = 13),

    ssgsea = list(mode  = "padj_go_prefilter",     # GO-only terms selected BEFORE BH
                  title = "ssGSEA — Significant Pathways: Placebo (PR vs PD)"),

    tf = list(mode = "padj_with_fallback",         # padj < 0.05, else raw p for plot
              subtitle = "DoRothEA A/B/C | Top TFs filtered by padj | positive score = more active in PR")
  ),

  # -------------------------------------------------------------------- 24h ---
  h24 = list(
    treatment   = "24h",
    tag         = "24h",
    prefix      = "24_RES",
    title_label = "24 Hours (PR vs PD)",
    short_label = "24h",
    output_dir  = paste0(OUTPUT_ROOT, "long_24_pr_pd/"),
    de_csv      = "24h_PR_vs_PD_DESeq2_results.csv",

    pca  = list(title = "PCA: 24h PR vs PD", style = "plain", save = TRUE),

    gsea = list(n_perm   = 10000,
                eps      = NULL,                   # NULL = clusterProfiler default
                sig_mode = "nes_desc",             # filter p.adjust, order by NES
                top_mode = "head_tail",            # first/last 5 of NES-sorted list
                individual_plots = TRUE,           # also one bar plot per collection
                bar_theme = "large24",
                combined_title    = "GSEA — significant pathways: 24h (PR vs PD)",
                combined_subtitle = "ordered by NES, p.adj < 0.05",
                annot_sizes = c(title = 20, subtitle = 14),
                width = 18, height = 11),

    ssgsea = list(mode  = "rawp_individual",       # plots use raw Wilcoxon p < 0.05
                  title = "ssGSEA Differential Enrichment Pathway: 24h"),

    tf = list(mode = "padj_strict",                # plot only if any TF padj < 0.05
              subtitle = "DoRothEA A/B/C | positive score = more active in PR | padj < 0.05")
  ),

  # ------------------------------------------------------------------ 6 wks ---
  w6 = list(
    treatment   = "6wks",
    tag         = "6w",
    prefix      = "6w_RES",
    title_label = "6 Weeks (PR vs PD)",
    short_label = "6wks",
    output_dir  = paste0(OUTPUT_ROOT, "long__6w_pr_pd/"),   # double '_' kept as in original
    de_csv      = "6wks_PR_vs_PD_DESeq2_results.csv",

    pca  = list(title = "PCA: 6wks PR vs PD", style = "bold", save = FALSE),

    gsea = list(n_perm   = 100000,
                eps      = 0,
                sig_mode = "padj",                 # uses clusterProfiler p.adjust
                top_mode = "slice",
                individual_plots = FALSE,
                bar_theme = "y12",
                combined_title    = "GSEA — Significant Pathways: 6 Weeks (PR vs PD)",
                combined_subtitle = "Top 5 significant activated & suppressed pathways/terms (p.adj < 0.05), ordered by NES",
                annot_sizes = c(title = 18, subtitle = 13),
                width = 18, height = 11),

    ssgsea = list(mode  = "padj_go_postfilter",    # BH over all C5, then GO-only kept
                  title = "ssGSEA — Significant Pathways: 6 Weeks (PR vs PD)"),

    tf = list(mode = "raw_p",                      # raw p < 0.05 (no BH)
              subtitle = "DoRothEA A/B/C | positive score = more active in PR")
  )
)


#===============================================================================
# 2. LOAD AND ALIGN SHARED INPUTS
#    Counts and metadata are identical for all comparisons, so load once.
#===============================================================================

# ---- 2.1 Raw counts (genes x samples) ----
GEX <- read.csv(paste0(COUNTS_DIR, COUNTS_FILE), row.names = 1, check.names = FALSE)

# ---- 2.2 Sample metadata (samples x variables) ----
TRT <- read.csv(paste0(COUNTS_DIR, COLDATA_FILE), row.names = 1, check.names = FALSE)

cat("GEX dimensions:", dim(GEX), "\n")
cat("TRT dimensions:", dim(TRT), "\n")

# ---- 2.3 Keep samples present in both tables; order counts to match metadata ----
# DESeq2 requires colnames(counts) to be identical to, and in the same order as,
# rownames(colData).
TRT <- TRT[rownames(TRT) %in% colnames(GEX), ]
GEX <- GEX[, rownames(TRT)]
cat("TRT after filtering:", nrow(TRT), "\n")


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

cat("Hallmark sets:", length(unique(h_t2g$gs_name)), "\n")
cat("KEGG sets:",     length(unique(kegg_t2g$gs_name)), "\n")
cat("GO sets:",       length(unique(c5_t2g$gs_name)), "\n")

# Per-collection settings used by GSEA and ssGSEA
#   file_tag    : name used in GSEA output files
#   wrap        : character width for wrapping pathway labels on plots
#   strip       : prefix removed from labels in 24h ssGSEA plots
#   indiv_width : width (in) of 24h single-collection GSEA bar plots
GENESET_COLLECTIONS <- list(
  Hallmark = list(t2g = h_t2g,    file_tag = "Hallmarks", wrap = 22,
                  strip = "^HALLMARK_", indiv_width = 11),
  KEGG     = list(t2g = kegg_t2g, file_tag = "KEGG",      wrap = 22,
                  strip = "^KEGG_",     indiv_width = 11),
  GO       = list(t2g = c5_t2g,   file_tag = "C5_GO",     wrap = 16,
                  strip = NULL,         indiv_width = 12)
)

# ---- 3.2 DoRothEA regulons (high / medium confidence: A, B, C) ----
data("dorothea_hs", package = "dorothea")

DOROTHEA_NET <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B", "C")) %>%
  dplyr::rename(source = tf) %>%
  dplyr::mutate(weight = as.numeric(mor))   # mode of regulation: +1 activates, -1 represses

cat("DoRothEA network rows (A/B/C confidence):", nrow(DOROTHEA_NET), "\n")


#===============================================================================
# 4. HELPER FUNCTIONS
#===============================================================================

# -----------------------------------------------------------------------------
# 4.1 Plot styling
# -----------------------------------------------------------------------------

# Bold theme for horizontal pathway bar plots. Sizes are passed in so that each
# comparison keeps the exact font sizes of its original script.
bold_bar_theme <- function(title, axis_title, x_text, y_text, legend_title,
                           legend_text, y_lineheight = NULL, subtitle = NULL) {
  th <- theme(
    plot.title   = element_text(hjust = 0.5, face = "bold", size = title),
    axis.title   = element_text(size = axis_title, face = "bold"),
    axis.text.x  = element_text(size = x_text, face = "bold"),
    axis.text.y  = element_text(size = y_text, face = "bold", lineheight = y_lineheight),
    legend.title = element_text(size = legend_title, face = "bold"),
    legend.text  = element_text(size = legend_text, face = "bold")
  )
  if (!is.null(subtitle)) {
    th <- th + theme(plot.subtitle = element_text(hjust = 0.5, size = subtitle,
                                                  face = "bold", colour = "black"))
  }
  th
}

BAR_THEMES <- list(
  y13      = bold_bar_theme(16, 13, 12, 13, 12, 11, y_lineheight = 0.85),  # Placebo GSEA
  y12      = bold_bar_theme(16, 13, 12, 12, 12, 11, y_lineheight = 0.85),  # 6w GSEA; Placebo/6w ssGSEA
  large24  = bold_bar_theme(18, 14, 14, 12, 13, 12),                       # 24h GSEA
  ss24     = bold_bar_theme(16, 13, 13, 13, 12, 11, subtitle = 13)         # 24h ssGSEA
)

# Tidy MSigDB names for display: drop collection prefix, underscores -> spaces,
# title case, and wrap long names over several lines.
clean_labels <- function(df, width = 22) {
  if (nrow(df) == 0) return(df)
  df$Description <- df$Description %>%
    str_remove("^(HALLMARK_|KEGG_|GOBP_|GOCC_|GOMF_)") %>%
    str_replace_all("_", " ") %>%
    str_to_title() %>%
    str_wrap(width = width)
  df
}

# Put a list of panels side by side with a shared legend and an overall title.
assemble_panels <- function(plot_list, title, subtitle, title_size, subtitle_size) {
  wrap_plots(plot_list, ncol = length(plot_list)) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title    = title,
      subtitle = subtitle,
      theme    = theme(
        plot.title    = element_text(size = title_size,    face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = subtitle_size, face = "bold", hjust = 0.5)
      )
    ) &
    theme(legend.position = "bottom")
}


# -----------------------------------------------------------------------------
# 4.2 Sample selection and differential expression
# -----------------------------------------------------------------------------

# Keep one timepoint, then only PR / PD (drops SD and missing labels).
# PD is the reference level, so positive log2FC = higher in PR.
subset_pr_pd <- function(TRT, treatment) {
  trt_tp <- TRT[TRT$Treatment == treatment, ]
  cat(treatment, "samples total:", nrow(trt_tp), "\n")
  print(table(trt_tp[[RESPONSE_COL]], useNA = "always"))

  trt_prpd <- trt_tp[trt_tp[[RESPONSE_COL]] %in% c("PR", "PD"), ]
  trt_prpd[[RESPONSE_COL]] <- factor(trt_prpd[[RESPONSE_COL]], levels = c("PD", "PR"))

  cat(treatment, "samples with PR/PD label:", nrow(trt_prpd), "\n")
  print(table(trt_prpd[[RESPONSE_COL]]))
  trt_prpd
}

# Build, pre-filter and fit the DESeq2 model.
# Unpaired design (~ Response): PR vs PD is a between-model comparison.
fit_deseq2 <- function(counts, coldata, cfg) {
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData   = coldata,
                                design    = ~ Response_3_6wk)

  # Remove near-empty genes: speeds fitting and reduces multiple-testing burden
  dds <- dds[rowSums(counts(dds)) >= MIN_TOTAL_CNT, ]
  cat("Genes going into DESeq2:", nrow(dds), "\n")
  cat("Samples:", ncol(dds), "\n")

  set.seed(SEED)
  dds <- DESeq(dds, fitType = "glmGamPoi")   # glmGamPoi: fast, stable dispersion fit

  # Save the fitted object straight away so it can be reloaded without refitting
  if (!dir.exists(PIPELINES_DIR)) dir.create(PIPELINES_DIR, recursive = TRUE)
  saveRDS(dds, paste0(PIPELINES_DIR, "dds_", cfg$tag, "_res_fitted.rds"))
  cat("Fitted object saved to disk.\n")
  print(resultsNames(dds))
  dds
}

# Mean Cook's distance per sample: higher values = sample drives more
# outlier calls. Printed to console only (no file output).
report_cooks <- function(dds, coldata) {
  cooks_per_sample <- apply(assays(dds)[["cooks"]], 2, mean, na.rm = TRUE)
  cooks_df <- data.frame(SampleID  = names(cooks_per_sample),
                         MeanCooks = cooks_per_sample,
                         Response  = coldata[names(cooks_per_sample), RESPONSE_COL])
  print(cooks_df[order(-cooks_df$MeanCooks), ])
}

# Variance-stabilising transform (design-aware, blind = FALSE) + PCA.
plot_pca <- function(dds, cfg) {
  vsd <- vst(dds, blind = FALSE)
  p <- plotPCA(vsd, intgroup = RESPONSE_COL) + labs(title = cfg$pca$title)

  p <- if (cfg$pca$style == "bold") {
    p + theme_classic(base_size = 14) +
      theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  } else {
    p + theme_classic()
  }
  print(p)

  if (cfg$pca$save) {
    ggsave(paste0(cfg$output_dir, cfg$prefix, "_PCA_", cfg$short_label, "_PR_vs_PD.png"),
           plot = p, width = 8, height = 6, dpi = 400)
    cat("Saved PCA plot to output directory.\n")
  }
  invisible(p)
}


# -----------------------------------------------------------------------------
# 4.3 Volcano plot
# -----------------------------------------------------------------------------

# Five gene classes:
#   UP / DOWN in PR       : padj < cutoff AND |log2FC| > cutoff
#   High Change, Low Sig  : large effect but not significant
#   Low Change, High Sig  : significant but small effect
#   NS                    : neither
make_volcano <- function(res, prefix, title_label, output_dir,
                         padj_thresh = PADJ_CUTOFF, fc_thresh = LFC_CUTOFF, n_label = 5) {

  # ---- Classify genes (rows with NA padj, e.g. independent-filtered, dropped) ----
  res_filtered <- as.data.frame(res) %>%
    dplyr::filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    dplyr::mutate(
      gene = rownames(.),
      diffexpressed = dplyr::case_when(
        padj <  padj_thresh & log2FoldChange >  fc_thresh       ~ "UP in PR",
        padj <  padj_thresh & log2FoldChange < -fc_thresh       ~ "DOWN in PR",
        padj >= padj_thresh & abs(log2FoldChange) >  fc_thresh  ~ "High Change, Low Sig",
        padj <  padj_thresh & abs(log2FoldChange) <= fc_thresh  ~ "Low Change, High Sig",
        TRUE ~ "NS"
      ),
      diffexpressed = factor(diffexpressed,
                             levels = c("UP in PR", "DOWN in PR", "High Change, Low Sig",
                                        "Low Change, High Sig", "NS"))
    )

  # ---- Counts of significant genes (annotated on the plot) ----
  total_up   <- sum(res_filtered$diffexpressed == "UP in PR")
  total_down <- sum(res_filtered$diffexpressed == "DOWN in PR")
  cat(prefix, "— UP:", total_up, "| DOWN:", total_down, "\n")

  # ---- Genes to label: most significant first, then largest |log2FC| ----
  top_up   <- res_filtered %>% dplyr::filter(diffexpressed == "UP in PR") %>%
    dplyr::arrange(padj, dplyr::desc(log2FoldChange)) %>% dplyr::slice_head(n = n_label)
  top_down <- res_filtered %>% dplyr::filter(diffexpressed == "DOWN in PR") %>%
    dplyr::arrange(padj, log2FoldChange) %>% dplyr::slice_head(n = n_label)
  top_genes <- dplyr::bind_rows(top_up, top_down)

  # ---- Symmetric x-axis and headroom on y ----
  ymax <- max(-log10(res_filtered$padj), na.rm = TRUE) + 1
  xlim <- ceiling(max(abs(res_filtered$log2FoldChange), na.rm = TRUE))

  # ---- Plot ----
  p <- ggplot(res_filtered,
              aes(x = log2FoldChange, y = -log10(padj),
                  colour = diffexpressed, size = diffexpressed, alpha = diffexpressed)) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major = element_line(colour = "grey95", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle    = element_text(colour = "grey40", size = 10, hjust = 0.5),
      legend.position  = "right",
      legend.title     = element_blank(),
      legend.text      = element_text(size = 11)
    ) +
    # threshold guide lines
    geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = c(-fc_thresh, fc_thresh), linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    geom_point() +
    # colour-blind-safe palette (Okabe-Ito orange / blue)
    scale_colour_manual(values = c("UP in PR"             = "#D55E00",
                                   "DOWN in PR"           = "#0072B2",
                                   "High Change, Low Sig" = "#B0BEC5",
                                   "Low Change, High Sig" = "#8E8E8E",
                                   "NS"                   = "grey85"), drop = FALSE) +
    scale_size_manual(values  = c("UP in PR" = 2, "DOWN in PR" = 2,
                                  "High Change, Low Sig" = 1.2, "Low Change, High Sig" = 1.2,
                                  "NS" = 0.8), drop = FALSE) +
    scale_alpha_manual(values = c("UP in PR" = 0.9, "DOWN in PR" = 0.9,
                                  "High Change, Low Sig" = 0.5, "Low Change, High Sig" = 0.5,
                                  "NS" = 0.3), drop = FALSE) +
    geom_text_repel(data = top_genes, aes(label = gene),
                    size = 3.2, fontface = "bold", colour = "black",
                    max.overlaps = 30, box.padding = 0.4, point.padding = 0.2,
                    segment.colour = "grey50", segment.size = 0.3,
                    show.legend = FALSE, inherit.aes = TRUE) +
    annotate("text", x =  xlim * 0.85, y = ymax * 0.97, label = paste0("UP: ", total_up),
             colour = "#D55E00", fontface = "bold", size = 5) +
    annotate("text", x = -xlim * 0.85, y = ymax * 0.97, label = paste0("DOWN: ", total_down),
             colour = "#0072B2", fontface = "bold", size = 5) +
    coord_cartesian(xlim = c(-xlim, xlim), ylim = c(0, ymax)) +
    guides(size = "none", alpha = "none",
           colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(title    = paste("Volcano Plot —", title_label),
         subtitle = paste0("padj < ", padj_thresh, " & |log2FC| > ", fc_thresh,
                           " | Top ", n_label, " UP/DOWN labelled"),
         x = "log2 Fold Change (PR / PD)",
         y = expression(-log[10]~"(adjusted p-value)"))

  ggsave(paste0(output_dir, prefix, "_VolcanoPlot.png"), p, width = 11, height = 8, dpi = 500)
  cat("Saved volcano plot:", prefix, "\n")
  invisible(p)
}


# -----------------------------------------------------------------------------
# 4.4 Pre-ranked GSEA
# -----------------------------------------------------------------------------

# Rank genes by the DESeq2 Wald statistic (sign = direction, size = evidence).
make_ranked_list <- function(res) {
  res_df <- as.data.frame(res)
  res_df <- res_df[!is.na(res_df$stat), ]
  gene_list <- sort(setNames(res_df$stat, rownames(res_df)), decreasing = TRUE)
  cat("Genes in ranked list:", length(gene_list), "\n")
  gene_list
}

# Run clusterProfiler::GSEA for one collection. pvalueCutoff = 1 returns all
# tested sets so that filtering is done explicitly downstream.
# eps = NULL leaves the clusterProfiler default untouched (24h script).
run_gsea <- function(gene_list, t2g, n_perm, eps = NULL) {
  args <- list(geneList = gene_list, exponent = 1, minGSSize = 10, maxGSSize = 10000,
               pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = t2g,
               verbose = TRUE, seed = TRUE, nPermSimple = n_perm)
  if (!is.null(eps)) args$eps <- eps
  as.data.frame(do.call(GSEA, args))
}

# Significant pathways, sorted. Three variants, matching the original scripts:
#   recompute_padj : add BH 'padj' column, filter, sort padj up then NES down
#   padj           : filter on clusterProfiler 'p.adjust', same sort
#   nes_desc       : filter on 'p.adjust', sort by NES only (high -> low)
select_sig_gsea <- function(df, mode, cutoff = PADJ_CUTOFF) {
  switch(mode,
    recompute_padj = df %>%
      dplyr::mutate(padj = p.adjust(pvalue, method = "BH")) %>%
      dplyr::filter(!is.na(padj), padj < cutoff) %>%
      dplyr::arrange(padj, dplyr::desc(NES)),
    padj = df %>%
      dplyr::filter(!is.na(p.adjust), p.adjust < cutoff) %>%
      dplyr::arrange(p.adjust, dplyr::desc(NES)),
    nes_desc = {
      sig <- df[df$p.adjust < cutoff, ]
      sig[order(-sig$NES), ]
    }
  )
}

# Pick pathways to plot and tag their direction.
#   slice     : top n positive NES + top n most negative NES
#   head_tail : first n and last n rows of the NES-sorted table
select_top_gsea <- function(sig_df, mode, n = N_TOP_PATHWAYS) {
  if (nrow(sig_df) == 0) return(sig_df)
  top <- switch(mode,
    slice = dplyr::bind_rows(
      sig_df %>% dplyr::filter(NES > 0) %>% dplyr::slice_max(NES, n = n, with_ties = FALSE),
      sig_df %>% dplyr::filter(NES < 0) %>% dplyr::slice_min(NES, n = n, with_ties = FALSE)
    ),
    head_tail = unique(rbind(head(sig_df, n), tail(sig_df, n)))  # unique(): avoid duplicates if < 2n
  )
  top$Direction <- ifelse(top$NES > 0, "Activated", "Suppressed")
  top
}

# Horizontal NES bar plot for one collection (used in the combined figure).
make_gsea_bar <- function(plot_df, title, wrap, bar_theme) {
  if (nrow(plot_df) == 0) return(NULL)
  plot_df <- clean_labels(plot_df, width = wrap)
  ggplot(plot_df, aes(reorder(Description, NES), NES, fill = Direction)) +
    geom_col(width = 0.7) + coord_flip() +
    scale_fill_manual(values = c("Activated" = "#2E7D32", "Suppressed" = "#C62828")) +
    labs(x = NULL, y = "Normalized Enrichment Score", title = title) +
    theme_classic(base_size = 14) + bar_theme
}

# Single-collection bar plot with raw pathway names (24h script only).
save_gsea_bar_individual <- function(plot_df, coll_name, coll, cfg) {
  p <- ggplot(plot_df, aes(reorder(Description, NES), NES, fill = Direction)) +
    geom_col(width = 0.7) + coord_flip() +
    scale_fill_manual(values = c("Activated" = "#2E7D32", "Suppressed" = "#C62828")) +
    labs(x = "Pathway", y = "Normalized Enrichment Score",
         title = paste0("Top 5 Positive & Negative Significant ", coll_name, " Pathways")) +
    theme_classic(base_size = 14) +
    theme(plot.title      = element_text(hjust = 0.5, face = "bold"),
          axis.title      = element_text(size = 12),
          legend.position = "right")
  ggsave(paste0(cfg$output_dir, cfg$prefix, "_", coll$file_tag, "_GSEA_barplot.png"),
         p, width = coll$indiv_width, height = 5, dpi = 400)
  cat("Saved", coll_name, "barplot.\n")
}

# Full GSEA block: run -> save all -> save significant -> plot.
run_gsea_block <- function(res, cfg) {
  g         <- cfg$gsea
  gene_list <- make_ranked_list(res)
  panels    <- list()

  for (coll_name in names(GENESET_COLLECTIONS)) {
    coll <- GENESET_COLLECTIONS[[coll_name]]

    # ---- Run and save every tested pathway ----
    gsea_df <- run_gsea(gene_list, coll$t2g, n_perm = g$n_perm, eps = g$eps)
    write.csv(gsea_df, paste0(cfg$output_dir, cfg$prefix, "_", coll$file_tag, "_GSEA_RESULTS.csv"))

    # ---- Significant pathways, sorted ----
    sig_df <- select_sig_gsea(gsea_df, g$sig_mode)
    write.csv(sig_df, paste0(cfg$output_dir, cfg$prefix, "_", coll$file_tag, "_GSEA_Sorted.csv"),
              row.names = (g$sig_mode == "nes_desc"))   # 24h tables kept row names
    cat(coll_name, "tested:", nrow(gsea_df), "| Significant:", nrow(sig_df), "\n")

    # ---- Top pathways per direction ----
    top_df <- select_top_gsea(sig_df, g$top_mode)
    if (g$individual_plots && nrow(top_df) > 0) {
      save_gsea_bar_individual(top_df, coll_name, coll, cfg)
    }
    panels[[coll_name]] <- make_gsea_bar(top_df, coll_name, coll$wrap, BAR_THEMES[[g$bar_theme]])
  }

  # ---- Combined Hallmark | KEGG | GO figure ----
  panels <- Filter(Negate(is.null), panels)
  if (length(panels) > 0) {
    combined <- assemble_panels(panels, g$combined_title, g$combined_subtitle,
                                g$annot_sizes[["title"]], g$annot_sizes[["subtitle"]])
    ggsave(paste0(cfg$output_dir, cfg$prefix, "_Combined_GSEA_barplots_horizontal.png"),
           combined, width = g$width, height = g$height, dpi = 400)
    cat("Saved combined", cfg$short_label, "GSEA barplot.\n")
  } else {
    cat("No significant pathways crossed padj < 0.05 across Hallmark, KEGG, or GO.\n")
  }
}


# -----------------------------------------------------------------------------
# 4.5 ssGSEA (single-sample pathway scores) + Wilcoxon PR vs PD
# -----------------------------------------------------------------------------

GO_TERM_REGEX <- "^GOBP_|^GOCC_|^GOMF_"   # true GO terms (drops HPO sets in C5)

# Score each sample for each pathway, save the score matrix, then test PR vs PD
# per pathway (two-sided Wilcoxon rank-sum). Returns unadjusted statistics.
# wilcox_exact = NULL is the wilcox.test default (24h script).
score_and_test_ssgsea <- function(expr_mat, t2g, collection_label, cfg, coldata,
                                  wilcox_exact) {
  gene_sets <- split(t2g$gene_symbol, t2g$gs_name)

  set.seed(SEED)
  ssgsea_param <- ssgseaParam(exprData = expr_mat, geneSets = gene_sets,
                              normalize = TRUE, minSize = 10, maxSize = 10000)
  ssgsea_df <- as.data.frame(gsva(ssgsea_param, verbose = TRUE))
  write.csv(ssgsea_df, paste0(cfg$output_dir, cfg$prefix, "_ssGSEA_", collection_label, "_scores.csv"))

  # Long format: one row per pathway x sample, with the response label attached
  ssgsea_long <- ssgsea_df %>%
    rownames_to_column("Pathway") %>%
    pivot_longer(-Pathway, names_to = "SampleID", values_to = "ssGSEA_score") %>%
    dplyr::left_join(coldata %>% rownames_to_column("SampleID") %>%
                       dplyr::select(SampleID, all_of(RESPONSE_COL)),
                     by = "SampleID")

  # delta = mean(PR) - mean(PD); positive = pathway higher in responders
  ssgsea_long %>%
    dplyr::group_by(Pathway) %>%
    dplyr::summarise(
      mean_PR = mean(ssGSEA_score[.data[[RESPONSE_COL]] == "PR"], na.rm = TRUE),
      mean_PD = mean(ssGSEA_score[.data[[RESPONSE_COL]] == "PD"], na.rm = TRUE),
      delta   = mean_PR - mean_PD,
      p_value = wilcox.test(ssGSEA_score[.data[[RESPONSE_COL]] == "PR"],
                            ssGSEA_score[.data[[RESPONSE_COL]] == "PD"],
                            exact = wilcox_exact)$p.value,
      .groups = "drop"
    )
}

# Multiple-testing correction and output tables. Variants, as in the originals:
#   padj_go_prefilter  (Placebo): GO-only terms selected, then BH; sorted table
#   padj_go_postfilter (6w)     : BH over all C5 sets, then GO-only kept
#   rawp_individual    (24h)    : BH over all sets, 'sig' column, *_Sig_Sorted.csv
finalise_ssgsea_stats <- function(stats, collection_label, cfg) {
  mode     <- cfg$ssgsea$mode
  out_stem <- paste0(cfg$output_dir, cfg$prefix, "_ssGSEA_", collection_label)
  is_go    <- collection_label == "GO"

  if (mode == "padj_go_prefilter") {
    if (is_go) stats <- stats %>% dplyr::filter(grepl(GO_TERM_REGEX, Pathway))
    stats <- stats %>%
      dplyr::mutate(padj = p.adjust(p_value, method = "BH")) %>%
      dplyr::arrange(padj, dplyr::desc(delta))

  } else if (mode == "padj_go_postfilter") {
    stats <- stats %>% dplyr::mutate(padj = p.adjust(p_value, method = "BH"))
    if (is_go) stats <- stats %>% dplyr::filter(grepl(GO_TERM_REGEX, Pathway))

  } else if (mode == "rawp_individual") {
    stats <- stats %>%
      dplyr::mutate(padj = p.adjust(p_value, method = "BH"), sig = padj < PADJ_CUTOFF) %>%
      dplyr::arrange(padj)
  }

  # ---- Full table (every pathway tested) ----
  write.csv(stats, paste0(out_stem, "_pathway_stats.csv"), row.names = FALSE)

  # ---- Significant-only table ----
  if (mode == "rawp_individual") {
    sig <- stats %>% dplyr::filter(sig) %>% dplyr::arrange(padj)
    write.csv(sig, paste0(out_stem, "_Sig_Sorted.csv"), row.names = FALSE)
    cat(collection_label, "— pathways with raw p<0.05:", sum(stats$p_value < 0.05), "\n")
  } else {
    sig <- stats %>%
      dplyr::filter(!is.na(padj), padj < PADJ_CUTOFF) %>%
      dplyr::arrange(padj, dplyr::desc(delta))
    write.csv(sig, paste0(out_stem, "_Sorted.csv"), row.names = FALSE)
  }

  cat(collection_label, "— tested:", nrow(stats),
      "| padj <", PADJ_CUTOFF, ":", nrow(sig), "\n")
  list(stats = stats, sig = sig)
}

# Top n pathways per direction by delta (Placebo / 6w).
select_top_delta <- function(sig_df, n = N_TOP_PATHWAYS) {
  dplyr::bind_rows(
    sig_df %>% dplyr::filter(delta > 0) %>% dplyr::slice_max(delta, n = n, with_ties = FALSE),
    sig_df %>% dplyr::filter(delta < 0) %>% dplyr::slice_min(delta, n = n, with_ties = FALSE)
  ) %>%
    dplyr::mutate(Direction   = ifelse(delta > 0, "Higher in PR", "Higher in PD"),
                  Description = Pathway)   # clean_labels() works on 'Description'
}

# Delta bar plot in the same style as the GSEA panels (Placebo / 6w).
make_ssgsea_bar <- function(plot_df, title, wrap, bar_theme) {
  if (nrow(plot_df) == 0) return(NULL)
  plot_df <- clean_labels(plot_df, width = wrap)
  ggplot(plot_df, aes(reorder(Description, delta), delta, fill = Direction)) +
    geom_col(width = 0.7) + coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    scale_fill_manual(values = c("Higher in PR" = "#2E7D32", "Higher in PD" = "#C62828")) +
    labs(x = NULL, y = "Delta ssGSEA Score (PR - PD)", title = title) +
    theme_classic(base_size = 14) + bar_theme
}

# 24h-style panel: raw Wilcoxon p < 0.05, first/last 5 by delta, saved
# individually and returned for the combined figure.
make_ssgsea_bar_24h <- function(stats, collection_label, coll, cfg) {
  plot_data <- stats %>% dplyr::filter(p_value < 0.05)
  if (collection_label == "GO") {
    plot_data <- plot_data %>% dplyr::filter(grepl(GO_TERM_REGEX, Pathway))
  }

  plot_data <- plot_data %>% dplyr::arrange(dplyr::desc(delta))
  plot_data <- unique(rbind(head(plot_data, N_TOP_PATHWAYS), tail(plot_data, N_TOP_PATHWAYS)))

  if (nrow(plot_data) == 0) {
    cat("No pathways passed p<0.05 for", collection_label, "barplot.\n")
    return(NULL)
  }

  # Readable labels, ordered by delta
  plot_data <- plot_data %>%
    dplyr::mutate(Pathway_Clean = if (!is.null(coll$strip)) gsub(coll$strip, "", Pathway) else Pathway,
                  Pathway_Clean = gsub(GO_TERM_REGEX, "", Pathway_Clean),
                  Pathway_Clean = gsub("_", " ", Pathway_Clean),
                  Pathway_Clean = tools::toTitleCase(tolower(Pathway_Clean))) %>%
    dplyr::arrange(delta) %>%
    dplyr::mutate(Pathway_Clean = factor(Pathway_Clean, levels = Pathway_Clean))

  # Size the figure and labels to the number of pathways
  n_paths     <- nrow(plot_data)
  plot_height <- max(5, n_paths * 0.35)
  label_size  <- if (n_paths > 20) 7.5 else if (n_paths > 10) 9 else 10

  p <- ggplot(plot_data, aes(x = delta, y = Pathway_Clean, fill = delta > 0)) +
    geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
    scale_fill_manual(values = c("TRUE" = "#2E7D32", "FALSE" = "#C62828"),
                      labels = c("TRUE" = "Enriched in PR", "FALSE" = "Enriched in PD")) +
    theme_bw(base_size = 12) +
    labs(title    = paste0(collection_label, " Pathway: ", cfg$short_label, " "),
         subtitle = paste("Top", n_paths, "most significant pathways | ssGSEA Wilcoxon p < 0.05"),
         x = "Delta Enrichment Score (Mean PR - Mean PD)", y = NULL, fill = "Enrichment") +
    theme(plot.title       = element_text(face = "bold", hjust = 0.5, size = 13),
          plot.subtitle    = element_text(hjust = 0.5, size = 9, colour = "black"),
          axis.text.y      = element_text(size = label_size),
          axis.text.x      = element_text(size = 9),
          legend.position  = "bottom",
          plot.margin      = margin(t = 10, r = 20, b = 10, l = 10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  ggsave(paste0(cfg$output_dir, cfg$prefix, "_ssGSEA_", collection_label, "_Barplot.png"),
         plot = p, width = 10, height = plot_height, dpi = 400, limitsize = FALSE)
  cat("Saved", collection_label, "ssGSEA barplot (", n_paths, "pathways ).\n")

  # Restyle for the combined figure: bold theme + wrapped y labels
  p + BAR_THEMES$ss24 +
    scale_y_discrete(labels = function(x) str_wrap(x, width = coll$wrap))
}

# Full ssGSEA block: score -> test -> correct -> save -> plot.
run_ssgsea_block <- function(dds, coldata, cfg) {
  s       <- cfg$ssgsea
  is_24h  <- s$mode == "rawp_individual"
  w_exact <- if (is_24h) NULL else FALSE   # 24h used the wilcox.test default

  # log2(normalised counts + 1) as ssGSEA input
  expr_mat <- log2(counts(dds, normalized = TRUE) + 1)
  panels   <- list()

  for (coll_name in names(GENESET_COLLECTIONS)) {
    coll  <- GENESET_COLLECTIONS[[coll_name]]
    stats <- score_and_test_ssgsea(expr_mat, coll$t2g, coll_name, cfg, coldata, w_exact)
    fin   <- finalise_ssgsea_stats(stats, coll_name, cfg)

    panels[[coll_name]] <- if (is_24h) {
      make_ssgsea_bar_24h(fin$stats, coll_name, coll, cfg)
    } else {
      make_ssgsea_bar(select_top_delta(fin$sig), coll_name, coll$wrap, BAR_THEMES$y12)
    }
  }

  # ---- Combined Hallmark | KEGG | GO figure ----
  panels <- Filter(Negate(is.null), panels)
  if (length(panels) == 0) {
    cat("No pathways crossed the plotting threshold across Hallmark, KEGG, or GO.\n")
    return(invisible(NULL))
  }

  if (is_24h) {
    combined <- assemble_panels(panels, s$title, "Wilcoxon p < 0.05, mean PR vs mean PD", 20, 14)
    ggsave(paste0(cfg$output_dir, cfg$prefix, "_Combined_ssGSEA_barplots.png"),
           plot = combined, width = 22, height = 13, dpi = 400, limitsize = FALSE)
  } else {
    combined <- assemble_panels(panels, s$title,
                                paste0("Top 5 higher in PR & top 5 higher in PD by delta | Wilcoxon p.adj < ",
                                       PADJ_CUTOFF),
                                18, 13)
    ggsave(paste0(cfg$output_dir, cfg$prefix, "_Combined_ssGSEA_barplots_horizontal.png"),
           combined, width = 18, height = 11, dpi = 400)
  }
  cat("Saved combined", cfg$short_label, "ssGSEA barplot.\n")
}


# -----------------------------------------------------------------------------
# 4.6 Transcription factor activity (DoRothEA + decoupleR ULM)
# -----------------------------------------------------------------------------

# ULM fits, for each TF, a linear model of the gene-level Wald statistic on
# the TF's signed regulon. Positive score = TF more active in PR.
run_tf_block <- function(res, cfg) {
  mode <- cfg$tf$mode

  # ---- Gene-level statistic as a 1-column matrix (genes x contrast) ----
  stat_vec <- res$stat
  names(stat_vec) <- rownames(res)
  stat_vec <- stat_vec[!is.na(stat_vec)]
  mat <- matrix(stat_vec, ncol = 1, dimnames = list(names(stat_vec), cfg$prefix))
  cat("Genes in matrix:", nrow(mat), "\n")

  # ---- Infer TF activities ----
  tf_res <- decoupleR::run_ulm(mat = mat, net = DOROTHEA_NET,
                               .source = "source", .target = "target",
                               .mor = "mor", minsize = 5)

  # ---- Significance (BH unless mode = raw_p, as in the 6w script) ----
  if (mode == "raw_p") {
    n_sig <- sum(tf_res$p_value < 0.05)
    cat("TFs tested:", nrow(tf_res), "| Significant (p<0.05):", n_sig, "\n")
  } else {
    tf_res <- tf_res %>% dplyr::mutate(padj = p.adjust(p_value, method = "BH"))
    n_sig  <- sum(tf_res$padj < PADJ_CUTOFF)
    cat("TFs tested:", nrow(tf_res), "| Significant (padj < 0.05):", n_sig, "\n")
  }
  write.csv(tf_res, paste0(cfg$output_dir, cfg$prefix, "_TF_activity_ULM.csv"), row.names = FALSE)

  # ---- TFs to plot ----
  tf_sig <- if (mode == "raw_p") {
    tf_res %>% dplyr::filter(p_value < 0.05) %>% dplyr::arrange(dplyr::desc(score))
  } else {
    tf_res %>% dplyr::filter(padj < PADJ_CUTOFF) %>% dplyr::arrange(dplyr::desc(score))
  }

  if (mode == "padj_with_fallback" && nrow(tf_sig) == 0) {
    warning("No TFs reached padj < 0.05! Falling back to raw p_value < 0.05 for plotting.")
    tf_sig <- tf_res %>% dplyr::filter(p_value < 0.05) %>% dplyr::arrange(dplyr::desc(score))
  }

  if (mode == "padj_strict" && nrow(tf_sig) == 0) {
    cat("No significant TFs found at padj < 0.05. Barplot was not generated.\n")
    return(invisible(tf_res))
  }

  # Top n most positive and n most negative scores
  n_show  <- min(N_TOP_TFS, nrow(tf_sig))
  tf_plot <- dplyr::bind_rows(head(tf_sig, n_show), tail(tf_sig, n_show)) %>%
    dplyr::distinct(source, .keep_all = TRUE) %>%
    dplyr::mutate(direction = ifelse(score > 0, "Active in PR", "Active in PD"))

  # ---- Plot (24h used larger fonts and wrapped labels) ----
  if (mode == "padj_strict") {
    tf_plot <- tf_plot %>%
      dplyr::mutate(label = str_replace_all(source, "^(TF_|HALLMARK_|KEGG_|GO_)", ""),
                    label = str_replace_all(label, "_", " "),
                    label = str_wrap(label, width = 15),
                    label = factor(label, levels = label[order(score)]))
    title    <- paste("TF activity (ULM) —", cfg$short_label)
    tf_theme <- theme_classic(base_size = 14) +
      theme(plot.title      = element_text(face = "bold", size = 18, hjust = 0.5),
            plot.subtitle   = element_text(size = 13, hjust = 0.5),
            axis.title      = element_text(face = "bold", size = 15),
            axis.text.x     = element_text(size = 14),
            axis.text.y     = element_text(size = 14, lineheight = 0.85),
            legend.text     = element_text(size = 13),
            legend.position = "top")
  } else {
    tf_plot  <- tf_plot %>% dplyr::mutate(label = factor(source, levels = source[order(score)]))
    title    <- paste("TF activity (ULM) —", cfg$prefix)
    tf_theme <- theme_classic(base_size = 11) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5),
            legend.position = "top")
  }

  p <- ggplot(tf_plot, aes(x = label, y = score, fill = direction)) +
    geom_col(width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c("Active in PR" = "#4DBBD5", "Active in PD" = "#F8766D")) +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey30") +
    labs(title = title, subtitle = cfg$tf$subtitle,
         x = "Transcription factor", y = "Activity score (ULM)", fill = NULL) +
    tf_theme

  ggsave(paste0(cfg$output_dir, cfg$prefix, "_TF_activity_barplot.png"),
         p, width = 10, height = 8, dpi = 400)
  cat("Saved TF activity barplot.\n")
  invisible(tf_res)
}


#===============================================================================
# 5. MAIN DRIVER — one full PR vs PD analysis for a single timepoint
#===============================================================================

run_response_comparison <- function(cfg, GEX, TRT) {

  cat("\n", strrep("=", 78), "\n", "  ", cfg$title_label, "\n", strrep("=", 78), "\n", sep = "")
  if (!dir.exists(cfg$output_dir)) dir.create(cfg$output_dir, recursive = TRUE)

  # ---- Step 1. Sample selection ----
  coldata <- subset_pr_pd(TRT, cfg$treatment)
  counts  <- GEX[, rownames(coldata)]          # align counts to selected samples

  # ---- Step 2. Differential expression ----
  dds <- fit_deseq2(counts, coldata, cfg)
  res <- results(dds, contrast = c(RESPONSE_COL, "PR", "PD"))   # log2(PR / PD)
  summary(res)
  write.csv(as.data.frame(res), paste0(PIPELINES_DIR, cfg$de_csv))
  cat("Saved", cfg$short_label, "PR vs PD DESeq2 results.\n")

  # ---- Step 3. Sample-level QC ----
  report_cooks(dds, coldata)
  plot_pca(dds, cfg)

  # ---- Step 4. Volcano plot ----
  make_volcano(res, prefix = cfg$prefix, title_label = cfg$title_label,
               output_dir = cfg$output_dir)

  # ---- Step 5. Pre-ranked GSEA ----
  run_gsea_block(res, cfg)

  # ---- Step 6. ssGSEA ----
  run_ssgsea_block(dds, coldata, cfg)

  # ---- Step 7. TF activity ----
  run_tf_block(res, cfg)

  invisible(list(dds = dds, res = res))
}


#===============================================================================
# 6. RUN ALL COMPARISONS
#    To run only one, e.g.:  run_response_comparison(COMPARISONS$h24, GEX, TRT)
#===============================================================================

results_all <- lapply(COMPARISONS, run_response_comparison, GEX = GEX, TRT = TRT)


#===============================================================================
# 7. REPRODUCIBILITY
#===============================================================================

sessionInfo()
