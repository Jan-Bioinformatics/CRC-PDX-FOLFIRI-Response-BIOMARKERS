################################################################################
#
#  Script      : Long_Interaction_Placebo_vs_6wks.R
#  Project     : Transcriptomic profiling of FOLFIRI / 5-FU response in
#                colorectal cancer (CRC) patient-derived xenograft (PDX) models
#  Analysis    : Treatment x Response interaction — does the transcriptional
#                change from Placebo to 6 weeks of treatment differ between
#                responders (PR) and non-responders (PD)?
#  Author      : Janmita Kaverimane Umesh (MSc Bioinformatics, QUB)
#
#  Workflow
#  ---------------------------------------------------------------------------
#    1. Case-level PR / PD lookup and label-consistency check
#    2. Matched cohorts — PDX models sampled at both timepoints of a pair
#    3. Waterfall plots of tumour volume change for each matched cohort
#    4. DESeq2 interaction model (~ Response + Treatment + Response:Treatment)
#    5. HBA2 sanity check (mouse-blood / haemoglobin contamination signal)
#    6. Volcano + MA plots of the interaction term
#    7. Pre-ranked GSEA on the interaction Wald statistic
#    8. ssGSEA — Kruskal-Wallis across 4 groups + paired interaction test
#    9. TF activity — DoRothEA (A/B/C) regulons + decoupleR ULM
#   10. Dual consensus — GSEA NES vs ssGSEA difference-in-differences (ΔΔ)
#
#  Interpretation of the interaction term (Response_3_6wkPR.Treatment6wks)
#    log2FC_interaction = (6wks - Placebo) in PR  -  (6wks - Placebo) in PD
#    Positive = gene rises more (or falls less) with treatment in responders.
#
#  Inputs  (0_data/PDX_Longitudinal/counts/)
#    - CRC_graft_raw_counts_matrix_12022024.csv  : raw counts, genes x samples
#    - Complete_GraftCRC_colData_19022024.csv    : sample metadata
#
#  Outputs
#    - 2_pipelines/.../LONG_INTERACTION/ : fitted DESeq2 object + full DE table
#    - 3_output/.../LONG_INTERACTION/    : plots and enrichment tables
#
#  Notes
#    - All thresholds and parameters are unchanged from the original script.
#    - Placebo vs 6wks is the best-powered pair and is the only one modelled;
#      Placebo vs 24h and 24h vs 6wks are built for waterfall plots only.
#
################################################################################


#===============================================================================
# 0. SESSION SETUP
#===============================================================================

rm(list = ls(all.names = TRUE))     # start from a clean environment
gc()                                # release memory from previous sessions

# ---- 0.1 One-off package installation (set TRUE on a new machine) ----
INSTALL_PACKAGES <- FALSE

if (INSTALL_PACKAGES) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(c("DESeq2", "clusterProfiler", "org.Hs.eg.db", "enrichplot",
                         "AnnotationDbi", "ComplexHeatmap", "GSVA", "dorothea",
                         "decoupleR", "GO.db"), update = FALSE, ask = FALSE)
  install.packages(c("tidyverse", "ggrepel", "pheatmap", "circlize", "gridExtra",
                     "knitr", "readxl", "msigdbr", "patchwork"))
}

# ---- 0.2 Libraries ----
# Core DE
library(DESeq2)          # differential expression (NB GLM) with interaction terms
# Enrichment
library(clusterProfiler) # pre-ranked GSEA
library(msigdbr)         # MSigDB gene-set collections
library(GSVA)            # single-sample GSEA (ssGSEA)
library(dorothea)        # TF -> target regulons
library(decoupleR)       # TF activity inference (ULM)
# Data handling and plotting
library(tidyverse)       # dplyr, tidyr, tibble, stringr, ggplot2
library(ggrepel)         # non-overlapping text labels
library(patchwork)       # multi-panel figure assembly

# ---- 0.3 Project paths ----
setwd("C:/Users/40496110/Downloads/PDX_PROJECT")

PIPELINES_DIR <- "./2_pipelines/PDX_Longitudinal/LONG_INTERACTION/"
OUTPUT_DIR    <- "./3_output/PDX_Longitudinal/LONG_INTERACTION/"
COUNTS_DIR    <- "./0_data/PDX_Longitudinal/counts/"

COUNTS_FILE   <- "CRC_graft_raw_counts_matrix_12022024.csv"
COLDATA_FILE  <- "Complete_GraftCRC_colData_19022024.csv"

dir.create(PIPELINES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_DIR,    recursive = TRUE, showWarnings = FALSE)


#===============================================================================
# 1. ANALYSIS CONFIGURATION
#===============================================================================

PREFIX         <- "Interaction_Placebo_6wks"   # prefix for most output files

# ---- 1.1 Metadata columns and levels ----
RESPONSE_COL   <- "Response_3_6wk"   # response label measured at 3-6 weeks
CASE_COL       <- "Case_ID"          # PDX model identifier (links matched samples)
TRT_LEVELS     <- c("Placebo", "24h", "6wks")
TRT_BASE       <- "Placebo"          # baseline arm of the modelled pair
TRT_LATE       <- "6wks"             # treated arm of the modelled pair

# ---- 1.2 DESeq2 ----
MIN_TOTAL_CNT  <- 10                 # pre-filter: keep genes with >= 10 reads in total
INTERACTION_TERM <- "Response_3_6wkPR.Treatment6wks"
SEED           <- 123

# ---- 1.3 Thresholds ----
PADJ_CUTOFF    <- 0.05               # BH-adjusted significance threshold (all tests)
LFC_CUTOFF     <- 1                  # |log2FC| > 1 (2-fold) for interaction volcano
RESP_VOL_LINE  <- 35                 # reference line on waterfall plots (% volume change)
N_TOP_PATHWAYS <- 5                  # pathways shown per direction in bar plots
N_TOP_TFS      <- 10                 # TFs shown per direction in TF bar plot
N_LABEL_CONS   <- 10                 # pathways labelled on consensus scatter

# ---- 1.4 ssGSEA options ----
SS_SIG_COL     <- "int_padj"   # which test defines significance in ssGSEA bar plots:
                               #   "int_padj" = paired Wilcoxon on per-model change
                               #   "kw_padj"  = Kruskal-Wallis across the 4 groups
CENTER_SS      <- TRUE         # ΔΔ: centre each sample on its median ssGSEA score
SS_DD_TEST     <- "wilcox"     # ΔΔ test: "wilcox" (as in Methods) or "t"


#===============================================================================
# 2. LOAD AND ALIGN INPUTS
#===============================================================================

# ---- 2.1 Raw counts (genes x samples) and sample metadata ----
GEX <- read.csv(paste0(COUNTS_DIR, COUNTS_FILE),  row.names = 1, check.names = FALSE)
TRT <- read.csv(paste0(COUNTS_DIR, COLDATA_FILE), row.names = 1, check.names = FALSE)

# ---- 2.2 Keep samples present in both tables; order counts to match metadata ----
TRT <- TRT[rownames(TRT) %in% colnames(GEX), ]
GEX <- GEX[, rownames(TRT)]

# ---- 2.3 Factor coding ----
TRT$Treatment <- factor(TRT$Treatment, levels = TRT_LEVELS)   # Placebo = reference
TRT$Case_ID   <- as.factor(TRT$Case_ID)

cat("GEX:", dim(GEX), "| TRT:", dim(TRT), "\n")


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
#   gsea_tag    : name used in GSEA output files
#   wrap        : character width for wrapping pathway labels on plots
#   indiv_width : width (in) of single-collection bar plots
#   go_filter   : keep only GOBP/GOCC/GOMF terms in the consensus step
GENESET_COLLECTIONS <- list(
  Hallmark = list(t2g = h_t2g,    gsea_tag = "Hallmarks", wrap = 22, indiv_width = 11, go_filter = FALSE),
  KEGG     = list(t2g = kegg_t2g, gsea_tag = "KEGG",      wrap = 22, indiv_width = 11, go_filter = FALSE),
  GO       = list(t2g = c5_t2g,   gsea_tag = "C5_GO",     wrap = 16, indiv_width = 12, go_filter = TRUE)
)

# ---- 3.2 DoRothEA regulons (high / medium confidence: A, B, C) ----
data("dorothea_hs", package = "dorothea")

DOROTHEA_NET <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B", "C")) %>%
  dplyr::rename(source = tf)            # 'mor' = mode of regulation (+1 / -1)

cat("DoRothEA network rows (A/B/C):", nrow(DOROTHEA_NET), "\n")


#===============================================================================
# 4. HELPER FUNCTIONS
#===============================================================================

# -----------------------------------------------------------------------------
# 4.1 Plot styling and labels
# -----------------------------------------------------------------------------

# Bold theme for horizontal pathway bar plots (GSEA uses y = 15, ssGSEA y = 14)
bold_bar_theme <- function(y_text) {
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title   = element_text(size = 14, face = "bold"),
    axis.text.x  = element_text(size = 14, face = "bold"),
    axis.text.y  = element_text(size = y_text, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text  = element_text(size = 12, face = "bold")
  )
}
THEME_GSEA   <- bold_bar_theme(y_text = 15)
THEME_SSGSEA <- bold_bar_theme(y_text = 14)

# Tidy MSigDB names for bar plots: drop collection prefix, underscores -> spaces,
# title case, and wrap long names over several lines.
clean_labels <- function(df, width = 22) {
  if (is.null(df)) return(NULL)
  df$Description <- df$Description %>%
    str_remove("^HALLMARK_") %>% str_remove("^KEGG_") %>%
    str_remove("^GOBP_") %>% str_remove("^GOCC_") %>% str_remove("^GOMF_") %>%
    str_replace_all("_", " ") %>% str_to_title() %>% str_wrap(width = width)
  df
}

# Same tidy-up without wrapping (consensus scatter labels)
tidy_names <- function(x) {
  x %>% str_remove("^(HALLMARK|KEGG|GOBP|GOCC|GOMF)_") %>%
    str_replace_all("_", " ") %>% str_to_title()
}

# Placeholder panel when a collection has no significant pathways, so the
# combined figure always shows all three collections.
empty_panel <- function(title, bar_theme) {
  ggplot() + theme_void() + bar_theme + labs(title = title) +
    annotate("text", x = 0, y = 0, label = "No significant pathways\n(p.adj < 0.05)",
             size = 5, fontface = "bold")
}

# Horizontal bar plot of a signed effect (NES or ssGSEA interaction).
# Green = higher in PR interaction, red = lower.
make_direction_bar <- function(df, value_col, title, y_label, bar_theme) {
  if (is.null(df)) return(empty_panel(title, bar_theme))
  ggplot(df, aes(reorder(Description, .data[[value_col]]), .data[[value_col]], fill = Direction)) +
    geom_col(width = 0.7) + coord_flip() +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    scale_fill_manual(values = c("Activated" = "#2E7D32", "Suppressed" = "#C62828"),
                      labels = c("Activated"  = "Higher in PR interaction",
                                 "Suppressed" = "Lower in PR interaction"),
                      drop = FALSE) +
    labs(x = NULL, y = y_label, title = title, fill = "Direction") +
    theme_classic(base_size = 14) + bar_theme
}

# Three panels side by side, shared legend at the bottom, overall title.
assemble_panels <- function(plot_list, title, subtitle) {
  wrap_plots(plot_list, ncol = length(plot_list)) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title    = title,
      subtitle = subtitle,
      theme    = theme(
        plot.title    = element_text(size = 20, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, face = "bold", hjust = 0.5)
      )
    ) &
    theme(legend.position = "bottom")
}


# -----------------------------------------------------------------------------
# 4.2 Cohort construction
# -----------------------------------------------------------------------------

# Matched cohort for one pair of timepoints:
#   - keep PDX models (Case_ID) sampled at BOTH level_A and level_B
#   - one sample per model per timepoint (first if replicated)
#   - attach the model-level PR / PD label (models without one are dropped)
build_matched_pairs_with_response <- function(TRT, GEX, level_A, level_B, case_response) {
  sub_meta <- TRT %>% as.data.frame() %>% rownames_to_column("SampleID") %>%
    dplyr::filter(Treatment %in% c(level_A, level_B)) %>%
    dplyr::select(-all_of(RESPONSE_COL))          # replaced by case-level label below

  ids_with_both <- sub_meta %>% dplyr::group_by(Case_ID) %>%
    dplyr::summarise(has_A = level_A %in% Treatment,
                     has_B = level_B %in% Treatment, .groups = "drop") %>%
    dplyr::filter(has_A & has_B) %>% dplyr::pull(Case_ID)

  sub_meta <- sub_meta %>% dplyr::filter(Case_ID %in% ids_with_both) %>%
    dplyr::group_by(Case_ID, Treatment) %>% dplyr::slice(1) %>% dplyr::ungroup() %>%
    dplyr::inner_join(case_response, by = "Case_ID") %>%
    column_to_rownames("SampleID")

  # Reference levels: first timepoint and PD
  sub_meta$Treatment        <- factor(sub_meta$Treatment, levels = c(level_A, level_B))
  sub_meta[[RESPONSE_COL]]  <- factor(sub_meta[[RESPONSE_COL]], levels = c("PD", "PR"))
  sub_meta$Case_ID          <- factor(sub_meta$Case_ID)

  list(meta = sub_meta, counts = GEX[, rownames(sub_meta)])
}

# Waterfall plot: one bar per PDX model, ordered by % tumour volume change.
plot_response_waterfall <- function(meta_df, comparison_label, output_dir,
                                    resp_threshold = RESP_VOL_LINE) {
  plot_data <- meta_df %>% as.data.frame() %>%
    dplyr::distinct(Case_ID, .keep_all = TRUE) %>%          # one bar per model
    dplyr::arrange(dplyr::desc(percentage_vol_change_3_6wk)) %>%
    dplyr::mutate(Case_ID = factor(Case_ID, levels = Case_ID))
  n_pr <- sum(plot_data[[RESPONSE_COL]] == "PR")
  n_pd <- sum(plot_data[[RESPONSE_COL]] == "PD")

  p <- ggplot(plot_data, aes(x = Case_ID, y = percentage_vol_change_3_6wk,
                             fill = .data[[RESPONSE_COL]])) +
    geom_col(width = 1, color = "white", linewidth = 0.05) +
    geom_hline(yintercept = resp_threshold, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    scale_fill_manual(values = c("PD" = "#F8766D", "PR" = "#4DBBD5")) +
    labs(title    = paste("Waterfall Plot:", comparison_label),
         subtitle = paste0("n = ", nrow(plot_data), " | PR(n=", n_pr, ") / PD(n=", n_pd, ")"),
         x = "PDX models (ordered by response)", y = "% tumour volume change (3-6wk avg)",
         fill = "Response") +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "top",
          plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
          plot.subtitle = element_text(size = 11, hjust = 0.5))

  ggsave(paste0(output_dir, gsub(" ", "_", comparison_label), "_waterfall.png"),
         p, width = 8, height = 5, dpi = 400)
  cat("Saved waterfall for", comparison_label, "\n")
}


# -----------------------------------------------------------------------------
# 4.3 Pre-ranked GSEA
# -----------------------------------------------------------------------------

# pvalueCutoff = 1 returns every tested set; filtering is done explicitly.
run_gsea <- function(gene_list, t2g) {
  GSEA(gene_list, exponent = 1, minGSSize = 10, maxGSSize = 10000,
       pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = t2g,
       verbose = TRUE, seed = TRUE, nPermSimple = 10000)
}

# Save significant pathways (sorted by NES) and return the top 5 positive and
# top 5 negative NES for plotting. NULL if none are significant.
get_gsea_plot_data <- function(df, gsea_tag) {
  sig <- df[df$p.adjust < PADJ_CUTOFF, ]
  sig <- sig[order(-sig$NES), ]
  write.csv(sig, paste0(OUTPUT_DIR, PREFIX, "_", gsea_tag, "_GSEA_Sorted.csv"))
  if (nrow(sig) == 0) return(NULL)

  pos <- head(sig[sig$NES > 0, ], N_TOP_PATHWAYS)                                       # highest NES
  neg <- head(sig[sig$NES < 0, ][order(sig$NES[sig$NES < 0]), ], N_TOP_PATHWAYS)        # lowest NES
  out <- rbind(pos, neg)
  out$Direction <- ifelse(out$NES > 0, "Activated", "Suppressed")
  out
}


# -----------------------------------------------------------------------------
# 4.4 ssGSEA
# -----------------------------------------------------------------------------

# Score every sample for every pathway, then test three ways:
#   (a) Kruskal-Wallis across the four Treatment x Response groups
#   (b) Interaction effect on group means: (6wk PR - Placebo PR) - (6wk PD - Placebo PD)
#   (c) Paired interaction test: per-model change (6wk - Placebo), PR vs PD (Wilcoxon)
run_ssgsea_interaction <- function(expr_mat, t2g, collection_label, meta_df) {
  gene_sets <- split(t2g$gene_symbol, t2g$gs_name)

  set.seed(SEED)
  ssgsea_param <- ssgseaParam(exprData = expr_mat, geneSets = gene_sets,
                              normalize = TRUE, minSize = 10, maxSize = 10000)
  ssgsea_df <- as.data.frame(gsva(ssgsea_param, verbose = TRUE))
  write.csv(ssgsea_df, paste0(OUTPUT_DIR, PREFIX, "_ssGSEA_", collection_label, "_scores.csv"))

  # Long format: one row per pathway x sample, with group labels attached
  ssgsea_long <- ssgsea_df %>% rownames_to_column("Pathway") %>%
    pivot_longer(-Pathway, names_to = "SampleID", values_to = "score") %>%
    dplyr::left_join(meta_df %>% dplyr::select(SampleID, Group, Treatment, Response,
                                               all_of(CASE_COL)),
                     by = "SampleID")

  # ---- (a) Kruskal-Wallis across the four groups ----
  kw <- ssgsea_long %>% dplyr::group_by(Pathway) %>%
    dplyr::summarise(kw_p = kruskal.test(score ~ Group)$p.value, .groups = "drop") %>%
    dplyr::mutate(kw_padj = p.adjust(kw_p, method = "BH"))

  # ---- (b) Interaction effect from group means ----
  means <- ssgsea_long %>% dplyr::group_by(Pathway, Treatment, Response) %>%
    dplyr::summarise(m = mean(score), .groups = "drop") %>%
    dplyr::mutate(key = paste(Treatment, Response, sep = "_")) %>%
    dplyr::select(Pathway, key, m) %>%
    pivot_wider(names_from = key, values_from = m)
  means$delta_PR    <- means[[paste(TRT_LATE, "PR", sep = "_")]] - means[[paste(TRT_BASE, "PR", sep = "_")]]
  means$delta_PD    <- means[[paste(TRT_LATE, "PD", sep = "_")]] - means[[paste(TRT_BASE, "PD", sep = "_")]]
  means$interaction <- means$delta_PR - means$delta_PD

  # ---- (c) Paired interaction test on per-model change ----
  paired <- ssgsea_long %>%
    dplyr::filter(Treatment %in% c(TRT_BASE, TRT_LATE)) %>%
    dplyr::group_by(Pathway, .data[[CASE_COL]], Response, Treatment) %>%
    dplyr::summarise(score = mean(score), .groups = "drop") %>%
    pivot_wider(names_from = Treatment, values_from = score) %>%
    dplyr::filter(!is.na(.data[[TRT_BASE]]), !is.na(.data[[TRT_LATE]])) %>%
    dplyr::mutate(change = .data[[TRT_LATE]] - .data[[TRT_BASE]]) %>%
    dplyr::group_by(Pathway) %>%
    dplyr::summarise(int_p = tryCatch(wilcox.test(change ~ Response)$p.value,
                                      error = function(e) NA_real_),
                     n_PR = sum(Response == "PR"), n_PD = sum(Response == "PD"),
                     .groups = "drop") %>%
    dplyr::mutate(int_padj = p.adjust(int_p, method = "BH"))

  # ---- Combine and save ----
  stats <- kw %>%
    dplyr::left_join(means %>% dplyr::select(Pathway, delta_PR, delta_PD, interaction), by = "Pathway") %>%
    dplyr::left_join(paired, by = "Pathway") %>%
    dplyr::arrange(kw_p)

  write.csv(stats, paste0(OUTPUT_DIR, PREFIX, "_ssGSEA_", collection_label, "_pathway_stats.csv"),
            row.names = FALSE)
  cat(collection_label, "| Kruskal padj<0.05:", sum(stats$kw_padj < PADJ_CUTOFF, na.rm = TRUE),
      "| Paired interaction padj<0.05:", sum(stats$int_padj < PADJ_CUTOFF, na.rm = TRUE),
      "| pairs PR/PD:", stats$n_PR[1], "/", stats$n_PD[1], "\n")

  list(stats = stats, scores = ssgsea_df)
}

# Top 5 positive / top 5 negative interaction effects among significant pathways.
get_ss_plot_data <- function(stats, sig_col = SS_SIG_COL) {
  sig <- stats[!is.na(stats[[sig_col]]) & stats[[sig_col]] < PADJ_CUTOFF, ]
  if (nrow(sig) == 0) return(NULL)

  pos <- head(sig[sig$interaction > 0, ][order(-sig$interaction[sig$interaction > 0]), ], N_TOP_PATHWAYS)
  neg <- head(sig[sig$interaction < 0, ][order( sig$interaction[sig$interaction < 0]), ], N_TOP_PATHWAYS)
  out <- rbind(pos, neg)
  out$Direction   <- ifelse(out$interaction > 0, "Activated", "Suppressed")
  out$Description <- out$Pathway                 # clean_labels() works on 'Description'
  out
}


# -----------------------------------------------------------------------------
# 4.5 Difference-in-differences (ΔΔ) and dual consensus
# -----------------------------------------------------------------------------

# ssGSEA ΔΔ = (late - base)_PR - (late - base)_PD from group means, plus a
# per-model test of PR vs PD change. Optional per-sample centring removes a
# global sample-level shift so it is not mistaken for pathway-specific change.
compute_ssgsea_dd <- function(ssgsea_scores_df, meta_df, center = CENTER_SS, test = SS_DD_TEST) {

  scores <- as.matrix(ssgsea_scores_df)
  if (center) scores <- sweep(scores, 2, apply(scores, 2, median, na.rm = TRUE), "-")

  scores_long <- as.data.frame(scores) %>% rownames_to_column("Pathway") %>%
    pivot_longer(-Pathway, names_to = "SampleID", values_to = "Score") %>%
    dplyr::left_join(meta_df %>% dplyr::select(SampleID, Treatment,
                                               Response = all_of(RESPONSE_COL),
                                               Case     = all_of(CASE_COL)),
                     by = "SampleID") %>%
    dplyr::mutate(Treatment = as.character(Treatment), Response = as.character(Response))

  # ---- Group means and ΔΔ ----
  group_means <- scores_long %>%
    dplyr::group_by(Pathway, Treatment, Response) %>%
    dplyr::summarise(mean_score = mean(Score, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(Group = paste(Treatment, Response, sep = "_")) %>%
    dplyr::select(Pathway, Group, mean_score) %>%
    pivot_wider(names_from = Group, values_from = mean_score)
  group_means$delta_PR    <- group_means[[paste0(TRT_LATE, "_PR")]] - group_means[[paste0(TRT_BASE, "_PR")]]
  group_means$delta_PD    <- group_means[[paste0(TRT_LATE, "_PD")]] - group_means[[paste0(TRT_BASE, "_PD")]]
  group_means$delta_delta <- group_means$delta_PR - group_means$delta_PD

  # ---- Per-model change (late - base); replicates averaged so each model counts once ----
  case_deltas <- scores_long %>%
    dplyr::filter(Treatment %in% c(TRT_BASE, TRT_LATE)) %>%
    dplyr::group_by(Pathway, Case, Response, Treatment) %>%
    dplyr::summarise(Score = mean(Score, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Treatment, values_from = Score) %>%
    dplyr::filter(!is.na(.data[[TRT_BASE]]), !is.na(.data[[TRT_LATE]])) %>%
    dplyr::mutate(case_delta = .data[[TRT_LATE]] - .data[[TRT_BASE]])

  # ---- PR vs PD test on per-model change ----
  pval_df <- case_deltas %>% dplyr::group_by(Pathway) %>%
    dplyr::summarise(
      p_value = tryCatch({
        x <- case_delta[Response == "PR"]; y <- case_delta[Response == "PD"]
        if (test == "wilcox") wilcox.test(x, y)$p.value else t.test(x, y)$p.value
      }, error = function(e) NA_real_),
      n_PR = sum(Response == "PR"), n_PD = sum(Response == "PD"),
      .groups = "drop") %>%
    dplyr::mutate(padj = p.adjust(p_value, method = "BH"))

  dplyr::full_join(group_means, pval_df, by = "Pathway")
}

# Classify each pathway by agreement between GSEA (NES) and ssGSEA (ΔΔ):
#   Consensus up / down    : significant in both, same sign
#   Conflicting direction  : significant in both, opposite sign
#   Significant in one     : significant in only one method
#   Not significant        : neither
plot_interaction_dual_consensus <- function(gsea_df, ssgsea_dd_df, collection_label, outdir,
                                            gsea_padj_cutoff = PADJ_CUTOFF,
                                            ssgsea_padj_cutoff = PADJ_CUTOFF,
                                            n_label = N_LABEL_CONS, go_filter = FALSE) {

  gsea_sub <- gsea_df %>% dplyr::select(ID, Description, NES, gsea_padj = p.adjust)
  if (go_filter) gsea_sub <- gsea_sub %>% dplyr::filter(grepl("^GOBP_|^GOCC_|^GOMF_", ID))

  merged <- dplyr::inner_join(gsea_sub, ssgsea_dd_df, by = c("ID" = "Pathway"))
  if (nrow(merged) == 0) { cat("No shared pathways for", collection_label, "\n"); return(NULL) }

  # ---- Consensus classification ----
  sig_gsea   <- !is.na(merged$gsea_padj) & merged$gsea_padj < gsea_padj_cutoff
  sig_ssgsea <- !is.na(merged$padj)      & merged$padj      < ssgsea_padj_cutoff

  lab_up   <- "Consensus: higher in PR interaction"
  lab_down <- "Consensus: lower in PR interaction"
  merged$consensus <- "Not significant"
  merged$consensus[xor(sig_gsea, sig_ssgsea)] <- "Significant in one method"
  merged$consensus[sig_gsea & sig_ssgsea & sign(merged$NES) != sign(merged$delta_delta)] <- "Conflicting direction"
  merged$consensus[sig_gsea & sig_ssgsea & merged$NES > 0 & merged$delta_delta > 0] <- lab_up
  merged$consensus[sig_gsea & sig_ssgsea & merged$NES < 0 & merged$delta_delta < 0] <- lab_down
  merged$consensus <- factor(merged$consensus,
                             levels = c(lab_up, lab_down, "Significant in one method",
                                        "Conflicting direction", "Not significant"))
  merged$Label <- tidy_names(merged$ID)

  cat(collection_label, ": Up =", sum(merged$consensus == lab_up),
      "| Down =", sum(merged$consensus == lab_down),
      "| One method =", sum(merged$consensus == "Significant in one method"),
      "| Conflicting =", sum(merged$consensus == "Conflicting direction"),
      "| NS =", sum(merged$consensus == "Not significant"), "\n")

  write.csv(merged, paste0(outdir, "InteractionDualConsensus_", collection_label, "_table.csv"),
            row.names = FALSE)

  # ---- Plot: grey points first so coloured points sit on top ----
  merged     <- merged[order(merged$consensus != "Not significant"), ]
  label_pool <- merged[merged$consensus %in% c(lab_up, lab_down), ]
  label_pool <- label_pool[order(-abs(label_pool$delta_delta)), ]   # label largest |ΔΔ|
  top_labels <- head(label_pool, n_label)

  p <- ggplot(merged, aes(NES, delta_delta, colour = consensus)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_point(alpha = 0.8, size = 2.6) +
    geom_text_repel(data = top_labels, aes(label = Label), size = 4.2, colour = "black",
                    fontface = "bold", max.overlaps = Inf, show.legend = FALSE,
                    box.padding = 0.5, force = 2, min.segment.length = 0,
                    segment.colour = "grey50", bg.color = "white", bg.r = 0.15) +
    scale_colour_manual(values = c(setNames("#2E7D32", lab_up),
                                   setNames("#C62828", lab_down),
                                   "Significant in one method" = "#F4B942",
                                   "Conflicting direction"     = "#9B59B6",
                                   "Not significant"           = "grey80"),
                        name = NULL, drop = FALSE) +
    labs(title    = paste0("Dual-GSEA interaction consensus — ", collection_label),
         subtitle = paste0("GSEA interaction NES vs ssGSEA ΔΔ [(", TRT_LATE, " − ", TRT_BASE, ")PR − (",
                           TRT_LATE, " − ", TRT_BASE, ")PD]",
                           if (CENTER_SS) " | sample-centred ssGSEA" else ""),
         x = "GSEA interaction NES",
         y = "ssGSEA ΔΔ (difference-in-differences)") +
    theme_bw(base_size = 15) +
    theme(plot.title      = element_text(face = "bold", size = 18),
          plot.subtitle   = element_text(size = 12, colour = "grey30"),
          axis.title      = element_text(face = "bold"),
          legend.position = "bottom",
          legend.text     = element_text(size = 12, face = "bold")) +
    guides(colour = guide_legend(nrow = 2, override.aes = list(size = 4))) +
    coord_cartesian(clip = "off")

  ggsave(paste0(outdir, "InteractionDualConsensus_", collection_label, ".png"),
         plot = p, width = 10, height = 8.5, dpi = 300)
  invisible(merged)
}


#===============================================================================
# 5. STEP 1 — CASE-LEVEL PR / PD LOOKUP
#    Response is a property of the PDX model, so every sample from one
#    Case_ID should carry the same label. Check, then build a lookup table.
#===============================================================================

case_response_check <- TRT %>% as.data.frame() %>% rownames_to_column("SampleID") %>%
  dplyr::filter(.data[[RESPONSE_COL]] %in% c("PR", "PD")) %>%
  dplyr::group_by(Case_ID) %>%
  dplyr::summarise(n_distinct_response = dplyr::n_distinct(.data[[RESPONSE_COL]]), .groups = "drop")
cat("Case_IDs with inconsistent PR/PD labels:", sum(case_response_check$n_distinct_response > 1), "\n")

case_response <- TRT %>% as.data.frame() %>% rownames_to_column("SampleID") %>%
  dplyr::filter(.data[[RESPONSE_COL]] %in% c("PR", "PD")) %>%
  dplyr::distinct(Case_ID, Response_3_6wk)
cat("Total cases with known PR/PD:", nrow(case_response), "\n")


#===============================================================================
# 6. STEP 2 — MATCHED COHORTS
#===============================================================================

pair_p24_int <- build_matched_pairs_with_response(TRT, GEX, "Placebo", "24h",  case_response)
pair_p6_int  <- build_matched_pairs_with_response(TRT, GEX, "Placebo", "6wks", case_response)
pair_246_int <- build_matched_pairs_with_response(TRT, GEX, "24h",     "6wks", case_response)

# Cell sizes (Treatment x Response) — determines power of the interaction model
cat("\nPlacebo vs 24h cell sizes:\n");  print(table(pair_p24_int$meta$Treatment, pair_p24_int$meta[[RESPONSE_COL]]))
cat("\nPlacebo vs 6wks cell sizes:\n"); print(table(pair_p6_int$meta$Treatment,  pair_p6_int$meta[[RESPONSE_COL]]))
cat("\n24h vs 6wks cell sizes:\n");     print(table(pair_246_int$meta$Treatment, pair_246_int$meta[[RESPONSE_COL]]))
# NOTE: 24h vs 6wks has only 2 PR cases — underpowered for interaction
#       modelling; treat with caution / consider excluding.


#===============================================================================
# 7. STEP 3 — WATERFALL PLOTS
#===============================================================================

plot_response_waterfall(pair_p24_int$meta, "Placebo vs 24h - matched PR-PD cases",  OUTPUT_DIR)
plot_response_waterfall(pair_p6_int$meta,  "Placebo vs 6wks - matched PR-PD cases", OUTPUT_DIR)
plot_response_waterfall(pair_246_int$meta, "24h vs 6wks - matched PR-PD cases",     OUTPUT_DIR)


#===============================================================================
# 8. STEP 4 — DESeq2 INTERACTION MODEL (Placebo vs 6wks)
#===============================================================================

# ---- 8.1 Build and pre-filter ----
# Main effects + interaction. The interaction coefficient tests whether the
# Placebo -> 6wks change differs between PR and PD.
dds_int_p6 <- DESeqDataSetFromMatrix(
  countData = pair_p6_int$counts,
  colData   = pair_p6_int$meta,
  design    = ~ Response_3_6wk + Treatment + Response_3_6wk:Treatment
)
dds_int_p6 <- dds_int_p6[rowSums(counts(dds_int_p6)) >= MIN_TOTAL_CNT, ]
cat("Genes going into DESeq2:", nrow(dds_int_p6), "| Samples:", ncol(dds_int_p6), "\n")

# ---- 8.2 Fit (default parametric dispersion fit) and save ----
set.seed(SEED)
dds_int_p6 <- DESeq(dds_int_p6)
saveRDS(dds_int_p6, paste0(PIPELINES_DIR, "dds_interaction_Placebo_6wks_fitted.rds"))
print(resultsNames(dds_int_p6))   # confirm the interaction term name before extracting

# ---- 8.3 Extract interaction results ----
res_interaction_p6 <- results(dds_int_p6, name = INTERACTION_TERM)
summary(res_interaction_p6)
write.csv(as.data.frame(res_interaction_p6),
          paste0(PIPELINES_DIR, "Interaction_Placebo_vs_6wks_ResponsePRvsPD.csv"))

# ---- 8.4 Significant genes ----
sig_interaction_p6_05 <- as.data.frame(res_interaction_p6) %>%
  dplyr::filter(!is.na(padj), padj < PADJ_CUTOFF) %>% dplyr::arrange(padj)
cat("Genes significant at padj < 0.05:", nrow(sig_interaction_p6_05), "\n")
write.csv(sig_interaction_p6_05,
          paste0(OUTPUT_DIR, "Interaction_Placebo_6wks_Significant_Genes_padj05.csv"),
          row.names = FALSE)


#===============================================================================
# 9. STEP 5 — HBA2 SANITY CHECK
#    Haemoglobin genes can reflect blood content rather than tumour biology.
#    Check whether HBA2 detection is driven by a few samples in one group.
#===============================================================================

if ("HBA2" %in% rownames(dds_int_p6)) {
  hba2_counts <- counts(dds_int_p6, normalized = TRUE)["HBA2", ]
  hba2_df <- data.frame(SampleID = names(hba2_counts), NormCount = hba2_counts) %>%
    dplyr::left_join(as.data.frame(colData(dds_int_p6)) %>% rownames_to_column("SampleID"),
                     by = "SampleID") %>%
    dplyr::mutate(Group = paste(Treatment, .data[[RESPONSE_COL]], sep = ":"))

  print(hba2_df %>% dplyr::group_by(Group) %>%
          dplyr::summarise(n            = dplyr::n(),
                           n_detected   = sum(NormCount > 0),
                           pct_detected = round(100 * mean(NormCount > 0), 1),
                           median_count = median(NormCount)))
}


#===============================================================================
# 10. STEP 6 — VOLCANO AND MA PLOTS (interaction term)
#===============================================================================

# ---- 10.1 Volcano ----
res_int_df <- as.data.frame(res_interaction_p6) %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::mutate(gene   = rownames(.),
                status = dplyr::case_when(
                  padj < PADJ_CUTOFF & log2FoldChange >  LFC_CUTOFF ~ "UP",
                  padj < PADJ_CUTOFF & log2FoldChange < -LFC_CUTOFF ~ "DOWN",
                  padj < PADJ_CUTOFF                                ~ "Significant (small effect)",
                  TRUE                                              ~ "NS"))

n_up        <- sum(res_int_df$status == "UP")
n_down      <- sum(res_int_df$status == "DOWN")
n_sig_small <- sum(res_int_df$status == "Significant (small effect)")
label_genes <- res_int_df %>% dplyr::filter(padj < PADJ_CUTOFF)    # label all significant genes

p_volcano_int <- ggplot(res_int_df, aes(x = log2FoldChange, y = -log10(padj), color = status)) +
  geom_point(data = dplyr::filter(res_int_df, status == "NS"), alpha = 0.25, size = 0.8) +
  geom_point(data = dplyr::filter(res_int_df, status != "NS"), alpha = 0.85, size = 2) +
  geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_color_manual(
    values = c("NS" = "grey75", "UP" = "#C62828", "DOWN" = "#2E7D32",
               "Significant (small effect)" = "#F9A825"),
    labels = c("NS"   = "Not significant",
               "UP"   = paste0("Up in PR interaction (n=", n_up, ")"),
               "DOWN" = paste0("Down in PR interaction (n=", n_down, ")"),
               "Significant (small effect)" = paste0("Significant, |log2FC|<", LFC_CUTOFF,
                                                     " (n=", n_sig_small, ")")),
    name = NULL) +
  geom_text_repel(data = label_genes, aes(label = gene), size = 3.5, fontface = "italic",
                  color = "black", max.overlaps = 20, box.padding = 0.4, segment.size = 0.3,
                  show.legend = FALSE) +
  labs(title    = "Volcano plot - Treatment x Response interaction",
       subtitle = "Placebo vs 6wks | tests whether the treatment effect differs between PR and PD",
       x = "log2 Fold Change (interaction term)", y = "-log10(adjusted p-value)") +
  theme_classic(base_size = 13) +
  theme(plot.title      = element_text(face = "bold", hjust = 0.5, size = 15),
        plot.subtitle   = element_text(hjust = 0.5, size = 10, color = "grey40"),
        legend.position = "bottom")

ggsave(paste0(OUTPUT_DIR, PREFIX, "_Interaction_VolcanoPlot.png"),
       p_volcano_int, width = 10, height = 8, dpi = 400)

# ---- 10.2 MA plot (significance by padj only, no fold-change cut-off) ----
ma_data <- as.data.frame(res_interaction_p6) %>%
  dplyr::filter(!is.na(padj), !is.na(log2FoldChange), baseMean > 0) %>%
  dplyr::mutate(status = dplyr::case_when(
    padj < PADJ_CUTOFF & log2FoldChange > 0 ~ "UP",
    padj < PADJ_CUTOFF & log2FoldChange < 0 ~ "DOWN",
    TRUE ~ "NS"))
n_up_ma   <- sum(ma_data$status == "UP")
n_down_ma <- sum(ma_data$status == "DOWN")

p_ma_int <- ggplot(ma_data, aes(x = log10(baseMean), y = log2FoldChange, color = status)) +
  geom_point(data = dplyr::filter(ma_data, status == "NS"), alpha = 0.3,  size = 0.7) +
  geom_point(data = dplyr::filter(ma_data, status != "NS"), alpha = 0.85, size = 1.5) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  scale_color_manual(values = c("NS" = "grey70", "UP" = "#D55E00", "DOWN" = "#0072B2"),
                     labels = c("NS"   = "Not significant",
                                "UP"   = paste0("Up (n=", n_up_ma, ")"),
                                "DOWN" = paste0("Down (n=", n_down_ma, ")")),
                     name = NULL) +
  labs(title = "MA plot - Treatment x Response interaction", subtitle = "Placebo vs 6wks | padj < 0.05",
       x = "Mean expression [log10(baseMean)]", y = "log2 Fold Change (interaction term)") +
  theme_classic(base_size = 13) +
  theme(plot.title      = element_text(face = "bold", hjust = 0.5, size = 15),
        plot.subtitle   = element_text(hjust = 0.5, size = 10, color = "grey40"),
        legend.position = "bottom")

ggsave(paste0(OUTPUT_DIR, PREFIX, "_Interaction_MAPlot.png"),
       p_ma_int, width = 9, height = 7, dpi = 400)


#===============================================================================
# 11. STEP 7 — PRE-RANKED GSEA ON THE INTERACTION STATISTIC
#===============================================================================

# ---- 11.1 Ranked gene list (Wald statistic of the interaction term) ----
res_int_ranked  <- as.data.frame(res_interaction_p6)
res_int_ranked  <- res_int_ranked[!is.na(res_int_ranked$stat), ]
geneList_int_p6 <- sort(setNames(res_int_ranked$stat, rownames(res_int_ranked)), decreasing = TRUE)
cat("Genes in ranked list:", length(geneList_int_p6), "\n")

# ---- 11.2 Run, save, plot per collection ----
gsea_results <- list()   # full GSEA tables, reused in the consensus step
gsea_panels  <- list()

for (coll_name in names(GENESET_COLLECTIONS)) {
  coll <- GENESET_COLLECTIONS[[coll_name]]

  # Run and save every tested pathway
  gsea_df <- as.data.frame(run_gsea(geneList_int_p6, coll$t2g))
  write.csv(gsea_df, paste0(OUTPUT_DIR, PREFIX, "_", coll$gsea_tag, "_GSEA_RESULTS.csv"))
  cat(coll_name, "tested:", nrow(gsea_df), "| Significant:", sum(gsea_df$p.adjust < PADJ_CUTOFF), "\n")
  gsea_results[[coll_name]] <- gsea_df

  # Significant table + top pathways -> bar plot
  plot_df <- clean_labels(get_gsea_plot_data(gsea_df, coll$gsea_tag), width = coll$wrap)
  p <- make_direction_bar(plot_df, "NES", coll_name, "Normalized Enrichment Score", THEME_GSEA)

  ggsave(paste0(OUTPUT_DIR, PREFIX, "_", coll$gsea_tag, "_GSEA_barplot.png"),
         p, width = coll$indiv_width, height = 6, dpi = 400)
  gsea_panels[[coll_name]] <- p
}

# ---- 11.3 Combined Hallmark | KEGG | GO figure ----
combined_int_plot <- assemble_panels(
  gsea_panels,
  title    = "GSEA — Treatment × Response interaction (Placebo vs 6 weeks)",
  subtitle = "Top 5 positive and top 5 negative significant pathways, ordered by NES (p.adj < 0.05)"
)
ggsave(paste0(OUTPUT_DIR, PREFIX, "_Combined_GSEA_barplots_horizontal.png"),
       combined_int_plot, width = 18, height = 14, dpi = 400)
cat("Saved combined interaction GSEA barplot.\n")


#===============================================================================
# 12. STEP 8 — ssGSEA (single-sample pathway scores)
#===============================================================================

# ---- 12.1 Input: log2(normalised counts + 1) and sample groups ----
expr_mat_int_p6 <- log2(counts(dds_int_p6, normalized = TRUE) + 1)

meta_int_p6 <- as.data.frame(colData(dds_int_p6)) %>%
  rownames_to_column("SampleID") %>%
  dplyr::mutate(Treatment = as.character(Treatment),
                Response  = as.character(.data[[RESPONSE_COL]]),
                Group     = paste(Treatment, Response, sep = "_"))   # e.g. "6wks_PR"
cat("Groups:", paste(unique(meta_int_p6$Group), collapse = ", "), "\n")

# ---- 12.2 Score, test, save, plot per collection ----
ssgsea_results <- list()   # scores reused in the ΔΔ step
ssgsea_panels  <- list()

for (coll_name in names(GENESET_COLLECTIONS)) {
  coll <- GENESET_COLLECTIONS[[coll_name]]

  ssgsea_results[[coll_name]] <- run_ssgsea_interaction(expr_mat_int_p6, coll$t2g, coll_name, meta_int_p6)

  plot_df <- clean_labels(get_ss_plot_data(ssgsea_results[[coll_name]]$stats), width = coll$wrap)
  p <- make_direction_bar(plot_df, "interaction", coll_name,
                          "Interaction effect\n(ΔPR − ΔPD, ssGSEA score)", THEME_SSGSEA)

  ggsave(paste0(OUTPUT_DIR, PREFIX, "_ssGSEA_", coll_name, "_barplot.png"),
         p, width = coll$indiv_width, height = 6, dpi = 400)
  ssgsea_panels[[coll_name]] <- p
}

# ---- 12.3 Combined figure ----
combined_ss_int <- assemble_panels(
  ssgsea_panels,
  title    = "ssGSEA — Treatment × Response interaction (Placebo vs 6 weeks)",
  subtitle = ifelse(SS_SIG_COL == "int_padj",
                    "Top 5 up and top 5 down; paired Wilcoxon on per-model change, p.adj < 0.05",
                    "Top 5 up and top 5 down; Kruskal-Wallis across groups, p.adj < 0.05")
)
ggsave(paste0(OUTPUT_DIR, PREFIX, "_ssGSEA_Combined_barplots_horizontal.png"),
       combined_ss_int, width = 18, height = 13, dpi = 400)
cat("Saved combined interaction ssGSEA barplot.\n")


#===============================================================================
# 13. STEP 9 — TRANSCRIPTION FACTOR ACTIVITY (DoRothEA + ULM)
#     Positive score = TF targets shift more towards activation in PR.
#===============================================================================

# ---- 13.1 Interaction statistic as a 1-column matrix (genes x contrast) ----
stat_vec <- res_interaction_p6$stat
names(stat_vec) <- rownames(res_interaction_p6)
stat_vec <- stat_vec[!is.na(stat_vec)]
mat_tf <- matrix(stat_vec, ncol = 1, dimnames = list(names(stat_vec), PREFIX))
cat("Genes in matrix:", nrow(mat_tf), "\n")

# ---- 13.2 Infer TF activities + BH correction ----
tf_int_p6 <- decoupleR::run_ulm(mat = mat_tf, net = DOROTHEA_NET,
                                .source = "source", .target = "target",
                                .mor = "mor", minsize = 5) %>%
  dplyr::mutate(padj = p.adjust(p_value, method = "BH"))

cat("TFs tested:", nrow(tf_int_p6),
    "| raw p<0.05:", sum(tf_int_p6$p_value < 0.05),
    "| padj<0.05:",  sum(tf_int_p6$padj < PADJ_CUTOFF), "\n")
write.csv(tf_int_p6, paste0(OUTPUT_DIR, PREFIX, "_TF_activity_ULM.csv"), row.names = FALSE)

# ---- 13.3 Bar plot: top 10 positive + top 10 negative (padj < 0.05) ----
tf_sig <- tf_int_p6 %>% dplyr::filter(padj < PADJ_CUTOFF)

if (nrow(tf_sig) > 0) {
  tf_pos <- tf_sig %>% dplyr::filter(score > 0) %>% dplyr::arrange(dplyr::desc(score)) %>% head(N_TOP_TFS)
  tf_neg <- tf_sig %>% dplyr::filter(score < 0) %>% dplyr::arrange(score)              %>% head(N_TOP_TFS)

  tf_plot <- dplyr::bind_rows(tf_pos, tf_neg) %>%
    dplyr::mutate(
      source_clean = str_replace_all(source, "_", " ") %>% str_wrap(width = 15),
      direction    = ifelse(score > 0, "Higher in PR interaction", "Lower in PR interaction"),
      source_clean = factor(source_clean, levels = source_clean[order(score)])
    )

  p_tf_int <- ggplot(tf_plot, aes(x = source_clean, y = score, fill = direction)) +
    geom_col(width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c("Higher in PR interaction" = "#4DBBD5",
                                 "Lower in PR interaction"  = "#F8766D")) +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey30") +
    labs(title    = "TF activity (ULM) — Treatment × Response interaction",
         subtitle = "Placebo vs 6 weeks | DoRothEA A/B/C | positive = greater increase in PR | padj < 0.05",
         x = "Transcription factor", y = "Interaction activity score (ULM)", fill = NULL) +
    theme_classic(base_size = 14) +
    theme(plot.title      = element_text(face = "bold", size = 18, hjust = 0.5),
          plot.subtitle   = element_text(size = 13, hjust = 0.5),
          axis.title      = element_text(face = "bold", size = 15),
          axis.text.x     = element_text(size = 14),
          axis.text.y     = element_text(size = 14, face = "bold", lineheight = 0.85),
          legend.text     = element_text(size = 13, face = "bold"),
          legend.position = "top")

  ggsave(paste0(OUTPUT_DIR, PREFIX, "_TF_activity_barplot.png"),
         p_tf_int, width = 10, height = 8, dpi = 400)
  cat("Saved interaction TF activity barplot (padj < 0.05).\n")
} else {
  cat("No significant TFs at padj < 0.05. Barplot not generated.\n")
}


#===============================================================================
# 14. STEP 10 — DUAL CONSENSUS (GSEA NES vs ssGSEA ΔΔ)
#     A pathway is a robust interaction hit when both the rank-based GSEA and
#     the sample-level ssGSEA ΔΔ are significant and agree in direction.
#===============================================================================

consensus_results <- list()

for (coll_name in names(GENESET_COLLECTIONS)) {
  coll  <- GENESET_COLLECTIONS[[coll_name]]
  dd_df <- compute_ssgsea_dd(ssgsea_results[[coll_name]]$scores, meta_int_p6)
  consensus_results[[coll_name]] <- plot_interaction_dual_consensus(
    gsea_results[[coll_name]], dd_df, coll_name, OUTPUT_DIR, go_filter = coll$go_filter
  )
}


#===============================================================================
# 15. REPRODUCIBILITY
#===============================================================================

sessionInfo()
