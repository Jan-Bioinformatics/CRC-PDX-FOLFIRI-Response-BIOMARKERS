rm(list = ls(all.names = TRUE))     # Remove all objects
gc()                                # Trigger garbage collection


################################################################################
#                   DESeq2 Pipeline/ Placebo TIME — RESPONSE
#                         PR vs. PD
################################################################################


library(ggrepel)
library(pheatmap)
library(DESeq2)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(AnnotationDbi)
library(circlize)
library(ComplexHeatmap)
library(gridExtra)
library(knitr)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readxl)
library(msigdbr)
library(readxl)


setwd("C:/Users/40496110/Downloads/PDX_PROJECT")
pipelines     <- "./2_pipelines/PDX_Longitudinal/LONG_RESPONSE/"
output        <- "./3_output/PDX_Longitudinal/LONG_RESPONSE/long_placebo_pr_pd/"
clinicals     <- "./0_data/PDX_Longitudinal/counts/"
annotations   <- "./0_data/PDX_Longitudinal/annotations/"



################################################################################
#  BLOCK 1 — PR vs PD
################################################################################

week   <- "Placebo"
prefix <- "Placebo_RES"

# COUNTS DATA
GEX <- read.csv(paste0(clinicals, "CRC_graft_raw_counts_matrix_12022024.csv"), row.names = 1,
                check.names = FALSE)

# TREATMENT DATA
TRT <- read.csv(paste0(clinicals, "Complete_GraftCRC_colData_19022024.csv"), row.names = 1,
                check.names = FALSE)

# ---- 1. Basic checks ----
cat("GEX dimensions:", dim(GEX), "\n")
cat("TRT dimensions:", dim(TRT), "\n")

# ---- 2. Keep only rows with count data ----
TRT <- TRT[rownames(TRT) %in% colnames(GEX), ]
GEX <- GEX[, rownames(TRT)]
cat("TRT after filtering:", nrow(TRT), "\n")

# ---- 3. Subset to Placebo samples only ----
TRT_placebo <- TRT[TRT$Treatment == "Placebo", ]
cat("Placebo samples total:", nrow(TRT_placebo), "\n")

table(TRT_placebo$Response_3_6wk, useNA = "always")

# ---- 4. Keep only PR and PD ----
TRT_placebo_PRPD <- TRT_placebo[TRT_placebo$Response_3_6wk %in% c("PR", "PD"), ]
TRT_placebo_PRPD$Response_3_6wk <- factor(TRT_placebo_PRPD$Response_3_6wk, levels = c("PD", "PR"))

cat("Placebo samples with PR/PD label:", nrow(TRT_placebo_PRPD), "\n")
table(TRT_placebo_PRPD$Response_3_6wk)

# ---- 5. Align counts ----
GEX_placebo_PRPD <- GEX[, rownames(TRT_placebo_PRPD)]

# ---- 6. Build DESeq2 object (unpaired, same as 24h/6wks response comparisons) ----
dds_placebo_res <- DESeqDataSetFromMatrix(
  countData = GEX_placebo_PRPD,
  colData   = TRT_placebo_PRPD,
  design    = ~ Response_3_6wk
)
dds_placebo_res <- dds_placebo_res[rowSums(counts(dds_placebo_res)) >= 10, ]

cat("Genes going into DESeq2:", nrow(dds_placebo_res), "\n")
cat("Samples:", ncol(dds_placebo_res), "\n")

# ---- 7. Run DESeq2 ----
set.seed(123)
dds_placebo_res <- DESeq(dds_placebo_res, fitType = "glmGamPoi")

# Ensure the directory exists before saving
if (!dir.exists(pipelines)) {
  dir.create(pipelines, recursive = TRUE)
}
saveRDS(dds_placebo_res, paste0(pipelines, "dds_placebo_res_fitted.rds"))



saveRDS(dds_placebo_res, paste0(pipelines, "dds_placebo_res_fitted.rds"))
cat("Fitted object saved to disk.\n")

resultsNames(dds_placebo_res)


#==================================================

res_placebo_res <- results(dds_placebo_res, contrast = c("Response_3_6wk", "PR", "PD"))
summary(res_placebo_res)

write.csv(as.data.frame(res_placebo_res), paste0(pipelines, "Placebo_PR_vs_PD_DESeq2_results.csv"))

vsd_placebo_res <- vst(dds_placebo_res, blind = FALSE)
p_pca_placebo_res <- plotPCA(vsd_placebo_res, intgroup = "Response_3_6wk") +
  labs(title = "PCA: Placebo PR vs PD") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
print(p_pca_placebo_res)


#================================================================================
#        Standard GSEA (Hallmark / KEGG / GO) — Placebo, PR vs PD
#================================================================================

library(msigdbr)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)

h_t2g <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, gene_symbol)

kegg_t2g <- msigdbr(species = "Homo sapiens", collection = "C2",
                    subcollection = "CP:KEGG_LEGACY") %>%
  dplyr::select(gs_name, gene_symbol)

c5_t2g <- msigdbr(species = "Homo sapiens", collection = "C5") %>%
  dplyr::select(gs_name, gene_symbol)

cat("Hallmark sets:", length(unique(h_t2g$gs_name)), "\n")
cat("KEGG sets:", length(unique(kegg_t2g$gs_name)), "\n")
cat("GO sets:", length(unique(c5_t2g$gs_name)), "\n")

# ---- 1. Ranked gene list ----
res_placebo_res_df <- as.data.frame(res_placebo_res)
res_placebo_res_df <- res_placebo_res_df[!is.na(res_placebo_res_df$stat), ]

geneList_placebo_res <- res_placebo_res_df$stat
names(geneList_placebo_res) <- rownames(res_placebo_res_df)
geneList_placebo_res <- sort(geneList_placebo_res, decreasing = TRUE)

cat("Genes in ranked list:", length(geneList_placebo_res), "\n")

# ---- 2. Run GSEA (same settings as 6-week) ----
run_gsea <- function(geneList, t2g) {
  as.data.frame(GSEA(geneList, exponent = 1, minGSSize = 10, maxGSSize = 10000,
                     pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = t2g,
                     verbose = TRUE, seed = TRUE, nPermSimple = 100000, eps = 0))
}

h_df_placebo_res    <- run_gsea(geneList_placebo_res, h_t2g)
kegg_df_placebo_res <- run_gsea(geneList_placebo_res, kegg_t2g)
c5_df_placebo_res   <- run_gsea(geneList_placebo_res, c5_t2g)

write.csv(h_df_placebo_res,    paste0(output, prefix, "_Hallmarks_GSEA_RESULTS.csv"))
write.csv(kegg_df_placebo_res, paste0(output, prefix, "_KEGG_GSEA_RESULTS.csv"))
write.csv(c5_df_placebo_res,   paste0(output, prefix, "_C5_GO_GSEA_RESULTS.csv"))

#================================================================================
#        padj, Sorting, Top-5 Selection & Combined Plot
#================================================================================

# ---- 3. Calculate padj (BH) → filter < 0.05 → sort padj ↑, then NES ↓ ----
sort_sig <- function(df, cutoff = 0.05) {
  df %>%
    dplyr::mutate(padj = p.adjust(pvalue, method = "BH")) %>%
    dplyr::filter(!is.na(padj), padj < cutoff) %>%
    dplyr::arrange(padj, dplyr::desc(NES))
}

h_sig_placebo_res    <- sort_sig(h_df_placebo_res)
kegg_sig_placebo_res <- sort_sig(kegg_df_placebo_res)
c5_sig_placebo_res   <- sort_sig(c5_df_placebo_res)

write.csv(h_sig_placebo_res,    paste0(output, prefix, "_Hallmarks_GSEA_Sorted.csv"), row.names = FALSE)
write.csv(kegg_sig_placebo_res, paste0(output, prefix, "_KEGG_GSEA_Sorted.csv"),      row.names = FALSE)
write.csv(c5_sig_placebo_res,   paste0(output, prefix, "_C5_GO_GSEA_Sorted.csv"),     row.names = FALSE)

cat("Hallmark significant:", nrow(h_sig_placebo_res),
    "| KEGG:", nrow(kegg_sig_placebo_res),
    "| GO:", nrow(c5_sig_placebo_res), "\n")

# ---- 4. Top 5 positive & 5 most negative NES ----
top_pos_neg <- function(sig_df, n = 5) {
  pos <- sig_df %>% dplyr::filter(NES > 0) %>% dplyr::slice_max(NES, n = n, with_ties = FALSE)
  neg <- sig_df %>% dplyr::filter(NES < 0) %>% dplyr::slice_min(NES, n = n, with_ties = FALSE)
  dplyr::bind_rows(pos, neg) %>%
    dplyr::mutate(Direction = ifelse(NES > 0, "Activated", "Suppressed"))
}

# ---- 5. Plot styling (same as 6-week) ----
bold_large_theme <- theme(
  plot.title    = element_text(hjust = 0.5, face = "bold", size = 16),
  axis.title    = element_text(size = 13, face = "bold"),
  axis.text.x   = element_text(size = 12, face = "bold"),
  axis.text.y   = element_text(size = 13, face = "bold", lineheight = 0.85),
  legend.title  = element_text(size = 12, face = "bold"),
  legend.text   = element_text(size = 11, face = "bold")
)

clean_labels <- function(df, width = 22) {
  if (nrow(df) == 0) return(df)
  df$Description <- df$Description %>%
    str_remove("^(HALLMARK_|KEGG_|GOBP_|GOCC_|GOMF_)") %>%
    str_replace_all("_", " ") %>%
    str_to_title() %>%
    str_wrap(width = width)
  df
}

make_gsea_bar <- function(plot_df, title, wrap = 22) {
  if (nrow(plot_df) == 0) return(NULL)
  plot_df <- clean_labels(plot_df, width = wrap)
  ggplot(plot_df, aes(reorder(Description, NES), NES, fill = Direction)) +
    geom_col(width = 0.7) + coord_flip() +
    scale_fill_manual(values = c("Activated" = "#2E7D32", "Suppressed" = "#C62828")) +
    labs(x = NULL, y = "Normalized Enrichment Score", title = title) +
    theme_classic(base_size = 14) + bold_large_theme
}

p1 <- make_gsea_bar(top_pos_neg(h_sig_placebo_res),    "Hallmark", wrap = 22)
p2 <- make_gsea_bar(top_pos_neg(kegg_sig_placebo_res), "KEGG",     wrap = 22)
p3 <- make_gsea_bar(top_pos_neg(c5_sig_placebo_res),   "GO",       wrap = 16)

# ---- 6. Combined plot ----
plot_list <- Filter(Negate(is.null), list(p1, p2, p3))

if (length(plot_list) > 0) {
  combined_placebo_plot <- wrap_plots(plot_list, ncol = length(plot_list)) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title    = "GSEA — Significant Pathways: Placebo (PR vs PD)",
      subtitle = "Top 5 activated & top 5 suppressed by NES | p.adj < 0.05",
      theme = theme(
        plot.title    = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 13, face = "bold", hjust = 0.5)
      )
    ) &
    theme(legend.position = "bottom")
  
  ggsave(paste0(output, prefix, "_Combined_GSEA_barplots_horizontal.png"),
         combined_placebo_plot, width = 18, height = 13, dpi = 400)
  cat("Saved unified horizontal Placebo GSEA barplot assembly.\n")
} else {
  cat("No significant pathways crossed padj < 0.05 across Hallmark, KEGG, or GO.\n")
}
#===============================================================================
#                     ssGSEA pathway — Placebo, PR vs PD
#===============================================================================

library(GSVA)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)
library(stringr)
library(patchwork)

# ---- Settings ----
padj_cutoff <- 0.05   # significance threshold (BH-adjusted)
fs          <- 12    # base font size — change only this

normCount_placebo_res <- counts(dds_placebo_res, normalized = TRUE)
expr_mat_placebo_res  <- log2(normCount_placebo_res + 1)

# ---- 1. Run ssGSEA + Wilcoxon, calculate padj, save full & sorted tables ----
run_ssgsea_collection <- function(expr_mat, t2g, collection_label, prefix, output_dir,
                                  TRT_sub, response_col, cutoff = 0.05) {
  
  gene_sets <- split(t2g$gene_symbol, t2g$gs_name)
  
  set.seed(123)
  ssgsea_param <- ssgseaParam(exprData = expr_mat, geneSets = gene_sets,
                              normalize = TRUE, minSize = 10, maxSize = 10000)
  ssgsea_df <- as.data.frame(gsva(ssgsea_param, verbose = TRUE))
  write.csv(ssgsea_df, paste0(output_dir, prefix, "_ssGSEA_", collection_label, "_scores.csv"))
  
  ssgsea_long <- ssgsea_df %>%
    rownames_to_column("Pathway") %>%
    pivot_longer(-Pathway, names_to = "SampleID", values_to = "ssGSEA_score") %>%
    left_join(TRT_sub %>% rownames_to_column("SampleID") %>%
                dplyr::select(SampleID, all_of(response_col)),
              by = "SampleID")
  
  pathway_stats <- ssgsea_long %>%
    group_by(Pathway) %>%
    summarise(
      mean_PR = mean(ssGSEA_score[.data[[response_col]] == "PR"], na.rm = TRUE),
      mean_PD = mean(ssGSEA_score[.data[[response_col]] == "PD"], na.rm = TRUE),
      delta   = mean_PR - mean_PD,
      p_value = wilcox.test(ssGSEA_score[.data[[response_col]] == "PR"],
                            ssGSEA_score[.data[[response_col]] == "PD"],
                            exact = FALSE)$p.value,
      .groups = "drop"
    )
  
  # GO: keep only BP/CC/MF (drops HPO sets in C5) BEFORE padj
  if (collection_label == "GO") {
    pathway_stats <- pathway_stats %>% filter(grepl("^GOBP_|^GOCC_|^GOMF_", Pathway))
  }
  
  pathway_stats <- pathway_stats %>%
    mutate(padj = p.adjust(p_value, method = "BH")) %>%
    arrange(padj, desc(delta))
  
  write.csv(pathway_stats,
            paste0(output_dir, prefix, "_ssGSEA_", collection_label, "_pathway_stats.csv"),
            row.names = FALSE)
  
  # Filter padj < cutoff → sort padj ↑, then delta ↓
  sig_sorted <- pathway_stats %>%
    dplyr::filter(!is.na(padj), padj < cutoff) %>%
    dplyr::arrange(padj, dplyr::desc(delta))
  
  write.csv(sig_sorted,
            paste0(output_dir, prefix, "_ssGSEA_", collection_label, "_Sorted.csv"),
            row.names = FALSE)
  
  cat(collection_label, "— tested:", nrow(pathway_stats),
      "| padj <", cutoff, ":", nrow(sig_sorted),
      "| min padj:", signif(min(pathway_stats$padj, na.rm = TRUE), 3), "\n")
  
  sig_sorted
}

# ---- 2. Top 5 positive & 5 most negative delta ----
top_pos_neg_delta <- function(sig_df, n = 5) {
  pos <- sig_df %>% dplyr::filter(delta > 0) %>% dplyr::slice_max(delta, n = n, with_ties = FALSE)
  neg <- sig_df %>% dplyr::filter(delta < 0) %>% dplyr::slice_min(delta, n = n, with_ties = FALSE)
  dplyr::bind_rows(pos, neg) %>%
    dplyr::mutate(Direction   = ifelse(delta > 0, "Higher in PR", "Higher in PD"),
                  Description = Pathway)
}

# ---- 3. Styling ----
bold_large_theme <- theme(
  plot.title    = element_text(hjust = 0.5, face = "bold", size = fs + 4),
  axis.title    = element_text(size = fs + 1, face = "bold"),
  axis.text.x   = element_text(size = fs, face = "bold"),
  axis.text.y   = element_text(size = fs, face = "bold", lineheight = 0.85),
  legend.title  = element_text(size = fs, face = "bold"),
  legend.text   = element_text(size = fs - 1, face = "bold")
)

clean_labels <- function(df, width = 22) {
  if (nrow(df) == 0) return(df)
  df$Description <- df$Description %>%
    str_remove("^(HALLMARK_|KEGG_|GOBP_|GOCC_|GOMF_)") %>%
    str_replace_all("_", " ") %>%
    str_to_title() %>%
    str_wrap(width = width)
  df
}

make_ssgsea_bar <- function(plot_df, title, wrap = 22) {
  if (nrow(plot_df) == 0) return(NULL)
  plot_df <- clean_labels(plot_df, width = wrap)
  ggplot(plot_df, aes(reorder(Description, delta), delta, fill = Direction)) +
    geom_col(width = 0.7) + coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    scale_fill_manual(values = c("Higher in PR" = "#2E7D32", "Higher in PD" = "#C62828")) +
    labs(x = NULL, y = "Delta ssGSEA Score (PR - PD)", title = title) +
    theme_classic(base_size = 14) + bold_large_theme
}

# ---- 4. Run ----
ss_hallmark_placebo <- run_ssgsea_collection(expr_mat_placebo_res, h_t2g,    "Hallmark", prefix, output,
                                             TRT_placebo_PRPD, "Response_3_6wk", cutoff = padj_cutoff)
ss_kegg_placebo     <- run_ssgsea_collection(expr_mat_placebo_res, kegg_t2g, "KEGG",     prefix, output,
                                             TRT_placebo_PRPD, "Response_3_6wk", cutoff = padj_cutoff)
ss_go_placebo       <- run_ssgsea_collection(expr_mat_placebo_res, c5_t2g,   "GO",       prefix, output,
                                             TRT_placebo_PRPD, "Response_3_6wk", cutoff = padj_cutoff)

s1 <- make_ssgsea_bar(top_pos_neg_delta(ss_hallmark_placebo), "Hallmark", wrap = 22)
s2 <- make_ssgsea_bar(top_pos_neg_delta(ss_kegg_placebo),     "KEGG",     wrap = 22)
s3 <- make_ssgsea_bar(top_pos_neg_delta(ss_go_placebo),       "GO",       wrap = 16)

# ---- 5. Combined plot ----
ss_plot_list <- Filter(Negate(is.null), list(s1, s2, s3))

if (length(ss_plot_list) > 0) {
  combined_ss_placebo_plot <- wrap_plots(ss_plot_list, ncol = length(ss_plot_list)) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title    = "ssGSEA — Significant Pathways: Placebo (PR vs PD)",
      subtitle = paste0("Top 5 higher in PR & top 5 higher in PD by delta | Wilcoxon p.adj < ", padj_cutoff),
      theme = theme(
        plot.title    = element_text(size = fs + 6, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = fs + 1, face = "bold", hjust = 0.5)
      )
    ) &
    theme(legend.position = "bottom")
  
  ggsave(paste0(output, prefix, "_Combined_ssGSEA_barplots_horizontal.png"),
         combined_ss_placebo_plot, width = 18, height = 11, dpi = 400)
  cat("Saved unified horizontal Placebo ssGSEA barplot assembly.\n")
} else {
  cat("No pathways crossed padj <", padj_cutoff, "across Hallmark, KEGG, or GO.\n")
}



#-----------------------------------------------------------------------------------
#                            TF Activity (DoRothEA)
#-----------------------------------------------------------------------------------

library(dorothea)
library(dplyr)
library(ggplot2)

data("dorothea_hs", package = "dorothea")

net <- dorothea_hs %>%
  filter(confidence %in% c("A", "B", "C")) %>%
  rename(source = tf) %>%
  mutate(weight = as.numeric(mor))

cat("DoRothEA network rows (A/B/C confidence):", nrow(net), "\n")

stat_vec <- res_placebo_res$stat
names(stat_vec) <- rownames(res_placebo_res)
stat_vec <- stat_vec[!is.na(stat_vec)]

mat <- matrix(stat_vec, ncol = 1, dimnames = list(names(stat_vec), prefix))
cat("Genes in matrix:", nrow(mat), "\n")

tf_placebo_res <- decoupleR::run_ulm(
  mat     = mat,
  net     = net,
  .source = "source",
  .target = "target",
  .mor    = "mor",
  minsize = 5
)

# --- CALCULATE PADJ VALUE HERE ---
tf_placebo_res <- tf_placebo_res %>%
  mutate(padj = p.adjust(p_value, method = "BH"))

# Count significant hits based on padj
n_sig <- tf_placebo_res %>% filter(padj < 0.05) %>% nrow()
cat("TFs tested:", nrow(tf_placebo_res), "| Significant (padj < 0.05):", n_sig, "\n")

write.csv(tf_placebo_res, paste0(output, prefix, "_TF_activity_ULM.csv"), row.names = FALSE)

# --- FILTER BY PADJ INSTEAD OF P VALUE ---
tf_sig <- tf_placebo_res %>% filter(padj < 0.05) %>% arrange(desc(score))

# Fallback: If no TFs pass padj < 0.05, take the top 20 nominal p-value TFs to avoid an empty plot
if (nrow(tf_sig) == 0) {
  warning("No TFs reached padj < 0.05! Falling back to raw p_value < 0.05 for plotting.")
  tf_sig <- tf_placebo_res %>% filter(p_value < 0.05) %>% arrange(desc(score))
}

n_show <- min(10, nrow(tf_sig))

tf_plot <- bind_rows(tf_sig %>% head(n_show), tf_sig %>% tail(n_show)) %>%
  distinct(source, .keep_all = TRUE) %>%
  mutate(
    direction = ifelse(score > 0, "Active in PR", "Active in PD"),
    source    = factor(source, levels = source[order(score)])
  )

p_tf_placebo_res <- ggplot(tf_plot, aes(x = source, y = score, fill = direction)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Active in PR" = "#4DBBD5", "Active in PD" = "#F8766D")) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey30") +
  labs(title = paste("TF activity (ULM) —", prefix),
       subtitle = "DoRothEA A/B/C | Top TFs filtered by padj | positive score = more active in PR",
       x = "Transcription factor", y = "Activity score (ULM)", fill = NULL) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "top")

ggsave(paste0(output, prefix, "_TF_activity_barplot.png"), p_tf_placebo_res, width = 10, height = 8, dpi = 400)
cat("Saved TF activity barplot.\n")



#================================================================================
#        Volcano Plot — reusable for all comparisons
#================================================================================

library(ggplot2)
library(ggrepel)
library(dplyr)

make_volcano <- function(res, prefix, title_label, output_dir,
                         padj_thresh = 0.05, fc_thresh = 0.585, n_label = 5) {
  
  # 1. Classify genes into five categories
  res_filtered <- as.data.frame(res) %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
      gene = rownames(.),
      diffexpressed = case_when(
        padj <  padj_thresh & log2FoldChange >  fc_thresh ~ "UP in PR",
        padj <  padj_thresh & log2FoldChange < -fc_thresh ~ "DOWN in PR",
        padj >= padj_thresh & abs(log2FoldChange) > fc_thresh ~ "High Change, Low Sig",
        padj <  padj_thresh & abs(log2FoldChange) <= fc_thresh ~ "Low Change, High Sig",
        TRUE ~ "NS"
      ),
      diffexpressed = factor(diffexpressed,
                             levels = c("UP in PR", "DOWN in PR", "High Change, Low Sig",
                                        "Low Change, High Sig", "NS"))
    )
  
  # 2. Totals
  total_up   <- sum(res_filtered$diffexpressed == "UP in PR")
  total_down <- sum(res_filtered$diffexpressed == "DOWN in PR")
  cat(prefix, "— UP:", total_up, "| DOWN:", total_down, "\n")
  
  # 3. Top genes to label (padj ↑, then strongest log2FC)
  top_up   <- res_filtered %>% filter(diffexpressed == "UP in PR") %>%
    arrange(padj, desc(log2FoldChange)) %>% slice_head(n = n_label)
  top_down <- res_filtered %>% filter(diffexpressed == "DOWN in PR") %>%
    arrange(padj, log2FoldChange) %>% slice_head(n = n_label)
  top_genes <- bind_rows(top_up, top_down)
  
  # 4. Axis limits
  ymax <- max(-log10(res_filtered$padj), na.rm = TRUE) + 1
  xlim <- ceiling(max(abs(res_filtered$log2FoldChange), na.rm = TRUE))
  
  # 5. Plot
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
    geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = c(-fc_thresh, fc_thresh), linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    geom_point() +
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

#================================================================================
#        Run for each comparison
#================================================================================

p_volc_placebo <- make_volcano(res_placebo_res, prefix = "Placebo_RES",
                               title_label = "Placebo (PR vs PD)", output_dir = output)

# In the 6-week script:
# p_volc_6w  <- make_volcano(res_6w_res,  prefix = prefix,
#                            title_label = "6 Weeks (PR vs PD)", output_dir = output)

# In the 24h script (adjust object name if different):
# p_volc_24h <- make_volcano(res_24h_res, prefix = prefix,
#                            title_label = "24 Hours (PR vs PD)", output_dir = output)