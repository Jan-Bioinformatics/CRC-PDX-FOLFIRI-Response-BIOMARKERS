rm(list = ls(all.names = TRUE))     # Remove all objects
gc()                                # Trigger garbage collection


################################################################################
#                   DESeq2 Pipeline/ 24h TIME — RESPONSE
#                         PR vs. PD
################################################################################


library(ggrepel)
library(pheatmap)
library(DESeq2)
library(tidyverse)
library(clusterProfiler)  # For performing GSEA and enrichment analysis
library(org.Hs.eg.db)
library(enrichplot)
library(AnnotationDbi)
library(circlize)
library(ComplexHeatmap)  # Advanced and customizable heatmaps
library(gridExtra)       # Arrange multiple plots in a grid
library(knitr)           # Data table formatting
library(tidyverse)       # Data manipulation and visualization
library(ggplot2)
library(dplyr)
library(readxl)
library(msigdbr)         # For retrieving gene sets from MSigDB
library(readxl)


setwd("C:/Users/40496110/Downloads/PDX_PROJECT")



pipelines     <- "./2_pipelines/PDX_Longitudinal/LONG_RESPONSE/"
output        <- "./3_output/PDX_Longitudinal/LONG_RESPONSE/long__6w_pr_pd/"
clinicals     <- "./0_data/PDX_Longitudinal/counts/"
annotations   <- "./0_data/PDX_Longitudinal/annotations/"



################################################################################
#  BLOCK 1 — PR vs PD
################################################################################


week   <- "6w"
prefix <- "6w_RES"

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

# ---- 3. Subset to 6wks samples only ----
TRT_6w <- TRT[TRT$Treatment == "6wks", ]
cat("6wks samples total:", nrow(TRT_6w), "\n")

table(TRT_6w$Response_3_6wk, useNA = "always")

# ---- 4. Keep only PR and PD ----
TRT_6w_PRPD <- TRT_6w[TRT_6w$Response_3_6wk %in% c("PR", "PD"), ]
TRT_6w_PRPD$Response_3_6wk <- factor(TRT_6w_PRPD$Response_3_6wk, levels = c("PD", "PR"))

cat("6wks samples with PR/PD label:", nrow(TRT_6w_PRPD), "\n")
table(TRT_6w_PRPD$Response_3_6wk)


# ---- 5. Align counts to filtered metadata ----
GEX_6w_PRPD <- GEX[, rownames(TRT_6w_PRPD)]

# ---- 6. Build DESeq2 object (unpaired, same as 24h/Placebo response comparisons) ----
dds_6w_res <- DESeqDataSetFromMatrix(
  countData = GEX_6w_PRPD,
  colData   = TRT_6w_PRPD,
  design    = ~ Response_3_6wk
)
dds_6w_res <- dds_6w_res[rowSums(counts(dds_6w_res)) >= 10, ]

cat("Genes going into DESeq2:", nrow(dds_6w_res), "\n")
cat("Samples:", ncol(dds_6w_res), "\n")

# ---- 7. Run DESeq2 ----
set.seed(123)
dds_6w_res <- DESeq(dds_6w_res, fitType = "glmGamPoi")

if (!dir.exists(pipelines)) {
  dir.create(pipelines, recursive = TRUE)
}

saveRDS(dds_6w_res, paste0(pipelines, "dds_6w_res_fitted.rds"))


saveRDS(dds_6w_res, paste0(pipelines, "dds_6w_res_fitted.rds"))
cat("Fitted object saved to disk.\n")

resultsNames(dds_6w_res)


#==================================================

res_6w_res <- results(dds_6w_res, contrast = c("Response_3_6wk", "PR", "PD"))
summary(res_6w_res)

write.csv(as.data.frame(res_6w_res), paste0(pipelines, "6wks_PR_vs_PD_DESeq2_results.csv"))

vsd_6w_res <- vst(dds_6w_res, blind = FALSE)
p_pca_6w_res <- plotPCA(vsd_6w_res, intgroup = "Response_3_6wk") +
  labs(title = "PCA: 6wks PR vs PD") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
print(p_pca_6w_res)



res_6w_res <- results(dds_6w_res, contrast = c("Response_3_6wk", "PR", "PD"))
summary(res_6w_res)

write.csv(as.data.frame(res_6w_res), paste0(pipelines, "6wks_PR_vs_PD_DESeq2_results.csv"))

vsd_6w_res <- vst(dds_6w_res, blind = FALSE)
p_pca_6w_res <- plotPCA(vsd_6w_res, intgroup = "Response_3_6wk") +
  labs(title = "PCA: 6wks PR vs PD") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
print(p_pca_6w_res)


#================================================================================
#        GSEA — 6 weeks, PR vs PD (Updated Layout, Typography & Warnings Fix)
#================================================================================

library(msigdbr)
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

# ---- 1. Ranked gene list from res_6w_res ----
res_6w_res_df <- as.data.frame(res_6w_res)
res_6w_res_df <- res_6w_res_df[!is.na(res_6w_res_df$stat), ]

geneList_6w_res <- res_6w_res_df$stat
names(geneList_6w_res) <- rownames(res_6w_res_df)
geneList_6w_res <- sort(geneList_6w_res, decreasing = TRUE)

cat("Genes in ranked list:", length(geneList_6w_res), "\n")

# ---- 2. Run GSEA Calculations ----
h_gsea_6w_res <- GSEA(geneList_6w_res, exponent = 1, minGSSize = 10, maxGSSize = 10000,
                      pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = h_t2g,
                      verbose = TRUE, seed = TRUE, nPermSimple = 100000, eps = 0)
h_df_6w_res <- as.data.frame(h_gsea_6w_res)
write.csv(h_df_6w_res, paste0(output, prefix, "_Hallmarks_GSEA_RESULTS.csv"))

kegg_gsea_6w_res <- GSEA(geneList_6w_res, exponent = 1, minGSSize = 10, maxGSSize = 10000,
                         pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = kegg_t2g,
                         verbose = TRUE, seed = TRUE, nPermSimple = 100000, eps = 0)
kegg_df_6w_res <- as.data.frame(kegg_gsea_6w_res)
write.csv(kegg_df_6w_res, paste0(output, prefix, "_KEGG_GSEA_RESULTS.csv"))

c5_gsea_6w_res <- GSEA(geneList_6w_res, exponent = 1, minGSSize = 10, maxGSSize = 10000,
                       pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = c5_t2g,
                       verbose = TRUE, seed = TRUE, nPermSimple = 100000, eps = 0)
c5_df_6w_res <- as.data.frame(c5_gsea_6w_res)
write.csv(c5_df_6w_res, paste0(output, prefix, "_C5_GO_GSEA_RESULTS.csv"))


# ================================================================================
#        Plot Layout Customizations (Label Formatting & Unified Grid Styling)
# ================================================================================

bold_large_theme <- theme(
  plot.title    = element_text(hjust = 0.5, face = "bold", size = 16),
  axis.title    = element_text(size = 13, face = "bold"),
  axis.text.x   = element_text(size = 12, face = "bold"),
  axis.text.y   = element_text(size = 12, face = "bold", lineheight = 0.85),
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

# ---- 1. Filter p.adj < 0.05, sort by p.adj (small → large), then NES (large → small) ----
sort_sig <- function(df, cutoff = 0.05) {
  df %>%
    dplyr::filter(!is.na(p.adjust), p.adjust < cutoff) %>%
    dplyr::arrange(p.adjust, dplyr::desc(NES))
}

# ---- 2. From the sorted list, take the top 5 positive and top 5 negative NES ----
# Because the list is already sorted by p.adj, slice_head() picks the
# 5 MOST SIGNIFICANT pathways in each direction.
top_pos_neg <- function(sig_df, n = 5) {
  pos <- sig_df %>% dplyr::filter(NES > 0) %>%
    dplyr::slice_max(NES, n = n, with_ties = FALSE)   # 5 highest positive NES
  neg <- sig_df %>% dplyr::filter(NES < 0) %>%
    dplyr::slice_min(NES, n = n, with_ties = FALSE)   # 5 most negative NES
  dplyr::bind_rows(pos, neg) %>%
    dplyr::mutate(Direction = ifelse(NES > 0, "Activated", "Suppressed"))
}

# ---- 3. Bar plot ----
make_gsea_bar <- function(plot_df, title, wrap = 22) {
  if (nrow(plot_df) == 0) return(NULL)
  plot_df <- clean_labels(plot_df, width = wrap)
  ggplot(plot_df, aes(reorder(Description, NES), NES, fill = Direction)) +
    geom_col(width = 0.7) + coord_flip() +
    scale_fill_manual(values = c("Activated" = "#2E7D32", "Suppressed" = "#C62828")) +
    labs(x = NULL, y = "Normalized Enrichment Score", title = title) +
    theme_classic(base_size = 14) + bold_large_theme
}

# ---- Run for each collection ----
h_sig_6w_res    <- sort_sig(h_df_6w_res)
kegg_sig_6w_res <- sort_sig(kegg_df_6w_res)
c5_sig_6w_res   <- sort_sig(c5_df_6w_res)

write.csv(h_sig_6w_res,    paste0(output, prefix, "_Hallmarks_GSEA_Sorted.csv"), row.names = FALSE)
write.csv(kegg_sig_6w_res, paste0(output, prefix, "_KEGG_GSEA_Sorted.csv"),      row.names = FALSE)
write.csv(c5_sig_6w_res,   paste0(output, prefix, "_C5_GO_GSEA_Sorted.csv"),     row.names = FALSE)

p1 <- make_gsea_bar(top_pos_neg(h_sig_6w_res),    "Hallmark", wrap = 22)
p2 <- make_gsea_bar(top_pos_neg(kegg_sig_6w_res), "KEGG",     wrap = 22)
p3 <- make_gsea_bar(top_pos_neg(c5_sig_6w_res),   "GO",       wrap = 16)

# ---- Assemble (base R replacement for purrr::compact) ----
plot_list <- Filter(Negate(is.null), list(p1, p2, p3))

if (length(plot_list) > 0) {
  combined_6w_res_plot <- wrap_plots(plot_list, ncol = length(plot_list)) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title    = "GSEA — Significant Pathways: 6 Weeks (PR vs PD)",
      subtitle = "Top 5 significant activated & suppressed pathways/terms (p.adj < 0.05), ordered by NES",
      theme = theme(
        plot.title    = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 13, face = "bold", hjust = 0.5)
      )
    ) &
    theme(legend.position = "bottom")
  
  ggsave(paste0(output, prefix, "_Combined_GSEA_barplots_horizontal.png"),
         combined_6w_res_plot, width = 18, height = 11, dpi = 400)
  cat("Saved unified horizontal 6wks GSEA barplot assembly.\n")
} else {
  cat("No significant pathways crossed padj < 0.05 across Hallmark, KEGG, or GO.\n")
}





#===============================================================================
#                            ssGSEA pathway
#===============================================================================

library(GSVA)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)
library(stringr)
library(patchwork)

normCount_6w_res <- counts(dds_6w_res, normalized = TRUE)
expr_mat_6w_res  <- log2(normCount_6w_res + 1)

# ---- 1. Run ssGSEA + Wilcoxon stats, save full & sorted tables ----
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
    ) %>%
    mutate(padj = p.adjust(p_value, method = "BH"))
  
  # GO: keep only BP/CC/MF (drops HPO sets in C5)
  if (collection_label == "GO") {
    pathway_stats <- pathway_stats %>% filter(grepl("^GOBP_|^GOCC_|^GOMF_", Pathway))
  }
  
  write.csv(pathway_stats,
            paste0(output_dir, prefix, "_ssGSEA_", collection_label, "_pathway_stats.csv"),
            row.names = FALSE)
  
  # Filter padj < cutoff, sort by padj (small → large), then delta (large → small)
  sig_sorted <- pathway_stats %>%
    dplyr::filter(!is.na(padj), padj < cutoff) %>%
    dplyr::arrange(padj, dplyr::desc(delta))
  
  write.csv(sig_sorted,
            paste0(output_dir, prefix, "_ssGSEA_", collection_label, "_Sorted.csv"),
            row.names = FALSE)
  cat(collection_label, "— pathways with padj <", cutoff, ":", nrow(sig_sorted), "\n")
  
  sig_sorted
}

# ---- 2. Top 5 positive & 5 most negative delta ----
top_pos_neg_delta <- function(sig_df, n = 5) {
  pos <- sig_df %>% dplyr::filter(delta > 0) %>% dplyr::slice_max(delta, n = n, with_ties = FALSE)
  neg <- sig_df %>% dplyr::filter(delta < 0) %>% dplyr::slice_min(delta, n = n, with_ties = FALSE)
  dplyr::bind_rows(pos, neg) %>%
    dplyr::mutate(Direction   = ifelse(delta > 0, "Higher in PR", "Higher in PD"),
                  Description = Pathway)          # clean_labels() expects this column
}

# ---- 3. Bar plot (same style/colours as GSEA) ----
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

# ---- Run ----
ss_hallmark_6w <- run_ssgsea_collection(expr_mat_6w_res, h_t2g,    "Hallmark", prefix, output,
                                        TRT_6w_PRPD, "Response_3_6wk")
ss_kegg_6w     <- run_ssgsea_collection(expr_mat_6w_res, kegg_t2g, "KEGG",     prefix, output,
                                        TRT_6w_PRPD, "Response_3_6wk")
ss_go_6w       <- run_ssgsea_collection(expr_mat_6w_res, c5_t2g,   "GO",       prefix, output,
                                        TRT_6w_PRPD, "Response_3_6wk")

s1 <- make_ssgsea_bar(top_pos_neg_delta(ss_hallmark_6w), "Hallmark", wrap = 22)
s2 <- make_ssgsea_bar(top_pos_neg_delta(ss_kegg_6w),     "KEGG",     wrap = 22)
s3 <- make_ssgsea_bar(top_pos_neg_delta(ss_go_6w),       "GO",       wrap = 16)

ss_plot_list <- Filter(Negate(is.null), list(s1, s2, s3))

if (length(ss_plot_list) > 0) {
  combined_ss_plot <- wrap_plots(ss_plot_list, ncol = length(ss_plot_list)) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title    = "ssGSEA — Significant Pathways: 6 Weeks (PR vs PD)",
      subtitle = "Top 5 higher in PR & top 5 higher in PD by delta | Wilcoxon p.adj < 0.05",
      theme = theme(
        plot.title    = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 13, face = "bold", hjust = 0.5)
      )
    ) &
    theme(legend.position = "bottom")
  
  ggsave(paste0(output, prefix, "_Combined_ssGSEA_barplots_horizontal.png"),
         combined_ss_plot, width = 18, height = 11, dpi = 400)
  cat("Saved unified horizontal 6wks ssGSEA barplot assembly.\n")
} else {
  cat("No pathways crossed padj < 0.05 across Hallmark, KEGG, or GO.\n")
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

stat_vec <- res_6w_res$stat
names(stat_vec) <- rownames(res_6w_res)
stat_vec <- stat_vec[!is.na(stat_vec)]

mat <- matrix(stat_vec, ncol = 1, dimnames = list(names(stat_vec), prefix))
cat("Genes in matrix:", nrow(mat), "\n")

tf_6w_res <- decoupleR::run_ulm(
  mat     = mat,
  net     = net,
  .source = "source",
  .target = "target",
  .mor    = "mor",
  minsize = 5
)

n_sig <- tf_6w_res %>% filter(p_value < 0.05) %>% nrow()
cat("TFs tested:", nrow(tf_6w_res), "| Significant (p<0.05):", n_sig, "\n")

write.csv(tf_6w_res, paste0(output, prefix, "_TF_activity_ULM.csv"), row.names = FALSE)

tf_sig <- tf_6w_res %>% filter(p_value < 0.05) %>% arrange(desc(score))
n_show <- min(10, nrow(tf_sig))

tf_plot <- bind_rows(tf_sig %>% head(n_show), tf_sig %>% tail(n_show)) %>%
  distinct(source, .keep_all = TRUE) %>%
  mutate(
    direction = ifelse(score > 0, "Active in PR", "Active in PD"),
    source    = factor(source, levels = source[order(score)])
  )

p_tf_6w_res <- ggplot(tf_plot, aes(x = source, y = score, fill = direction)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Active in PR" = "#4DBBD5", "Active in PD" = "#F8766D")) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey30") +
  labs(title = paste("TF activity (ULM) —", prefix),
       subtitle = "DoRothEA A/B/C | positive score = more active in PR",
       x = "Transcription factor", y = "Activity score (ULM)", fill = NULL) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "top")

ggsave(paste0(output, prefix, "_TF_activity_barplot.png"), p_tf_6w_res, width = 10, height = 8, dpi = 400)
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

# p_volc_placebo <- make_volcano(res_placebo_res, prefix = "Placebo_RES",
#                               title_label = "Placebo (PR vs PD)", output_dir = output)

# In the 6-week script:
 p_volc_6w  <- make_volcano(res_6w_res,  prefix = prefix,
                            title_label = "6 Weeks (PR vs PD)", output_dir = output)

# In the 24h script (adjust object name if different):
#p_volc_24h <- make_volcano(res_24h_res, prefix = prefix,
#                           title_label = "24 Hours (PR vs PD)", output_dir = output)
