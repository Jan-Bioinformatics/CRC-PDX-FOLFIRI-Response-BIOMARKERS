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
output        <- "./3_output/PDX_Longitudinal/LONG_RESPONSE/long_24_pr_pd/"
clinicals     <- "./0_data/PDX_Longitudinal/counts/"
annotations   <- "./0_data/PDX_Longitudinal/annotations/"




################################################################################
#  BLOCK 1 — PR vs PD
################################################################################


week   <- "24h"
prefix <- "24_RES"

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

# ---- 3. Subset to 24h samples only ----
TRT_24h <- TRT[TRT$Treatment == "24h", ]
cat("24h samples total:", nrow(TRT_24h), "\n")

# ---- 4. Check Response_3_6wk distribution within 24h samples ----
table(TRT_24h$Response_3_6wk, useNA = "always")

# ---- 5. Keep only PR and PD (drop SD, NA) ----
TRT_24h_PRPD <- TRT_24h[TRT_24h$Response_3_6wk %in% c("PR", "PD"), ]
cat("24h samples with PR/PD label:", nrow(TRT_24h_PRPD), "\n")
table(TRT_24h_PRPD$Response_3_6wk)

# Set Response as a factor, PD as reference (so positive log2FC = higher in PR)
TRT_24h_PRPD$Response_3_6wk <- factor(TRT_24h_PRPD$Response_3_6wk, levels = c("PD", "PR"))

# ---- 6. Align GEX to filtered samples ----
GEX_24h_PRPD <- GEX[, rownames(TRT_24h_PRPD)]

# ---- 7. Build DESeq2 object ----
# Note: unpaired design here (~ Response only) since PR/PD is a between-patient
# comparison, not within-patient like the timepoint comparisons were.
dds_24h_res <- DESeqDataSetFromMatrix(
  countData = GEX_24h_PRPD,
  colData   = TRT_24h_PRPD,
  design    = ~ Response_3_6wk
)

dds_24h_res <- dds_24h_res[rowSums(counts(dds_24h_res)) >= 10, ]

cat("Genes going into DESeq2:", nrow(dds_24h_res), "\n")
cat("Samples:", ncol(dds_24h_res), "\n")

# ---- 8. Run DESeq2 ----
set.seed(123)
dds_24h_res <- DESeq(dds_24h_res, fitType = "glmGamPoi")


if (!dir.exists(pipelines)) {
  dir.create(pipelines, recursive = TRUE)
}

saveRDS(dds_24h_res, paste0(pipelines, "dds_24h_res_fitted.rds"))


# Save immediately, before touching anything else
saveRDS(dds_24h_res, paste0(pipelines, "dds_24h_res_fitted.rds"))
cat("Fitted object saved to disk.\n")

resultsNames(dds_24h_res)

res_24h_res <- results(dds_24h_res, contrast = c("Response_3_6wk", "PR", "PD"))
summary(res_24h_res)

write.csv(as.data.frame(res_24h_res), paste0(pipelines, "24h_PR_vs_PD_DESeq2_results.csv"))
cat("Saved 24h PR vs PD results.\n")


#===========================================================

# ---- Check Cook's distances per sample ----
cooks_mat <- assays(dds_24h_res)[["cooks"]]

# Average Cook's distance per sample (higher = more influential/outlier-prone)
cooks_per_sample <- apply(cooks_mat, 2, mean, na.rm = TRUE)
cooks_per_sample_df <- data.frame(
  SampleID = names(cooks_per_sample),
  MeanCooks = cooks_per_sample,
  Response = TRT_24h_PRPD[names(cooks_per_sample), "Response_3_6wk"]
)

cooks_per_sample_df[order(-cooks_per_sample_df$MeanCooks), ]

vsd_24h_res <- vst(dds_24h_res, blind = FALSE)
p_pca_24h_res <- plotPCA(vsd_24h_res, intgroup = "Response_3_6wk") +
  labs(title = "PCA: 24h PR vs PD") + theme_classic()
print(p_pca_24h_res)

# Save the PCA plot to your output folder
ggsave(
  filename = paste0(output, prefix, "_PCA_24h_PR_vs_PD.png"), 
  plot     = p_pca_24h_res, 
  width    = 8, 
  height   = 6, 
  dpi      = 400
)
cat("Saved PCA plot to output directory.\n")
#================================================================================
#================================================================================

library(msigdbr)
library(dplyr)

# Hallmark
h_t2g <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, gene_symbol)

# KEGG
kegg_t2g <- msigdbr(species = "Homo sapiens", collection = "C2",
                    subcollection = "CP:KEGG_LEGACY") %>%
  dplyr::select(gs_name, gene_symbol)

# GO (C5)
c5_t2g <- msigdbr(species = "Homo sapiens", collection = "C5") %>%
  dplyr::select(gs_name, gene_symbol)

cat("Hallmark sets:", length(unique(h_t2g$gs_name)), "\n")
cat("KEGG sets:", length(unique(kegg_t2g$gs_name)), "\n")
cat("GO sets:", length(unique(c5_t2g$gs_name)), "\n")

################################################################################
#  24h PR vs PD — Standard GSEA (Hallmark / KEGG / GO)
################################################################################

# ---- 1. Ranked gene list from res_24h_res ----
res_24h_res_df <- as.data.frame(res_24h_res)
res_24h_res_df <- res_24h_res_df[!is.na(res_24h_res_df$stat), ]

geneList_24h_res <- res_24h_res_df$stat
names(geneList_24h_res) <- rownames(res_24h_res_df)
geneList_24h_res <- sort(geneList_24h_res, decreasing = TRUE)

cat("Genes in ranked list:", length(geneList_24h_res), "\n")

# ---- 2. Hallmark GSEA ----
h_gsea_24h_res <- GSEA(geneList_24h_res, exponent = 1, minGSSize = 10, maxGSSize = 10000,
                       pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = h_t2g,
                       verbose = TRUE, seed = TRUE, nPermSimple = 10000)
h_df_24h_res <- as.data.frame(h_gsea_24h_res)
write.csv(h_df_24h_res, paste0(output, prefix, "_Hallmarks_GSEA_RESULTS.csv"))
cat("Hallmark tested:", nrow(h_df_24h_res), "| Significant:", sum(h_df_24h_res$p.adjust < 0.05), "\n")

h_sig_24h_res <- h_df_24h_res[h_df_24h_res$p.adjust < 0.05, ]
h_sig_24h_res <- h_sig_24h_res[order(-h_sig_24h_res$NES), ]
write.csv(h_sig_24h_res, paste0(output, prefix, "_Hallmarks_GSEA_Sorted.csv"))

if (nrow(h_sig_24h_res) > 0) {
  # Extract top 5 positive and top 5 negative rows
  h_top5_pos <- head(h_sig_24h_res, 5)
  h_top5_neg <- tail(h_sig_24h_res, 5)
  h_plot_data <- unique(rbind(h_top5_pos, h_top5_neg)) # unique() avoids duplicates if < 10 paths total
  
  # Create a direction column for plotting color
  h_plot_data$Direction <- ifelse(h_plot_data$NES > 0, "Activated", "Suppressed")
  
  p_hallmark_24h_res <- ggplot(h_plot_data, aes(reorder(Description, NES), NES, fill = Direction)) +
    geom_col(width = 0.7) + coord_flip() +
    scale_fill_manual(values = c("Activated" = "#2E7D32", "Suppressed" = "#C62828")) +
    labs(x = "Pathway", y = "Normalized Enrichment Score",
         title = "Top 5 Positive & Negative Significant Hallmark Pathways") +
    theme_classic(base_size = 14) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title = element_text(size = 12),
          legend.position = "right")
  
  ggsave(paste0(output, prefix, "_Hallmarks_GSEA_barplot.png"),
         p_hallmark_24h_res, width = 11, height = 5, dpi = 400)
  cat("Saved Hallmark barplot.\n")
}

# ---- 3. KEGG GSEA ----
kegg_gsea_24h_res <- GSEA(geneList_24h_res, exponent = 1, minGSSize = 10, maxGSSize = 10000,
                          pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = kegg_t2g,
                          verbose = TRUE, seed = TRUE, nPermSimple = 10000)
kegg_df_24h_res <- as.data.frame(kegg_gsea_24h_res)
write.csv(kegg_df_24h_res, paste0(output, prefix, "_KEGG_GSEA_RESULTS.csv"))
cat("KEGG tested:", nrow(kegg_df_24h_res), "| Significant:", sum(kegg_df_24h_res$p.adjust < 0.05), "\n")

kegg_sig_24h_res <- kegg_df_24h_res[kegg_df_24h_res$p.adjust < 0.05, ]
kegg_sig_24h_res <- kegg_sig_24h_res[order(-kegg_sig_24h_res$NES), ]
write.csv(kegg_sig_24h_res, paste0(output, prefix, "_KEGG_GSEA_Sorted.csv"))

if (nrow(kegg_sig_24h_res) > 0) {
  # Extract top 5 positive and top 5 negative rows
  kegg_top5_pos <- head(kegg_sig_24h_res, 5)
  kegg_top5_neg <- tail(kegg_sig_24h_res, 5)
  kegg_plot_data <- unique(rbind(kegg_top5_pos, kegg_top5_neg))
  
  kegg_plot_data$Direction <- ifelse(kegg_plot_data$NES > 0, "Activated", "Suppressed")
  
  p_kegg_24h_res <- ggplot(kegg_plot_data, aes(reorder(Description, NES), NES, fill = Direction)) +
    geom_col(width = 0.7) + coord_flip() +
    scale_fill_manual(values = c("Activated" = "#2E7D32", "Suppressed" = "#C62828")) +
    labs(x = "Pathway", y = "Normalized Enrichment Score",
         title = "Top 5 Positive & Negative Significant KEGG Pathways") +
    theme_classic(base_size = 14) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title = element_text(size = 12),
          legend.position = "right")
  
  ggsave(paste0(output, prefix, "_KEGG_GSEA_barplot.png"),
         p_kegg_24h_res, width = 11, height = 5, dpi = 400)
  cat("Saved KEGG barplot.\n")
}

# ---- 4. GO (C5) GSEA ----
c5_gsea_24h_res <- GSEA(geneList_24h_res, exponent = 1, minGSSize = 10, maxGSSize = 10000,
                        pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = c5_t2g,
                        verbose = TRUE, seed = TRUE, nPermSimple = 10000)
c5_df_24h_res <- as.data.frame(c5_gsea_24h_res)
write.csv(c5_df_24h_res, paste0(output, prefix, "_C5_GO_GSEA_RESULTS.csv"))
cat("GO tested:", nrow(c5_df_24h_res), "| Significant:", sum(c5_df_24h_res$p.adjust < 0.05), "\n")

c5_sig_24h_res <- c5_df_24h_res[c5_df_24h_res$p.adjust < 0.05, ]
c5_sig_24h_res <- c5_sig_24h_res[order(-c5_sig_24h_res$NES), ]
write.csv(c5_sig_24h_res, paste0(output, prefix, "_C5_GO_GSEA_Sorted.csv"))

if (nrow(c5_sig_24h_res) > 0) {
  # Extract top 5 positive and top 5 negative rows
  c5_top5_pos <- head(c5_sig_24h_res, 5)
  c5_top5_neg <- tail(c5_sig_24h_res, 5) # Fixed assignment operator here
  c5_plot_data <- unique(rbind(c5_top5_pos, c5_top5_neg))
  
  c5_plot_data$Direction <- ifelse(c5_plot_data$NES > 0, "Activated", "Suppressed")
  
  p_go_24h_res <- ggplot(c5_plot_data, aes(reorder(Description, NES), NES, fill = Direction)) +
    geom_col(width = 0.7) + coord_flip() +
    scale_fill_manual(values = c("Activated" = "#2E7D32", "Suppressed" = "#C62828")) +
    labs(x = "Pathway", y = "Normalized Enrichment Score",
         title = "Top 5 Positive & Negative Significant GO Pathways") +
    theme_classic(base_size = 14) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title = element_text(size = 12),
          legend.position = "right")
  
  ggsave(paste0(output, prefix, "_C5_GO_GSEA_barplot.png"),
         p_go_24h_res, width = 12, height = 5, dpi = 400)
  cat("Saved GO barplot.\n")
}


library(patchwork)
library(ggplot2)
library(stringr)

# --- Bold, large-font theme ---
bold_large_theme <- theme(
  plot.title   = element_text(hjust = 0.5, face = "bold", size = 18),
  axis.title   = element_text(size = 14, face = "bold"),
  axis.text.x  = element_text(size = 14, face = "bold"),
  axis.text.y  = element_text(size = 12, face = "bold"),
  legend.title = element_text(size = 13, face = "bold"),
  legend.text  = element_text(size = 12, face = "bold")
)

# --- Helper: strip redundant prefix, tidy underscores, wrap onto multiple lines ---
clean_labels <- function(df, width = 22) {
  df$Description <- df$Description %>%
    str_remove("^HALLMARK_") %>%
    str_remove("^KEGG_") %>%
    str_remove("^GOBP_") %>%
    str_remove("^GOCC_") %>%
    str_remove("^GOMF_") %>%
    str_replace_all("_", " ") %>%
    str_to_title() %>%
    str_wrap(width = width)
  df
}

# Hallmark/KEGG names are shorter -> width 22; GO names run longer -> wrap tighter (width 16)
h_plot_data    <- clean_labels(h_plot_data,    width = 22)
kegg_plot_data <- clean_labels(kegg_plot_data, width = 22)
c5_plot_data   <- clean_labels(c5_plot_data,   width = 16)

# --- Rebuild each plot with cleaned/wrapped labels + bold theme ---
p1 <- ggplot(h_plot_data, aes(reorder(Description, NES), NES, fill = Direction)) +
  geom_col(width = 0.7) + coord_flip() +
  scale_fill_manual(values = c("Activated" = "#2E7D32", "Suppressed" = "#C62828")) +
  labs(x = NULL, y = "Normalized Enrichment Score", title = "Hallmark") +
  theme_classic(base_size = 14) + bold_large_theme

p2 <- ggplot(kegg_plot_data, aes(reorder(Description, NES), NES, fill = Direction)) +
  geom_col(width = 0.7) + coord_flip() +
  scale_fill_manual(values = c("Activated" = "#2E7D32", "Suppressed" = "#C62828")) +
  labs(x = NULL, y = "Normalized Enrichment Score", title = "KEGG") +
  theme_classic(base_size = 14) + bold_large_theme

p3 <- ggplot(c5_plot_data, aes(reorder(Description, NES), NES, fill = Direction)) +
  geom_col(width = 0.7) + coord_flip() +
  scale_fill_manual(values = c("Activated" = "#2E7D32", "Suppressed" = "#C62828")) +
  labs(x = NULL, y = "Normalized Enrichment Score", title = "GO") +
  theme_classic(base_size = 14) + bold_large_theme

# --- Combine side-by-side, shared legend at bottom, bold overall title ---
combined_24h_res_plot <- (p1 | p2 | p3) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "GSEA — significant pathways: 24h (PR vs PD)",
    subtitle = "ordered by NES, p.adj < 0.05",
    theme = theme(
      plot.title    = element_text(size = 20, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  ) &
  theme(legend.position = "bottom")

ggsave(
  filename = paste0(output, prefix, "_Combined_GSEA_barplots_horizontal.png"),
  plot     = combined_24h_res_plot,
  width    = 18,
  height   = 11,   # slightly taller to accommodate 2-line GO labels
  dpi      = 400
)

cat("Saved horizontal combined 24h PR vs PD GSEA barplot.\n")

#===============================================================================
#                            ssGSEA pathway
#===============================================================================

library(GSVA)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)
library(patchwork)
library(stringr)

# ---- Normalized counts, log-transformed (reuse across all 3 collections) ----
normCount_24h_res <- counts(dds_24h_res, normalized = TRUE)
expr_mat_24h_res   <- log2(normCount_24h_res + 1)

# ---- Reusable function: run ssGSEA + Wilcoxon test + barplot for one collection ----
run_ssgsea_collection <- function(expr_mat, t2g, collection_label, prefix, output_dir,
                                  TRT_sub, response_col, name_strip = NULL) {
  
  gene_sets <- split(t2g$gene_symbol, t2g$gs_name)
  
  set.seed(123)
  ssgsea_param <- ssgseaParam(
    exprData  = expr_mat,
    geneSets  = gene_sets,
    normalize = TRUE,
    minSize   = 10,
    maxSize   = 10000
  )
  ssgsea_result <- gsva(ssgsea_param, verbose = TRUE)
  ssgsea_df <- as.data.frame(ssgsea_result)
  write.csv(ssgsea_df, paste0(output_dir, prefix, "_ssGSEA_", collection_label, "_scores.csv"))
  
  ssgsea_long <- ssgsea_df %>%
    rownames_to_column("Pathway") %>%
    pivot_longer(-Pathway, names_to = "SampleID", values_to = "ssGSEA_score") %>%
    left_join(
      TRT_sub %>% rownames_to_column("SampleID") %>% dplyr::select(SampleID, all_of(response_col)),
      by = "SampleID"
    )
  
  pathway_stats <- ssgsea_long %>%
    group_by(Pathway) %>%
    summarise(
      mean_PR = mean(ssGSEA_score[.data[[response_col]] == "PR"], na.rm = TRUE),
      mean_PD = mean(ssGSEA_score[.data[[response_col]] == "PD"], na.rm = TRUE),
      delta   = mean_PR - mean_PD,
      p_value = wilcox.test(
        ssGSEA_score[.data[[response_col]] == "PR"],
        ssGSEA_score[.data[[response_col]] == "PD"]
      )$p.value,
      .groups = "drop"
    ) %>%
    mutate(padj = p.adjust(p_value, method = "BH"), sig = padj < 0.05) %>%
    arrange(padj)
  
  # ---- Full stats table (every pathway tested) ----
  write.csv(pathway_stats, paste0(output_dir, prefix, "_ssGSEA_", collection_label, "_pathway_stats.csv"),
            row.names = FALSE)
  cat(collection_label, "— pathways with raw p<0.05:", sum(pathway_stats$p_value < 0.05), "\n")
  
  # ---- Significant-only sorted table, for the thesis (mirrors GSEA's *_Sorted.csv) ----
  sig_stats <- pathway_stats %>% filter(sig) %>% arrange(padj)
  write.csv(sig_stats, paste0(output_dir, prefix, "_ssGSEA_", collection_label, "_Sig_Sorted.csv"),
            row.names = FALSE)
  
  # ---- Filter to significant pathways for plotting ----
  plot_data <- pathway_stats %>% filter(p_value < 0.05)
  
  # ---- for GO collection, keep only true GO terms (drop HP/other mixed-in terms) ----
  if (collection_label == "GO") {
    plot_data <- plot_data %>% filter(grepl("^GOBP_|^GOCC_|^GOMF_", Pathway))
  }
  
  # ---- Keep only top 5 positive (PR-enriched) and top 5 negative (PD-enriched) by delta ----
  plot_data <- plot_data %>% arrange(desc(delta))
  plot_top5_pos <- head(plot_data, 5)
  plot_top5_neg <- tail(plot_data, 5)
  plot_data <- unique(rbind(plot_top5_pos, plot_top5_neg))
  
  # ---- Clean pathway names for display ----
  plot_data <- plot_data %>%
    mutate(Pathway_Clean = if (!is.null(name_strip)) gsub(name_strip, "", Pathway) else Pathway,
           Pathway_Clean = gsub("^GOBP_|^GOCC_|^GOMF_", "", Pathway_Clean),
           Pathway_Clean = gsub("_", " ", Pathway_Clean),
           Pathway_Clean = tools::toTitleCase(tolower(Pathway_Clean))) %>%
    arrange(delta) %>%
    mutate(Pathway_Clean = factor(Pathway_Clean, levels = Pathway_Clean))
  
  if (nrow(plot_data) > 0) {
    n_paths <- nrow(plot_data)
    
    # ---- Dynamic sizing so labels get room to breathe ----
    plot_height <- max(5, n_paths * 0.35)
    label_size  <- if (n_paths > 20) 7.5 else if (n_paths > 10) 9 else 10
    
    p <- ggplot(plot_data, aes(x = delta, y = Pathway_Clean, fill = delta > 0)) +
      geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.3) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
      scale_fill_manual(values = c("TRUE" = "#2E7D32", "FALSE" = "#C62828"),
                        labels = c("TRUE" = "Enriched in PR", "FALSE" = "Enriched in PD")) +
      theme_bw(base_size = 12) +
      labs(title = paste0(collection_label, " Pathway: 24h "),
           subtitle = paste("Top", n_paths, "most significant pathways | ssGSEA Wilcoxon p < 0.05"),
           x = "Delta Enrichment Score (Mean PR - Mean PD)", y = NULL, fill = "Enrichment") +
      theme(
        plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
        plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "black"),
        axis.text.y   = element_text(size = label_size),
        axis.text.x   = element_text(size = 9),
        legend.position = "bottom",
        plot.margin   = margin(t = 10, r = 20, b = 10, l = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    # ---- Stash the plot object globally so it can be combined later ----
    assign(paste0("p_ssgsea_", tolower(collection_label)), p, envir = .GlobalEnv)
    
    ggsave(paste0(output_dir, prefix, "_ssGSEA_", collection_label, "_Barplot.png"),
           plot = p, width = 10, height = plot_height, dpi = 400, limitsize = FALSE)
    cat("Saved", collection_label, "ssGSEA barplot (", n_paths, "pathways ).\n")
  } else {
    cat("No pathways passed p<0.05 for", collection_label, "barplot.\n")
  }
  
  return(pathway_stats)
}

# ---- Run for all three collections ----
ss_hallmark_24h <- run_ssgsea_collection(expr_mat_24h_res, h_t2g, "Hallmark", prefix, output,
                                         TRT_24h_PRPD, "Response_3_6wk", name_strip = "^HALLMARK_")

ss_kegg_24h <- run_ssgsea_collection(expr_mat_24h_res, kegg_t2g, "KEGG", prefix, output,
                                     TRT_24h_PRPD, "Response_3_6wk", name_strip = "^KEGG_")

ss_go_24h <- run_ssgsea_collection(expr_mat_24h_res, c5_t2g, "GO", prefix, output,
                                   TRT_24h_PRPD, "Response_3_6wk")


#===============================================================================
#                    Combine the three ssGSEA barplots into one figure
#===============================================================================

# ---- Bold, large-font theme ----
bold_large_theme <- theme(
  plot.title    = element_text(face = "bold", hjust = 0.5, size = 16),
  plot.subtitle = element_text(hjust = 0.5, size = 13, face = "bold", colour = "black"),
  axis.title    = element_text(size = 13, face = "bold"),
  axis.text.x   = element_text(size = 13, face = "bold"),
  axis.text.y   = element_text(size = 13, face = "bold"),
  legend.title  = element_text(size = 12, face = "bold"),
  legend.text   = element_text(size = 11, face = "bold")
)

# ---- Wrap y-axis labels without touching the underlying data ----
p1 <- p_ssgsea_hallmark + bold_large_theme +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 22))

p2 <- p_ssgsea_kegg + bold_large_theme +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 22))

p3 <- p_ssgsea_go + bold_large_theme +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 16))   # GO terms run longer

# ---- Combine side-by-side, shared legend, bold overall title ----
combined_ssgsea_plot <- (p1 | p2 | p3) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = paste0("ssGSEA Differential Enrichment Pathway: 24h"),
    subtitle = "Wilcoxon p < 0.05, mean PR vs mean PD",
    theme = theme(
      plot.title    = element_text(size = 20, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  ) &
  theme(legend.position = "bottom")

ggsave(
  paste0(output, prefix, "_Combined_ssGSEA_barplots.png"),
  plot   = combined_ssgsea_plot,
  width  = 22,
  height = 13,
  dpi    = 400,
  limitsize = FALSE
)

cat("Saved combined ssGSEA barplot.\n")

#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
# library(decoupleR)
library(dorothea)
library(dplyr)
library(ggplot2)

# ---- 1. Load DoRothEA network ----
data("dorothea_hs", package = "dorothea")

net <- dorothea_hs %>%
  filter(confidence %in% c("A", "B", "C")) %>%
  rename(source = tf) %>%
  mutate(weight = as.numeric(mor))

cat("DoRothEA network rows (A/B/C confidence):", nrow(net), "\n")

# ---- 2. Build stat matrix from res_24h_res ----
stat_vec <- res_24h_res$stat
names(stat_vec) <- rownames(res_24h_res)
stat_vec <- stat_vec[!is.na(stat_vec)]

mat <- matrix(stat_vec, ncol = 1, dimnames = list(names(stat_vec), prefix))
cat("Genes in matrix:", nrow(mat), "\n")

# ---- 3. Run ULM (TF activity inference) ----
tf_24h_res <- decoupleR::run_ulm(
  mat     = mat,
  net     = net,
  .source = "source",
  .target = "target",
  .mor    = "mor",
  minsize = 5
)

# Calculate adjusted p-values using Benjamini-Hochberg method
tf_24h_res <- tf_24h_res %>% 
  mutate(padj = p.adjust(p_value, method = "BH"))

n_sig <- tf_24h_res %>% filter(padj < 0.05) %>% nrow()
cat("TFs tested:", nrow(tf_24h_res), "| Significant (padj<0.05):", n_sig, "\n")

write.csv(tf_24h_res, paste0(output, prefix, "_TF_activity_ULM.csv"), row.names = FALSE)

# ---- 4. Barplot: top 10 active in PR, top 10 active in PD ----
# Filter strictly by padj instead of raw p_value
tf_sig <- tf_24h_res %>% filter(padj < 0.05) %>% arrange(desc(score))

library(ggplot2)
library(stringr)

if (nrow(tf_sig) > 0) {
  n_show <- min(10, nrow(tf_sig))
  
  tf_plot <- bind_rows(tf_sig %>% head(n_show), tf_sig %>% tail(n_show)) %>%
    distinct(source, .keep_all = TRUE) %>%
    mutate(
      # Clean up names: remove common tags, change underscores to spaces
      source_clean = str_replace_all(source, "^(TF_|HALLMARK_|KEGG_|GO_)", ""),
      source_clean = str_replace_all(source_clean, "_", " "),
      # Wrap text into 2 lines if it exceeds ~15 characters
      source_clean = str_wrap(source_clean, width = 15),
      
      direction = ifelse(score > 0, "Active in PR", "Active in PD"),
      source_clean = factor(source_clean, levels = source_clean[order(score)])
    )
  
  p_tf_24h_res <- ggplot(tf_plot, aes(x = source_clean, y = score, fill = direction)) +
    geom_col(width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c("Active in PR" = "#4DBBD5", "Active in PD" = "#F8766D")) +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey30") +
    labs(title = "TF activity (ULM) — 24h",
         subtitle = "DoRothEA A/B/C | positive score = more active in PR | padj < 0.05",
         x = "Transcription factor", y = "Activity score (ULM)", fill = NULL) +
    theme_classic(base_size = 14) + # Increased overall baseline font size
    theme(plot.title      = element_text(face = "bold", size = 18, hjust = 0.5),
          plot.subtitle   = element_text(size = 13, hjust = 0.5),
          axis.title      = element_text(face = "bold", size = 15),
          axis.text.x     = element_text(size = 14),
          axis.text.y     = element_text(size = 14, lineheight = 0.85), # Adjusted line height for wrapped text
          legend.text     = element_text(size = 13),
          legend.position = "top")
  
  ggsave(paste0(output, prefix, "_TF_activity_barplot.png"), p_tf_24h_res, width = 10, height = 8, dpi = 400)
  cat("Saved TF activity barplot based on padj < 0.05.\n")
} else {
  cat("No significant TFs found at padj < 0.05. Barplot was not generated.\n")
}




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
# p_volc_6w  <- make_volcano(res_6w_res,  prefix = prefix,
#                            title_label = "6 Weeks (PR vs PD)", output_dir = output)

# In the 24h script (adjust object name if different):
 p_volc_24h <- make_volcano(res_24h_res, prefix = prefix,
                            title_label = "24 Hours (PR vs PD)", output_dir = output)