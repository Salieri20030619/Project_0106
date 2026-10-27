
cat("\n========================================\n")
cat("20. Abca1 KO 下游基因 GO + KEGG 富集分析\n")
cat("========================================\n\n")

rm(list = ls()); gc()
options(bitmapType = "cairo")

BASE_DIR   <- ""
INPUT_DIR  <- file.path(BASE_DIR, "19_scTenifoldKnk/Results")
OUTPUT_DIR <- file.path(BASE_DIR, "20_downstream_enrichment")

setwd(BASE_DIR)

res.path <- file.path(OUTPUT_DIR, "Results")
fig.path <- file.path(OUTPUT_DIR, "Figures")

if (dir.exists(OUTPUT_DIR)) {
  cat("清理旧结果目录...\n")
  unlink(OUTPUT_DIR, recursive = TRUE)
}

for (p in c(OUTPUT_DIR, res.path, fig.path)) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

BIOMARKER <- "Abca1"
INPUT_FILE <- file.path(INPUT_DIR, "Abca1_significant_P005.csv")
TOP_N <- 10

cat("\n========== 加载 R 包 ==========\n")

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(DOSE)
})

cat("✓ clusterProfiler:", as.character(packageVersion("clusterProfiler")), "\n")
cat("✓ org.Mm.eg.db:", as.character(packageVersion("org.Mm.eg.db")), "\n")

t_total <- Sys.time()

cat("\n========== 1. 加载 Abca1 KO 下游基因 ==========\n")

gene_df <- read.csv(INPUT_FILE, stringsAsFactors = FALSE)
cat("输入文件:", INPUT_FILE, "\n")
cat("总基因数:", nrow(gene_df), "\n")

genes_downstream <- setdiff(gene_df$gene, BIOMARKER)
cat("去除 Abca1 自身后下游基因数:", length(genes_downstream), "\n")
cat("前 20 个下游基因:", paste(head(genes_downstream, 20), collapse = ", "), "\n")

cat("\n========== 2. 基因名转换 ==========\n")

gene_entrez <- suppressWarnings(
  bitr(genes_downstream,
       fromType = "SYMBOL",
       toType   = "ENTREZID",
       OrgDb    = org.Mm.eg.db)
)
cat("成功转换:", nrow(gene_entrez), "/", length(genes_downstream), "\n")

if (nrow(gene_entrez) < length(genes_downstream)) {
  missing <- setdiff(genes_downstream, gene_entrez$SYMBOL)
  cat("未匹配基因 (", length(missing), "):", paste(missing, collapse = ", "), "\n")
}

entrez_ids <- gene_entrez$ENTREZID

cat("\n========== 3. GO 富集分析 ==========\n")

go_results <- list()

for (ont in c("BP", "CC", "MF")) {
  cat("\n--- GO", ont, "---\n")

  go_res <- enrichGO(
    gene          = entrez_ids,
    OrgDb         = org.Mm.eg.db,
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )

  n_sig <- ifelse(is.null(go_res), 0, nrow(as.data.frame(go_res)))
  cat("adj.P<0.05 显著通路:", n_sig, "\n")

  if (n_sig < 3) {
    cat("通路数太少, 按方案放宽到 P<0.05\n")
    go_res <- enrichGO(
      gene          = entrez_ids,
      OrgDb         = org.Mm.eg.db,
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 1,
      readable      = TRUE
    )
    n_sig <- ifelse(is.null(go_res), 0, nrow(as.data.frame(go_res)))
    cat("P<0.05 显著通路:", n_sig, "\n")
  }

  go_results[[ont]] <- go_res

  if (n_sig > 0) {
    go_df <- as.data.frame(go_res)
    write.csv(go_df,
              file.path(res.path, paste0(BIOMARKER, "_GO_", ont, ".csv")),
              row.names = FALSE)
    cat("\nTop 5 通路:\n")
    print(head(go_df[, c("Description", "pvalue", "p.adjust", "Count")], 5))
  }
}

cat("\n========== 4. KEGG 富集分析 ==========\n")

kegg_res <- tryCatch({
  enrichKEGG(
    gene          = entrez_ids,
    organism      = "mmu",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
}, error = function(e) {
  cat("KEGG 富集出错 (可能网络问题):", conditionMessage(e), "\n")
  return(NULL)
})

n_kegg <- ifelse(is.null(kegg_res), 0, nrow(as.data.frame(kegg_res)))
cat("adj.P<0.05 KEGG 通路:", n_kegg, "\n")

if (!is.null(kegg_res) && n_kegg < 3) {
  cat("通路数太少, 按方案放宽到 P<0.05\n")
  kegg_res <- tryCatch({
    enrichKEGG(
      gene          = entrez_ids,
      organism      = "mmu",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 1
    )
  }, error = function(e) NULL)
  n_kegg <- ifelse(is.null(kegg_res), 0, nrow(as.data.frame(kegg_res)))
  cat("P<0.05 KEGG 通路:", n_kegg, "\n")
}

if (!is.null(kegg_res) && n_kegg > 0) {
  kegg_res <- setReadable(kegg_res, OrgDb = org.Mm.eg.db,
                           keyType = "ENTREZID")
  kegg_df <- as.data.frame(kegg_res)
  write.csv(kegg_df,
            file.path(res.path, paste0(BIOMARKER, "_KEGG.csv")),
            row.names = FALSE)
  cat("\nTop 5 KEGG 通路:\n")
  print(head(kegg_df[, c("Description", "pvalue", "p.adjust", "Count")], 5))
}

cat("\n========== 5. 画图 ==========\n")

plot_go_combined <- function(go_results) {
  df_list <- list()
  for (ont in c("BP", "CC", "MF")) {
    go_res <- go_results[[ont]]
    if (!is.null(go_res) && nrow(as.data.frame(go_res)) > 0) {
      tmp <- as.data.frame(go_res)
      tmp <- tmp[order(tmp$pvalue), ]
      top_n <- min(TOP_N, nrow(tmp))
      tmp <- tmp[1:top_n, ]
      tmp$Ontology <- ont
      df_list[[ont]] <- tmp
    }
  }

  if (length(df_list) == 0) return(NULL)

  go_combined <- do.call(rbind, df_list)
  go_combined$Description <- factor(
    go_combined$Description,
    levels = rev(unique(go_combined$Description))
  )
  go_combined$logP <- -log10(go_combined$pvalue)

  ont_colors <- c("BP" = "#4DBBD5", "CC" = "#00A087", "MF" = "#F39B7F")

  ggplot(go_combined,
         aes(x = logP, y = Description, fill = Ontology)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = ont_colors) +
    facet_grid(Ontology ~ ., scales = "free_y", space = "free_y") +
    labs(title = paste0("Fig. 17 (", BIOMARKER,
                        "): GO enrichment (Top ", TOP_N, " per ontology)"),
         x = "-log10(P-value)", y = "") +
    theme_bw(base_size = 11) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 13),
      axis.text.y     = element_text(size = 9),
      strip.background = element_rect(fill = "grey90"),
      strip.text      = element_text(face = "bold", size = 11),
      legend.position = "none"
    )
}

p_go <- plot_go_combined(go_results)
if (!is.null(p_go)) {
  ggsave(file.path(fig.path, "Fig17a_Abca1_GO.pdf"),
         p_go, width = 11, height = 9)
  ggsave(file.path(fig.path, "Fig17a_Abca1_GO.png"),
         p_go, width = 11, height = 9, dpi = 300, bg = "white")
  cat("✓ Fig 17a (GO) 保存\n")
} else {
  cat("⚠ 无显著 GO 通路\n")
}

if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
  kegg_df <- as.data.frame(kegg_res)
  kegg_df <- kegg_df[order(kegg_df$pvalue), ]
  top_n <- min(TOP_N, nrow(kegg_df))
  kegg_df <- kegg_df[1:top_n, ]
  kegg_df$Description <- factor(
    kegg_df$Description,
    levels = rev(kegg_df$Description)
  )
  kegg_df$logP <- -log10(kegg_df$pvalue)

  p_kegg <- ggplot(kegg_df,
                   aes(x = logP, y = Description, fill = logP)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_gradient(low = "#FBDDDD", high = "#E64B35") +
    labs(title = paste0("Fig. 17 (", BIOMARKER,
                        "): KEGG enrichment (Top ", top_n, ")"),
         x = "-log10(P-value)", y = "", fill = "-log10(P)") +
    theme_bw(base_size = 11) +
    theme(
      plot.title  = element_text(hjust = 0.5, face = "bold", size = 13),
      axis.text.y = element_text(size = 10)
    )

  ggsave(file.path(fig.path, "Fig17b_Abca1_KEGG.pdf"),
         p_kegg, width = 10, height = max(5, top_n * 0.4))
  ggsave(file.path(fig.path, "Fig17b_Abca1_KEGG.png"),
         p_kegg, width = 10, height = max(5, top_n * 0.4),
         dpi = 300, bg = "white")
  cat("✓ Fig 17b (KEGG) 保存\n")
} else {
  cat("⚠ 无显著 KEGG 通路\n")
}

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("       第 20 部分 GO + KEGG 富集分析 完成\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

total_min <- as.numeric(difftime(Sys.time(), t_total, units = "mins"))
cat(sprintf("\n总耗时: %.1f 分钟\n", total_min))

summary_df <- data.frame(
  Analysis  = c("GO_BP", "GO_CC", "GO_MF", "KEGG"),
  N_Pathways = c(
    ifelse(is.null(go_results$BP), 0, nrow(as.data.frame(go_results$BP))),
    ifelse(is.null(go_results$CC), 0, nrow(as.data.frame(go_results$CC))),
    ifelse(is.null(go_results$MF), 0, nrow(as.data.frame(go_results$MF))),
    n_kegg
  )
)
write.csv(summary_df,
          file.path(res.path, "enrichment_summary.csv"),
          row.names = FALSE)

cat("\n富集结果汇总:\n")
print(summary_df)

cat("\n输出文件:\n")
cat("  Results/\n")
cat("    Abca1_GO_{BP,CC,MF}.csv     GO 三本体结果表\n")
cat("    Abca1_KEGG.csv              KEGG 结果表\n")
cat("    enrichment_summary.csv      汇总\n")
cat("  Figures/\n")
cat("    Fig17a_Abca1_GO.pdf/png     GO 组合图 (BP/CC/MF Top 10)\n")
cat("    Fig17b_Abca1_KEGG.pdf/png   KEGG Top 10\n")

cat("\n✅ 第 20 部分完成!\n")

sink(file.path(res.path, "session_info.txt"))
sessionInfo()
sink()
