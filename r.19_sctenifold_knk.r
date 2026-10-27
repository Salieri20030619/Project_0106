
rm(list = ls()); gc()

library(dplyr)
library(Seurat)
library(scTenifoldKnk)
library(ggplot2)
library(ggrepel)

options(bitmapType = "cairo")

setwd("~")

OUT_DIR <- "./19_scTenifoldKnk"
RES_DIR <- file.path(OUT_DIR, "Results")
FIG_DIR <- file.path(OUT_DIR, "Figures")
for (p in c(OUT_DIR, RES_DIR, FIG_DIR)) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE)
}

cat("========== 加载数据 ==========\n")

INPUT_RDS <- "./seurat_obj_annotated_clean.rds"

if (!file.exists(INPUT_RDS)) {
  stop("请先从 km112 传 seurat_obj_annotated_clean.rds 到当前目录")
}

seu_obj <- readRDS(INPUT_RDS)
cat("细胞数:", ncol(seu_obj), "\n")
cat("基因数:", nrow(seu_obj), "\n")
cat("细胞类型:\n"); print(table(seu_obj$cell_type))

run_scTenifoldKnk <- function(
    seu_obj,
    genes_to_analyze,
    nfeatures = 5000,
    qc = FALSE,
    prefix = "scTenifoldKnk",
    save_dir = "."
) {

  cat("工作目录:", getwd(), "\n")
  cat("细胞数:", ncol(seu_obj), "基因数:", nrow(seu_obj), "\n")

  seu <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = nfeatures)
  hvg <- VariableFeatures(seu)
  cat("HVG 数:", length(hvg), "\n")

  results <- list()

  for (gene in genes_to_analyze) {
    cat("\n", paste(rep("=", 60), collapse = ""), "\n")
    cat(sprintf("==========  KO: %s  ==========\n", gene))
    cat(paste(rep("=", 60), collapse = ""), "\n")

    if (!gene %in% rownames(seu)) {
      cat("跳过: 基因不存在\n")
      results[[gene]] <- "Gene not found"
      next
    }

    expr_pct <- 100 * sum(
      LayerData(seu, layer = "counts")[gene, ] > 0
    ) / ncol(seu)
    cat(sprintf("  %s 全数据集表达率: %.2f%%\n", gene, expr_pct))

    ct_expr <- sapply(levels(factor(seu_obj$cell_type)), function(ct) {
      cells_ct <- colnames(seu_obj)[seu_obj$cell_type == ct]
      if (length(cells_ct) == 0) return(0)
      100 * sum(LayerData(seu, layer = "counts")[gene, cells_ct] > 0) /
        length(cells_ct)
    })
    cat("  各细胞类型表达率 (%):\n")
    for (ct in names(ct_expr)) {
      cat(sprintf("    %-15s %5.1f%%\n", ct, ct_expr[ct]))
    }

    genes_use <- unique(c(gene, hvg))
    genes_use <- genes_use[genes_use %in% rownames(seu)]

    countMatrix <- as.matrix(LayerData(seu, layer = "counts")[genes_use, ])
    cat("矩阵维度:", dim(countMatrix), "\n")

    if (which(rownames(countMatrix) == gene) != 1) {
      new_order <- c(gene, setdiff(rownames(countMatrix), gene))
      countMatrix <- countMatrix[new_order, ]
    }
    cat("目标基因位置:", which(rownames(countMatrix) == gene), "\n")

    t_start <- Sys.time()
    cat("开始时间:", format(t_start), "\n")

    tryCatch({
      res <- scTenifoldKnk(
        countMatrix = countMatrix,
        gKO = gene,
        qc = qc
      )

      t_end <- Sys.time()
      elapsed <- as.numeric(t_end - t_start, units = "hours")
      cat(sprintf("结束时间: %s\n", format(t_end)))
      cat(sprintf("耗时: %.2f 小时\n", elapsed))

      fname <- file.path(save_dir,
                          sprintf("%s_%dHVG_%s.rds", prefix, nfeatures, gene))
      saveRDS(res, fname)
      cat("OK 已保存:", fname, "\n")
      results[[gene]] <- res

    }, error = function(e) {
      cat("FAIL:", e$message, "\n")
      results[[gene]] <- e$message
    })
  }

  return(results)
}

process_ko_result <- function(KO_gene, result, RES_DIR, FIG_DIR) {

  cat("\n", paste(rep("-", 60), collapse = ""), "\n")
  cat(sprintf("处理 %s 结果\n", KO_gene))
  cat(paste(rep("-", 60), collapse = ""), "\n")

  if (is.character(result)) {
    cat(sprintf("⚠ %s 失败, 跳过: %s\n", KO_gene, result))
    return(NULL)
  }

  diffReg <- result$diffRegulation

  cat("diffRegulation 维度:", dim(diffReg), "\n")
  cat("前 10 行:\n"); print(head(diffReg, 10))

  write.csv(diffReg,
            file.path(RES_DIR, sprintf("%s_diffRegulation.csv", KO_gene)),
            row.names = FALSE)

  sig_p <- diffReg[diffReg$p.value < 0.05, ]
  sig_adj <- diffReg[diffReg$p.adj < 0.05, ]
  cat(sprintf("\n显著基因 (%s):\n", KO_gene))
  cat(sprintf("  P < 0.05:    %d 个\n", nrow(sig_p)))
  cat(sprintf("  adj.P<0.05:  %d 个\n", nrow(sig_adj)))

  write.csv(sig_p,
            file.path(RES_DIR, sprintf("%s_significant_P005.csv", KO_gene)),
            row.names = FALSE)
  write.csv(sig_adj,
            file.path(RES_DIR, sprintf("%s_significant_adjP005.csv", KO_gene)),
            row.names = FALSE)

  top20_FC <- head(diffReg[order(-diffReg$FC), ], 21)
  top20_FC <- top20_FC[top20_FC$gene != KO_gene, ]
  top20_FC <- head(top20_FC, 20)
  write.csv(top20_FC,
            file.path(RES_DIR, sprintf("%s_top20_byFC.csv", KO_gene)),
            row.names = FALSE)

  top20_Z <- head(diffReg[order(-abs(diffReg$Z)), ], 21)
  top20_Z <- top20_Z[top20_Z$gene != KO_gene, ]
  top20_Z <- head(top20_Z, 20)
  write.csv(top20_Z,
            file.path(RES_DIR, sprintf("%s_top20_byZ.csv", KO_gene)),
            row.names = FALSE)

  cat(sprintf("\n%s Top 20 by FC:\n", KO_gene)); print(top20_FC)

  cat(sprintf("\n画 Fig.16a (%s)...\n", KO_gene))

  if (nrow(top20_FC) > 0) {
    p_bar <- ggplot(top20_FC, aes(x = reorder(gene, FC), y = FC)) +
      geom_bar(stat = "identity", fill = "#8e1a7f") +
      coord_flip() +
      labs(
        title    = sprintf("Top 20 differentially regulated genes after %s knockout",
                            KO_gene),
        subtitle = "Ranked by FC (perturbation magnitude)",
        x        = "Gene",
        y        = "FC"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
        axis.text.y   = element_text(size = 10)
      )

    ggsave(file.path(FIG_DIR, sprintf("Fig16a_%s_KO_Top20_FC.pdf", KO_gene)),
           p_bar, width = 8, height = 7)
    ggsave(file.path(FIG_DIR, sprintf("Fig16a_%s_KO_Top20_FC.png", KO_gene)),
           p_bar, width = 8, height = 7, dpi = 300, bg = "white")
    cat(sprintf("OK Fig.16a (%s) 保存\n", KO_gene))
  } else {
    cat(sprintf("⚠ Top 20 表为空, 跳过 Fig.16a (%s)\n", KO_gene))
  }

  cat(sprintf("\n画 Fig.16b (%s)...\n", KO_gene))

  df <- diffReg
  df$log_pval <- -log10(df$p.value)

  label_candidates <- df[abs(df$Z) > 2 & df$p.value < 0.05, ]
  label_candidates <- label_candidates[order(-abs(label_candidates$Z)), ]
  label_genes <- head(label_candidates, 30)

  df$category <- "Not significant"
  df$category[df$p.value < 0.05 & abs(df$Z) > 2] <- "Significant (P<0.05, |Z|>2)"
  df$category[df$p.adj < 0.05] <- "Highly significant (adj.P<0.05)"
  df$category <- factor(df$category,
                         levels = c("Not significant",
                                    "Significant (P<0.05, |Z|>2)",
                                    "Highly significant (adj.P<0.05)"))

  p_volcano <- ggplot(df, aes(x = Z, y = log_pval, color = category)) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_color_manual(values = c(
      "Not significant"                 = "gray70",
      "Significant (P<0.05, |Z|>2)"     = "#3182bd",
      "Highly significant (adj.P<0.05)" = "#e41a1c"
    )) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "gray50") +
    geom_text_repel(
      data = label_genes,
      aes(label = gene),
      size = 3,
      color = "black",
      max.overlaps = 30,
      box.padding = 0.4,
      segment.color = "gray40"
    ) +
    labs(
      title    = sprintf("Perturbed genes after %s knockout", KO_gene),
      subtitle = sprintf("Total: %d | P<0.05: %d | adj.P<0.05: %d",
                         nrow(df), nrow(sig_p), nrow(sig_adj)),
      x     = "Z-score",
      y     = "-log10(p-value)",
      color = ""
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.subtitle   = element_text(hjust = 0.5, size = 10, color = "gray40"),
      legend.position = "bottom"
    )

  ggsave(file.path(FIG_DIR, sprintf("Fig16b_%s_KO_Volcano.pdf", KO_gene)),
         p_volcano, width = 9, height = 7)
  ggsave(file.path(FIG_DIR, sprintf("Fig16b_%s_KO_Volcano.png", KO_gene)),
         p_volcano, width = 9, height = 7, dpi = 300, bg = "white")
  cat(sprintf("OK Fig.16b (%s) 保存\n", KO_gene))

  return(data.frame(
    gene         = KO_gene,
    n_total      = nrow(diffReg),
    n_p005       = nrow(sig_p),
    n_adjp005    = nrow(sig_adj),
    top1_byFC    = ifelse(nrow(top20_FC) > 0, top20_FC$gene[1], NA),
    top1_byZ     = ifelse(nrow(top20_Z) > 0, top20_Z$gene[1], NA),
    stringsAsFactors = FALSE
  ))
}

cat("\n========== 开始敲除分析 (3 个生物标志物) ==========\n")

genes_to_ko <- c("Mfge8", "Abca1", "C5ar1")

for (g in genes_to_ko) {
  if (!g %in% rownames(seu_obj)) {
    cat(sprintf("⚠ 基因 %s 不在数据中, 检查命名\n", g))
  }
}

res_list <- run_scTenifoldKnk(
  seu_obj          = seu_obj,
  genes_to_analyze = genes_to_ko,
  nfeatures        = 5000,
  qc               = FALSE,
  prefix           = "scTenifoldKnk",
  save_dir         = RES_DIR
)

rm(seu_obj); gc()

cat("\n========== 整理所有 KO 结果 + 出图 ==========\n")

summary_list <- list()
for (KO_gene in genes_to_ko) {
  result <- res_list[[KO_gene]]
  summary_list[[KO_gene]] <- process_ko_result(KO_gene, result, RES_DIR, FIG_DIR)
}

summary_df <- do.call(rbind, summary_list)
write.csv(summary_df,
          file.path(RES_DIR, "all_KO_summary.csv"),
          row.names = FALSE)

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("    第十九部分 完成\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("敲除基因:", paste(genes_to_ko, collapse = ", "), "\n")
cat("HVG 数量: 5000 (+ 1 目标基因)\n")
cat("\n汇总表:\n")
print(summary_df)

cat("\n输出文件 (每个基因一套):\n")
for (KO_gene in genes_to_ko) {
  cat(sprintf("  Results/scTenifoldKnk_5000HVG_%s.rds\n", KO_gene))
  cat(sprintf("  Results/%s_diffRegulation.csv\n", KO_gene))
  cat(sprintf("  Results/%s_significant_P005.csv\n", KO_gene))
  cat(sprintf("  Results/%s_significant_adjP005.csv\n", KO_gene))
  cat(sprintf("  Results/%s_top20_byFC.csv\n", KO_gene))
  cat(sprintf("  Results/%s_top20_byZ.csv\n", KO_gene))
  cat(sprintf("  Figures/Fig16a_%s_KO_Top20_FC.pdf/png\n", KO_gene))
  cat(sprintf("  Figures/Fig16b_%s_KO_Volcano.pdf/png\n", KO_gene))
}
cat("  Results/all_KO_summary.csv (全局汇总表)\n")

cat("\nOK 第十九部分 完成\n")

sink(file.path(RES_DIR, "session_info.txt"))
sessionInfo()
sink()
