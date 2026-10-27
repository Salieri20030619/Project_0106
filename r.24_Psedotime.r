
cat("\n========================================\n")
cat("24. 拟时序分析 (Slingshot + tradeSeq)\n")
cat("========================================\n\n")

rm(list = ls()); gc()
options(bitmapType = "cairo")

setwd("")

output_dir <- "./24_Pseudotime"
res.path   <- file.path(output_dir, "Results")
fig.path   <- file.path(output_dir, "Figures")

if (dir.exists(output_dir)) {
  cat("清理旧结果目录...\n")
  unlink(output_dir, recursive = TRUE)
}

for (p in c(output_dir, res.path, fig.path)) {
  dir.create(p, recursive = TRUE)
}

BIOMARKERS    <- c("C5ar1", "Abca1", "Mfge8")
BBB_CELLS     <- c("Pericyte", "Endo")
N_HVG         <- 500
N_KNOTS       <- 6
TOP_GENES_N   <- 30

cat("\n========== 加载 R 包 ==========\n")

suppressPackageStartupMessages({
  library(Seurat)
  library(slingshot)
  library(tradeSeq)
  library(SingleCellExperiment)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(pheatmap)
  library(viridis)
  library(mgcv)
  library(BiocParallel)
  library(patchwork)
})
cat("✓ R 包加载完成\n")
cat("  slingshot:", as.character(packageVersion("slingshot")), "\n")
cat("  tradeSeq:",  as.character(packageVersion("tradeSeq")), "\n")

register(MulticoreParam(workers = 4))

t_total <- Sys.time()

cat("\n========== 1. 加载 Seurat 并抽取 BBB 簇 ==========\n")

INPUT_RDS <- "./17_annotation/seurat_obj_annotated_clean.rds"
seu_full <- readRDS(INPUT_RDS)
cat("原始数据:", ncol(seu_full), "细胞 x", nrow(seu_full), "基因\n")

bm_found <- BIOMARKERS[BIOMARKERS %in% rownames(seu_full)]
cat("Biomarker 匹配:", paste(bm_found, collapse = ", "), "\n")
if (length(bm_found) == 0) {
  BIOMARKERS <- toupper(BIOMARKERS)
  bm_found <- BIOMARKERS[BIOMARKERS %in% rownames(seu_full)]
  cat("  重试 (全大写) 匹配:", paste(bm_found, collapse = ", "), "\n")
}
BIOMARKERS <- bm_found

bbb_idx <- seu_full$cell_type %in% BBB_CELLS
seu <- seu_full[, bbb_idx]
cat("\nBBB 簇 (Pericyte + Endo):\n")
print(table(seu$cell_type, seu$condition))

rm(seu_full); gc()

analyze_one_group <- function(seu_grp, group_name) {
  cat("\n", paste(rep("#", 75), collapse = ""), "\n")
  cat("### 开始", group_name, "组分析\n")
  cat(paste(rep("#", 75), collapse = ""), "\n")

  t_grp <- Sys.time()
  n_cells <- ncol(seu_grp)
  cat("细胞数:", n_cells, "\n")
  cat("细胞类型分布:\n")
  print(table(seu_grp$cell_type))

  cat("\n[", group_name, "] 步骤 1/5: PCA + UMAP...\n")
  seu_grp <- NormalizeData(seu_grp, verbose = FALSE)
  seu_grp <- FindVariableFeatures(seu_grp, nfeatures = 2000, verbose = FALSE)
  seu_grp <- ScaleData(seu_grp, verbose = FALSE)
  seu_grp <- RunPCA(seu_grp, npcs = 30, verbose = FALSE)
  n_neighbors <- min(30, n_cells - 1)
  seu_grp <- RunUMAP(seu_grp, dims = 1:20,
                      n.neighbors = n_neighbors, verbose = FALSE)
  cat("  ✓ 降维完成\n")

  cat("\n[", group_name, "] 步骤 2/5: Slingshot 推断拟时序...\n")
  sce <- as.SingleCellExperiment(seu_grp)
  colData(sce)$cell_type <- factor(seu_grp$cell_type)
  umap_coord <- Embeddings(seu_grp, reduction = "umap")
  reducedDim(sce, "UMAP") <- umap_coord

  sce <- slingshot(
    sce,
    clusterLabels = "cell_type",
    reducedDim    = "UMAP",
    start.clus    = NULL
  )

  pseudotime_df <- slingPseudotime(sce)
  n_lineages <- ncol(pseudotime_df)
  cat("  ✓ 推断出", n_lineages, "条 lineage\n")
  cat("  pseudotime 范围:\n")
  print(apply(pseudotime_df, 2, function(x) range(x, na.rm = TRUE)))

  saveRDS(sce, file.path(res.path, paste0("sce_", group_name, ".rds")))

  cat("\n[", group_name, "] 步骤 3/5: 画 Fig 21a 轨迹图...\n")

  ct_colors <- c("Endo" = "#4DBBD5", "Pericyte" = "#E64B35")

  curves_list <- slingCurves(sce)
  curves_df <- do.call(rbind, lapply(seq_along(curves_list), function(i) {
    crv <- curves_list[[i]]
    pts <- crv$s[crv$ord, , drop = FALSE]
    data.frame(
      UMAP_1  = pts[, 1],
      UMAP_2  = pts[, 2],
      Lineage = paste0("Lineage", i)
    )
  }))

  plot_df <- data.frame(
    UMAP_1   = umap_coord[, 1],
    UMAP_2   = umap_coord[, 2],
    CellType = as.character(seu_grp$cell_type),
    Pseudotime = pseudotime_df[, 1]
  )

  p_ct <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = CellType), size = 0.7, alpha = 0.7) +
    scale_color_manual(values = ct_colors) +
    geom_path(
      data = curves_df,
      aes(x = UMAP_1, y = UMAP_2, group = Lineage),
      color = "black", linewidth = 1.2, inherit.aes = FALSE
    ) +
    labs(title = paste0("Fig. 21a (", group_name, "): Trajectory by cell type"),
         x = "UMAP 1", y = "UMAP 2", color = "Cell type") +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "right")

  p_pt <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = Pseudotime), size = 0.7, alpha = 0.8) +
    scale_color_viridis_c(option = "viridis", na.value = "grey80") +
    geom_path(
      data = curves_df,
      aes(x = UMAP_1, y = UMAP_2, group = Lineage),
      color = "black", linewidth = 1.2, inherit.aes = FALSE
    ) +
    labs(title = paste0("Fig. 21a (", group_name, "): Pseudotime"),
         x = "UMAP 1", y = "UMAP 2", color = "Pseudotime") +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "right")

  p_21a <- p_ct | p_pt

  ggsave(file.path(fig.path, paste0("Fig21a_Trajectory_", group_name, ".pdf")),
         p_21a, width = 14, height = 6)
  ggsave(file.path(fig.path, paste0("Fig21a_Trajectory_", group_name, ".png")),
         p_21a, width = 14, height = 6, dpi = 300, bg = "white")
  cat("  ✓ Fig 21a 保存\n")

  cat("\n[", group_name, "] 步骤 4/5: tradeSeq fitGAM ...\n")
  cat("  (基因数:", N_HVG, ", 细胞数:", n_cells, ", 预计 15-30 分钟)\n")
  t_gam <- Sys.time()

  hvg <- VariableFeatures(seu_grp)[1:min(N_HVG, length(VariableFeatures(seu_grp)))]
  hvg <- unique(c(hvg, BIOMARKERS[BIOMARKERS %in% rownames(seu_grp)]))
  cat("  分析基因数:", length(hvg), "\n")

  counts_mat <- GetAssayData(seu_grp, assay = "RNA", layer = "counts")[hvg, ]

  sce_gam <- fitGAM(
    counts     = as.matrix(counts_mat),
    sds        = SlingshotDataSet(sce),
    nknots     = N_KNOTS,
    verbose    = TRUE,
    parallel   = TRUE,
    BPPARAM    = MulticoreParam(workers = 4)
  )

  elapsed_gam <- as.numeric(difftime(Sys.time(), t_gam, units = "mins"))
  cat(sprintf("  ✓ fitGAM 完成, 耗时 %.1f 分钟\n", elapsed_gam))

  assoc_res <- associationTest(sce_gam)
  assoc_res <- assoc_res[order(assoc_res$waldStat, decreasing = TRUE), ]
  assoc_res$gene <- rownames(assoc_res)
  assoc_res$group <- group_name
  write.csv(assoc_res,
            file.path(res.path, paste0("tradeSeq_assoc_", group_name, ".csv")),
            row.names = FALSE)

  cat("  ✓ 显著基因 (pvalue < 0.05):", sum(assoc_res$pvalue < 0.05, na.rm = TRUE), "\n")
  cat("  Top 10 基因:\n")
  print(head(assoc_res[, c("gene", "waldStat", "pvalue")], 10))

  saveRDS(sce_gam, file.path(res.path, paste0("sce_gam_", group_name, ".rds")))

  cat("\n[", group_name, "] 步骤 5/5: 画 Fig 21b 基因表达曲线...\n")

  valid_bm <- BIOMARKERS[BIOMARKERS %in% rownames(sce_gam)]
  if (length(valid_bm) > 0) {
    plots_bm <- list()
    for (bm in valid_bm) {
      p <- plotSmoothers(sce_gam, counts = counts_mat,
                         gene = bm, alpha = 0.5) +
        labs(title = bm) +
        theme_bw(base_size = 11) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))
      plots_bm[[bm]] <- p
    }

    p_21b_bm <- wrap_plots(plots_bm, ncol = length(valid_bm)) +
      plot_annotation(
        title = paste0("Fig. 21b (", group_name, "): Biomarkers along pseudotime"),
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
      )

    ggsave(file.path(fig.path, paste0("Fig21b_Biomarkers_", group_name, ".pdf")),
           p_21b_bm, width = 4 * length(valid_bm), height = 4)
    ggsave(file.path(fig.path, paste0("Fig21b_Biomarkers_", group_name, ".png")),
           p_21b_bm, width = 4 * length(valid_bm), height = 4,
           dpi = 300, bg = "white")
    cat("  ✓ Fig 21b (biomarkers) 保存\n")
  }

  top_genes <- head(assoc_res$gene[!is.na(assoc_res$pvalue)], 6)
  plots_top <- list()
  for (gn in top_genes) {
    if (gn %in% rownames(sce_gam)) {
      p <- plotSmoothers(sce_gam, counts = counts_mat,
                         gene = gn, alpha = 0.5) +
        labs(title = gn) +
        theme_bw(base_size = 11) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
      plots_top[[gn]] <- p
    }
  }

  if (length(plots_top) > 0) {
    p_21b_top <- wrap_plots(plots_top, ncol = 3) +
      plot_annotation(
        title = paste0("Fig. 21b (", group_name, "): Top ",
                       length(plots_top), " dynamic genes"),
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
      )

    ggsave(file.path(fig.path, paste0("Fig21b_TopGenes_", group_name, ".pdf")),
           p_21b_top, width = 12, height = 7)
    ggsave(file.path(fig.path, paste0("Fig21b_TopGenes_", group_name, ".png")),
           p_21b_top, width = 12, height = 7, dpi = 300, bg = "white")
    cat("  ✓ Fig 21b (Top genes) 保存\n")
  }

  cat("\n[", group_name, "] 附加: 画热图...\n")
  top30_genes <- head(assoc_res$gene[!is.na(assoc_res$pvalue) & assoc_res$pvalue < 0.05],
                       TOP_GENES_N)

  if (length(top30_genes) >= 5) {
    pt_vec <- pseudotime_df[, 1]
    cells_ordered <- order(pt_vec, na.last = NA)

    expr_mat <- GetAssayData(seu_grp, assay = "RNA", layer = "data")[top30_genes, cells_ordered]
    expr_mat_scaled <- t(scale(t(as.matrix(expr_mat))))
    expr_mat_scaled[is.na(expr_mat_scaled)] <- 0
    expr_mat_scaled[expr_mat_scaled >  2] <-  2
    expr_mat_scaled[expr_mat_scaled < -2] <- -2

    anno_col <- data.frame(
      CellType   = as.character(seu_grp$cell_type)[cells_ordered],
      Pseudotime = pt_vec[cells_ordered]
    )
    rownames(anno_col) <- colnames(expr_mat_scaled)

    png(file.path(fig.path, paste0("Fig21c_Heatmap_", group_name, ".png")),
        width = 12, height = 8, units = "in", res = 300, bg = "white")
    pheatmap(
      expr_mat_scaled,
      cluster_rows      = TRUE,
      cluster_cols      = FALSE,
      show_colnames     = FALSE,
      annotation_col    = anno_col,
      annotation_colors = list(CellType = ct_colors),
      color             = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
      fontsize_row      = 9,
      main              = paste0("Fig. 21c (", group_name,
                                  "): Top ", TOP_GENES_N, " dynamic genes")
    )
    dev.off()
    cat("  ✓ Fig 21c 热图保存\n")
  }

  elapsed_grp <- as.numeric(difftime(Sys.time(), t_grp, units = "mins"))
  cat(sprintf("\n### %s 组完成, 耗时 %.1f 分钟\n", group_name, elapsed_grp))

  return(list(
    group      = group_name,
    n_cells    = n_cells,
    n_lineages = n_lineages,
    assoc_res  = assoc_res,
    elapsed    = elapsed_grp
  ))
}

results_list <- list()

for (grp in c("Sham", "Stroke")) {
  cells_grp <- colnames(seu)[seu$condition == grp]
  seu_grp <- seu[, cells_grp]
  results_list[[grp]] <- analyze_one_group(seu_grp, grp)
  gc()
}

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("        Sham vs Stroke 对比\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

sham_top   <- head(results_list$Sham$assoc_res, 50)
stroke_top <- head(results_list$Stroke$assoc_res, 50)

cat("Sham Top 50 和 Stroke Top 50 的重叠基因:\n")
overlap <- intersect(sham_top$gene, stroke_top$gene)
cat("重叠数:", length(overlap), "/50\n")
cat("重叠基因:", paste(head(overlap, 20), collapse = ", "),
    ifelse(length(overlap) > 20, "...", ""), "\n")

cat("\nSham 特异 Top 20:", paste(head(setdiff(sham_top$gene, stroke_top$gene), 20),
                                     collapse = ", "), "\n")
cat("\nStroke 特异 Top 20:", paste(head(setdiff(stroke_top$gene, sham_top$gene), 20),
                                       collapse = ", "), "\n")

compare_df <- data.frame(
  metric     = c("n_cells", "n_lineages", "n_sig_pvalue0.05"),
  Sham       = c(results_list$Sham$n_cells,
                  results_list$Sham$n_lineages,
                  sum(results_list$Sham$assoc_res$pvalue < 0.05, na.rm = TRUE)),
  Stroke     = c(results_list$Stroke$n_cells,
                  results_list$Stroke$n_lineages,
                  sum(results_list$Stroke$assoc_res$pvalue < 0.05, na.rm = TRUE))
)
write.csv(compare_df,
          file.path(res.path, "sham_vs_stroke_summary.csv"),
          row.names = FALSE)
print(compare_df)

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("       第 24 部分 Slingshot + tradeSeq 完成\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

total_min <- as.numeric(difftime(Sys.time(), t_total, units = "mins"))
cat(sprintf("\n总耗时: %.1f 分钟\n", total_min))

cat("\n输出文件:\n")
cat("  Results/\n")
cat("    sce_{Sham,Stroke}.rds              Slingshot 对象\n")
cat("    sce_gam_{Sham,Stroke}.rds          tradeSeq GAM 对象\n")
cat("    tradeSeq_assoc_{Sham,Stroke}.csv   沿 pseudotime 差异基因\n")
cat("    sham_vs_stroke_summary.csv         两组对比\n")
cat("  Figures/\n")
cat("    Fig21a_Trajectory_{Sham,Stroke}.pdf/png\n")
cat("    Fig21b_Biomarkers_{Sham,Stroke}.pdf/png\n")
cat("    Fig21b_TopGenes_{Sham,Stroke}.pdf/png\n")
cat("    Fig21c_Heatmap_{Sham,Stroke}.png\n")

cat("\n✅ 第 24 部分完成!\n")

sink(file.path(res.path, "session_info.txt"))
sessionInfo()
sink()
