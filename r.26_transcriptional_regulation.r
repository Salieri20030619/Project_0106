
cat("\n========================================\n")
cat("26. SCENIC 转录调控网络分析\n")
cat("========================================\n\n")

rm(list = ls()); gc()
options(bitmapType = "cairo")

BASE_DIR   <- ""
OUTPUT_DIR <- file.path(BASE_DIR, "26_transcriptional_regulation")
DB_DIR     <- "SCENIC"

setwd(BASE_DIR)

N_CELLS_PER_TYPE <- 500
N_TREES          <- 1000
N_TOP_TF         <- 30

res.path <- file.path(OUTPUT_DIR, "Results")
fig.path <- file.path(OUTPUT_DIR, "Figures")
int.path <- file.path(OUTPUT_DIR, "int")
out.path <- file.path(OUTPUT_DIR, "output")

cat("========== 清理旧结果 (注释体系已切换, 必须重跑全流程) ==========\n")

old_files <- c(
  file.path(OUTPUT_DIR, "02.scenic_regulons.csv"),
  file.path(OUTPUT_DIR, "03.Fig23_tf_heatmap.pdf"),
  file.path(OUTPUT_DIR, "03.Fig23_tf_heatmap.png"),
  file.path(OUTPUT_DIR, "04.Fig23_tf_dotplot.pdf"),
  file.path(OUTPUT_DIR, "04.Fig23_tf_dotplot.png"),
  file.path(OUTPUT_DIR, "05.tf_umap.png"),
  file.path(OUTPUT_DIR, "06.auc_scores.csv"),
  file.path(OUTPUT_DIR, "09.analysis_summary.txt"),
  file.path(OUTPUT_DIR, "10.output_files_list.csv"),
  file.path(OUTPUT_DIR, "analysis_summary.txt"),
  file.path(OUTPUT_DIR, "exprMatrix_filtered.rds"),
  file.path(OUTPUT_DIR, "output_files_list.csv"),
  file.path(OUTPUT_DIR, "scenic_analysis.log"),
  file.path(OUTPUT_DIR, "scenicSettings.rds")
)
for (f in old_files) {
  if (file.exists(f)) file.remove(f)
}

for (sub in c("SCENIC", "int", "output", "Data", "Figures", "Results", "SCENIC_db")) {
  sd <- file.path(OUTPUT_DIR, sub)
  if (dir.exists(sd)) {
    unlink(sd, recursive = TRUE)
    cat(sprintf("  ✓ 已清空 %s/\n", sub))
  }
}
cat("✓ 旧结果已彻底清空\n")

for (p in c(res.path, fig.path, int.path, out.path)) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

setwd(OUTPUT_DIR)

cat("\n========== 加载 R 包 ==========\n")
suppressPackageStartupMessages({
  library(Seurat)
  library(SCENIC)
  library(AUCell)
  library(GENIE3)
  library(RcisTarget)
  library(arrow)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(tidyr)
  library(dplyr)
  library(doParallel)
  library(foreach)
  library(SingleCellExperiment)
})

cat("✓ SCENIC:", as.character(packageVersion("SCENIC")), "\n")

nCores <- max(1, parallel::detectCores() - 4)
if (nCores > 16) nCores <- 16
cat("使用核心数:", nCores, "\n")
registerDoParallel(cores = nCores)

set.seed(12345)
t_total <- Sys.time()

cat("\n========== 加载 motifAnnotations ==========\n")

motif_tbl <- file.path(DB_DIR, "motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl")
if (!file.exists(motif_tbl)) stop("❌ motif 注释缺失: ", motif_tbl)

if (exists("importAnnotations", where = asNamespace("RcisTarget"))) {
  motifAnnotations_mgi <- RcisTarget::importAnnotations(motif_tbl)
  cat("✓ 用 importAnnotations 加载\n")
} else {
  motifAnnotations_mgi <- data.table::fread(motif_tbl)
  cat("✓ 用 fread 加载\n")
}

cat("  维度:", dim(motifAnnotations_mgi), "\n")
cat("  列名:", paste(colnames(motifAnnotations_mgi), collapse = ", "), "\n")

assign("motifAnnotations_mgi", motifAnnotations_mgi, envir = .GlobalEnv)

cat("\n========== 1. 加载 Seurat ==========\n")

INPUT_RDS <- file.path(BASE_DIR, "17_annotation/seurat_obj_annotated_clean.rds")
seu <- readRDS(INPUT_RDS)
cat("✓ 加载完成:", ncol(seu), "细胞 x", nrow(seu), "基因\n")

if (!"cell_type" %in% colnames(seu@meta.data)) stop("❌ 缺 cell_type")
if (!"condition" %in% colnames(seu@meta.data)) stop("❌ 缺 condition")

cat("\n当前细胞类型分布:\n")
print(table(seu$cell_type))

cat("\n========== 2. 降采样 (每类", N_CELLS_PER_TYPE, ") ==========\n")

set.seed(12345)
cells_keep <- c()
for (ct in unique(seu$cell_type)) {
  cells_ct <- colnames(seu)[seu$cell_type == ct]
  if (length(cells_ct) <= N_CELLS_PER_TYPE) {
    cells_keep <- c(cells_keep, cells_ct)
  } else {
    cells_keep <- c(cells_keep, sample(cells_ct, N_CELLS_PER_TYPE))
  }
}

seu_sub <- seu[, cells_keep]
cat("✓ 降采样完成:", ncol(seu_sub), "细胞\n")
cat("降采样后各细胞类型分布:\n")
print(table(seu_sub$cell_type))

rm(seu); gc()

cat("\n========== 3. 准备表达矩阵 ==========\n")

exprMat <- GetAssayData(seu_sub, assay = "RNA", layer = "counts")
exprMat <- as.matrix(exprMat)
cat("表达矩阵维度:", dim(exprMat), "\n")

cellInfo <- data.frame(
  CellType  = as.character(seu_sub$cell_type),
  Condition = as.character(seu_sub$condition),
  row.names = colnames(seu_sub)
)
saveRDS(cellInfo, file.path(int.path, "cellInfo.Rds"))

cat("\n========== 4. 初始化 SCENIC ==========\n")

dbs_files <- c(
  "10kb" = "mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
)

for (nm in names(dbs_files)) {
  f <- file.path(DB_DIR, dbs_files[nm])
  if (!file.exists(f)) stop("❌ 数据库缺失: ", f)
  cat(sprintf("  ✓ %s: %s (%.1f GB)\n", nm, dbs_files[nm],
              file.info(f)$size / 1024^3))
}

scenicOptions <- initializeScenic(
  org          = "mgi",
  dbDir        = DB_DIR,
  dbs          = dbs_files,
  datasetTitle = "Ischemic_Stroke_scRNA_v10_Yu2024Panel",
  nCores       = nCores
)

saveRDS(scenicOptions, file.path(int.path, "scenicOptions.Rds"))
cat("✓ SCENIC 配置完成\n")

cat("\n========== 5. 基因过滤 ==========\n")

genesKept <- geneFiltering(
  exprMat,
  scenicOptions    = scenicOptions,
  minCountsPerGene = 3 * 0.01 * ncol(exprMat),
  minSamples       = ncol(exprMat) * 0.01
)
cat("✓ 保留基因数:", length(genesKept), "\n")

exprMat_filtered <- exprMat[genesKept, ]
rm(exprMat); gc()

cat("\n========== 6. Spearman 相关性 ==========\n")
t1 <- Sys.time()

runCorrelation(exprMat_filtered, scenicOptions)

cat(sprintf("✓ 耗时 %.1f 分钟\n",
            as.numeric(difftime(Sys.time(), t1, units = "mins"))))

cat("\n========== 7. GENIE3 ==========\n")
cat("参数: RF, K=sqrt, nTrees=", N_TREES, "\n")
t_genie <- Sys.time()

exprMat_log <- log2(exprMat_filtered + 1)

runGenie3(
  exprMat_log,
  scenicOptions,
  nParts     = 10,
  nTrees     = N_TREES,
  treeMethod = "RF",
  K          = "sqrt"
)

elapsed_genie <- as.numeric(difftime(Sys.time(), t_genie, units = "mins"))
cat(sprintf("✓ GENIE3 耗时 %.1f 分钟\n", elapsed_genie))

cat("\n========== 8. regulons + AUCell ==========\n")
t2 <- Sys.time()

scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores  <- nCores
scenicOptions@settings$seed    <- 12345

runSCENIC_1_coexNetwork2modules(scenicOptions)
cat("  ✓ 共表达模块\n")

runSCENIC_2_createRegulons(
  scenicOptions,
  coexMethod = c("top5perTarget", "top10perTarget", "top50perTarget",
                 "top5", "top50")
)
cat("  ✓ RcisTarget\n")

runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
cat("  ✓ AUCell\n")

elapsed_reg <- as.numeric(difftime(Sys.time(), t2, units = "mins"))
cat(sprintf("\n✓ 耗时 %.1f 分钟\n", elapsed_reg))

cat("\n========== 9. 提取 AUC ==========\n")

regulonAUC_file <- getOutName(scenicOptions, "s3_aucellRegulonAUC")
regulonAUC <- readRDS(regulonAUC_file)

auc_mat <- getAUC(regulonAUC)
cat("  AUC 矩阵维度:", dim(auc_mat), "\n")

write.csv(auc_mat, file.path(res.path, "auc_scores.csv"))

cat("\n========== 10. 细胞类型平均 AUC ==========\n")

cellInfo_match <- cellInfo[colnames(auc_mat), , drop = FALSE]

regulons_ct <- sapply(unique(cellInfo_match$CellType), function(ct) {
  cells_ct <- rownames(cellInfo_match)[cellInfo_match$CellType == ct]
  rowMeans(auc_mat[, cells_ct, drop = FALSE])
})
cat("  维度:", dim(regulons_ct), "\n")
cat("  细胞类型 (应为 10 类):", paste(colnames(regulons_ct), collapse = ", "), "\n")

write.csv(regulons_ct,
          file.path(res.path, "tf_activity_by_celltype.csv"),
          row.names = TRUE)

regulons <- loadInt(scenicOptions, "regulons")
regulon_df <- data.frame(
  Regulon   = names(regulons),
  N_Targets = sapply(regulons, length),
  Targets   = sapply(regulons, function(x) paste(head(x, 20), collapse = ";"))
)
write.csv(regulon_df,
          file.path(res.path, "scenic_regulons.csv"),
          row.names = FALSE)
cat("✓ regulons:", nrow(regulon_df), "\n")

cat("\n========== 11. Fig 23a ==========\n")

tf_range <- apply(regulons_ct, 1, function(x) max(x) - min(x))
top_tf_idx <- order(tf_range, decreasing = TRUE)[1:min(N_TOP_TF, length(tf_range))]
top_tf_mat <- regulons_ct[top_tf_idx, , drop = FALSE]

heatmap_scaled <- t(scale(t(top_tf_mat)))
heatmap_scaled[is.na(heatmap_scaled)] <- 0
heatmap_scaled[heatmap_scaled >  2] <-  2
heatmap_scaled[heatmap_scaled < -2] <- -2

pretty_tf <- function(x) gsub("\\s*\\([0-9]+g\\)", "", x)
rownames(heatmap_scaled) <- pretty_tf(rownames(heatmap_scaled))

pheatmap(
  heatmap_scaled,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  cluster_rows = TRUE, cluster_cols = TRUE,
  show_rownames = TRUE, show_colnames = TRUE,
  fontsize_row = 9, fontsize_col = 11, angle_col = 45,
  main = paste0("Fig. 23a: Top ", N_TOP_TF, " TF regulon activity (AUC)"),
  filename = file.path(fig.path, "Fig23a_TF_AUC_Heatmap.pdf"),
  width = 10, height = 9
)

png(file.path(fig.path, "Fig23a_TF_AUC_Heatmap.png"),
    width = 10, height = 9, units = "in", res = 300, bg = "white")
pheatmap(
  heatmap_scaled,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  cluster_rows = TRUE, cluster_cols = TRUE,
  show_rownames = TRUE, show_colnames = TRUE,
  fontsize_row = 9, fontsize_col = 11, angle_col = 45,
  main = paste0("Fig. 23a: Top ", N_TOP_TF, " TF regulon activity (AUC)")
)
dev.off()
cat("✓ Fig 23a\n")

cat("\n========== 12. Fig 23b ==========\n")

top10_tf_names <- pretty_tf(rownames(heatmap_scaled)[1:min(10, nrow(heatmap_scaled))])
top10_tf_clean <- gsub("_extended$", "", top10_tf_names)
tf_in_data <- top10_tf_clean[top10_tf_clean %in% rownames(seu_sub)]
cat("Top TF:", paste(tf_in_data, collapse = ", "), "\n")

if (length(tf_in_data) > 0) {
  Idents(seu_sub) <- seu_sub$cell_type

  p_dot <- DotPlot(seu_sub,
                    features = tf_in_data,
                    cols = c("lightgrey", "#E64B35"),
                    dot.scale = 6) +
    RotatedAxis() +
    labs(title = "Fig. 23b: Top 10 TF gene expression",
         x = "", y = "Cell type") +
    theme(
      plot.title  = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave(file.path(fig.path, "Fig23b_TF_Expression_DotPlot.pdf"),
         p_dot, width = max(8, length(tf_in_data) * 0.8), height = 6)
  ggsave(file.path(fig.path, "Fig23b_TF_Expression_DotPlot.png"),
         p_dot, width = max(8, length(tf_in_data) * 0.8), height = 6,
         dpi = 300, bg = "white")
  cat("✓ Fig 23b\n")
}

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("       第 26 部分完成\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

total_min <- as.numeric(difftime(Sys.time(), t_total, units = "mins"))
cat(sprintf("\n总耗时: %.1f 分钟 (%.1f 小时)\n", total_min, total_min/60))

cat("\n关键数字:\n")
cat("  细胞数:", ncol(seu_sub), "\n")
cat("  细胞类型数:", length(unique(seu_sub$cell_type)), "\n")
cat("  基因数:", length(genesKept), "\n")
cat("  regulons:", nrow(regulon_df), "\n")
cat("  Top 10 TF:\n")
print(head(data.frame(
  Regulon = rownames(heatmap_scaled),
  Range   = round(tf_range[top_tf_idx], 3)
), 10))

cat("\n✅ 完成\n")

sink(file.path(res.path, "session_info.txt"))
sessionInfo()
sink()
