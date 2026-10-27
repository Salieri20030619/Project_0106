
cat("\n========================================\n")
cat("22. ReactomeGSA 富集分析\n")
cat("========================================\n\n")

rm(list = ls()); gc()
options(bitmapType = "cairo")

options(reactome_gsa.url = "http://gsa.reactome.org")
Sys.setenv(CURL_SSL_VERIFYPEER = "0")
Sys.setenv(CURL_SSL_VERIFYHOST = "0")
options(download.file.method = "libcurl")
suppressPackageStartupMessages(library(httr))
httr::set_config(httr::config(ssl_verifypeer = 0L, ssl_verifyhost = 0L))
cat("✓ API URL 设为 HTTP, SSL 验证禁用\n")

cat("\n========== 加载 R 包 ==========\n")
suppressPackageStartupMessages({
  library(Seurat)
  library(ReactomeGSA)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
})
cat("✓ ReactomeGSA 版本:", as.character(packageVersion("ReactomeGSA")), "\n")

setwd("")

output_dir <- "./22_ReactomeGSA"
res.path <- file.path(output_dir, "Results")
fig.path <- file.path(output_dir, "Figures")
for (p in c(output_dir, res.path, fig.path)) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE)
}

cat("\n========== 1. 加载 Seurat 对象 ==========\n")

INPUT_RDS <- "./17_annotation/seurat_obj_annotated_clean.rds"
if (!file.exists(INPUT_RDS)) {
  stop("❌ 输入文件不存在: ", INPUT_RDS)
}

t0 <- Sys.time()
seu <- readRDS(INPUT_RDS)
cat("✓ Seurat 对象加载完成\n")
cat("  细胞数:", ncol(seu), "\n")
cat("  基因数:", nrow(seu), "\n")

if (!"cell_type" %in% colnames(seu@meta.data)) {
  stop("❌ 没有 cell_type 字段")
}

Idents(seu) <- seu$cell_type

cat("  细胞类型分布:\n")
print(table(seu$cell_type))

cat("\n========== 2. 运行 ReactomeGSA ==========\n")
cat("(上传到 gsa.reactome.org API, ssGSEA 方法)\n")
cat("开始时间:", format(Sys.time()), "\n")

t1 <- Sys.time()

gsva_result <- tryCatch({
  analyse_sc_clusters(
    object  = seu,
    verbose = TRUE
  )
}, error = function(e) {
  cat("❌ 失败:", e$message, "\n")
  return(NULL)
})

t2 <- Sys.time()
elapsed <- as.numeric(difftime(t2, t1, units = "mins"))
cat(sprintf("\n耗时: %.1f 分钟\n", elapsed))

if (is.null(gsva_result)) {
  stop("❌ ReactomeGSA 分析失败, 可能是 Reactome 服务器暂时故障, 请稍后重试")
}

cat("✓ 分析完成\n")

saveRDS(gsva_result, file.path(res.path, "reactome_gsa_result.rds"))

cat("\n========== 3. 提取通路活性矩阵 ==========\n")

pathway_expression <- pathways(gsva_result)
cat("  通路数:", nrow(pathway_expression), "\n")
cat("  列数:", ncol(pathway_expression), "\n")

write.csv(pathway_expression,
          file.path(res.path, "reactome_pathway_matrix.csv"),
          row.names = TRUE)
cat("✓ 通路矩阵已保存\n")

cat("\n========== 4. 筛选 Top 10 最大差异通路 ==========\n")

name_col <- which(colnames(pathway_expression) == "Name")
if (length(name_col) > 0) {
  pathway_names <- pathway_expression$Name
  pathway_matrix <- as.matrix(pathway_expression[, -name_col, drop = FALSE])
  rownames(pathway_matrix) <- pathway_names
} else if (is.character(pathway_expression[, 1])) {
  pathway_matrix <- as.matrix(pathway_expression[, -1, drop = FALSE])
  rownames(pathway_matrix) <- pathway_expression[, 1]
} else {
  pathway_matrix <- as.matrix(pathway_expression)
}

storage.mode(pathway_matrix) <- "numeric"
valid_rows <- apply(pathway_matrix, 1, function(x) sum(!is.na(x)) >= 2)
pathway_matrix <- pathway_matrix[valid_rows, ]
cat("  有效通路数:", nrow(pathway_matrix), "\n")

pathway_range <- apply(pathway_matrix, 1, function(x) {
  max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
})
pathway_range[!is.finite(pathway_range)] <- 0

n_top <- min(10, sum(pathway_range > 0))
top10_idx <- order(pathway_range, decreasing = TRUE)[1:n_top]
top10_matrix <- pathway_matrix[top10_idx, , drop = FALSE]

cat(sprintf("\n✓ Top %d 最大差异通路:\n", n_top))
for (i in seq_along(top10_idx)) {
  cat(sprintf("  %2d. [range=%.3f] %s\n",
              i, pathway_range[top10_idx[i]],
              substr(rownames(top10_matrix)[i], 1, 80)))
}

top10_df <- data.frame(
  Pathway = rownames(top10_matrix),
  Range   = pathway_range[top10_idx],
  top10_matrix
)
write.csv(top10_df,
          file.path(res.path, "reactome_top10_pathways.csv"),
          row.names = FALSE)
cat("\n✓ Top 10 表已保存\n")

cat("\n========== 5. 画 Fig 19 热图 ==========\n")

clean_cols <- colnames(top10_matrix)
clean_cols <- gsub("\\.Seurat$", "", clean_cols)
clean_cols <- gsub("\\.[^.]+$", "", clean_cols)
colnames(top10_matrix) <- clean_cols

heatmap_scaled <- t(scale(t(top10_matrix)))
heatmap_scaled[is.na(heatmap_scaled)] <- 0
heatmap_scaled[heatmap_scaled >  3] <-  3
heatmap_scaled[heatmap_scaled < -3] <- -3

shorten_name <- function(x, maxlen = 60) {
  ifelse(nchar(x) > maxlen, paste0(substr(x, 1, maxlen - 3), "..."), x)
}
rownames(heatmap_scaled) <- shorten_name(rownames(heatmap_scaled))

pheatmap(
  heatmap_scaled,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  cluster_rows = TRUE, cluster_cols = TRUE,
  show_rownames = TRUE, show_colnames = TRUE,
  fontsize_row = 10, fontsize_col = 11, angle_col = 45,
  main = "Fig. 19: Top 10 differentially expressed Reactome pathways",
  filename = file.path(fig.path, "Fig19_Reactome_Heatmap.pdf"),
  width = 12, height = 7
)

png(file.path(fig.path, "Fig19_Reactome_Heatmap.png"),
    width = 12, height = 7, units = "in", res = 300, bg = "white")
pheatmap(
  heatmap_scaled,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  cluster_rows = TRUE, cluster_cols = TRUE,
  show_rownames = TRUE, show_colnames = TRUE,
  fontsize_row = 10, fontsize_col = 11, angle_col = 45,
  main = "Fig. 19: Top 10 differentially expressed Reactome pathways"
)
dev.off()

cat("✓ Fig 19 热图保存完成\n")

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("       第 22 部分 ReactomeGSA 完成\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

total_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("\n总耗时: %.1f 分钟\n", total_min))

cat("\n输出文件:\n")
cat("  Results/\n")
cat("    reactome_gsa_result.rds           完整 ReactomeGSA 对象\n")
cat("    reactome_pathway_matrix.csv       全部通路 × 11 细胞类型\n")
cat("    reactome_top10_pathways.csv       Top 10 差异通路\n")
cat("  Figures/\n")
cat("    Fig19_Reactome_Heatmap.pdf\n")
cat("    Fig19_Reactome_Heatmap.png\n")

cat("\n✅ 第 22 部分完成!\n")

sink(file.path(res.path, "session_info.txt"))
sessionInfo()
sink()
