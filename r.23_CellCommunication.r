
cat("\n========================================\n")
cat("23. 细胞通讯分析 (CellChat)\n")
cat("========================================\n\n")

rm(list = ls()); gc()
options(bitmapType = "cairo")

setwd("")

output_dir <- "./23_CellCommunication"
res.path   <- file.path(output_dir, "Results")
fig.path   <- file.path(output_dir, "Figures")
data.path  <- file.path(output_dir, "Data")

if (dir.exists(output_dir)) {
  cat("清理旧结果目录...\n")
  unlink(output_dir, recursive = TRUE)
}

for (p in c(output_dir, res.path, fig.path, data.path)) {
  dir.create(p, recursive = TRUE)
}

cat("\n========== 加载 R 包 ==========\n")

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

cat("✓ CellChat 版本:", as.character(packageVersion("CellChat")), "\n")

future::plan("multisession", workers = 4)
options(future.globals.maxSize = 20 * 1024^3)

t0 <- Sys.time()

cat("\n========== 1. 加载 Seurat 对象 ==========\n")

INPUT_RDS <- "./17_annotation/seurat_obj_annotated_clean.rds"
if (!file.exists(INPUT_RDS)) {
  stop("❌ 输入文件不存在: ", INPUT_RDS)
}

seu <- readRDS(INPUT_RDS)
cat("✓ 加载完成:", ncol(seu), "细胞 x", nrow(seu), "基因\n")

if (!"cell_type" %in% colnames(seu@meta.data)) stop("❌ 缺少 cell_type 字段")
if (!"condition" %in% colnames(seu@meta.data)) stop("❌ 缺少 condition 字段")

cat("细胞类型分布:\n"); print(table(seu$cell_type))
cat("分组分布:\n"); print(table(seu$condition))

cat("\n========== 2. 分组运行 CellChat (Sham / Stroke) ==========\n")

CellChatDB <- CellChatDB.mouse
cat("CellChatDB 配体-受体对数:", nrow(CellChatDB$interaction), "\n")

cellchat_list <- list()

for (grp in c("Sham", "Stroke")) {
  cat("\n========== 处理", grp, "组 ==========\n")
  t_grp <- Sys.time()

  cells_grp <- colnames(seu)[seu$condition == grp]
  cat("  细胞数:", length(cells_grp), "\n")

  data.input <- GetAssayData(seu, assay = "RNA", layer = "data")[, cells_grp, drop = FALSE]
  meta_grp <- data.frame(
    cell_type = factor(as.character(seu$cell_type[cells_grp])),
    row.names = cells_grp
  )

  cellchat <- createCellChat(
    object   = data.input,
    meta     = meta_grp,
    group.by = "cell_type"
  )
  cellchat@DB <- CellChatDB

  cat("  预处理: subsetData ...\n")
  cellchat <- subsetData(cellchat)

  cat("  识别过表达基因 ...\n")
  cellchat <- identifyOverExpressedGenes(cellchat)

  cat("  识别过表达相互作用 ...\n")
  cellchat <- identifyOverExpressedInteractions(cellchat)

  cat("  计算通讯概率 (triMean) ...\n")
  cellchat <- computeCommunProb(cellchat, type = "triMean")

  cat("  过滤通讯 (min.cells=10) ...\n")
  cellchat <- filterCommunication(cellchat, min.cells = 10)

  cat("  通路水平聚合 ...\n")
  cellchat <- computeCommunProbPathway(cellchat)

  cellchat <- aggregateNet(cellchat)

  t_grp_end <- Sys.time()
  elapsed_grp <- as.numeric(difftime(t_grp_end, t_grp, units = "mins"))
  cat(sprintf("\n  ✓ %s 完成, 耗时 %.1f 分钟\n", grp, elapsed_grp))
  cat("  通讯数:", sum(cellchat@net$count), "\n")
  cat("  通路数:", length(cellchat@netP$pathways), "\n")

  cellchat_list[[grp]] <- cellchat
  saveRDS(cellchat, file.path(data.path, paste0("cellchat_", grp, ".rds")))
}

cat("\n========== 3. Fig 20a 通讯网络图 ==========\n")

for (grp in names(cellchat_list)) {
  cellchat <- cellchat_list[[grp]]

  if (is.null(cellchat@net$count) || sum(cellchat@net$count) == 0) {
    cat("  ", grp, "无通讯, 跳过\n")
    next
  }

  vertex_weight <- as.numeric(table(cellchat@meta$cell_type))
  names(vertex_weight) <- names(table(cellchat@meta$cell_type))

  pdf(file.path(fig.path, paste0("Fig20a_Network_Count_", grp, ".pdf")),
      width = 10, height = 10)
  par(mar = c(1, 1, 3, 1))
  netVisual_circle(cellchat@net$count,
                   vertex.weight = vertex_weight,
                   weight.scale  = TRUE,
                   label.edge    = FALSE,
                   title.name    = paste0("Fig. 20a: ", grp, " - Number of interactions"))
  dev.off()

  png(file.path(fig.path, paste0("Fig20a_Network_Count_", grp, ".png")),
      width = 10, height = 10, units = "in", res = 300, bg = "white")
  par(mar = c(1, 1, 3, 1))
  netVisual_circle(cellchat@net$count,
                   vertex.weight = vertex_weight,
                   weight.scale  = TRUE,
                   label.edge    = FALSE,
                   title.name    = paste0("Fig. 20a: ", grp, " - Number of interactions"))
  dev.off()

  pdf(file.path(fig.path, paste0("Fig20a_Network_Weight_", grp, ".pdf")),
      width = 10, height = 10)
  par(mar = c(1, 1, 3, 1))
  netVisual_circle(cellchat@net$weight,
                   vertex.weight = vertex_weight,
                   weight.scale  = TRUE,
                   label.edge    = FALSE,
                   title.name    = paste0("Fig. 20a: ", grp, " - Interaction strength"))
  dev.off()

  png(file.path(fig.path, paste0("Fig20a_Network_Weight_", grp, ".png")),
      width = 10, height = 10, units = "in", res = 300, bg = "white")
  par(mar = c(1, 1, 3, 1))
  netVisual_circle(cellchat@net$weight,
                   vertex.weight = vertex_weight,
                   weight.scale  = TRUE,
                   label.edge    = FALSE,
                   title.name    = paste0("Fig. 20a: ", grp, " - Interaction strength"))
  dev.off()

  cat("  ✓", grp, "网络图保存\n")
}

cat("\n画 Sham vs Stroke 2x2 对比拼图...\n")
tryCatch({
  pdf(file.path(fig.path, "Fig20a_Combined_Networks.pdf"),
      width = 16, height = 16)
  par(mfrow = c(2, 2), mar = c(1, 1, 3, 1))
  for (grp in c("Sham", "Stroke")) {
    cellchat <- cellchat_list[[grp]]
    vertex_weight <- as.numeric(table(cellchat@meta$cell_type))
    names(vertex_weight) <- names(table(cellchat@meta$cell_type))
    netVisual_circle(cellchat@net$count, vertex.weight = vertex_weight,
                     weight.scale = TRUE, label.edge = FALSE,
                     title.name = paste0(grp, " - Number"))
    netVisual_circle(cellchat@net$weight, vertex.weight = vertex_weight,
                     weight.scale = TRUE, label.edge = FALSE,
                     title.name = paste0(grp, " - Strength"))
  }
  dev.off()

  png(file.path(fig.path, "Fig20a_Combined_Networks.png"),
      width = 16, height = 16, units = "in", res = 300, bg = "white")
  par(mfrow = c(2, 2), mar = c(1, 1, 3, 1))
  for (grp in c("Sham", "Stroke")) {
    cellchat <- cellchat_list[[grp]]
    vertex_weight <- as.numeric(table(cellchat@meta$cell_type))
    names(vertex_weight) <- names(table(cellchat@meta$cell_type))
    netVisual_circle(cellchat@net$count, vertex.weight = vertex_weight,
                     weight.scale = TRUE, label.edge = FALSE,
                     title.name = paste0(grp, " - Number"))
    netVisual_circle(cellchat@net$weight, vertex.weight = vertex_weight,
                     weight.scale = TRUE, label.edge = FALSE,
                     title.name = paste0(grp, " - Strength"))
  }
  dev.off()
  cat("  ✓ 对比拼图保存\n")
}, error = function(e) {
  cat("  ⚠ 拼图失败:", e$message, "\n")
})

cat("\n========== 4. Fig 20b 配体-受体气泡图 ==========\n")

for (grp in names(cellchat_list)) {
  cellchat <- cellchat_list[[grp]]

  if (is.null(cellchat@netP$pathways) || length(cellchat@netP$pathways) == 0) {
    cat("  ", grp, "无通路, 跳过\n")
    next
  }

  pathways_use <- cellchat@netP$pathways
  if (length(pathways_use) > 10) pathways_use <- pathways_use[1:10]

  tryCatch({
    p_bubble <- netVisual_bubble(
      cellchat,
      sources.use    = 1:length(unique(cellchat@meta$cell_type)),
      targets.use    = 1:length(unique(cellchat@meta$cell_type)),
      signaling      = pathways_use,
      remove.isolate = TRUE,
      title.name     = paste0("Fig. 20b: ", grp, " - Ligand-Receptor Interactions")
    )

    ggsave(file.path(fig.path, paste0("Fig20b_Bubble_", grp, ".pdf")),
           p_bubble, width = 14, height = 12)
    ggsave(file.path(fig.path, paste0("Fig20b_Bubble_", grp, ".png")),
           p_bubble, width = 14, height = 12, dpi = 300, bg = "white")
    cat("  ✓", grp, "气泡图保存\n")
  }, error = function(e) {
    cat("  ⚠", grp, "气泡图失败:", e$message, "\n")
  })
}

cat("\n========== 5. Sham vs Stroke 对比分析 ==========\n")

tryCatch({

  cellchat_merged <- mergeCellChat(cellchat_list,
                                     add.names = names(cellchat_list))
  saveRDS(cellchat_merged, file.path(data.path, "cellchat_comparison.rds"))

  gg_count <- compareInteractions(cellchat_merged,
                                    show.legend = FALSE, group = c(1, 2))
  gg_weight <- compareInteractions(cellchat_merged,
                                     show.legend = FALSE, group = c(1, 2),
                                     measure = "weight")

  p_compare <- gg_count + gg_weight
  ggsave(file.path(fig.path, "Fig20c_Compare_Total.pdf"),
         p_compare, width = 10, height = 5)
  ggsave(file.path(fig.path, "Fig20c_Compare_Total.png"),
         p_compare, width = 10, height = 5, dpi = 300, bg = "white")
  cat("  ✓ Sham vs Stroke 总数对比图保存\n")
}, error = function(e) {
  cat("  ⚠ 对比分析失败:", e$message, "\n")
})

cat("\n========== 6. 生成汇总表 ==========\n")

summary_df <- data.frame(
  group = names(cellchat_list),
  n_cells = sapply(cellchat_list, function(x) ncol(x@meta)),
  n_interactions = sapply(cellchat_list, function(x)
    ifelse(!is.null(x@net$count), sum(x@net$count), 0)),
  n_pathways = sapply(cellchat_list, function(x)
    ifelse(!is.null(x@netP$pathways), length(x@netP$pathways), 0)),
  stringsAsFactors = FALSE
)
write.csv(summary_df,
          file.path(res.path, "interaction_summary.csv"),
          row.names = FALSE)
print(summary_df)

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("       第 23 部分 CellChat 完成\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

total_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("\n总耗时: %.1f 分钟\n", total_min))

cat("\n输出文件:\n")
cat("  Data/        cellchat_Sham.rds, cellchat_Stroke.rds, cellchat_comparison.rds\n")
cat("  Results/     interaction_summary.csv\n")
cat("  Figures/     Fig20a_Network_Count_{Sham,Stroke}.pdf/png\n")
cat("               Fig20a_Network_Weight_{Sham,Stroke}.pdf/png\n")
cat("               Fig20a_Combined_Networks.pdf/png\n")
cat("               Fig20b_Bubble_{Sham,Stroke}.pdf/png\n")
cat("               Fig20c_Compare_Total.pdf/png\n")

cat("\n✅ 第 23 部分完成!\n")

sink(file.path(res.path, "session_info.txt"))
sessionInfo()
sink()
