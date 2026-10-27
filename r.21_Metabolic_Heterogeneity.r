
cat("\n========================================\n")
cat("21. 代谢异质性评估 - MEBOCOST\n")
cat("========================================\n\n")

rm(list = ls()); gc()
options(bitmapType = "cairo")

PROJECT_DIR <- ""
MEBOCOST_ENV <- "miniconda3/envs/mebocost"
MEBOCOST_CONF <- "MEBOCOST/mebocost.conf"

setwd(PROJECT_DIR)

output_dir <- "./21_Metabolic_Heterogeneity"
data.path  <- file.path(output_dir, "Data")
res.path   <- file.path(output_dir, "Results")
fig.path   <- file.path(output_dir, "Figures")

if (dir.exists(output_dir)) {
  cat("清理旧结果目录...\n")
  unlink(output_dir, recursive = TRUE)
}

for (p in c(output_dir, data.path, res.path, fig.path)) {
  dir.create(p, recursive = TRUE)
}

cat("\n========== 加载 R 包 ==========\n")
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
})
use_cellchat <- suppressPackageStartupMessages(
  requireNamespace("CellChat", quietly = TRUE))
if (!use_cellchat) {
  suppressPackageStartupMessages(library(igraph))
}
cat("✓ R 包加载完成 (CellChat:", use_cellchat, ")\n")

t0 <- Sys.time()

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("Step 1: Seurat → 中间数据格式\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

INPUT_RDS <- "./17_annotation/seurat_obj_annotated_clean.rds"
if (!file.exists(INPUT_RDS)) {
  stop("❌ 输入文件不存在: ", INPUT_RDS)
}

seu <- readRDS(INPUT_RDS)
cat("✓ Seurat 加载完成:", ncol(seu), "细胞 x", nrow(seu), "基因\n")

cat("细胞类型分布:\n"); print(table(seu$cell_type))
cat("分组分布:\n"); print(table(seu$condition))

cat("\n导出 counts.mtx...\n")
counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
writeMM(counts, file.path(data.path, "counts.mtx"))
writeLines(rownames(counts), file.path(data.path, "genes.tsv"))
writeLines(colnames(counts), file.path(data.path, "barcodes.tsv"))

meta <- seu@meta.data
meta$barcode <- rownames(meta)
keep_cols <- intersect(
  c("barcode", "orig.ident", "cell_type", "condition", "Age",
    "nCount_RNA", "nFeature_RNA"),
  colnames(meta)
)
write.csv(meta[, keep_cols],
          file.path(data.path, "metadata.csv"),
          row.names = FALSE)

rm(seu, counts); gc()
cat("✓ 中间数据导出完成 (counts.mtx, genes.tsv, barcodes.tsv, metadata.csv)\n")

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("Step 2: 运行 MEBOCOST (Python)\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

py_script <- file.path(data.path, "run_mebocost.py")

py_code <- sprintf('#!/usr/bin/env python
"""MEBOCOST 代谢物介导细胞通讯分析"""
import os, sys, time
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread
import anndata as ad
from mebocost import mebocost

DATA_DIR    = "%s"
RES_DIR     = "%s"
CONFIG_PATH = "%s"

print("=" * 60)
print("MEBOCOST 分析")
print("=" * 60)

t0 = time.time()

print("\\n[1/3] 加载数据并构建 AnnData ...")
counts_mtx = mmread(os.path.join(DATA_DIR, "counts.mtx")).tocsr()
genes    = pd.read_csv(os.path.join(DATA_DIR, "genes.tsv"),    header=None)[0].tolist()
barcodes = pd.read_csv(os.path.join(DATA_DIR, "barcodes.tsv"), header=None)[0].tolist()
meta = pd.read_csv(os.path.join(DATA_DIR, "metadata.csv"))
meta.index = meta["barcode"]

adata = ad.AnnData(
    X   = counts_mtx.T.tocsr(),
    obs = meta.loc[barcodes].copy(),
    var = pd.DataFrame(index=genes)
)
adata.var_names_make_unique()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.write_h5ad(os.path.join(RES_DIR, "mebocost_adata.h5ad"))
print("  AnnData shape:", adata.shape)

print("\\n[2/3] 运行 MEBOCOST 推断通讯 (thread=8, n_shuffle=1000) ...")

mebo_obj = mebocost.create_obj(
    adata           = adata,
    group_col       = "cell_type",
    condition_col   = None,
    met_est         = "mebocost",
    config_path     = CONFIG_PATH,
    species         = "mouse",
    cutoff_exp      = "auto",
    cutoff_met      = "auto",
    cutoff_prop     = 0.25,
    sensor_type     = "All",
    thread          = 8
)
mebo_obj._load_config_()
mebo_obj.estimator()

commu_res = mebo_obj.infer_commu(
    n_shuffle        = 1000,
    seed             = 12345,
    Return           = True,
    thread           = None,
    save_permuation  = False,
    min_cell_number  = 10,
    pval_method      = "permutation_test_fdr",
    pval_cutoff      = 0.05
)

print(f"  显著通讯对: {commu_res.shape[0]}")

print("\\n[3/3] 保存结果 ...")
commu_res.to_csv(os.path.join(RES_DIR, "mebocost_communication.csv"), index=False)
mebocost.save_obj(obj=mebo_obj, path=os.path.join(RES_DIR, "mebocost_obj.pk"))

sender_col   = next((c for c in commu_res.columns if "sender"   in c.lower()), "Sender")
receiver_col = next((c for c in commu_res.columns if "receiver" in c.lower()), "Receiver")
out_count = commu_res[sender_col].value_counts()
in_count  = commu_res[receiver_col].value_counts()

print("\\n作为发送方:")
print(out_count)
print("\\n作为接收方:")
print(in_count)

all_ct = sorted(set(out_count.index) | set(in_count.index))
sum_df = pd.DataFrame({
    "cell_type": all_ct,
    "outgoing": [out_count.get(ct, 0) for ct in all_ct],
    "incoming": [in_count.get(ct, 0) for ct in all_ct],
})
sum_df["total"] = sum_df["outgoing"] + sum_df["incoming"]
sum_df = sum_df.sort_values("total", ascending=False)
sum_df.to_csv(os.path.join(RES_DIR, "mebocost_cellcomm_summary.csv"), index=False)

print(f"\\n✓ MEBOCOST 完成, 耗时: {(time.time()-t0)/60:.1f} 分钟")
',
data.path, res.path, MEBOCOST_CONF)

writeLines(py_code, py_script)
cat("✓ Python 脚本已生成:", py_script, "\n")

cat("\n启动 Python (MEBOCOST 环境)...\n")
py_cmd <- file.path(MEBOCOST_ENV, "bin", "python")
if (!file.exists(py_cmd)) {
  stop("❌ Python 不存在: ", py_cmd)
}

py_exit <- system2(
  command = py_cmd,
  args    = c("-u", py_script),
  stdout  = "",
  stderr  = ""
)

if (py_exit != 0) {
  stop("❌ MEBOCOST 执行失败 (exit code ", py_exit, ")")
}

if (!file.exists(file.path(res.path, "mebocost_communication.csv"))) {
  stop("❌ MEBOCOST 结果文件缺失")
}

cat("\n✓ MEBOCOST 分析完成\n")

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("Step 3: 绘图 (Fig 18a + Fig 18b)\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

commu <- read.csv(file.path(res.path, "mebocost_communication.csv"),
                   stringsAsFactors = FALSE)
cat("✓ 通讯结果加载:", nrow(commu), "对\n")

sender_col <- "Sender"
receiver_col <- "Receiver"

cat("\n画 Fig 18a 柱状图...\n")

out_count <- commu %>%
  group_by(cell_type = !!sym(sender_col)) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(direction = "Outgoing (Sender)")

in_count <- commu %>%
  group_by(cell_type = !!sym(receiver_col)) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(direction = "Incoming (Receiver)")

comm_count <- bind_rows(out_count, in_count)
comm_count$direction <- factor(comm_count$direction,
                                 levels = c("Outgoing (Sender)", "Incoming (Receiver)"))

total_order <- comm_count %>%
  group_by(cell_type) %>%
  summarise(total = sum(n)) %>%
  arrange(desc(total))
comm_count$cell_type <- factor(comm_count$cell_type,
                                 levels = total_order$cell_type)

p_bar <- ggplot(comm_count, aes(x = cell_type, y = n, fill = direction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = c(
    "Outgoing (Sender)"  = "#E64B35",
    "Incoming (Receiver)" = "#4DBBD5"
  )) +
  labs(
    title = "Fig. 18a: Metabolite-mediated cell communication events",
    subtitle = sprintf("Total significant mCCC events: %d", nrow(commu)),
    x = "Cell type",
    y = "Number of communication events",
    fill = ""
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle   = element_text(hjust = 0.5, size = 11, color = "gray40"),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 11),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(fig.path, "Fig18a_Barplot_CommCounts.pdf"), p_bar,
       width = 10, height = 7)
ggsave(file.path(fig.path, "Fig18a_Barplot_CommCounts.png"), p_bar,
       width = 10, height = 7, dpi = 300, bg = "white")
cat("✓ Fig 18a 保存\n")

cat("\n画 Fig 18b 网络图...\n")

cell_types <- sort(unique(c(
  as.character(commu[[sender_col]]),
  as.character(commu[[receiver_col]])
)))
n_types <- length(cell_types)

count_mat <- matrix(0, nrow = n_types, ncol = n_types,
                     dimnames = list(cell_types, cell_types))
for (i in seq_len(nrow(commu))) {
  s <- as.character(commu[[sender_col]][i])
  r <- as.character(commu[[receiver_col]][i])
  count_mat[s, r] <- count_mat[s, r] + 1
}

write.csv(count_mat,
          file.path(res.path, "commu_count_matrix.csv"))

node_weights <- rowSums(count_mat) + colSums(count_mat)

if (use_cellchat) {
  library(CellChat)

  pdf(file.path(fig.path, "Fig18b_Network_CellChat.pdf"),
      width = 10, height = 10)
  par(mfrow = c(1, 1))
  netVisual_circle(
    count_mat,
    vertex.weight     = node_weights,
    weight.scale      = TRUE,
    label.edge        = FALSE,
    edge.width.max    = 8,
    title.name        = "Fig. 18b: Metabolite-mediated cell communication network"
  )
  dev.off()

  png(file.path(fig.path, "Fig18b_Network_CellChat.png"),
      width = 10, height = 10, units = "in", res = 300, bg = "white")
  par(mfrow = c(1, 1))
  netVisual_circle(
    count_mat,
    vertex.weight     = node_weights,
    weight.scale      = TRUE,
    label.edge        = FALSE,
    edge.width.max    = 8,
    title.name        = "Fig. 18b: Metabolite-mediated cell communication network"
  )
  dev.off()

  cat("✓ Fig 18b 保存 (CellChat)\n")

} else {

  edges_df <- data.frame()
  for (s in cell_types) for (r in cell_types) {
    if (s != r && count_mat[s, r] > 0) {
      edges_df <- rbind(edges_df, data.frame(
        from = s, to = r, count = count_mat[s, r]
      ))
    }
  }
  g <- graph_from_data_frame(edges_df, directed = TRUE, vertices = cell_types)
  V(g)$size <- 5 + scales::rescale(node_weights, to = c(0, 25))
  ct_colors <- setNames(
    colorRampPalette(brewer.pal(8, "Set2"))(length(cell_types)),
    cell_types
  )
  V(g)$color <- ct_colors[V(g)$name]
  V(g)$frame.color <- "white"
  V(g)$label.cex <- 0.9
  V(g)$label.color <- "black"
  E(g)$width <- scales::rescale(edges_df$count, to = c(0.5, 5))
  E(g)$arrow.size <- 0.4
  E(g)$curved <- 0.2

  set.seed(42)
  layout <- layout_in_circle(g)

  pdf(file.path(fig.path, "Fig18b_Network_CellChat.pdf"),
      width = 10, height = 10)
  par(mar = c(0, 0, 3, 0))
  plot(g, layout = layout,
       main = "Fig. 18b: Metabolite-mediated cell communication network")
  dev.off()

  png(file.path(fig.path, "Fig18b_Network_CellChat.png"),
      width = 10, height = 10, units = "in", res = 300, bg = "white")
  par(mar = c(0, 0, 3, 0))
  plot(g, layout = layout,
       main = "Fig. 18b: Metabolite-mediated cell communication network")
  dev.off()

  cat("✓ Fig 18b 保存 (igraph)\n")
}

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("       第 21 部分 MEBOCOST 完成\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

total_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("\n总耗时: %.1f 分钟\n", total_min))

cat("\n结果:\n")
cat(sprintf("  显著通讯对数: %d\n", nrow(commu)))
cat("\n输出文件:\n")
cat("  Data/        counts.mtx, genes.tsv, barcodes.tsv, metadata.csv\n")
cat("  Results/     mebocost_communication.csv\n")
cat("               mebocost_cellcomm_summary.csv\n")
cat("               commu_count_matrix.csv\n")
cat("               mebocost_obj.pk, mebocost_adata.h5ad\n")
cat("  Figures/     Fig18a_Barplot_CommCounts.pdf/png\n")
cat("               Fig18b_Network_CellChat.pdf/png\n")

cat("\n✅ 第 21 部分完成!\n")

sink(file.path(res.path, "session_info.txt"))
sessionInfo()
sink()
