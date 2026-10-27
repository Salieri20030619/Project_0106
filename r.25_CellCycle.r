
cat("\n========================================\n")
cat("25. 细胞周期鉴别分析\n")
cat("========================================\n\n")

rm(list = ls()); gc()
options(bitmapType = "cairo")

setwd("")

output_dir <- "./25_CellCycle"
res.path   <- file.path(output_dir, "Results")
fig.path   <- file.path(output_dir, "Figures")

if (dir.exists(output_dir)) {
  cat("清理旧结果目录...\n")
  unlink(output_dir, recursive = TRUE)
}

for (p in c(output_dir, res.path, fig.path)) {
  dir.create(p, recursive = TRUE)
}

cat("\n========== 加载 R 包 ==========\n")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
})
cat("✓ R 包加载完成\n")

t0 <- Sys.time()

cat("\n========== 1. 加载 Seurat 对象 ==========\n")

INPUT_RDS <- "./17_annotation/seurat_obj_annotated_clean.rds"
if (!file.exists(INPUT_RDS)) stop("❌ 输入文件不存在: ", INPUT_RDS)

seu <- readRDS(INPUT_RDS)
cat("✓ 加载完成:", ncol(seu), "细胞 x", nrow(seu), "基因\n")
cat("细胞类型分布:\n"); print(table(seu$cell_type))

cat("\n========== 2. CellCycleScoring ==========\n")

s.genes   <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

cat("S 期基因数:", length(s.genes), "\n")
cat("G2M 期基因数:", length(g2m.genes), "\n")

to_mouse <- function(x) {
  paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
}
s.genes_mouse   <- to_mouse(s.genes)
g2m.genes_mouse <- to_mouse(g2m.genes)

s.genes_use   <- intersect(s.genes_mouse,   rownames(seu))
g2m.genes_use <- intersect(g2m.genes_mouse, rownames(seu))

cat("匹配到的 S 期基因:",   length(s.genes_use), "/", length(s.genes_mouse), "\n")
cat("匹配到的 G2M 期基因:", length(g2m.genes_use), "/", length(g2m.genes_mouse), "\n")

seu <- CellCycleScoring(
  object       = seu,
  s.features   = s.genes_use,
  g2m.features = g2m.genes_use,
  set.ident    = FALSE
)

cat("✓ 细胞周期打分完成\n")
cat("\n周期分布:\n")
print(table(seu$Phase))
cat("\n分组 × 周期:\n")
if ("condition" %in% colnames(seu@meta.data)) {
  print(table(seu$Phase, seu$condition))
}

cc_scores <- data.frame(
  barcode   = colnames(seu),
  cell_type = as.character(seu$cell_type),
  condition = as.character(seu$condition),
  S.Score   = seu$S.Score,
  G2M.Score = seu$G2M.Score,
  Phase     = as.character(seu$Phase)
)
write.csv(cc_scores, file.path(res.path, "cell_cycle_scores.csv"), row.names = FALSE)

cat("\n========== 3. Fig 22a UMAP 按周期着色 ==========\n")

if (!"umap" %in% Reductions(seu)) {
  cat("⚠ 未找到 UMAP, 重新计算...\n")
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE)
  seu <- ScaleData(seu, verbose = FALSE)
  seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
  seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)
}

phase_colors <- c("G1" = "#4DBBD5", "S" = "#E64B35", "G2M" = "#00A087")

p_umap <- DimPlot(seu, reduction = "umap", group.by = "Phase",
                   cols = phase_colors, pt.size = 0.3) +
  labs(title = "Fig. 22a: Cell cycle phase distribution") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

ggsave(file.path(fig.path, "Fig22a_UMAP_CellCycle.pdf"), p_umap,
       width = 9, height = 7)
ggsave(file.path(fig.path, "Fig22a_UMAP_CellCycle.png"), p_umap,
       width = 9, height = 7, dpi = 300, bg = "white")
cat("✓ Fig 22a 保存\n")

cat("\n========== 4. Fig 22b 周期比例柱状图 ==========\n")

phase_df <- cc_scores %>%
  group_by(cell_type, Phase) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cell_type) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

phase_df$Phase <- factor(phase_df$Phase, levels = c("G1", "S", "G2M"))

ct_order <- cc_scores %>%
  group_by(cell_type) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  pull(cell_type)
phase_df$cell_type <- factor(phase_df$cell_type, levels = ct_order)

phase_wide <- phase_df %>%
  select(cell_type, Phase, percent) %>%
  pivot_wider(names_from = Phase, values_from = percent, values_fill = 0)
write.csv(phase_wide,
          file.path(res.path, "cycle_celltype_percent.csv"),
          row.names = FALSE)

p_bar <- ggplot(phase_df, aes(x = cell_type, y = percent, fill = Phase)) +
  geom_bar(stat = "identity", position = "stack", color = "white",
           linewidth = 0.2, width = 0.75) +
  scale_fill_manual(values = phase_colors) +
  labs(
    title = "Fig. 22b: Cell cycle phase composition by cell type",
    x = "Cell type",
    y = "Percentage (%)",
    fill = "Phase"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x   = element_text(angle = 45, hjust = 1, size = 11),
    legend.position = "top"
  )

ggsave(file.path(fig.path, "Fig22b_CellCycle_Percent_Bar.pdf"), p_bar,
       width = 10, height = 7)
ggsave(file.path(fig.path, "Fig22b_CellCycle_Percent_Bar.png"), p_bar,
       width = 10, height = 7, dpi = 300, bg = "white")
cat("✓ Fig 22b 保存\n")

cat("\n========== 5. Fig 22c 周期基因表达点图 ==========\n")

rep_s_genes <- intersect(c("Mcm6", "Pcna", "Tyms", "Fen1", "Mcm2", "Mcm4", "Mcm5"),
                          rownames(seu))
rep_g2m_genes <- intersect(c("Top2a", "Mki67", "Ccnb2", "Ccnb1", "Cdk1",
                               "Aurka", "Aurkb", "Birc5"),
                            rownames(seu))
rep_genes <- c(rep_s_genes, rep_g2m_genes)
cat("代表性周期基因数:", length(rep_genes), "\n")
cat("  S 期:",   paste(rep_s_genes, collapse = ", "), "\n")
cat("  G2M 期:", paste(rep_g2m_genes, collapse = ", "), "\n")

Idents(seu) <- factor(seu$cell_type, levels = ct_order)

p_dot <- DotPlot(seu,
                  features = rep_genes,
                  cols = c("lightgrey", "#E64B35"),
                  dot.scale = 6) +
  RotatedAxis() +
  labs(
    title = "Fig. 22c: Expression of representative cell cycle genes",
    x = "", y = "Cell type"
  ) +
  theme(
    plot.title  = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(fig.path, "Fig22c_CellCycle_Genes_DotPlot.pdf"), p_dot,
       width = 11, height = 7)
ggsave(file.path(fig.path, "Fig22c_CellCycle_Genes_DotPlot.png"), p_dot,
       width = 11, height = 7, dpi = 300, bg = "white")
cat("✓ Fig 22c 保存\n")

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("       第 25 部分 细胞周期 完成\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

total_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("\n总耗时: %.1f 分钟\n", total_min))

cat("\n周期分布汇总:\n")
print(table(seu$Phase))
cat("\n各细胞类型周期比例 (前 5 行):\n")
print(head(phase_wide, 5))

cat("\n输出文件:\n")
cat("  Results/\n")
cat("    cell_cycle_scores.csv\n")
cat("    cycle_celltype_percent.csv\n")
cat("  Figures/\n")
cat("    Fig22a_UMAP_CellCycle.pdf/png\n")
cat("    Fig22b_CellCycle_Percent_Bar.pdf/png\n")
cat("    Fig22c_CellCycle_Genes_DotPlot.pdf/png\n")

cat("\n✅ 第 25 部分完成!\n")

sink(file.path(res.path, "session_info.txt"))
sessionInfo()
sink()
