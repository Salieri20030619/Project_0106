
cat("\n========================================\n")
cat("04. 脾脏组织免疫浸润分析 (22种细胞)\n")
cat("========================================\n\n")

rm(list = ls())

setwd("")

spleen_file <- "./00_rawdata/spleen_ready_mouse.rds"

candidate_dir <- "./02_CandidateGene_Intersect/Results"
candidate_file <- file.path(candidate_dir, "final_candidate_genes.txt")
candidate_detail_file <- file.path(candidate_dir, "candidate_genes_detail.csv")

output_dir <- "./04_Spleen_Immune_Infiltration_Mouse"
fig_dir <- file.path(output_dir, "figures")
res_dir <- file.path(output_dir, "results")

for (dir_path in c(output_dir, fig_dir, res_dir)) {
  if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)
}

cat("候选基因目录:", candidate_dir, "\n")
cat("输出目录:", output_dir, "\n\n")

cat("加载R包...\n")
library(GSVA)
library(GSEABase)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(corrplot)
library(ggpubr)
library(tibble)
cat("✅ R包加载完成\n\n")

cat("========== 1. 读取数据 ==========\n")
if (file.exists(spleen_file)) {
  spleen <- readRDS(spleen_file)
  spleen_expr <- spleen$expression
  spleen_group <- spleen$group
  cat("✅ 小鼠脾脏组织数据读取成功\n")
  cat("  基因数:", nrow(spleen_expr), "\n")
  cat("  样本数:", ncol(spleen_expr), "\n")
  cat("  分组统计:\n")
  print(table(spleen_group$group))
} else {
  stop("❌ 小鼠脾脏组织数据不存在: ", spleen_file)
}

cat("\n========== 2. 提取候选基因 ==========\n")

candidate_genes <- character(0)
candidate_up <- character(0)
candidate_down <- character(0)

if (dir.exists(candidate_dir)) {
  if (file.exists(candidate_file)) {
    candidate_genes <- readLines(candidate_file)
    candidate_genes <- candidate_genes[candidate_genes != "" & !is.na(candidate_genes)]
    cat("✅ 从 final_candidate_genes.txt 读取到", length(candidate_genes), "个候选基因\n")
  }
  if (file.exists(candidate_detail_file)) {
    candidate_detail <- read.csv(candidate_detail_file, stringsAsFactors = FALSE)
    if (nrow(candidate_detail) > 0) {
      candidate_genes <- candidate_detail$Gene
      candidate_up <- candidate_detail$Gene[candidate_detail$Regulation == "Up"]
      candidate_down <- candidate_detail$Gene[candidate_detail$Regulation == "Down"]
      cat("  上调基因:", length(candidate_up), "  下调基因:", length(candidate_down), "\n")
    }
  }
}

if (length(candidate_genes) == 0) {
  cat("⚠️ 使用备用候选基因列表\n")
  candidate_genes <- c("C5AR1","PDK4","CRISPLD2","ADM","CYP1B1","CEBPD",
                       "ABCA1","IRS2","FCGR2A","HSPA1A","KLF4","CD14",
                       "EMILIN2","FKBP5","LRG1","HP","LCN2","LGALS3","MFGE8")
  candidate_up <- setdiff(candidate_genes, "MFGE8")
  candidate_down <- "MFGE8"
}

cat("\n候选基因:", paste(candidate_genes, collapse = ", "), "\n\n")

cat("========== 3. 基因名转换 ==========\n")
rownames(spleen_expr) <- toupper(rownames(spleen_expr))
candidate_genes <- toupper(candidate_genes)
candidate_up <- toupper(candidate_up)
candidate_down <- toupper(candidate_down)

available_genes <- candidate_genes[candidate_genes %in% rownames(spleen_expr)]
missing_genes <- setdiff(candidate_genes, rownames(spleen_expr))

cat("✅ 可用候选基因:", length(available_genes), "/", length(candidate_genes), "\n")
if (length(missing_genes) > 0) {
  cat("⚠️ 缺失:", paste(missing_genes, collapse = ", "), "\n")
}

cat("\n========== 4. 准备22种免疫细胞标记基因集 ==========\n")

immune_genes <- list(
  "B_cells_naive" = c("MS4A1","CD19","CD79A","CD79B","BANK1","TCL1A","FCER2","IGHD"),
  "B_cells_memory" = c("CD27","CD19","CD79A","MS4A1","TNFRSF13B","AIM2"),
  "Plasma_cells" = c("MZB1","JCHAIN","XBP1","DERL3","CD38","TNFRSF17","SDC1","IGHG1"),
  "T_cells_CD8" = c("CD8A","CD8B","GZMA","GZMB","PRF1","CD3D","CD3E","NKG7"),
  "T_cells_CD4_naive" = c("CD4","CCR7","SELL","LEF1","TCF7","IL7R","CD3D","CD3E"),
  "T_cells_CD4_memory_resting" = c("CD4","IL7R","CD3D","CD3E","CCR7","S100A4","CD40LG"),
  "T_cells_CD4_memory_activated" = c("CD4","CD3D","CD3E","CD69","IFNG","IL2RA","CXCR3"),
  "T_cells_follicular_helper" = c("CXCR5","ICOS","PDCD1","BCL6","CD4","IL21","CXCL13"),
  "T_cells_regulatory_Tregs" = c("FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18","TIGIT"),
  "T_cells_gamma_delta" = c("TRGC1","TRGC2","TRDC","TRGV9","CD3D","CD3E","NKG7"),
  "NK_cells_resting" = c("NCAM1","KLRD1","NCR1","KLRC1","XCL1","XCL2"),
  "NK_cells_activated" = c("NCAM1","NCR1","KLRD1","GZMB","PRF1","FCGR3A","FGFBP2","SPON2"),
  "Monocytes" = c("CD14","FCGR3A","CSF1R","LYZ","S100A8","S100A9","VCAN","CCR2"),
  "Macrophages_M0" = c("CD68","CD163","MSR1","MRC1","CSF1R","ITGAM"),
  "Macrophages_M1" = c("CD68","NOS2","TNF","IL6","CXCL10","CCL2","FCGR1A","CD86"),
  "Macrophages_M2" = c("CD163","MRC1","ARG1","IL10","TGFB1","CCL22","MSR1","STAB1"),
  "Dendritic_cells_resting" = c("CD1C","ITGAX","FLT3","IRF8","ZBTB46","CLEC10A","CLEC9A","CD1A"),
  "Dendritic_cells_activated" = c("CD83","CCR7","LAMP3","CD80","CD86","CCL19","CCL22","CD40"),
  "Mast_cells_resting" = c("KIT","MS4A2","FCER1A","CPA3","CMA1","TPSAB1","TPSB2"),
  "Mast_cells_activated" = c("KIT","TPSAB1","TPSB2","CMA1","HDC","IL4","TNF"),
  "Eosinophils" = c("SIGLEC8","IL5RA","CCR3","PRG2","EPX","RNASE2","RNASE3","CLC"),
  "Neutrophils" = c("FCGR3B","CSF3R","FPR1","S100A8","S100A9","MPO","ELANE","MMP9","LCN2","CXCR2")
)

filtered_genes <- list()
for (cell_type in names(immune_genes)) {
  genes <- immune_genes[[cell_type]]
  available <- genes[genes %in% rownames(spleen_expr)]
  if (length(available) >= 3) {
    filtered_genes[[cell_type]] <- available
    cat(sprintf("  ✅ %-32s %d/%d 可用\n", cell_type, length(available), length(genes)))
  } else {
    cat(sprintf("  ⚠️ %-32s 仅 %d/%d 可用，跳过\n", cell_type, length(available), length(genes)))
  }
}

cat("\n✅ 最终用于分析的免疫细胞类型:", length(filtered_genes), "种\n")

if (length(filtered_genes) == 0) {
  stop("❌ 没有可用的免疫细胞标记基因集")
}

cat("\n========== 5. 准备表达矩阵 ==========\n")
expr_matrix <- as.matrix(spleen_expr)
expr_matrix <- expr_matrix[rowSums(expr_matrix) > 0, ]
expr_matrix <- log2(expr_matrix + 1)
cat("  表达矩阵维度:", dim(expr_matrix), "\n")

cat("\n========== 6. 运行GSVA-ssGSEA分析 ==========\n")
gene_sets <- lapply(names(filtered_genes), function(name) {
  GeneSet(setName = name, geneIds = filtered_genes[[name]])
})
gene_sets <- GeneSetCollection(gene_sets)

ssgsea_param <- ssgseaParam(
  exprData = expr_matrix,
  geneSets = gene_sets,
  alpha = 0.25,
  normalize = TRUE
)
ssgsea_result <- gsva(ssgsea_param, verbose = TRUE)

immune_scores <- as.data.frame(t(ssgsea_result))
immune_scores$Group <- spleen_group$group
immune_scores$Sample <- rownames(immune_scores)

cat("\n✅ 免疫浸润评分计算完成\n")
cat("  免疫细胞类型数:", ncol(immune_scores) - 2, "\n")
cat("  样本数:", nrow(immune_scores), "\n")

write.csv(immune_scores, file.path(res_dir, "immune_infiltration_scores.csv"), row.names = FALSE)

immune_scores[is.na(immune_scores)] <- 0
immune_cols <- colnames(immune_scores)[!colnames(immune_scores) %in% c("Group","Sample")]

cat("\n========== 8. 差异分析 ==========\n")

control_scores <- immune_scores[immune_scores$Group == "control", immune_cols, drop = FALSE]
maco_scores <- immune_scores[immune_scores$Group == "MACO", immune_cols, drop = FALSE]

cat("Control组样本数:", nrow(control_scores), "\n")
cat("MACO组样本数:", nrow(maco_scores), "\n")

control_means <- colMeans(control_scores, na.rm = TRUE)
maco_means <- colMeans(maco_scores, na.rm = TRUE)

diff_results <- data.frame(
  CellType = immune_cols,
  Control_Mean = control_means,
  MACO_Mean = maco_means,
  Difference = maco_means - control_means,
  Pvalue = NA,
  stringsAsFactors = FALSE
)

for (i in seq_along(immune_cols)) {
  test <- wilcox.test(control_scores[, i], maco_scores[, i], exact = FALSE)
  diff_results$Pvalue[i] <- test$p.value
}

diff_results$FDR <- p.adjust(diff_results$Pvalue, method = "BH")
diff_results$Significant <- ifelse(
  diff_results$Pvalue < 0.05,
  ifelse(diff_results$Difference > 0, "Up", "Down"),
  "NS"
)
diff_results <- diff_results[order(diff_results$Pvalue), ]

cat("\n差异分析结果:\n")
print(diff_results)

sig_cells <- diff_results$CellType[diff_results$Pvalue < 0.05]
cat("\n显著差异的免疫细胞 (p<0.05):", length(sig_cells), "种\n")
if (length(sig_cells) > 0) print(sig_cells)

write.csv(diff_results, file.path(res_dir, "differential_immune_cells.csv"), row.names = FALSE)

cat("\n========== 9. 开始绘图 ==========\n")

heat_colors <- colorRampPalette(c("#4DBBD5","white","#E64B35"))(100)
bar_colors <- c("control" = "#4DAF4A", "MACO" = "#E41A1C")

cat("\n1. Fig.2a 免疫浸润热图...\n")
heatmap_matrix <- t(as.matrix(immune_scores[, immune_cols, drop = FALSE]))
colnames(heatmap_matrix) <- immune_scores$Sample

annotation_col <- data.frame(Group = immune_scores$Group)
rownames(annotation_col) <- immune_scores$Sample
ann_colors <- list(Group = c(control = "#4DAF4A", MACO = "#E41A1C"))

heatmap_args <- list(
  mat = heatmap_matrix,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  main = "Immune Cell Infiltration - Mouse Spleen",
  color = heat_colors,
  scale = "row",
  show_rownames = TRUE, show_colnames = TRUE,
  fontsize = 12,
  fontsize_row = 12,
  fontsize_col = 11,
  cellwidth = 28,
  cellheight = 22,
  angle_col = 45,
  cluster_rows = TRUE, cluster_cols = TRUE,
  border_color = "gray85"
)

pdf(file.path(fig_dir, "Fig2a_immune_heatmap.pdf"), width = 14, height = 11)
do.call(pheatmap, heatmap_args)
dev.off()

png(file.path(fig_dir, "Fig2a_immune_heatmap.png"),
    width = 14, height = 11, units = "in", res = 300)
do.call(pheatmap, heatmap_args)
dev.off()
cat("  ✅ Fig2a_immune_heatmap.pdf/png\n")

cat("\n2. Fig.2a 免疫浸润柱状图...\n")
plot_data <- immune_scores %>%
  dplyr::select(-Sample) %>%
  tidyr::pivot_longer(cols = -Group, names_to = "CellType", values_to = "Score")

plot_summary <- plot_data %>%
  dplyr::group_by(CellType, Group) %>%
  dplyr::summarise(
    Mean = mean(Score, na.rm = TRUE),
    SD = sd(Score, na.rm = TRUE),
    .groups = "drop"
  )

control_means_order <- plot_summary %>%
  dplyr::filter(Group == "control") %>%
  dplyr::arrange(desc(Mean)) %>%
  dplyr::pull(CellType)

plot_summary$CellType <- factor(plot_summary$CellType, levels = control_means_order)

p_bar <- ggplot(plot_summary, aes(x = CellType, y = Mean, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.75) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                position = position_dodge(0.9), width = 0.25, linewidth = 0.5) +
  scale_fill_manual(values = bar_colors) +
  labs(title = "Immune Cell Infiltration - Mouse Spleen",
       x = NULL, y = "ssGSEA Score") +
  theme_bw(base_size = 13) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.text.y  = element_text(size = 12),
        axis.title.y = element_text(size = 13),
        plot.title   = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.position = "top",
        legend.text  = element_text(size = 12),
        legend.title = element_text(size = 12),
        panel.grid.minor = element_blank(),
        plot.margin = margin(10, 15, 10, 15))

ggsave(file.path(fig_dir, "Fig2a_immune_barplot.pdf"), p_bar, width = 15, height = 8)
ggsave(file.path(fig_dir, "Fig2a_immune_barplot.png"), p_bar, width = 15, height = 8, dpi = 300)
cat("  ✅ Fig2a_immune_barplot.pdf/png\n")

if (length(sig_cells) > 0) {
  cat("\n3. Fig.2b 差异免疫细胞箱线图...\n")
  plot_data_subset <- plot_data[plot_data$CellType %in% sig_cells, ]
  plot_data_subset$CellType <- factor(plot_data_subset$CellType, levels = sig_cells)

  p_box <- ggplot(plot_data_subset, aes(x = CellType, y = Score, fill = Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.1), size = 1, alpha = 0.5) +
    scale_fill_manual(values = bar_colors) +
    labs(title = "Differential Immune Cell Infiltration",
         x = "", y = "ssGSEA Score") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    stat_compare_means(aes(group = Group), label = "p.signif",
                       method = "wilcox.test", hide.ns = TRUE)

  fig_w <- max(8, length(sig_cells) * 0.7)
  ggsave(file.path(fig_dir, "Fig2b_differential_immune_cells.pdf"), p_box,
         width = fig_w, height = 7)
  ggsave(file.path(fig_dir, "Fig2b_differential_immune_cells.png"), p_box,
         width = fig_w, height = 7, dpi = 300)
  cat("  ✅ Fig2b_differential_immune_cells.pdf/png\n")
} else {
  cat("\n3. 无显著差异细胞，跳过箱线图\n")
}

cat("\n4. Fig.2c 免疫细胞相关性热图...\n")
immune_cor <- cor(immune_scores[, immune_cols], method = "spearman")
write.csv(immune_cor, file.path(res_dir, "immune_cell_correlation.csv"))

pdf(file.path(fig_dir, "Fig2c_immune_cell_correlation.pdf"), width = 11, height = 10)
corrplot(immune_cor, method = "color", type = "upper", order = "hclust",
         tl.col = "black", tl.srt = 45, tl.cex = 0.75,
         col = colorRampPalette(c("#4DBBD5","white","#E64B35"))(100),
         title = "Immune Cell Correlation - Mouse Spleen",
         mar = c(0, 0, 2, 0))
dev.off()

png(file.path(fig_dir, "Fig2c_immune_cell_correlation.png"), width = 11, height = 10, units = "in", res = 300)
corrplot(immune_cor, method = "color", type = "upper", order = "hclust",
         tl.col = "black", tl.srt = 45, tl.cex = 0.75,
         col = colorRampPalette(c("#4DBBD5","white","#E64B35"))(100),
         title = "Immune Cell Correlation - Mouse Spleen",
         mar = c(0, 0, 2, 0))
dev.off()
cat("  ✅ Fig2c_immune_cell_correlation.pdf/png\n")

cat("\n5. Fig.2d 候选基因与免疫细胞相关性...\n")
if (length(available_genes) > 0 && length(immune_cols) > 0) {

  gene_expr <- as.data.frame(t(spleen_expr[available_genes, , drop = FALSE]))
  colnames(gene_expr) <- available_genes
  gene_expr$Sample <- rownames(gene_expr)

  merge_data <- merge(gene_expr, immune_scores[, c("Sample", immune_cols)], by = "Sample")

  gene_immune_cor <- data.frame()
  for (gene in available_genes) {
    for (cell in immune_cols) {
      test <- suppressWarnings(
        cor.test(merge_data[[gene]], merge_data[[cell]], method = "spearman", exact = FALSE)
      )
      gene_immune_cor <- rbind(gene_immune_cor, data.frame(
        Gene = gene, CellType = cell,
        Correlation = test$estimate, Pvalue = test$p.value,
        stringsAsFactors = FALSE
      ))
    }
  }

  gene_immune_cor$FDR <- p.adjust(gene_immune_cor$Pvalue, method = "BH")
  gene_immune_cor$Significant <- ifelse(
    abs(gene_immune_cor$Correlation) > 0.3 & gene_immune_cor$Pvalue < 0.05,
    ifelse(gene_immune_cor$Correlation > 0, "Positive", "Negative"),
    "NS"
  )
  write.csv(gene_immune_cor, file.path(res_dir, "gene_immune_correlation.csv"), row.names = FALSE)

  if (length(unique(gene_immune_cor$Gene)) >= 2 && length(unique(gene_immune_cor$CellType)) >= 2) {
    cor_matrix <- gene_immune_cor %>%
      dplyr::select(Gene, CellType, Correlation) %>%
      tidyr::pivot_wider(names_from = CellType, values_from = Correlation, values_fill = 0) %>%
      tibble::column_to_rownames("Gene") %>%
      as.matrix()
    mode(cor_matrix) <- "numeric"

    annotation_row <- data.frame(
      Regulation = ifelse(rownames(cor_matrix) %in% candidate_up, "Up",
                   ifelse(rownames(cor_matrix) %in% candidate_down, "Down", "Other"))
    )
    rownames(annotation_row) <- rownames(cor_matrix)
    ann_row_colors <- list(Regulation = c(Up = "#E64B35", Down = "#4DBBD5", Other = "gray80"))

    n_genes <- nrow(cor_matrix)
    n_cells <- ncol(cor_matrix)
    plot_width <- max(12, n_cells * 0.55)
    plot_height <- max(8, n_genes * 0.4)

    pdf(file.path(fig_dir, "Fig2d_gene_immune_correlation.pdf"),
        width = plot_width, height = plot_height)
    pheatmap(cor_matrix,
             display_numbers = TRUE, number_format = "%.2f", fontsize_number = 7,
             color = colorRampPalette(c("#4DBBD5","white","#E64B35"))(100),
             breaks = seq(-1, 1, length.out = 101),
             main = "Gene-Immune Cell Correlation",
             fontsize_row = 10, fontsize_col = 8, angle_col = 45,
             annotation_row = annotation_row,
             annotation_colors = ann_row_colors)
    dev.off()

    png(file.path(fig_dir, "Fig2d_gene_immune_correlation.png"),
        width = plot_width, height = plot_height, units = "in", res = 300)
    pheatmap(cor_matrix,
             display_numbers = TRUE, number_format = "%.2f", fontsize_number = 7,
             color = colorRampPalette(c("#4DBBD5","white","#E64B35"))(100),
             breaks = seq(-1, 1, length.out = 101),
             main = "Gene-Immune Cell Correlation",
             fontsize_row = 10, fontsize_col = 8, angle_col = 45,
             annotation_row = annotation_row,
             annotation_colors = ann_row_colors)
    dev.off()
    cat("  ✅ Fig2d_gene_immune_correlation.pdf/png\n")
  }
}

cat("\n========== 10. 保存分析结果 ==========\n")
saveRDS(list(
  immune_scores = immune_scores,
  diff_results = diff_results,
  immune_cor = immune_cor,
  gene_immune_cor = if (exists("gene_immune_cor")) gene_immune_cor else NULL,
  filtered_genes = filtered_genes,
  available_genes = available_genes,
  candidate_genes = candidate_genes,
  candidate_up = candidate_up,
  candidate_down = candidate_down
), file = file.path(res_dir, "immune_infiltration_results.rds"))

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("结果汇总\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat(sprintf("\n免疫细胞类型: %d 种\n", length(immune_cols)))
cat(sprintf("显著差异细胞 (p<0.05): %d 种\n\n", length(sig_cells)))

if (length(sig_cells) > 0) {
  for (cell in sig_cells) {
    idx <- which(diff_results$CellType == cell)
    cat(sprintf("  %s: Ctrl=%.4f MACO=%.4f Δ=%.4f p=%.4f %s\n",
                cell, diff_results$Control_Mean[idx], diff_results$MACO_Mean[idx],
                diff_results$Difference[idx], diff_results$Pvalue[idx],
                diff_results$Significant[idx]))
  }
}

cat("\n📁 生成文件:\n")
for (f in sort(list.files(fig_dir))) {
  size <- round(file.size(file.path(fig_dir, f)) / 1024, 1)
  cat(sprintf("  ✅ %-45s %8.1f KB\n", f, size))
}

cat("\n✅ 完成！\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
