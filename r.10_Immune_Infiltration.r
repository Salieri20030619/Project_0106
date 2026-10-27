
cat("\n========================================\n")
cat("10. 免疫浸润分析\n")
cat("========================================\n\n")

rm(list = ls())

options(bitmapType = "cairo")

setwd("")

output_dir <- "./10_Immune_Infiltration"
res.path   <- file.path(output_dir, "Results")
fig.path   <- file.path(output_dir, "Figures")
cor.path   <- file.path(output_dir, "Figures", "cor")

for (path in c(res.path, fig.path, cor.path)) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

cat("\n========== 加载R包 ==========\n")

pkgs <- c(
  "tidyverse", "GSVA", "ggplot2", "pheatmap",
  "RColorBrewer", "psych", "rstatix", "ggpubr",
  "ggbeeswarm", "corrplot", "igraph", "tibble", "grid"
)

for (pkg in pkgs) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org", quiet = TRUE)
    library(pkg, character.only = TRUE)
  }
  cat("✅ 已加载:", pkg, "\n")
}

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate

COLOR_CTRL   <- "#45a9b8"
COLOR_STROKE <- "#f76a56"
COLOR_LOW    <- "#355783"
COLOR_HIGH   <- "#d9b02c"
ANN_COLORS   <- list(group = c(Control = COLOR_CTRL, Stroke = COLOR_STROKE))

cat("\n========== 读取生物标志物 ==========\n")

biomarkers_file <- "./07_Expression_Validation/Results/selected_biomarkers.txt"
if (!file.exists(biomarkers_file)) {
  stop("❌ selected_biomarkers.txt 不存在，请先运行第七部分")
}
selected_biomarkers <- readLines(biomarkers_file)
selected_biomarkers <- selected_biomarkers[nzchar(selected_biomarkers) & !is.na(selected_biomarkers)]
cat("生物标志物:", paste(selected_biomarkers, collapse = ", "), "\n")

cat("\n========== 读取训练集数据 ==========\n")

train_file <- "./00_rawdata/GSE16561_for_DEG.rds"
if (!file.exists(train_file)) stop("❌ 训练集文件不存在:", train_file)

train_raw   <- readRDS(train_file)
Train_expr  <- train_raw$expression
Train_group <- train_raw$group

if (any(grepl("GSM", Train_group[[1]]))) {
  rownames(Train_group) <- Train_group[[1]]
} else if ("sample" %in% colnames(Train_group)) {
  rownames(Train_group) <- Train_group$sample
}

Train_group$group_label <- ifelse(tolower(Train_group$group) == "control", "Control", "Stroke")
cat("分组: Control =", sum(Train_group$group_label == "Control"),
    ", Stroke =", sum(Train_group$group_label == "Stroke"), "\n")

merge_dup <- function(mat) {
  genes <- rownames(mat)
  if (!any(duplicated(genes))) return(mat)
  cat("发现重复基因", sum(duplicated(genes)), "个，取均值合并\n")
  ug <- unique(genes)
  new_mat <- do.call(rbind, lapply(ug, function(g) {
    idx <- which(genes == g)
    if (length(idx) == 1) mat[idx, ] else colMeans(mat[idx, , drop = FALSE])
  }))
  rownames(new_mat) <- ug
  new_mat
}

Train_expr <- merge_dup(Train_expr)
rownames(Train_expr) <- toupper(rownames(Train_expr))

common_samples <- intersect(rownames(Train_group), colnames(Train_expr))
Train_expr  <- Train_expr[, common_samples]
Train_group <- Train_group[common_samples, ]
cat("对齐后表达矩阵:", dim(Train_expr), "\n")

group_df <- data.frame(
  sample = rownames(Train_group),
  group  = Train_group$group_label,
  stringsAsFactors = FALSE
)

cat("\n========== 读取 mmc3 标准基因集 ==========\n")

mmc3_paths <- c(
  "ssGSEA/mmc3.txt",
  "./mmc3.txt"
)
mmc3_file <- Filter(file.exists, mmc3_paths)[1]
if (is.na(mmc3_file)) stop("❌ 找不到 mmc3.txt，请确认路径")

gene_set  <- read.table(mmc3_file, header = TRUE, sep = "\t")

colnames(gene_set) <- c("gene", "cell_type", "immunity")
gene_list <- split(gene_set$gene, gene_set$cell_type)
cat("mmc3 包含细胞类型数:", length(gene_list), "\n")

gene_list_filtered <- Filter(function(g) sum(g %in% rownames(Train_expr)) >= 3, gene_list)
gene_list_filtered <- lapply(gene_list_filtered, function(g) g[g %in% rownames(Train_expr)])
cat("过滤后可用细胞类型:", length(gene_list_filtered), "\n")

if (length(gene_list_filtered) == 0) stop("❌ 没有可用的免疫细胞基因集")

cat("\n========== 运行 ssGSEA（GSVA 包）==========\n")

set.seed(123)

scores <- tryCatch({
  param <- ssgseaParam(as.matrix(Train_expr), gene_list_filtered, normalize = TRUE)
  gsva(param, verbose = FALSE)
}, error = function(e) {
  message("使用旧版 GSVA API...")
  gsva(as.matrix(Train_expr), gene_list_filtered,
       method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE)
})

write.csv(scores, file.path(res.path, "ssgsea_scores.csv"))
cat("✅ ssGSEA 计算完成，维度:", dim(scores), "\n")

violin_dat <- as.data.frame(t(scores))
violin_dat$sample <- rownames(violin_dat)
violin_dat <- merge(violin_dat, group_df, by = "sample")

violin_long <- violin_dat %>%
  pivot_longer(
    cols      = -c("sample", "group"),
    names_to  = "ImmuneCell",
    values_to = "Score"
  )

cat("\n========== 差异分析 ==========\n")

wilcox_res <- violin_long %>%
  group_by(ImmuneCell) %>%
  wilcox_test(Score ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p")

DE_wilcox <- wilcox_res %>% filter(p < 0.05)

write.csv(wilcox_res,  file.path(res.path, "wilcox_res.csv"),    row.names = FALSE)
write.csv(DE_wilcox,   file.path(res.path, "DE_wilcox_res.csv"), row.names = FALSE)
cat("显著差异细胞类型数:", nrow(DE_wilcox), "\n")

cat("\n========== 绘制热图 ==========\n")

train_group_ord <- group_df[order(group_df$group), ]
rownames(train_group_ord) <- train_group_ord$sample

annotation_col <- data.frame(
  group = train_group_ord$group,
  row.names = train_group_ord$sample
)

heatmap_mat <- scores[, rownames(train_group_ord)]
heatmap_mat[] <- apply(heatmap_mat, 2, as.numeric)

pheatmap_args <- list(
  mat = heatmap_mat,
  color             = colorRampPalette(c(COLOR_LOW, "white", COLOR_HIGH))(50),
  annotation_col    = annotation_col,
  annotation_colors = ANN_COLORS,
  clustering_method = "ward.D2",
  show_rownames     = TRUE,
  show_colnames     = FALSE,
  fontsize_col      = 5,
  cluster_cols      = FALSE,
  cluster_rows      = TRUE,
  silent            = TRUE
)

ph_obj <- do.call(pheatmap, pheatmap_args)

cairo_pdf(file.path(fig.path, "01.heatmap.pdf"),
          width = 10, height = 8, onefile = TRUE)
grid::grid.newpage()
grid::grid.draw(ph_obj$gtable)
dev.off()

png(file.path(fig.path, "01.heatmap.png"),
    width = 10, height = 8, units = "in", res = 600)
grid::grid.newpage()
grid::grid.draw(ph_obj$gtable)
dev.off()

cat("✅ 热图已保存\n")

cat("\n========== 绘制半小提琴图 ==========\n")

if (nrow(DE_wilcox) > 0) {
  plot_data <- violin_long %>%
    filter(ImmuneCell %in% DE_wilcox$ImmuneCell)

  p_half <- ggplot(plot_data, aes(x = ImmuneCell, y = Score, fill = group)) +
    geom_violin(width = 1.2, trim = FALSE, alpha = 0.4, colour = NA) +
    geom_quasirandom(width = 0.15, alpha = 0.4, size = 1.2) +
    stat_summary(fun.data = mean_cl_boot, geom = "crossbar",
                 width = 0.2, colour = "black", fatten = 1.5) +
    stat_pwc(aes(group = group),
             method = "wilcox.test",
             paired = FALSE,
             p.adjust.method = "BH",
             label = "p.signif",
             tip.length = 0,
             bracket.size = 0,
             vjust = -1) +
    scale_fill_manual(values = c(Control = COLOR_CTRL, Stroke = COLOR_STROKE), name = "Group") +
    coord_flip() +
    labs(x = NULL, y = "ssGSEA score") +
    theme_bw() +
    theme(
      plot.title          = element_text(size = 16, face = "bold", family = "Times", hjust = 0.5),
      axis.text.y         = element_text(size = 14, face = "bold", family = "Times"),
      panel.grid.major.y  = element_blank(),
      panel.grid.minor    = element_blank(),
      legend.position     = "top"
    )

  ggsave(file.path(fig.path, "02.immune_cell_halfviolin.pdf"), p_half, width = 7, height = 9)
  ggsave(file.path(fig.path, "02.immune_cell_halfviolin.png"), p_half, width = 7, height = 9, dpi = 600)
  cat("✅ 半小提琴图已保存\n")
} else {
  cat("⚠️ 无显著差异细胞类型，跳过半小提琴图\n")
}

cat("\n========== 绘制模板风格箱线图（全细胞 + 显著性 stars）==========\n")

sig_label_df <- wilcox_res %>%
  dplyr::mutate(
    star = dplyr::case_when(
      p <  0.001 ~ "***",
      p <  0.01  ~ "**",
      p <  0.05  ~ "*",
      TRUE       ~ "ns"
    )
  ) %>%
  dplyr::select(ImmuneCell, p, star)

y_max_per_cell <- violin_long %>%
  dplyr::group_by(ImmuneCell) %>%
  dplyr::summarise(y_top = max(Score, na.rm = TRUE), .groups = "drop")

sig_anno_df <- dplyr::left_join(sig_label_df, y_max_per_cell, by = "ImmuneCell")
sig_anno_df$y_star <- sig_anno_df$y_top * 1.08 + 0.02

cell_order <- violin_long %>%
  dplyr::group_by(ImmuneCell) %>%
  dplyr::summarise(m = mean(Score, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(m)) %>%
  dplyr::pull(ImmuneCell)

violin_long$ImmuneCell <- factor(violin_long$ImmuneCell, levels = cell_order)
sig_anno_df$ImmuneCell <- factor(sig_anno_df$ImmuneCell, levels = cell_order)

y_axis_max <- max(sig_anno_df$y_star, na.rm = TRUE) * 1.05

p_box <- ggplot(violin_long, aes(x = ImmuneCell, y = Score, fill = group)) +
  geom_boxplot(
    position = position_dodge(0.8),
    width = 0.65,
    outlier.size = 0.6,
    outlier.alpha = 0.5,
    linewidth = 0.35
  ) +

  geom_text(
    data = sig_anno_df,
    aes(x = ImmuneCell, y = y_star, label = star),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold",
    color = "black"
  ) +
  scale_fill_manual(
    values = c(Control = COLOR_CTRL, Stroke = COLOR_STROKE),
    name   = "Group"
  ) +
  scale_y_continuous(limits = c(NA, y_axis_max), expand = expansion(mult = c(0.02, 0.05))) +
  labs(
    title = "Immune Cell Infiltration - GSE16561 (Control vs Stroke)",
    subtitle = "Wilcoxon test: *** p<0.001, ** p<0.01, * p<0.05, ns p\u22650.05",
    x = NULL,
    y = "ssGSEA score"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle    = element_text(hjust = 0.5, size = 10, color = "gray40"),
    axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y      = element_text(size = 10),
    axis.title.y     = element_text(size = 12),
    legend.position  = "top",
    legend.title     = element_text(size = 11),
    legend.text      = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

n_cells <- length(cell_order)
fig_w <- max(10, n_cells * 0.45 + 2)

ggsave(file.path(fig.path, "02b.immune_cell_boxplot_template.pdf"),
       p_box, width = fig_w, height = 7)
ggsave(file.path(fig.path, "02b.immune_cell_boxplot_template.png"),
       p_box, width = fig_w, height = 7, dpi = 300)
cat(sprintf("✅ 模板风格箱线图已保存 (%.0fx7 英寸, %d 种细胞)\n", fig_w, n_cells))

cat("\n========== 相关性分析 ==========\n")

available_biomarkers <- intersect(selected_biomarkers, rownames(Train_expr))
if (length(available_biomarkers) == 0) stop("❌ 生物标志物均不在表达矩阵中")

keyGs_exp <- Train_expr[available_biomarkers, common_samples, drop = FALSE]
keyGs_exp <- apply(keyGs_exp, 2, as.numeric)
rownames(keyGs_exp) <- available_biomarkers

ssgsea_result_t <- as.data.frame(t(scores))
colnames(ssgsea_result_t) <- gsub("\\.", " ", colnames(ssgsea_result_t))

DE_cells <- if (nrow(DE_wilcox) > 0) unique(DE_wilcox$ImmuneCell) else colnames(ssgsea_result_t)
DE_cells  <- gsub("\\.", " ", DE_cells)
DEcells_spearman <- ssgsea_result_t[, intersect(DE_cells, colnames(ssgsea_result_t)), drop = FALSE]

immuscore <- function(gene) {
  y <- as.numeric(keyGs_exp[gene, ])
  do.call(rbind, lapply(colnames(DEcells_spearman), function(x) {
    dd <- corr.test(as.numeric(DEcells_spearman[, x]), y, method = "spearman", adjust = "fdr")
    data.frame(gene = gene, immune_cells = x, cor = dd$r, p.value = dd$p)
  }))
}

KGs_cor <- do.call(rbind, lapply(available_biomarkers, immuscore))
KGs_cor$sig <- ifelse(KGs_cor$p.value < 0.001, "***",
               ifelse(KGs_cor$p.value < 0.01,  "**",
               ifelse(KGs_cor$p.value < 0.05,  "*",  " ")))

write.csv(KGs_cor, file.path(cor.path, "05_gene_immunecell_cor.csv"), quote = FALSE, row.names = FALSE)

draw_cor_heatmap <- function(filename, w = 15, h = 5, is_pdf = TRUE) {
  p <- ggplot(KGs_cor, aes(x = immune_cells, y = gene, fill = cor)) +
    geom_tile() +
    scale_fill_gradient2(low  = COLOR_LOW,
                         mid  = "white",
                         high = COLOR_HIGH,
                         limits = c(-1, 1),
                         breaks = c(-1, -0.5, 0, 0.5, 1)) +
    labs(x = NULL, y = NULL,
         fill  = paste0(" * p < 0.05\n\n** p < 0.01\n\n*** p < 0.001\n\nCorrelation"),
         title = "Immune Cell") +
    geom_text(aes(label = sig), color = "black", size = 5) +
    theme_bw(base_size = 15) +
    theme(
      plot.title    = element_text(hjust = 0.5, size = 10, face = "bold"),
      axis.text.x   = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
      axis.text.y   = element_text(size = 10, face = "bold")
    )
  if (is_pdf) {
    pdf(filename, width = w, height = h); print(p); dev.off()
  } else {
    png(filename, height = h, width = w, units = "in", res = 600); print(p); dev.off()
  }
}

draw_cor_heatmap(file.path(cor.path, "05_gene_immunecell_cor.pdf"), is_pdf = TRUE)
draw_cor_heatmap(file.path(cor.path, "05_gene_immunecell_cor.png"), is_pdf = FALSE)
cat("✅ 基因-免疫细胞相关性热图已保存\n")

cat("\n========== 免疫细胞相关性矩阵 ==========\n")

hub_t    <- as.data.frame(t(keyGs_exp)); hub_t$sample <- rownames(hub_t)
de_cells <- DEcells_spearman; de_cells$sample <- rownames(de_cells)
cor_data <- merge(de_cells, hub_t, by = "sample")
rownames(cor_data) <- cor_data$sample
cor_data <- cor_data[, colnames(de_cells)[colnames(de_cells) != "sample"]]

write.csv(cor_data, file.path(cor.path, "06_cor_data.csv"))

Immcell_cor <- corr.test(cor_data, method = "spearman", adjust = "BH")
Immcell     <- Immcell_cor$r
p_imm       <- Immcell_cor$p

write.csv(Immcell, file.path(cor.path, "07_ImmuneCell_deg_cor.csv"))
write.csv(p_imm,   file.path(cor.path, "07_ImmuneCell_deg_pvalue.csv"))

col_grad <- colorRampPalette(c(COLOR_LOW, "#FF9999", COLOR_HIGH), space = "Lab")

draw_corrplot <- function(filename, w = 8, h = 10, is_pdf = TRUE) {
  if (is_pdf) pdf(filename, width = w, height = h, family = "Times")
  else        png(filename, width = w, height = h, units = "in", res = 600)

  corrplot(t(Immcell),
           p.mat       = p_imm,
           method      = "pie",
           type        = "upper",
           tl.pos      = "lt", tl.cex = 1, tl.col = "black",
           tl.srt      = 45, tl.offset = 0.5,
           insig       = "label_sig",
           sig.level   = c(0.001, 0.01, 0.05),
           pch.cex     = 1,
           col         = col_grad(20))
  corrplot(t(Immcell),
           type        = "lower",
           add         = TRUE,
           method      = "number",
           tl.pos      = "n",
           cl.pos      = "n",
           diag        = FALSE,
           number.cex  = 0.7,
           col         = col_grad(20))
  dev.off()
}

draw_corrplot(file.path(cor.path, "08.Immcell_cor.pdf"), is_pdf = TRUE)
draw_corrplot(file.path(cor.path, "08.Immcell_cor.png"), is_pdf = FALSE)
cat("✅ 免疫细胞相关性矩阵已保存\n")

save.image(file.path(res.path, "ssGSEA_all.RData"))

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("             第十部分：免疫浸润分析完成\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")
cat("生物标志物:", paste(available_biomarkers, collapse = ", "), "\n")
cat("免疫细胞类型（mmc3）:", length(gene_list_filtered), "\n")
cat("显著差异细胞:", nrow(DE_wilcox), "\n\n")
cat("输出目录:", output_dir, "\n")
cat("  Results/ssgsea_scores.csv\n")
cat("  Results/wilcox_res.csv / DE_wilcox_res.csv\n")
cat("  Figures/01.heatmap.pdf/png\n")
cat("  Figures/02.immune_cell_halfviolin.pdf/png\n")
cat("  Figures/cor/05_gene_immunecell_cor.pdf/png\n")
cat("  Figures/cor/08.Immcell_cor.pdf/png\n")
cat("\n✅ 第十部分完成！\n")

sink(file.path(res.path, "session_info.txt"))
sessionInfo()
sink()
