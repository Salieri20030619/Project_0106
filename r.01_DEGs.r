
cat("\n========================================\n")
cat("01. 差异表达基因分析\n")
cat("========================================\n\n")

rm(list = ls())

cat("\n========== 初始化 ==========\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

setwd("")

output_dir <- "./01_DEGs"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

fig.path <- file.path(output_dir, "Figures")
res.path <- file.path(output_dir, "Results")
code.path <- file.path(output_dir, "Codes")

for (path in c(fig.path, res.path, code.path)) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

cat("输出目录:", output_dir, "\n")
cat("图片目录:", fig.path, "\n")
cat("结果目录:", res.path, "\n\n")

cat("\n========== 加载R包 ==========\n")

packages <- c(
  "tidyverse", "ggplot2", "ggrepel", "pheatmap", "RColorBrewer",
  "DESeq2", "limma", "edgeR", "dplyr", "tidyr"
)

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg %in% c("DESeq2", "limma", "edgeR")) {
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg, repos = "http://cran.us.r-project.org")
    }
    library(pkg, character.only = TRUE)
  }
}

cat("✅ 所有R包加载完成\n")

plot_volcano_limma <- function(deg_df, sig_genes, title, filename, threshold_text) {

  deg_df$logFC_plot <- deg_df$logFC
  deg_df$p_plot <- deg_df$P.Value
  deg_df$neg_log10_p <- -log10(deg_df$p_plot + 1e-300)

  sig_up <- sig_genes$Gene[sig_genes$logFC > 0]
  sig_down <- sig_genes$Gene[sig_genes$logFC < 0]

  top10_up <- head(sig_up, 10)
  top10_down <- head(sig_down, 10)
  top10_genes <- c(top10_up, top10_down)
  deg_df$Label <- ifelse(deg_df$Gene %in% top10_genes, deg_df$Gene, "")

  deg_df$Regulation <- "Not Significant"
  deg_df$Regulation[deg_df$Gene %in% sig_up] <- "Up-regulated"
  deg_df$Regulation[deg_df$Gene %in% sig_down] <- "Down-regulated"

  up_count <- length(sig_up)
  down_count <- length(sig_down)

  p <- ggplot(deg_df, aes(x = logFC_plot, y = neg_log10_p, color = Regulation)) +
    geom_point(size = 1.5, alpha = 0.6) +
    scale_color_manual(
      values = c("Up-regulated" = "#E64B35", "Down-regulated" = "#4DBBD5", "Not Significant" = "gray70"),
      labels = c(
        "Down-regulated" = paste0("Down-regulated (Blue) (n=", down_count, ")"),
        "Not Significant" = "Not Significant (Gray)",
        "Up-regulated"   = paste0("Up-regulated (Red) (n=", up_count, ")")
      )
    ) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray50", alpha = 0.8) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    geom_text_repel(aes(label = Label), size = 3, max.overlaps = 15, box.padding = 0.5) +
    theme_minimal() +
    labs(title = paste0(title, "\nThreshold: ", threshold_text),
         x = "log2(Fold Change)", y = "-log10(p-value)", color = "Regulation") +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          legend.position = "bottom")

  ggsave(file.path(fig.path, filename), p, width = 8, height = 7)
  png_filename <- gsub("\\.pdf$", ".png", filename)
  ggsave(file.path(fig.path, png_filename), p, width = 8, height = 7, dpi = 300)
  cat("  ✅ 已保存:", filename, "和", png_filename, "\n")
  return(p)
}

plot_volcano_deseq2 <- function(deg_df, sig_genes, title, filename, threshold_text) {

  deg_df$logFC_plot <- deg_df$log2FoldChange
  deg_df$p_plot <- deg_df$pvalue
  deg_df$neg_log10_p <- -log10(deg_df$p_plot + 1e-300)

  sig_up <- sig_genes$Gene[sig_genes$log2FoldChange > 0]
  sig_down <- sig_genes$Gene[sig_genes$log2FoldChange < 0]

  top10_up <- head(sig_up, 10)
  top10_down <- head(sig_down, 10)
  top10_genes <- c(top10_up, top10_down)
  deg_df$Label <- ifelse(deg_df$Gene %in% top10_genes, deg_df$Gene, "")

  deg_df$Regulation <- "Not Significant"
  deg_df$Regulation[deg_df$Gene %in% sig_up] <- "Up-regulated"
  deg_df$Regulation[deg_df$Gene %in% sig_down] <- "Down-regulated"

  up_count <- length(sig_up)
  down_count <- length(sig_down)

  p <- ggplot(deg_df, aes(x = logFC_plot, y = neg_log10_p, color = Regulation)) +
    geom_point(size = 1.5, alpha = 0.6) +
    scale_color_manual(
      values = c("Up-regulated" = "#E64B35", "Down-regulated" = "#4DBBD5", "Not Significant" = "gray70"),
      labels = c(
        "Down-regulated" = paste0("Down-regulated (Blue) (n=", down_count, ")"),
        "Not Significant" = "Not Significant (Gray)",
        "Up-regulated"   = paste0("Up-regulated (Red) (n=", up_count, ")")
      )
    ) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray50", alpha = 0.8) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    geom_text_repel(aes(label = Label), size = 3, max.overlaps = 15, box.padding = 0.5) +
    theme_minimal() +
    labs(title = paste0(title, "\nThreshold: ", threshold_text),
         x = "log2(Fold Change)", y = "-log10(p-value)", color = "Regulation") +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          legend.position = "bottom")

  ggsave(file.path(fig.path, filename), p, width = 8, height = 7)
  png_filename <- gsub("\\.pdf$", ".png", filename)
  ggsave(file.path(fig.path, png_filename), p, width = 8, height = 7, dpi = 300)
  cat("  ✅ 已保存:", filename, "和", png_filename, "\n")
  return(p)
}

plot_deg_heatmap <- function(expr_matrix, group_info, sig_genes, title, filename,
                             top_n = 10) {

  p_col  <- if ("P.Value"        %in% colnames(sig_genes)) "P.Value"
             else if ("pvalue"   %in% colnames(sig_genes)) "pvalue"
             else stop("sig_genes 中找不到 P.Value 或 pvalue 列")

  fc_col <- if ("logFC"              %in% colnames(sig_genes)) "logFC"
             else if ("log2FoldChange" %in% colnames(sig_genes)) "log2FoldChange"
             else stop("sig_genes 中找不到 logFC 或 log2FoldChange 列")

  up_genes   <- sig_genes[sig_genes[[fc_col]] > 0, ]
  down_genes <- sig_genes[sig_genes[[fc_col]] < 0, ]

  up_genes   <- up_genes[order(up_genes[[p_col]]),   ][seq_len(min(top_n, nrow(up_genes))),   ]
  down_genes <- down_genes[order(down_genes[[p_col]]),][seq_len(min(top_n, nrow(down_genes))),]

  selected <- c(up_genes$Gene, down_genes$Gene)
  genes_plot <- intersect(selected, rownames(expr_matrix))

  cat(sprintf("  热图基因: 上调 %d 个, 下调 %d 个（共 %d 个）\n",
              sum(genes_plot %in% up_genes$Gene),
              sum(genes_plot %in% down_genes$Gene),
              length(genes_plot)))

  if (length(genes_plot) < 2) {
    cat("  ⚠️ 可用基因少于2个，跳过热图\n")
    return(NULL)
  }

  expr_plot <- expr_matrix[genes_plot, , drop = FALSE]

  if ("sample" %in% colnames(group_info)) {
    sample_order <- group_info %>% arrange(group) %>% pull(sample)
  } else {
    sample_order <- rownames(group_info)
  }
  sample_order <- intersect(sample_order, colnames(expr_plot))
  if (length(sample_order) == 0) {
    cat("  ⚠️ 没有匹配的样本，跳过热图\n")
    return(NULL)
  }
  expr_plot <- expr_plot[, sample_order, drop = FALSE]

  if ("sample" %in% colnames(group_info)) {
    annotation_col <- data.frame(
      Group = factor(group_info$group, levels = unique(group_info$group)),
      row.names = group_info$sample
    )
  } else {
    annotation_col <- data.frame(
      Group = factor(group_info$group, levels = unique(group_info$group)),
      row.names = rownames(group_info)
    )
  }
  annotation_col <- annotation_col[sample_order, , drop = FALSE]

  group_levels <- unique(group_info$group)
  ann_colors <- NULL
  if (length(group_levels) == 2) {
    if ("control" %in% group_levels && "stroke" %in% group_levels) {
      ann_colors <- list(Group = c(control = "#4DBBD5", stroke = "#E64B35"))
    } else if ("control" %in% group_levels && "MACO" %in% group_levels) {
      ann_colors <- list(Group = c(control = "#4DBBD5", MACO = "#E64B35"))
    }
  }

  row_anno <- data.frame(
    Direction = factor(
      ifelse(genes_plot %in% up_genes$Gene, "Up", "Down"),
      levels = c("Up", "Down")
    ),
    row.names = genes_plot
  )
  ann_colors$Direction <- c(Up = "#E64B35", Down = "#4DBBD5")

  expr_scaled <- t(scale(t(expr_plot)))
  expr_scaled[is.na(expr_scaled)] <- 0
  expr_scaled <- pmin(pmax(expr_scaled, -2), 2)

  expr_scaled <- expr_scaled[genes_plot, , drop = FALSE]

  n_genes   <- nrow(expr_scaled)
  n_samples <- ncol(expr_scaled)
  plot_width  <- max(7, n_samples * 0.35 + 3)
  plot_height <- max(5, n_genes  * 0.38 + 2)

  heatmap_args <- list(
    mat               = expr_scaled,
    annotation_col    = annotation_col,
    annotation_row    = row_anno,
    annotation_colors = ann_colors,
    show_rownames     = TRUE,
    show_colnames     = FALSE,
    fontsize_row      = 10,
    fontsize          = 10,
    color             = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
    breaks            = seq(-2, 2, length.out = 101),
    main              = paste0(title, "  (Top ", top_n, " Up & Down by p-value)"),
    cluster_rows      = FALSE,
    cluster_cols      = FALSE,
    gaps_row          = sum(genes_plot %in% up_genes$Gene),
    border_color      = "white",
    cellwidth         = 18,
    cellheight        = 18,
    treeheight_row    = 0,
    treeheight_col    = 0,
    legend_breaks     = c(-2, -1, 0, 1, 2),
    legend_labels     = c("-2", "-1", "0", "1", "2")
  )

  do.call(pdf, list(file = file.path(fig.path, filename),
                    width = plot_width, height = plot_height))
  do.call(pheatmap, heatmap_args)
  dev.off()

  png_filename <- gsub("\\.pdf$", ".png", filename)
  do.call(png, list(filename = file.path(fig.path, png_filename),
                    width = plot_width, height = plot_height,
                    units = "in", res = 300))
  do.call(pheatmap, heatmap_args)
  dev.off()

  cat("  ✅ 已保存:", filename, "和", png_filename, "\n")
}

run_deseq2_analysis <- function(expr_matrix, group_info, group_col, dataset_name) {

  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("分析:", dataset_name, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")

  common_samples <- intersect(colnames(expr_matrix), group_info$sample)
  cat("共同样本数:", length(common_samples), "\n")

  if (length(common_samples) == 0) {
    cat("❌ 没有共同样本\n")
    return(NULL)
  }

  expr_clean <- expr_matrix[, common_samples, drop = FALSE]
  group_clean <- group_info[group_info$sample %in% common_samples, , drop = FALSE]
  group_clean <- group_clean[match(colnames(expr_clean), group_clean$sample), , drop = FALSE]

  cat("清理后样本数:", ncol(expr_clean), "\n")
  print(table(group_clean[[group_col]]))

  expr_round <- round(expr_clean)
  group_clean$group <- factor(group_clean[[group_col]], levels = c("control", "MACO"))

  dds <- DESeqDataSetFromMatrix(
    countData = expr_round,
    colData = data.frame(group = group_clean$group, row.names = colnames(expr_round)),
    design = ~ group
  )

  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep, ]
  cat("过滤后基因数:", nrow(dds), "\n")

  dds <- DESeq(dds)
  res <- results(dds, contrast = c("group", "MACO", "control"))

  deg_df <- as.data.frame(res)
  deg_df$Gene <- rownames(deg_df)
  deg_df <- deg_df[!is.na(deg_df$pvalue), ]

  deg_sig <- deg_df[deg_df$pvalue < 0.05 & abs(deg_df$log2FoldChange) > 0.5, ]

  cat("\n✅ 筛选结果:\n")
  cat("  显著基因数:", nrow(deg_sig), "\n")
  cat("  上调基因:", sum(deg_sig$log2FoldChange > 0, na.rm = TRUE), "\n")
  cat("  下调基因:", sum(deg_sig$log2FoldChange < 0, na.rm = TRUE), "\n")

  return(list(
    all = deg_df,
    sig = deg_sig,
    expr = expr_clean,
    group = group_clean
  ))
}

cat("\n========== 1. 读取数据 ==========\n")

train_file <- "./00_rawdata/GSE16561_for_DEG.rds"
if (file.exists(train_file)) {
  train_data <- readRDS(train_file)
  Train_expr <- train_data$expression
  Train_group <- train_data$group

  cat("✅ GSE16561训练集读取成功\n")
  cat("  表达矩阵维度:", dim(Train_expr), "\n")
  cat("  分组统计:\n")
  print(table(Train_group$group))

  if (any(grepl("GSM", Train_group[[1]]))) {
    rownames(Train_group) <- Train_group[[1]]
    cat("  ✅ 使用第一列作为行名\n")
  } else if ("sample" %in% colnames(Train_group) && any(grepl("GSM", Train_group$sample))) {
    rownames(Train_group) <- Train_group$sample
    cat("  ✅ 使用sample列作为行名\n")
  }

  Train_group$group <- ifelse(tolower(Train_group$group) == "control", "control", "stroke")

} else {
  stop("❌ 训练集文件不存在:", train_file)
}

sample_info_file <- "./00_rawdata/mouse_sample_info.csv"
if (file.exists(sample_info_file)) {
  mouse_sample_info <- read.csv(sample_info_file, stringsAsFactors = FALSE)
  cat("\n✅ 小鼠样本信息读取成功\n")
  cat("  样本数:", nrow(mouse_sample_info), "\n")
} else {
  stop("❌ 小鼠样本信息文件不存在:", sample_info_file)
}

brain_file <- "./00_rawdata/brain_for_DEG.rds"
if (file.exists(brain_file)) {
  brain_data <- readRDS(brain_file)
  Brain_expr <- brain_data$expression
  cat("\n✅ 脑组织数据读取成功\n")
  cat("  表达矩阵维度:", dim(Brain_expr), "\n")
} else {
  stop("❌ 脑组织文件不存在:", brain_file)
}

spleen_file <- "./00_rawdata/spleen_for_DEG.rds"
if (file.exists(spleen_file)) {
  spleen_data <- readRDS(spleen_file)
  Spleen_expr <- spleen_data$expression
  cat("\n✅ 脾脏组织数据读取成功\n")
  cat("  表达矩阵维度:", dim(Spleen_expr), "\n")
} else {
  stop("❌ 脾脏组织文件不存在:", spleen_file)
}

cat("\n========== 2. 处理小鼠数据分组 ==========\n")

brain_samples <- mouse_sample_info[mouse_sample_info$tissue == "brain", ]
cat("\n脑组织样本数:", nrow(brain_samples), "\n")
print(table(brain_samples$strain, brain_samples$group))

spleen_samples <- mouse_sample_info[mouse_sample_info$tissue == "spleen", ]
cat("\n脾脏组织样本数:", nrow(spleen_samples), "\n")
print(table(spleen_samples$strain, spleen_samples$group))

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("3. DEG1: GSE16561差异分析 (limma)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("使用阈值: P < 0.05 & |logFC| > 0.5\n")

common_samples <- intersect(rownames(Train_group), colnames(Train_expr))
cat("共同样本数:", length(common_samples), "\n")

if (length(common_samples) == 0) {
  cat("❌ 没有共同样本，无法进行差异分析\n")
  cat("  请检查 Train_group 的行名是否正确\n")
  stop("样本名不匹配")
}

Train_expr <- Train_expr[, common_samples, drop = FALSE]
Train_group <- Train_group[common_samples, , drop = FALSE]

Train_group$group <- factor(Train_group$group, levels = c("control", "stroke"))
design <- model.matrix(~ group, data = Train_group)

fit <- lmFit(Train_expr, design)
fit <- eBayes(fit)

DEGs1 <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
DEGs1$Gene <- rownames(DEGs1)

DEG1_sig <- DEGs1[DEGs1$P.Value < 0.05 & abs(DEGs1$logFC) > 0.5, ]
threshold_used1 <- "P<0.05 & |logFC|>0.5"

cat("\n✅ DEG1筛选结果:\n")
cat("  显著基因数:", nrow(DEG1_sig), "\n")
cat("  上调基因:", sum(DEG1_sig$logFC > 0, na.rm = TRUE), "\n")
cat("  下调基因:", sum(DEG1_sig$logFC < 0, na.rm = TRUE), "\n")

write.csv(DEGs1, file.path(res.path, "DEG1_GSE16561_all.csv"), row.names = FALSE)
write.csv(DEG1_sig, file.path(res.path, "DEG1_GSE16561_significant.csv"), row.names = FALSE)
writeLines(DEG1_sig$Gene, file.path(res.path, "DEG1_GSE16561_genes.txt"))

cat("\n绘制 Fig.1a: GSE16561火山图\n")
p1a <- plot_volcano_limma(DEGs1, DEG1_sig, "DEG1: GSE16561 (Training Set)",
                          "Fig1a_DEG1_Volcano_final.pdf", threshold_used1)

cat("\n绘制 Fig.1a热图\n")
plot_deg_heatmap(Train_expr, Train_group, DEG1_sig,
                 "DEG1: GSE16561 - DEGs", "Fig1a_DEG1_Heatmap_final.pdf")

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("4. DEG2: 脑组织 C57品系 (control vs MACO)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

c57_brain_samples <- brain_samples[brain_samples$strain == "C57", ]
c57_brain_samples <- c57_brain_samples[c57_brain_samples$group %in% c("control", "MACO"), ]

cat("C57脑组织样本数:", nrow(c57_brain_samples), "\n")
print(table(c57_brain_samples$group))

if (nrow(c57_brain_samples) >= 4) {
  c57_brain_expr <- Brain_expr[, c57_brain_samples$sample_id, drop = FALSE]

  c57_brain_group <- data.frame(
    sample = c57_brain_samples$sample_id,
    group = c57_brain_samples$group,
    stringsAsFactors = FALSE
  )

  result2 <- run_deseq2_analysis(c57_brain_expr, c57_brain_group, "group", "DEG2: Brain_C57")

  if (!is.null(result2)) {
    DEG2_all <- result2$all
    DEG2_sig <- result2$sig
    Brain_expr_clean <- result2$expr
    Brain_group_clean <- result2$group

    write.csv(DEG2_all, file.path(res.path, "DEG2_Brain_C57_all.csv"), row.names = FALSE)
    write.csv(DEG2_sig, file.path(res.path, "DEG2_Brain_C57_significant.csv"), row.names = FALSE)
    writeLines(DEG2_sig$Gene, file.path(res.path, "DEG2_Brain_C57_genes.txt"))

    cat("\n绘制 Fig.1b: 脑组织C57火山图\n")
    p1b <- plot_volcano_deseq2(DEG2_all, DEG2_sig, "DEG2: Brain Tissue (C57, MACO vs Control)",
                                "Fig1b_DEG2_Brain_C57_Volcano_final.pdf", "p<0.05 & |logFC|>0.5")

    cat("\n绘制 Fig.1b热图\n")
    plot_deg_heatmap(Brain_expr_clean, Brain_group_clean, DEG2_sig,
                     "DEG2: Brain C57 - DEGs", "Fig1b_DEG2_Brain_C57_Heatmap_final.pdf")
  }
} else {
  cat("⚠️ C57脑组织样本不足，跳过分析\n")
}

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("5. DEG3: 脑组织 IL-21品系 (control vs MACO)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

il21_brain_samples <- brain_samples[brain_samples$strain == "IL-21", ]
il21_brain_samples <- il21_brain_samples[il21_brain_samples$group %in% c("control", "MACO"), ]

cat("IL-21脑组织样本数:", nrow(il21_brain_samples), "\n")
print(table(il21_brain_samples$group))

if (nrow(il21_brain_samples) >= 4) {
  il21_brain_expr <- Brain_expr[, il21_brain_samples$sample_id, drop = FALSE]

  il21_brain_group <- data.frame(
    sample = il21_brain_samples$sample_id,
    group = il21_brain_samples$group,
    stringsAsFactors = FALSE
  )

  result3 <- run_deseq2_analysis(il21_brain_expr, il21_brain_group, "group", "DEG3: Brain_IL21")

  if (!is.null(result3)) {
    DEG3_all <- result3$all
    DEG3_sig <- result3$sig
    Brain_expr_clean <- result3$expr
    Brain_group_clean <- result3$group

    write.csv(DEG3_all, file.path(res.path, "DEG3_Brain_IL21_all.csv"), row.names = FALSE)
    write.csv(DEG3_sig, file.path(res.path, "DEG3_Brain_IL21_significant.csv"), row.names = FALSE)
    writeLines(DEG3_sig$Gene, file.path(res.path, "DEG3_Brain_IL21_genes.txt"))

    cat("\n绘制 Fig.1c: 脑组织IL-21火山图\n")
    p1c <- plot_volcano_deseq2(DEG3_all, DEG3_sig, "DEG3: Brain Tissue (IL-21, MACO vs Control)",
                                "Fig1c_DEG3_Brain_IL21_Volcano_final.pdf", "p<0.05 & |logFC|>0.5")

    cat("\n绘制 Fig.1c热图\n")
    plot_deg_heatmap(Brain_expr_clean, Brain_group_clean, DEG3_sig,
                     "DEG3: Brain IL-21 - DEGs", "Fig1c_DEG3_Brain_IL21_Heatmap_final.pdf")
  }
} else {
  cat("⚠️ IL-21脑组织样本不足，跳过分析\n")
}

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("6. DEG4: 脾脏组织 C57品系 (control vs MACO)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

c57_spleen_samples <- spleen_samples[spleen_samples$strain == "C57", ]
c57_spleen_samples <- c57_spleen_samples[c57_spleen_samples$group %in% c("control", "MACO"), ]

cat("C57脾脏组织样本数:", nrow(c57_spleen_samples), "\n")
print(table(c57_spleen_samples$group))

if (nrow(c57_spleen_samples) >= 4) {
  c57_spleen_expr <- Spleen_expr[, c57_spleen_samples$sample_id, drop = FALSE]

  c57_spleen_group <- data.frame(
    sample = c57_spleen_samples$sample_id,
    group = c57_spleen_samples$group,
    stringsAsFactors = FALSE
  )

  result4 <- run_deseq2_analysis(c57_spleen_expr, c57_spleen_group, "group", "DEG4: Spleen_C57")

  if (!is.null(result4)) {
    DEG4_all <- result4$all
    DEG4_sig <- result4$sig
    Spleen_expr_clean <- result4$expr
    Spleen_group_clean <- result4$group

    write.csv(DEG4_all, file.path(res.path, "DEG4_Spleen_C57_all.csv"), row.names = FALSE)
    write.csv(DEG4_sig, file.path(res.path, "DEG4_Spleen_C57_significant.csv"), row.names = FALSE)
    writeLines(DEG4_sig$Gene, file.path(res.path, "DEG4_Spleen_C57_genes.txt"))

    cat("\n绘制 Fig.1d: 脾脏组织C57火山图\n")
    p1d <- plot_volcano_deseq2(DEG4_all, DEG4_sig, "DEG4: Spleen Tissue (C57, MACO vs Control)",
                                "Fig1d_DEG4_Spleen_C57_Volcano_final.pdf", "p<0.05 & |logFC|>0.5")

    cat("\n绘制 Fig.1d热图\n")
    plot_deg_heatmap(Spleen_expr_clean, Spleen_group_clean, DEG4_sig,
                     "DEG4: Spleen C57 - DEGs", "Fig1d_DEG4_Spleen_C57_Heatmap_final.pdf")
  }
} else {
  cat("⚠️ C57脾脏组织样本不足，跳过分析\n")
}

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("7. DEG5: 脾脏组织 IL-21品系 (control vs MACO)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

il21_spleen_samples <- spleen_samples[spleen_samples$strain == "IL-21", ]
il21_spleen_samples <- il21_spleen_samples[il21_spleen_samples$group %in% c("control", "MACO"), ]

cat("IL-21脾脏组织样本数:", nrow(il21_spleen_samples), "\n")
print(table(il21_spleen_samples$group))

if (nrow(il21_spleen_samples) >= 4) {
  il21_spleen_expr <- Spleen_expr[, il21_spleen_samples$sample_id, drop = FALSE]

  il21_spleen_group <- data.frame(
    sample = il21_spleen_samples$sample_id,
    group = il21_spleen_samples$group,
    stringsAsFactors = FALSE
  )

  result5 <- run_deseq2_analysis(il21_spleen_expr, il21_spleen_group, "group", "DEG5: Spleen_IL21")

  if (!is.null(result5)) {
    DEG5_all <- result5$all
    DEG5_sig <- result5$sig
    Spleen_expr_clean <- result5$expr
    Spleen_group_clean <- result5$group

    write.csv(DEG5_all, file.path(res.path, "DEG5_Spleen_IL21_all.csv"), row.names = FALSE)
    write.csv(DEG5_sig, file.path(res.path, "DEG5_Spleen_IL21_significant.csv"), row.names = FALSE)
    writeLines(DEG5_sig$Gene, file.path(res.path, "DEG5_Spleen_IL21_genes.txt"))

    cat("\n绘制 Fig.1e: 脾脏组织IL-21火山图\n")
    p1e <- plot_volcano_deseq2(DEG5_all, DEG5_sig, "DEG5: Spleen Tissue (IL-21, MACO vs Control)",
                                "Fig1e_DEG5_Spleen_IL21_Volcano_final.pdf", "p<0.05 & |logFC|>0.5")

    cat("\n绘制 Fig.1e热图\n")
    plot_deg_heatmap(Spleen_expr_clean, Spleen_group_clean, DEG5_sig,
                     "DEG5: Spleen IL-21 - DEGs", "Fig1e_DEG5_Spleen_IL21_Heatmap_final.pdf")
  }
} else {
  cat("⚠️ IL-21脾脏组织样本不足，跳过分析\n")
}

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("8. 保存结果汇总\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

deg_stats <- data.frame(
  Dataset = c(
    "DEG1_GSE16561",
    ifelse(exists("DEG2_sig"), "DEG2_Brain_C57", NA),
    ifelse(exists("DEG3_sig"), "DEG3_Brain_IL21", NA),
    ifelse(exists("DEG4_sig"), "DEG4_Spleen_C57", NA),
    ifelse(exists("DEG5_sig"), "DEG5_Spleen_IL21", NA)
  ),
  Significant_Genes = c(
    nrow(DEG1_sig),
    ifelse(exists("DEG2_sig"), nrow(DEG2_sig), NA),
    ifelse(exists("DEG3_sig"), nrow(DEG3_sig), NA),
    ifelse(exists("DEG4_sig"), nrow(DEG4_sig), NA),
    ifelse(exists("DEG5_sig"), nrow(DEG5_sig), NA)
  ),
  Upregulated = c(
    sum(DEG1_sig$logFC > 0, na.rm = TRUE),
    ifelse(exists("DEG2_sig"), sum(DEG2_sig$log2FoldChange > 0, na.rm = TRUE), NA),
    ifelse(exists("DEG3_sig"), sum(DEG3_sig$log2FoldChange > 0, na.rm = TRUE), NA),
    ifelse(exists("DEG4_sig"), sum(DEG4_sig$log2FoldChange > 0, na.rm = TRUE), NA),
    ifelse(exists("DEG5_sig"), sum(DEG5_sig$log2FoldChange > 0, na.rm = TRUE), NA)
  ),
  Downregulated = c(
    sum(DEG1_sig$logFC < 0, na.rm = TRUE),
    ifelse(exists("DEG2_sig"), sum(DEG2_sig$log2FoldChange < 0, na.rm = TRUE), NA),
    ifelse(exists("DEG3_sig"), sum(DEG3_sig$log2FoldChange < 0, na.rm = TRUE), NA),
    ifelse(exists("DEG4_sig"), sum(DEG4_sig$log2FoldChange < 0, na.rm = TRUE), NA),
    ifelse(exists("DEG5_sig"), sum(DEG5_sig$log2FoldChange < 0, na.rm = TRUE), NA)
  ),
  Threshold = rep("P<0.05 & |logFC|>0.5", 5)
)

deg_stats <- deg_stats[!is.na(deg_stats$Dataset), ]

write.csv(deg_stats, file.path(res.path, "DEGs_summary.csv"), row.names = FALSE)

cat("\n📊 差异分析结果汇总:\n")
print(deg_stats)

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("             差异表达分析完成！\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\n📁 输出目录:", output_dir, "\n")
cat("  Results/ - 差异分析结果\n")
cat("  Figures/ - 火山图和热图\n")
cat("\n🎯 下一步分析:\n")
cat("  1. 运行 02_CandidateGene_Intersect 获取候选基因\n")
cat("  2. 运行 03_GO_KEGG_PPI 进行富集分析和PPI网络构建\n")

cat("\n✅ 所有处理完成！\n")
