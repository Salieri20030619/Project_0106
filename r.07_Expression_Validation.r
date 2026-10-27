
rm(list = ls())

options(bitmapType = "cairo")

setwd("")

output_dir <- "./07_Expression_Validation"
res.path <- file.path(output_dir, "Results")
fig.path <- file.path(output_dir, "Figures")

for (path in c(output_dir, res.path, fig.path)) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

cat("\n========== 加载 R 包 ==========\n")

packages <- c("ggplot2", "pheatmap", "RColorBrewer", "gridExtra",
              "ggpubr", "caret", "pROC", "dplyr", "tidyr", "mboost",
              "glmnet", "randomForestSRC", "gbm", "e1071", "MASS", "xgboost")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org", quiet = TRUE)
    library(pkg, character.only = TRUE)
  }
}

select <- dplyr::select
filter <- dplyr::filter
arrange <- dplyr::arrange
mutate <- dplyr::mutate

cat("✅ R 包加载完成\n")

cat("\n========== Step 1: 加载第五部分结果 ==========\n")

ml_results_path <- "./05_Machine_Learning/Results"

required_files <- c(
  file.path(ml_results_path, "best_model_final.rds"),
  file.path(ml_results_path, "best_features_final.txt"),
  file.path(ml_results_path, "all_results_final.rds")
)
for (f in required_files) {
  if (!file.exists(f)) stop("❌ 缺失文件: ", f)
}

best_features <- readLines(file.path(ml_results_path, "best_features_final.txt"))
best_features <- best_features[nchar(best_features) > 0]

all_results <- readRDS(file.path(ml_results_path, "all_results_final.rds"))
best_model_name <- all_results$best_model
best_validation_auc <- all_results$best_validation_auc
best_cv_auc <- all_results$best_cv_auc
best_model_fit <- readRDS(file.path(ml_results_path, "best_model_final.rds"))

cat("最佳模型:", best_model_name, "\n")
cat("交叉验证 AUC:", round(best_cv_auc, 4), "\n")
cat("验证集 AUC:", round(best_validation_auc, 4), "\n")
cat("特征基因数:", length(best_features), "\n")
cat("特征基因:", paste(best_features, collapse = ", "), "\n")

cat("\n========== Step 2: 加载已处理数据 ==========\n")

processed_data_file <- "./05_Machine_Learning/Intermediate/processed_data.rds"
if (!file.exists(processed_data_file)) {
  stop("❌ processed_data.rds 不存在")
}

processed <- readRDS(processed_data_file)
X_train <- processed$Train_set
X_test  <- processed$Test_set
y_train <- processed$Train_class$group
y_test  <- processed$Test_class$group
common_features <- colnames(X_train)

cat("训练集:", dim(X_train), "(0/1:", sum(y_train == 0), "/", sum(y_train == 1), ")\n")
cat("验证集:", dim(X_test),  "(0/1:", sum(y_test  == 0), "/", sum(y_test  == 1), ")\n")

cat("\n========== Step 3: 基于 SHAP 重要性取前 5 基因 ==========\n")

shap_file <- "./06_SHAP_Analysis/Results/shap_importance.csv"

if (!file.exists(shap_file)) {
  stop("❌ SHAP 重要性文件不存在: ", shap_file,
       "\n请先运行第六部分 SHAP 分析")
}

shap_df <- read.csv(shap_file, stringsAsFactors = FALSE)
cat("SHAP 重要性表 (来自第六部分):\n")
print(shap_df)

shap_col_feature <- intersect(c("Feature", "feature", "Gene", "gene"), colnames(shap_df))[1]
shap_col_value <- intersect(c("Mean_SHAP", "mean_shap", "Mean_abs_SHAP",
                               "mean_abs_shap", "SHAP", "shap"),
                             colnames(shap_df))[1]

if (is.na(shap_col_feature) || is.na(shap_col_value)) {
  stop("❌ SHAP 表列名无法识别, 应包含 Feature 和 Mean_SHAP")
}

shap_df <- shap_df[order(shap_df[[shap_col_value]], decreasing = TRUE), ]

shap_df_filtered <- shap_df[toupper(shap_df[[shap_col_feature]]) %in%
                              toupper(common_features), ]

cat(sprintf("\nSHAP 表中匹配的特征基因: %d / %d\n",
            nrow(shap_df_filtered), nrow(shap_df)))

top_n <- min(5, nrow(shap_df_filtered))
top5_genes_raw <- shap_df_filtered[[shap_col_feature]][seq_len(top_n)]

top5_genes <- common_features[toupper(common_features) %in% toupper(top5_genes_raw)]

cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
cat("SHAP 重要性 Top", length(top5_genes), "基因 (作为生物标志物候选):\n")
for (i in seq_along(top5_genes)) {
  g <- top5_genes[i]
  shap_val <- shap_df_filtered[[shap_col_value]][
    toupper(shap_df_filtered[[shap_col_feature]]) == toupper(g)][1]
  cat(sprintf("  %d. %-15s (SHAP = %.4f)\n", i, g, shap_val))
}
cat(paste(rep("=", 60), collapse = ""), "\n", sep = "")

candidate_genes <- top5_genes

cat("\n========== Step 4: AUC + 表达量验证 ==========\n")

results_df <- data.frame()

for (gene in candidate_genes) {

  train_control <- X_train[y_train == 0, gene]
  train_stroke  <- X_train[y_train == 1, gene]
  train_p <- wilcox.test(train_control, train_stroke)$p.value
  train_median_control <- median(train_control)
  train_median_stroke  <- median(train_stroke)
  train_trend <- ifelse(train_median_stroke > train_median_control, "上调", "下调")

  test_control <- X_test[y_test == 0, gene]
  test_stroke  <- X_test[y_test == 1, gene]
  test_p <- wilcox.test(test_control, test_stroke)$p.value
  test_median_control <- median(test_control)
  test_median_stroke  <- median(test_stroke)
  test_trend <- ifelse(test_median_stroke > test_median_control, "上调", "下调")

  roc_train_g <- roc(y_train, X_train[, gene], quiet = TRUE)
  train_auc <- round(roc_train_g$auc, 4)

  roc_test_g <- roc(y_test, X_test[, gene], quiet = TRUE)
  test_auc <- round(roc_test_g$auc, 4)

  shap_rank <- which(toupper(shap_df_filtered[[shap_col_feature]]) == toupper(gene))[1]
  shap_value <- shap_df_filtered[[shap_col_value]][shap_rank]

  pass_filter <- (train_auc > 0.7) & (test_auc > 0.7) &
                 (train_p   < 0.05) & (test_p   < 0.05) &
                 (train_trend == test_trend)

  results_df <- rbind(results_df, data.frame(
    Gene = gene,
    SHAP_Rank = shap_rank,
    SHAP_Value = round(shap_value, 4),
    Train_AUC = train_auc,
    Test_AUC = test_auc,
    Train_p = train_p,
    Test_p = test_p,
    Train_Trend = train_trend,
    Test_Trend = test_trend,
    Pass_Filter = pass_filter,
    stringsAsFactors = FALSE
  ))
}

results_df <- results_df %>% arrange(SHAP_Rank)

cat("\nSHAP 前 5 基因的验证结果:\n")
print(results_df)

selected_genes <- results_df$Gene[results_df$Pass_Filter]

cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("生物标志物筛选结果\n")
cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("候选: SHAP Top 5 基因\n")
cat("筛选标准: 训练/验证集 AUC>0.7 且 p<0.05 且 趋势一致\n\n")

for (i in seq_len(nrow(results_df))) {
  status <- ifelse(results_df$Pass_Filter[i], "✅ 通过", "❌ 未通过")
  cat(sprintf("  [%d] %-12s | SHAP=%.4f | TrainAUC=%.3f TestAUC=%.3f | TrainP=%.2e TestP=%.2e | %s/%s | %s\n",
              i, results_df$Gene[i], results_df$SHAP_Value[i],
              results_df$Train_AUC[i], results_df$Test_AUC[i],
              results_df$Train_p[i],   results_df$Test_p[i],
              results_df$Train_Trend[i], results_df$Test_Trend[i], status))
}

cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat("✅ 最终生物标志物: ", length(selected_genes), " 个\n", sep = "")
if (length(selected_genes) > 0) {
  cat("   ", paste(selected_genes, collapse = ", "), "\n", sep = "")
} else {
  cat("   (无基因通过筛选, 请检查数据或调整阈值)\n")
}
cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")

if (length(selected_genes) == 0) {
  warning("❌ SHAP Top 5 基因全部未通过筛选!")
  cat("\n警告: 0 个基因通过筛选, 后续单细胞分析将无法进行\n")
  cat("建议:\n")
  cat("  1. 检查训练集/验证集的批次效应是否过大\n")
  cat("  2. 考虑放宽筛选阈值 (如 AUC > 0.6)\n")
  cat("  3. 确认 SHAP 前 5 基因在两数据集中都可量化\n")
}

if (length(selected_genes) >= 1 && length(selected_genes) < 2) {
  cat("\n⚠️ 仅 1 个生物标志物通过筛选, 后续分析功效有限\n")
}

cat("\n========== Step 5: 加载预测函数 ==========\n")

ml_r_file <- "ML.R"
if (!file.exists(ml_r_file)) {
  stop("❌ ML.R 不存在")
}
source(ml_r_file)
cat("✅ ML.R 加载\n")

predict_function <- function(model, newdata) {
  tryCatch({
    if (is.data.frame(newdata)) newdata <- as.matrix(newdata)
    newdata[is.na(newdata)] <- 0
    newdata <- as.data.frame(newdata)
    CalPredictScore(fit = model, new_data = newdata)
  }, error = function(e) {
    warning(sprintf("predict_function 失败 [%s]: %s",
                    paste(class(model), collapse = "/"), e$message))
    rep(NA_real_, nrow(newdata))
  })
}

cat("\n========== Step 6: 绘制 ROC 曲线 ==========\n")

train_pred <- predict_function(best_model_fit, X_train)
test_pred  <- predict_function(best_model_fit, X_test)

if (any(is.na(train_pred)) || any(is.na(test_pred))) {
  stop("预测返回 NA, 请检查 ML.R")
}

roc_train <- roc(y_train, train_pred, quiet = TRUE)
roc_test  <- roc(y_test,  test_pred,  quiet = TRUE)

train_auc <- round(roc_train$auc, 4)
test_auc  <- round(roc_test$auc,  4)

draw_roc_main <- function() {
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "1 - Specificity", ylab = "Sensitivity",
       main = paste("ROC Curves -", best_model_name), cex.main = 1.2)
  abline(0, 1, col = "gray", lty = 2, lwd = 2)
  lines(1 - roc_train$specificities, roc_train$sensitivities,
        col = "#4DBBD5", lwd = 3)
  lines(1 - roc_test$specificities, roc_test$sensitivities,
        col = "#E64B35", lwd = 3)
  legend("bottomright",
         legend = c(paste("Training AUC =", train_auc),
                    paste("Validation AUC =", test_auc)),
         col = c("#4DBBD5", "#E64B35"), lwd = 3, bty = "n")
}

pdf(file.path(fig.path, "ROC_Curves.pdf"), width = 8, height = 7)
draw_roc_main()
dev.off()

png(file.path(fig.path, "ROC_Curves.png"),
    width = 8, height = 7, units = "in", res = 300)
draw_roc_main()
dev.off()

cat("✅ ROC 曲线已保存\n")

cat("\n========== Step 7: 表达水平箱线图 ==========\n")

plot_genes <- candidate_genes

train_expr_df <- as.data.frame(X_train[, plot_genes, drop = FALSE])
train_expr_df$Group <- ifelse(y_train == 1, "Stroke", "Control")
train_expr_df$Dataset <- "Training"

test_expr_df <- as.data.frame(X_test[, plot_genes, drop = FALSE])
test_expr_df$Group <- ifelse(y_test == 1, "Stroke", "Control")
test_expr_df$Dataset <- "Validation"

all_expr_df <- rbind(train_expr_df, test_expr_df)

plot_long <- all_expr_df %>%
  pivot_longer(cols = all_of(plot_genes),
               names_to = "Gene", values_to = "Expression")

plot_long <- plot_long %>%
  left_join(results_df[, c("Gene", "Train_p", "Test_p", "Pass_Filter")],
            by = "Gene") %>%
  mutate(Signif = case_when(
    Dataset == "Training"   & Train_p < 0.001 ~ "***",
    Dataset == "Training"   & Train_p < 0.01  ~ "**",
    Dataset == "Training"   & Train_p < 0.05  ~ "*",
    Dataset == "Validation" & Test_p  < 0.001 ~ "***",
    Dataset == "Validation" & Test_p  < 0.01  ~ "**",
    Dataset == "Validation" & Test_p  < 0.05  ~ "*",
    TRUE ~ "ns"
  ))

plot_long$Gene <- factor(plot_long$Gene, levels = plot_genes)

max_vals <- plot_long %>%
  group_by(Gene, Dataset) %>%
  summarise(max_exp = max(Expression, na.rm = TRUE), .groups = "drop")

signif_labels <- plot_long %>%
  select(Gene, Dataset, Signif) %>%
  distinct()

p_box <- ggplot(plot_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1.0) +
  scale_fill_manual(values = c(Control = "#45a9b8", Stroke = "#f76a56")) +
  facet_grid(Gene ~ Dataset, scales = "free_y") +
  labs(title = "SHAP Top 5 Biomarkers - Expression Validation",
       subtitle = "Wilcoxon test: * p<0.05, ** p<0.01, *** p<0.001",
       x = "", y = "Normalized Expression") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9),
        strip.background = element_rect(fill = "lightgray"),
        strip.text = element_text(size = 11, face = "bold"),
        legend.position = "top")

p_box_final <- p_box +
  geom_text(data = merge(max_vals, signif_labels,
                          by = c("Gene", "Dataset")),
            aes(x = 1.5, y = max_exp + 0.15, label = Signif),
            inherit.aes = FALSE, size = 5, color = "black")

n_genes <- length(plot_genes)
plot_height <- min(max(6, n_genes * 2.5), 49)

ggsave(file.path(fig.path, "Expression_Boxplot.pdf"),
       p_box_final, width = 10, height = plot_height, limitsize = FALSE)
ggsave(file.path(fig.path, "Expression_Boxplot.png"),
       p_box_final, width = 10, height = plot_height, dpi = 300, limitsize = FALSE)

cat("✅ 表达水平箱线图已保存\n")

cat("\n========== Step 8: 单基因 ROC 曲线 ==========\n")

n_roc_genes <- length(candidate_genes)
n_roc_cols  <- min(2, n_roc_genes)
n_roc_rows  <- ceiling(n_roc_genes / n_roc_cols)
roc_width   <- n_roc_cols * 5
roc_height  <- min(n_roc_rows * 4, 49)

draw_single_roc <- function() {
  par(mfrow = c(n_roc_rows, n_roc_cols), mar = c(4, 4, 3, 2))
  for (gene in candidate_genes) {
    rt <- roc(y_train, X_train[, gene], quiet = TRUE)
    rv <- roc(y_test,  X_test[, gene],  quiet = TRUE)

    is_pass <- results_df$Pass_Filter[results_df$Gene == gene]
    title_str <- paste0("ROC - ", gene, ifelse(is_pass, " ✅", " ✗"))

    plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
         xlab = "1 - Specificity", ylab = "Sensitivity", main = title_str)
    abline(0, 1, col = "gray", lty = 2)
    lines(1 - rt$specificities, rt$sensitivities, col = "#4DBBD5", lwd = 2)
    lines(1 - rv$specificities, rv$sensitivities, col = "#E64B35", lwd = 2)
    legend("bottomright",
           legend = c(paste("Train AUC =", round(rt$auc, 3)),
                      paste("Test AUC =",  round(rv$auc,  3))),
           col = c("#4DBBD5", "#E64B35"), lwd = 2, bty = "n", cex = 0.8)
  }
}

pdf(file.path(fig.path, "Single_Gene_ROC_Curves.pdf"),
    width = roc_width, height = roc_height)
draw_single_roc()
dev.off()

png(file.path(fig.path, "Single_Gene_ROC_Curves.png"),
    width = roc_width, height = roc_height, units = "in", res = 150)
draw_single_roc()
dev.off()

cat("✅ 单基因 ROC 已保存\n")

cat("\n========== Step 9: 特征相关性热图 ==========\n")

combined_expr <- rbind(X_train[, candidate_genes, drop = FALSE],
                       X_test[, candidate_genes, drop = FALSE])
cor_matrix <- cor(combined_expr, method = "spearman")
cor_matrix[is.na(cor_matrix)] <- 0

col_fun <- colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100)

pdf(file.path(fig.path, "Feature_Correlation_Heatmap.pdf"), width = 6, height = 5)
pheatmap(cor_matrix,
         color = col_fun,
         main = "Spearman correlation among SHAP Top 5",
         display_numbers = TRUE,
         number_format = "%.2f",
         fontsize_number = 11,
         border_color = "grey")
dev.off()

png(file.path(fig.path, "Feature_Correlation_Heatmap.png"),
    width = 6, height = 5, units = "in", res = 300)
pheatmap(cor_matrix,
         color = col_fun,
         main = "Spearman correlation among SHAP Top 5",
         display_numbers = TRUE,
         number_format = "%.2f",
         fontsize_number = 11,
         border_color = "grey")
dev.off()

cat("✅ 相关性热图已保存\n")

cat("\n========== Step 10: 保存结果 ==========\n")

write.csv(results_df, file.path(res.path, "biomarker_screening_results.csv"),
          row.names = FALSE)

writeLines(selected_genes, file.path(res.path, "selected_biomarkers.txt"))

writeLines(candidate_genes,
           file.path(res.path, "shap_top5_candidates.txt"))

expr_data <- list(
  training   = data.frame(Sample = rownames(X_train), Group = y_train, X_train),
  validation = data.frame(Sample = rownames(X_test),  Group = y_test,  X_test)
)
saveRDS(expr_data, file.path(res.path, "expression_data.rds"))

prediction_data <- data.frame(
  Sample = rownames(X_test),
  Group = y_test,
  Predicted_Score = test_pred
)
write.csv(prediction_data, file.path(res.path, "validation_predictions.csv"),
          row.names = FALSE)

report_lines <- c(
  "# 第七部分 报告",
  paste("生成时间:", Sys.time()),
  "",
  "## 改进说明 (vs 旧版)",
  "- 去除硬编码优先候选基因 (ABCA1/LRG1/CYP1B1/HSPA1A/ADM)",
  "- 直接基于 SHAP 重要性取 Top 5 作为候选生物标志物",
  "- 对 5 个候选做 AUC + 表达量 + 趋势一致性筛选",
  "- 通过几个用几个, 不强制要求 5 个全通过",
  "",
  "## 模型信息",
  paste("- 最佳模型:", best_model_name),
  paste("- 训练集 AUC:", train_auc),
  paste("- 验证集 AUC:", test_auc),
  paste("- 总特征数:", length(common_features)),
  "",
  "## SHAP Top 5 候选基因"
)

for (i in seq_along(candidate_genes)) {
  g <- candidate_genes[i]
  shap_val <- results_df$SHAP_Value[results_df$Gene == g]
  report_lines <- c(report_lines,
    sprintf("- [%d] %s (SHAP = %.4f)", i, g, shap_val))
}

report_lines <- c(report_lines,
  "",
  "## 筛选标准",
  "- 训练集 AUC > 0.7",
  "- 验证集 AUC > 0.7",
  "- 训练集 Wilcoxon p < 0.05",
  "- 验证集 Wilcoxon p < 0.05",
  "- 训练集和验证集表达趋势一致",
  "",
  "## 候选基因详细结果",
  "",
  "| Gene | SHAP | TrainAUC | TestAUC | TrainP | TestP | Trend | Pass |",
  "|---|---|---|---|---|---|---|---|"
)

for (i in seq_len(nrow(results_df))) {
  r <- results_df[i, ]
  report_lines <- c(report_lines,
    sprintf("| %s | %.4f | %.3f | %.3f | %.2e | %.2e | %s/%s | %s |",
            r$Gene, r$SHAP_Value, r$Train_AUC, r$Test_AUC,
            r$Train_p, r$Test_p, r$Train_Trend, r$Test_Trend,
            ifelse(r$Pass_Filter, "✅", "❌")))
}

report_lines <- c(report_lines,
  "",
  paste("## 最终生物标志物: ", length(selected_genes), " 个", sep = ""),
  "",
  if (length(selected_genes) > 0) paste("- ", selected_genes, collapse = "\n") else "- 无 (全部未通过筛选)",
  "",
  "## 输出文件",
  "",
  "### Results/",
  "- biomarker_screening_results.csv  (5 个候选的完整结果)",
  "- selected_biomarkers.txt          (通过筛选的生物标志物)",
  "- shap_top5_candidates.txt         (SHAP 前 5 候选, 包括未通过的)",
  "- expression_data.rds              (训练+验证集表达矩阵)",
  "- validation_predictions.csv       (验证集预测评分)",
  "",
  "### Figures/",
  "- ROC_Curves.pdf/png                  (模型整体 ROC)",
  "- Expression_Boxplot.pdf/png          (5 个候选的表达箱线图)",
  "- Single_Gene_ROC_Curves.pdf/png      (5 个候选的单独 ROC)",
  "- Feature_Correlation_Heatmap.pdf/png (5 个候选间相关性)"
)

writeLines(report_lines, file.path(res.path, "Expression_Validation_Report.md"))

cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
cat("             第七部分 完成\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n", sep = "")

cat("候选: SHAP Top", length(candidate_genes), "基因\n")
for (i in seq_along(candidate_genes)) {
  cat(sprintf("  [%d] %s\n", i, candidate_genes[i]))
}

cat(sprintf("\n通过筛选: %d / %d 个\n",
            length(selected_genes), length(candidate_genes)))
if (length(selected_genes) > 0) {
  cat("最终生物标志物:", paste(selected_genes, collapse = ", "), "\n")
}

cat("\n输出目录:", output_dir, "\n")
cat("\n✅ 第七部分完成\n")

sink(file.path(res.path, "session_info.txt"))
sessionInfo()
sink()
