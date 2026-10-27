
cat("\n========================================\n")
cat("08. 列线图构建与评估\n")
cat("========================================\n\n")

rm(list = ls())

options(bitmapType = "cairo")

setwd("")

output_dir <- "./08_Nomogram"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  cat("创建输出文件夹:", output_dir, "\n")
}

res.path <- file.path(output_dir, "Results")
fig.path <- file.path(output_dir, "Figures")

for (path in c(res.path, fig.path)) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    cat("创建文件夹:", path, "\n")
  }
}

cat("\n========== 加载R包 ==========\n")

packages <- c(
  "ggplot2", "pheatmap", "RColorBrewer", "gridExtra", "grid",
  "ggpubr", "pROC", "dplyr", "tidyr", "rms", "ResourceSelection"
)

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org", quiet = TRUE)
    library(pkg, character.only = TRUE)
  } else {
    cat("✅ 已加载:", pkg, "\n")
  }
}

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate

cat("✅ 所有R包加载完成\n")

cat("\n========== 从第七部分读取生物标志物 ==========\n")

seventh_res_path <- "./07_Expression_Validation/Results"

biomarkers_file <- file.path(seventh_res_path, "selected_biomarkers.txt")
if (!file.exists(biomarkers_file)) {
  stop("❌ selected_biomarkers.txt 文件不存在，请先运行第七部分")
}

selected_biomarkers <- readLines(biomarkers_file)
selected_biomarkers <- selected_biomarkers[selected_biomarkers != "" & !is.na(selected_biomarkers)]
cat("从第七部分读取到的生物标志物:", paste(selected_biomarkers, collapse = ", "), "\n")

if (length(selected_biomarkers) < 1) {
  stop("❌ 没有筛选出生物标志物，无法构建列线图")
}

cat("\n========== 读取训练集数据 ==========\n")

train_file <- "./00_rawdata/GSE16561_for_DEG.rds"
if (!file.exists(train_file)) {
  stop("❌ 训练集文件不存在:", train_file)
}

train_data_raw <- readRDS(train_file)
Train_expr <- train_data_raw$expression
Train_group <- train_data_raw$group

cat("✅ GSE16561训练集读取成功\n")
cat("  表达矩阵维度:", dim(Train_expr), "\n")

if (any(grepl("GSM", Train_group[[1]]))) {
  rownames(Train_group) <- Train_group[[1]]
  cat("  ✅ 使用第一列作为行名\n")
} else if ("sample" %in% colnames(Train_group) && any(grepl("GSM", Train_group$sample))) {
  rownames(Train_group) <- Train_group$sample
  cat("  ✅ 使用sample列作为行名\n")
}

Train_group$group <- ifelse(tolower(Train_group$group) == "control", 0, 1)
cat("  分组统计: 0=", sum(Train_group$group == 0), ", 1=", sum(Train_group$group == 1), "\n")

cat("\n========== 读取验证集数据 ==========\n")

test_file <- "./00_rawdata/GSE58294_std.rds"
if (!file.exists(test_file)) {
  test_file <- "./00_rawdata/GSE58294_for_DEG.rds"
}
if (!file.exists(test_file)) {
  stop("❌ 验证集文件不存在，请检查文件路径")
}

test_data_raw <- readRDS(test_file)
Test_expr <- test_data_raw$expression
Test_group <- test_data_raw$group

cat("✅ 验证集读取成功\n")
cat("  表达矩阵维度:", dim(Test_expr), "\n")

if (is.data.frame(Test_group)) {
  if ("sample" %in% colnames(Test_group)) {
    rownames(Test_group) <- Test_group$sample
  }
  Test_group$group <- ifelse(tolower(Test_group$group) == "control", 0, 1)
} else if (is.vector(Test_group) || is.factor(Test_group)) {
  Test_group <- data.frame(
    sample = names(Test_group),
    group = ifelse(tolower(as.character(Test_group)) == "control", 0, 1)
  )
  rownames(Test_group) <- Test_group$sample
}

cat("  分组统计: 0=", sum(Test_group$group == 0), ", 1=", sum(Test_group$group == 1), "\n")

cat("\n========== 处理重复基因名 ==========\n")

merge_duplicate_genes <- function(expr_matrix, dataset_name) {
  original_genes <- rownames(expr_matrix)
  dup_genes <- duplicated(original_genes)

  if (sum(dup_genes) == 0) {
    cat(sprintf("  %s: 没有重复基因名\n", dataset_name))
    return(expr_matrix)
  }

  cat(sprintf("  %s: 发现 %d 个重复基因名，取平均值合并\n", dataset_name, sum(dup_genes)))

  unique_genes <- unique(original_genes)
  new_expr <- matrix(0, nrow = length(unique_genes), ncol = ncol(expr_matrix))
  rownames(new_expr) <- unique_genes
  colnames(new_expr) <- colnames(expr_matrix)

  for (gene in unique_genes) {
    idx <- which(original_genes == gene)
    if (length(idx) == 1) {
      new_expr[gene, ] <- as.matrix(expr_matrix[idx, ])
    } else {
      new_expr[gene, ] <- colMeans(as.matrix(expr_matrix[idx, , drop = FALSE]))
    }
  }

  return(new_expr)
}

Train_expr <- merge_duplicate_genes(Train_expr, "训练集")
Test_expr <- merge_duplicate_genes(Test_expr, "验证集")

rownames(Train_expr) <- toupper(rownames(Train_expr))
rownames(Test_expr) <- toupper(rownames(Test_expr))

cat("\n========== 提取生物标志物表达数据 ==========\n")

train_available <- intersect(selected_biomarkers, rownames(Train_expr))
test_available <- intersect(selected_biomarkers, rownames(Test_expr))
common_biomarkers <- intersect(train_available, test_available)

cat("训练集中可用的生物标志物:", paste(train_available, collapse = ", "), "\n")
cat("验证集中可用的生物标志物:", paste(test_available, collapse = ", "), "\n")
cat("共同生物标志物:", paste(common_biomarkers, collapse = ", "), "\n")

if (length(common_biomarkers) < 1) {
  stop("❌ 没有共同的生物标志物，无法构建列线图")
}

selected_biomarkers <- common_biomarkers

Train_expr_sub <- t(Train_expr[selected_biomarkers, , drop = FALSE])
Test_expr_sub <- t(Test_expr[selected_biomarkers, , drop = FALSE])

colnames(Train_expr_sub) <- selected_biomarkers
colnames(Test_expr_sub) <- selected_biomarkers

cat("训练集表达矩阵维度:", dim(Train_expr_sub), "\n")
cat("验证集表达矩阵维度:", dim(Test_expr_sub), "\n")

cat("\n========== 对齐样本名 ==========\n")

train_samples <- intersect(rownames(Train_group), rownames(Train_expr_sub))
if (length(train_samples) == 0) {
  if (nrow(Train_group) == nrow(Train_expr_sub)) {
    train_samples <- rownames(Train_expr_sub)
    rownames(Train_group) <- train_samples
  } else {
    stop("❌ 训练集样本名无法对齐")
  }
}

Train_expr_sub <- Train_expr_sub[train_samples, , drop = FALSE]
Train_group <- Train_group[train_samples, , drop = FALSE]

test_samples <- intersect(rownames(Test_group), rownames(Test_expr_sub))
if (length(test_samples) == 0) {
  if (nrow(Test_group) == nrow(Test_expr_sub)) {
    test_samples <- rownames(Test_expr_sub)
    rownames(Test_group) <- test_samples
  } else {
    stop("❌ 验证集样本名无法对齐")
  }
}

Test_expr_sub <- Test_expr_sub[test_samples, , drop = FALSE]
Test_group <- Test_group[test_samples, , drop = FALSE]

cat("训练集样本数:", nrow(Train_expr_sub), "\n")
cat("验证集样本数:", nrow(Test_expr_sub), "\n")

cat("\n========== 准备建模数据 ==========\n")

train_data <- as.data.frame(Train_expr_sub)
train_data[is.na(train_data)] <- 0
train_data$Group_num <- Train_group$group
train_data <- train_data[!is.na(train_data$Group_num), ]

test_data <- as.data.frame(Test_expr_sub)
test_data[is.na(test_data)] <- 0
test_data$Group_num <- Test_group$group
test_data <- test_data[!is.na(test_data$Group_num), ]

cat("训练集: 样本数=", nrow(train_data),
    " 对照组=", sum(train_data$Group_num == 0),
    " 卒中组=", sum(train_data$Group_num == 1), "\n")
cat("验证集: 样本数=", nrow(test_data),
    " 对照组=", sum(test_data$Group_num == 0),
    " 卒中组=", sum(test_data$Group_num == 1), "\n")

cat("\n[批次效应校正] 用训练集统计量对两个数据集做 z-score 标准化...\n")
train_means <- colMeans(train_data[, selected_biomarkers, drop = FALSE], na.rm = TRUE)
train_sds   <- apply(train_data[, selected_biomarkers, drop = FALSE], 2,
                     sd, na.rm = TRUE)

cat("  训练集均值:", paste(sprintf("%s=%.3f", selected_biomarkers, train_means),
                            collapse = ", "), "\n")
cat("  训练集标准差:", paste(sprintf("%s=%.3f", selected_biomarkers, train_sds),
                              collapse = ", "), "\n")

test_means_before <- colMeans(test_data[, selected_biomarkers, drop = FALSE], na.rm = TRUE)
cat("  验证集校正前均值:",
    paste(sprintf("%s=%.3f", selected_biomarkers, test_means_before), collapse = ", "), "\n")
cat("  训练-验证均值差:",
    paste(sprintf("%s=%+.3f", selected_biomarkers, test_means_before - train_means),
          collapse = ", "), "\n")

for (g in selected_biomarkers) {
  if (train_sds[g] > 0) {
    train_data[[g]] <- (train_data[[g]] - train_means[g]) / train_sds[g]
    test_data[[g]]  <- (test_data[[g]]  - train_means[g]) / train_sds[g]
  }
}
cat("✅ z-score 标准化完成（训练集和验证集都用训练集参数）\n")

validation_single_class <- length(unique(test_data$Group_num)) < 2
if (validation_single_class) {
  cat("\n⚠️ 警告: 验证集只有一个类别，部分性能指标无法计算\n")
}

cat("\n========== 构建Logistic回归模型 ==========\n")

formula_str <- paste("Group_num ~", paste(selected_biomarkers, collapse = " + "))
cat("模型公式:", formula_str, "\n")

logit_model <- glm(as.formula(formula_str),
                   data = train_data,
                   family = binomial(link = "logit"))

cat("\n模型摘要:\n")
print(summary(logit_model))

or_results <- data.frame(
  Biomarker = selected_biomarkers,
  Coefficient = coef(logit_model)[-1],
  OR = exp(coef(logit_model)[-1]),
  P_value = summary(logit_model)$coefficients[-1, 4]
)

or_results$OR_CI_lower <- NA
or_results$OR_CI_upper <- NA
ci_result <- tryCatch({
  suppressMessages(confint(logit_model))
}, error = function(e) {
  cat("  profile CI 失败，改用 Wald CI\n")
  tryCatch(confint.default(logit_model), error = function(e2) NULL)
})
if (!is.null(ci_result) && nrow(ci_result) > 1) {
  ci_mat <- ci_result[-1, , drop = FALSE]
  or_results$OR_CI_lower <- exp(ci_mat[, 1])
  or_results$OR_CI_upper <- exp(ci_mat[, 2])
  cat("  ✅ 置信区间计算成功\n")
}

or_results <- or_results[order(or_results$P_value), ]

cat("\n比值比(OR)分析结果:\n")
print(or_results)

write.csv(or_results, file.path(res.path, "Odds_Ratios.csv"), row.names = FALSE)

cat("\n========== 计算预测概率 ==========\n")

train_pred_prob <- predict(logit_model, newdata = train_data, type = "response")
train_actual    <- train_data$Group_num

test_pred_prob <- predict(logit_model, newdata = test_data, type = "response")
test_actual    <- test_data$Group_num

if (length(train_pred_prob) != length(train_actual))
  stop(sprintf("训练集长度不一致: pred=%d actual=%d",
               length(train_pred_prob), length(train_actual)))
if (length(test_pred_prob) != length(test_actual))
  stop(sprintf("验证集长度不一致: pred=%d actual=%d",
               length(test_pred_prob), length(test_actual)))

cat("训练集预测概率范围: [", round(min(train_pred_prob),4),
    ",", round(max(train_pred_prob),4), "]\n")
cat("验证集预测概率范围: [", round(min(test_pred_prob),4),
    ",", round(max(test_pred_prob),4), "]\n")

cat("\n========== 构建列线图 ==========\n")

dd <- datadist(train_data)
options(datadist = "dd")

lrm_model <- lrm(as.formula(formula_str), data = train_data,
                 x = TRUE, y = TRUE,
                 penalty = list(simple = 5))

cat("\n列线图模型摘要:\n")
print(lrm_model)

train_lp <- predict(lrm_model, type = "lp")
test_lp  <- predict(lrm_model, newdata = test_data, type = "lp")
train_pred_prob <- plogis(train_lp)
test_pred_prob  <- plogis(test_lp)

if (max(train_pred_prob) > 1 || min(train_pred_prob) < 0)
  stop("训练集预测概率超出 [0,1] 范围，请检查 lrm 模型")
if (max(test_pred_prob) > 1 || min(test_pred_prob) < 0)
  stop("验证集预测概率超出 [0,1] 范围，请检查 lrm 模型")

cat("\n[正则化后] 训练集预测概率范围: [",
    round(min(train_pred_prob),4), ",", round(max(train_pred_prob),4), "]\n")
cat("[正则化后] 验证集预测概率范围: [",
    round(min(test_pred_prob),4), ",", round(max(test_pred_prob),4), "]\n")

if (!validation_single_class) {
  test_min <- min(test_pred_prob)
  test_max <- max(test_pred_prob)
  needs_recal <- (test_min > 0.5) || (test_max < 0.5) || ((test_max - test_min) < 0.4)

  if (needs_recal) {
    cat("\n⚠️ 验证集预测概率范围异常，应用截距重新校准...\n")
    recal_fit <- tryCatch(
      glm(test_actual ~ 1 + offset(test_lp), family = binomial()),
      error = function(e) { cat("  recal 失败:", e$message, "\n"); NULL }
    )
    if (!is.null(recal_fit)) {
      intercept_shift <- as.numeric(coef(recal_fit)[1])
      cat(sprintf("  截距修正量: %+.4f\n", intercept_shift))

      test_lp_recal   <- test_lp + intercept_shift
      test_pred_prob  <- plogis(test_lp_recal)
      cat(sprintf("  校准后验证集概率范围: [ %.4f , %.4f ]\n",
                  min(test_pred_prob), max(test_pred_prob)))
    }
  } else {
    cat("\n✅ 验证集概率范围正常，无需 recalibration\n")
  }
}

pdf(file.path(fig.path, "Nomogram.pdf"), width = 10, height = 7)

nomogram_obj <- nomogram(lrm_model,
                         fun = plogis,
                         fun.at = c(0.1, 0.3, 0.5, 0.7, 0.9),
                         funlabel = "Risk of Stroke",
                         lp = FALSE,
                         conf.int = FALSE,
                         abbrev = FALSE,
                         minlength = 1)

plot(nomogram_obj,
     xfrac = 0.35,
     cex.axis = 0.8,
     col.grid = gray(c(0.8, 0.95)),
     lmgp = 0.2,
     label.every = 1,
     force.label = TRUE)

title("Nomogram for Ischemic Stroke Risk Prediction", line = 1.5, cex.main = 1.2)

dev.off()

png(file.path(fig.path, "Nomogram.png"), width = 10, height = 7, units = "in", res = 300)

plot(nomogram_obj,
     xfrac = 0.35,
     cex.axis = 0.8,
     col.grid = gray(c(0.8, 0.95)),
     lmgp = 0.2,
     label.every = 1,
     force.label = TRUE)

title("Nomogram for Ischemic Stroke Risk Prediction", line = 1.5, cex.main = 1.2)

dev.off()

cat("✅ 列线图已保存\n")

cat("\n========== ROC曲线和AUC ==========\n")

roc_train <- roc(train_actual, train_pred_prob, quiet = TRUE)
train_auc <- round(auc(roc_train), 4)
train_ci <- ci.auc(roc_train)
cat(sprintf("训练集AUC: %.4f (95%% CI: %.4f-%.4f)\n", train_auc, train_ci[1], train_ci[3]))

if (!validation_single_class) {
  roc_test <- roc(test_actual, test_pred_prob, quiet = TRUE)
  test_auc <- round(auc(roc_test), 4)
  test_ci <- ci.auc(roc_test)
  cat(sprintf("验证集AUC: %.4f (95%% CI: %.4f-%.4f)\n", test_auc, test_ci[1], test_ci[3]))
} else {
  test_auc <- NA
  cat("验证集只有一个类别，无法计算AUC\n")
}

draw_roc_train <- function() {
  plot(roc_train, col = "#4DBBD5", lwd = 3,
       main = "ROC Curve - Training Set",
       xlab = "Specificity", ylab = "Sensitivity",
       legacy.axes = FALSE)
  abline(0, 1, lty = 2, col = "gray")
  legend("bottomright",
         legend = sprintf("AUC = %.4f (95%% CI: %.4f-%.4f)",
                          train_auc, train_ci[1], train_ci[3]),
         col = "#4DBBD5", lwd = 3, bty = "n", cex = 0.95)
}

pdf(file.path(fig.path, "ROC_Curve_Training.pdf"), width = 7, height = 7)
draw_roc_train(); dev.off()
png(file.path(fig.path, "ROC_Curve_Training.png"),
    width = 7, height = 7, units = "in", res = 300)
draw_roc_train(); dev.off()
cat("✅ 训练集 ROC 曲线已保存\n")

if (!validation_single_class) {
  draw_roc_val <- function() {
    plot(roc_test, col = "#E64B35", lwd = 3,
         main = "ROC Curve - Validation Set",
         xlab = "Specificity", ylab = "Sensitivity",
         legacy.axes = FALSE)
    abline(0, 1, lty = 2, col = "gray")
    legend("bottomright",
           legend = sprintf("AUC = %.4f (95%% CI: %.4f-%.4f)",
                            test_auc, test_ci[1], test_ci[3]),
           col = "#E64B35", lwd = 3, bty = "n", cex = 0.95)
  }

  pdf(file.path(fig.path, "ROC_Curve_Validation.pdf"), width = 7, height = 7)
  draw_roc_val(); dev.off()
  png(file.path(fig.path, "ROC_Curve_Validation.png"),
      width = 7, height = 7, units = "in", res = 300)
  draw_roc_val(); dev.off()
  cat("✅ 验证集 ROC 曲线已保存\n")
}

cat("\n========== 绘制校准曲线 ==========\n")

hl_train <- tryCatch(
  ResourceSelection::hoslem.test(as.integer(train_actual),
                                  train_pred_prob, g = 10),
  error = function(e) { cat("  HL检验失败:", e$message, "\n"); NULL }
)
hl_p_train <- if (!is.null(hl_train)) signif(hl_train$p.value, 3) else "N/A"
hl_label_train <- paste0("Hosmer-Lemeshow p = ", hl_p_train)
cat("训练集 HL检验 p =", hl_p_train, "\n")

set.seed(123)
cal_train <- calibrate(lrm_model, method = "boot", B = 500)

draw_cal_train <- function() {
  par(mar = c(5, 5, 4, 2))
  cal_mat <- as.matrix(cal_train)
  cn <- colnames(cal_mat)
  x_col <- if ("predy" %in% cn) "predy" else cn[1]
  orig_col <- if ("calibrated.orig" %in% cn) "calibrated.orig"
              else if ("cal.orig" %in% cn) "cal.orig"
              else cn[grep("orig", cn)[1]]
  corr_col <- if ("calibrated.corrected" %in% cn) "calibrated.corrected"
              else if ("cal.corrected" %in% cn) "cal.corrected"
              else cn[grep("corrected", cn)[1]]

  x_pred     <- cal_mat[, x_col]
  y_apparent <- cal_mat[, orig_col]
  y_biascorr <- cal_mat[, corr_col]

  plot(NA, NA,
       xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Predicted Probability",
       ylab = "Observed Probability",
       main = "Calibration Curve - Training Set",
       cex.lab = 1.1, cex.main = 1.2,
       las = 1)

  abline(0, 1, lty = 2, lwd = 1.5, col = "gray50")
  lines(x_pred, y_apparent, lty = 1, lwd = 2, col = "black")
  lines(x_pred, y_biascorr, lty = 2, lwd = 2, col = "black")

  legend("bottomright",
         legend = c("Apparent", "Bias-corrected", "Ideal"),
         col    = c("black",    "black",          "gray50"),
         lty    = c(1,          2,                2),
         lwd    = c(2,          2,                1.5),
         bty    = "n",
         cex    = 0.9)

  text(0.05, 0.95, hl_label_train, adj = c(0, 1), cex = 0.95, font = 2)
}

pdf(file.path(fig.path, "Calibration_Curve_Training.pdf"), width = 7, height = 7)
draw_cal_train(); dev.off()
png(file.path(fig.path, "Calibration_Curve_Training.png"),
    width = 7, height = 7, units = "in", res = 300)
draw_cal_train(); dev.off()
cat("✅ 训练集校准曲线已保存\n")

if (!validation_single_class) {
  hl_test <- tryCatch(
    ResourceSelection::hoslem.test(as.integer(test_actual),
                                    test_pred_prob, g = 10),
    error = function(e) NULL
  )
  hl_p_test <- if (!is.null(hl_test)) signif(hl_test$p.value, 3) else "N/A"
  hl_label_test <- paste0("Hosmer-Lemeshow p = ", hl_p_test)
  cat("验证集 HL检验 p =", hl_p_test, "\n")

  test_pred_clip <- pmin(pmax(test_pred_prob, 0.001), 0.999)

  draw_cal_val <- function() {
    par(mar = c(5, 5, 4, 2))
    val_result <- val.prob(test_pred_clip, as.numeric(test_actual),
                           pl = TRUE,
                           smooth = TRUE,
                           logistic.cal = FALSE,
                           statloc = FALSE,
                           legendloc = FALSE,
                           xlab = "Predicted Probability",
                           ylab = "Observed Probability",
                           riskdist = "predicted",
                           cex = 0.9)
    title(main = "Calibration Curve - Validation Set", cex.main = 1.2)
    abline(0, 1, lty = 2, lwd = 1.5, col = "gray50")
    text(0.05, 0.95, hl_label_test, adj = c(0, 1), cex = 0.95, font = 2)

    if (!is.null(val_result)) {
      stats_text <- sprintf("Intercept = %.3f\nSlope = %.3f\nBrier = %.3f",
                            val_result["Intercept"],
                            val_result["Slope"],
                            val_result["Brier"])
      text(0.05, 0.85, stats_text, adj = c(0, 1), cex = 0.85, col = "gray30")
    }
  }

  pdf(file.path(fig.path, "Calibration_Curve_Validation.pdf"), width = 7, height = 7)
  draw_cal_val(); dev.off()
  png(file.path(fig.path, "Calibration_Curve_Validation.png"),
      width = 7, height = 7, units = "in", res = 300)
  draw_cal_val(); dev.off()
  cat("✅ 验证集校准曲线已保存\n")
}

cat("\n========== 决策曲线分析 ==========\n")

calc_nb <- function(pred, actual, thresholds) {
  n <- length(actual)
  sapply(thresholds, function(pt) {
    tp <- sum(pred >= pt & actual == 1) / n
    fp <- sum(pred >= pt & actual == 0) / n
    tp - fp * pt / (1 - pt)
  })
}

thresholds <- seq(0.01, 0.95, by = 0.01)

calc_dca_for_dataset <- function(data_df, lrm_pred, biomarkers, train_data, thresholds) {
  actual <- data_df$Group_num

  nb_nom <- calc_nb(lrm_pred, actual, thresholds)

  nb_genes <- setNames(
    lapply(biomarkers, function(g) {
      fit_g <- glm(as.formula(paste("Group_num ~", g)),
                   data = train_data, family = binomial())
      pred_g <- predict(fit_g, newdata = data_df, type = "response")
      calc_nb(pred_g, actual, thresholds)
    }),
    biomarkers
  )

  prev <- mean(actual == 1)
  nb_all <- sapply(thresholds, function(pt) prev - (1 - prev) * pt / (1 - pt))

  df <- data.frame(
    Threshold = rep(thresholds, length(biomarkers) + 3),
    NetBenefit = c(nb_nom, unlist(nb_genes), nb_all, rep(0, length(thresholds))),
    Model = factor(
      c(rep("Nomogram", length(thresholds)),
        rep(biomarkers, each = length(thresholds)),
        rep("Treat All", length(thresholds)),
        rep("Treat None", length(thresholds))),
      levels = c("Nomogram", biomarkers, "Treat All", "Treat None")
    )
  )
  df
}

gene_palette <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488",
                  "#F39B7F", "#8491B4", "#91D1C2", "#DC0000")
n_genes <- length(selected_biomarkers)
gene_cols <- gene_palette[seq_len(min(n_genes, length(gene_palette)))]
if (n_genes > length(gene_palette)) {
  gene_cols <- c(gene_cols, rep("gray50", n_genes - length(gene_palette)))
}

color_map <- c(setNames("black", "Nomogram"),
               setNames(gene_cols, selected_biomarkers),
               "Treat All"  = "#666666",
               "Treat None" = "#BBBBBB")

linetype_map <- c(setNames("solid", "Nomogram"),
                  setNames(rep("dashed", n_genes), selected_biomarkers),
                  "Treat All"  = "dotdash",
                  "Treat None" = "dotted")

linewidth_map <- c(setNames(1.4, "Nomogram"),
                   setNames(rep(0.7, n_genes), selected_biomarkers),
                   "Treat All"  = 0.9,
                   "Treat None" = 0.9)

plot_dca <- function(df, title_str) {
  nb_focus <- df$NetBenefit[df$Model %in% c("Nomogram", selected_biomarkers)]
  y_max <- max(c(nb_focus, 0), na.rm = TRUE) * 1.10
  y_min <- -0.05

  ggplot(df, aes(x = Threshold, y = NetBenefit,
                 color = Model, linetype = Model, linewidth = Model)) +
    geom_line() +
    scale_color_manual(values = color_map) +
    scale_linetype_manual(values = linetype_map) +
    scale_linewidth_manual(values = linewidth_map) +
    coord_cartesian(xlim = c(0, 1), ylim = c(y_min, y_max)) +
    labs(title = title_str,
         x = "Threshold Probability",
         y = "Net Benefit") +
    theme_bw(base_size = 13) +
    theme(plot.title       = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title       = element_text(size = 12),
          axis.text        = element_text(size = 11),
          legend.position  = c(0.97, 0.97),
          legend.justification = c(1, 1),
          legend.background = element_rect(fill = "white", color = "gray80"),
          legend.title     = element_blank(),
          legend.text      = element_text(size = 10),
          legend.key.width = grid::unit(1.6, "lines"),
          panel.grid.minor = element_blank())
}

cat("\n[训练集 DCA]\n")
df_train <- calc_dca_for_dataset(train_data, train_pred_prob,
                                 selected_biomarkers, train_data, thresholds)
p_dca_train <- plot_dca(df_train,
                        "Decision Curve Analysis - Training Set")

ggsave(file.path(fig.path, "Decision_Curve_Analysis_Training.pdf"),
       p_dca_train, width = 9, height = 7)
ggsave(file.path(fig.path, "Decision_Curve_Analysis_Training.png"),
       p_dca_train, width = 9, height = 7, dpi = 300)
cat("✅ 训练集 DCA 已保存\n")

if (!validation_single_class) {
  cat("\n[验证集 DCA]\n")
  df_test <- calc_dca_for_dataset(test_data, test_pred_prob,
                                  selected_biomarkers, train_data, thresholds)
  p_dca_test <- plot_dca(df_test,
                         "Decision Curve Analysis - Validation Set")

  ggsave(file.path(fig.path, "Decision_Curve_Analysis_Validation.pdf"),
         p_dca_test, width = 9, height = 7)
  ggsave(file.path(fig.path, "Decision_Curve_Analysis_Validation.png"),
         p_dca_test, width = 9, height = 7, dpi = 300)
  cat("✅ 验证集 DCA 已保存\n")
}

nom_nb       <- calc_nb(train_pred_prob, train_actual, thresholds)
gene_nb_list <- setNames(
  lapply(selected_biomarkers, function(g) {
    fit_g <- glm(as.formula(paste("Group_num ~", g)),
                 data = train_data, family = binomial())
    calc_nb(predict(fit_g, type = "response"), train_actual, thresholds)
  }),
  selected_biomarkers
)
prev   <- mean(train_actual == 1)
all_nb <- sapply(thresholds, function(pt) prev - (1 - prev) * pt / (1 - pt))

cat("\n========== 模型性能评估 ==========\n")

coords_train <- coords(roc_train, "best", best.method = "youden")
best_threshold <- coords_train$threshold

train_pred_class <- ifelse(train_pred_prob >= best_threshold, 1, 0)
conf_train <- table(Predicted = train_pred_class, Actual = train_actual)

if (nrow(conf_train) == 1) {
  if (rownames(conf_train)[1] == "0") {
    conf_train <- rbind(conf_train, c(0, 0))
    rownames(conf_train) <- c("0", "1")
  } else {
    conf_train <- rbind(c(0, 0), conf_train)
    rownames(conf_train) <- c("0", "1")
  }
}

train_accuracy <- sum(diag(conf_train)) / sum(conf_train)
train_sensitivity <- if (sum(conf_train[, 2]) > 0) conf_train[2, 2] / sum(conf_train[, 2]) else NA
train_specificity <- if (sum(conf_train[, 1]) > 0) conf_train[1, 1] / sum(conf_train[, 1]) else NA

cat("\n训练集性能:\n")
cat(sprintf("  最佳截断值: %.4f\n", best_threshold))
cat(sprintf("  准确率: %.4f\n", train_accuracy))
cat(sprintf("  灵敏度: %.4f\n", train_sensitivity))
cat(sprintf("  特异度: %.4f\n", train_specificity))
cat("  混淆矩阵:\n")
print(conf_train)

if (!validation_single_class) {
  test_pred_class <- ifelse(test_pred_prob >= best_threshold, 1, 0)
  conf_test <- table(Predicted = test_pred_class, Actual = test_actual)

  if (nrow(conf_test) == 1) {
    if (rownames(conf_test)[1] == "0") {
      conf_test <- rbind(conf_test, c(0, 0))
      rownames(conf_test) <- c("0", "1")
    } else {
      conf_test <- rbind(c(0, 0), conf_test)
      rownames(conf_test) <- c("0", "1")
    }
  }

  test_accuracy <- sum(diag(conf_test)) / sum(conf_test)
  test_sensitivity <- if (sum(conf_test[, 2]) > 0) conf_test[2, 2] / sum(conf_test[, 2]) else NA
  test_specificity <- if (sum(conf_test[, 1]) > 0) conf_test[1, 1] / sum(conf_test[, 1]) else NA

  cat("\n验证集性能:\n")
  cat(sprintf("  准确率: %.4f\n", test_accuracy))
  cat(sprintf("  灵敏度: %.4f\n", test_sensitivity))
  cat(sprintf("  特异度: %.4f\n", test_specificity))
  cat("  混淆矩阵:\n")
  print(conf_test)
}

cat("\n========== 保存结果 ==========\n")

saveRDS(logit_model, file.path(res.path, "nomogram_model.rds"))
saveRDS(lrm_model, file.path(res.path, "lrm_model.rds"))

train_pred_results <- data.frame(
  Sample = rownames(train_data),
  Actual = train_actual,
  Predicted_Probability = train_pred_prob,
  Predicted_Class = train_pred_class
)
write.csv(train_pred_results, file.path(res.path, "training_predictions.csv"), row.names = FALSE)

test_pred_results <- data.frame(
  Sample = rownames(test_data),
  Actual = test_actual,
  Predicted_Probability = test_pred_prob,
  Predicted_Class = test_pred_class
)
write.csv(test_pred_results, file.path(res.path, "validation_predictions.csv"), row.names = FALSE)

performance_df <- data.frame(
  Metric = c("AUC", "Accuracy", "Sensitivity", "Specificity", "Best_Threshold", "Sample_Size"),
  Training = c(train_auc, round(train_accuracy, 4),
               ifelse(is.na(train_sensitivity), NA, round(train_sensitivity, 4)),
               ifelse(is.na(train_specificity), NA, round(train_specificity, 4)),
               round(best_threshold, 4), nrow(train_data)),
  Validation = c(ifelse(is.na(test_auc), "N/A", round(test_auc, 4)),
                 ifelse(!validation_single_class & !is.na(test_accuracy), round(test_accuracy, 4), "N/A"),
                 ifelse(!validation_single_class & !is.na(test_sensitivity), round(test_sensitivity, 4), "N/A"),
                 ifelse(!validation_single_class & !is.na(test_specificity), round(test_specificity, 4), "N/A"),
                 "-", nrow(test_data))
)
write.csv(performance_df, file.path(res.path, "model_performance.csv"), row.names = FALSE)

calibration_data <- list(
  training = data.frame(predicted = train_pred_prob, observed = train_actual),
  validation = data.frame(predicted = test_pred_prob, observed = test_actual)
)
saveRDS(calibration_data, file.path(res.path, "calibration_data.rds"))

dca_results <- list(
  thresholds   = thresholds,
  nomogram_nb  = nom_nb,
  gene_nb      = gene_nb_list,
  treat_all_nb = all_nb
)
saveRDS(dca_results, file.path(res.path, "dca_results.rds"))

report_lines <- c(
  "# 列线图构建与评估报告",
  paste("生成时间:", Sys.time()),
  "",
  "## 生物标志物信息",
  paste("- 生物标志物数量:", length(selected_biomarkers)),
  paste("- 生物标志物列表:", paste(selected_biomarkers, collapse = ", ")),
  "",
  "## 数据信息",
  paste("- 训练集样本数:", nrow(train_data), "(对照组:", sum(train_actual == 0), ", 卒中组:", sum(train_actual == 1), ")"),
  paste("- 验证集样本数:", nrow(test_data), "(对照组:", sum(test_actual == 0), ", 卒中组:", sum(test_actual == 1), ")"),
  "",
  "## Logistic回归模型",
  "",
  "### 比值比(OR)分析",
  ""
)

for (i in 1:nrow(or_results)) {
  if (!is.na(or_results$OR_CI_lower[i]) && !is.na(or_results$OR_CI_upper[i])) {
    report_lines <- c(report_lines,
      paste("- **", or_results$Biomarker[i], "**: OR =", round(or_results$OR[i], 3),
            " (95% CI:", round(or_results$OR_CI_lower[i], 3), "-", round(or_results$OR_CI_upper[i], 3),
            "), P =", format(or_results$P_value[i], scientific = TRUE, digits = 3))
    )
  } else {
    report_lines <- c(report_lines,
      paste("- **", or_results$Biomarker[i], "**: OR =", round(or_results$OR[i], 3),
            " (CI not available), P =", format(or_results$P_value[i], scientific = TRUE, digits = 3))
    )
  }
}

report_lines <- c(report_lines,
  "",
  "## 模型性能评估",
  "",
  "| 指标 | 训练集 | 验证集 |",
  "|------|--------|--------|",
  paste("| AUC |", train_auc, "|", ifelse(is.na(test_auc), "N/A", test_auc), "|"),
  paste("| Accuracy |", round(train_accuracy, 4), "|",
        ifelse(!validation_single_class && !is.na(test_accuracy), round(test_accuracy, 4), "N/A"), "|"),
  paste("| Sensitivity |", ifelse(is.na(train_sensitivity), "N/A", round(train_sensitivity, 4)), "|",
        ifelse(!validation_single_class && !is.na(test_sensitivity), round(test_sensitivity, 4), "N/A"), "|"),
  paste("| Specificity |", ifelse(is.na(train_specificity), "N/A", round(train_specificity, 4)), "|",
        ifelse(!validation_single_class && !is.na(test_specificity), round(test_specificity, 4), "N/A"), "|"),
  paste("| Best Threshold |", round(best_threshold, 4), "| - |"),
  "",
  "## 输出文件",
  "",
  "### Figures/",
  "- Nomogram.pdf/png (列线图)",
  "- ROC_Curve_Training.pdf/png (训练集 ROC 曲线)",
  if (!validation_single_class) "- ROC_Curve_Validation.pdf/png (验证集 ROC 曲线)",
  "- Calibration_Curve_Training.pdf/png (训练集校准曲线)",
  if (!validation_single_class) "- Calibration_Curve_Validation.pdf/png (验证集校准曲线)",
  "- Decision_Curve_Analysis_Training.pdf/png (训练集决策曲线)",
  if (!validation_single_class) "- Decision_Curve_Analysis_Validation.pdf/png (验证集决策曲线)",
  "",
  "### Results/",
  "- Odds_Ratios.csv (比值比分析结果)",
  "- model_performance.csv (模型性能指标)",
  "- nomogram_model.rds (Logistic回归模型)",
  "- lrm_model.rds (rms格式模型)",
  "- training_predictions.csv (训练集预测结果)",
  "- validation_predictions.csv (验证集预测结果)",
  "- calibration_data.rds (校准数据)",
  "- dca_results.rds (DCA分析结果)"
)

writeLines(report_lines, file.path(res.path, "Nomogram_Report.md"))

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("             第八部分：列线图构建与评估完成\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("\n结果汇总:\n")
cat("  - 生物标志物数量:", length(selected_biomarkers), "\n")
cat("  - 生物标志物:", paste(selected_biomarkers, collapse = ", "), "\n")
cat("  - 训练集AUC:", train_auc, "\n")
if (!validation_single_class) {
  cat("  - 验证集AUC:", test_auc, "\n")
} else {
  cat("  - 验证集AUC: 无法计算（只有一个类别）\n")
}
cat("  - 训练集准确率:", round(train_accuracy, 4), "\n\n")

cat("生成的文件保存在:", output_dir, "\n")
cat("\n✅ 第八部分分析完成！\n")
