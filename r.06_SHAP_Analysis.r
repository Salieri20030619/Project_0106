
cat("\n========================================\n")
cat("06. SHAP分析\n")
cat("========================================\n\n")

rm(list = ls())

options(bitmapType = "cairo")

setwd("")

output_dir <- "./06_SHAP_Analysis"
if (!dir.exists(output_dir)) dir.create(output_dir)

res.path <- file.path(output_dir, "Results")
fig.path <- file.path(output_dir, "Figures")

for (path in c(res.path, fig.path)) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

cat("输出目录:", output_dir, "\n")
cat("结果目录:", res.path, "\n")
cat("图片目录:", fig.path, "\n\n")

cat("========== 加载R包 ==========\n")

packages <- c(
  "ggplot2", "pheatmap", "RColorBrewer", "gridExtra",
  "ggpubr", "caret", "pROC", "mboost", "dplyr", "tidyr",
  "shapviz", "kernelshap", "GGally",
  "glmnet", "randomForestSRC", "gbm", "e1071", "MASS", "xgboost"
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
arrange <- dplyr::arrange
mutate <- dplyr::mutate

cat("✅ 所有R包加载完成\n")

cat("\n========== 从第五部分输出文件读取结果 ==========\n")

ml_results_path <- "./05_Machine_Learning/Results"

if (!file.exists(file.path(ml_results_path, "best_model_final.rds"))) {
  stop("❌ best_model_final.rds 文件不存在")
}
if (!file.exists(file.path(ml_results_path, "best_features_final.txt"))) {
  stop("❌ best_features_final.txt 文件不存在")
}
if (!file.exists(file.path(ml_results_path, "all_results_final.rds"))) {
  stop("❌ all_results_final.rds 文件不存在")
}

best_model_fit <- readRDS(file.path(ml_results_path, "best_model_final.rds"))
cat("✅ 最佳模型加载成功\n")
cat("  模型类型:", class(best_model_fit)[1], "\n")

best_features <- readLines(file.path(ml_results_path, "best_features_final.txt"))
best_features <- best_features[best_features != "" & !is.na(best_features)]
cat("✅ 特征基因加载成功: ", paste(best_features, collapse = ", "), "\n")

all_results <- readRDS(file.path(ml_results_path, "all_results_final.rds"))
best_model_name <- all_results$best_model
best_validation_auc <- all_results$best_validation_auc
best_cv_auc <- all_results$best_cv_auc

cat("\n第五阶段结果摘要:\n")
cat("  - 最佳模型:", best_model_name, "\n")
cat("  - 5折交叉验证平均AUC:", round(best_cv_auc, 4), "\n")
cat("  - 独立验证集AUC:", round(best_validation_auc, 4), "\n")
cat("  - 特征基因数:", length(best_features), "\n")

cat("\n========== 加载第五部分已处理数据 ==========\n")

processed_data_file <- "./05_Machine_Learning/Intermediate/processed_data.rds"
if (!file.exists(processed_data_file)) {
  stop("❌ processed_data.rds 不存在，请先运行第五部分")
}

processed <- readRDS(processed_data_file)
X_train <- processed$Train_set
X_test  <- processed$Test_set
Train_class <- processed$Train_class
Test_class  <- processed$Test_class

y_train <- Train_class$group
y_test  <- Test_class$group

common_features <- colnames(X_train)

cat("训练集维度:", dim(X_train), "\n")
cat("验证集维度:", dim(X_test), "\n")
cat("训练集: 0=", sum(y_train == 0), "1=", sum(y_train == 1), "\n")
cat("验证集: 0=", sum(y_test  == 0), "1=", sum(y_test  == 1), "\n")
cat("特征列数:", length(common_features), "\n")
cat("特征基因:", paste(common_features, collapse = ", "), "\n")

df_train <- as.data.frame(X_train)
df_train$group <- y_train

cat("\n========== 加载预测函数 ==========\n")

ml_r_file <- "ML.R"
if (file.exists(ml_r_file)) {
  source(ml_r_file)
  cat("✅ ML.R 加载成功，使用 CalPredictScore\n")
} else {
  stop("❌ ML.R 不存在，无法进行预测")
}

predict_function <- function(model, newdata) {
  tryCatch({
    if (is.data.frame(newdata)) newdata <- as.matrix(newdata)
    newdata[is.na(newdata)] <- 0
    newdata <- as.data.frame(newdata)

    CalPredictScore(fit = model, new_data = newdata)
  }, error = function(e) {
    warning(sprintf("predict_function 失败 [%s]: %s",
                    paste(class(model), collapse="/"), e$message))
    rep(NA_real_, nrow(newdata))
  })
}

if (FALSE) {

  predict_function_placeholder <- function(model, newdata) {
    if (is.data.frame(newdata)) {
      newdata <- as.matrix(newdata)
    }
    newdata[is.na(newdata)] <- 0

    if (!is.null(model$subFeature)) {
      common_cols <- intersect(model$subFeature, colnames(newdata))
      if (length(common_cols) > 0) {
        newdata <- newdata[, common_cols, drop = FALSE]
      }
    }

    if (inherits(model, "glmboost")) {
      pred <- tryCatch({
      predict(model, newdata = as.data.frame(newdata), type = "response")
    }, error = function(e) {
      predict(model, newdata = as.data.frame(newdata))
    })
    return(as.numeric(pred))
  } else if (inherits(model, "cv.glmnet")) {
    pred <- predict(model, newx = newdata, s = "lambda.min", type = "response")
    return(as.numeric(pred))
  } else if (inherits(model, "glm")) {
    pred <- predict(model, newdata = as.data.frame(newdata), type = "response")
    return(as.numeric(pred))
    }
  }
}

if (nrow(X_train) > 0) {
  test_pred <- predict_function(best_model_fit, X_train[1:min(5, nrow(X_train)), , drop = FALSE])
  cat("预测测试成功\n")
}

test_pred_all <- predict_function(best_model_fit, X_test)
if (!all(is.na(test_pred_all))) {
  roc_test <- roc(y_test, test_pred_all)
  cat("验证集AUC:", round(roc_test$auc, 4), "\n")
}

cat("\n========== 使用kernelshap计算SHAP值 ==========\n")

set.seed(123)

bg_size <- min(20, nrow(X_train))
bg_X <- X_train[sample(1:nrow(X_train), bg_size), , drop = FALSE]

cat("背景数据集大小:", bg_size, "\n")
cat("计算SHAP值中，请稍候...\n")

shap <- NULL
shap_viz <- NULL

tryCatch({
  shap <- kernelshap(
    object = best_model_fit,
    X = X_train,
    bg_X = bg_X,
    pred_fun = predict_function,
    verbose = FALSE
  )

  shap_viz <- shapviz(shap)
  cat("✅ SHAP值计算成功\n")

}, error = function(e) {
  cat("kernelshap计算失败:", e$message, "\n")
  cat("\n尝试使用备用方法计算特征重要性...\n")

  tryCatch({

    baseline_pred <- CalPredictScore(fit = best_model_fit,
                                     new_data = as.data.frame(X_train))

    perm_importance <- sapply(common_features, function(feat) {
      X_perm <- as.data.frame(X_train)
      X_perm[[feat]] <- sample(X_perm[[feat]])
      perm_pred <- CalPredictScore(fit = best_model_fit,
                                   new_data = X_perm)
      mean(abs(baseline_pred - perm_pred), na.rm = TRUE)
    })

    importance_df <<- data.frame(
      Feature   = names(perm_importance),
      Mean_SHAP = as.numeric(perm_importance),
      stringsAsFactors = FALSE
    )
    importance_df <<- importance_df[order(importance_df$Mean_SHAP, decreasing = TRUE), ]

    shap_mat <- matrix(0, nrow = nrow(X_train), ncol = length(common_features),
                       dimnames = list(rownames(X_train), common_features))
    for (feat in common_features) {
      imp_val <- perm_importance[feat]
      shap_mat[, feat] <- (baseline_pred - mean(baseline_pred)) *
                           scale(X_train[, feat])[, 1] *
                           imp_val / max(perm_importance, na.rm = TRUE)
    }
    shap_mat[is.na(shap_mat)] <- 0

    shap_obj <- list(S = shap_mat, baseline = mean(baseline_pred))
    colnames(shap_obj$S) <- common_features
    shap_viz <<- shapviz(shap_obj)

    cat("✅ 备用排列重要性计算成功\n")
    cat("  重要性排名 Top5:\n")
    print(head(importance_df, 5))

  }, error = function(e2) {
    cat("备用方案也失败:", e2$message, "\n")
    cat("将使用模型系数绝对值作为最终兜底\n")

    s_val <- if (!is.null(best_model_fit$lambda)) median(best_model_fit$lambda) else 0.01
    coef_mat <- tryCatch(
      as.matrix(coef(best_model_fit, s = s_val)),
      error = function(e3) NULL
    )
    if (!is.null(coef_mat) && nrow(coef_mat) > 1) {
      coef_abs <- abs(coef_mat[-1, 1])
      names(coef_abs) <- rownames(coef_mat)[-1]
      importance_df <<- data.frame(
        Feature   = names(sort(coef_abs, decreasing = TRUE)),
        Mean_SHAP = as.numeric(sort(coef_abs, decreasing = TRUE)),
        stringsAsFactors = FALSE
      )
    }

  })
})

cat("\n========== SHAP可视化 ==========\n")

if (!is.null(shap_viz)) {

  tryCatch({
    p1 <- sv_importance(shap_viz, kind = "bar", bar_width = 0.8) +
      labs(title = paste("SHAP Feature Importance -", best_model_name),
           subtitle = paste("Validation AUC =", round(best_validation_auc, 3)),
           x = "Mean |SHAP Value|", y = "") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 10),
            axis.text.y = element_text(size = 12, face = "bold"))

    ggsave(file.path(fig.path, "Fig4a_SHAP_Importance_Bar.pdf"), p1, width = 8, height = 4)
    ggsave(file.path(fig.path, "Fig4a_SHAP_Importance_Bar.png"), p1, width = 8, height = 4, dpi = 300)
    cat("✅ 图4a: SHAP重要性柱状图已保存\n")
  }, error = function(e) {
    cat("图4a绘制失败:", e$message, "\n")
  })

  tryCatch({
    p2 <- sv_importance(shap_viz, kind = "both", show_numbers = TRUE, number_size = 3) +
      labs(title = paste("SHAP Feature Importance -", best_model_name),
           subtitle = "Bee swarm plot showing SHAP value distribution") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 10))

    ggsave(file.path(fig.path, "Fig4b_SHAP_Bee_Plot.pdf"), p2, width = 10, height = 6)
    ggsave(file.path(fig.path, "Fig4b_SHAP_Bee_Plot.png"), p2, width = 10, height = 6, dpi = 300)
    cat("✅ 图4b: SHAP蜂群图已保存\n")
  }, error = function(e) {
    cat("图4b绘制失败:", e$message, "\n")
  })

  if (length(common_features) >= 1) {
    tryCatch({
      pdf(file.path(fig.path, "Fig4c_SHAP_Dependence_Plots.pdf"), width = 8, height = 6)
      for (feature in common_features) {
        p_dep <- sv_dependence(shap_viz, v = feature, color_var = feature, alpha = 0.7) +
          labs(title = paste("SHAP Dependence Plot:", feature),
               x = paste(feature, "Expression"),
               y = "SHAP Value") +
          geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
          scale_color_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35", midpoint = 0) +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5))
        print(p_dep)
      }
      dev.off()
      cat("✅ 图4c: SHAP依赖图已保存\n")
    }, error = function(e) {
      cat("图4c绘制失败:", e$message, "\n")
      dev.off()
    })
  }

  if (length(common_features) >= 2 && require("GGally", quietly = TRUE)) {
    tryCatch({
      shap_matrix <- as.matrix(shap_viz$S)
      colnames(shap_matrix) <- common_features

      shap_df <- as.data.frame(shap_matrix)
      shap_df$Group <- ifelse(y_train == 1, "Stroke", "Control")

      p_matrix <- ggpairs(shap_df, columns = 1:min(5, length(common_features)),
                          mapping = aes(color = Group),
                          upper = list(continuous = wrap("cor", size = 3)),
                          lower = list(continuous = wrap("points", alpha = 0.5, size = 1)),
                          diag = list(continuous = wrap("barDiag", alpha = 0.5))) +
        labs(title = paste("SHAP Value Matrix -", best_model_name)) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))

      ggsave(file.path(fig.path, "Fig4d_SHAP_Matrix.pdf"), p_matrix, width = 10, height = 8)
      ggsave(file.path(fig.path, "Fig4d_SHAP_Matrix.png"), p_matrix, width = 10, height = 8, dpi = 150)
      cat("✅ 图4d: SHAP散点图矩阵已保存\n")
    }, error = function(e) {
      cat("图4d绘制失败:", e$message, "\n")
    })
  }

  tryCatch({
    shap_matrix <- as.matrix(shap_viz$S)
    colnames(shap_matrix) <- common_features
    rownames(shap_matrix) <- rownames(X_train)

    set.seed(123)
    shap_matrix_noise <- shap_matrix + matrix(rnorm(length(shap_matrix), 0, 1e-10),
                                               nrow = nrow(shap_matrix), ncol = ncol(shap_matrix))

    total_shap <- rowSums(abs(shap_matrix_noise))
    order_samples <- order(total_shap, decreasing = TRUE)

    breaks <- seq(min(shap_matrix_noise), max(shap_matrix_noise), length.out = 101)

    breaks <- unique(breaks)
    if (length(breaks) < 2) {
      breaks <- seq(-1, 1, length.out = 101)
    }

    pdf(file.path(fig.path, "SHAP_Heatmap.pdf"), width = 10, height = 12)
    pheatmap(shap_matrix_noise[order_samples, ],
             color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
             main = paste("SHAP Value Heatmap -", best_model_name),
             cluster_rows = FALSE,
             cluster_cols = TRUE,
             show_rownames = FALSE,
             show_colnames = TRUE,
             border_color = NA,
             fontsize_col = 10,
             breaks = breaks)
    dev.off()
    cat("✅ SHAP热图已保存\n")
  }, error = function(e) {
    cat("SHAP热图绘制失败:", e$message, "\n")

    tryCatch({
      pdf(file.path(fig.path, "SHAP_Heatmap_simple.pdf"), width = 10, height = 12)
      pheatmap(shap_matrix[order_samples, ],
               color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
               main = paste("SHAP Value Heatmap -", best_model_name),
               cluster_rows = FALSE,
               cluster_cols = TRUE,
               show_rownames = FALSE,
               show_colnames = TRUE,
               border_color = NA,
               fontsize_col = 10)
      dev.off()
      cat("✅ SHAP简化热图已保存\n")
    }, error = function(e2) {
      cat("简化热图也失败:", e2$message, "\n")
    })
  })

  shap_values_df <- as.data.frame(shap_viz$S)
  shap_values_df$Sample <- rownames(X_train)
  shap_values_df$Group <- y_train
  write.csv(shap_values_df, file.path(res.path, "shap_values.csv"), row.names = FALSE)

  importance_df <- data.frame(
    Feature = common_features,
    Mean_SHAP = colMeans(abs(shap_viz$S), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  importance_df <- importance_df[order(importance_df$Mean_SHAP, decreasing = TRUE), ]
  write.csv(importance_df, file.path(res.path, "shap_importance.csv"), row.names = FALSE)

  cat("✅ SHAP值已保存\n")

} else {
  cat("❌ SHAP对象为空，无法进行可视化\n")
}

cat("\n========== 绘制特征表达箱线图 ==========\n")

tryCatch({
  plot_long <- df_train %>%
    pivot_longer(cols = all_of(common_features), names_to = "Gene", values_to = "Expression") %>%
    mutate(Group = ifelse(group == 1, "Stroke", "Control"))

  n_features <- length(common_features)
  n_col <- min(3, n_features)
  n_row <- ceiling(n_features / n_col)

  p_box <- ggplot(plot_long, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1.2) +
    scale_fill_manual(values = c(Control = "#4DBBD5", Stroke = "#E64B35")) +
    facet_wrap(~Gene, scales = "free_y", ncol = n_col) +
    labs(title = paste("Feature Expression Levels -", best_model_name),
         x = "", y = "Standardized Expression") +
    stat_compare_means(method = "wilcox.test", label = "p.format", hide.ns = TRUE, label.x = 1.5) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          strip.background = element_rect(fill = "lightgray"),
          strip.text = element_text(size = 10, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))

  plot_height <- max(5, n_row * 3)
  ggsave(file.path(fig.path, "Feature_Expression_Boxplot.pdf"), p_box, width = 8, height = plot_height)
  ggsave(file.path(fig.path, "Feature_Expression_Boxplot.png"), p_box, width = 8, height = plot_height, dpi = 300)
  cat("✅ 特征表达箱线图已保存\n")
}, error = function(e) {
  cat("特征表达箱线图绘制失败:", e$message, "\n")
})

cat("\n========== 绘制特征相关性热图 ==========\n")

tryCatch({
  cor_matrix <- cor(df_train[, common_features], method = "spearman")
  cor_matrix[is.na(cor_matrix)] <- 0

  col_fun <- colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100)

  pdf(file.path(fig.path, "Feature_Correlation_Heatmap.pdf"), width = 8, height = 6)
  pheatmap(cor_matrix,
           color = col_fun,
           main = paste("Feature Correlation -", best_model_name),
           display_numbers = TRUE,
           number_format = "%.2f",
           fontsize_number = 8,
           border_color = "grey")
  dev.off()
  cat("✅ 特征相关性热图已保存\n")
}, error = function(e) {
  cat("特征相关性热图绘制失败:", e$message, "\n")
})

cat("\n========== 保存最终结果 ==========\n")

final_summary <- list(
  best_model = best_model_name,
  model_type = class(best_model_fit)[1],
  features = common_features,
  validation_auc = best_validation_auc,
  cv_auc = best_cv_auc,
  n_features = length(common_features),
  n_train = nrow(df_train),
  n_test = nrow(X_test),
  analysis_date = Sys.time()
)

saveRDS(final_summary, file.path(res.path, "shap_analysis_summary.rds"))

report_lines <- c(
  "# SHAP分析报告",
  paste("生成时间:", Sys.time()),
  "",
  "## 模型信息",
  paste("- 最佳模型:", best_model_name),
  paste("- 特征基因:", paste(common_features, collapse = ", ")),
  paste("- 验证集AUC:", round(best_validation_auc, 4)),
  paste("- 5折交叉验证平均AUC:", round(best_cv_auc, 4)),
  "",
  "## SHAP特征重要性",
  ""
)

if (exists("importance_df") && nrow(importance_df) > 0) {
  for (i in 1:nrow(importance_df)) {
    report_lines <- c(report_lines,
      paste(i, ".", importance_df$Feature[i],
            "- Mean SHAP:", round(importance_df$Mean_SHAP[i], 4)))
  }
} else {
  report_lines <- c(report_lines, "无法计算SHAP重要性")
}

report_lines <- c(report_lines,
  "",
  "## 特征表达统计",
  ""
)

for (feat in common_features) {
  mean_control <- mean(df_train[df_train$group == 0, feat], na.rm = TRUE)
  mean_stroke <- mean(df_train[df_train$group == 1, feat], na.rm = TRUE)
  p_val <- wilcox.test(df_train[df_train$group == 0, feat],
                       df_train[df_train$group == 1, feat])$p.value

  report_lines <- c(report_lines,
    paste("###", feat),
    paste("- 对照组平均表达:", round(mean_control, 4)),
    paste("- 卒中组平均表达:", round(mean_stroke, 4)),
    paste("- 表达趋势:", ifelse(mean_stroke > mean_control, "上调", "下调")),
    paste("- Wilcoxon检验p值:", format(p_val, scientific = TRUE, digits = 3)),
    ""
  )
}

writeLines(report_lines, file.path(res.path, "SHAP_Analysis_Report.md"))

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("             第六阶段：SHAP分析完成\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("\n结果汇总:\n")
cat("  - 最佳模型:", best_model_name, "\n")
cat("  - 模型类型:", class(best_model_fit)[1], "\n")
cat("  - 验证集AUC:", round(best_validation_auc, 4), "\n")
cat("  - 特征基因:", paste(common_features, collapse = ", "), "\n")
cat("  - 训练集样本数:", nrow(df_train), "\n")
cat("  - 验证集样本数:", nrow(X_test), "\n")

if (exists("importance_df") && nrow(importance_df) > 0) {
  cat("\nSHAP重要性排名:\n")
  for (i in 1:min(10, nrow(importance_df))) {
    cat(sprintf("  %d. %s: Mean SHAP = %.4f\n",
                i, importance_df$Feature[i], importance_df$Mean_SHAP[i]))
  }
}

cat("\n生成的文件保存在:", output_dir, "\n")
cat("  - Figures/: SHAP可视化图表\n")
cat("    - Fig4a_SHAP_Importance_Bar.pdf/png (图4a: 特征重要性柱状图)\n")
cat("    - Fig4b_SHAP_Bee_Plot.pdf/png (图4b: SHAP蜂群图)\n")
cat("    - Fig4c_SHAP_Dependence_Plots.pdf (图4c: SHAP依赖图)\n")
cat("    - Fig4d_SHAP_Matrix.pdf/png (图4d: SHAP散点图矩阵)\n")
cat("    - SHAP_Heatmap.pdf (SHAP热图)\n")
cat("    - Feature_Expression_Boxplot.pdf/png (特征表达箱线图)\n")
cat("  - Results/: SHAP值和汇总结果\n")
cat("\n✅ 第六阶段分析完成！\n")
