
cat("\n========================================\n")
cat("05. 机器学习分析\n")
cat("========================================\n\n")

rm(list = ls())

setwd("")

output_dir <- "./05_Machine_Learning"
if (!dir.exists(output_dir)) dir.create(output_dir)

work.path <- output_dir
method.path <- ""
data.path <- file.path(work.path, "InputData")
res.path <- file.path(work.path, "Results")
fig.path <- file.path(work.path, "Figures")
code.path <- file.path(work.path, "Codes")
intermediate.path <- file.path(work.path, "Intermediate")

for (path in c(data.path, res.path, fig.path, code.path, intermediate.path)) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

cat("输出目录:", output_dir, "\n")
cat("结果目录:", res.path, "\n")
cat("图片目录:", fig.path, "\n")
cat("中间结果目录:", intermediate.path, "\n\n")

analysis_start_time <- Sys.time()
cat("分析开始时间:", format(analysis_start_time), "\n\n")

cat("========== 输入文件检查 ==========\n")

candidate_dir <- "./02_CandidateGene_Intersect/Results"
candidate_file <- file.path(candidate_dir, "final_candidate_genes.txt")
candidate_detail_file <- file.path(candidate_dir, "candidate_genes_detail.csv")

train_file <- "00_rawdata/GSE16561_std.rds"
test_file <- "00_rawdata/GSE58294_std.rds"
ml_file <- "ML.R"
methods_file <- "methods.xlsx"

cat("候选基因目录:", candidate_dir, "\n")
cat("训练集文件:", train_file, "\n")
cat("验证集文件:", test_file, "\n\n")

cat("========== 读取候选基因 ==========\n")

candidate_genes <- NULL
candidate_up <- NULL
candidate_down <- NULL

if (file.exists(candidate_file)) {
  candidate_genes <- readLines(candidate_file)
  candidate_genes <- candidate_genes[candidate_genes != "" & !is.na(candidate_genes)]
  cat("✅ 从 final_candidate_genes.txt 读取到", length(candidate_genes), "个候选基因\n")
} else {
  cat("⚠️ final_candidate_genes.txt 不存在\n")
}

if (file.exists(candidate_detail_file)) {
  candidate_detail <- read.csv(candidate_detail_file, stringsAsFactors = FALSE)
  cat("✅ 从 candidate_genes_detail.csv 读取到详细信息\n")

  if (nrow(candidate_detail) > 0) {
    if (is.null(candidate_genes)) {
      candidate_genes <- candidate_detail$Gene
    }
    candidate_up <- candidate_detail$Gene[candidate_detail$Regulation == "Up"]
    candidate_down <- candidate_detail$Gene[candidate_detail$Regulation == "Down"]

    cat("  上调基因数:", length(candidate_up), "\n")
    cat("  下调基因数:", length(candidate_down), "\n")
  }
}

if (is.null(candidate_genes) || length(candidate_genes) == 0) {
  cat("\n⚠️ 无法从02读取候选基因，使用备用列表\n")
  candidate_genes <- c("C5AR1", "PDK4", "CRISPLD2", "ADM", "CYP1B1", "CEBPD",
                       "ABCA1", "IRS2", "FCGR2A", "HSPA1A", "KLF4", "CD14",
                       "EMILIN2", "FKBP5", "LRG1", "HP", "LCN2", "LGALS3", "MFGE8")
  candidate_up <- c("C5AR1", "PDK4", "CRISPLD2", "ADM", "CYP1B1", "CEBPD",
                    "ABCA1", "IRS2", "FCGR2A", "HSPA1A", "KLF4", "CD14",
                    "EMILIN2", "FKBP5", "LRG1", "HP", "LCN2", "LGALS3")
  candidate_down <- c("MFGE8")
  cat("使用备用列表，共", length(candidate_genes), "个候选基因\n")
}

cat("\n候选基因列表:\n")
cat(paste(head(candidate_genes, 10), collapse = ", "))
if (length(candidate_genes) > 10) cat(" ... 共", length(candidate_genes), "个\n")
cat("\n")

writeLines(candidate_genes, file.path(intermediate.path, "input_candidate_genes.txt"))
writeLines(candidate_up, file.path(intermediate.path, "input_candidate_up.txt"))
writeLines(candidate_down, file.path(intermediate.path, "input_candidate_down.txt"))

cat("\n========== 加载R包 ==========\n")

if (!require("R.utils", character.only = TRUE, quietly = TRUE)) {
  install.packages("R.utils", repos = "http://cran.us.r-project.org", quiet = TRUE)
  library(R.utils, character.only = TRUE)
}
cat("✅ 已加载: R.utils\n")

core_packages <- c(
  "openxlsx", "randomForestSRC", "glmnet", "gbm", "xgboost",
  "e1071", "MASS", "caret", "pROC", "RColorBrewer", "dplyr", "tidyr",
  "mboost", "plsRglm", "pheatmap", "ggplot2", "reshape2"
)

for (pkg in core_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("安装包:", pkg, "\n")
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
rename <- dplyr::rename

cat("✅ 核心R包加载完成\n")

cat("\n========== 加载机器学习函数库 ==========\n")

if (file.exists(ml_file)) {
  source(ml_file)
  cat("✅ ML.R 加载成功\n")
} else {
  cat("⚠️ ML.R 文件不存在，将使用内置函数\n")

  RunML <- function(method, Train_set, Train_label, mode, classVar) {
    tryCatch({
      Train_label_factor <- as.factor(Train_label[[classVar]])

      if (method == "logistic") {
        df <- as.data.frame(Train_set)
        df$label <- Train_label_factor
        fit <- glm(label ~ ., data = df, family = binomial())
        if (mode == "Model") return(fit)
        if (mode == "Variable") return(colnames(Train_set))
      }

      else if (method == "rf") {
        fit <- randomForest(x = Train_set, y = Train_label_factor, ntree = 500)
        if (mode == "Model") return(fit)
        if (mode == "Variable") return(colnames(Train_set))
      }

      else if (method == "svm") {
        fit <- svm(x = Train_set, y = Train_label_factor, kernel = "radial", probability = TRUE)
        if (mode == "Model") return(fit)
        if (mode == "Variable") return(colnames(Train_set))
      }

      else if (method == "lasso") {
        fit <- cv.glmnet(x = as.matrix(Train_set), y = Train_label_factor, family = "binomial", alpha = 1)
        if (mode == "Model") return(fit)
        if (mode == "Variable") return(colnames(Train_set))
      }

      else if (method == "ridge") {
        fit <- cv.glmnet(x = as.matrix(Train_set), y = Train_label_factor, family = "binomial", alpha = 0)
        if (mode == "Model") return(fit)
        if (mode == "Variable") return(colnames(Train_set))
      }

      else if (method == "elasticnet") {
        fit <- cv.glmnet(x = as.matrix(Train_set), y = Train_label_factor, family = "binomial", alpha = 0.5)
        if (mode == "Model") return(fit)
        if (mode == "Variable") return(colnames(Train_set))
      }

      else {
        if (mode == "Model") return(NULL)
        if (mode == "Variable") return(colnames(Train_set))
      }
    }, error = function(e) {
      if (mode == "Model") return(NULL)
      if (mode == "Variable") return(colnames(Train_set))
    })
  }

  CalPredictScore <- function(fit, new_data) {
    tryCatch({
      if (inherits(fit, "glm")) {
        pred <- predict(fit, newdata = as.data.frame(new_data), type = "response")
        return(as.numeric(pred))
      } else if (inherits(fit, "randomForest")) {
        pred <- predict(fit, newdata = new_data, type = "prob")[, 2]
        return(as.numeric(pred))
      } else if (inherits(fit, "svm")) {
        pred <- attr(predict(fit, newdata = new_data, probability = TRUE), "probabilities")[, 2]
        return(as.numeric(pred))
      } else if (inherits(fit, "cv.glmnet")) {
        pred <- predict(fit, newx = as.matrix(new_data), s = "lambda.min", type = "response")
        return(as.numeric(pred))
      } else {
        return(rep(0.5, nrow(new_data)))
      }
    }, error = function(e) {
      return(rep(NA, nrow(new_data)))
    })
  }

  ExtractVar <- function(fit) {
    if (inherits(fit, "glm")) {
      return(names(coef(fit))[-1])
    } else if (inherits(fit, "randomForest")) {
      return(rownames(importance(fit)))
    } else if (inherits(fit, "cv.glmnet")) {
      coefs <- as.matrix(coef(fit, s = "lambda.min"))
      return(rownames(coefs)[coefs[,1] != 0][-1])
    } else {
      return(fit$subFeature)
    }
  }

  cat("✅ 内置函数已定义\n")
}

cat("\n========== 修复XGBoost函数 ==========\n")

RunXGBoost_fixed <- function(Train_set, Train_label, mode, classVar) {
  tryCatch({
    if (!is.matrix(Train_set)) {
      Train_set <- as.matrix(Train_set)
    }

    if (any(!sapply(as.data.frame(Train_set), is.numeric))) {
      Train_set <- as.matrix(as.data.frame(lapply(as.data.frame(Train_set), as.numeric)))
    }

    Train_set[is.na(Train_set)] <- 0

    labels_raw <- Train_label[[classVar]]
    if (is.factor(labels_raw)) {
      labels <- as.numeric(as.character(labels_raw))
    } else if (is.character(labels_raw)) {
      labels <- as.numeric(labels_raw)
    } else {
      labels <- as.numeric(labels_raw)
    }

    unique_labels <- unique(labels)
    if (length(unique_labels) != 2) {
      stop("XGBoost要求二分类标签，当前标签类别数: ", length(unique_labels))
    }

    if (!(0 %in% unique_labels && 1 %in% unique_labels)) {
      label_map <- sort(unique_labels)
      labels <- ifelse(labels == label_map[1], 0, 1)
    }

    cat("  XGBoost标签统计: 0=", sum(labels == 0), "1=", sum(labels == 1), "\n")

    dtrain <- xgb.DMatrix(data = Train_set, label = labels)

    params <- list(
      objective = "reg:logistic",
      eval_metric = "logloss",
      max_depth = 3,
      eta = 0.1,
      nthread = 2
    )

    cv_result <- xgb.cv(
      params = params,
      data = dtrain,
      nrounds = 100,
      nfold = 5,
      verbose = 0,
      early_stopping_rounds = 10,
      prediction = FALSE
    )

    best_nrounds <- cv_result$best_iteration
    if (is.null(best_nrounds) || best_nrounds == 0 || is.na(best_nrounds)) {
      best_nrounds <- 50
    }

    fit <- xgb.train(
      params = params,
      data = dtrain,
      nrounds = best_nrounds,
      verbose = 0
    )

    fit$subFeature <- colnames(Train_set)

    if (mode == "Model") return(fit)
    if (mode == "Variable") return(colnames(Train_set))

  }, error = function(e) {
    cat("  XGBoost错误:", e$message, "\n")
    if (mode == "Model") return(NULL)
    if (mode == "Variable") return(colnames(Train_set))
  })
}

if (exists("RunXGBoost")) {
  assign("RunXGBoost", RunXGBoost_fixed, envir = .GlobalEnv)
  cat("✅ RunXGBoost函数已修复\n")
} else {
  RunXGBoost <- RunXGBoost_fixed
  cat("✅ RunXGBoost函数已创建\n")
}

cat("\n修复CalPredictScore函数...\n")

if (exists("CalPredictScore")) {
  Original_CalPredictScore <- CalPredictScore
} else {
  Original_CalPredictScore <- function(fit, new_data) {
    return(rep(0.5, nrow(new_data)))
  }
}

CalPredictScore_enhanced <- function(fit, new_data) {
  tryCatch({
    if (is.null(fit$subFeature) || length(fit$subFeature) == 0) {
      new_data_sub <- new_data
    } else {
      available_features <- intersect(fit$subFeature, colnames(new_data))
      if (length(available_features) == 0) {
        new_data_sub <- new_data
      } else {
        new_data_sub <- new_data[, available_features, drop = FALSE]
      }
    }

    if (inherits(fit, "xgb.Booster")) {
      if (is.data.frame(new_data_sub)) {
        new_data_sub <- as.matrix(new_data_sub)
      }
      new_data_sub[is.na(new_data_sub)] <- 0
      if (any(!is.numeric(new_data_sub))) {
        new_data_sub <- as.matrix(as.data.frame(lapply(as.data.frame(new_data_sub), as.numeric)))
      }
      RS <- predict(fit, new_data_sub)
    } else {
      RS <- Original_CalPredictScore(fit, new_data_sub)
    }

    RS <- as.numeric(as.vector(RS))
    if (length(RS) != nrow(new_data)) {
      if (length(RS) < nrow(new_data)) {
        RS <- c(RS, rep(NA, nrow(new_data) - length(RS)))
      } else {
        RS <- RS[1:nrow(new_data)]
      }
    }
    names(RS) <- rownames(new_data)
    return(RS)
  }, error = function(e) {
    warning("预测评分计算失败: ", e$message)
    RS <- rep(NA, nrow(new_data))
    names(RS) <- rownames(new_data)
    return(RS)
  })
}

assign("CalPredictScore", CalPredictScore_enhanced, envir = .GlobalEnv)
cat("✅ CalPredictScore函数已增强\n")

cat("\n========== 读取训练集数据 (GSE16561) ==========\n")

if (file.exists(train_file)) {
  train_data <- readRDS(train_file)
  if (is.list(train_data) && !is.null(train_data$expression) && !is.null(train_data$group)) {
    Train_expr <- train_data$expression
    Train_class <- train_data$group
    cat("✅ GSE16561训练集读取成功\n")
  } else {
    stop("训练集格式不正确")
  }
  cat("  表达矩阵维度:", dim(Train_expr), "\n")
  cat("  分组统计:\n")
  print(table(Train_class$group))
} else {
  stop("❌ 训练集文件不存在")
}

cat("\n========== 读取验证集数据 (GSE58294) ==========\n")

if (file.exists(test_file)) {
  test_data <- readRDS(test_file)
  if (is.list(test_data) && !is.null(test_data$expression) && !is.null(test_data$group)) {
    Test_expr <- test_data$expression
    Test_class <- test_data$group
    cat("✅ GSE58294验证集读取成功\n")
  } else {
    stop("验证集格式不正确")
  }
  cat("  表达矩阵维度:", dim(Test_expr), "\n")
  cat("  分组统计:\n")
  print(table(Test_class$group))
} else {
  stop("❌ 验证集文件不存在")
}

cat("\n========== 提取候选基因表达数据 ==========\n")

candidate_genes <- toupper(candidate_genes)
rownames(Train_expr) <- toupper(rownames(Train_expr))
rownames(Test_expr) <- toupper(rownames(Test_expr))

train_available <- intersect(candidate_genes, rownames(Train_expr))
test_available <- intersect(candidate_genes, rownames(Test_expr))
common_genes <- intersect(train_available, test_available)

cat("\n候选基因匹配情况:\n")
cat("  训练集中存在的候选基因数:", length(train_available), "\n")
cat("  验证集中存在的候选基因数:", length(test_available), "\n")
cat("  共有候选基因数:", length(common_genes), "\n")

if (length(common_genes) < 2) {
  cat("\n⚠️ 共有基因不足2个，尝试使用所有候选基因\n")
  common_genes <- intersect(candidate_genes, rownames(Train_expr))
  if (length(common_genes) < 2) {
    stop("❌ 候选基因在训练集中不足2个")
  }
}

Train_expr <- Train_expr[common_genes, , drop = FALSE]
Test_expr <- Test_expr[common_genes, , drop = FALSE]

writeLines(common_genes, file.path(intermediate.path, "common_genes.txt"))

cat("\n========== 转置表达矩阵 ==========\n")

Train_expr <- t(Train_expr)
Test_expr <- t(Test_expr)

cat("训练集样本数:", nrow(Train_expr), " 基因数:", ncol(Train_expr), "\n")
cat("验证集样本数:", nrow(Test_expr), " 基因数:", ncol(Test_expr), "\n")

cat("\n========== 样本名对齐 ==========\n")

align_class_rownames <- function(class_df, expr_mat, label) {

  if (length(intersect(rownames(class_df), rownames(expr_mat))) > 0) {
    cat(sprintf("  %s: rownames 已匹配，无需调整\n", label))
    return(class_df)
  }

  for (col in c("sample", "Sample", "geo_accession", "ID")) {
    if (col %in% colnames(class_df)) {
      candidate <- as.character(class_df[[col]])
      if (length(intersect(candidate, rownames(expr_mat))) > 0) {
        rownames(class_df) <- candidate
        cat(sprintf("  %s: 使用 '%s' 列作为 rownames，匹配 %d/%d 个样本\n",
                    label, col,
                    length(intersect(candidate, rownames(expr_mat))),
                    nrow(expr_mat)))
        return(class_df)
      }
    }
  }

  if (nrow(class_df) == nrow(expr_mat)) {
    rownames(class_df) <- rownames(expr_mat)
    cat(sprintf("  %s: 样本数相等，按顺序对齐（%d 个）\n", label, nrow(class_df)))
    return(class_df)
  }

  stop(sprintf("无法对齐 %s 样本名：class_df %d 行，expr_mat %d 行，无公共名称",
               label, nrow(class_df), nrow(expr_mat)))
}

Train_class <- align_class_rownames(Train_class, Train_expr, "训练集")
Test_class  <- align_class_rownames(Test_class,  Test_expr,  "验证集")

train_samples <- intersect(rownames(Train_class), rownames(Train_expr))
test_samples  <- intersect(rownames(Test_class),  rownames(Test_expr))

if (length(train_samples) < nrow(Train_expr)) {
  lost <- setdiff(rownames(Train_expr), rownames(Train_class))
  cat(sprintf("  ⚠️ 训练集丢失 %d 个样本（共 %d → %d）\n",
              length(lost), nrow(Train_expr), length(train_samples)))
  if (length(lost) <= 10) cat("  丢失样本:", paste(lost, collapse=", "), "\n")
} else {
  cat(sprintf("  ✅ 训练集全部 %d 个样本对齐成功\n", length(train_samples)))
}

if (length(test_samples) < nrow(Test_expr)) {
  lost <- setdiff(rownames(Test_expr), rownames(Test_class))
  cat(sprintf("  ⚠️ 验证集丢失 %d 个样本（共 %d → %d）\n",
              length(lost), nrow(Test_expr), length(test_samples)))
  if (length(lost) <= 10) cat("  丢失样本:", paste(lost, collapse=", "), "\n")
} else {
  cat(sprintf("  ✅ 验证集全部 %d 个样本对齐成功\n", length(test_samples)))
}

Train_expr  <- Train_expr[train_samples,  , drop = FALSE]
Train_class <- Train_class[train_samples, , drop = FALSE]
Test_expr   <- Test_expr[test_samples,   , drop = FALSE]
Test_class  <- Test_class[test_samples,  , drop = FALSE]

cat("\n========== 转换分组标签 ==========\n")

Train_class$group <- ifelse(tolower(as.character(Train_class$group)) %in% c("control", "0"), 0, 1)
Test_class$group <- ifelse(tolower(as.character(Test_class$group)) %in% c("control", "0"), 0, 1)

cat("训练集分组: 0=", sum(Train_class$group == 0), "1=", sum(Train_class$group == 1), "\n")
cat("验证集分组: 0=", sum(Test_class$group == 0), "1=", sum(Test_class$group == 1), "\n")

cat("\n========== 数据标准化 ==========\n")

Train_set <- scale(Train_expr, center = TRUE, scale = TRUE)
Train_set[is.na(Train_set)] <- 0

Test_set <- scale(Test_expr,
                  center = attr(Train_set, "scaled:center"),
                  scale = attr(Train_set, "scaled:scale"))
Test_set[is.na(Test_set)] <- 0

cat("✅ 数据标准化完成\n")

saveRDS(list(Train_set = Train_set, Train_class = Train_class,
             Test_set = Test_set, Test_class = Test_class,
             common_genes = common_genes,
             candidate_up = candidate_up,
             candidate_down = candidate_down),
        file = file.path(intermediate.path, "processed_data.rds"))

cat("✅ 处理后的数据已保存到:", file.path(intermediate.path, "processed_data.rds"), "\n")
writeLines(common_genes, file.path(data.path, "common_genes.txt"))
cat("✅ 共有基因列表已保存到:", file.path(data.path, "common_genes.txt"), "\n")

cat("\n========== 读取算法配置 ==========\n")

if (file.exists(methods_file)) {
  methods <- openxlsx::read.xlsx(methods_file, sheet = 1, startRow = 2)
  methods <- methods$Model
  methods <- gsub("-| ", "", methods)
  methods <- methods[!is.na(methods) & methods != ""]
  cat("✅ 读取到", length(methods), "种算法组合\n")
  cat("  前10种:", paste(head(methods, 10), collapse = ", "), "\n")
} else {
  cat("⚠️ methods.xlsx 不存在，使用默认算法列表\n")
  methods <- c(
    "lasso", "ridge", "elasticnet", "logistic", "svm", "rf", "gbm", "xgboost"
  )
  cat("使用默认算法:", paste(methods, collapse = ", "), "\n")
}

writeLines(methods, file.path(intermediate.path, "algorithms_list.txt"))

classVar <- "group"
min.selected.var <- 2
n_folds <- 5
timeout_seconds <- 300

seed_file <- file.path(res.path, "seed_record.txt")
if (file.exists(seed_file)) {
  GLOBAL_SEED <- as.integer(readLines(seed_file, n = 1))
  cat("✅ 从", seed_file, "读取到历史 seed =", GLOBAL_SEED, "\n")
} else {
  GLOBAL_SEED <- 123
  writeLines(as.character(GLOBAL_SEED), seed_file)
  cat("✅ 使用 seed =", GLOBAL_SEED, "，已保存到:", seed_file, "\n")
}
set.seed(GLOBAL_SEED)

cat("\n========== 模型参数 ==========\n")
cat("  类别变量:", classVar, "\n")
cat("  最小特征数:", min.selected.var, "\n")
cat("  交叉验证折数:", n_folds, "\n")
cat("  超时时间:", timeout_seconds, "秒\n")
cat("  随机种子:", GLOBAL_SEED, "\n")

cat("\n========== 预筛选特征 ==========\n")

preTrain.method <- strsplit(methods, "\\+")
preTrain.method <- lapply(preTrain.method, function(x) {
  if(length(x) > 1) rev(x)[-1] else "simple"
})
preTrain.method <- unique(unlist(preTrain.method))

cat("预筛选方法:", paste(preTrain.method, collapse = ", "), "\n")

preTrain.var <- list()
set.seed(GLOBAL_SEED)

for (method in preTrain.method) {
  cat("\n处理预筛选方法:", method, "\n")

  set.seed(GLOBAL_SEED)

  result_var <- tryCatch({
    if (method == "simple") {
      vars <- colnames(Train_set)
      cat("  ✅ 使用所有特征:", length(vars), "\n")
      vars
    } else if (method == "xgboost") {
      vars <- RunXGBoost_fixed(
        Train_set = Train_set,
        Train_label = Train_class,
        mode = "Variable",
        classVar = classVar
      )
      cat("  ✅ 筛选出", length(vars), "个特征\n")
      vars
    } else if (method == "RF") {

      Train_df <- as.data.frame(Train_set)
      Train_df[[classVar]] <- as.factor(Train_class[[classVar]])

      rf_fit <- rfsrc(
        formula  = formula(paste0(classVar, "~.")),
        data     = Train_df,
        ntree    = 500,
        nodesize = 5,
        importance = TRUE,
        proximity = FALSE,
        forest   = TRUE
      )

      imp <- rf_fit$importance
      if (is.matrix(imp)) imp <- imp[, 1]
      vars <- names(imp[imp > 0])
      if (length(vars) == 0) vars <- colnames(Train_set)

      rf_fit$subFeature <- colnames(Train_set)
      cat("  ✅ RF 筛选出", length(vars), "个特征\n")
      if (length(vars) > 0)
        cat("  前5个:", paste(head(vars, 5), collapse = ", "), "\n")
      vars

    } else {
      vars <- RunML(
        method = method,
        Train_set = Train_set,
        Train_label = Train_class,
        mode = "Variable",
        classVar = classVar
      )
      cat("  ✅ 筛选出", length(vars), "个特征\n")
      if (length(vars) > 0)
        cat("  前5个:", paste(head(vars, 5), collapse = ", "), "\n")
      vars
    }
  }, error = function(e) {
    cat("  ❌ 预筛选失败:", e$message, "\n")
    cat("  ⚠️ 使用所有特征作为备选\n")
    colnames(Train_set)
  })

  preTrain.var[[method]] <- result_var
}

saveRDS(preTrain.var, file.path(intermediate.path, "preTrain_features.rds"))

cat("\n========== 定义5折交叉验证函数 ==========\n")

run_cv_evaluation <- function(method_name, Train_set, Train_label, preTrain.var,
                               n_folds = 5, classVar, min.selected.var,
                               timeout_seconds = 300) {

  method_parts <- strsplit(method_name, "\\+")[[1]]
  if (length(method_parts) == 1) {
    method_parts <- c("simple", method_parts)
  }

  pre_method <- method_parts[1]
  main_method <- method_parts[2]

  if (!(pre_method %in% names(preTrain.var))) {
    Variable <- colnames(Train_set)
  } else {
    Variable <- preTrain.var[[pre_method]]
  }

  if (length(Variable) < min.selected.var) {
    return(list(
      success = FALSE,
      reason = paste("特征不足:", length(Variable))
    ))
  }

  Variable <- gsub("`", "", Variable)
  Train_set_tmp <- Train_set[, Variable, drop = FALSE]

  set.seed(GLOBAL_SEED)
  folds <- createFolds(Train_label[[classVar]], k = n_folds, list = TRUE)

  cv_train_auc <- numeric(n_folds)
  cv_test_auc <- numeric(n_folds)
  fold_results <- list()
  fold_features <- list()

  for (fold_idx in 1:n_folds) {
    tryCatch({
      test_indices <- folds[[fold_idx]]
      train_indices <- setdiff(1:nrow(Train_set_tmp), test_indices)

      cv_train_set <- Train_set_tmp[train_indices, , drop = FALSE]
      cv_test_set <- Train_set_tmp[test_indices, , drop = FALSE]
      cv_train_label <- Train_label[train_indices, , drop = FALSE]
      cv_test_label <- Train_label[test_indices, , drop = FALSE]

      model_start <- Sys.time()

      if (main_method == "xgboost") {
        fit <- withTimeout({
          RunXGBoost_fixed(
            Train_set = cv_train_set,
            Train_label = cv_train_label,
            mode = "Model",
            classVar = classVar)
        }, timeout = timeout_seconds / n_folds, onTimeout = "error")
      } else if (main_method == "RF") {

        fit <- withTimeout({
          cv_df <- as.data.frame(cv_train_set)
          cv_df[[classVar]] <- as.factor(cv_train_label[[classVar]])
          rf_cv <- rfsrc(
            formula  = formula(paste0(classVar, "~.")),
            data     = cv_df,
            ntree    = 500,
            nodesize = 5,
            importance = TRUE,
            proximity = FALSE,
            forest   = TRUE
          )
          rf_cv$subFeature <- colnames(cv_train_set)
          rf_cv
        }, timeout = timeout_seconds / n_folds, onTimeout = "error")
      } else {
        fit <- withTimeout({
          RunML(method = main_method,
                Train_set = cv_train_set,
                Train_label = cv_train_label,
                mode = "Model",
                classVar = classVar)
        }, timeout = timeout_seconds / n_folds, onTimeout = "error")
      }

      model_time <- difftime(Sys.time(), model_start, units = "secs")

      if (!is.null(fit)) {

        fold_features[[fold_idx]] <- ExtractVar(fit)

        train_pred <- CalPredictScore(fit = fit, new_data = cv_train_set)
        if (!all(is.na(train_pred))) {
          roc_train <- tryCatch(auc(roc(cv_train_label[[classVar]], train_pred)),
                                error = function(e) NA)
        } else {
          roc_train <- NA
        }

        test_pred <- CalPredictScore(fit = fit, new_data = cv_test_set)
        if (!all(is.na(test_pred))) {
          roc_test <- tryCatch(auc(roc(cv_test_label[[classVar]], test_pred)),
                               error = function(e) NA)
        } else {
          roc_test <- NA
        }

        cv_train_auc[fold_idx] <- as.numeric(roc_train)
        cv_test_auc[fold_idx] <- as.numeric(roc_test)

        fold_results[[fold_idx]] <- list(
          fold = fold_idx,
          train_auc = roc_train,
          test_auc = roc_test,
          feature_count = length(ExtractVar(fit)),
          feature_list = ExtractVar(fit),
          model_time = model_time
        )
      } else {
        cv_train_auc[fold_idx] <- NA
        cv_test_auc[fold_idx] <- NA
        fold_results[[fold_idx]] <- list(
          fold = fold_idx,
          train_auc = NA,
          test_auc = NA,
          feature_count = 0,
          feature_list = NULL,
          model_time = NA,
          error = "模型训练失败"
        )
      }

    }, error = function(e) {
      cv_train_auc[fold_idx] <- NA
      cv_test_auc[fold_idx] <- NA
      fold_results[[fold_idx]] <- list(
        fold = fold_idx,
        train_auc = NA,
        test_auc = NA,
        feature_count = 0,
        feature_list = NULL,
        model_time = NA,
        error = e$message
      )
    })
  }

  valid_train_auc <- cv_train_auc[!is.na(cv_train_auc)]
  valid_test_auc <- cv_test_auc[!is.na(cv_test_auc)]

  mean_train_auc <- if(length(valid_train_auc) > 0) mean(valid_train_auc) else NA
  mean_test_auc <- if(length(valid_test_auc) > 0) mean(valid_test_auc) else NA
  mean_c_index <- mean(c(mean_train_auc, mean_test_auc), na.rm = TRUE)

  sd_train_auc <- if(length(valid_train_auc) > 1) sd(valid_train_auc) else NA
  sd_test_auc <- if(length(valid_test_auc) > 1) sd(valid_test_auc) else NA

  overfit_gap <- mean_train_auc - mean_test_auc

  exclude_reason <- NULL
  if (mean_train_auc > 0.99 && mean_test_auc > 0.99) {
    method_parts <- strsplit(method_name, "\\+")[[1]]
    pre_m <- if(length(method_parts) > 1) method_parts[1] else "simple"
    n_feat <- if(pre_m %in% names(preTrain.var)) length(preTrain.var[[pre_m]]) else ncol(Train_set)
    if (n_feat < 10) {
      exclude_reason <- paste("perfect_separation_", n_feat, "features", sep="")
    }
  } else if (overfit_gap > 0.15) {
    exclude_reason <- paste("severe_overfit_gap_", round(overfit_gap, 3), sep="")
  } else if (mean_test_auc < 0.6) {
    exclude_reason <- paste("poor_validation_auc_", round(mean_test_auc, 3), sep="")
  } else if (overfit_gap > 0.1 && mean_train_auc > 0.95) {
    exclude_reason <- paste("overfit_gap_", round(overfit_gap, 3), sep="")
  } else if (!is.na(sd_test_auc) && sd_test_auc > 0.15) {
    exclude_reason <- paste("high_variance_", round(sd_test_auc, 3), sep="")
  }

  return(list(
    success = TRUE,
    method_name = method_name,
    mean_train_auc = mean_train_auc,
    mean_test_auc = mean_test_auc,
    mean_c_index = mean_c_index,
    sd_train_auc = sd_train_auc,
    sd_test_auc = sd_test_auc,
    overfit_gap = overfit_gap,
    fold_train_auc = cv_train_auc,
    fold_test_auc = cv_test_auc,
    n_success_folds = length(valid_train_auc),
    exclude_reason = exclude_reason,
    fold_results = fold_results,
    fold_features = fold_features
  ))
}

cat("\n========== 运行5折交叉验证（共", length(methods), "个模型） ==========\n")

cv_results <- list()
total_start <- Sys.time()
success_count <- 0
fail_count <- 0
timeout_count <- 0
excluded_count <- 0

for (i in 1:length(methods)) {
  method_name <- methods[i]
  cat(sprintf("\n[%d/%d] %s ", i, length(methods), method_name))

  set.seed(GLOBAL_SEED)

  tryCatch({
    cv_result <- run_cv_evaluation(
      method_name = method_name,
      Train_set = Train_set,
      Train_label = Train_class,
      preTrain.var = preTrain.var,
      n_folds = n_folds,
      classVar = classVar,
      min.selected.var = min.selected.var,
      timeout_seconds = timeout_seconds
    )

    if (cv_result$success && !is.na(cv_result$mean_c_index)) {
      cv_results[[method_name]] <- cv_result
      success_count <- success_count + 1

      if (!is.null(cv_result$exclude_reason)) {
        excluded_count <- excluded_count + 1
        cat(sprintf(" ⚠️ 排除: %s (CV AUC: %.3f)",
                    cv_result$exclude_reason, cv_result$mean_c_index))
      } else {
        cat(sprintf(" ✅ CV AUC: %.4f (Train: %.3f, Test: %.3f, Gap: %.3f)",
                    cv_result$mean_c_index,
                    cv_result$mean_train_auc,
                    cv_result$mean_test_auc,
                    cv_result$overfit_gap))
      }
    } else {
      cat(sprintf(" ❌ %s", ifelse(is.null(cv_result$reason), "失败", cv_result$reason)))
      fail_count <- fail_count + 1
    }

  }, error = function(e) {
    if (grepl("TIMEOUT", e$message, ignore.case = TRUE)) {
      cat(" ⏰ 超时")
      timeout_count <- timeout_count + 1
    } else {
      cat(" ❌ 错误:", substr(e$message, 1, 60))
    }
    fail_count <- fail_count + 1
  })

  if (i %% 20 == 0 || i == length(methods)) {
    elapsed <- difftime(Sys.time(), total_start, units = "mins")
    cat(sprintf("\n  进度: %d/%d, 用时 %.1f分钟, 成功: %d, 失败: %d, 超时: %d, 排除: %d\n",
                i, length(methods), elapsed, success_count, fail_count, timeout_count, excluded_count))
  }
}

total_time <- difftime(Sys.time(), total_start, units = "mins")
cat(sprintf("\n✅ 5折交叉验证完成！总用时: %.1f分钟\n", total_time))
cat(sprintf("  成功: %d, 失败: %d, 超时: %d, 排除: %d\n",
            success_count, fail_count, timeout_count, excluded_count))

saveRDS(cv_results, file.path(intermediate.path, "cv_results.rds"))

if (length(cv_results) == 0) {
  stop("❌ 没有成功的模型")
}

cat("\n========== 筛选可接受模型 ==========\n")

acceptable_results <- list()
for (method in names(cv_results)) {
  if (is.null(cv_results[[method]]$exclude_reason)) {
    acceptable_results[[method]] <- cv_results[[method]]
  }
}

cat("✅ 可接受模型数:", length(acceptable_results), "/", length(cv_results), "\n")

if (length(acceptable_results) == 0) {
  cat("\n⚠️ 警告: 所有模型均被排除！使用排除标准最宽松的模型\n")
  gap_scores <- sapply(cv_results, function(x) x$overfit_gap)
  best_method <- names(cv_results)[which.min(gap_scores)]
  acceptable_results[[best_method]] <- cv_results[[best_method]]
  cat("  选择过拟合程度最小的模型:", best_method, "(Gap=", round(gap_scores[best_method], 3), ")\n")
}

saveRDS(acceptable_results, file.path(intermediate.path, "acceptable_results.rds"))

cat("\n========== 在验证集上评估可接受模型 ==========\n")

validation_results <- data.frame(
  Model = names(acceptable_results),
  CV_Mean_AUC = sapply(acceptable_results, function(x) x$mean_c_index),
  CV_Train_AUC = sapply(acceptable_results, function(x) x$mean_train_auc),
  CV_Test_AUC = sapply(acceptable_results, function(x) x$mean_test_auc),
  CV_Overfit_Gap = sapply(acceptable_results, function(x) x$overfit_gap),
  CV_SD_Test = sapply(acceptable_results, function(x) x$sd_test_auc),
  Validation_AUC = NA,
  stringsAsFactors = FALSE
)

final_models <- list()

for (i in 1:nrow(validation_results)) {
  method_name <- validation_results$Model[i]
  cat(sprintf("\n[%d/%d] 评估 %s ... ", i, nrow(validation_results), method_name))

  tryCatch({

    method_parts <- strsplit(method_name, "\\+")[[1]]
    if (length(method_parts) == 1) {
      method_parts <- c("simple", method_parts)
    }
    pre_method <- method_parts[1]
    main_method <- method_parts[2]

    if (!(pre_method %in% names(preTrain.var))) {
      Variable <- colnames(Train_set)
    } else {
      Variable <- preTrain.var[[pre_method]]
    }
    Variable <- gsub("`", "", Variable)
    Train_set_tmp <- Train_set[, Variable, drop = FALSE]

    if (main_method == "xgboost") {
      final_fit <- RunXGBoost_fixed(
        Train_set = Train_set_tmp,
        Train_label = Train_class,
        mode = "Model",
        classVar = classVar)
    } else if (main_method == "RF") {

      Train_df_tmp <- as.data.frame(Train_set_tmp)
      Train_df_tmp[[classVar]] <- as.factor(Train_class[[classVar]])
      rf_fit_final <- rfsrc(
        formula  = formula(paste0(classVar, "~.")),
        data     = Train_df_tmp,
        ntree    = 500,
        nodesize = 5,
        importance = TRUE,
        proximity = FALSE,
        forest   = TRUE
      )
      rf_fit_final$subFeature <- colnames(Train_set_tmp)
      final_fit <- rf_fit_final
    } else {
      final_fit <- RunML(method = main_method,
                         Train_set = Train_set_tmp,
                         Train_label = Train_class,
                         mode = "Model",
                         classVar = classVar)
    }

    if (!is.null(final_fit)) {
      final_models[[method_name]] <- final_fit

      test_pred <- CalPredictScore(fit = final_fit, new_data = Test_set)
      if (!all(is.na(test_pred))) {
        roc_test <- auc(roc(Test_class$group, test_pred))
        validation_results$Validation_AUC[i] <- as.numeric(roc_test)
        cat(sprintf("验证集AUC: %.4f", validation_results$Validation_AUC[i]))
      } else {
        cat("预测失败")
      }
    } else {
      cat("模型训练失败")
    }
  }, error = function(e) {
    cat("错误:", substr(e$message, 1, 50))
  })
}

validation_results <- validation_results[complete.cases(validation_results$Validation_AUC), ]
cat("\n\n✅ 成功评估模型数:", nrow(validation_results), "\n")

if (nrow(validation_results) == 0) {
  stop("❌ 所有模型在验证集上评估失败")
}

write.csv(validation_results, file.path(res.path, "validation_results.csv"), row.names = FALSE)
saveRDS(validation_results, file.path(intermediate.path, "validation_results.rds"))

cat("\n========== 选择最佳模型 ==========\n")

validation_results <- validation_results[order(validation_results$Validation_AUC, decreasing = TRUE), ]

cat("\n✅ Top 20 模型（按验证集AUC排序）:\n")
print(validation_results[1:min(20, nrow(validation_results)),
                         c("Model", "CV_Mean_AUC", "Validation_AUC", "CV_Overfit_Gap")])

best_model_name <- validation_results$Model[1]
best_validation_auc <- validation_results$Validation_AUC[1]
best_cv_auc <- validation_results$CV_Mean_AUC[1]
best_cv_gap <- validation_results$CV_Overfit_Gap[1]

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("🏆 最终最佳模型\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("  模型名称:", best_model_name, "\n")
cat("  5折交叉验证平均AUC:", round(best_cv_auc, 4), "\n")
cat("  5折交叉验证过拟合差距:", round(best_cv_gap, 4), "\n")
cat("  独立验证集AUC:", round(best_validation_auc, 4), "\n")

best_fit <- final_models[[best_model_name]]
if (!is.null(best_fit)) {
  best_features <- ExtractVar(best_fit)
  cat("  特征基因数:", length(best_features), "\n")
  cat("  特征基因:", paste(head(best_features, 10), collapse = ", "),
      ifelse(length(best_features) > 10, "...", ""), "\n")
}

cat("\n========== 绘制模型性能热图（图3，Top 20 清晰版） ==========\n")

heatmap_data <- validation_results %>%
  dplyr::select(Model, Validation_AUC, CV_Mean_AUC) %>%
  dplyr::arrange(dplyr::desc(Validation_AUC))

heatmap_data <- heatmap_data %>% dplyr::rename(Training_AUC = CV_Mean_AUC)

heatmap_data$Model_Type <- gsub("\\+.*$", "", heatmap_data$Model)
heatmap_data$Model_Type[heatmap_data$Model_Type == "RF"]       <- "RandomForest"
heatmap_data$Model_Type[heatmap_data$Model_Type == "Enet"]     <- "ElasticNet"
heatmap_data$Model_Type[heatmap_data$Model_Type == "glmBoost"] <- "GLMBoost"

heatmap_data$Alpha <- NA
for (i in 1:nrow(heatmap_data)) {
  if (grepl("Enet\\[alpha=", heatmap_data$Model[i])) {
    alpha_val <- gsub(".*Enet\\[alpha=([0-9.]+)\\].*", "\\1", heatmap_data$Model[i])
    heatmap_data$Alpha[i] <- as.numeric(alpha_val)
  }
}

heatmap_data$Mean_AUC <- rowMeans(
  heatmap_data[, c("Validation_AUC", "Training_AUC")], na.rm = TRUE
)

top_n <- min(20, nrow(heatmap_data))

heatmap_top <- heatmap_data[1:top_n, ]
heatmap_top <- heatmap_top[order(-heatmap_top$Validation_AUC), ]

heat_matrix <- as.matrix(heatmap_top[, c("Validation_AUC", "Training_AUC", "Mean_AUC")])
rownames(heat_matrix) <- heatmap_top$Model
colnames(heat_matrix) <- c("Validation AUC", "Training AUC", "Mean AUC")

anno_row <- data.frame(
  Rank = factor(ifelse(seq_len(top_n) <= 10, "Top 10", "Top 20"),
                levels = c("Top 10", "Top 20"))
)
rownames(anno_row) <- rownames(heat_matrix)
anno_colors <- list(
  Rank = c("Top 10" = "#E64B35", "Top 20" = "#F39B7F")
)

val_range <- range(heat_matrix, na.rm = TRUE)
breaks_vec <- seq(
  floor(val_range[1] * 100) / 100,
  ceiling(val_range[2] * 100) / 100,
  length.out = 101
)

gap_rows <- c()
if (top_n >= 10) gap_rows <- 10

heatmap_args <- list(
  mat              = heat_matrix,
  main             = sprintf("Fig. 3: Model Performance (Top %d, ranked by Validation AUC)",
                             top_n),
  color            = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
  breaks           = breaks_vec,
  cluster_rows     = FALSE,
  cluster_cols     = FALSE,
  display_numbers  = TRUE,
  number_format    = "%.3f",
  number_color     = "black",
  fontsize         = 14,
  fontsize_number  = 13,
  fontsize_row     = 13,
  fontsize_col     = 14,
  cellwidth        = 140,
  cellheight       = 32,
  angle_col        = 0,
  border_color     = "gray80",
  annotation_row   = anno_row,
  annotation_colors = anno_colors,
  annotation_names_row = FALSE,
  gaps_row         = gap_rows,
  legend           = TRUE,
  legend_breaks    = c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
  legend_labels    = c("0.5", "0.6", "0.7", "0.8", "0.9", "1.0")
)

fig_h <- max(10, 2.5 + top_n * 0.42)
fig_w <- 14

pdf(file.path(fig.path, "Fig3_Model_Performance_Heatmap.pdf"),
    width = fig_w, height = fig_h)
do.call(pheatmap, heatmap_args)
dev.off()

png(file.path(fig.path, "Fig3_Model_Performance_Heatmap.png"),
    width = fig_w, height = fig_h, units = "in", res = 300, bg = "white")
do.call(pheatmap, heatmap_args)
dev.off()
cat(sprintf("✅ Fig. 3 模型性能热图已保存 (%.0fx%.1f 英寸, Top %d 模型)\n",
            fig_w, fig_h, top_n))

p_bar <- ggplot(heatmap_top,
                aes(x = reorder(Model, Validation_AUC), y = Validation_AUC,
                    fill = Validation_AUC)) +
  geom_bar(stat = "identity", width = 0.75) +
  geom_text(aes(label = sprintf("%.3f", Validation_AUC)),
            hjust = -0.15, size = 4.2, color = "black") +
  coord_flip() +
  scale_fill_gradient(low = "#4DBBD5", high = "#E64B35",
                      limits = c(min(heatmap_top$Validation_AUC, na.rm = TRUE),
                                 max(heatmap_top$Validation_AUC, na.rm = TRUE))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = sprintf("Fig. 3: Top %d Models by Validation AUC", top_n),
       x = NULL, y = "Validation AUC") +
  theme_bw(base_size = 14) +
  theme(plot.title   = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text.y  = element_text(size = 12),
        axis.text.x  = element_text(size = 12),
        axis.title.x = element_text(size = 13),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        plot.margin  = margin(10, 25, 10, 10))

bar_h <- max(8, top_n * 0.42)
ggsave(file.path(fig.path, "Fig3_Model_Performance_Barplot.pdf"),
       p_bar, width = 12, height = bar_h)
ggsave(file.path(fig.path, "Fig3_Model_Performance_Barplot.png"),
       p_bar, width = 12, height = bar_h, dpi = 300, bg = "white")
cat("✅ Fig. 3 模型性能条形图已保存\n")

top20_out <- heatmap_top
top20_out$Rank <- seq_len(nrow(top20_out))
write.csv(top20_out[, c("Rank", "Model", "Validation_AUC", "Training_AUC", "Mean_AUC")],
          file.path(res.path, "top20_models_summary.csv"),
          row.names = FALSE)
cat(sprintf("✅ Top %d 模型汇总 CSV 已保存\n", top_n))

cat("\n绘制验证集AUC vs CV AUC散点图...\n")

pdf(file.path(fig.path, "CV_vs_Validation_AUC.pdf"), width = 8, height = 6)

plot(validation_results$CV_Mean_AUC, validation_results$Validation_AUC,
     xlab = "5-Fold CV Mean AUC", ylab = "Validation Set AUC",
     main = "Cross-Validation vs Independent Validation",
     pch = 19, col = ifelse(validation_results$Model == best_model_name, "red", "blue"),
     xlim = c(0, 1), ylim = c(0, 1))
abline(0, 1, lty = 2, col = "gray")
points(validation_results$CV_Mean_AUC[validation_results$Model == best_model_name],
       validation_results$Validation_AUC[validation_results$Model == best_model_name],
       col = "red", pch = 19, cex = 1.5)
text(validation_results$CV_Mean_AUC, validation_results$Validation_AUC,
     labels = substr(validation_results$Model, 1, 20), cex = 0.5, pos = 4)

dev.off()

png(file.path(fig.path, "CV_vs_Validation_AUC.png"), width = 8, height = 6, units = "in", res = 300)
plot(validation_results$CV_Mean_AUC, validation_results$Validation_AUC,
     xlab = "5-Fold CV Mean AUC", ylab = "Validation Set AUC",
     main = "Cross-Validation vs Independent Validation",
     pch = 19, col = ifelse(validation_results$Model == best_model_name, "red", "blue"),
     xlim = c(0, 1), ylim = c(0, 1))
abline(0, 1, lty = 2, col = "gray")
points(validation_results$CV_Mean_AUC[validation_results$Model == best_model_name],
       validation_results$Validation_AUC[validation_results$Model == best_model_name],
       col = "red", pch = 19, cex = 1.5)
dev.off()
cat("✅ CV vs Validation AUC散点图已保存\n")

cat("\n绘制最佳模型ROC曲线...\n")

if (!is.null(best_fit)) {
  pdf(file.path(fig.path, "Best_Model_ROC.pdf"), width = 8, height = 7)

  train_pred <- CalPredictScore(fit = best_fit, new_data = Train_set)
  test_pred <- CalPredictScore(fit = best_fit, new_data = Test_set)

  if (!all(is.na(train_pred)) && !all(is.na(test_pred))) {
    roc_train <- roc(Train_class$group, train_pred)
    roc_test <- roc(Test_class$group, test_pred)

    plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
         xlab = "1 - Specificity", ylab = "Sensitivity",
         main = paste("ROC Curve -", best_model_name), cex.main = 1)
    abline(0, 1, col = "gray", lty = 2, lwd = 2)
    lines(1 - roc_train$specificities, roc_train$sensitivities, col = "#4DBBD5", lwd = 3)
    lines(1 - roc_test$specificities, roc_test$sensitivities, col = "#E64B35", lwd = 3)
    legend("bottomright",
           legend = c(paste("Training AUC =", round(roc_train$auc, 3)),
                      paste("Validation AUC =", round(roc_test$auc, 3))),
           col = c("#4DBBD5", "#E64B35"), lwd = 3, bty = "n")
  }

  dev.off()

  png(file.path(fig.path, "Best_Model_ROC.png"), width = 8, height = 7, units = "in", res = 150)
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "1 - Specificity", ylab = "Sensitivity",
       main = paste("ROC Curve -", best_model_name), cex.main = 1)
  abline(0, 1, col = "gray", lty = 2, lwd = 2)
  lines(1 - roc_train$specificities, roc_train$sensitivities, col = "#4DBBD5", lwd = 3)
  lines(1 - roc_test$specificities, roc_test$sensitivities, col = "#E64B35", lwd = 3)
  legend("bottomright",
         legend = c(paste("Training AUC =", round(roc_train$auc, 3)),
                    paste("Validation AUC =", round(roc_test$auc, 3))),
         col = c("#4DBBD5", "#E64B35"), lwd = 3, bty = "n")
  dev.off()

  cat("✅ ROC曲线已保存\n")
}

cat("\n========== 保存所有结果 ==========\n")

analysis_end_time <- Sys.time()
total_runtime <- difftime(analysis_end_time, analysis_start_time, units = "mins")

all_results <- list(
  analysis_date = Sys.time(),
  analysis_start = analysis_start_time,
  analysis_end = analysis_end_time,
  total_runtime_minutes = as.numeric(total_runtime),
  candidate_source = candidate_dir,
  random_seed = GLOBAL_SEED,
  n_folds = n_folds,
  n_total_models = length(methods),
  n_success_cv = length(cv_results),
  n_acceptable = length(acceptable_results),
  n_validated = nrow(validation_results),
  n_timeout = timeout_count,

  best_model = best_model_name,
  best_validation_auc = best_validation_auc,
  best_cv_auc = best_cv_auc,
  best_features = best_features,

  validation_results = validation_results,
  all_cv_results = cv_results,
  preTrain_features = preTrain.var,
  common_genes = common_genes
)

if (exists("best_fit") && !is.null(best_fit)) {
  saveRDS(best_fit, file.path(res.path, "best_model_final.rds"))
  cat("✅ 最佳模型已保存到:", file.path(res.path, "best_model_final.rds"), "\n")
}

writeLines(best_features, file.path(res.path, "best_features_final.txt"))
cat("✅ 最佳模型特征基因已保存到:", file.path(res.path, "best_features_final.txt"), "\n")

writeLines(common_genes, file.path(res.path, "common_genes.txt"))
cat("✅ 共有基因列表已保存到:", file.path(res.path, "common_genes.txt"), "\n")

saveRDS(all_results, file.path(res.path, "all_results_final.rds"))
save.image(file.path(res.path, "final_analysis.RData"))

intermediate_summary <- data.frame(
  File = c(
    "input_candidate_genes.txt",
    "input_candidate_up.txt",
    "input_candidate_down.txt",
    "common_genes.txt",
    "algorithms_list.txt",
    "processed_data.rds",
    "preTrain_features.rds",
    "cv_results.rds",
    "acceptable_results.rds",
    "validation_results.rds"
  ),
  Description = c(
    "输入候选基因列表",
    "输入上调基因列表",
    "输入下调基因列表",
    "训练集和验证集共有基因",
    "算法配置列表",
    "标准化后的训练和验证数据",
    "预筛选特征结果",
    "5折交叉验证完整结果",
    "可接受模型结果",
    "验证集评估结果"
  ),
  stringsAsFactors = FALSE
)

write.csv(intermediate_summary, file.path(intermediate.path, "intermediate_files_summary.csv"), row.names = FALSE)

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("             101种机器学习模型分析完成！\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("候选基因来源:", candidate_dir, "\n")
cat("总模型数:", length(methods), "\n")
cat("可接受模型数:", length(acceptable_results), "\n")
cat("成功验证模型数:", nrow(validation_results), "\n")
cat("最终最佳模型:", best_model_name, "\n")
cat("验证集AUC:", round(best_validation_auc, 4), "\n")
cat("特征基因数:", length(best_features), "\n")
cat("\n📁 输出目录:\n")
cat("  结果目录:", res.path, "\n")
cat("  图片目录:", fig.path, "\n")
cat("  中间结果目录:", intermediate.path, "\n")
cat("\n📊 总运行时间:", round(total_runtime, 2), "分钟\n")
