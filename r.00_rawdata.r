
rm(list = ls())

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("                    初始化\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

setwd("")

if (!dir.exists("./00_rawdata")) {
  dir.create("./00_rawdata")
  cat("创建00_rawdata文件夹\n")
}

setwd("./00_rawdata")

cat("\n加载R包...\n")
packages <- c("GEOquery", "limma", "dplyr", "tidyr", "BiocManager")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    if (pkg == "BiocManager") {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    } else {
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
    library(pkg, character.only = TRUE)
  } else {
    cat("✅ 已加载:", pkg, "\n")
  }
}

if (!require("homologene", quietly = TRUE)) {
  install.packages("homologene", repos = "https://cloud.r-project.org")
  library(homologene)
} else {
  cat("✅ 已加载: homologene\n")
}

log_file <- "processing_log.txt"
cat("数据处理日志\n", file = log_file)
cat(paste("开始时间:", Sys.time(), "\n"), file = log_file, append = TRUE)

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("              1. 处理 GSE16561 (训练集)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

if (file.exists("GSE16561_ready.rds")) {
  cat("✅ GSE16561_ready.rds 已存在，跳过处理\n")
  gse16561 <- readRDS("GSE16561_ready.rds")
} else {

  cat("使用GEOquery下载GSE16561...\n")
  gset16561 <- getGEO("GSE16561", destdir = ".", GSEMatrix = TRUE)
  gset16561 <- gset16561[[1]]

  expr16561 <- exprs(gset16561)
  cat("表达矩阵维度:", dim(expr16561), "\n")

  qx <- as.numeric(quantile(expr16561, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
  LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (LogC) {
    expr16561[which(expr16561 <= 0)] <- 0
    expr16561 <- log2(expr16561 + 1)
    cat("  log2转换完成\n")
  } else {
    cat("  无需log2转换\n")
  }

  expr16561 <- as.data.frame(expr16561)

  cat("获取GPL6883平台注释...\n")
  gpl <- getGEO(gset16561@annotation, destdir = ".")
  gpl <- Table(gpl)

  probe2symbol <- dplyr::select(gpl, 'ID', 'Symbol')
  colnames(probe2symbol) <- c('ID', 'symbol')
  probe2symbol <- filter(probe2symbol, symbol != '' & !is.na(symbol))
  probe2symbol <- probe2symbol[!grepl('^---$', probe2symbol$symbol), ]
  probe2symbol$symbol <- gsub(' ', '', probe2symbol$symbol)

  expr16561$ID <- rownames(expr16561)
  dat16561 <- merge(expr16561, probe2symbol, by = 'ID') %>%
    dplyr::select(-ID) %>%
    dplyr::select(symbol, everything())

  sample_cols <- grep('GSM', colnames(dat16561), value = TRUE)

  dat16561 <- dat16561 %>%
    group_by(symbol) %>%
    summarise(across(all_of(sample_cols), ~ mean(.x, na.rm = TRUE))) %>%
    as.data.frame()

  rownames(dat16561) <- dat16561$symbol
  dat16561 <- dat16561[, -1, drop = FALSE]

  pd16561 <- pData(gset16561)
  group16561 <- data.frame(
    sample = pd16561$geo_accession,
    title = pd16561$title,
    stringsAsFactors = FALSE
  )

  group16561$group <- ifelse(
    grepl("Control", group16561$title, ignore.case = TRUE),
    "control",
    ifelse(grepl("Stroke", group16561$title, ignore.case = TRUE),
           "stroke", NA)
  )

  cat("\nGSE16561分组统计:\n")
  print(table(group16561$group))

  group16561 <- group16561[order(group16561$group), ]
  dat16561 <- dat16561[, group16561$sample]

  gse16561 <- list(expression = dat16561, group = group16561)
  saveRDS(gse16561, "GSE16561_ready.rds")
  write.csv(dat16561, "GSE16561_expression.csv")
  write.csv(group16561, "GSE16561_group.csv", row.names = FALSE)

  cat("✅ GSE16561处理完成\n")
  cat("  样本数:", ncol(dat16561), "\n")
  cat("  基因数:", nrow(dat16561), "\n")
}

cat(paste(Sys.time(), "- GSE16561完成\n"), file = log_file, append = TRUE)

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("              2. 处理 GSE58294 (验证集)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

if (file.exists("GSE58294_ready.rds")) {
  cat("✅ GSE58294_ready.rds 已存在，跳过处理\n")
  gse58294 <- readRDS("GSE58294_ready.rds")
} else {

  cat("使用GEOquery下载GSE58294...\n")
  gset58294 <- getGEO("GSE58294", destdir = ".", GSEMatrix = TRUE)
  gset58294 <- gset58294[[1]]

  expr58294 <- exprs(gset58294)
  cat("表达矩阵维度:", dim(expr58294), "\n")

  qx <- as.numeric(quantile(expr58294, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
  LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (LogC) {
    expr58294[which(expr58294 <= 0)] <- 0
    expr58294 <- log2(expr58294 + 1)
    cat("  log2转换完成\n")
  } else {
    cat("  无需log2转换\n")
  }

  expr58294 <- as.data.frame(expr58294)

  cat("获取GPL570平台注释...\n")
  gpl570 <- getGEO("GPL570", destdir = ".")
  gpl570_table <- Table(gpl570)

  probe2symbol_df <- dplyr::select(gpl570_table, ID, `Gene Symbol`)
  colnames(probe2symbol_df) <- c("ID", "symbol")

  probe2symbol_df <- probe2symbol_df %>%
    filter(symbol != "" & !is.na(symbol)) %>%
    filter(!grepl("^---$", symbol))

  probe2symbol_df$symbol <- gsub(" ///.*", "", probe2symbol_df$symbol)
  probe2symbol_df$symbol <- gsub(" ", "", probe2symbol_df$symbol)
  probe2symbol_df$ID     <- as.character(probe2symbol_df$ID)

  expr58294$ID <- rownames(expr58294)
  dat58294 <- merge(expr58294, probe2symbol_df, by = 'ID') %>%
    dplyr::select(-ID) %>%
    dplyr::select(symbol, everything())

  sample_cols <- grep('GSM', colnames(dat58294), value = TRUE)

  dat58294 <- dat58294 %>%
    group_by(symbol) %>%
    summarise(across(all_of(sample_cols), ~ mean(.x, na.rm = TRUE))) %>%
    as.data.frame()

  rownames(dat58294) <- dat58294$symbol
  dat58294 <- dat58294[, -1, drop = FALSE]

  pd58294 <- pData(gset58294)
  group58294 <- data.frame(
    sample = pd58294$geo_accession,
    title = pd58294$title,
    stringsAsFactors = FALSE
  )

  group58294$group <- ifelse(
    grepl("Control", group58294$title, ignore.case = TRUE),
    "control",
    ifelse(grepl("Stroke", group58294$title, ignore.case = TRUE),
           "stroke", NA)
  )

  cat("\nGSE58294分组统计:\n")
  print(table(group58294$group))

  group58294 <- group58294[order(group58294$group), ]
  dat58294 <- dat58294[, group58294$sample]

  gse58294 <- list(expression = dat58294, group = group58294)
  saveRDS(gse58294, "GSE58294_ready.rds")
  write.csv(dat58294, "GSE58294_expression.csv")
  write.csv(group58294, "GSE58294_group.csv", row.names = FALSE)

  cat("✅ GSE58294处理完成\n")
  cat("  样本数:", ncol(dat58294), "\n")
  cat("  基因数:", nrow(dat58294), "\n")
}

cat(paste(Sys.time(), "- GSE58294完成\n"), file = log_file, append = TRUE)

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("              3. 处理小鼠自测数据\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat(paste(Sys.time(), "- 开始处理小鼠自测数据\n"), file = log_file, append = TRUE)

if (!file.exists("../mRNA_counts_all.txt")) {
  cat("⚠️ mRNA_counts_all.txt 不存在，跳过小鼠数据处理\n")
  has_mouse_data <- FALSE
  brain_mouse <- NULL
  spleen_mouse <- NULL
} else {
  has_mouse_data <- TRUE

  if (file.exists("brain_ready_mouse.rds") && file.exists("spleen_ready_mouse.rds")) {
    cat("✅ 小鼠数据已存在，跳过处理\n")
    brain_mouse <- readRDS("brain_ready_mouse.rds")
    spleen_mouse <- readRDS("spleen_ready_mouse.rds")
  } else {
    cat("使用命名规则进行分组...\n")

    counts_all <- read.delim("../mRNA_counts_all.txt", row.names = 1, check.names = FALSE)
    cat("表达矩阵维度:", dim(counts_all), "\n")
    cat("样本名:\n")
    print(colnames(counts_all))

    sample_info <- data.frame(
      sample_id = colnames(counts_all),
      tissue = NA,
      group = NA,
      strain = NA,
      stringsAsFactors = FALSE
    )

    cat("\n根据命名规则进行分组...\n")
    for (i in 1:nrow(sample_info)) {
      sample_id <- sample_info$sample_id[i]

      if (grepl("101$|102$|103$", sample_id)) {
        sample_info$tissue[i] <- "brain"
      } else if (grepl("201$|202$|203$", sample_id)) {
        sample_info$tissue[i] <- "spleen"
      }

      if (grepl("CM", sample_id)) {
        sample_info$strain[i] <- "C57"
        sample_info$group[i] <- "MACO"
      } else if (grepl("CN", sample_id)) {
        sample_info$strain[i] <- "C57"
        sample_info$group[i] <- "control"
      } else if (grepl("IN", sample_id)) {
        sample_info$strain[i] <- "IL-21"
        sample_info$group[i] <- "control"
      } else if (grepl("IM", sample_id)) {
        sample_info$strain[i] <- "IL-21"
        sample_info$group[i] <- "MACO"
      }
    }

    cat("\n分组统计:\n")
    print(table(sample_info$group, useNA = "ifany"))
    cat("\n组织类型统计:\n")
    print(table(sample_info$tissue, useNA = "ifany"))

    sample_info <- sample_info[!is.na(sample_info$group) & !is.na(sample_info$tissue), ]
    cat("\n过滤后有效样本数:", nrow(sample_info), "\n")

    if (nrow(sample_info) == 0) {
      cat("❌ 没有有效样本，请检查命名规则\n")
      has_mouse_data <- FALSE
      brain_mouse <- NULL
      spleen_mouse <- NULL
    } else {
      counts_filtered <- counts_all[, sample_info$sample_id]

      write.csv(sample_info, "mouse_sample_info.csv", row.names = FALSE)

      brain_samples <- sample_info$sample_id[sample_info$tissue == "brain"]
      if (length(brain_samples) > 0) {
        brain_counts <- counts_filtered[, brain_samples, drop = FALSE]
        brain_info <- sample_info[sample_info$tissue == "brain", ]
        colnames(brain_info)[1] <- "sample"

        brain_mouse <- list(expression = brain_counts, group = brain_info)
        saveRDS(brain_mouse, "brain_ready_mouse.rds")
        cat("✅ 脑组织小鼠数据已保存\n")
        cat("  脑组织样本数:", ncol(brain_counts), "\n")
      } else {
        brain_mouse <- NULL
      }

      spleen_samples <- sample_info$sample_id[sample_info$tissue == "spleen"]
      if (length(spleen_samples) > 0) {
        spleen_counts <- counts_filtered[, spleen_samples, drop = FALSE]
        spleen_info <- sample_info[sample_info$tissue == "spleen", ]
        colnames(spleen_info)[1] <- "sample"

        spleen_mouse <- list(expression = spleen_counts, group = spleen_info)
        saveRDS(spleen_mouse, "spleen_ready_mouse.rds")
        cat("✅ 脾脏组织小鼠数据已保存\n")
        cat("  脾脏组织样本数:", ncol(spleen_counts), "\n")
      } else {
        spleen_mouse <- NULL
      }
    }
  }
}

cat(paste(Sys.time(), "- 小鼠数据处理完成\n"), file = log_file, append = TRUE)

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("              4. 小鼠基因转人基因\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat(paste(Sys.time(), "- 开始基因转换\n"), file = log_file, append = TRUE)

if (!has_mouse_data || (!exists("brain_mouse") && !exists("spleen_mouse"))) {
  cat("\n⚠️ 没有小鼠数据，跳过基因转换\n")
  brain_human <- NULL
  spleen_human <- NULL
} else {

  if (!exists("brain_mouse") && file.exists("brain_ready_mouse.rds")) {
    brain_mouse <- readRDS("brain_ready_mouse.rds")
  }
  if (!exists("spleen_mouse") && file.exists("spleen_ready_mouse.rds")) {
    spleen_mouse <- readRDS("spleen_ready_mouse.rds")
  }

  if (exists("brain_mouse") && !is.null(brain_mouse)) {
    mouse_brain_genes <- rownames(brain_mouse$expression)
    cat("脑组织小鼠基因数:", length(mouse_brain_genes), "\n")

    cat("\n使用homologene进行基因转换...\n")
    brain_homolog <- homologene(mouse_brain_genes, inTax = 10090, outTax = 9606)
    cat("  homologene找到", nrow(brain_homolog), "个脑组织映射关系\n")

    brain_all_map <- data.frame(
      mouseGene = brain_homolog[,1],
      humanGene = brain_homolog[,2],
      method = "homologene",
      stringsAsFactors = FALSE
    )

    brain_all_map <- distinct(brain_all_map)
    brain_human_genes <- unique(brain_all_map$humanGene)
    cat("脑组织最终人类基因数:", length(brain_human_genes), "\n")

    write.csv(brain_all_map, "brain_mouse_human_map.csv", row.names = FALSE)

    cat("\n构建脑组织人类基因表达矩阵...\n")
    brain_expr_human <- matrix(0, nrow = length(brain_human_genes), ncol = ncol(brain_mouse$expression))
    rownames(brain_expr_human) <- brain_human_genes
    colnames(brain_expr_human) <- colnames(brain_mouse$expression)

    for (human_gene in brain_human_genes) {
      mouse_genes <- brain_all_map$mouseGene[brain_all_map$humanGene == human_gene]
      mouse_genes <- mouse_genes[mouse_genes %in% rownames(brain_mouse$expression)]
      if (length(mouse_genes) >= 1) {
        brain_expr_human[human_gene, ] <- colMeans(brain_mouse$expression[mouse_genes, , drop = FALSE])
      }
    }

    brain_human <- list(expression = brain_expr_human, group = brain_mouse$group)
    saveRDS(brain_human, "brain_ready_human.rds")
    cat("✅ 脑组织基因转换完成，人类基因数:", nrow(brain_expr_human), "\n")
  } else {
    brain_human <- NULL
  }

  if (exists("spleen_mouse") && !is.null(spleen_mouse)) {
    mouse_spleen_genes <- rownames(spleen_mouse$expression)
    cat("脾脏组织小鼠基因数:", length(mouse_spleen_genes), "\n")

    spleen_homolog <- homologene(mouse_spleen_genes, inTax = 10090, outTax = 9606)
    cat("  homologene找到", nrow(spleen_homolog), "个脾脏组织映射关系\n")

    spleen_all_map <- data.frame(
      mouseGene = spleen_homolog[,1],
      humanGene = spleen_homolog[,2],
      method = "homologene",
      stringsAsFactors = FALSE
    )

    spleen_all_map <- distinct(spleen_all_map)
    spleen_human_genes <- unique(spleen_all_map$humanGene)
    cat("脾脏组织最终人类基因数:", length(spleen_human_genes), "\n")

    write.csv(spleen_all_map, "spleen_mouse_human_map.csv", row.names = FALSE)

    cat("\n构建脾脏组织人类基因表达矩阵...\n")
    spleen_expr_human <- matrix(0, nrow = length(spleen_human_genes), ncol = ncol(spleen_mouse$expression))
    rownames(spleen_expr_human) <- spleen_human_genes
    colnames(spleen_expr_human) <- colnames(spleen_mouse$expression)

    for (human_gene in spleen_human_genes) {
      mouse_genes <- spleen_all_map$mouseGene[spleen_all_map$humanGene == human_gene]
      mouse_genes <- mouse_genes[mouse_genes %in% rownames(spleen_mouse$expression)]
      if (length(mouse_genes) >= 1) {
        spleen_expr_human[human_gene, ] <- colMeans(spleen_mouse$expression[mouse_genes, , drop = FALSE])
      }
    }

    spleen_human <- list(expression = spleen_expr_human, group = spleen_mouse$group)
    saveRDS(spleen_human, "spleen_ready_human.rds")
    cat("✅ 脾脏组织基因转换完成，人类基因数:", nrow(spleen_expr_human), "\n")
  } else {
    spleen_human <- NULL
  }
}

cat(paste(Sys.time(), "- 基因转换完成\n"), file = log_file, append = TRUE)

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("              5. 基因名标准化\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat(paste(Sys.time(), "- 开始基因名标准化\n"), file = log_file, append = TRUE)

standardize_gene_names <- function(expr_matrix, dataset_name) {
  cat("\n处理", dataset_name, "...\n")
  original_genes <- rownames(expr_matrix)
  cat("  原始基因数:", length(original_genes), "\n")

  upper_genes <- toupper(original_genes)

  if (any(duplicated(upper_genes))) {
    dup_count <- sum(duplicated(upper_genes))
    cat("  发现", dup_count, "个重复基因名\n")

    unique_genes <- unique(upper_genes)
    new_expr <- matrix(0, nrow = length(unique_genes), ncol = ncol(expr_matrix))
    rownames(new_expr) <- unique_genes
    colnames(new_expr) <- colnames(expr_matrix)

    for (gene in unique_genes) {
      idx <- which(upper_genes == gene)
      if (length(idx) == 1) {
        new_expr[gene, ] <- as.numeric(expr_matrix[idx, ])
      } else {
        new_expr[gene, ] <- colMeans(expr_matrix[idx, , drop = FALSE])
      }
    }
    cat("  合并后基因数:", nrow(new_expr), "\n")
    return(new_expr)
  } else {
    rownames(expr_matrix) <- upper_genes
    return(expr_matrix)
  }
}

if (exists("gse16561")) {
  gse16561$expression <- standardize_gene_names(gse16561$expression, "GSE16561")
  saveRDS(gse16561, "GSE16561_std.rds")
  cat("✅ GSE16561标准化完成\n")
}

if (exists("gse58294")) {
  gse58294$expression <- standardize_gene_names(gse58294$expression, "GSE58294")
  saveRDS(gse58294, "GSE58294_std.rds")
  cat("✅ GSE58294标准化完成\n")
}

if (exists("brain_human") && !is.null(brain_human)) {
  brain_human$expression <- standardize_gene_names(brain_human$expression, "Brain Human")
  saveRDS(brain_human, "brain_ready_human_std.rds")
  cat("✅ 脑组织人类基因标准化完成\n")
}

if (exists("spleen_human") && !is.null(spleen_human)) {
  spleen_human$expression <- standardize_gene_names(spleen_human$expression, "Spleen Human")
  saveRDS(spleen_human, "spleen_ready_human_std.rds")
  cat("✅ 脾脏组织人类基因标准化完成\n")
}

cat("\n✅ 基因名标准化完成\n")
cat(paste(Sys.time(), "- 基因名标准化完成\n"), file = log_file, append = TRUE)

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("              6. 准备DEG分析数据\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat(paste(Sys.time(), "- 开始准备\n"), file = log_file, append = TRUE)

if (file.exists("GSE16561_std.rds")) {
  file.copy("GSE16561_std.rds", "GSE16561_for_DEG.rds", overwrite = TRUE)
  cat("✅ GSE16561_for_DEG.rds 已准备\n")
}

if (file.exists("GSE58294_std.rds")) {
  file.copy("GSE58294_std.rds", "GSE58294_for_DEG.rds", overwrite = TRUE)
  cat("✅ GSE58294_for_DEG.rds 已准备\n")
}

if (file.exists("brain_ready_human_std.rds")) {
  file.copy("brain_ready_human_std.rds", "brain_for_DEG.rds", overwrite = TRUE)
  cat("✅ brain_for_DEG.rds 已准备\n")
}

if (file.exists("spleen_ready_human_std.rds")) {
  file.copy("spleen_ready_human_std.rds", "spleen_for_DEG.rds", overwrite = TRUE)
  cat("✅ spleen_for_DEG.rds 已准备\n")
}

cat("\n✅ DEG分析数据准备完成\n")

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("                    处理完成报告\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("📊 数据集统计:\n")
if (file.exists("GSE16561_for_DEG.rds")) {
  gse16561_final <- readRDS("GSE16561_for_DEG.rds")
  cat("  GSE16561 (训练集):", nrow(gse16561_final$expression), "基因,",
      ncol(gse16561_final$expression), "样本\n")
}
if (file.exists("GSE58294_for_DEG.rds")) {
  gse58294_final <- readRDS("GSE58294_for_DEG.rds")
  cat("  GSE58294 (验证集):", nrow(gse58294_final$expression), "基因,",
      ncol(gse58294_final$expression), "样本\n")
}
if (file.exists("brain_for_DEG.rds")) {
  brain_final <- readRDS("brain_for_DEG.rds")
  cat("  脑组织(人):", nrow(brain_final$expression), "基因,",
      ncol(brain_final$expression), "样本\n")
}
if (file.exists("spleen_for_DEG.rds")) {
  spleen_final <- readRDS("spleen_for_DEG.rds")
  cat("  脾脏组织(人):", nrow(spleen_final$expression), "基因,",
      ncol(spleen_final$expression), "样本\n")
}

cat("\n📁 输出文件:\n")
cat("  - GSE16561_std.rds: 训练集标准化数据\n")
cat("  - GSE58294_std.rds: 验证集标准化数据\n")
cat("  - brain_ready_human_std.rds: 脑组织人类基因数据\n")
cat("  - spleen_ready_human_std.rds: 脾脏组织人类基因数据\n")
cat("  - GSE16561_for_DEG.rds: 训练集DEG分析数据\n")
cat("  - GSE58294_for_DEG.rds: 验证集DEG分析数据\n")
cat("  - brain_for_DEG.rds: 脑组织DEG分析数据\n")
cat("  - spleen_for_DEG.rds: 脾脏组织DEG分析数据\n")
cat("  - mouse_sample_info.csv: 小鼠样本信息\n")
cat("  - brain_mouse_human_map.csv: 脑组织基因映射关系\n")
cat("  - spleen_mouse_human_map.csv: 脾脏组织基因映射关系\n")

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("                    全部处理完成！\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\n📝 重要说明:\n")
cat("  - 训练集和验证集未进行批次校正，保持独立\n")
cat("  - 分组信息根据title列正确构建\n")
cat("  - GSE58294使用GPL570平台文件进行探针注释\n")
cat("  - 小鼠基因使用homologene转为人基因\n")
cat("  - 去重方式统一为取平均值\n")
cat("\n✅ 所有处理完成！\n")

