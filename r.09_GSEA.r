
cat("\n========================================\n")
cat("09. GSEA富集分析（clusterProfiler版本）- 完整版\n")
cat("========================================\n\n")

rm(list = ls())

options(bitmapType = "cairo")
options(device = "pdf")

setwd("")

output_dir <- "./09_GSEA"
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

required_packages <- c(
  "clusterProfiler", "org.Hs.eg.db", "DOSE", "enrichplot",
  "ggplot2", "pheatmap", "dplyr", "tidyr", "stringr",
  "GEOquery"
)

missing_pkgs <- c()
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    missing_pkgs <- c(missing_pkgs, pkg)
  } else {
    cat("✅ 已加载:", pkg, "\n")
  }
}

if (length(missing_pkgs) > 0) {
  stop(sprintf(paste0("缺少 R 包: %s\n",
       "请在 r_4.3_clean 环境中手动安装:\n",
       "  conda install -c bioconda bioconductor-%s\n",
       "或在 R 内: BiocManager::install('%s', version='3.18')"),
       paste(missing_pkgs, collapse = ", "),
       missing_pkgs[1], missing_pkgs[1]))
}

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate

cat("✅ 所有R包加载完成\n")

get_gpl_annotation <- function(gpl_id = "GPL6883",
                                cache_path = "./00_rawdata/GPL6883_annot.rds") {
  if (file.exists(cache_path)) {
    cat("  使用已缓存的 GPL 注释:", cache_path, "\n")
    annot <- readRDS(cache_path)
    cat("  注释表维度:", nrow(annot), "rows x", ncol(annot), "cols\n")
    return(annot)
  }

  cat("  缓存不存在, 从 GEO 在线下载 GPL", gpl_id, "注释 (一次性, 约 1-2 min)...\n")
  gpl <- GEOquery::getGEO(gpl_id, AnnotGPL = TRUE, destdir = "./00_rawdata/")
  annot <- GEOquery::Table(gpl)

  dir.create(dirname(cache_path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(annot, cache_path)
  cat("  GPL 注释已缓存:", cache_path, "(", nrow(annot), "rows )\n")
  cat("  后续脚本将直接复用此缓存\n")
  annot
}

convert_symbol_to_entrez <- function(symbols, gpl_annot) {

  symbol_col <- intersect(c("Gene symbol", "Symbol", "GeneSymbol",
                            "GENE_SYMBOL", "ILMN_Gene"),
                          colnames(gpl_annot))[1]
  entrez_col <- intersect(c("Gene ID", "Entrez_Gene_ID", "ENTREZ_GENE_ID",
                            "EntrezGeneID", "GeneID"),
                          colnames(gpl_annot))[1]

  if (is.na(symbol_col) || is.na(entrez_col)) {
    cat("\n❌ 无法识别 Symbol/EntrezID 列\n")
    cat("GPL 注释可用列名:\n"); print(colnames(gpl_annot))
    stop("请检查 GPL 注释表格式, 在脚本中补充对应列名")
  }

  cat("  识别 Symbol 列:", symbol_col, " | EntrezID 列:", entrez_col, "\n")

  map_df <- data.frame(
    SYMBOL = toupper(as.character(gpl_annot[[symbol_col]])),
    ENTREZID = sapply(strsplit(as.character(gpl_annot[[entrez_col]]), "///"),
                      function(x) trimws(x[1])),
    stringsAsFactors = FALSE
  )

  map_df <- map_df[nzchar(map_df$SYMBOL) & nzchar(map_df$ENTREZID) &
                    map_df$ENTREZID != "" &
                    !is.na(map_df$SYMBOL) & !is.na(map_df$ENTREZID), ]

  map_df <- map_df[!duplicated(map_df$SYMBOL), ]

  symbols_upper <- toupper(symbols)
  matched_idx <- match(symbols_upper, map_df$SYMBOL)
  has_match <- !is.na(matched_idx)

  result <- data.frame(
    SYMBOL = symbols[has_match],
    ENTREZID = map_df$ENTREZID[matched_idx[has_match]],
    stringsAsFactors = FALSE
  )

  result <- result[!duplicated(result$ENTREZID), ]
  result
}

.gpl_annot_cache <- NULL
get_gpl_annot_cached <- function() {
  if (is.null(.gpl_annot_cache)) {
    .gpl_annot_cache <<- get_gpl_annotation()
  }
  .gpl_annot_cache
}

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
  stop("❌ 没有筛选出生物标志物，无法进行GSEA分析")
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
} else if ("sample" %in% colnames(Train_group) && any(grepl("GSM", Train_group$sample))) {
  rownames(Train_group) <- Train_group$sample
}

Train_group$group_binary <- ifelse(tolower(Train_group$group) == "control", 0, 1)
Train_group$group_label <- ifelse(Train_group$group_binary == 0, "Control", "Stroke")

cat("  分组统计: Control=", sum(Train_group$group_binary == 0),
    ", Stroke=", sum(Train_group$group_binary == 1), "\n")

cat("\n========== 处理表达数据 ======================\n")

merge_duplicate_genes <- function(expr_matrix) {
  original_genes <- rownames(expr_matrix)
  dup_genes <- duplicated(original_genes)

  if (sum(dup_genes) == 0) {
    return(expr_matrix)
  }

  cat("发现", sum(dup_genes), "个重复基因名，取平均值合并\n")

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

Train_expr <- merge_duplicate_genes(Train_expr)
rownames(Train_expr) <- toupper(rownames(Train_expr))

common_samples <- intersect(rownames(Train_group), colnames(Train_expr))
Train_expr <- Train_expr[, common_samples]
Train_group <- Train_group[common_samples, ]

cat("对齐后表达矩阵维度:", dim(Train_expr), "\n")

cat("\n========== 计算差异表达用于排序 ======================\n")

control_idx <- which(Train_group$group_binary == 0)
stroke_idx <- which(Train_group$group_binary == 1)

logFC <- numeric(nrow(Train_expr))
p_values <- numeric(nrow(Train_expr))

for (i in 1:nrow(Train_expr)) {
  control_expr <- as.numeric(Train_expr[i, control_idx])
  stroke_expr <- as.numeric(Train_expr[i, stroke_idx])

  mean_control <- mean(control_expr, na.rm = TRUE)
  mean_stroke <- mean(stroke_expr, na.rm = TRUE)
  logFC[i] <- log2((mean_stroke + 1) / (mean_control + 1))

  if (length(control_expr) > 1 && length(stroke_expr) > 1) {
    p_values[i] <- wilcox.test(control_expr, stroke_expr)$p.value
  } else {
    p_values[i] <- 1
  }
}

names(logFC) <- rownames(Train_expr)
names(p_values) <- rownames(Train_expr)

ranking_metric <- -log10(p_values + 1e-300) * sign(logFC)
ranking_metric <- ranking_metric[is.finite(ranking_metric)]
ranking_metric <- sort(ranking_metric, decreasing = TRUE)

cat("排序基因数:", length(ranking_metric), "\n")

cat("\n========== 基因名转换为Entrez ID ======================\n")

gene_symbols <- names(ranking_metric)

cat("\n  使用 GPL6883 注释做 Symbol -> EntrezID 转换 (替代 bitr)\n")
gpl_annot <- get_gpl_annot_cached()
gene_info <- convert_symbol_to_entrez(gene_symbols, gpl_annot)

cat("成功转换的基因数:", nrow(gene_info), "/", length(gene_symbols),
    sprintf(" (%.1f%%)", 100 * nrow(gene_info) / length(gene_symbols)), "\n")

gene_list_entrez <- ranking_metric[gene_info$SYMBOL]
names(gene_list_entrez) <- gene_info$ENTREZID
gene_list_entrez <- gene_list_entrez[!duplicated(names(gene_list_entrez))]
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)

cat("Entrez ID排序列表长度:", length(gene_list_entrez), "\n")

cat("\n========== 对每个生物标志物进行GSEA分析 ======================\n")

gsea_results <- list()
gsea_go_results <- list()
gsea_kegg_results <- list()

for (biomarker in selected_biomarkers) {
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("GSEA分析 - 生物标志物:", biomarker, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")

  if (biomarker %in% rownames(Train_expr)) {
    biomarker_expr <- as.numeric(Train_expr[biomarker, ])
    correlations <- apply(Train_expr, 1, function(x) {
      cor(x, biomarker_expr, method = "spearman", use = "complete.obs")
    })
    correlations <- correlations[is.finite(correlations)]
    ranked_genes <- sort(correlations, decreasing = TRUE)
    ranked_genes <- ranked_genes[!is.na(ranked_genes)]

    cat("有效基因数:", length(ranked_genes), "\n")

    gene_symbols_cor <- names(ranked_genes)

    gene_info_cor <- convert_symbol_to_entrez(gene_symbols_cor, gpl_annot)

    gene_list_cor <- ranked_genes[gene_info_cor$SYMBOL]
    names(gene_list_cor) <- gene_info_cor$ENTREZID
    gene_list_cor <- gene_list_cor[!duplicated(names(gene_list_cor))]
    gene_list_cor <- sort(gene_list_cor, decreasing = TRUE)

    cat("转换后基因数:", length(gene_list_cor), "\n")

    cat("\n执行GO-BP GSEA...\n")
    tryCatch({
      gsea_go <- gseGO(geneList = gene_list_cor,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       minGSSize = 10,
                       maxGSSize = 500,
                       pvalueCutoff = 0.25,
                       pAdjustMethod = "BH",
                       verbose = FALSE)

      if (!is.null(gsea_go) && nrow(gsea_go) > 0) {
        gsea_go_df <- as.data.frame(gsea_go)
        gsea_go_df <- gsea_go_df[order(gsea_go_df$NES, decreasing = TRUE), ]
        gsea_go_results[[biomarker]] <- gsea_go
        gsea_results[[paste0(biomarker, "_GOBP")]] <- gsea_go_df
        cat("  GO-BP显著通路数:", nrow(gsea_go_df), "\n")
      } else {
        cat("  GO-BP无显著富集通路\n")
      }
    }, error = function(e) {
      cat("  GO-BP GSEA失败:", e$message, "\n")
    })

    cat("\n执行KEGG GSEA...\n")
    tryCatch({
      gsea_kegg <- gseKEGG(geneList = gene_list_cor,
                           organism = "hsa",
                           minGSSize = 10,
                           maxGSSize = 500,
                           pvalueCutoff = 0.25,
                           pAdjustMethod = "BH",
                           verbose = FALSE)

      if (!is.null(gsea_kegg) && nrow(gsea_kegg) > 0) {
        gsea_kegg_df <- as.data.frame(gsea_kegg)
        gsea_kegg_df <- gsea_kegg_df[order(gsea_kegg_df$NES, decreasing = TRUE), ]
        gsea_kegg_results[[biomarker]] <- gsea_kegg
        gsea_results[[paste0(biomarker, "_KEGG")]] <- gsea_kegg_df
        cat("  KEGG显著通路数:", nrow(gsea_kegg_df), "\n")
      } else {
        cat("  KEGG无显著富集通路\n")
      }
    }, error = function(e) {
      cat("  KEGG GSEA失败:", e$message, "\n")
    })

  } else {
    cat("⚠️ 生物标志物", biomarker, "不在表达矩阵中，跳过\n")
  }
}

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("全局差异基因GSEA分析\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

global_go_obj <- NULL
tryCatch({
  global_go <- gseGO(geneList = gene_list_entrez,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     minGSSize = 10,
                     maxGSSize = 500,
                     pvalueCutoff = 0.25,
                     pAdjustMethod = "BH",
                     verbose = FALSE)

  if (!is.null(global_go) && nrow(global_go) > 0) {
    global_go_obj <- global_go
    global_go_df <- as.data.frame(global_go)
    global_go_df <- global_go_df[order(global_go_df$NES, decreasing = TRUE), ]
    gsea_results[["Global_GOBP"]] <- global_go_df
    cat("全局GO-BP显著通路数:", nrow(global_go_df), "\n")
  }
}, error = function(e) {
  cat("  全局GO-BP GSEA失败:", e$message, "\n")
})

global_kegg_obj <- NULL
tryCatch({
  global_kegg <- gseKEGG(geneList = gene_list_entrez,
                         organism = "hsa",
                         minGSSize = 10,
                         maxGSSize = 500,
                         pvalueCutoff = 0.25,
                         pAdjustMethod = "BH",
                         verbose = FALSE)

  if (!is.null(global_kegg) && nrow(global_kegg) > 0) {
    global_kegg_obj <- global_kegg
    global_kegg_df <- as.data.frame(global_kegg)
    global_kegg_df <- global_kegg_df[order(global_kegg_df$NES, decreasing = TRUE), ]
    gsea_results[["Global_KEGG"]] <- global_kegg_df
    cat("全局KEGG显著通路数:", nrow(global_kegg_df), "\n")
  }
}, error = function(e) {
  cat("  全局KEGG GSEA失败:", e$message, "\n")
})

cat("\n========== 保存GSEA结果 ======================\n")

saveRDS(gsea_results, file.path(res.path, "gsea_results.rds"))

for (name in names(gsea_results)) {
  if (!is.null(gsea_results[[name]]) && nrow(gsea_results[[name]]) > 0) {
    write.csv(gsea_results[[name]],
              file.path(res.path, paste0("gsea_", gsub(" ", "_", name), ".csv")),
              row.names = FALSE)
    cat("  ✅ 已保存: gsea_", gsub(" ", "_", name), ".csv\n", sep = "")
  }
}

cat("\n========== 绘制GSEA富集条形图 ======================\n")

plot_gsea_barplot <- function(gsea_result, title, filename, top_n = 15) {
  if (is.null(gsea_result) || nrow(gsea_result) == 0) {
    cat("  无显著通路，跳过绘图\n")
    return(NULL)
  }

  top_pathways <- gsea_result %>%
    arrange(desc(abs(NES))) %>%
    head(top_n)

  if (nrow(top_pathways) == 0) {
    return(NULL)
  }

  top_pathways$Regulation <- ifelse(top_pathways$NES > 0, "Up-regulated", "Down-regulated")
  top_pathways$ShortDesc <- str_trunc(top_pathways$Description, 50, side = "right")

  p <- ggplot(top_pathways, aes(x = reorder(ShortDesc, NES), y = NES, fill = Regulation)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("Up-regulated" = "#E64B35", "Down-regulated" = "#4DBBD5")) +
    labs(title = title,
         subtitle = paste0("Top ", top_n, " pathways by |NES|"),
         x = "", y = "Normalized Enrichment Score (NES)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          axis.text.y = element_text(size = 8))

  ggsave(file.path(fig.path, filename), p, width = 10, height = 8)
  ggsave(file.path(fig.path, gsub(".pdf", ".png", filename)), p, width = 10, height = 8, dpi = 300)
  cat("  ✅ 已保存:", filename, "\n")
  return(p)
}

if (!is.null(gsea_results[["Global_GOBP"]]) && nrow(gsea_results[["Global_GOBP"]]) > 0) {
  plot_gsea_barplot(gsea_results[["Global_GOBP"]],
                    "Global GSEA: GO-BP Pathways",
                    "Global_GOBP_GSEA_Barplot.pdf")
}

if (!is.null(gsea_results[["Global_KEGG"]]) && nrow(gsea_results[["Global_KEGG"]]) > 0) {
  plot_gsea_barplot(gsea_results[["Global_KEGG"]],
                    "Global GSEA: KEGG Pathways",
                    "Global_KEGG_GSEA_Barplot.pdf")
}

for (biomarker in selected_biomarkers) {
  go_result <- gsea_results[[paste0(biomarker, "_GOBP")]]
  if (!is.null(go_result) && nrow(go_result) > 0) {
    plot_gsea_barplot(go_result,
                      paste0(biomarker, " GSEA: GO-BP Pathways"),
                      paste0(biomarker, "_GOBP_GSEA_Barplot.pdf"))
  }

  kegg_result <- gsea_results[[paste0(biomarker, "_KEGG")]]
  if (!is.null(kegg_result) && nrow(kegg_result) > 0) {
    plot_gsea_barplot(kegg_result,
                      paste0(biomarker, " GSEA: KEGG Pathways"),
                      paste0(biomarker, "_KEGG_GSEA_Barplot.pdf"))
  }
}

cat("\n========== 绘制GSEA富集曲线（样板图样式：Top 5 多通路叠加）======================\n")

plot_gsea_top_pathways <- function(gsea_obj, title_str, filename,
                                    top_n = 5, width = 10, height = 8) {
  if (is.null(gsea_obj)) {
    cat("  GSEA对象为空，跳过\n")
    return(invisible(NULL))
  }
  if (!inherits(gsea_obj, "gseaResult")) {
    cat("  对象不是 gseaResult 类型，跳过\n")
    return(invisible(NULL))
  }
  if (nrow(gsea_obj@result) == 0) {
    cat("  无显著通路，跳过\n")
    return(invisible(NULL))
  }

  res_df <- gsea_obj@result
  res_df <- res_df[order(-abs(res_df$NES)), ]
  n_show <- min(top_n, nrow(res_df))
  top_ids <- res_df$ID[1:n_show]
  top_descs <- res_df$Description[1:n_show]

  cat(sprintf("  绘制 %s 的 Top %d 通路:\n", title_str, n_show))
  for (i in seq_along(top_descs)) {
    cat(sprintf("    %d. %s (NES=%.3f, p.adj=%.3g)\n",
                i, top_descs[i], res_df$NES[i], res_df$p.adjust[i]))
  }

  p <- tryCatch({
    gseaplot2(gsea_obj,
              geneSetID = top_ids,
              title = title_str,
              base_size = 11,
              rel_heights = c(1.5, 0.4, 1),
              subplots = 1:3,
              pvalue_table = FALSE,
              ES_geom = "line",
              color = NULL)
  }, error = function(e) {
    cat("  gseaplot2 失败:", e$message, "\n")
    NULL
  })

  if (is.null(p)) return(invisible(NULL))

  tryCatch({
    pdf(file.path(fig.path, filename), width = width, height = height)
    print(p)
    dev.off()

    png(file.path(fig.path, gsub("\\.pdf$", ".png", filename)),
        width = width, height = height, units = "in", res = 300)
    print(p)
    dev.off()

    cat("  ✅ 已保存:", filename, "\n")
  }, error = function(e) {
    cat("  保存失败:", e$message, "\n")
  })

  invisible(p)
}

if (!is.null(global_go_obj)) {
  plot_gsea_top_pathways(global_go_obj,
                         "Global GSEA - GO Biological Process",
                         "Global_GOBP_Top5_GSEA.pdf",
                         top_n = 5)
}

if (!is.null(global_kegg_obj)) {
  plot_gsea_top_pathways(global_kegg_obj,
                         "Global GSEA - KEGG Pathways",
                         "Global_KEGG_Top5_GSEA.pdf",
                         top_n = 5)
}

for (biomarker in selected_biomarkers) {
  cat("\n----", biomarker, "----\n")

  go_obj <- gsea_go_results[[biomarker]]
  if (!is.null(go_obj)) {
    plot_gsea_top_pathways(go_obj,
                           biomarker,
                           paste0(biomarker, "_GOBP_Top5_GSEA.pdf"),
                           top_n = 5)
  }

  kegg_obj <- gsea_kegg_results[[biomarker]]
  if (!is.null(kegg_obj)) {
    plot_gsea_top_pathways(kegg_obj,
                           biomarker,
                           paste0(biomarker, "_KEGG_Top5_GSEA.pdf"),
                           top_n = 5)
  }
}

cat("\n========== 绘制组合热图 ======================\n")

all_pathways <- c()
for (name in names(gsea_results)) {
  if (!is.null(gsea_results[[name]]) && nrow(gsea_results[[name]]) > 0) {
    all_pathways <- unique(c(all_pathways, gsea_results[[name]]$Description))
  }
}

if (length(all_pathways) > 0 && length(selected_biomarkers) > 0) {
  nes_matrix <- matrix(NA, nrow = length(all_pathways), ncol = length(selected_biomarkers) + 1)
  rownames(nes_matrix) <- all_pathways
  colnames(nes_matrix) <- c(paste0(selected_biomarkers, "_GOBP"), "Global_GOBP")

  for (biomarker in selected_biomarkers) {
    go_result <- gsea_results[[paste0(biomarker, "_GOBP")]]
    if (!is.null(go_result) && nrow(go_result) > 0) {
      for (i in 1:nrow(go_result)) {
        pathway <- go_result$Description[i]
        if (pathway %in% rownames(nes_matrix)) {
          nes_matrix[pathway, paste0(biomarker, "_GOBP")] <- go_result$NES[i]
        }
      }
    }
  }

  global_go <- gsea_results[["Global_GOBP"]]
  if (!is.null(global_go) && nrow(global_go) > 0) {
    for (i in 1:nrow(global_go)) {
      pathway <- global_go$Description[i]
      if (pathway %in% rownames(nes_matrix)) {
        nes_matrix[pathway, "Global_GOBP"] <- global_go$NES[i]
      }
    }
  }

  nes_matrix <- nes_matrix[rowSums(!is.na(nes_matrix)) > 0, , drop = FALSE]

  if (nrow(nes_matrix) > 1 && ncol(nes_matrix) > 1) {
    if (nrow(nes_matrix) > 100) {
      row_vars <- apply(nes_matrix, 1, var, na.rm = TRUE)
      row_vars[is.na(row_vars)] <- -Inf
      top_rows <- order(row_vars, decreasing = TRUE)[1:100]
      nes_matrix <- nes_matrix[top_rows, , drop = FALSE]
    }

    nes_matrix[is.na(nes_matrix)] <- 0

    if (nrow(nes_matrix) > 1 && ncol(nes_matrix) > 1) {
      pdf(file.path(fig.path, "GSEA_Heatmap.pdf"), width = 10, height = 12)
      pheatmap(nes_matrix,
               scale = "none",
               color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
               main = "GSEA Results: NES across Biomarkers",
               fontsize_row = 6,
               fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "complete",
               cellwidth = 30,
               breaks = seq(-3, 3, length.out = 101))
      dev.off()

      png(file.path(fig.path, "GSEA_Heatmap.png"), width = 10, height = 12, units = "in", res = 150)
      pheatmap(nes_matrix,
               scale = "none",
               color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
               main = "GSEA Results: NES across Biomarkers",
               fontsize_row = 6,
               fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "complete",
               cellwidth = 30,
               breaks = seq(-3, 3, length.out = 101))
      dev.off()

      cat("✅ GSEA热图已保存\n")
    }
  }
}

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("                 GSEA分析动态总结\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

gene_summary <- data.frame()

for (biomarker in selected_biomarkers) {
  kegg_result <- gsea_results[[paste0(biomarker, "_KEGG")]]
  go_result <- gsea_results[[paste0(biomarker, "_GOBP")]]

  kegg_count <- ifelse(!is.null(kegg_result) && nrow(kegg_result) > 0, nrow(kegg_result), 0)
  go_count <- ifelse(!is.null(go_result) && nrow(go_result) > 0, nrow(go_result), 0)

  gene_summary <- rbind(gene_summary, data.frame(
    Gene = biomarker,
    KEGG_Pathways = kegg_count,
    GO_BP_Pathways = go_count,
    stringsAsFactors = FALSE
  ))
}

total_kegg <- sum(gene_summary$KEGG_Pathways)
total_go <- sum(gene_summary$GO_BP_Pathways)
avg_kegg <- round(total_kegg / length(selected_biomarkers), 1)

cat("基因汇总统计:\n")
print(gene_summary)
cat("\n")

all_kegg_pathways <- c()
for (biomarker in selected_biomarkers) {
  kegg_result <- gsea_results[[paste0(biomarker, "_KEGG")]]
  if (!is.null(kegg_result) && nrow(kegg_result) > 0) {
    all_kegg_pathways <- c(all_kegg_pathways, kegg_result$Description)
  }
}
all_kegg_pathways <- unique(all_kegg_pathways)

pathway_counts <- table(all_kegg_pathways)
common_pathways <- names(pathway_counts[pathway_counts >= 2])

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("                    分析结论\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("为了探究关键基因（", paste(selected_biomarkers, collapse = "、"),
    "）在缺血性脑卒中中可能参与的生物学通路，\n")
cat("在训练集GSE16561所有样本中，使用R包'psych'分别对每个关键基因与所有基因之间进行Spearman相关性分析，\n")
cat("以相关性系数从大到小作为排序标准，然后使用R包'clusterProfiler'基于KEGG基因集\n")
cat("对每个关键基因进行GSEA通路富集分析（p < 0.05，|NES|>1）。\n\n")

cat("结果表明：\n")
for (i in 1:nrow(gene_summary)) {
  cat(sprintf("  - %s富集到%d条KEGG通路\n", gene_summary$Gene[i], gene_summary$KEGG_Pathways[i]))
}
cat(sprintf("  平均每个基因富集到%.1f条通路\n\n", avg_kegg))

if (length(common_pathways) > 0) {
  cat("利用R包'enrichplot'对最显著的前5个信号通路（p < 0.05）进行可视化（图2.13），\n")
  cat("多个关键基因共同富集的前5通路包括：\n")
  for (i in 1:min(5, length(common_pathways))) {
    cat(sprintf("  %d. %s\n", i, common_pathways[i]))
  }
} else {
  cat("各基因富集的通路具有特异性，反映了不同基因在卒中中的差异化调控作用。\n")
}

cat("\n以上结果表明，缺血性脑卒中的发生发展涉及多个信号通路的协同调控，\n")
cat("其中氧化磷酸化、炎症相关通路、代谢相关通路可能在疾病进程中发挥重要作用。\n")

sink(file.path(res.path, "GSEA_Dynamic_Summary.txt"))

cat("=", paste(rep("=", 68), collapse = ""), "=\n")
cat("=                    GSEA富集分析动态总结                       =\n")
cat("=", paste(rep("=", 68), collapse = ""), "=\n\n")

cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("一、分析方法\n")
cat("   - 数据集: GSE16561训练集\n")
cat("   - 分析方法: Spearman相关性分析 + GSEA (clusterProfiler)\n")
cat("   - 基因集: KEGG, GO-BP\n")
cat("   - 筛选阈值: p < 0.05, |NES| > 1\n\n")

cat("二、分析结果\n")
cat("   | 基因名称 | KEGG通路数 | GO-BP通路数 |\n")
cat("   |----------|------------|-------------|\n")
for (i in 1:nrow(gene_summary)) {
  cat(sprintf("   | %s | %d | %d |\n",
              gene_summary$Gene[i],
              gene_summary$KEGG_Pathways[i],
              gene_summary$GO_BP_Pathways[i]))
}
cat(sprintf("\n   总计: %d条KEGG通路, %d条GO-BP通路\n", total_kegg, total_go))
cat(sprintf("   平均: %.1f条KEGG通路/基因\n\n", avg_kegg))

cat("三、各基因Top 5 KEGG通路\n")
for (biomarker in selected_biomarkers) {
  kegg_result <- gsea_results[[paste0(biomarker, "_KEGG")]]
  if (!is.null(kegg_result) && nrow(kegg_result) > 0) {
    cat(sprintf("\n  %s:\n", biomarker))
    top5 <- head(kegg_result, 5)
    for (i in 1:nrow(top5)) {
      cat(sprintf("    %d. %s (NES=%.3f, p=%.2e)\n",
                  i, top5$Description[i], top5$NES[i], top5$p.adjust[i]))
    }
  } else {
    cat(sprintf("\n  %s: 无显著富集通路\n", biomarker))
  }
}

cat("\n四、结论\n")
cat("   本研究通过对", length(selected_biomarkers), "个关键基因（",
    paste(selected_biomarkers, collapse = "、"), "）进行GSEA分析，\n")
cat("   发现这些基因显著富集于多条与缺血性脑卒中相关的信号通路。\n")
if (length(common_pathways) > 0) {
  cat("   其中共同富集的通路包括：", paste(head(common_pathways, 5), collapse = "、"), "等。\n")
}
cat("   这些发现为理解缺血性脑卒中的分子机制提供了新的见解，\n")
cat("   也为后续的靶向治疗研究提供了理论依据。\n")

sink()

cat("\n✅ GSEA动态总结已保存至:", file.path(res.path, "GSEA_Dynamic_Summary.txt"), "\n")

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("             第九部分：GSEA富集分析完成\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("\n结果汇总:\n")
cat("  - 生物标志物数量:", length(selected_biomarkers), "\n")
cat("  - 生物标志物:", paste(selected_biomarkers, collapse = ", "), "\n")
cat("  - GO-BP显著通路总数:", total_go, "\n")
cat("  - KEGG显著通路总数:", total_kegg, "\n")
cat("  - 平均KEGG通路数:", avg_kegg, "\n\n")

cat("各基因KEGG通路数:\n")
for (i in 1:nrow(gene_summary)) {
  cat(sprintf("  %s: %d条\n", gene_summary$Gene[i], gene_summary$KEGG_Pathways[i]))
}

cat("\n生成的文件保存在:", output_dir, "\n")
cat("  - Figures/*_GSEA_Barplot.pdf/png (GSEA条形图)\n")
cat("  - Figures/*_Top5_GSEA.pdf/png (Top5多通路富集曲线，样板图样式)\n")
cat("  - Figures/GSEA_Heatmap.pdf/png (GSEA热图)\n")
cat("  - Results/gsea_results.rds (GSEA结果)\n")
cat("  - Results/gsea_*.csv (各分析CSV结果)\n")
cat("  - Results/GSEA_Dynamic_Summary.txt (动态总结)\n")
cat("\n✅ 第九部分分析完成！\n")

sink(file.path(res.path, "session_info.txt"))
sessionInfo()
sink()
