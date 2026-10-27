
cat("\n========================================\n")
cat("14a. 药物预测 - CellMiner\n")
cat("========================================\n\n")

rm(list = ls())

options(bitmapType = "cairo")
options(device = "pdf")

setwd("")

output_dir <- "./14_Drug_Prediction"
res.path <- file.path(output_dir, "Results")
fig.path <- file.path(output_dir, "Figures")

for (path in c(output_dir, res.path, fig.path)) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

cat("\n========== 加载 R 包 ==========\n")

packages <- c("ggplot2", "dplyr", "tidyr", "stringr", "igraph",
              "ggraph", "tidygraph", "ggrepel", "RColorBrewer", "pheatmap")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org", quiet = TRUE)
    library(pkg, character.only = TRUE)
  }
  cat("✅", pkg, "\n")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "http://cran.us.r-project.org", quiet = TRUE)
}

cellminer_available <- FALSE
tryCatch({
  if (!require("rcellminer", quietly = TRUE)) {
    cat("安装 rcellminer (Bioconductor)...\n")
    BiocManager::install("rcellminer", update = FALSE, ask = FALSE, quiet = TRUE)
    library(rcellminer)
  }
  if (!require("rcellminerData", quietly = TRUE)) {
    BiocManager::install("rcellminerData", update = FALSE, ask = FALSE, quiet = TRUE)
    library(rcellminerData)
  }

  if (!require("Biobase", quietly = TRUE)) {
    BiocManager::install("Biobase", update = FALSE, ask = FALSE, quiet = TRUE)
    library(Biobase)
  }
  cellminer_available <- TRUE
  cat("✅ rcellminer + rcellminerData + Biobase\n")
}, error = function(e) {
  cat("❌ rcellminer 加载失败:", e$message, "\n")
  cat("   请手动安装: BiocManager::install('rcellminer')\n")
})

if (!cellminer_available) {
  stop("❌ 无 rcellminer, 无法继续 CellMiner 分析")
}

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate

cat("\n========== 读取生物标志物 ==========\n")

biomarkers_file <- "./07_Expression_Validation/Results/selected_biomarkers.txt"
if (!file.exists(biomarkers_file)) stop("❌ selected_biomarkers.txt 不存在")

biomarker_genes <- readLines(biomarkers_file)
biomarker_genes <- biomarker_genes[nzchar(biomarker_genes) & !is.na(biomarker_genes)]
cat("生物标志物:", paste(biomarker_genes, collapse = ", "), "\n")

cat("\n========== Step 1: 加载 NCI-60 表达 + 药物数据 ==========\n")

safe_get_matrix <- function(obj) {
  if (is.matrix(obj)) return(obj)
  if (inherits(obj, "ExpressionSet")) return(Biobase::exprs(obj))
  if (inherits(obj, "eSetNCI60"))     return(Biobase::exprs(obj))
  if (is.data.frame(obj))             return(as.matrix(obj))

  tryCatch(as.matrix(obj), error = function(e) {
    stop("无法将对象转为 matrix: class=", paste(class(obj), collapse = ","))
  })
}

drug_act_raw <- getAct(rcellminerData::drugData)
nci60_drug_act <- safe_get_matrix(drug_act_raw)
cat("药物活性对象 class:", paste(class(drug_act_raw), collapse = ","), "\n")

mol_data <- getAllFeatureData(rcellminerData::molData)
cat("molData 包含的特征类型:", paste(names(mol_data), collapse = ", "), "\n")

exp_slot <- if ("exp" %in% names(mol_data)) {
  "exp"
} else if ("xai" %in% names(mol_data)) {
  "xai"
} else {
  names(mol_data)[1]
}
cat("使用基因表达 slot:", exp_slot, "\n")

gene_exp_raw <- mol_data[[exp_slot]]
nci60_gene_exp <- safe_get_matrix(gene_exp_raw)
cat("基因表达对象 class:", paste(class(gene_exp_raw), collapse = ","), "\n")

cat("NCI-60 基因表达矩阵: ", nrow(nci60_gene_exp), "基因 x ",
    ncol(nci60_gene_exp), "细胞系\n")
cat("NCI-60 药物活性矩阵: ", nrow(nci60_drug_act), "药物 x ",
    ncol(nci60_drug_act), "细胞系\n")

sample_rownames <- head(rownames(nci60_gene_exp), 5)
cat("基因表达矩阵行名示例:", paste(sample_rownames, collapse = ", "), "\n")
if (any(grepl("^exp_", sample_rownames))) {
  cat("⚠️ 行名带有 'exp_' 前缀, 自动去除\n")
  rownames(nci60_gene_exp) <- gsub("^exp_", "", rownames(nci60_gene_exp))
}

drug_annot_raw <- getFeatureAnnot(rcellminerData::drugData)[["drug"]]
cat("药物注释对象 class:", paste(class(drug_annot_raw), collapse = ","), "\n")
drug_annot <- if (is.data.frame(drug_annot_raw)) {
  drug_annot_raw
} else {
  as.data.frame(drug_annot_raw)
}
cat("药物注释字段:", paste(colnames(drug_annot), collapse = ", "), "\n")

cat("\n========== Step 2: 筛选 FDA / Clinical trial 药物 ==========\n")

fda_status_col <- intersect(c("FDA_STATUS", "FDA.Status", "FDA_Status"),
                            colnames(drug_annot))[1]
if (is.na(fda_status_col)) {
  cat("⚠️ 找不到 FDA_STATUS 字段, 使用全部药物\n")
  fda_drugs_idx <- seq_len(nrow(drug_annot))
} else {
  cat("使用字段:", fda_status_col, "\n")
  cat("FDA_STATUS 分布:\n")
  print(table(drug_annot[[fda_status_col]], useNA = "ifany"))

  keep <- grepl("FDA|Clinical", drug_annot[[fda_status_col]], ignore.case = TRUE)
  fda_drugs_idx <- which(keep)
}

cat(sprintf("筛选后药物数: %d / %d\n", length(fda_drugs_idx), nrow(drug_annot)))

fda_drug_ids <- rownames(drug_annot)[fda_drugs_idx]
common_drugs <- intersect(fda_drug_ids, rownames(nci60_drug_act))
nci60_drug_act_fda <- nci60_drug_act[common_drugs, , drop = FALSE]

cat("最终用于相关性分析的药物数:", nrow(nci60_drug_act_fda), "\n")

cat("\n========== Step 3: 计算 Spearman 相关性 ==========\n")

available_genes <- intersect(biomarker_genes, rownames(nci60_gene_exp))
missing_genes <- setdiff(biomarker_genes, rownames(nci60_gene_exp))
cat("可用生物标志物:", paste(available_genes, collapse = ", "), "\n")
if (length(missing_genes) > 0) {
  cat("⚠️ 缺失生物标志物:", paste(missing_genes, collapse = ", "), "\n")
}

if (length(available_genes) == 0) {
  stop("❌ 没有生物标志物在 NCI-60 表达矩阵中")
}

common_cells <- intersect(colnames(nci60_gene_exp), colnames(nci60_drug_act_fda))
cat("共同细胞系数:", length(common_cells), "\n")

gene_expr_sub <- nci60_gene_exp[available_genes, common_cells, drop = FALSE]
drug_act_sub  <- nci60_drug_act_fda[, common_cells, drop = FALSE]

cat("\n开始计算相关性 (基因 × 药物 =",
    length(available_genes), "×", nrow(drug_act_sub), ")...\n")

cor_results <- data.frame()
for (gene in available_genes) {
  expr_vec <- as.numeric(gene_expr_sub[gene, ])

  for (i in seq_len(nrow(drug_act_sub))) {
    drug_id <- rownames(drug_act_sub)[i]
    act_vec <- as.numeric(drug_act_sub[i, ])

    valid_idx <- which(!is.na(expr_vec) & !is.na(act_vec))
    if (length(valid_idx) < 10) next

    test <- tryCatch(
      cor.test(expr_vec[valid_idx], act_vec[valid_idx],
               method = "spearman", exact = FALSE),
      error = function(e) NULL
    )
    if (is.null(test)) next

    cor_results <- rbind(cor_results, data.frame(
      Gene = gene,
      Drug_NSC = drug_id,
      Cor = as.numeric(test$estimate),
      P_value = test$p.value,
      N = length(valid_idx),
      stringsAsFactors = FALSE
    ))
  }
  cat(sprintf("  %s: 计算完成 (累计 %d 条)\n", gene, nrow(cor_results)))
}

cat(sprintf("\n相关性计算完成: %d 个 gene-drug 对\n", nrow(cor_results)))

cat("\n========== Step 4: 筛选显著相关 (|r|>0.4 & p<0.05) ==========\n")

sig_results <- cor_results %>%
  filter(abs(Cor) > 0.4, P_value < 0.05) %>%
  arrange(Gene, desc(abs(Cor)))

cat(sprintf("显著药物-基因关系对: %d\n", nrow(sig_results)))

name_col <- intersect(c("NAME", "Name", "DRUG_NAME"), colnames(drug_annot))[1]
moa_col  <- intersect(c("MOA", "Mechanism", "MECHANISM"), colnames(drug_annot))[1]

sig_results$Drug_Name <- if (!is.na(name_col)) {
  drug_annot[sig_results$Drug_NSC, name_col]
} else {
  sig_results$Drug_NSC
}
sig_results$MOA <- if (!is.na(moa_col)) {
  drug_annot[sig_results$Drug_NSC, moa_col]
} else {
  NA
}

sig_results$Drug_Name[is.na(sig_results$Drug_Name) | sig_results$Drug_Name == ""] <-
  sig_results$Drug_NSC[is.na(sig_results$Drug_Name) | sig_results$Drug_Name == ""]

gene_drug_count <- sig_results %>%
  group_by(Gene) %>%
  summarise(N_drugs = n(),
            Mean_abs_cor = mean(abs(Cor)),
            .groups = "drop")
cat("\n每个生物标志物的显著药物数:\n")
print(gene_drug_count)

write.csv(cor_results,
          file.path(res.path, "all_gene_drug_correlations.csv"),
          row.names = FALSE)
write.csv(sig_results,
          file.path(res.path, "significant_drug_gene_pairs.csv"),
          row.names = FALSE)

cat("\n========== Step 5: 构建 mRNA-药物网络 ==========\n")

if (nrow(sig_results) == 0) {
  cat("⚠️ 没有显著结果, 跳过网络绘制\n")
} else {

  edges <- data.frame(
    from = sig_results$Gene,
    to   = sig_results$Drug_Name,
    cor  = sig_results$Cor,
    stringsAsFactors = FALSE
  )

  all_nodes <- unique(c(edges$from, edges$to))
  node_df <- data.frame(
    name = all_nodes,
    type = ifelse(all_nodes %in% biomarker_genes, "Gene", "Drug"),
    stringsAsFactors = FALSE
  )

  cat(sprintf("网络规模: %d 个节点 (%d 基因 + %d 药物), %d 边\n",
              nrow(node_df),
              sum(node_df$type == "Gene"),
              sum(node_df$type == "Drug"),
              nrow(edges)))

  write.csv(edges, file.path(res.path, "network_edges.csv"), row.names = FALSE)
  write.csv(node_df, file.path(res.path, "network_nodes.csv"), row.names = FALSE)

  g <- graph_from_data_frame(edges, vertices = node_df, directed = FALSE)
  tg <- as_tbl_graph(g)

  set.seed(123)

  p_net <- ggraph(tg, layout = "fr") +
    geom_edge_link(aes(edge_alpha = abs(cor),
                       edge_color = cor),
                   edge_width = 0.6) +
    scale_edge_color_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35",
                               midpoint = 0,
                               name = "Spearman r",
                               limits = c(-1, 1)) +
    scale_edge_alpha(range = c(0.3, 0.9), guide = "none") +
    geom_node_point(aes(color = type, size = type)) +
    scale_color_manual(values = c(Gene = "#E64B35", Drug = "#4DBBD5"),
                       name = "Node type") +
    scale_size_manual(values = c(Gene = 10, Drug = 4), guide = "none") +
    geom_node_text(aes(label = ifelse(type == "Gene", name, "")),
                   color = "black", size = 5, fontface = "bold",
                   nudge_y = 0.05, show.legend = FALSE) +
    geom_node_text(aes(label = ifelse(type == "Drug", name, "")),
                   color = "gray20", size = 2.6, fontface = "plain",
                   repel = TRUE,
                   max.overlaps = Inf,
                   box.padding = 0.3,
                   point.padding = 0.3,
                   segment.size = 0.15,
                   segment.color = "gray70",
                   show.legend = FALSE) +
    labs(title = "mRNA-Drug Regulatory Network (CellMiner)",
         subtitle = sprintf("Spearman correlation between %d biomarkers and FDA/Clinical trial drugs (|r|>0.4, p<0.05)",
                            length(available_genes))) +
    theme_void(base_size = 13) +
    theme(plot.background  = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          plot.title       = element_text(hjust = 0.5, size = 15, face = "bold"),
          plot.subtitle    = element_text(hjust = 0.5, size = 10, color = "gray40"),
          legend.position  = "right",
          legend.title     = element_text(size = 11, face = "bold"),
          legend.text      = element_text(size = 10),
          plot.margin      = margin(15, 25, 15, 25))

  ggsave(file.path(fig.path, "mRNA_Drug_Network.pdf"),
         p_net, width = 14, height = 10)
  ggsave(file.path(fig.path, "mRNA_Drug_Network.png"),
         p_net, width = 14, height = 10, dpi = 300)
  cat("✅ mRNA-药物网络图已保存\n")

  igraph::write_graph(g, file.path(res.path, "network.graphml"), format = "graphml")
  write.csv(edges, file.path(res.path, "network_for_cytoscape_edges.csv"), row.names = FALSE)
  write.csv(node_df, file.path(res.path, "network_for_cytoscape_nodes.csv"), row.names = FALSE)
  cat("✅ Cytoscape 文件已导出\n")
}

cat("\n========== Step 6: Top 药物相关性热图 ==========\n")

if (nrow(sig_results) > 0) {

  top_pairs <- sig_results %>%
    group_by(Gene) %>%
    slice_max(abs(Cor), n = 10) %>%
    ungroup()

  unique_drugs <- unique(top_pairs$Drug_Name)
  cor_mat <- matrix(NA,
                    nrow = length(available_genes),
                    ncol = length(unique_drugs))
  rownames(cor_mat) <- available_genes
  colnames(cor_mat) <- unique_drugs

  for (i in seq_len(nrow(top_pairs))) {
    cor_mat[top_pairs$Gene[i], top_pairs$Drug_Name[i]] <- top_pairs$Cor[i]
  }

  cor_mat[is.na(cor_mat)] <- 0

  draw_heatmap <- function() {
    pheatmap(cor_mat,
             color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
             breaks = seq(-1, 1, length.out = 101),
             main = "Top Gene-Drug Correlations (Spearman, CellMiner)",
             display_numbers = TRUE,
             number_format = "%.2f",
             fontsize_number = 7,
             fontsize_row = 12,
             fontsize_col = 8,
             cluster_rows = FALSE,
             cluster_cols = TRUE,
             border_color = "gray85",
             angle_col = 45,
             cellheight = 26)
  }

  pdf(file.path(fig.path, "Gene_Drug_Correlation_Heatmap.pdf"),
      width = max(10, ncol(cor_mat) * 0.25 + 4), height = 5)
  draw_heatmap()
  dev.off()

  png(file.path(fig.path, "Gene_Drug_Correlation_Heatmap.png"),
      width = max(10, ncol(cor_mat) * 0.25 + 4), height = 5,
      units = "in", res = 300)
  draw_heatmap()
  dev.off()
  cat("✅ 相关性热图已保存\n")
}

cat("\n========== Step 7: 候选化合物 Top 表 + SMILES (供 ADMETlab) ==========\n")

if (nrow(sig_results) > 0) {

  top10 <- sig_results %>%
    arrange(desc(abs(Cor))) %>%
    head(10) %>%
    select(Drug_Name, Drug_NSC, Gene, Cor, P_value, MOA, N)
  cat("\n按 |r| 排序的 Top 10 药物-基因对:\n")
  print(top10)
  write.csv(top10, file.path(res.path, "Top10_drug_gene_pairs.csv"),
            row.names = FALSE)

  cat("\n从 CellMiner 提取 Top 候选药物的 SMILES (供 ADMETlab 2.0)...\n")

  smiles_col <- intersect(c("SMILES", "smiles", "Smiles"), colnames(drug_annot))[1]
  if (is.na(smiles_col)) {
    cat("⚠️ drug_annot 没有 SMILES 字段, 跳过 SMILES 导出\n")
    cat("   你需要去 PubChem 手动查找候选药物的 SMILES\n")
  } else {

    unique_drugs <- sig_results %>%
      group_by(Drug_NSC, Drug_Name) %>%
      summarise(max_abs_cor = max(abs(Cor)), .groups = "drop") %>%
      arrange(desc(max_abs_cor))

    unique_drugs$SMILES <- drug_annot[unique_drugs$Drug_NSC, smiles_col]
    unique_drugs <- unique_drugs[!is.na(unique_drugs$SMILES) &
                                 nchar(unique_drugs$SMILES) > 0, ]

    cat(sprintf("提取到 %d 个独特候选药物的 SMILES\n", nrow(unique_drugs)))

    write.csv(unique_drugs,
              file.path(res.path, "candidate_drugs_with_SMILES.csv"),
              row.names = FALSE)

    writeLines(unique_drugs$SMILES,
               file.path(res.path, "SMILES_for_ADMETlab.txt"))

    write.table(
      data.frame(SMILES = unique_drugs$SMILES,
                 Name   = unique_drugs$Drug_Name),
      file.path(res.path, "SMILES_with_names_for_ADMETlab.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )

    cat("✅ SMILES 文件已导出:\n")
    cat("   - candidate_drugs_with_SMILES.csv  (完整表, R 可读)\n")
    cat("   - SMILES_for_ADMETlab.txt          (纯 SMILES, 复制粘贴到 ADMETlab)\n")
    cat("   - SMILES_with_names_for_ADMETlab.tsv  (带药物名, 用于追溯)\n")
  }
}

cat("\n========== 生成报告 ==========\n")

report_lines <- c(
  "# CellMiner 药物预测报告",
  paste("生成时间:", Sys.time()),
  "",
  "## 方法",
  "使用 CellMiner 数据库 (NCI-60 癌细胞系) 计算 5 个生物标志物表达水平与 FDA 批准 / 临床试验药物 IC50 之间的 Spearman 相关性, 筛选 |r| > 0.4 且 p < 0.05 的药物-基因关系对。",
  "",
  "## 数据来源",
  paste("- R 包: rcellminer, rcellminerData (Bioconductor)"),
  paste("- 数据库: NCI-60 (https://discover.nci.nih.gov/cellminer/)"),
  paste("- 细胞系数:", length(common_cells)),
  paste("- 候选药物 (FDA / Clinical trial):", nrow(nci60_drug_act_fda)),
  paste("- 生物标志物:", paste(biomarker_genes, collapse = ", ")),
  paste("- NCI-60 中可用的生物标志物:", paste(available_genes, collapse = ", ")),
  "",
  "## 筛选结果",
  paste("- 总相关性计算数:", nrow(cor_results)),
  paste("- 显著药物-基因对 (|r|>0.4 & p<0.05):", nrow(sig_results)),
  ""
)

if (nrow(sig_results) > 0) {
  report_lines <- c(report_lines, "## 各基因显著药物数")
  for (i in seq_len(nrow(gene_drug_count))) {
    report_lines <- c(report_lines,
      sprintf("- **%s**: %d 个药物 (平均 |r| = %.3f)",
              gene_drug_count$Gene[i],
              gene_drug_count$N_drugs[i],
              gene_drug_count$Mean_abs_cor[i]))
  }
  report_lines <- c(report_lines, "")

  report_lines <- c(report_lines, "## Top 10 药物-基因对")
  for (i in seq_len(min(10, nrow(sig_results)))) {
    report_lines <- c(report_lines,
      sprintf("- %s ↔ %s (r = %.3f, p = %.2e)",
              sig_results$Drug_Name[i],
              sig_results$Gene[i],
              sig_results$Cor[i],
              sig_results$P_value[i]))
  }
}

report_lines <- c(report_lines,
  "",
  "## 输出文件",
  "",
  "### Results/",
  "- all_gene_drug_correlations.csv  (全部相关性)",
  "- significant_drug_gene_pairs.csv  (显著筛选后)",
  "- Top10_drug_gene_pairs.csv        (Top 10 候选)",
  "- network_edges.csv / network_nodes.csv",
  "- network.graphml                  (Cytoscape)",
  "",
  "### Figures/",
  "- mRNA_Drug_Network.pdf/png        (mRNA-药物调控网络)",
  "- Gene_Drug_Correlation_Heatmap.pdf/png  (Top 相关性热图)",
  "",
  "## 下一步 (第十五部分)",
  "1. 从 PubChem 下载 Top 候选药物的 3D-SDF",
  "2. 使用 CB-Dock 与生物标志物蛋白结构进行分子对接",
  "3. 使用 PyMoL 可视化对接复合物"
)

writeLines(report_lines, file.path(res.path, "Drug_Prediction_Report.md"))

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("    第十四部分: CellMiner 药物预测 完成\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("方法: CellMiner NCI-60 细胞系\n")
cat("筛选: FDA/Clinical trial 药物 + Spearman |r|>0.4 & p<0.05\n\n")

cat(sprintf("生物标志物: %d 个 (%s)\n",
            length(biomarker_genes), paste(biomarker_genes, collapse = ", ")))
cat(sprintf("可用于 NCI-60 的生物标志物: %d 个\n", length(available_genes)))
cat(sprintf("候选药物数: %d (FDA / Clinical trial)\n", nrow(nci60_drug_act_fda)))
cat(sprintf("显著药物-基因对: %d\n\n", nrow(sig_results)))

if (nrow(sig_results) > 0) {
  cat("各基因显著药物数:\n")
  for (i in seq_len(nrow(gene_drug_count))) {
    cat(sprintf("  %-8s : %3d 个 (mean |r| = %.3f)\n",
                gene_drug_count$Gene[i],
                gene_drug_count$N_drugs[i],
                gene_drug_count$Mean_abs_cor[i]))
  }
}

cat("\n输出文件:\n")
cat("  Figures/mRNA_Drug_Network.pdf/png\n")
cat("  Figures/Gene_Drug_Correlation_Heatmap.pdf/png\n")
cat("  Results/significant_drug_gene_pairs.csv\n")
cat("  Results/Top10_drug_gene_pairs.csv\n")
cat("  Results/candidate_drugs_with_SMILES.csv\n")
cat("  Results/SMILES_for_ADMETlab.txt   ← 用这个去 ADMETlab 跑 ADMET\n")
cat("  Results/network.graphml  (Cytoscape)\n")

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("    下一步: ADMET 筛选 (手动 + 自动两步)\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("\n1. 打开浏览器: https://admetmesh.scbdd.com/service/screening/index\n")
cat("2. 把 Results/SMILES_for_ADMETlab.txt 的内容粘贴到输入框\n")
cat("3. 选 'All endpoints', 点 Submit, 等待计算\n")
cat("4. 计算完成后, 点 'Download', 下载 CSV 结果\n")
cat("5. 把下载的 CSV 上传到服务器:\n")
cat("   scp admet_result.csv km112:14_Drug_Prediction/Results/\n")
cat("6. 跑 r.14_Drug_Prediction_ADMET.r 完成最终筛选\n")

cat("\n✅ 第十四部分 (CellMiner) 完成\n")

sink(file.path(res.path, "session_info.txt"))
sessionInfo()
sink()
