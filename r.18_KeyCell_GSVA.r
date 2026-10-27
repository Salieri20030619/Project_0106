
cat("\n========================================\n")
cat("18. 筛选关键细胞 + GSVA 通路富集\n")
cat("========================================\n\n")

rm(list = ls())
options(bitmapType = "cairo")
options(device = "pdf")

setwd("")

output_dir <- "./18_KeyCell_GSVA"
res.path   <- file.path(output_dir, "Results")
fig.path   <- file.path(output_dir, "Figures")
for (p in c(output_dir, res.path, fig.path)) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE)
}

LOCAL_KEGG_GMT <- "./msigdb_cache/kegg_mouse.gmt"

cat("========== 加载 R 包 ==========\n")

required_packages <- c(
  "Seurat", "dplyr", "tidyr", "ggplot2", "pheatmap",
  "RColorBrewer", "patchwork", "ggrepel", "stringr",
  "GSVA", "GSEABase", "limma", "msigdbr"
)

missing_pkgs <- c()
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    missing_pkgs <- c(missing_pkgs, pkg)
  } else {
    cat("OK", pkg, "\n")
  }
}

if (length(missing_pkgs) > 0) {
  cat("\n缺少 R 包:", paste(missing_pkgs, collapse = ", "), "\n")
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(missing_pkgs, update = FALSE, ask = FALSE)
  for (pkg in missing_pkgs) {
    library(pkg, character.only = TRUE)
  }
}

cat(sprintf("msigdbr 版本: %s\n", packageVersion("msigdbr")))
cat(sprintf("GSVA 版本: %s\n", packageVersion("GSVA")))

cat("\n========== Step 1: 加载数据 ==========\n")

seurat_file <- "./17_annotation/seurat_obj_annotated_clean.rds"
if (!file.exists(seurat_file)) stop("seurat_obj_annotated_clean.rds 不存在, 请先运行 17 部分注释")

obj <- readRDS(seurat_file)
cat("细胞数:", ncol(obj), "\n")
cat("基因数:", nrow(obj), "\n")
cat("细胞类型:\n")
print(table(obj$cell_type))

if (!"condition" %in% colnames(obj@meta.data)) {
  cat("\ncondition 字段缺失, 从 orig.ident 推断...\n")
  gsm_to_group <- function(gsm) {
    num <- as.numeric(gsub("GSM", "", gsm))
    if (is.na(num)) return(NA)
    if (num <= 7060820) return("Sham") else return("Stroke")
  }
  obj$condition <- sapply(obj$orig.ident, gsm_to_group)
  cat("  OK condition 已补充\n")
}

cat("\n实验分组:\n")
print(table(obj$condition))

if ("Age" %in% colnames(obj@meta.data)) {
  cat("\nAge 分布:\n"); print(table(obj$Age))
} else {
  stop("数据缺失 Age 字段, 无法做 Young 子集分析")
}

biomarker_file <- "./07_Expression_Validation/Results/selected_biomarkers.txt"
if (file.exists(biomarker_file)) {
  biomarkers_human <- readLines(biomarker_file)
  biomarkers_human <- biomarkers_human[nchar(biomarkers_human) > 0]
  cat("\nOK 从第七部分读取生物标志物:", paste(biomarkers_human, collapse = ", "), "\n")
} else {
  biomarkers_human <- c("C5AR1", "ABCA1", "MFGE8")
}

to_mouse_symbol <- function(human_sym) {
  paste0(substr(human_sym, 1, 1),
         tolower(substr(human_sym, 2, nchar(human_sym))))
}
biomarkers_mouse <- to_mouse_symbol(biomarkers_human)
cat("人 -> 小鼠基因名:\n")
for (i in seq_along(biomarkers_human)) {
  cat(sprintf("  %s -> %s\n", biomarkers_human[i], biomarkers_mouse[i]))
}

available_genes <- rownames(obj)
biomarkers_found <- biomarkers_mouse[biomarkers_mouse %in% available_genes]
biomarkers_missing <- setdiff(biomarkers_mouse, available_genes)

cat(sprintf("\n在数据中找到 %d / %d 个生物标志物: %s\n",
            length(biomarkers_found), length(biomarkers_mouse),
            paste(biomarkers_found, collapse = ", ")))
if (length(biomarkers_missing) > 0) {
  cat("未找到:", paste(biomarkers_missing, collapse = ", "), "\n")
  for (m in biomarkers_missing) {
    found <- available_genes[toupper(available_genes) == toupper(m)]
    if (length(found) > 0) {
      cat(sprintf("  %s -> %s (大小写兜底)\n", m, found[1]))
      biomarkers_found <- c(biomarkers_found, found[1])
    }
  }
}

if (length(biomarkers_found) == 0) stop("没有生物标志物在数据中找到")

cat(sprintf("\n最终使用 %d 个生物标志物: %s\n",
            length(biomarkers_found), paste(biomarkers_found, collapse = ", ")))

score_cell_types <- function(seu, biomarkers, tag = "Full") {
  cat(sprintf("\n----- 打分 (%s 数据) -----\n", tag))
  cat(sprintf("细胞数: %d\n", ncol(seu)))

  Idents(seu) <- "cell_type"
  cell_types <- sort(unique(as.character(seu$cell_type)))

  avg_expr <- AverageExpression(seu, features = biomarkers,
                                 group.by = "cell_type",
                                 assays = "RNA", slot = "data")$RNA
  avg_expr_df <- as.data.frame(avg_expr)

  expr_mat <- GetAssayData(seu, assay = "RNA", layer = "data")
  pct_expr_df <- matrix(0, nrow = length(biomarkers),
                         ncol = length(cell_types),
                         dimnames = list(biomarkers, cell_types))
  for (g in biomarkers) {
    for (ct in cell_types) {
      cells <- WhichCells(seu, idents = ct)
      if (length(cells) == 0) next
      expr_vals <- expr_mat[g, cells]
      pct_expr_df[g, ct] <- mean(expr_vals > 0) * 100
    }
  }

  log2fc_df <- matrix(NA, nrow = length(biomarkers),
                       ncol = length(cell_types),
                       dimnames = list(biomarkers, cell_types))
  pvalue_df <- matrix(NA, nrow = length(biomarkers),
                       ncol = length(cell_types),
                       dimnames = list(biomarkers, cell_types))

  for (ct in cell_types) {
    ct_cells <- WhichCells(seu, idents = ct)
    if (length(ct_cells) < 20) next
    ct_obj <- subset(seu, cells = ct_cells)
    ct_sham   <- WhichCells(ct_obj, expression = condition == "Sham")
    ct_stroke <- WhichCells(ct_obj, expression = condition == "Stroke")
    if (length(ct_sham) < 5 || length(ct_stroke) < 5) next
    for (g in biomarkers) {
      sham_expr   <- expr_mat[g, ct_sham]
      stroke_expr <- expr_mat[g, ct_stroke]
      log2fc_df[g, ct] <- log2(mean(expm1(stroke_expr)) + 1) -
                           log2(mean(expm1(sham_expr))   + 1)
      pvalue_df[g, ct] <- tryCatch(
        wilcox.test(sham_expr, stroke_expr)$p.value,
        error = function(e) NA
      )
    }
  }

  common_ct <- intersect(colnames(avg_expr_df),
                          intersect(colnames(pct_expr_df), colnames(log2fc_df)))
  avg_expr_df <- avg_expr_df[, common_ct, drop = FALSE]
  pct_expr_df <- pct_expr_df[, common_ct, drop = FALSE]
  log2fc_df   <- log2fc_df[,   common_ct, drop = FALSE]

  score_mat <- matrix(0, nrow = length(biomarkers), ncol = length(common_ct),
                       dimnames = list(biomarkers, common_ct))
  for (g in biomarkers) {
    score_mat[g, ] <- as.numeric(avg_expr_df[g, ]) *
                      (as.numeric(pct_expr_df[g, ]) / 100) *
                      abs(as.numeric(log2fc_df[g, ]))
  }
  score_mat[is.na(score_mat)] <- 0

  total_score <- colSums(score_mat, na.rm = TRUE)
  ranked <- sort(total_score, decreasing = TRUE)

  rank_df <- data.frame(
    Cell_Type = names(ranked),
    Total_Score = round(as.numeric(ranked), 4),
    Rank = seq_along(ranked)
  )
  for (g in biomarkers) {
    rank_df[[sprintf("%s_score", g)]] <- round(as.numeric(score_mat[g, names(ranked)]), 4)
  }

  cat("Ranking:\n"); print(rank_df)

  list(rank_df = rank_df, score_mat = score_mat,
       avg_expr = avg_expr_df, pct_expr = pct_expr_df,
       log2fc = log2fc_df, pvalue = pvalue_df)
}

cat("\n========== Step 2: 全数据 表达分析 + 打分 ==========\n")

full_scoring <- score_cell_types(obj, biomarkers_found, tag = "Full")

write.csv(full_scoring$avg_expr,
          file.path(res.path, "biomarker_avg_expression_by_celltype.csv"))
write.csv(full_scoring$pct_expr,
          file.path(res.path, "biomarker_pct_expressing_by_celltype.csv"))
write.csv(full_scoring$log2fc,
          file.path(res.path, "biomarker_log2FC_stroke_vs_sham.csv"))
write.csv(full_scoring$pvalue,
          file.path(res.path, "biomarker_wilcox_pvalue.csv"))
write.csv(full_scoring$rank_df,
          file.path(res.path, "keycell_ranking.csv"), row.names = FALSE)

key_full <- full_scoring$rank_df$Cell_Type[1]
cat(sprintf("\n>>> 全数据 Top 1: %s\n", key_full))

cat("\n========== Step 3: Young 子集 表达分析 + 打分 ==========\n")

obj_young <- subset(obj, Age == "Young")
cat(sprintf("Young 子集: %d 细胞 (Sham %d, Stroke %d)\n",
            ncol(obj_young),
            sum(obj_young$condition == "Sham"),
            sum(obj_young$condition == "Stroke")))

young_scoring <- score_cell_types(obj_young, biomarkers_found, tag = "Young")

write.csv(young_scoring$rank_df,
          file.path(res.path, "keycell_ranking_YOUNG.csv"), row.names = FALSE)

key_young <- young_scoring$rank_df$Cell_Type[1]
cat(sprintf("\n>>> Young 子集 Top 1: %s\n", key_young))

cat("\n========== Step 4: 合并全数据和 Young 的 Top 1 ==========\n")

key_cells <- unique(c(key_full, key_young))

cat("关键细胞:\n")
cat(sprintf("  · 全数据 Top 1: %s\n", key_full))
cat(sprintf("  · Young 子集 Top 1: %s\n", key_young))
if (length(key_cells) == 1) {
  cat(sprintf("  · 两者一致, 最终 1 个关键细胞: %s\n", key_cells))
} else {
  cat(sprintf("  · 合并后最终 %d 个关键细胞: %s\n",
              length(key_cells), paste(key_cells, collapse = ", ")))
}

writeLines(key_cells, file.path(res.path, "key_cells.txt"))

compare_df <- merge(
  full_scoring$rank_df[, c("Cell_Type", "Total_Score", "Rank")],
  young_scoring$rank_df[, c("Cell_Type", "Total_Score", "Rank")],
  by = "Cell_Type", suffixes = c("_Full", "_Young"), all = TRUE
)
compare_df <- compare_df[order(compare_df$Rank_Full), ]
write.csv(compare_df,
          file.path(res.path, "comparison_Full_vs_Young.csv"), row.names = FALSE)
cat("\n全数据 vs Young 对比表:\n"); print(compare_df)

cat("\n========== Step 4b: 新旧注释体系关键细胞对比 (主管要求) ==========\n")

find_old_ranking <- function() {
  candidates <- c(
    Sys.glob("./18_KeyCell_GSVA_BACKUP_*/Results/keycell_ranking.csv"),
    Sys.glob("./18_KeyCell_GSVA_BACKUP*/Results/keycell_ranking.csv"),
    Sys.glob("./18_BACKUP*/Results/keycell_ranking.csv")
  )
  candidates_y <- c(
    Sys.glob("./18_KeyCell_GSVA_BACKUP_*/Results/keycell_ranking_YOUNG.csv"),
    Sys.glob("./18_KeyCell_GSVA_BACKUP*/Results/keycell_ranking_YOUNG.csv"),
    Sys.glob("./18_BACKUP*/Results/keycell_ranking_YOUNG.csv")
  )
  if (length(candidates) > 0 && length(candidates_y) > 0) {
    return(list(full = candidates[1], young = candidates_y[1]))
  }
  NULL
}

build_compare_table <- function() {
  old_paths <- find_old_ranking()
  if (is.null(old_paths)) {
    cat("  未找到原 18 部分备份 (./18_KeyCell_GSVA_BACKUP_*),\n")
    cat("  跳过旧 vs 新对比表生成. 如需对比, 请确保原版备份在该位置.\n")
    return(NULL)
  }

  cat(sprintf("  从备份读取原版 ranking:\n    %s\n    %s\n",
              old_paths$full, old_paths$young))

  v17_full  <- read.csv(old_paths$full,  stringsAsFactors = FALSE)
  v17_young <- read.csv(old_paths$young, stringsAsFactors = FALSE)

  cat(sprintf("\n  旧体系 ranking: %d 个细胞类型\n",   nrow(v17_full)))
  cat(sprintf("  新体系 ranking: %d 个细胞类型\n",  nrow(full_scoring$rank_df)))

  v17_key_full  <- v17_full$Cell_Type[1]
  v17_key_young <- v17_young$Cell_Type[1]
  v17_keys <- unique(c(v17_key_full, v17_key_young))

  v172_keys <- key_cells

  summary_table <- data.frame(
    项目 = c(
      "细胞类型总数",
      "全数据 Top 1 关键细胞",
      "全数据 Top 1 Score",
      "Young 子集 Top 1 关键细胞",
      "Young 子集 Top 1 Score",
      "最终关键细胞 (Full ∪ Young)",
      "最终关键细胞数目"
    ),
    `旧体系 (原 11 类)` = c(
      as.character(nrow(v17_full)),
      v17_key_full,
      sprintf("%.4f", v17_full$Total_Score[1]),
      v17_key_young,
      sprintf("%.4f", v17_young$Total_Score[1]),
      paste(v17_keys, collapse = ", "),
      as.character(length(v17_keys))
    ),
    `新体系` = c(
      as.character(nrow(full_scoring$rank_df)),
      key_full,
      sprintf("%.4f", full_scoring$rank_df$Total_Score[1]),
      key_young,
      sprintf("%.4f", young_scoring$rank_df$Total_Score[1]),
      paste(v172_keys, collapse = ", "),
      as.character(length(v172_keys))
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  cat("\n----- 主管要求的旧 vs 新关键细胞对比表 -----\n")
  print(summary_table, row.names = FALSE)

  write.csv(summary_table,
            file.path(res.path, "OLD_vs_NEW_keycell_summary.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  cat(sprintf("\n  OK 主表已保存: %s\n",
              file.path(res.path, "OLD_vs_NEW_keycell_summary.csv")))

  v17_full_sub  <- v17_full[,  c("Cell_Type", "Total_Score", "Rank")]
  v172_full_sub <- full_scoring$rank_df[, c("Cell_Type", "Total_Score", "Rank")]

  detailed <- merge(v17_full_sub, v172_full_sub,
                     by = "Cell_Type",
                     suffixes = c("_OLD_Full", "_NEW_Full"), all = TRUE)
  detailed <- detailed[order(pmin(detailed$Rank_OLD_Full,
                                   detailed$Rank_NEW_Full,
                                   na.rm = TRUE)), ]

  v17_young_sub  <- v17_young[,  c("Cell_Type", "Total_Score", "Rank")]
  v172_young_sub <- young_scoring$rank_df[, c("Cell_Type", "Total_Score", "Rank")]
  detailed_young <- merge(v17_young_sub, v172_young_sub,
                           by = "Cell_Type",
                           suffixes = c("_OLD_Young", "_NEW_Young"), all = TRUE)
  detailed_young <- detailed_young[order(pmin(detailed_young$Rank_OLD_Young,
                                              detailed_young$Rank_NEW_Young,
                                              na.rm = TRUE)), ]

  detailed_all <- merge(detailed, detailed_young, by = "Cell_Type", all = TRUE)
  write.csv(detailed_all,
            file.path(res.path, "OLD_vs_NEW_ranking_detail.csv"),
            row.names = FALSE)
  cat(sprintf("  OK 详细对比表已保存: %s\n",
              file.path(res.path, "OLD_vs_NEW_ranking_detail.csv")))

  cat("\n----- 关键细胞是否变化判定 -----\n")
  changed_full  <- v17_key_full  != key_full
  changed_young <- v17_key_young != key_young
  changed_overall <- !setequal(v17_keys, v172_keys)

  cat(sprintf("  全数据 Top 1: %s   (旧: %s, 新: %s)\n",
              ifelse(changed_full, "已改变", "未改变"), v17_key_full, key_full))
  cat(sprintf("  Young Top 1:  %s   (旧: %s, 新: %s)\n",
              ifelse(changed_young, "已改变", "未改变"), v17_key_young, key_young))
  cat(sprintf("  关键细胞集合: %s   (旧: %s; 新: %s)\n",
              ifelse(changed_overall, "已改变", "未改变"),
              paste(v17_keys, collapse=", "),
              paste(v172_keys, collapse=", ")))

  if (changed_overall) {
    cat("\n  >>> 关键细胞集合已变, 新增的关键细胞需补做 GSVA / scTenifoldKnk\n")
  } else {
    cat("\n  >>> 关键细胞集合未变, 针对关键细胞的下游分析可直接复用原结果\n")
  }

  invisible(summary_table)
}

build_compare_table()

cat("\n========== Step 5: Fig.15a 细胞比例图 ==========\n")

cell_prop <- obj@meta.data %>%
  dplyr::group_by(condition, cell_type) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(pct = n / sum(n) * 100)

all_cts <- sort(unique(cell_prop$cell_type))
ct_colors <- setNames(colorRampPalette(brewer.pal(8, "Set2"))(length(all_cts)),
                       all_cts)

p_prop <- ggplot(cell_prop,
                  aes(x = condition, y = pct, fill = cell_type)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.2) +
  scale_fill_manual(values = ct_colors) +
  labs(title = "Fig. 15a: Cell type proportions (Sham vs Stroke)",
       subtitle = sprintf("Key cells: %s", paste(key_cells, collapse = ", ")),
       x = "", y = "Proportion (%)", fill = "Cell type") +
  theme_bw(base_size = 12) +
  theme(plot.title    = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
        legend.text   = element_text(size = 9))

ggsave(file.path(fig.path, "Fig15a_CellProportions.pdf"), p_prop,
       width = 9, height = 7)
ggsave(file.path(fig.path, "Fig15a_CellProportions.png"), p_prop,
       width = 9, height = 7, dpi = 300, bg = "white")
cat("OK Fig15a 保存\n")

p_dot <- DotPlot(obj, features = biomarkers_found,
                  group.by = "cell_type",
                  cols = c("lightgrey", "#E64B35"),
                  scale = TRUE) +
  RotatedAxis() +
  labs(title = "Biomarker expression across cell types",
       subtitle = sprintf("Key cells: %s", paste(key_cells, collapse = ", "))) +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"))
ggsave(file.path(fig.path, "Fig15a_Biomarker_DotPlot.pdf"), p_dot,
       width = 8, height = 6)
ggsave(file.path(fig.path, "Fig15a_Biomarker_DotPlot.png"), p_dot,
       width = 8, height = 6, dpi = 300, bg = "white")
cat("OK DotPlot 保存\n")

cat("\n========== Step 6: GSVA + LIMMA 通路富集 ==========\n")

BUILTIN_PATHWAYS <- list(
  "Complement_and_Coagulation_Cascades" = c("C1qa", "C1qb", "C1qc", "C3", "C4a", "C4b", "C5", "C5ar1", "C5ar2", "C6", "C7", "C8a", "C8b", "C8g", "C9", "Cfb", "Cfh", "Cfi", "Cfp", "F10", "F13a1", "F2", "F3", "F5", "F7", "F8", "Masp1", "Masp2", "Plau", "Plaur", "Plg", "Serpina1a", "Serping1"),
  "Innate_Immune_System" = c("Aim2", "Casp1", "Ccl2", "Ccl3", "Ccl4", "Cxcl1", "Cxcl2", "Gsdmd", "Il18", "Il1a", "Il1b", "Il6", "Irak4", "Mefv", "Myd88", "Nfkb1", "Nfkbia", "Nlrp3", "Nod1", "Nod2", "Pycard", "Rela", "Tlr2", "Tlr4", "Tlr9", "Tnf", "Traf6"),
  "Neutrophil_Degranulation" = c("Camp", "Cd14", "Cd63", "Ctsg", "Cybb", "Defa1", "Elane", "Fcgr2b", "Fcgr3", "Itgam", "Itgax", "Lcn2", "Ltf", "Lyz1", "Lyz2", "Mmp8", "Mmp9", "Mpo", "Ncf1", "Ncf2", "Ncf4", "Ngp", "Prtn3", "Rac2", "S100a8", "S100a9", "Sell"),
  "Phagosome" = c("Atp6v0a1", "Atp6v1a", "Axl", "Cd14", "Cd36", "Cd47", "Cd9", "Ctsb", "Ctsd", "Ctsl", "Ctss", "Fcgr1", "Fcgr2b", "Fcgr3", "Gas6", "Itgam", "Itgax", "Itgb2", "Lamp1", "Lamp2", "Marco", "Mertk", "Mrc1", "Msr1", "Rab5a", "Rab7", "Sirpa", "Tyro3"),
  "NLR_Signaling_Pathway" = c("Aim2", "Atg16l1", "Atg7", "Casp1", "Casp11", "Casp4", "Gsdmd", "Hmgb1", "Il18", "Il1b", "Map3k7", "Mefv", "Nlrc4", "Nlrp1a", "Nlrp1b", "Nlrp3", "Nod1", "Nod2", "Pycard", "Ripk2", "Traf6", "Txnip"),
  "Chemokine_Signaling" = c("Ccl2", "Ccl3", "Ccl4", "Ccl5", "Ccl7", "Ccr1", "Ccr2", "Ccr5", "Cdc42", "Cxcl1", "Cxcl10", "Cxcl12", "Cxcl2", "Cxcr2", "Cxcr3", "Cxcr4", "Gnb1", "Gng2", "Jak2", "Nfkb1", "Pik3cg", "Rac1", "Rhoa", "Stat1", "Stat3"),
  "Cholesterol_Biosynthesis" = c("Cyp51", "Dhcr24", "Dhcr7", "Ebp", "Fdft1", "Fdps", "Ggps1", "Hmgcr", "Hmgcs1", "Hsd17b7", "Idi1", "Insig1", "Insig2", "Lss", "Msmo1", "Mvd", "Mvk", "Nsdhl", "Sc5d", "Scap", "Sqle", "Srebf1", "Srebf2", "Tm7sf2"),
  "Cholesterol_Efflux_and_HDL" = c("Abca1", "Abcg1", "Acat1", "Apoa1", "Apoa2", "Apoa4", "Apoe", "Cetp", "Lcat", "Ldlr", "Lrp1", "Ncr1", "Nr0b2", "Nr1h2", "Nr1h3", "Osbpl1a", "Pcsk9", "Pltp", "Ppara", "Scarb1", "Soat1", "Stard3"),
  "PPAR_Signaling" = c("Acaa1a", "Acox1", "Angptl4", "Apoa1", "Apoa2", "Apoc3", "Cd36", "Cpt1a", "Cpt1b", "Ehhadh", "Fabp1", "Fabp4", "Lpl", "Ncoa1", "Ncoa2", "Ncor1", "Ncor2", "Pck1", "Plin1", "Ppara", "Ppard", "Pparg", "Rxra", "Scd1"),
  "Fatty_Acid_Metabolism" = c("Acaa1b", "Acaa2", "Acaca", "Acacb", "Acadl", "Acadm", "Acads", "Acadvl", "Acsl1", "Acsl3", "Acsl4", "Acsl5", "Cpt1a", "Cpt2", "Crot", "Ech1", "Echs1", "Fasn", "Hadha", "Hadhb", "Hmgcs1", "Hmgcs2", "Scd1", "Scd2"),
  "Lipid_Particle_Organization" = c("Abca1", "Abcg1", "Abcg5", "Abcg8", "Apoa1", "Apoa2", "Apoa4", "Apoa5", "Apob", "Apoc1", "Apoc2", "Apoc3", "Apoc4", "Apoe", "Cav1", "Cav2", "Dgat1", "Dgat2", "Lcat", "Plin1", "Plin2", "Plin3", "Plin5", "Pltp"),
  "Sterol_Regulatory_Element_Binding" = c("Acaca", "Acly", "Elovl6", "Fasn", "Gpam", "Hmgcr", "Insig1", "Insig2", "Ldlr", "Mbtps1", "Mbtps2", "Nr1h2", "Nr1h3", "Pnpla3", "Scap", "Scd1", "Srebf1", "Srebf2"),
  "Efferocytosis_and_Apoptotic_Clearance" = c("Abca1", "Axl", "Bai1", "Cd14", "Cd36", "Dock1", "Elmo1", "Gas6", "Itga5", "Itgav", "Itgb3", "Itgb5", "Lrp1", "Marco", "Mertk", "Mfge8", "Msr1", "Pros1", "Rac1", "Rock1", "Scarb1", "Stab1", "Stab2", "Timd4", "Tyro3"),
  "Apoptosis" = c("Aifm1", "Apaf1", "Bad", "Bak1", "Bax", "Bcl2", "Bcl2l1", "Bid", "Bim", "Birc3", "Casp10", "Casp3", "Casp6", "Casp7", "Casp8", "Casp9", "Cycs", "Diablo", "Fadd", "Fas", "Fasl", "Parp1", "Puma", "Tnfrsf1a", "Tp53", "Tradd", "Xiap"),
  "ECM_Receptor_Interaction" = c("Cd36", "Cd44", "Col1a1", "Col1a2", "Col3a1", "Col4a1", "Fn1", "Itga5", "Itgav", "Itgb1", "Itgb3", "Itgb5", "Lamb1", "Lamc1", "Sdc1", "Sdc4", "Spp1", "Thbs1", "Thbs2", "Tnc", "Vtn"),
  "Oxidative_Phosphorylation" = c("Atp5a1", "Atp5b", "Atp5c1", "Atp5g1", "Atp5o", "Atp6v0a1", "Cox4i1", "Cox5a", "Cox5b", "Cox6a1", "Cycs", "Ndufa1", "Ndufa2", "Ndufb1", "Ndufb2", "Sdha", "Sdhb", "Sdhc", "Sdhd", "Uqcrc1", "Uqcrc2", "Uqcrq"),
  "Glycolysis_Gluconeogenesis" = c("Aldoa", "Aldob", "Aldoc", "Eno1", "Eno2", "Fbp1", "G6pc", "Gapdh", "Gpi1", "Hk1", "Hk2", "Hk3", "Ldha", "Ldhb", "Pck1", "Pck2", "Pcx", "Pfkl", "Pfkm", "Pfkp", "Pgam1", "Pgk1", "Pgk2", "Pklr", "Pkm", "Tpi1"),
  "HIF1_Signaling" = c("Arnt", "Bnip3", "Crebbp", "Edn1", "Egln1", "Egln2", "Egln3", "Eno1", "Ep300", "Epas1", "Epo", "Hif1a", "Hk2", "Ldha", "Nos2", "Pdk1", "Pfkl", "Pfkm", "Serpine1", "Slc2a1", "Tgfb3", "Vegfa", "Vhl"),
  "Oxidative_Stress_Response" = c("Cat", "Gclc", "Gclm", "Gpx1", "Gpx3", "Gpx4", "Hmox1", "Hmox2", "Keap1", "Nfe2l2", "Nqo1", "Nrf1", "Ppargc1a", "Prdx1", "Prdx2", "Prdx3", "Prdx4", "Sod1", "Sod2", "Sod3", "Txn1", "Txn2", "Txnrd1"),
  "Ferroptosis" = c("Acsl4", "Aifm2", "Alox15", "Atg5", "Atg7", "Cybb", "Fpn1", "Fth1", "Ftl1", "Gclc", "Gclm", "Gpx4", "Hmox1", "Lpcat3", "Map1lc3b", "Ncoa4", "Nfe2l2", "Nox1", "Slc3a2", "Slc7a11", "Tfrc"),
  "Regulation_of_Autophagy" = c("Atg10", "Atg12", "Atg16l1", "Atg3", "Atg4b", "Atg5", "Atg7", "Becn1", "Gabarap", "Gabarapl1", "Lamp1", "Map1lc3a", "Map1lc3b", "Mtor", "Nbr1", "Pik3c3", "Rptor", "Rraga", "Sqstm1", "Ulk1", "Ulk2", "Wipi1", "Wipi2"),
  "Cellular_Senescence" = c("Ccna2", "Ccnb1", "Ccnd1", "Ccne1", "Cdkn1a", "Cdkn2a", "Cdkn2b", "Cxcl1", "Cxcl2", "E2f1", "Igfbp3", "Igfbp7", "Il1b", "Il6", "Mdm2", "Mmp13", "Mmp3", "Mmp9", "Nfkb1", "Nfkb2", "Plau", "Rb1", "Rela", "Serpine1", "Sirt1", "Tp53"),
  "VEGF_Signaling" = c("Akt1", "Flt1", "Flt4", "Hras", "Kdr", "Mapk1", "Mapk3", "Nfat5", "Nos3", "Nrp1", "Nrp2", "Pik3ca", "Pik3cb", "Pik3cg", "Pla2g4a", "Ptk2", "Pxn", "Sphk1", "Src", "Vegfa", "Vegfb", "Vegfc", "Vegfd"),
  "Angiogenesis" = c("Angpt1", "Angpt2", "Dll4", "Efna1", "Epha2", "Flt1", "Jag1", "Kdr", "Notch1", "Nrp1", "Nrp2", "Pdgfb", "Pdgfrb", "Smad2", "Smad3", "Smad4", "Tek", "Tgfb1", "Tgfb2", "Tgfbr1", "Tgfbr2", "Tie1", "Vegfa", "Vegfb", "Vegfc"),
  "Leukocyte_Transendothelial_Migration" = c("Cd99", "Cdh5", "Cldn5", "Ezr", "Icam1", "Icam2", "Itga4", "Itga5", "Itgal", "Itgam", "Itgax", "Itgb1", "Itgb2", "Jam1", "Jam2", "Jam3", "Msn", "Ocln", "Pecam1", "Sele", "Sell", "Selp", "Tjp1", "Vcam1", "Vcl"),
  "TNF_Signaling" = c("Casp8", "Ccl2", "Chuk", "Cxcl1", "Fadd", "Fos", "Ikbkb", "Ikbkg", "Jun", "Mapk14", "Mapk8", "Mapk9", "Mlkl", "Nfkb1", "Nfkbia", "Rela", "Ripk1", "Ripk3", "Tnf", "Tnfrsf1a", "Tnfrsf1b", "Tradd", "Traf2", "Traf5"),
  "Regulation_of_Actin_Cytoskeleton" = c("Actb", "Actg1", "Arhgef6", "Arpc1a", "Arpc2", "Arpc3", "Arpc4", "Cdc42", "Cfl1", "Cfl2", "Limk1", "Mapk1", "Pak1", "Pak2", "Pak4", "Pip5k1a", "Rac1", "Rhoa", "Rhoc", "Rock1", "Rock2", "Was", "Wasf2", "Wasl"),
  "Toll_Like_Receptor_Signaling" = c("Chuk", "Ikbkb", "Ikbke", "Ikbkg", "Irak1", "Irak4", "Map3k7", "Myd88", "Nfkb1", "Rela", "Tab1", "Tab2", "Tbk1", "Ticam1", "Ticam2", "Tlr1", "Tlr2", "Tlr3", "Tlr4", "Tlr5", "Tlr6", "Tlr7", "Tlr8", "Tlr9", "Traf3", "Traf6"),
  "Interferon_Signaling" = c("Ifit1", "Ifit2", "Ifit3", "Ifna1", "Ifnar1", "Ifnar2", "Ifnb1", "Ifng", "Ifngr1", "Ifngr2", "Irf1", "Irf3", "Irf7", "Irf9", "Isg15", "Isg20", "Jak1", "Jak2", "Mx1", "Mx2", "Oas1a", "Oas2", "Oas3", "Stat1", "Stat2", "Tyk2"),
  "IL6_JAK_STAT3" = c("Bcl2", "Ccnd1", "Fos", "Grb2", "Hras", "Il6", "Il6ra", "Il6st", "Jak1", "Jak2", "Jun", "Mapk1", "Mapk3", "Myc", "Ptpn1", "Ptpn11", "Ptpn2", "Socs1", "Socs3", "Sos1", "Stat3", "Tyk2", "Vegfa"),
  "MAPK_Signaling" = c("Atf2", "Elk1", "Fos", "Jun", "Map2k1", "Map2k2", "Map2k3", "Map2k4", "Map2k6", "Map2k7", "Map3k1", "Map3k2", "Map3k5", "Map3k7", "Mapk1", "Mapk10", "Mapk11", "Mapk12", "Mapk14", "Mapk3", "Mapk8", "Mapk9", "Myc"),
  "DNA_Damage_Response" = c("Apex1", "Atm", "Atr", "Brca1", "Brca2", "Chek1", "Chek2", "Fen1", "Lig1", "Lig3", "Lig4", "Mdm2", "Mgmt", "Mre11a", "Nbn", "Neil1", "Ogg1", "Parp1", "Pcna", "Rad50", "Rad51", "Tp53", "Xrcc1", "Xrcc2", "Xrcc3", "Xrcc4"),
  "Cell_Cycle" = c("Ccna1", "Ccna2", "Ccnb1", "Ccnb2", "Ccnd1", "Ccnd2", "Ccne1", "Ccne2", "Cdk1", "Cdk2", "Cdk4", "Cdk6", "Cdkn1a", "Cdkn1b", "Cdkn2a", "Cdkn2b", "E2f1", "E2f2", "Mcm2", "Mcm3", "Mcm4", "Mcm5", "Mcm6", "Mcm7", "Myc", "Rb1"),
  "Adipocytokine_Signaling" = c("Acaca", "Adipoq", "Cpt1a", "Ikbkg", "Il6", "Irs1", "Irs2", "Lep", "Lepr", "Mapk8", "Mapk9", "Nfkb1", "Pik3r1", "Ppargc1a", "Prkaa1", "Prkaa2", "Slc2a4", "Socs3", "Stat3", "Tnf"),
  "Insulin_Signaling" = c("Akt1", "Akt2", "Akt3", "Foxo1", "Foxo3", "Glut4", "Gsk3a", "Gsk3b", "Hk2", "Insr", "Irs1", "Irs2", "Mapk1", "Mapk3", "Mtor", "Pdk1", "Pik3ca", "Pik3r1", "Ppp1r3a", "Prkci", "Pygl", "Rptor", "Slc2a1", "Slc2a4"),
  "T_Cell_Receptor_Signaling" = c("Cd247", "Cd3d", "Cd3e", "Cd3g", "Cdc42", "Fyn", "Grb2", "Hras", "Il17a", "Il2", "Il4", "Lat", "Lck", "Mapk1", "Mapk3", "Nck1", "Nfat5", "Nfatc1", "Nfatc2", "Pak1", "Plcg1", "Rac1", "Rhoa", "Sos1", "Vav1", "Zap70"),
  "Pyruvate_Metabolism" = c("Acot11", "Acss1", "Acss2", "Acyp1", "Dlat", "Dld", "Gpi1", "Hagh", "Ldha", "Ldhb", "Ldhc", "Me1", "Me2", "Pck1", "Pck2", "Pcx", "Pdha1", "Pdhb", "Pdhx", "Pklr", "Pkm"),
  "Citrate_Cycle_TCA_Cycle" = c("Aco1", "Aco2", "Cs", "Dld", "Dlst", "Fh1", "Idh1", "Idh2", "Idh3a", "Idh3b", "Idh3g", "Mdh1", "Mdh2", "Ogdh", "Ogdhl", "Pcx", "Sdha", "Sdhb", "Sdhc", "Sdhd", "Sucla2", "Suclg1", "Suclg2"),
  "Arachidonic_Acid_Metabolism" = c("Alox12", "Alox15", "Alox5", "Cyp2c44", "Cyp4a10", "Cyp4a12a", "Ephx2", "Pla2g4a", "Pla2g4b", "Pla2g4c", "Pla2g5", "Pla2g6", "Pld3", "Ptges", "Ptges2", "Ptges3", "Ptgs1", "Ptgs2"),
  "Glutathione_Metabolism" = c("Anpep", "Gclc", "Gclm", "Ggt1", "Ggt5", "Gpx1", "Gpx2", "Gpx3", "Gpx4", "Gpx5", "Gpx7", "Gpx8", "Gsr", "Gss", "Gsta1", "Gsta2", "Gstm1", "Gstm2", "Gstp1", "Gstp2", "Mgst1", "Mgst2", "Mgst3"),
  "Sphingolipid_Metabolism" = c("Asah1", "Asah2", "Cerk", "Cers2", "Cers4", "Cers5", "Cers6", "Degs1", "Degs2", "Sgms1", "Sgms2", "Sgpp1", "Sgpp2", "Smpd1", "Smpd2", "Smpd3", "Smpd4", "Sphk1", "Sphk2", "Sptlc1", "Sptlc2", "Sptssa"),
  "Glycerophospholipid_Metabolism" = c("Cept1", "Chpt1", "Etnk1", "Gpd1", "Gpd2", "Lpin1", "Lpin2", "Lpin3", "Pcyt1a", "Pcyt2", "Pemt", "Plcb1", "Plcb2", "Plcb3", "Plcg1", "Plcg2", "Plch1", "Pld1", "Pld2", "Pld3", "Pnpla2", "Pnpla3"),
  "mTOR_Signaling" = c("Akt1", "Atg13", "Deptor", "Eif4ebp1", "Lamtor1", "Lamtor2", "Mlst8", "Mtor", "Rheb", "Rps6", "Rps6ka1", "Rps6kb1", "Rptor", "Rraga", "Rragb", "Rragc", "Rragd", "Sesn1", "Sesn2", "Tsc1", "Tsc2", "Ulk1"),
  "P53_Signaling" = c("Bak1", "Bax", "Bbc3", "Bcl2", "Bcl2l1", "Bid", "Cdkn1a", "Cdkn2a", "Fas", "Gadd45a", "Gadd45b", "Gadd45g", "Mdm2", "Mdm4", "Pmaip1", "Pten", "Sfn", "Tp53", "Trp53bp1", "Trp53inp1", "Zmat3"),
  "Oncostatin_M_Signaling" = c("Grb2", "Hras", "Il6st", "Jak1", "Jak2", "Lifr", "Mapk1", "Mapk3", "Osm", "Osmr", "Ptpn11", "Shp2", "Socs3", "Sos1", "Stat1", "Stat3", "Stat5a", "Stat5b", "Tyk2"),
  "TGFbeta_Signaling" = c("Acvr1", "Acvr2a", "Bmp2", "Bmp4", "Bmpr1a", "Bmpr2", "Ppp1ca", "Ppp1cb", "Ski", "Skil", "Smad2", "Smad3", "Smad4", "Smad6", "Smad7", "Tgfb1", "Tgfb2", "Tgfb3", "Tgfbr1", "Tgfbr2", "Tgfbr3", "Tgfbrap1"),
  "Notch_Signaling" = c("Adam10", "Adam17", "Dll1", "Dll3", "Dll4", "Dtx1", "Hes1", "Hes5", "Hey1", "Hey2", "Jag1", "Jag2", "Maml1", "Maml2", "Ncstn", "Notch1", "Notch2", "Notch3", "Notch4", "Psen1", "Psen2", "Rbpj"),
  "Wnt_Signaling" = c("Apc", "Axin1", "Ccnd1", "Ctnnb1", "Dvl1", "Dvl2", "Dvl3", "Fzd1", "Fzd2", "Fzd4", "Fzd7", "Gsk3b", "Lef1", "Lrp5", "Lrp6", "Myc", "Tcf7", "Wif1", "Wnt1", "Wnt3", "Wnt3a", "Wnt5a", "Wnt7a", "Wnt7b")
)

cat(sprintf("内嵌通路集: %d 条通路, 平均 %.1f 基因/通路\n",
            length(BUILTIN_PATHWAYS),
            mean(sapply(BUILTIN_PATHWAYS, length))))

pathway_list <- BUILTIN_PATHWAYS
pathway_list <- pathway_list[sapply(pathway_list, length) >= 5 &
                              sapply(pathway_list, length) <= 500]
cat(sprintf("过滤后保留 %d 条通路 (基因数 5-500)\n", length(pathway_list)))

all_pathway_genes <- unique(unlist(pathway_list))
n_match <- sum(all_pathway_genes %in% rownames(obj))
cat(sprintf("基因匹配率: %d / %d (%.1f%%) 在数据中\n",
            n_match, length(all_pathway_genes),
            100 * n_match / length(all_pathway_genes)))

run_gsva_for_celltype <- function(seu, cell_type_name, pathway_list,
                                    tag = "", fig_prefix = "Fig15b") {
  cat(sprintf("\n----- GSVA: %s (%s 数据) -----\n", cell_type_name, tag))

  key_obj <- subset(seu, cell_type == cell_type_name)
  cat(sprintf("  细胞数: %d\n", ncol(key_obj)))

  if (ncol(key_obj) < 30) {
    cat("  细胞数 < 30, 跳过 GSVA\n")
    return(NULL)
  }

  MAX_CELLS_GSVA <- 3000
  if (ncol(key_obj) > MAX_CELLS_GSVA) {
    cat(sprintf("  下采样 %d -> %d\n", ncol(key_obj), MAX_CELLS_GSVA))
    set.seed(42)
    sham_cells   <- WhichCells(key_obj, expression = condition == "Sham")
    stroke_cells <- WhichCells(key_obj, expression = condition == "Stroke")
    n_sham   <- min(round(MAX_CELLS_GSVA * length(sham_cells) /
                           (length(sham_cells) + length(stroke_cells))),
                     length(sham_cells))
    n_stroke <- MAX_CELLS_GSVA - n_sham
    sampled <- c(sample(sham_cells, n_sham),
                  sample(stroke_cells, min(n_stroke, length(stroke_cells))))
    key_obj <- subset(key_obj, cells = sampled)
  }

  key_expr <- as.matrix(GetAssayData(key_obj, assay = "RNA", layer = "data"))
  key_expr <- key_expr[rowSums(key_expr > 0) >= 10, , drop = FALSE]
  cat(sprintf("  GSVA 输入: %d 基因 x %d 细胞\n", nrow(key_expr), ncol(key_expr)))
  cat(sprintf("  Sham: %d, Stroke: %d\n",
              sum(key_obj$condition == "Sham"),
              sum(key_obj$condition == "Stroke")))

  t_start <- Sys.time()
  gsva_scores <- tryCatch({
    param <- gsvaParam(exprData = key_expr, geneSets = pathway_list,
                        kcdf = "Gaussian", minSize = 5, maxSize = 500)
    gsva(param, verbose = FALSE)
  }, error = function(e) {
    cat("  新版 API 失败:", e$message, "\n"); NULL
  })

  if (is.null(gsva_scores)) {
    gsva_scores <- tryCatch({
      gsva(expr = key_expr, gset.idx.list = pathway_list,
           method = "gsva", kcdf = "Gaussian",
           min.sz = 5, max.sz = 500, verbose = FALSE, parallel.sz = 4)
    }, error = function(e) { cat("  旧版 API 失败:", e$message, "\n"); NULL })
  }

  if (is.null(gsva_scores)) { cat("  GSVA 失败\n"); return(NULL) }

  cat(sprintf("  GSVA 完成, %.1f 秒, 结果 %d 通路 x %d 细胞\n",
              as.numeric(Sys.time() - t_start, units = "secs"),
              nrow(gsva_scores), ncol(gsva_scores)))

  condition_vec <- factor(key_obj$condition[colnames(gsva_scores)],
                            levels = c("Sham", "Stroke"))
  design <- model.matrix(~0 + condition_vec)
  colnames(design) <- levels(condition_vec)

  fit <- lmFit(gsva_scores, design)
  contrast_matrix <- makeContrasts(StrokeVsSham = Stroke - Sham, levels = design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)

  diff_result <- topTable(fit2, coef = "StrokeVsSham", number = Inf, sort.by = "t")
  diff_result$Pathway <- rownames(diff_result)

  significant <- diff_result[abs(diff_result$t) >= 2 & diff_result$P.Value < 0.05, ]
  cat(sprintf("  显著通路 (|t|>=2, P<0.05): %d 条\n", nrow(significant)))

  if (nrow(significant) < 5) {
    significant <- diff_result[diff_result$P.Value < 0.05, ]
    cat(sprintf("  放宽 P<0.05: %d 条\n", nrow(significant)))
  }
  if (nrow(significant) < 5) {
    diff_result_sorted <- diff_result[order(abs(diff_result$t), decreasing = TRUE), ]
    significant <- head(diff_result_sorted, 10)
    cat("  使用 Top 10 by |t|\n")
  }

  significant <- significant[order(abs(significant$t), decreasing = TRUE), ]
  top10 <- head(significant, 10)

  cat("  Top 10 通路:\n")
  for (i in seq_len(nrow(top10))) {
    cat(sprintf("    %d. %s (t=%.2f, P=%.3g)\n",
                i, top10$Pathway[i], top10$t[i], top10$P.Value[i]))
  }

  ct_clean <- gsub("[^A-Za-z0-9]", "_", cell_type_name)
  tag_clean <- if (nchar(tag) > 0) paste0("_", tag) else ""

  write.csv(diff_result,
            file.path(res.path, sprintf("gsva_limma_all_%s%s.csv", ct_clean, tag_clean)),
            row.names = FALSE)
  write.csv(top10,
            file.path(res.path, sprintf("gsva_top10_%s%s.csv", ct_clean, tag_clean)),
            row.names = FALSE)
  saveRDS(gsva_scores,
          file.path(res.path, sprintf("gsva_scores_%s%s.rds", ct_clean, tag_clean)))

  if (nrow(top10) > 0) {
    top_pathways <- top10$Pathway
    heatmap_data <- gsva_scores[top_pathways, , drop = FALSE]

    condition_ordered <- order(condition_vec)
    heatmap_data <- heatmap_data[, condition_ordered, drop = FALSE]
    condition_sorted <- condition_vec[condition_ordered]

    annotation_col <- data.frame(
      Condition = condition_sorted,
      row.names = colnames(heatmap_data)
    )
    ann_colors <- list(Condition = c(Sham = "#45a9b8", Stroke = "#f76a56"))

    clean_names <- gsub("_", " ", top_pathways)
    clean_names <- stringr::str_trunc(clean_names, 50)
    rownames(heatmap_data) <- clean_names

    fig_base <- sprintf("%s_%s%s_Heatmap", fig_prefix, ct_clean, tag_clean)
    main_title <- sprintf("%s: %s (%s) - Top pathways",
                           fig_prefix, cell_type_name,
                           ifelse(nchar(tag) > 0, tag, "Full"))

    pdf(file.path(fig.path, paste0(fig_base, ".pdf")), width = 12, height = 6)
    pheatmap(heatmap_data, scale = "row",
             color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
             annotation_col = annotation_col, annotation_colors = ann_colors,
             cluster_cols = FALSE, cluster_rows = TRUE,
             show_colnames = FALSE, fontsize_row = 10,
             main = main_title, border_color = NA)
    dev.off()

    png(file.path(fig.path, paste0(fig_base, ".png")),
        width = 12, height = 6, units = "in", res = 300)
    pheatmap(heatmap_data, scale = "row",
             color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
             annotation_col = annotation_col, annotation_colors = ann_colors,
             cluster_cols = FALSE, cluster_rows = TRUE,
             show_colnames = FALSE, fontsize_row = 10,
             main = main_title, border_color = NA)
    dev.off()

    cat(sprintf("  OK 热图保存: %s.pdf/png\n", fig_base))
  }

  list(gsva_scores = gsva_scores, diff_result = diff_result, top10 = top10)
}

cat(sprintf("\n>>> 关键细胞共 %d 个: %s\n",
            length(key_cells), paste(key_cells, collapse = ", ")))

gsva_results <- list()

cat(sprintf("\n>>> [Full] 关键细胞 GSVA: %s\n", key_full))
gsva_results[[paste0(key_full, "_Full")]] <-
  run_gsva_for_celltype(obj, key_full, pathway_list,
                          tag = "Full", fig_prefix = "Fig15b")

if (key_young != key_full) {
  cat(sprintf("\n>>> [Young] 关键细胞 GSVA: %s\n", key_young))
  gsva_results[[paste0(key_young, "_Young")]] <-
    run_gsva_for_celltype(obj_young, key_young, pathway_list,
                            tag = "Young", fig_prefix = "Fig15c")
} else {
  cat(sprintf("\n>>> Young Top 1 与全数据一致 (%s), 无需重复\n", key_young))
}

remaining <- setdiff(key_cells, c(key_full, key_young))
if (length(remaining) > 0) {
  for (extra_ct in remaining) {
    cat(sprintf("\n>>> [Full, 额外] 关键细胞 GSVA: %s\n", extra_ct))
    gsva_results[[paste0(extra_ct, "_Full_extra")]] <-
      run_gsva_for_celltype(obj, extra_ct, pathway_list,
                              tag = "Full", fig_prefix = "Fig15d")
  }
}

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("    第十八部分 完成\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("关键细胞:\n")
cat(sprintf("  · 全数据 Top 1: %s (%d 个细胞)\n",
            key_full, sum(obj$cell_type == key_full)))
cat(sprintf("  · Young 子集 Top 1: %s (%d 个细胞, Young)\n",
            key_young, sum(obj_young$cell_type == key_young)))
if (length(key_cells) == 1) {
  cat("  · 两者一致, 最终 1 个关键细胞\n")
} else {
  cat(sprintf("  · 合并后最终 %d 个关键细胞: %s\n",
              length(key_cells), paste(key_cells, collapse = ", ")))
}

cat("\nGSVA 完成的关键细胞:\n")
for (k in names(gsva_results)) {
  if (!is.null(gsva_results[[k]])) {
    cat(sprintf("  · %s: Top 10 通路已生成\n", k))
  }
}

cat("\n输出文件:\n")
cat("  Results/\n")
cat("    keycell_ranking.csv               (全数据 ranking)\n")
cat("    keycell_ranking_YOUNG.csv         (Young ranking)\n")
cat("    comparison_Full_vs_Young.csv      (新体系内部对比)\n")
cat("    OLD_vs_NEW_keycell_summary.csv    (主管要求: 旧 vs 新关键细胞对比)\n")
cat("    OLD_vs_NEW_ranking_detail.csv     (详细 ranking 对比)\n")
cat("    key_cells.txt                     (最终关键细胞)\n")
cat("    biomarker_*.csv                   (4 个: avg/pct/log2FC/pvalue)\n")
cat("    gsva_top10_*.csv                  (每个关键细胞的 Top 10 通路)\n")
cat("    gsva_limma_all_*.csv              (每个关键细胞的全部通路结果)\n")
cat("    gsva_scores_*.rds                 (GSVA 原始打分矩阵)\n")
cat("  Figures/\n")
cat("    Fig15a_CellProportions.pdf/png\n")
cat("    Fig15a_Biomarker_DotPlot.pdf/png\n")
cat("    Fig15b/c/d_*_Heatmap.pdf/png      (每个关键细胞的 Top 10 通路热图)\n")

cat("\nOK 第十八部分 完成\n")

sink(file.path(res.path, "session_info.txt"))
sessionInfo()
sink()
