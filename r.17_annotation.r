
cat("\n========================================\n")
cat("17. 单细胞细胞类型注释\n")
cat("========================================\n\n")

rm(list = ls())
options(bitmapType = "cairo")
options(device = "pdf")

setwd("")

output_dir <- "./17_annotation"
res.path   <- file.path(output_dir, "Results")
fig.path   <- file.path(output_dir, "Figures")
for (p in c(output_dir, res.path, fig.path)) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE)
}

cat("========== 加载 R 包 ==========\n")

required_packages <- c(
  "Seurat", "dplyr", "tidyr", "ggplot2", "pheatmap",
  "RColorBrewer", "patchwork", "SingleR", "SummarizedExperiment",
  "celldex"
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
  stop(sprintf("缺少 R 包: %s", paste(missing_pkgs, collapse=", ")))
}

cat("\n========== 读取 Seurat 对象 ==========\n")

seurat_file <- "./16_scRNAseq/Data/seurat_obj_clustered.rds"
if (!file.exists(seurat_file)) stop("seurat 对象不存在: ", seurat_file)

obj <- readRDS(seurat_file)
cat("细胞数:", ncol(obj), "\n")
cat("基因数:", nrow(obj), "\n")
cat("Clusters:", length(unique(obj$seurat_clusters)), "\n")

Idents(obj) <- "seurat_clusters"
cluster_levels <- sort(unique(obj$seurat_clusters))

cat("\n========== 定义 marker 库 ==========\n")

marker_libs <- list(
  Endo         = c("Cldn5", "Flt1", "Vwf"),
  Micro        = c("C1qc", "Csf1r", "P2ry12", "Tmem119"),
  Astro        = c("Slc1a2", "Aqp4", "Aldoc", "Aldh1l1"),
  CP           = c("Clic6"),
  SMC          = c("Tagln", "Acta2"),
  Pericyte     = c("Pdgfrb", "Rgs5", "Notch3"),
  VLMC         = c("Dcn", "Col1a1", "Col1a2"),
  Macro        = c("Lyz2", "Cd74", "Mrc1"),
  Mono         = c("Ly6c2", "Ccr2", "Plac8"),
  Oligo        = c("Mog", "Mbp"),
  Neut         = c("S100a8", "S100a9"),
  NK_T         = c("Cd3d", "Nkg7", "Gzmb"),
  B_cell       = c("Cd19", "Ms4a1", "Cd79a"),
  DC           = c("H2-Ab1", "H2-Aa", "Flt3"),
  Epend        = c("Dnah9", "Cfap54"),
  Prolif_Micro = c("Mki67", "Top2a", "Ccnb1")
)

dotplot_marker_order <- c(

  "Cldn5", "Flt1",

  "Rgs5",

  "P2ry12", "Tmem119",

  "Mki67", "Top2a",

  "Mrc1",

  "Ly6c2", "Ccr2", "Plac8",

  "Flt3", "H2-Ab1",

  "S100a8", "S100a9",

  "Cd3d", "Nkg7",

  "Cd19", "Ms4a1", "Cd79a"
)

celltype_display_order <- c(
  "Endo", "Pericyte", "SMC", "VLMC",
  "Astro", "CP", "Epend", "Oligo",
  "Micro", "Prolif_Micro", "Macro", "Mono", "DC",
  "Neut", "NK_T", "B_cell"
)

cat("已定义", length(marker_libs), "个细胞类型 marker 库\n")
cat("总 marker 数(去重):", length(unique(unlist(marker_libs))), "\n")

marker_libs_filtered <- lapply(marker_libs, function(genes) {
  intersect(genes, rownames(obj))
})
marker_libs_filtered <- marker_libs_filtered[sapply(marker_libs_filtered, length) >= 1]
cat("过滤后可用库数:", length(marker_libs_filtered), "\n")

cat("\n各 marker 库可用基因数:\n")
for (ct in names(marker_libs_filtered)) {
  n_total <- length(marker_libs[[ct]])
  n_avail <- length(marker_libs_filtered[[ct]])
  cat(sprintf("  %-13s %2d / %2d 可用  [%s]\n", ct, n_avail, n_total,
              paste(marker_libs_filtered[[ct]], collapse = ", ")))
}

cat("\n========== 写出 marker 文献溯源表 ==========\n")

marker_source_file <- "./17_annotation/marker_gene_sources.csv"

if (file.exists(marker_source_file)) {
  marker_src_df <- read.csv(marker_source_file, stringsAsFactors = FALSE)
  cat("  从", marker_source_file, "读入", nrow(marker_src_df), "个细胞类型的溯源信息\n")

  marker_src_df$N_Markers_Total <- sapply(marker_src_df$CellType, function(ct) {
    if (ct %in% names(marker_libs)) length(marker_libs[[ct]]) else NA
  })
  marker_src_df$N_Markers_Available <- sapply(marker_src_df$CellType, function(ct) {
    if (ct %in% names(marker_libs_filtered)) length(marker_libs_filtered[[ct]]) else 0
  })
  marker_src_df$Available_Markers <- sapply(marker_src_df$CellType, function(ct) {
    if (ct %in% names(marker_libs_filtered))
      paste(marker_libs_filtered[[ct]], collapse = "; ")
    else ""
  })

  write.csv(marker_src_df,
            file.path(res.path, "marker_source_table.csv"),
            row.names = FALSE)
  cat("OK 已保存 Results/marker_source_table.csv\n")
} else {
  cat("  未找到", marker_source_file, ", 生成简化版溯源表\n")
  fallback_df <- data.frame(
    CellType = names(marker_libs),
    Markers  = sapply(marker_libs, function(g) paste(g, collapse = "; ")),
    N_Markers_Total = sapply(marker_libs, length),
    N_Markers_Available = sapply(names(marker_libs), function(ct) {
      if (ct %in% names(marker_libs_filtered)) length(marker_libs_filtered[[ct]]) else 0
    }),
    Available_Markers = sapply(names(marker_libs), function(ct) {
      if (ct %in% names(marker_libs_filtered))
        paste(marker_libs_filtered[[ct]], collapse = "; ")
      else ""
    }),
    stringsAsFactors = FALSE
  )
  write.csv(fallback_df, file.path(res.path, "marker_source_table.csv"), row.names = FALSE)
  cat("OK 已保存 Results/marker_source_table.csv\n")
}

cat("\n========== Step 1: Marker-based 打分 ==========\n")

markers_cache_primary <- file.path(res.path, "all_cluster_markers.rds")

if (file.exists(markers_cache_primary)) {
  cat("使用本目录缓存:", markers_cache_primary, "\n")
  all_markers <- readRDS(markers_cache_primary)
} else {
  cat("运行 FindAllMarkers (耗时 10-30 min)...\n")
  all_markers <- FindAllMarkers(obj, only.pos = TRUE,
                                min.pct = 0.25, logfc.threshold = 0.25,
                                max.cells.per.ident = 500)
  saveRDS(all_markers, markers_cache_primary)
}

cat("总 DE 基因:", nrow(all_markers), "\n")

score_cluster_library <- function(markers_df, library_genes) {
  lib_markers <- markers_df[markers_df$gene %in% library_genes, , drop = FALSE]
  if (nrow(lib_markers) == 0) return(0)
  weighted_sum <- sum(lib_markers$avg_log2FC * lib_markers$pct.1)
  weighted_sum / sqrt(length(library_genes))
}

score_matrix <- matrix(0,
                       nrow = length(marker_libs_filtered),
                       ncol = length(cluster_levels),
                       dimnames = list(names(marker_libs_filtered),
                                        paste0("C", cluster_levels)))

for (cl in cluster_levels) {
  cl_markers <- all_markers[all_markers$cluster == cl, ]
  for (lib_name in names(marker_libs_filtered)) {
    score_matrix[lib_name, paste0("C", cl)] <-
      score_cluster_library(cl_markers, marker_libs_filtered[[lib_name]])
  }
}

write.csv(score_matrix, file.path(res.path, "marker_score_matrix.csv"))

cat("\n----- Marker-based 打分 (每个 cluster Top 3) -----\n")
marker_top_calls <- data.frame(
  Cluster = paste0("C", cluster_levels),
  Top1 = NA_character_, Top1_score = NA_real_,
  Top2 = NA_character_, Top2_score = NA_real_,
  Top3 = NA_character_, Top3_score = NA_real_,
  stringsAsFactors = FALSE
)
for (i in seq_along(cluster_levels)) {
  cl_col <- paste0("C", cluster_levels[i])
  scores <- score_matrix[, cl_col]
  ord <- order(scores, decreasing = TRUE)
  marker_top_calls$Top1[i] <- names(scores)[ord[1]]
  marker_top_calls$Top1_score[i] <- round(scores[ord[1]], 3)
  if (length(ord) >= 2) {
    marker_top_calls$Top2[i] <- names(scores)[ord[2]]
    marker_top_calls$Top2_score[i] <- round(scores[ord[2]], 3)
  }
  if (length(ord) >= 3) {
    marker_top_calls$Top3[i] <- names(scores)[ord[3]]
    marker_top_calls$Top3_score[i] <- round(scores[ord[3]], 3)
  }
}
print(marker_top_calls)
write.csv(marker_top_calls,
          file.path(res.path, "marker_top_calls.csv"), row.names = FALSE)

cat("\n========== Step 2: SingleR 双 reference 注释 ==========\n")

expr_for_singleR <- GetAssayData(obj, assay = "RNA", layer = "data")

mouse_ref <- tryCatch(celldex::MouseRNAseqData(),
                      error = function(e) { cat("  MouseRNAseqData 失败:", e$message, "\n"); NULL })
immgen_ref <- tryCatch(celldex::ImmGenData(),
                       error = function(e) { cat("  ImmGenData 失败:", e$message, "\n"); NULL })

run_singleR_per_cluster <- function(ref_data, ref_name, expr, cluster_vec) {
  if (is.null(ref_data)) return(NULL)
  cat(sprintf("  运行 SingleR (%s)...\n", ref_name))
  pred <- tryCatch(
    SingleR(test = expr, ref = ref_data, labels = ref_data$label.main,
            clusters = cluster_vec, fine.tune = TRUE),
    error = function(e) { cat(sprintf("  SingleR %s 失败: %s\n", ref_name, e$message)); NULL }
  )
  if (is.null(pred)) return(NULL)
  data.frame(Cluster = rownames(pred),
             Label = pred$labels,
             Pruned_Label = pred$pruned.labels,
             stringsAsFactors = FALSE)
}

singleR_mouse  <- run_singleR_per_cluster(mouse_ref, "MouseRNAseqData",
                                           expr_for_singleR,
                                           as.character(obj$seurat_clusters))
singleR_immgen <- run_singleR_per_cluster(immgen_ref, "ImmGenData",
                                           expr_for_singleR,
                                           as.character(obj$seurat_clusters))

if (!is.null(singleR_mouse))
  write.csv(singleR_mouse,
            file.path(res.path, "singleR_MouseRNAseqData.csv"), row.names = FALSE)
if (!is.null(singleR_immgen))
  write.csv(singleR_immgen,
            file.path(res.path, "singleR_ImmGenData.csv"), row.names = FALSE)

cat("\n========== Step 3: 三方法交叉验证 + Final annotation ==========\n")

final_calls <- data.frame(
  Cluster = paste0("C", cluster_levels),
  N_cells = sapply(cluster_levels, function(cl) sum(obj$seurat_clusters == cl)),
  Marker_Top1  = marker_top_calls$Top1,
  Marker_Score = marker_top_calls$Top1_score,
  SingleR_Mouse  = NA_character_,
  SingleR_ImmGen = NA_character_,
  Final_Call = NA_character_,
  Confidence = NA_character_,
  stringsAsFactors = FALSE
)

if (!is.null(singleR_mouse)) {
  m_map <- setNames(singleR_mouse$Label, paste0("C", singleR_mouse$Cluster))
  final_calls$SingleR_Mouse <- m_map[final_calls$Cluster]
}
if (!is.null(singleR_immgen)) {
  i_map <- setNames(singleR_immgen$Label, paste0("C", singleR_immgen$Cluster))
  final_calls$SingleR_ImmGen <- i_map[final_calls$Cluster]
}

normalize_name <- function(x) {
  if (is.na(x) || nchar(x) == 0) return("")
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[[:space:]_-]+", "", x)
  if (x == "endo" || grepl("endothelial|endothelium", x)) return("endothelial")
  if (x %in% c("micro", "prolifmicro") || grepl("microglia", x)) return("microglia")
  if (x == "astro" || grepl("astrocyte", x)) return("astrocyte")
  if (x == "cp" || grepl("choroid|epithel", x)) return("cp_epithelial")
  if (x == "smc" || grepl("smoothmuscle|smc", x)) return("smc")
  if (x %in% c("macro", "mono") || grepl("macrophage|monocyte", x)) return("myeloid")
  if (x == "oligo" || grepl("oligodendrocyte|oligo", x)) return("oligo")
  if (x == "neut" || grepl("neutrophil|granulo", x)) return("neutrophil")
  if (x == "nkt" || grepl("nkcell|tcell|tcells|lymphocyte|nk", x)) return("lymphocyte")
  if (x == "bcell" || grepl("bcell|bcells", x)) return("bcell")
  if (x == "dc" || grepl("dendritic", x)) return("dc")
  if (x == "epend" || grepl("ependym", x)) return("ependymal")
  if (x == "pericyte" || grepl("pericyte", x)) return("pericyte")
  if (x == "vlmc" || grepl("fibroblast|stromal|leptomenin", x)) return("fibroblast")
  x
}

for (i in seq_len(nrow(final_calls))) {
  n_m    <- normalize_name(final_calls$Marker_Top1[i])
  n_sr_m <- normalize_name(final_calls$SingleR_Mouse[i])
  n_sr_i <- normalize_name(final_calls$SingleR_ImmGen[i])

  agreements <- sum(c(n_m == n_sr_m, n_m == n_sr_i, n_sr_m == n_sr_i), na.rm = TRUE)

  final_calls$Final_Call[i] <- final_calls$Marker_Top1[i]

  if (n_m == n_sr_m && n_m == n_sr_i && nchar(n_m) > 0) {
    final_calls$Confidence[i] <- "High (all 3 agree)"
  } else if (agreements >= 1) {
    final_calls$Confidence[i] <- "Medium (2/3 agree)"
  } else {
    final_calls$Confidence[i] <- "Low (inconsistent)"
  }

  cl <- cluster_levels[i]
  cells_in_cl <- WhichCells(obj, idents = cl)
  if ("percent.mt" %in% colnames(obj@meta.data)) {
    mean_mt <- mean(obj$percent.mt[cells_in_cl], na.rm = TRUE)
    if (mean_mt > 10) {
      final_calls$Final_Call[i] <- "Low_Quality"
      final_calls$Confidence[i] <- sprintf("QC-fail (mt=%.1f%%)", mean_mt)
    }
  }
}

cat("\n----- Final annotation -----\n")
print(final_calls)
write.csv(final_calls,
          file.path(res.path, "final_annotation.csv"), row.names = FALSE)

cell_to_call <- setNames(final_calls$Final_Call, paste0("C", cluster_levels))
obj$cell_type <- unname(cell_to_call[paste0("C", obj$seurat_clusters)])
obj$cell_type_confidence <- unname(
  setNames(final_calls$Confidence, paste0("C", cluster_levels))[
    paste0("C", obj$seurat_clusters)])

obj$is_high_quality <- obj$cell_type != "Low_Quality"

cat("\n----- 注释分布 (含 QC 标记) -----\n")
print(table(obj$cell_type))
cat(sprintf("\n高质量细胞: %d / %d (%.1f%%)\n",
            sum(obj$is_high_quality), ncol(obj),
            100 * sum(obj$is_high_quality) / ncol(obj)))

gsm_to_group <- function(gsm) {
  num <- as.numeric(gsub("GSM", "", gsm))
  if (is.na(num)) return(NA)
  if (num <= 7060820) return("Sham") else return("Stroke")
}
obj$condition <- sapply(obj$orig.ident, gsm_to_group)

cat("\n========== 永久剔除 Low_Quality 细胞 ==========\n")

n_before <- ncol(obj)
obj_clean <- subset(obj, subset = is_high_quality == TRUE)
n_after <- ncol(obj_clean)

cat(sprintf("剔除前: %d -> 剔除后: %d (剔除 %d, %.1f%%)\n",
            n_before, n_after, n_before - n_after,
            100 * (n_before - n_after) / n_before))

saveRDS(obj_clean, file.path(output_dir, "seurat_obj_annotated_clean.rds"))
cat("OK 高质量注释对象已保存: seurat_obj_annotated_clean.rds\n")

cat("\n========== Step 4: 可视化 ==========\n")

p_umap_clean <- DimPlot(obj_clean, reduction = "umap", group.by = "cell_type",
                        label = TRUE, label.size = 4, repel = TRUE) +
  ggtitle("UMAP - Cell types") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 8))
ggsave(file.path(fig.path, "UMAP_celltype_clean.pdf"), p_umap_clean,
       width = 10, height = 8)
ggsave(file.path(fig.path, "UMAP_celltype_clean.png"), p_umap_clean,
       width = 10, height = 8, dpi = 200)

p_umap_cluster <- DimPlot(obj_clean, reduction = "umap",
                           group.by = "seurat_clusters",
                           label = TRUE, label.size = 5) +
  ggtitle("UMAP - Seurat clusters (high-quality cells)") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(file.path(fig.path, "UMAP_clusters.pdf"), p_umap_cluster,
       width = 8, height = 7)
ggsave(file.path(fig.path, "UMAP_clusters.png"), p_umap_cluster,
       width = 8, height = 7, dpi = 200)
cat("  OK UMAP 图已保存\n")

pdf(file.path(fig.path, "Marker_Score_Heatmap.pdf"), width = 11, height = 9)
pheatmap(score_matrix,
         scale = "column",
         color = colorRampPalette(c("#F7FCF5", "#74C476", "#00441B"))(100),
         cluster_rows = TRUE, cluster_cols = FALSE,
         main = "Marker-based cell type score",
         fontsize_row = 11, fontsize_col = 11,
         cellwidth = 22, cellheight = 20)
dev.off()
png(file.path(fig.path, "Marker_Score_Heatmap.png"),
    width = 11, height = 9, units = "in", res = 200)
pheatmap(score_matrix,
         scale = "column",
         color = colorRampPalette(c("#F7FCF5", "#74C476", "#00441B"))(100),
         cluster_rows = TRUE, cluster_cols = FALSE,
         main = "Marker-based cell type score",
         fontsize_row = 11, fontsize_col = 11,
         cellwidth = 22, cellheight = 20)
dev.off()
cat("  OK Marker score 热图已保存\n")

cat("\n========== 气泡图 (主管范式: 绿色 + Identity/Features, 精简 marker) ==========\n")

dotplot_markers <- intersect(dotplot_marker_order, rownames(obj_clean))
cat(sprintf("气泡图 marker (精简 panel 中可用): %d / %d\n",
            length(dotplot_markers), length(dotplot_marker_order)))

dot_width <- max(10, length(dotplot_markers) * 0.45 + 3.0)

present_celltypes <- intersect(celltype_display_order, unique(obj_clean$cell_type))
obj_clean$cell_type <- factor(obj_clean$cell_type, levels = present_celltypes)
Idents(obj_clean) <- "cell_type"

p_dot_after <- DotPlot(obj_clean,
                        features = dotplot_markers,
                        cols = c("#F7FCF5", "#00441B"),
                        dot.scale = 8) +
  labs(title = "Canonical markers by annotated cell type (after annotation)",
       x = "Features", y = "Identity") +
  theme_bw(base_size = 14) +
  theme(plot.title    = element_text(hjust = 0.5, size = 15, face = "bold"),
        axis.text.x   = element_text(angle = 45, hjust = 1, vjust = 1,
                                     size = 12, face = "italic"),
        axis.text.y   = element_text(size = 13),
        axis.title    = element_text(size = 13, face = "bold"),
        legend.title  = element_text(size = 11),
        legend.text   = element_text(size = 10),
        panel.grid.minor = element_blank(),
        plot.margin = margin(10, 15, 10, 10))

n_celltypes_clean <- length(unique(obj_clean$cell_type))
dot_height_after <- max(5, n_celltypes_clean * 0.45 + 2.5)

ggsave(file.path(fig.path, "DotPlot_after_annotation.pdf"),
       p_dot_after, width = dot_width, height = dot_height_after,
       limitsize = FALSE)
ggsave(file.path(fig.path, "DotPlot_after_annotation.png"),
       p_dot_after, width = dot_width, height = dot_height_after,
       dpi = 200, bg = "white", limitsize = FALSE)
cat(sprintf("  OK 注释后气泡图已保存 (%.1f x %.1f 英寸)\n",
            dot_width, dot_height_after))

Idents(obj_clean) <- "seurat_clusters"

p_dot_before <- DotPlot(obj_clean,
                         features = dotplot_markers,
                         cols = c("#F7FCF5", "#00441B"),
                         dot.scale = 8) +
  labs(title = "Canonical markers by Seurat cluster (before annotation)",
       x = "Features", y = "Identity") +
  theme_bw(base_size = 14) +
  theme(plot.title    = element_text(hjust = 0.5, size = 15, face = "bold"),
        axis.text.x   = element_text(angle = 45, hjust = 1, vjust = 1,
                                     size = 12, face = "italic"),
        axis.text.y   = element_text(size = 13),
        axis.title    = element_text(size = 13, face = "bold"),
        legend.title  = element_text(size = 11),
        legend.text   = element_text(size = 10),
        panel.grid.minor = element_blank(),
        plot.margin = margin(10, 15, 10, 10))

n_clusters_clean <- length(unique(obj_clean$seurat_clusters))
dot_height_before <- max(5, n_clusters_clean * 0.45 + 2.5)

ggsave(file.path(fig.path, "DotPlot_before_annotation.pdf"),
       p_dot_before, width = dot_width, height = dot_height_before,
       limitsize = FALSE)
ggsave(file.path(fig.path, "DotPlot_before_annotation.png"),
       p_dot_before, width = dot_width, height = dot_height_before,
       dpi = 200, bg = "white", limitsize = FALSE)
cat(sprintf("  OK 注释前气泡图已保存 (%.1f x %.1f 英寸)\n",
            dot_width, dot_height_before))

p_dot_combined <- p_dot_before / p_dot_after +
  plot_annotation(
    title = "Cell annotation: before vs after",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 17, face = "bold"))
  )

combined_h <- dot_height_before + dot_height_after + 1
ggsave(file.path(fig.path, "DotPlot_before_vs_after.pdf"),
       p_dot_combined, width = dot_width, height = combined_h,
       limitsize = FALSE)
ggsave(file.path(fig.path, "DotPlot_before_vs_after.png"),
       p_dot_combined, width = dot_width, height = combined_h,
       dpi = 200, bg = "white", limitsize = FALSE)
cat(sprintf("  OK 注释前后合并对比图已保存 (%.1f x %.1f 英寸)\n",
            dot_width, combined_h))

cat("\n========== Sham vs Stroke 细胞比例 ==========\n")

cell_prop <- obj_clean@meta.data %>%
  group_by(condition, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(pct = n / sum(n) * 100)

p_prop <- ggplot(cell_prop,
                  aes(x = condition, y = pct, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Greens"))(
    length(unique(cell_prop$cell_type)))) +
  labs(title = "Cell type proportions: Sham vs Stroke",
       x = "", y = "Proportion (%)") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.text = element_text(size = 9))
ggsave(file.path(fig.path, "CellProportions_Sham_vs_Stroke.pdf"), p_prop,
       width = 9, height = 7)
ggsave(file.path(fig.path, "CellProportions_Sham_vs_Stroke.png"), p_prop,
       width = 9, height = 7, dpi = 200)
cat("  OK 细胞比例图已保存\n")

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("    第 17 部分 单细胞细胞类型注释 完成\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("Cluster 数量:", length(cluster_levels), "\n")
cat("\n所有 cluster 注释 (含 Low_Quality 标记):\n")
print(table(final_calls$Final_Call))

cat(sprintf("\n>>> 永久剔除 Low_Quality 细胞: %d -> %d (保留 %.1f%%) <<<\n",
            n_before, n_after, 100 * n_after / n_before))

cat("\n高质量细胞类型分布:\n")
print(table(obj_clean$cell_type))

cat("\nConfidence 分布:\n")
print(table(final_calls$Confidence))

cat("\n输出文件:\n")
cat("  Results/marker_source_table.csv   marker + 文献 PMID 溯源表\n")
cat("  Results/final_annotation.csv\n")
cat("  Results/marker_top_calls.csv\n")
cat("  Results/marker_score_matrix.csv\n")
cat("  Results/singleR_*.csv\n")
cat("  seurat_obj_annotated_clean.rds (已剔除 Low_Quality)\n")
cat("\n  Figures/UMAP_*.pdf/png\n")
cat("  Figures/Marker_Score_Heatmap.pdf/png\n")
cat("  Figures/DotPlot_before_annotation.pdf/png  (主管范式: 绿色, 精简 marker)\n")
cat("  Figures/DotPlot_after_annotation.pdf/png   (主管范式: 绿色, 精简 marker)\n")
cat("  Figures/DotPlot_before_vs_after.pdf/png\n")
cat("  Figures/CellProportions_Sham_vs_Stroke.pdf/png\n")

cat("\nOK 第 17 部分分析完成\n")

sink(file.path(res.path, "session_info.txt"))
sessionInfo()
sink()
