
cat("\n========================================\n")
cat("16. 单细胞测序数据处理\n")
cat("========================================\n\n")

rm(list = ls())

options(bitmapType = "cairo")
options(device = "pdf")

setwd("")

output_dir <- "./16_scRNAseq"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  cat("创建输出文件夹:", output_dir, "\n")
}

res.path <- file.path(output_dir, "Results")
fig.path <- file.path(output_dir, "Figures")
data.path <- file.path(output_dir, "Data")

for (path in c(res.path, fig.path, data.path)) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    cat("创建文件夹:", path, "\n")
  }
}

cat("\n========== 加载R包 ==========\n")

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")
  library(BiocManager)
}

seurat_packages <- c("Seurat", "SeuratObject")
for (pkg in seurat_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
    library(pkg, character.only = TRUE)
  } else {
    cat("✅ 已加载:", pkg, "\n")
  }
}

packages <- c(
  "ggplot2", "pheatmap", "RColorBrewer", "dplyr",
  "tidyr", "stringr", "Matrix", "patchwork"
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
if (file.exists(biomarkers_file)) {
  selected_biomarkers <- readLines(biomarkers_file)
  selected_biomarkers <- selected_biomarkers[selected_biomarkers != "" & !is.na(selected_biomarkers)]
  cat("从第七部分读取到的生物标志物:", paste(selected_biomarkers, collapse = ", "), "\n")
} else {
  selected_biomarkers <- c("ABCA1", "LRG1", "CYP1B1", "HSPA1A", "ADM")
  cat("使用默认生物标志物:", paste(selected_biomarkers, collapse = ", "), "\n")
}

cat("\n========== 读取单细胞数据 ======================\n")

all_files <- list.files(output_dir, pattern = "*_counts.csv.gz", full.names = TRUE)

brain_counts_files <- all_files[!grepl("Blood", all_files)]

if (length(brain_counts_files) == 0) {

  all_files_current <- list.files(".", pattern = "*_counts.csv.gz", full.names = TRUE)
  brain_counts_files <- all_files_current[!grepl("Blood", all_files_current)]
  cat("在当前目录找到", length(brain_counts_files), "个脑组织counts文件\n")
}

if (length(brain_counts_files) == 0) {
  stop("❌ 未找到任何脑组织counts文件，请检查文件路径")
}

cat("找到", length(brain_counts_files), "个脑组织counts文件\n")

all_counts <- list()
sample_info <- data.frame()

for (file in brain_counts_files) {

  basename <- gsub(".csv.gz", "", basename(file))
  sample_id <- gsub("_Brain_.*", "", basename)

  cat("读取样本:", sample_id, "\n")

  tryCatch({

    counts <- read.csv(file, row.names = 1, check.names = FALSE)

    colnames(counts) <- paste0(sample_id, "_", colnames(counts))

    all_counts[[sample_id]] <- counts

    temp_info <- data.frame(
      cell = colnames(counts),
      sample = sample_id,
      stringsAsFactors = FALSE
    )
    sample_info <- rbind(sample_info, temp_info)

    cat("  ", ncol(counts), "个细胞,", nrow(counts), "个基因\n")

  }, error = function(e) {
    cat("  ⚠️ 读取失败:", e$message, "\n")
  })
}

cat("\n成功读取", length(all_counts), "个样本\n")

cat("\n========== 合并数据 ======================\n")

if (length(all_counts) > 0) {

  all_genes <- unique(unlist(lapply(all_counts, rownames)))
  cat("总基因数:", length(all_genes), "\n")

  combined_counts <- Matrix(0, nrow = length(all_genes), ncol = sum(sapply(all_counts, ncol)))
  rownames(combined_counts) <- all_genes

  col_idx <- 1
  for (sample_id in names(all_counts)) {
    counts <- all_counts[[sample_id]]
    n_cells <- ncol(counts)

    common_genes <- intersect(rownames(counts), all_genes)
    combined_counts[common_genes, col_idx:(col_idx + n_cells - 1)] <- as.matrix(counts[common_genes, ])

    col_idx <- col_idx + n_cells
  }

  colnames(combined_counts) <- sample_info$cell

  cat("合并矩阵维度:", dim(combined_counts), "\n")

  seurat_obj <- CreateSeuratObject(
    counts = combined_counts,
    project = "GSE225948_Brain",
    min.cells = 3,
    min.features = 200
  )

  seurat_obj$sample <- sample_info$sample[match(colnames(seurat_obj), sample_info$cell)]

  cat("Seurat对象创建成功\n")
  cat("  细胞数:", ncol(seurat_obj), "\n")
  cat("  基因数:", nrow(seurat_obj), "\n")
} else {
  stop("❌ 没有成功读取任何样本")
}

cat("\n========== 质量控制 ======================\n")

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")

cat("\n质控前:\n")
cat("  细胞数:", ncol(seurat_obj), "\n")
cat("  基因数(中位数):", median(seurat_obj$nFeature_RNA), "\n")
cat("  UMI数(中位数):", median(seurat_obj$nCount_RNA), "\n")
cat("  线粒体比例(中位数):", median(seurat_obj$percent.mt), "%\n")

p1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  ggtitle("QC: nCount vs % Mitochondrial")
p2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  ggtitle("QC: nCount vs nFeature")
p_qc <- p1 + p2

ggsave(file.path(fig.path, "Fig14a_QC_Scatter.pdf"), p_qc, width = 12, height = 5)
ggsave(file.path(fig.path, "Fig14a_QC_Scatter.png"), p_qc, width = 12, height = 5, dpi = 300)

seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA > 200 &
           nFeature_RNA < 4000 &
           percent.mt < 15
)

cat("\n质控后:\n")
cat("  细胞数:", ncol(seurat_obj), "\n")
cat("  基因数(中位数):", median(seurat_obj$nFeature_RNA), "\n")
cat("  UMI数(中位数):", median(seurat_obj$nCount_RNA), "\n")
cat("  线粒体比例(中位数):", median(seurat_obj$percent.mt), "%\n")

saveRDS(seurat_obj, file.path(data.path, "seurat_obj_qced.rds"))

cat("\n========== 标准化和归一化 ======================\n")

seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

cat("\n========== 识别高变基因 ======================\n")

seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)

top10 <- head(VariableFeatures(seurat_obj), 10)

p_var <- VariableFeaturePlot(seurat_obj)
p_var_label <- LabelPoints(plot = p_var, points = top10, repel = TRUE)

ggsave(file.path(fig.path, "Fig14b_Variable_Genes.pdf"), p_var_label, width = 10, height = 6)
ggsave(file.path(fig.path, "Fig14b_Variable_Genes.png"), p_var_label, width = 10, height = 6, dpi = 300)

cat("高变基因数:", length(VariableFeatures(seurat_obj)), "\n")

cat("\n========== 数据缩放 ======================\n")

all_genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all_genes)

cat("\n========== PCA降维 ======================\n")

seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), npcs = 50)

p_pca <- DimPlot(seurat_obj, reduction = "pca", group.by = "sample") +
  ggtitle("PCA: Cells by Sample")

ggsave(file.path(fig.path, "Fig14c_PCA_Plot.pdf"), p_pca, width = 10, height = 8)
ggsave(file.path(fig.path, "Fig14c_PCA_Plot.png"), p_pca, width = 10, height = 8, dpi = 300)

cat("\n========== JackStraw方法确定PC数量 ======================\n")

cat("运行JackStraw分析（这可能需要几分钟）...\n")

tryCatch({
  seurat_obj <- JackStraw(seurat_obj, num.replicate = 100, dims = 50)
  seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:50)

  p_jackstraw <- JackStrawPlot(seurat_obj, dims = 1:30)
  ggsave(file.path(fig.path, "Fig14c_JackStraw_Plot.pdf"), p_jackstraw, width = 12, height = 8)
  ggsave(file.path(fig.path, "Fig14c_JackStraw_Plot.png"), p_jackstraw, width = 12, height = 8, dpi = 300)

}, error = function(e) {
  cat("JackStraw分析失败:", e$message, "\n")
})

p_elbow <- ElbowPlot(seurat_obj, ndims = 50)
ggsave(file.path(fig.path, "Fig14c_Elbow_Plot.pdf"), p_elbow, width = 8, height = 6)
ggsave(file.path(fig.path, "Fig14c_Elbow_Plot.png"), p_elbow, width = 8, height = 6, dpi = 300)

elbow_data <- p_elbow$data
stdev <- elbow_data$stdev

diffs <- diff(stdev)
diffs2 <- diff(diffs)
elbow_point <- which.min(diffs2) + 2

if (is.na(elbow_point) || elbow_point < 2) elbow_point <- 10
if (elbow_point > 30) elbow_point <- 20

cat("\n碎石图拐点PC数:", elbow_point, "\n")

cat("\n========== 聚类分析 ======================\n")

pcs_use <- min(elbow_point, 30)

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:pcs_use)

seurat_obj <- FindClusters(seurat_obj, resolution = 0.8)

seurat_obj <- RunUMAP(seurat_obj, dims = 1:pcs_use)

p_umap_clusters <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle(paste0("UMAP: ", length(unique(seurat_obj$seurat_clusters)), " Clusters"))

p_umap_sample <- DimPlot(seurat_obj, reduction = "umap", group.by = "sample", pt.size = 0.5) +
  ggtitle("UMAP: Colored by Sample")

p_umap_combined <- p_umap_clusters + p_umap_sample

ggsave(file.path(fig.path, "Fig14d_UMAP_Clusters.pdf"), p_umap_combined, width = 14, height = 6)
ggsave(file.path(fig.path, "Fig14d_UMAP_Clusters.png"), p_umap_combined, width = 14, height = 6, dpi = 300)

cat("聚类完成，聚类数:", length(unique(seurat_obj$seurat_clusters)), "\n")

cat("\n========== 保存结果 ======================\n")

saveRDS(seurat_obj, file.path(data.path, "seurat_obj_clustered.rds"))

cluster_info <- data.frame(
  cell = colnames(seurat_obj),
  cluster = seurat_obj$seurat_clusters,
  sample = seurat_obj$sample,
  nCount_RNA = seurat_obj$nCount_RNA,
  nFeature_RNA = seurat_obj$nFeature_RNA,
  percent.mt = seurat_obj$percent.mt,
  stringsAsFactors = FALSE
)
write.csv(cluster_info, file.path(res.path, "cluster_info.csv"), row.names = FALSE)

pc_info <- data.frame(
  Method = c("Elbow_Plot", "PCs_Used"),
  Value = c(elbow_point, pcs_use),
  stringsAsFactors = FALSE
)
write.csv(pc_info, file.path(res.path, "pc_determination.csv"), row.names = FALSE)

cat("\n========== 生成报告 ======================\n")

report_lines <- c(
  "# 单细胞测序数据处理报告",
  paste("生成时间:", Sys.time()),
  "",
  "## 数据信息",
  paste("- 数据集: GSE225948"),
  paste("- 成功读取样本数:", length(all_counts)),
  paste("- 原始细胞数:", nrow(sample_info)),
  "",
  "## 质量控制",
  paste("- 质控前细胞数:", nrow(sample_info)),
  paste("- 质控后细胞数:", ncol(seurat_obj)),
  paste("- 质控后基因数:", nrow(seurat_obj)),
  paste("- 高变基因数:", length(VariableFeatures(seurat_obj))),
  "",
  "## PCA主成分确定",
  paste("- 碎石图拐点PC数:", elbow_point),
  paste("- 实际使用PC数:", pcs_use),
  "",
  "## 降维和聚类",
  paste("- 聚类分辨率: 0.8"),
  paste("- 最终聚类数:", length(unique(seurat_obj$seurat_clusters))),
  "",
  "## 输出文件",
  "",
  "### Results/",
  "- cluster_info.csv (细胞聚类信息)",
  "- pc_determination.csv (PC确定信息)",
  "",
  "### Data/",
  "- seurat_obj_qced.rds (质控后Seurat对象)",
  "- seurat_obj_clustered.rds (聚类后Seurat对象，供第17部分使用)",
  "",
  "### Figures/",
  "- Fig14a_QC_Scatter.pdf/png (质量控制散点图)",
  "- Fig14b_Variable_Genes.pdf/png (高变基因)",
  "- Fig14c_PCA_Plot.pdf/png (PCA图)",
  "- Fig14c_JackStraw_Plot.pdf/png (JackStraw图)",
  "- Fig14c_Elbow_Plot.pdf/png (碎石图)",
  "- Fig14d_UMAP_Clusters.pdf/png (UMAP聚类图)"
)

writeLines(report_lines, file.path(res.path, "scRNAseq_Report.md"))

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("             第十六部分：单细胞测序数据处理完成\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("\n结果汇总:\n")
cat("  - 原始细胞数:", nrow(sample_info), "\n")
cat("  - 质控后细胞数:", ncol(seurat_obj), "\n")
cat("  - 聚类数:", length(unique(seurat_obj$seurat_clusters)), "\n")
cat("  - 使用PC数:", pcs_use, "\n\n")

cat("生成的文件保存在:", output_dir, "\n")
cat("  - Data/seurat_obj_qced.rds (质控后Seurat对象)\n")
cat("  - Data/seurat_obj_clustered.rds (聚类后Seurat对象)\n")
cat("  - Results/cluster_info.csv (细胞聚类信息)\n")
cat("  - Results/pc_determination.csv (PC确定信息)\n")
cat("  - Figures/Fig14a_QC_Scatter.pdf/png (质控图)\n")
cat("  - Figures/Fig14b_Variable_Genes.pdf/png (高变基因)\n")
cat("  - Figures/Fig14c_PCA_Plot.pdf/png (PCA图)\n")
cat("  - Figures/Fig14c_JackStraw_Plot.pdf/png (JackStraw图)\n")
cat("  - Figures/Fig14c_Elbow_Plot.pdf/png (碎石图)\n")
cat("  - Figures/Fig14d_UMAP_Clusters.pdf/png (UMAP聚类图)\n")
cat("\n✅ 第十六部分分析完成！\n")

sink(file.path(res.path, "session_info.txt"))
sessionInfo()
sink()
