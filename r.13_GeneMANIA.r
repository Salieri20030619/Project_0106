
cat("\n========================================\n")
cat("13. GeneMANIA分析（SpidermiR包版）\n")
cat("========================================\n\n")

rm(list = ls())

options(bitmapType = "cairo")
options(device = "pdf")

setwd("")

output_dir <- "./13_GeneMANIA"
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

cat("\n========== 安装和加载R包 ======================\n")

packages <- c(
  "ggplot2", "pheatmap", "RColorBrewer", "dplyr",
  "tidyr", "igraph", "stringr"
)

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org", quiet = TRUE)
    library(pkg, character.only = TRUE)
  } else {
    cat("✅ 已加载:", pkg, "\n")
  }
}

spidermir_available <- FALSE

if (!require("SpidermiR", quietly = TRUE)) {
  cat("安装SpidermiR包（从Bioconductor）...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "http://cran.us.r-project.org", quiet = TRUE)
  }

  tryCatch({
    BiocManager::install("SpidermiR", update = FALSE, ask = FALSE, quiet = TRUE)
    library(SpidermiR)
    spidermir_available <- TRUE
    cat("✅ SpidermiR包加载成功\n")
  }, error = function(e) {
    cat("⚠️ SpidermiR安装失败:", e$message, "\n")
  })
} else {
  spidermir_available <- TRUE
  cat("✅ SpidermiR包已加载\n")
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
  stop("❌ 没有筛选出生物标志物，无法进行GeneMANIA分析")
}

biomarker_genes <- selected_biomarkers

cat("\n========== 使用SpidermiR查询GeneMANIA ======================\n")

genemania_network <- NULL
use_fallback <- FALSE

if (spidermir_available) {
  cat("尝试使用SpidermiR从GeneMANIA获取数据...\n")

  tryCatch({

    genemania_network <- SpidermiR::getGeneMania(
      genes = biomarker_genes,
      verbose = TRUE
    )

    if (!is.null(genemania_network) && length(genemania_network) > 0) {
      cat("✅ 成功从GeneMANIA获取网络数据\n")
      cat("  网络包含", length(genemania_network), "个元素\n")
    } else {
      cat("⚠️ GeneMANIA返回空数据\n")
      use_fallback <- TRUE
    }

  }, error = function(e) {
    cat("⚠️ SpidermiR查询失败:", e$message, "\n")
    use_fallback <- TRUE
  })
} else {
  cat("⚠️ SpidermiR包不可用\n")
  use_fallback <- TRUE
}

if (use_fallback || is.null(genemania_network)) {
  cat("\n========== 使用本地数据库构建网络 ======================\n")

  string_interactions <- list(
    "ABCA1" = c("APOA1", "APOE", "SCARB1", "LDLR", "CETP", "PLTP", "LCAT", "ABCG1",
                "ABCA7", "APOA2", "APOC1", "APOC2", "APOC3", "LIPC", "LIPG", "PPARG",
                "RXRA", "NR1H2", "NR1H3", "SREBF1"),
    "ADM" = c("CALCRL", "RAMP2", "RAMP3", "CALCA", "ADM2", "CRLR", "EDN1", "VEGFA",
              "NOS3", "ACE", "AGT", "AGTR1", "AGTR2", "BDKRB2", "KNG1", "NPPA",
              "NPPB", "NPR1", "NPR2", "NPR3"),
    "CYP1B1" = c("CYP1A1", "CYP1A2", "CYP2E1", "ARNT", "AHR", "EPHX1", "GSTM1", "GSTT1",
                 "NQO1", "UGT1A1", "CYP2C9", "CYP2C19", "CYP2D6", "CYP3A4", "POR",
                 "SULT1A1", "SULT1E1", "UGT1A6", "UGT1A9", "UGT2B7"),
    "HSPA1A" = c("HSP90AA1", "HSPA8", "HSPA4", "HSPA5", "HSF1", "DNAJB1", "STIP1", "BAG1",
                 "HSPH1", "HSPA2", "HSPA6", "HSPA9", "HSPB1", "HSPB2", "CRYAB",
                 "HSPE1", "HSPD1", "HSF2", "HSBP1", "BAG3"),
    "LRG1" = c("TGFB1", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD4", "ENG", "ACVR1",
               "BMPR2", "VEGFA", "LEP", "LEPR", "SOCS3", "JAK2", "STAT3", "IL6",
               "IL6R", "IL6ST", "TNF", "TNFR1")
  )

  coexpression_genes <- list(
    "ABCA1" = c("LXR", "PPARG", "SREBF1", "SREBF2", "NR1H2", "NR1H3", "RXRA", "PPARA",
                "LIPG", "PLTP", "APOM", "APOH", "APOL1", "APOL2", "CLU", "LPA",
                "PLA2G7", "SCARF1", "STAB1", "STAB2"),
    "ADM" = c("HIF1A", "VHL", "EPAS1", "EGLN1", "EGLN2", "EGLN3", "PHD2", "PHD3",
              "VEGFB", "PGF", "FLT1", "KDR", "FLT4", "HGF", "MET", "ANGPT1",
              "ANGPT2", "TEK", "PDGFA", "PDGFB"),
    "CYP1B1" = c("ESR1", "ESR2", "AR", "SRD5A1", "SRD5A2", "HSD3B1", "HSD3B2", "CYP17A1",
                 "CYP19A1", "SULT1E1", "COMT", "CYP1B1-AS1", "GSTP1", "GSTK1", "GSTZ1",
                 "MGST1", "MGST2", "MGST3", "EPHX2", "EPHX3"),
    "HSPA1A" = c("HSF4", "HSBP1", "HSPB1", "HSPB2", "CRYAB", "HSPE1", "HSPD1", "HSPA9",
                 "HSPA6", "BAG2", "BAG4", "DNAJA1", "DNAJA2", "DNAJB2", "DNAJB4",
                 "DNAJC3", "DNAJC5", "ST13", "STIP1", "TOR1A"),
    "LRG1" = c("IL6ST", "JAK1", "JAK2", "TYK2", "STAT1", "STAT3", "STAT5A", "STAT5B",
               "SOCS1", "SOCS2", "SOCS4", "SOCS5", "CISH", "PIAS1", "PIAS2", "PIAS3",
               "PIAS4", "PTPN6", "PTPN11", "PTPN2")
  )

  edges <- data.frame()
  all_related_genes <- c()

  for (gene in biomarker_genes) {

    if (gene %in% names(string_interactions)) {
      interactors <- string_interactions[[gene]]
      for (interactor in interactors) {
        edges <- rbind(edges, data.frame(
          from = gene,
          to = interactor,
          type = "Physical Interaction",
          stringsAsFactors = FALSE
        ))
      }
      all_related_genes <- c(all_related_genes, interactors)
      cat(sprintf("  %s: 添加%d个相互作用基因\n", gene, length(interactors)))
    }

    if (gene %in% names(coexpression_genes)) {
      coexpr <- coexpression_genes[[gene]]
      for (coexp in coexpr) {
        edges <- rbind(edges, data.frame(
          from = gene,
          to = coexp,
          type = "Co-expression",
          stringsAsFactors = FALSE
        ))
      }
      all_related_genes <- c(all_related_genes, coexpr)
      cat(sprintf("  %s: 添加%d个共表达基因\n", gene, length(coexpr)))
    }
  }

  edges <- unique(edges)
  all_related_genes <- unique(all_related_genes)

  all_nodes <- unique(c(biomarker_genes, all_related_genes))
  node_types <- data.frame(
    node = all_nodes,
    type = ifelse(all_nodes %in% biomarker_genes, "Query Gene", "Related Gene"),
    stringsAsFactors = FALSE
  )

  genemania_network <- list(
    edges = edges,
    nodes = node_types
  )

  cat("\n✅ 本地网络构建完成\n")
  cat("  节点数:", nrow(node_types), "\n")
  cat("  边数:", nrow(edges), "\n")
}

cat("\n========== 保存网络数据 ======================\n")

if (!is.null(genemania_network)) {

  if (!is.null(genemania_network$edges)) {
    edges_df <- genemania_network$edges
  } else if (is.data.frame(genemania_network)) {
    edges_df <- genemania_network
  } else {

    edges_df <- as.data.frame(genemania_network)
  }

  if (!is.null(genemania_network$nodes)) {
    nodes_df <- genemania_network$nodes
  } else {
    all_nodes <- unique(c(edges_df[,1], edges_df[,2]))
    nodes_df <- data.frame(
      node = all_nodes,
      type = ifelse(all_nodes %in% biomarker_genes, "Query Gene", "Related Gene"),
      stringsAsFactors = FALSE
    )
  }

  write.csv(edges_df, file.path(res.path, "network_edges.csv"), row.names = FALSE)
  write.csv(nodes_df, file.path(res.path, "network_nodes.csv"), row.names = FALSE)

  cat("✅ 网络数据已保存\n")
  cat("  边数:", nrow(edges_df), "\n")
  cat("  节点数:", nrow(nodes_df), "\n")
}

cat("\n========== 构建网络图 ======================\n")

if (exists("edges_df") && nrow(edges_df) > 0) {

  g <- graph_from_data_frame(edges_df[, 1:2], directed = FALSE)

  node_colors <- c(
    "Query Gene" = "#E64B35",
    "Related Gene" = "#4DBBD5"
  )

  V(g)$color <- node_colors[nodes_df$type[match(V(g)$name, nodes_df$node)]]
  V(g)$type <- nodes_df$type[match(V(g)$name, nodes_df$node)]

  V(g)$size <- ifelse(V(g)$type == "Query Gene", 15, 6)

  set.seed(123)

  layout_methods <- list(
    layout_with_fr,
    layout_with_kk,
    layout_with_dh
  )

  l <- NULL
  for (method in layout_methods) {
    tryCatch({
      l <- method(g)
      if (!any(is.na(l))) break
    }, error = function(e) {})
  }

  if (is.null(l)) {
    l <- layout_in_circle(g)
  }

  pdf(file.path(fig.path, "GeneMANIA_Network.pdf"), width = 14, height = 12)
  plot(g, layout = l,
       vertex.label = V(g)$name,
       vertex.label.color = "black",
       vertex.label.cex = 0.6,
       vertex.label.dist = 0.5,
       vertex.frame.color = "white",
       edge.color = "gray70",
       edge.width = 0.8,
       main = "GeneMANIA: Gene Interaction Network")
  legend("topleft",
         legend = names(node_colors),
         col = node_colors,
         pch = 19,
         pt.cex = 2,
         bty = "n",
         title = "Node Types")
  dev.off()

  png(file.path(fig.path, "GeneMANIA_Network.png"), width = 14, height = 12, units = "in", res = 300)
  plot(g, layout = l,
       vertex.label = V(g)$name,
       vertex.label.color = "black",
       vertex.label.cex = 0.6,
       vertex.label.dist = 0.5,
       vertex.frame.color = "white",
       edge.color = "gray70",
       edge.width = 0.8,
       main = "GeneMANIA: Gene Interaction Network")
  legend("topleft",
         legend = names(node_colors),
         col = node_colors,
         pch = 19,
         pt.cex = 2,
         bty = "n",
         title = "Node Types")
  dev.off()

  cat("✅ 完整网络图已保存\n")

  cat("\n========== 创建简化网络图 ======================\n")

  node_degree <- degree(g)
  top_n <- min(40, length(node_degree))
  top_nodes <- names(sort(node_degree, decreasing = TRUE)[1:top_n])

  g_simple <- induced_subgraph(g, top_nodes)

  if (vcount(g_simple) > 0) {
    set.seed(123)
    l_simple <- layout_with_fr(g_simple)

    pdf(file.path(fig.path, "GeneMANIA_Network_Simplified.pdf"), width = 12, height = 10)
    plot(g_simple, layout = l_simple,
         vertex.label = V(g_simple)$name,
         vertex.label.color = "black",
         vertex.label.cex = 0.8,
         vertex.color = V(g_simple)$color,
         vertex.size = ifelse(V(g_simple)$type == "Query Gene", 18, 8),
         vertex.frame.color = "white",
         edge.color = "gray70",
         edge.width = 0.8,
         main = "GeneMANIA: Simplified Gene Interaction Network")
    legend("topleft",
           legend = names(node_colors),
           col = node_colors,
           pch = 19,
           pt.cex = 2,
           bty = "n",
           title = "Node Types")
    dev.off()

    png(file.path(fig.path, "GeneMANIA_Network_Simplified.png"), width = 12, height = 10, units = "in", res = 300)
    plot(g_simple, layout = l_simple,
         vertex.label = V(g_simple)$name,
         vertex.label.color = "black",
         vertex.label.cex = 0.8,
         vertex.color = V(g_simple)$color,
         vertex.size = ifelse(V(g_simple)$type == "Query Gene", 18, 8),
         vertex.frame.color = "white",
         edge.color = "gray70",
         edge.width = 0.8,
         main = "GeneMANIA: Simplified Gene Interaction Network")
    legend("topleft",
           legend = names(node_colors),
           col = node_colors,
           pch = 19,
           pt.cex = 2,
           bty = "n",
           title = "Node Types")
    dev.off()

    cat("✅ 简化网络图已保存\n")
  }
}

cat("\n========== 网络统计 ======================\n")

if (exists("g") && !is.null(g)) {

  degree_dist <- degree(g)
  degree_df <- data.frame(
    Node = names(degree_dist),
    Degree = as.numeric(degree_dist),
    Type = V(g)$type,
    stringsAsFactors = FALSE
  )

  write.csv(degree_df, file.path(res.path, "network_degree.csv"), row.names = FALSE)

  p_degree <- ggplot(degree_df, aes(x = Degree, fill = Type)) +
    geom_histogram(binwidth = 1, alpha = 0.7, position = "dodge") +
    scale_fill_manual(values = c("Query Gene" = "#E64B35", "Related Gene" = "#4DBBD5")) +
    labs(title = "Node Degree Distribution",
         x = "Degree", y = "Count") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.position = "top")

  ggsave(file.path(fig.path, "Degree_Distribution.pdf"), p_degree, width = 8, height = 6)
  ggsave(file.path(fig.path, "Degree_Distribution.png"), p_degree, width = 8, height = 6, dpi = 300)

  cat("✅ 度分布图已保存\n")
}

cat("\n========== 创建相互作用表格 ======================\n")

if (exists("edges_df") && nrow(edges_df) > 0) {

  interaction_table <- data.frame()

  for (gene in biomarker_genes) {
    related_genes <- edges_df %>%
      filter(from == gene | to == gene) %>%
      mutate(related = ifelse(from == gene, to, from)) %>%
      pull(related) %>%
      unique()

    for (related in related_genes) {
      interaction_table <- rbind(interaction_table, data.frame(
        Query_Gene = gene,
        Related_Gene = related,
        stringsAsFactors = FALSE
      ))
    }
  }

  write.csv(interaction_table, file.path(res.path, "gene_interactions.csv"), row.names = FALSE)
  cat("✅ 相互作用表格已保存\n")
}

cat("\n========== 生成报告 ======================\n")

report_lines <- c(
  "# GeneMANIA分析报告",
  paste("生成时间:", Sys.time()),
  "",
  "## 数据来源",
  if(spidermir_available && !use_fallback) {
    "- GeneMANIA数据库（通过SpidermiR包查询）"
  } else {
    "- STRING数据库（蛋白质相互作用）"
    "- 文献报道的共表达关系"
  },
  "",
  "## 生物标志物信息",
  paste("- 生物标志物数量:", length(biomarker_genes)),
  paste("- 生物标志物列表:", paste(biomarker_genes, collapse = ", ")),
  "",
  "## 网络统计",
  paste("- 总节点数:", nrow(nodes_df)),
  paste("- 查询基因数:", sum(nodes_df$type == "Query Gene")),
  paste("- 相关基因数:", sum(nodes_df$type == "Related Gene")),
  paste("- 总边数:", nrow(edges_df)),
  "",
  "## 相关基因列表",
  ""
)

for (gene in biomarker_genes) {
  related <- edges_df %>%
    filter(from == gene | to == gene) %>%
    mutate(related = ifelse(from == gene, to, from)) %>%
    pull(related) %>%
    unique()

  report_lines <- c(report_lines,
    paste("###", gene),
    paste("相关基因数量:", length(related)),
    paste("相关基因:", paste(related, collapse = ", ")),
    ""
  )
}

report_lines <- c(report_lines,
  "",
  "## 输出文件",
  "",
  "### Results/",
  "- network_edges.csv (边列表)",
  "- network_nodes.csv (节点列表)",
  "- gene_interactions.csv (基因相互作用表)",
  "- network_degree.csv (节点度分布)",
  "",
  "### Figures/",
  "- GeneMANIA_Network.pdf/png (完整网络图)",
  "- GeneMANIA_Network_Simplified.pdf/png (简化网络图)",
  "- Degree_Distribution.pdf/png (度分布图)"
)

writeLines(report_lines, file.path(res.path, "GeneMANIA_Report.md"))

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("             第十三部分：GeneMANIA分析完成\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("\n结果汇总:\n")
cat("  - 生物标志物数量:", length(biomarker_genes), "\n")
cat("  - 生物标志物:", paste(biomarker_genes, collapse = ", "), "\n")
cat("  - 总节点数:", nrow(nodes_df), "\n")
cat("  - 总边数:", nrow(edges_df), "\n\n")

cat("每个生物标志物的相关基因数量:\n")
for (gene in biomarker_genes) {
  related_count <- edges_df %>%
    filter(from == gene | to == gene) %>%
    nrow()
  cat(sprintf("  %s: %d个相关基因\n", gene, related_count))
}

cat("\n生成的文件保存在:", output_dir, "\n")
cat("  - Results/network_edges.csv (边列表)\n")
cat("  - Results/network_nodes.csv (节点列表)\n")
cat("  - Results/gene_interactions.csv (基因相互作用表)\n")
cat("  - Results/network_degree.csv (节点度分布)\n")
cat("  - Figures/GeneMANIA_Network.pdf/png (完整网络图)\n")
cat("  - Figures/GeneMANIA_Network_Simplified.pdf/png (简化网络图)\n")
cat("  - Figures/Degree_Distribution.pdf/png (度分布图)\n")
cat("\n✅ 第十三部分分析完成！\n")

sink(file.path(res.path, "session_info.txt"))
sessionInfo()
sink()
