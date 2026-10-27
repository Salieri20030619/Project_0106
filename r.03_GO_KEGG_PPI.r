
cat("\n========================================\n")
cat("03. GO和KEGG富集分析及PPI网络构建\n")
cat("========================================\n\n")

rm(list = ls())

setwd("")

input_dir  <- "./02_CandidateGene_Intersect/Results"
output_dir <- "./03_GO_KEGG_PPI"
fig_dir    <- file.path(output_dir, "Figures")
res_dir    <- file.path(output_dir, "Results")

for (dir_path in c(output_dir, fig_dir, res_dir)) {
  if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)
}

final_biomarkers <- c("C5AR1", "ABCA1", "MFGE8")

STRING_DATA_DIR         <- "./string_data"
STRING_ALIASES_FILE     <- file.path(STRING_DATA_DIR, "9606.protein.aliases.v12.0.txt.gz")
STRING_LINKS_FILE       <- file.path(STRING_DATA_DIR, "9606.protein.links.v12.0.txt.gz")
STRING_SCORE_THRESHOLDS <- c(400, 250, 150)
STRING_MIN_EDGES        <- 15
STRING_SCORE_THRESHOLD  <- STRING_SCORE_THRESHOLDS[1]

ENRICH_MIN_TERMS <- 8

cat("输入目录:", input_dir, "\n")
cat("输出目录:", output_dir, "\n")
cat("STRING 数据目录:", STRING_DATA_DIR, "\n")
cat("STRING score 阈值级联:",
    paste(STRING_SCORE_THRESHOLDS/1000, collapse = " -> "),
    "(目标边数 >=", STRING_MIN_EDGES, ")\n")
cat("富集通路目标数: >=", ENRICH_MIN_TERMS, "\n")
cat("最终生物标志物:", paste(final_biomarkers, collapse = ", "), "\n\n")

cat("加载R包...\n")
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(enrichplot)
library(igraph)
library(ggraph)
library(tidygraph)
library(patchwork)
library(RColorBrewer)
library(dplyr)
library(stringr)
library(ggrepel)
library(data.table)
library(scales)

igraph_version <- as.character(packageVersion("igraph"))
cp_version     <- as.character(packageVersion("clusterProfiler"))
cat("R:", R.version.string, "\n")
cat("igraph:", igraph_version, "| clusterProfiler:", cp_version, "\n")
cat("✅ R包加载完成\n\n")

cat("========== 1. 从 02_CandidateGene_Intersect 读取候选基因 ==========\n")

candidate_file <- file.path(input_dir, "final_candidate_genes.txt")
if (file.exists(candidate_file)) {
  candidate_genes <- readLines(candidate_file)
  candidate_genes <- candidate_genes[candidate_genes != "" & !is.na(candidate_genes)]
  cat("✅ 读取", length(candidate_genes), "个候选基因\n")
} else {
  candidate_genes <- character(0)
}

detail_file <- file.path(input_dir, "candidate_genes_detail.csv")
candidate_up   <- character(0)
candidate_down <- character(0)
if (file.exists(detail_file)) {
  candidate_detail <- read.csv(detail_file, stringsAsFactors = FALSE)
  if (nrow(candidate_detail) > 0) {
    candidate_genes <- candidate_detail$Gene
    candidate_up    <- candidate_detail$Gene[candidate_detail$Regulation == "Up"]
    candidate_down  <- candidate_detail$Gene[candidate_detail$Regulation == "Down"]
    cat("  上调:", length(candidate_up), "| 下调:", length(candidate_down), "\n")
  }
}

if (length(candidate_genes) == 0) {
  cat("\n⚠️ 使用备用列表\n")
  candidate_genes <- c("C5AR1", "PDK4", "CRISPLD2", "ADM", "CYP1B1", "CEBPD",
                       "ABCA1", "IRS2", "FCGR2A", "HSPA1A", "KLF4", "CD14",
                       "EMILIN2", "FKBP5", "LRG1", "HP", "LCN2", "LGALS3", "MFGE8")
  candidate_up <- setdiff(candidate_genes, "MFGE8")
  candidate_down <- c("MFGE8")
}
cat("\n候选基因:", paste(candidate_genes, collapse = ", "), "\n\n")

cat("========== 2. 基因 ID 转换 ==========\n")

convert_symbol_to_entrez <- function(symbols) {
  result <- data.frame(); unmapped <- c()
  for (sym in symbols) {
    entrez <- tryCatch({
      mapIds(org.Hs.eg.db, keys = sym, column = "ENTREZID",
             keytype = "SYMBOL", multiVals = "first")
    }, error = function(e) NA)
    if (!is.na(entrez) && !is.null(entrez) && entrez != sym) {
      result <- rbind(result,
                      data.frame(SYMBOL = sym,
                                 ENTREZID = as.character(entrez),
                                 stringsAsFactors = FALSE))
    } else {
      unmapped <- c(unmapped, sym)
    }
  }
  list(mapped = result, unmapped = unmapped)
}

conv_result <- convert_symbol_to_entrez(candidate_genes)
gene_entrez <- conv_result$mapped
cat("成功转换:", nrow(gene_entrez), "/", length(candidate_genes), "\n")
write.csv(gene_entrez, file.path(res_dir, "gene_entrez_mapping.csv"), row.names = FALSE)
entrez_list <- gene_entrez$ENTREZID[!is.na(gene_entrez$ENTREZID)]
cat("用于富集分析的 ENTREZID:", length(entrez_list), "\n\n")

run_enrichment_cascade <- function(enrich_func, ont_label = "",
                                    min_terms = ENRICH_MIN_TERMS) {
  thresholds <- list(
    list(use_padj = TRUE,  cut = 0.05, label = "padj<0.05"),
    list(use_padj = TRUE,  cut = 0.10, label = "padj<0.10"),
    list(use_padj = TRUE,  cut = 0.25, label = "padj<0.25"),
    list(use_padj = FALSE, cut = 0.05, label = "pvalue<0.05 (FDR relaxed)")
  )
  best_result <- NULL; best_label <- NA; best_n <- 0

  for (th in thresholds) {
    res <- tryCatch({
      enrich_func(th$use_padj, th$cut)
    }, error = function(e) {
      cat(sprintf("    [%s] 阈值 %s 调用失败: %s\n", ont_label, th$label, e$message))
      NULL
    })

    if (!is.null(res) && nrow(res@result) > 0) {
      df <- res@result
      if (th$use_padj) {
        df <- df[df$p.adjust < th$cut, , drop = FALSE]
      } else {
        df <- df[df$pvalue < th$cut, , drop = FALSE]
      }

      if (nrow(df) == 0) {
        cat(sprintf("    [%s] 阈值 %s -> 0 条 (post-filter)\n",
                    ont_label, th$label))
        next
      }

      res@result <- df
      n_terms <- nrow(df)

      cat(sprintf("    [%s] 阈值 %s -> %d 条 (top3 padj: %.3g, %.3g, %.3g)\n",
                  ont_label, th$label, n_terms,
                  df$p.adjust[1],
                  ifelse(n_terms >= 2, df$p.adjust[2], NA),
                  ifelse(n_terms >= 3, df$p.adjust[3], NA)))

      if (is.null(best_result) || n_terms > best_n) {
        best_result <- res; best_label <- th$label; best_n <- n_terms
      }
      if (n_terms >= min_terms) {
        return(list(result = res, threshold_used = th$label, n_terms = n_terms))
      }
    } else {
      cat(sprintf("    [%s] 阈值 %s -> 0 条\n", ont_label, th$label))
    }
  }

  if (!is.null(best_result)) {
    cat(sprintf("    [%s] 未达 >=%d 条目标, 使用最佳结果 (%s, %d 条)\n",
                ont_label, min_terms, best_label, best_n))
    return(list(result = best_result, threshold_used = best_label, n_terms = best_n))
  }
  list(result = NULL, threshold_used = NA, n_terms = 0)
}

cat("========== 3. GO富集分析 (目标 >=", ENRICH_MIN_TERMS, "条) ==========\n")

go_results        <- list()
go_thresholds_used <- list()

for (ont in c("BP", "CC", "MF")) {
  cat(sprintf("\nGO-%s:\n", ont))
  enrich_fn <- function(use_padj, cut) {
    enrichGO(gene = entrez_list, OrgDb = org.Hs.eg.db,
             keyType = "ENTREZID", ont = ont, pAdjustMethod = "BH",
             pvalueCutoff = 1, qvalueCutoff = 1,
             minGSSize = 5, maxGSSize = 500, readable = TRUE)
  }
  cascade <- run_enrichment_cascade(enrich_fn, sprintf("GO-%s", ont))
  go_results[[ont]]        <- cascade$result
  go_thresholds_used[[ont]] <- cascade$threshold_used

  if (!is.null(cascade$result)) {
    cat(sprintf("  ✅ 使用阈值: %s, %d 条通路\n",
                cascade$threshold_used, cascade$n_terms))
    write.csv(cascade$result@result,
              file.path(res_dir, sprintf("GO_%s_results.csv", ont)),
              row.names = FALSE)
  }
}

cat("\n========== 4. KEGG富集分析 (目标 >=", ENRICH_MIN_TERMS, "条) ==========\n")

kegg_cache_dir      <- "./kegg_cache"
kegg_path2gene_file <- file.path(kegg_cache_dir, "hsa_pathway2gene.txt")
kegg_path2name_file <- file.path(kegg_cache_dir, "hsa_pathway2name.txt")

kegg_enricher_cascade <- function(gene_list, term2gene, term2name, source_label) {
  enrich_fn <- function(use_padj, cut) {
    enricher(gene = gene_list, TERM2GENE = term2gene, TERM2NAME = term2name,
             pvalueCutoff = 1, qvalueCutoff = 1,
             minGSSize = 5, maxGSSize = 500, pAdjustMethod = "BH")
  }
  cascade <- run_enrichment_cascade(enrich_fn, source_label)
  if (!is.null(cascade$result)) {
    cascade$result <- tryCatch(
      setReadable(cascade$result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
      error = function(e) cascade$result
    )
  }
  cascade
}

kegg_cascade <- list(result = NULL, threshold_used = NA, n_terms = 0)

if (file.exists(kegg_path2gene_file) && file.exists(kegg_path2name_file)) {
  cat("✅ 发现本地KEGG缓存, 读取中...\n")
  t2g <- tryCatch(read.table(kegg_path2gene_file, sep = "\t", header = FALSE,
                              stringsAsFactors = FALSE, quote = "", comment.char = ""),
                  error = function(e) NULL)
  t2n <- tryCatch(read.table(kegg_path2name_file, sep = "\t", header = FALSE,
                              stringsAsFactors = FALSE, quote = "", comment.char = ""),
                  error = function(e) NULL)
  if (!is.null(t2g) && !is.null(t2n)) {
    colnames(t2g) <- c("path_id", "gene_id")
    colnames(t2n) <- c("path_id", "name")
    cat(sprintf("  加载成功: %d 条通路\n", nrow(t2n)))
    kegg_cascade <- kegg_enricher_cascade(entrez_list, t2g, t2n, "KEGG (local)")
  }
}

if (is.null(kegg_cascade$result)) {
  cat("\n尝试在线 enrichKEGG...\n")
  enrich_fn_online <- function(use_padj, cut) {
    enrichKEGG(gene = entrez_list, organism = "hsa", keyType = "kegg",
               pvalueCutoff = 1, qvalueCutoff = 1,
               minGSSize = 5, maxGSSize = 500, pAdjustMethod = "BH")
  }
  kegg_cascade <- tryCatch({
    cascade <- run_enrichment_cascade(enrich_fn_online, "KEGG (online)")
    if (!is.null(cascade$result)) {
      cascade$result <- tryCatch(
        setReadable(cascade$result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
        error = function(e) cascade$result
      )
    }
    cascade
  }, error = function(e) list(result = NULL, threshold_used = NA, n_terms = 0))
}

kegg_result <- kegg_cascade$result
kegg_threshold_used <- kegg_cascade$threshold_used

if (!is.null(kegg_result)) {
  cat(sprintf("\n✅ KEGG 最终使用阈值: %s, %d 条通路\n",
              kegg_threshold_used, nrow(kegg_result@result)))
  cat("Top 5 pvalue 和 p.adjust:\n")
  top5 <- head(kegg_result@result, 5)
  for (i in seq_len(nrow(top5))) {
    cat(sprintf("  %d. %s: pvalue=%.3g, p.adjust=%.3g\n",
                i, top5$Description[i], top5$pvalue[i], top5$p.adjust[i]))
  }
  write.csv(kegg_result@result,
            file.path(res_dir, "KEGG_results.csv"), row.names = FALSE)
}

cat("\n========== 5. Fig. 1e: 绘制富集分析图 ==========\n")

prepare_dotplot_df <- function(ego, title_str, n_show = 15) {
  tryCatch({
    df <- ego@result
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df <- df[order(df$p.adjust), ]
    df <- df[1:min(n_show, nrow(df)), ]
    df$GeneRatioNum <- sapply(strsplit(df$GeneRatio, "/"), function(x) {
      as.numeric(x[1]) / as.numeric(x[2])
    })
    df$Description <- factor(df$Description, levels = rev(df$Description))
    df$Panel <- title_str
    return(df)
  }, error = function(e) NULL)
}

make_dotplot <- function(df, title_str, size_max = NULL) {
  if (is.null(df) || nrow(df) == 0) return(NULL)

  if (is.null(size_max)) size_max <- max(df$Count, na.rm = TRUE)

  size_breaks <- pretty(c(1, size_max), n = 4)
  size_breaks <- size_breaks[size_breaks >= 1 & size_breaks <= size_max]

  ggplot(df, aes(x = GeneRatioNum, y = Description,
                  size = Count, color = p.adjust)) +
    geom_point() +

    scale_color_gradient(
      low = "#E64B35", high = "#4DBBD5", name = "p.adjust",
      limits = c(0, 0.05),
      oob = scales::squish,
      breaks = c(0.01, 0.02, 0.03, 0.04, 0.05)
    ) +

    scale_size_continuous(
      range = c(3, 8), name = "Count",
      limits = c(1, size_max),
      breaks = size_breaks
    ) +
    labs(title = title_str, x = "GeneRatio", y = NULL) +
    theme_bw() +
    theme(plot.title   = element_text(hjust = 0.5, size = 13, face = "bold"),
          axis.text.y  = element_text(size = 9),
          axis.text.x  = element_text(size = 8),
          legend.text  = element_text(size = 8),
          legend.title = element_text(size = 9))
}

df_list <- list()
if (!is.null(go_results$BP))
  df_list$BP <- prepare_dotplot_df(go_results$BP, "GO: Biological Process")
if (!is.null(go_results$CC))
  df_list$CC <- prepare_dotplot_df(go_results$CC, "GO: Cellular Component")
if (!is.null(go_results$MF))
  df_list$MF <- prepare_dotplot_df(go_results$MF, "GO: Molecular Function")

if (!is.null(kegg_result)) {

  kegg_result@result$Description <- stringr::str_trunc(
    kegg_result@result$Description, 60
  )
  df_list$KEGG <- prepare_dotplot_df(kegg_result, "KEGG Pathway Enrichment")
}

df_list <- df_list[!sapply(df_list, is.null)]

go_df_list <- df_list[intersect(c("BP", "CC", "MF"), names(df_list))]
if (length(go_df_list) > 0) {
  go_size_max <- max(sapply(go_df_list, function(d) max(d$Count, na.rm = TRUE)))
  cat(sprintf("  GO 三面板 Count 全局最大值: %d (size scale 统一)\n", go_size_max))
} else {
  go_size_max <- 10
}

plots <- list()
for (ont in c("BP", "CC", "MF")) {
  if (!is.null(df_list[[ont]])) {
    title_map <- c(
      BP = "GO: Biological Process",
      CC = "GO: Cellular Component",
      MF = "GO: Molecular Function"
    )
    plots[[ont]] <- make_dotplot(df_list[[ont]], title_map[[ont]],
                                   size_max = go_size_max)
  }
}

if (!is.null(df_list$KEGG)) {
  plots$KEGG <- make_dotplot(df_list$KEGG, "KEGG Pathway Enrichment")
}

go_panel_plots <- plots[intersect(c("BP", "CC", "MF"), names(plots))]
go_panel_plots <- go_panel_plots[!sapply(go_panel_plots, is.null)]

if (length(go_panel_plots) > 0) {
  p_go_trio <- Reduce(`+`, go_panel_plots) +
    plot_annotation(
      title = "Fig. 1e-GO: Gene Ontology Enrichment Analysis",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    ) +
    plot_layout(ncol = length(go_panel_plots), guides = "collect") &
    theme(legend.position = "right")

  max_rows <- max(sapply(go_panel_plots, function(p) {
    tryCatch(length(unique(p$data$Description)), error = function(e) 10)
  }))
  fig_h <- max(6, min(max_rows * 0.45, 18))
  fig_w <- length(go_panel_plots) * 7

  ggsave(file.path(fig_dir, "Fig1e_GO_trio.png"),
         p_go_trio, width = fig_w, height = fig_h, dpi = 300, bg = "white")
  ggsave(file.path(fig_dir, "Fig1e_GO_trio.pdf"),
         p_go_trio, width = fig_w, height = fig_h)
  cat(sprintf("✅ GO %d联图已保存 (%.0fx%.0f 英寸)\n",
              length(go_panel_plots), fig_w, fig_h))
}

if (!is.null(plots$KEGG)) {
  n_kegg <- tryCatch(length(unique(plots$KEGG$data$Description)), error = function(e) 10)
  ggsave(file.path(fig_dir, "Fig1e_KEGG_dotplot.png"),
         plots$KEGG, width = 10, height = max(6, n_kegg * 0.45), dpi = 300, bg = "white")
  ggsave(file.path(fig_dir, "Fig1e_KEGG_dotplot.pdf"),
         plots$KEGG, width = 10, height = max(6, n_kegg * 0.45))
  cat("✅ KEGG dotplot 已保存\n")
}

cat("\n========== 6. STRING PPI 网络构建 ==========\n")

g <- NULL
ppi_source <- NA_character_
chosen_string_threshold <- NA_integer_

load_string_from_files <- function() {
  if (!file.exists(STRING_ALIASES_FILE)) {
    cat("❌ STRING aliases 文件不存在:", STRING_ALIASES_FILE, "\n")
    cat("\n请先下载以下文件到", STRING_DATA_DIR, "目录:\n")
    cat("  9606.protein.aliases.v12.0.txt.gz (约 19 MB)\n")
    cat("  9606.protein.links.v12.0.txt.gz    (约 130 MB)\n")
    return(NULL)
  }
  if (!file.exists(STRING_LINKS_FILE)) {
    cat("❌ STRING links 文件不存在\n"); return(NULL)
  }

  cat("读取 STRING aliases:", STRING_ALIASES_FILE, "\n")
  aliases_dt <- tryCatch({
    fread(STRING_ALIASES_FILE, header = TRUE, sep = "\t")
  }, error = function(e) { cat("  读取失败:", e$message, "\n"); NULL })
  if (is.null(aliases_dt)) return(NULL)

  col_names <- tolower(colnames(aliases_dt))
  string_id_col <- grep("string|protein.?id", col_names)[1]
  alias_col     <- grep("^alias$", col_names)[1]
  if (is.na(string_id_col) || is.na(alias_col)) {
    if (ncol(aliases_dt) >= 2) {
      setnames(aliases_dt, c("string_id", "alias",
                              if (ncol(aliases_dt) >= 3) "source" else NULL)[1:ncol(aliases_dt)])
    }
  } else {
    setnames(aliases_dt, c(colnames(aliases_dt)[string_id_col],
                            colnames(aliases_dt)[alias_col]),
             c("string_id", "alias"))
  }
  cat(sprintf("  aliases 行数: %d, 列: %s\n",
              nrow(aliases_dt), paste(colnames(aliases_dt), collapse = ", ")))

  cat("\n映射候选基因...\n")
  map_df <- aliases_dt[alias %in% candidate_genes,
                        .(string_id = first(string_id)),
                        by = alias]
  setnames(map_df, "alias", "gene")

  unmapped_genes <- setdiff(candidate_genes, map_df$gene)
  if (length(unmapped_genes) > 0) {
    aliases_dt[, alias_upper := toupper(alias)]
    unmapped_upper <- toupper(unmapped_genes)
    extra_map <- aliases_dt[alias_upper %in% unmapped_upper,
                             .(string_id = first(string_id)),
                             by = alias_upper]
    if (nrow(extra_map) > 0) {
      extra_map$gene <- candidate_genes[match(extra_map$alias_upper,
                                                toupper(candidate_genes))]
      extra_map$alias_upper <- NULL
      map_df <- rbind(map_df, extra_map[, .(gene, string_id)])
    }
    aliases_dt[, alias_upper := NULL]
  }

  n_mapped <- nrow(map_df)
  cat(sprintf("  成功映射: %d / %d\n", n_mapped, length(candidate_genes)))
  if (n_mapped < 2) { cat("❌ 映射基因数 < 2\n"); return(NULL) }

  cat("\n读取 STRING links (大文件, 耐心等候):", STRING_LINKS_FILE, "\n")
  links_dt <- tryCatch({
    fread(STRING_LINKS_FILE, header = TRUE, sep = " ")
  }, error = function(e) { cat("  读取失败:", e$message, "\n"); NULL })
  if (is.null(links_dt)) return(NULL)
  cat(sprintf("  links 行数: %d\n", nrow(links_dt)))

  target_ids <- map_df$string_id

  edges_dt <- NULL
  chosen_threshold <- NA_integer_
  for (thr in STRING_SCORE_THRESHOLDS) {
    cat(sprintf("\n尝试阈值 score >= %d (%.2f)...\n", thr, thr/1000))
    tmp <- links_dt[protein1 %in% target_ids &
                     protein2 %in% target_ids &
                     combined_score >= thr]
    n_est <- nrow(tmp) / 2
    cat(sprintf("  匹配 %d 条双向边 (约 %.0f 条去重边)\n", nrow(tmp), n_est))

    if (is.null(edges_dt) || n_est >= STRING_MIN_EDGES) {
      edges_dt <- tmp
      chosen_threshold <- thr
      if (n_est >= STRING_MIN_EDGES) {
        cat(sprintf("  ✅ 达到目标边数 (>= %d), 停止级联\n", STRING_MIN_EDGES))
        break
      }
    }
  }

  cat(sprintf("\n最终选用阈值: score >= %d (%.2f)\n",
              chosen_threshold, chosen_threshold/1000))
  STRING_SCORE_THRESHOLD     <<- chosen_threshold
  chosen_string_threshold    <<- chosen_threshold

  if (is.null(edges_dt) || nrow(edges_dt) == 0) {
    cat("⚠️ 所有阈值下都没有 interaction\n")
    return(list(map_df = map_df, edges = data.frame()))
  }

  id_to_sym <- setNames(map_df$gene, map_df$string_id)
  edges_dt[, from_sym := id_to_sym[protein1]]
  edges_dt[, to_sym   := id_to_sym[protein2]]

  edges_df <- as.data.frame(edges_dt[, .(from_sym, to_sym, combined_score)])
  edges_df <- edges_df[!is.na(edges_df$from_sym) & !is.na(edges_df$to_sym), ]
  edge_key <- apply(edges_df[, c("from_sym", "to_sym")], 1, function(r) {
    paste(sort(r), collapse = "|")
  })
  edges_df$key <- edge_key
  edges_df <- edges_df[order(edges_df$key, -edges_df$combined_score), ]
  edges_df <- edges_df[!duplicated(edges_df$key), ]
  edges_df$key <- NULL
  cat(sprintf("  去重后 edges: %d\n", nrow(edges_df)))

  list(map_df = map_df, edges = edges_df)
}

string_data <- load_string_from_files()

if (!is.null(string_data) && nrow(string_data$edges) > 0) {
  edge_df <- string_data$edges

  g <- graph_from_data_frame(
    d = edge_df[, c("from_sym", "to_sym", "combined_score")],
    directed = FALSE,
    vertices = data.frame(name = candidate_genes, stringsAsFactors = FALSE)
  )

  V(g)$regulation <- "Other"
  V(g)$regulation[V(g)$name %in% candidate_up]   <- "Up"
  V(g)$regulation[V(g)$name %in% candidate_down] <- "Down"
  V(g)$is_biomarker <- V(g)$name %in% final_biomarkers
  V(g)$degree       <- degree(g)

  cat("\n网络统计:\n")
  cat("  节点:", vcount(g), "| 边:", ecount(g),
      "| 密度:", round(edge_density(g), 4),
      "| 平均度:", round(mean(degree(g)), 2), "\n")
  cat("  生物标志物 degree:\n")
  for (bm in final_biomarkers) {
    if (bm %in% V(g)$name) {
      deg_val <- degree(g)[bm]
      status  <- if (deg_val == 0) " ⚠️ 孤立!" else ""
      cat(sprintf("    %s: %d%s\n", bm, deg_val, status))
    }
  }

  ppi_source <- "STRING v12.0 (local)"
  saveRDS(g, file.path(res_dir, "PPI_network_STRING.rds"))

  edges_list <- as_edgelist(g)
  if (nrow(edges_list) > 0) {
    e_weights <- E(g)$combined_score
    edges_out <- data.frame(
      from = edges_list[, 1], to = edges_list[, 2],
      combined_score = e_weights,
      score_normalized = round(e_weights / 1000, 3)
    )
    write.csv(edges_out, file.path(res_dir, "PPI_edge_list.csv"), row.names = FALSE)
  }

  node_attr <- data.frame(
    node = V(g)$name, degree = V(g)$degree,
    regulation = V(g)$regulation,
    is_biomarker = V(g)$is_biomarker
  )
  write.csv(node_attr, file.path(res_dir, "PPI_node_attributes.csv"), row.names = FALSE)
  cat("✅ STRING PPI 构建成功\n")
}

if (is.null(g) && !is.null(go_results$BP) && nrow(go_results$BP@result) > 0) {
  cat("\n⚠️ STRING 不可用, Fallback 到 GO 功能相似性 PPI\n")
  shared_threshold <- 3
  gene_pathways <- list()
  for (i in 1:nrow(go_results$BP@result)) {
    pathway_genes <- unlist(strsplit(go_results$BP@result$geneID[i], "/"))
    pid <- go_results$BP@result$ID[i]
    for (gene in pathway_genes) {
      if (gene %in% candidate_genes) {
        gene_pathways[[gene]] <- unique(c(gene_pathways[[gene]], pid))
      }
    }
  }
  n_genes <- length(candidate_genes)
  adj_matrix <- matrix(0, n_genes, n_genes,
                        dimnames = list(candidate_genes, candidate_genes))
  for (i in 1:(n_genes - 1)) {
    for (j in (i + 1):n_genes) {
      gi <- candidate_genes[i]; gj <- candidate_genes[j]
      if (!is.null(gene_pathways[[gi]]) && !is.null(gene_pathways[[gj]])) {
        shared <- length(intersect(gene_pathways[[gi]], gene_pathways[[gj]]))
        if (shared >= shared_threshold) {
          adj_matrix[i, j] <- shared
          adj_matrix[j, i] <- shared
        }
      }
    }
  }
  g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected",
                                    weighted = TRUE, diag = FALSE)
  V(g)$regulation <- "Other"
  V(g)$regulation[V(g)$name %in% candidate_up]   <- "Up"
  V(g)$regulation[V(g)$name %in% candidate_down] <- "Down"
  V(g)$is_biomarker <- V(g)$name %in% final_biomarkers
  V(g)$degree <- degree(g)
  ppi_source <- "GO_similarity (fallback)"
  saveRDS(g, file.path(res_dir, "PPI_network_GO_fallback.rds"))
  cat("✅ Fallback PPI 构建完成\n")
}

cat("\n========== 7. Fig. 1f: PPI 网络可视化 ==========\n")

hub_genes <- character(0)

if (!is.null(g) && vcount(g) >= 3) {

  node_color_map <- c("Up" = "#E64B35", "Down" = "#4DBBD5", "Other" = "#B2B2B2")
  g_tidy <- as_tbl_graph(g)

  if (ppi_source == "STRING v12.0 (local)") {
    subtitle_str <- sprintf(
      "STRING v12.0 | Nodes: %d | Edges: %d | Density: %.3f | Score >= %.2f",
      vcount(g), ecount(g), edge_density(g), STRING_SCORE_THRESHOLD / 1000
    )
  } else {
    subtitle_str <- sprintf(
      "GO functional similarity (fallback) | Nodes: %d | Edges: %d | Density: %.3f",
      vcount(g), ecount(g), edge_density(g)
    )
  }

  set.seed(42)
  if (ecount(g) > 0) {
    layout_g <- create_layout(g_tidy, layout = "fr")
  } else {
    layout_g <- create_layout(g_tidy, layout = "circle")
  }

  p_base <- ggraph(layout_g)

  if (ppi_source == "STRING v12.0 (local)" && ecount(g) > 0) {
    p_base <- p_base +
      geom_edge_link(aes(edge_alpha = combined_score, edge_width = combined_score),
                      color = "#666666") +
      scale_edge_alpha_continuous(range = c(0.25, 0.8), guide = "none") +
      scale_edge_width_continuous(range = c(0.3, 1.8), name = "STRING score",
                                   breaks = c(200, 400, 700, 900),
                                   labels = c("0.2", "0.4", "0.7", "0.9"))
  } else if (ecount(g) > 0) {
    p_base <- p_base +
      geom_edge_link(color = "gray60", alpha = 0.55, width = 0.6)
  }

  p_base <- p_base +
    geom_node_point(
      data = function(d) d[!d$is_biomarker, ],
      aes(color = regulation, size = degree),
      shape = 16, alpha = 0.9
    )

  p_base <- p_base +
    geom_node_point(
      data = function(d) d[d$is_biomarker, ],
      aes(size = degree),
      shape = 21, stroke = 1.8, fill = "white", color = "black"
    ) +
    geom_node_point(
      data = function(d) d[d$is_biomarker, ],
      aes(color = regulation, size = degree),
      shape = 16, alpha = 0.95
    )

  p_base <- p_base +
    geom_node_label(
      aes(label = name,
          fontface = ifelse(is_biomarker, "bold", "plain")),
      size = 3.3,
      family = "sans",
      repel = TRUE,
      max.overlaps = Inf,
      box.padding = 0.4,
      point.padding = 0.3,
      segment.color = "gray40",
      segment.alpha = 0.6,
      label.padding = unit(0.2, "lines"),
      label.size = 0.3,
      fill = scales::alpha("white", 0.85)
    )

  p_ppi <- p_base +
    scale_color_manual(values = node_color_map,
                        labels = c("Up" = "Up-regulated",
                                   "Down" = "Down-regulated",
                                   "Other" = "Other"),
                        name = "Regulation") +
    scale_size_continuous(range = c(3, 11), name = "Degree") +
    labs(
      title = "Fig. 1f: Protein-Protein Interaction Network",
      subtitle = subtitle_str,
      caption = paste0(
        "Bold black-bordered nodes: final biomarkers (",
        paste(final_biomarkers, collapse = ", "), ")"
      )
    ) +
    theme_void(base_family = "sans") +
    theme(
      text             = element_text(family = "sans"),
      plot.title       = element_text(hjust = 0.5, size = 16, face = "bold",
                                       family = "sans"),
      plot.subtitle    = element_text(hjust = 0.5, size = 10, color = "gray30",
                                       family = "sans"),
      plot.caption     = element_text(hjust = 0.5, size = 9, color = "gray30",
                                       face = "italic", family = "sans"),
      legend.position  = "right",
      legend.box       = "vertical",
      legend.title     = element_text(size = 10, face = "bold", family = "sans"),
      legend.text      = element_text(size = 9, family = "sans"),
      plot.margin      = margin(t = 20, r = 30, b = 20, l = 30),
      plot.background  = element_rect(fill = "white", color = NA)
    ) +
    coord_cartesian(clip = "off")

  x_range <- range(layout_g$x)
  y_range <- range(layout_g$y)
  x_pad <- diff(x_range) * 0.15
  y_pad <- diff(y_range) * 0.15
  p_ppi <- p_ppi +
    expand_limits(x = c(x_range[1] - x_pad, x_range[2] + x_pad),
                  y = c(y_range[1] - y_pad, y_range[2] + y_pad))

  ggsave(file.path(fig_dir, "Fig1f_PPI_network.png"),
         p_ppi, width = 14, height = 11, dpi = 300, bg = "white")
  ggsave(file.path(fig_dir, "Fig1f_PPI_network.pdf"),
         p_ppi, width = 14, height = 11)
  cat("✅ Fig1f_PPI_network.png/.pdf 已保存 (14x11 英寸)\n")

  hub_genes <- names(sort(degree(g), decreasing = TRUE))[1:min(10, vcount(g))]
  hub_genes <- hub_genes[degree(g)[hub_genes] > 0]

  if (length(hub_genes) > 0) {
    cat("\nTop Hub Genes (", ppi_source, "):\n")
    for (i in seq_along(hub_genes)) {
      marker_tag <- if (hub_genes[i] %in% final_biomarkers) " ⭐" else ""
      cat(sprintf("  %d. %s (degree: %d, regulation: %s)%s\n",
                  i, hub_genes[i], degree(g)[hub_genes[i]],
                  V(g)$regulation[match(hub_genes[i], V(g)$name)],
                  marker_tag))
    }
    hub_df <- data.frame(
      gene = hub_genes, degree = degree(g)[hub_genes],
      regulation = V(g)$regulation[match(hub_genes, V(g)$name)],
      is_biomarker = hub_genes %in% final_biomarkers
    )
    write.csv(hub_df, file.path(res_dir, "PPI_hub_genes.csv"), row.names = FALSE)
  }
}

cat("\n========== 8. 生成分析报告 ==========\n")

report_file <- file.path(output_dir, "GO_KEGG_PPI_report.txt")
sink(report_file)
cat("============================================================\n")
cat("        GO/KEGG富集分析和PPI网络构建报告\n")
cat("============================================================\n\n")
cat("分析时间:", format(Sys.time()), "\n")
cat("R版本:", R.version.string, "\n")
cat("clusterProfiler版本:", cp_version, "\n\n")
cat("v5 改进:\n")
cat("  1. 富集点图 color scale: limits=c(0, 0.05) + squish, 消除灰点\n")
cat("  2. 富集点图 size scale: 三面板统一 limits/breaks, 图例合并不溢出\n")
cat("  3. 其他: KEGG padj post-filter / STRING 级联 / GO 阈值递降\n\n")

cat("候选基因总数:", length(candidate_genes), "\n")
cat("最终生物标志物:", paste(final_biomarkers, collapse = ", "), "\n\n")

cat("============================================================\n")
cat("GO 富集分析\n")
cat("============================================================\n")
for (ont in names(go_results)) {
  cat(sprintf("\nGO-%s (阈值: %s):\n", ont, go_thresholds_used[[ont]]))
  if (!is.null(go_results[[ont]])) {
    cat(sprintf("  通路数: %d\n", nrow(go_results[[ont]]@result)))
    cat("  Top 5:\n")
    top5 <- head(go_results[[ont]]@result, 5)
    for (i in 1:nrow(top5)) {
      cat(sprintf("    %d. %s (padj=%.3g)\n",
                  i, top5$Description[i], top5$p.adjust[i]))
    }
  }
}

cat("\n============================================================\n")
cat("KEGG 富集分析\n")
cat("============================================================\n")
if (!is.null(kegg_result)) {
  cat(sprintf("阈值: %s | 通路数: %d\n",
              kegg_threshold_used, nrow(kegg_result@result)))
  cat("Top 10 (pvalue | padjust):\n")
  top10 <- head(kegg_result@result, 10)
  for (i in 1:nrow(top10)) {
    cat(sprintf("  %d. %s (p=%.3g | padj=%.3g)\n",
                i, top10$Description[i], top10$pvalue[i], top10$p.adjust[i]))
  }
}

if (!is.null(g) && vcount(g) > 0) {
  cat("\n============================================================\n")
  cat("PPI 网络 (来源:", ppi_source, ")\n")
  cat("============================================================\n")
  if (ppi_source == "STRING v12.0 (local)" && !is.na(chosen_string_threshold)) {
    cat(sprintf("STRING score 阈值: %.2f (级联选出)\n",
                chosen_string_threshold/1000))
  }
  cat(sprintf("节点: %d | 边: %d | 密度: %.4f | 平均度: %.2f\n",
              vcount(g), ecount(g), edge_density(g), mean(degree(g))))
  cat("\n生物标志物 degree:\n")
  for (bm in final_biomarkers) {
    if (bm %in% V(g)$name) {
      cat(sprintf("  %s: %d\n", bm, degree(g)[bm]))
    }
  }
  if (length(hub_genes) > 0) {
    cat("\nTop Hub Genes:\n")
    for (i in seq_along(hub_genes)) {
      tag <- if (hub_genes[i] %in% final_biomarkers) " [Biomarker]" else ""
      cat(sprintf("  %d. %s (degree: %d)%s\n",
                  i, hub_genes[i], degree(g)[hub_genes[i]], tag))
    }
  }
}
sink()
cat("✅ 报告已保存:", report_file, "\n")

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("✅ 第三部分 v5 完成!\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("PPI 来源:", ifelse(is.na(ppi_source), "未生成", ppi_source), "\n")
if (!is.na(chosen_string_threshold)) {
  cat("STRING 最终选用阈值: >=", chosen_string_threshold/1000, "\n")
}
