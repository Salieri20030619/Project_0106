
cat("\n========================================\n")
cat("11. lncRNA-miRNA-mRNA 网络构建\n")
cat("========================================\n\n")

rm(list = ls())

options(bitmapType = "cairo")
options(device = "pdf")

setwd("")

output_dir <- "./11_LncRNA_miRNA_mRNA_Network"
res.path <- file.path(output_dir, "Results")
fig.path <- file.path(output_dir, "Figures")

for (path in c(output_dir, res.path, fig.path)) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

cat("\n========== 加载 R 包 ==========\n")

packages <- c(
  "tidyverse", "ggplot2", "igraph", "RColorBrewer",
  "ggraph", "ggrepel", "tidygraph", "grid"
)

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org", quiet = TRUE)
    library(pkg, character.only = TRUE)
  }
  cat("✅", pkg, "\n")
}

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "http://cran.us.r-project.org", quiet = TRUE)
  library(BiocManager)
}

multiMir_available <- FALSE
tryCatch({
  if (!require("multiMiR", quietly = TRUE)) {
    BiocManager::install("multiMiR", update = FALSE, ask = FALSE, quiet = TRUE)
    library(multiMiR)
  }
  multiMir_available <- TRUE
  cat("✅ multiMiR\n")
}, error = function(e) {
  cat("⚠️ multiMiR 加载失败:", e$message, "\n")
})

cat("\n========== 读取生物标志物 ==========\n")

biomarkers_file <- "./07_Expression_Validation/Results/selected_biomarkers.txt"
if (!file.exists(biomarkers_file)) {
  stop("❌ selected_biomarkers.txt 不存在")
}

biomarker_genes <- readLines(biomarkers_file)
biomarker_genes <- biomarker_genes[nzchar(biomarker_genes) & !is.na(biomarker_genes)]
cat("生物标志物:", paste(biomarker_genes, collapse = ", "), "\n")

builtin_mirna_database <- list(
  "ABCA1"  = c("hsa-miR-33a-5p", "hsa-miR-144-3p", "hsa-miR-26b-5p", "hsa-miR-758-3p", "hsa-miR-27a-3p"),
  "ADM"    = c("hsa-miR-200a-3p", "hsa-miR-146a-5p", "hsa-miR-34a-5p", "hsa-miR-155-5p", "hsa-miR-200b-3p"),
  "CYP1B1" = c("hsa-miR-27b-3p", "hsa-miR-27a-3p", "hsa-miR-19a-3p", "hsa-miR-133b", "hsa-miR-449a", "hsa-miR-17-5p", "hsa-miR-93-5p", "hsa-miR-335-5p"),
  "HSPA1A" = c("hsa-miR-320a", "hsa-miR-16-5p", "hsa-miR-34a-5p", "hsa-miR-21-5p", "hsa-miR-320b"),
  "LRG1"   = c("hsa-miR-146a-5p", "hsa-miR-155-5p", "hsa-miR-223-3p", "hsa-miR-21-5p", "hsa-miR-29b-3p")
)

builtin_lncrna_database <- list(
  "MALAT1"   = c("hsa-miR-146a-5p", "hsa-miR-155-5p", "hsa-miR-21-5p", "hsa-miR-223-3p",
                 "hsa-miR-33a-5p", "hsa-miR-34a-5p", "hsa-miR-16-5p", "hsa-miR-27b-3p",
                 "hsa-miR-200a-3p", "hsa-miR-320a", "hsa-miR-29b-3p", "hsa-miR-26b-5p",
                 "hsa-miR-144-3p", "hsa-miR-17-5p", "hsa-miR-93-5p", "hsa-miR-335-5p", "hsa-miR-423-5p"),
  "NEAT1"    = c("hsa-miR-146a-5p", "hsa-miR-155-5p", "hsa-miR-21-5p", "hsa-miR-223-3p",
                 "hsa-miR-33a-5p", "hsa-miR-34a-5p", "hsa-miR-16-5p", "hsa-miR-27a-3p",
                 "hsa-miR-200b-3p", "hsa-miR-320a", "hsa-miR-449a"),
  "XIST"     = c("hsa-miR-155-5p", "hsa-miR-21-5p", "hsa-miR-223-3p", "hsa-miR-146a-5p",
                 "hsa-miR-34a-5p", "hsa-miR-29b-3p"),
  "HOTAIR"   = c("hsa-miR-21-5p", "hsa-miR-146a-5p", "hsa-miR-223-3p", "hsa-miR-34a-5p",
                 "hsa-miR-27b-3p", "hsa-miR-200a-3p", "hsa-miR-19a-3p"),
  "GAS5"     = c("hsa-miR-21-5p", "hsa-miR-223-3p", "hsa-miR-146a-5p", "hsa-miR-16-5p",
                 "hsa-miR-34a-5p", "hsa-miR-155-5p"),
  "MEG3"     = c("hsa-miR-21-5p", "hsa-miR-223-3p", "hsa-miR-146a-5p", "hsa-miR-34a-5p",
                 "hsa-miR-29b-3p", "hsa-miR-27a-3p", "hsa-miR-133b"),
  "H19"      = c("hsa-miR-146a-5p", "hsa-miR-155-5p", "hsa-miR-21-5p", "hsa-miR-200a-3p",
                 "hsa-miR-200b-3p", "hsa-miR-758-3p"),
  "TUG1"     = c("hsa-miR-21-5p", "hsa-miR-146a-5p", "hsa-miR-223-3p", "hsa-miR-34a-5p",
                 "hsa-miR-27a-3p", "hsa-miR-320a", "hsa-miR-144-3p"),
  "SNHG1"    = c("hsa-miR-21-5p", "hsa-miR-146a-5p", "hsa-miR-223-3p", "hsa-miR-33a-5p",
                 "hsa-miR-16-5p", "hsa-miR-155-5p"),
  "PVT1"     = c("hsa-miR-146a-5p", "hsa-miR-155-5p", "hsa-miR-21-5p", "hsa-miR-223-3p",
                 "hsa-miR-33a-5p", "hsa-miR-27b-3p", "hsa-miR-200b-3p")
)

cat("\n========== Step 1: 查询每个 mRNA 的 miRNA (miRTarBase 优先, 主方案补足) ==========\n")

TOP_N_PER_GENE <- 5
mirna_results  <- list()

query_mirtarbase <- function(gene) {
  res <- tryCatch({
    get_multimir(org = "hsa", target = gene,
                 table = "mirtarbase", summary = TRUE)
  }, error = function(e) {
    cat("    miRTarBase 查询失败:", e$message, "\n"); NULL
  })
  if (is.null(res) || nrow(res@data) == 0) return(NULL)

  df <- res@data
  df <- df[!is.na(df$mature_mirna_id) & nzchar(df$mature_mirna_id), , drop = FALSE]
  if (nrow(df) == 0) return(NULL)

  agg <- do.call(rbind, lapply(unique(df$mature_mirna_id), function(m) {
    sub <- df[df$mature_mirna_id == m, , drop = FALSE]
    data.frame(
      mirna     = m,
      strong_n  = sum(!grepl("Weak", sub$support_type)),
      pmid_n    = length(unique(sub$pubmed_id)),
      record_n  = nrow(sub),
      stringsAsFactors = FALSE
    )
  }))
  agg <- agg[order(-agg$strong_n, -agg$pmid_n, -agg$record_n), ]
  agg
}

query_predicted_union <- function(gene) {
  ts_mirnas <- character(0)
  res_ts <- tryCatch({
    get_multimir(org = "hsa", target = gene, table = "targetscan", summary = TRUE)
  }, error = function(e) NULL)
  if (!is.null(res_ts) && nrow(res_ts@data) > 0) {
    ts_mirnas <- unique(res_ts@data$mature_mirna_id[
      !is.na(res_ts@data$mature_mirna_id) & nzchar(res_ts@data$mature_mirna_id)])
  }

  mirdb_df <- NULL
  res_mirdb <- tryCatch({
    get_multimir(org = "hsa", target = gene, table = "mirdb", summary = TRUE)
  }, error = function(e) NULL)
  if (!is.null(res_mirdb) && nrow(res_mirdb@data) > 0) {
    mirdb_df <- res_mirdb@data
    mirdb_df <- mirdb_df[!is.na(mirdb_df$mature_mirna_id) & nzchar(mirdb_df$mature_mirna_id), ]
  }

  if (length(ts_mirnas) == 0 && (is.null(mirdb_df) || nrow(mirdb_df) == 0)) {
    return(character(0))
  }

  union_set <- union(ts_mirnas,
                     if (!is.null(mirdb_df)) unique(mirdb_df$mature_mirna_id) else character(0))

  if (!is.null(mirdb_df) && nrow(mirdb_df) > 0) {
    mirdb_df$score_num <- suppressWarnings(as.numeric(mirdb_df$score))
    mirdb_df$score_num[is.na(mirdb_df$score_num)] <- 0
    score_by_mirna <- aggregate(score_num ~ mature_mirna_id, data = mirdb_df, FUN = max)
    rownames(score_by_mirna) <- score_by_mirna$mature_mirna_id
    scores <- setNames(rep(0, length(union_set)), union_set)
    in_mirdb <- intersect(union_set, score_by_mirna$mature_mirna_id)
    scores[in_mirdb] <- score_by_mirna[in_mirdb, "score_num"]
    union_set <- names(sort(scores, decreasing = TRUE))
  }
  union_set
}

for (gene in biomarker_genes) {
  cat("\n----", gene, "----\n")
  selected_mirnas <- character(0)
  source_tags     <- character(0)

  if (multiMir_available) {
    cat("  [1] 查询 miRTarBase (实验验证)...\n")
    mtb_agg <- query_mirtarbase(gene)
    if (!is.null(mtb_agg) && nrow(mtb_agg) > 0) {
      n_take <- min(TOP_N_PER_GENE, nrow(mtb_agg))
      take <- mtb_agg$mirna[1:n_take]
      selected_mirnas <- c(selected_mirnas, take)
      source_tags     <- c(source_tags, rep("miRTarBase", n_take))
      cat(sprintf("      miRTarBase: %d 个唯一 miRNA, 取 Top %d\n",
                  nrow(mtb_agg), n_take))
      cat("      选中:", paste(take, collapse = ", "), "\n")
    } else {
      cat("      miRTarBase: 无数据\n")
    }

    if (length(selected_mirnas) < TOP_N_PER_GENE) {
      need <- TOP_N_PER_GENE - length(selected_mirnas)
      cat(sprintf("  [2] 不足 %d, 从 TargetScan ∪ miRDB 补 %d 个...\n",
                  TOP_N_PER_GENE, need))
      pred <- query_predicted_union(gene)
      pred <- setdiff(pred, selected_mirnas)
      if (length(pred) > 0) {
        take <- pred[1:min(need, length(pred))]
        selected_mirnas <- c(selected_mirnas, take)
        source_tags     <- c(source_tags, rep("TargetScan∪miRDB", length(take)))
        cat(sprintf("      并集补 %d 个: %s\n",
                    length(take), paste(take, collapse = ", ")))
      } else {
        cat("      并集无可补 miRNA\n")
      }
    }

    if (length(selected_mirnas) < TOP_N_PER_GENE && gene %in% names(builtin_mirna_database)) {
      need <- TOP_N_PER_GENE - length(selected_mirnas)
      cat(sprintf("  [3] 仍不足, 从文献库补 %d 个...\n", need))
      lit <- setdiff(builtin_mirna_database[[gene]], selected_mirnas)
      if (length(lit) > 0) {
        take <- lit[1:min(need, length(lit))]
        selected_mirnas <- c(selected_mirnas, take)
        source_tags     <- c(source_tags, rep("Literature", length(take)))
        cat(sprintf("      文献库补 %d 个: %s\n",
                    length(take), paste(take, collapse = ", ")))
      }
    }
  } else {
    if (gene %in% names(builtin_mirna_database)) {
      selected_mirnas <- builtin_mirna_database[[gene]][1:min(TOP_N_PER_GENE,
                                                              length(builtin_mirna_database[[gene]]))]
      source_tags     <- rep("Literature (multiMiR unavailable)", length(selected_mirnas))
    }
  }

  src_summary <- paste(sprintf("%s(%d)",
                               names(table(source_tags)),
                               as.integer(table(source_tags))),
                       collapse = " + ")

  mirna_results[[gene]] <- list(
    gene        = gene,
    mirnas      = selected_mirnas,
    source_tags = source_tags,
    source      = src_summary
  )

  cat(sprintf("  ✅ %s 最终选 %d 个 miRNA, 来源: %s\n",
              gene, length(selected_mirnas), src_summary))
}

all_mirnas <- unique(unlist(lapply(mirna_results, function(x) x$mirnas)))
cat(sprintf("\n总计 %d 个独特 miRNA (3 个 mRNA 去重后)\n", length(all_mirnas)))

cat("\n========== Step 2: 查询 miRNA 上游 lncRNA (ENCORI/Starbase) ==========\n")

MAX_LNCRNA_PER_MIRNA <- 5

query_encori_lncrna <- function(mirna_id, timeout_sec = 15) {
  url <- sprintf(
    "https://rnasysu.com/encori/api/miRNATarget/?assembly=hg19&geneType=lncRNA&miRNA=%s&clipExpNum=1&degraExpNum=0&pancancerNum=0&programNum=0&program=None&target=All&cellType=all",
    URLencode(mirna_id, reserved = TRUE)
  )

  result <- tryCatch({
    tmp <- tempfile()
    suppressWarnings(
      utils::download.file(url, tmp, quiet = TRUE, method = "auto",
                           cacheOK = FALSE)
    )
    if (file.exists(tmp) && file.info(tmp)$size > 0) {
      lines <- readLines(tmp, warn = FALSE)
      unlink(tmp)
      lines
    } else {
      NULL
    }
  }, error = function(e) NULL,
     warning = function(w) NULL)

  if (is.null(result)) return(character(0))

  data_lines <- result[!grepl("^#", result) & nchar(result) > 0]
  if (length(data_lines) < 2) return(character(0))

  header <- strsplit(data_lines[1], "\t")[[1]]
  gene_col <- which(header == "geneName")
  if (length(gene_col) == 0) {
    gene_col <- which(header %in% c("Gene", "gene_name", "geneSymbol"))
  }
  if (length(gene_col) == 0) return(character(0))

  lncrnas <- character(0)
  for (i in 2:length(data_lines)) {
    fields <- strsplit(data_lines[i], "\t")[[1]]
    if (length(fields) >= gene_col) {
      lncrnas <- c(lncrnas, fields[gene_col])
    }
  }
  unique(lncrnas[nchar(lncrnas) > 0])
}

encori_available <- FALSE
if (length(all_mirnas) > 0) {
  cat("测试 ENCORI API...\n")
  test_result <- query_encori_lncrna(all_mirnas[1])
  if (length(test_result) > 0) {
    encori_available <- TRUE
    cat("  ✅ ENCORI API 可用\n")
  } else {
    cat("  ⚠️ ENCORI API 不可用 (网络/防火墙问题), 回退内置库\n")
  }
}

lncrna_results <- list()

for (mirna in all_mirnas) {
  matched <- character(0)
  source_label <- "none"

  if (encori_available) {
    enc_lncrnas <- query_encori_lncrna(mirna)
    if (length(enc_lncrnas) > 0) {
      matched <- enc_lncrnas
      source_label <- "ENCORI/Starbase API"
    }
    Sys.sleep(0.3)
  }

  if (length(matched) == 0) {
    for (lnc in names(builtin_lncrna_database)) {
      if (mirna %in% builtin_lncrna_database[[lnc]]) {
        matched <- c(matched, lnc)
      }
    }
    if (length(matched) > 0) {
      source_label <- "Builtin literature database"
    }
  }

  if (length(matched) > MAX_LNCRNA_PER_MIRNA) {
    matched <- matched[1:MAX_LNCRNA_PER_MIRNA]
  }
  if (length(matched) > 0) {
    lncrna_results[[mirna]] <- list(
      lncrnas = matched,
      source  = source_label
    )
  }
}

cat(sprintf("有 lncRNA 匹配的 miRNA 数量: %d / %d\n",
            length(lncrna_results), length(all_mirnas)))

all_lncrnas <- unique(unlist(lapply(lncrna_results, function(x) x$lncrnas)))
cat(sprintf("总计 %d 个独特 lncRNA\n", length(all_lncrnas)))

if (encori_available) {
  cat("数据源: ENCORI/Starbase API (实时) + 内置库回退\n")
} else {
  cat("数据源: 内置库 (ENCORI/Starbase 文献整理)\n")
}

cat("\n========== Step 3: 构建网络边列表 ==========\n")

edges <- data.frame(from = character(), to = character(),
                    type = character(), stringsAsFactors = FALSE)

for (gene in names(mirna_results)) {
  for (mirna in mirna_results[[gene]]$mirnas) {
    edges <- rbind(edges, data.frame(
      from = mirna, to = gene, type = "miRNA-mRNA",
      stringsAsFactors = FALSE
    ))
  }
}

for (mirna in names(lncrna_results)) {
  for (lnc in lncrna_results[[mirna]]$lncrnas) {
    edges <- rbind(edges, data.frame(
      from = lnc, to = mirna, type = "lncRNA-miRNA",
      stringsAsFactors = FALSE
    ))
  }
}

all_nodes <- unique(c(edges$from, edges$to))
node_df <- data.frame(
  name = all_nodes,
  type = ifelse(all_nodes %in% biomarker_genes, "mRNA",
                ifelse(all_nodes %in% all_mirnas, "miRNA", "lncRNA")),
  stringsAsFactors = FALSE
)

cat(sprintf("节点统计: mRNA=%d, miRNA=%d, lncRNA=%d, 总=%d\n",
            sum(node_df$type == "mRNA"),
            sum(node_df$type == "miRNA"),
            sum(node_df$type == "lncRNA"),
            nrow(node_df)))
cat(sprintf("边数: %d (miRNA-mRNA=%d, lncRNA-miRNA=%d)\n",
            nrow(edges),
            sum(edges$type == "miRNA-mRNA"),
            sum(edges$type == "lncRNA-miRNA")))

write.csv(edges, file.path(res.path, "network_edges.csv"), row.names = FALSE)
write.csv(node_df, file.path(res.path, "network_nodes.csv"), row.names = FALSE)

cat("\n========== Step 4: 计算分层布局 ==========\n")

LAYER_Y <- c(lncRNA = 2.4, miRNA = 1.2, mRNA = 0)

n_mrna <- length(biomarker_genes)
mrna_x <- setNames(seq(0, 1, length.out = max(n_mrna, 2)), biomarker_genes)

mirna_x <- setNames(numeric(0), character(0))
mrna_gap <- if (n_mrna > 1) (1 / (n_mrna - 1)) else 1.0
cluster_width <- mrna_gap * 0.85

for (gene in biomarker_genes) {
  gene_mirnas <- mirna_results[[gene]]$mirnas
  base_x <- mrna_x[gene]
  n <- length(gene_mirnas)
  if (n == 1) {
    offsets <- 0
  } else {
    offsets <- seq(-cluster_width / 2, cluster_width / 2, length.out = n)
  }
  for (i in seq_along(gene_mirnas)) {
    m <- gene_mirnas[i]
    if (m %in% names(mirna_x)) {
      mirna_x[m] <- (mirna_x[m] + (base_x + offsets[i])) / 2
    } else {
      mirna_x[m] <- base_x + offsets[i]
    }
  }
}

if (length(mirna_x) > 1) {
  ord <- order(mirna_x)
  sorted_x <- mirna_x[ord]
  min_gap_mirna <- 0.06
  for (i in 2:length(sorted_x)) {
    if (sorted_x[i] - sorted_x[i - 1] < min_gap_mirna) {
      sorted_x[i] <- sorted_x[i - 1] + min_gap_mirna
    }
  }
  mirna_x[ord] <- sorted_x
}

lncrna_x <- setNames(numeric(0), character(0))
for (lnc in all_lncrnas) {
  connected_mirnas <- names(lncrna_results)[
    sapply(lncrna_results, function(x) lnc %in% x$lncrnas)
  ]
  if (length(connected_mirnas) > 0) {
    xs <- mirna_x[connected_mirnas]
    xs <- xs[!is.na(xs)]
    if (length(xs) > 0) {
      lncrna_x[lnc] <- mean(xs)
    }
  }
}

if (length(lncrna_x) > 1) {
  ord <- order(lncrna_x)
  sorted_x <- lncrna_x[ord]
  min_gap <- 0.12
  for (i in 2:length(sorted_x)) {
    if (sorted_x[i] - sorted_x[i - 1] < min_gap) {
      sorted_x[i] <- sorted_x[i - 1] + min_gap
    }
  }
  total_span <- sorted_x[length(sorted_x)] - sorted_x[1]
  if (total_span > 0) {
    target_span <- 1.0
    sorted_x <- (sorted_x - sorted_x[1]) / total_span * target_span
  }
  lncrna_x[ord] <- sorted_x
}

node_df$x <- NA_real_
node_df$y <- NA_real_
for (i in seq_len(nrow(node_df))) {
  nm <- node_df$name[i]
  tp <- node_df$type[i]
  node_df$y[i] <- LAYER_Y[tp]
  if (tp == "mRNA") {
    node_df$x[i] <- mrna_x[nm]
  } else if (tp == "miRNA") {
    node_df$x[i] <- mirna_x[nm]
  } else if (tp == "lncRNA") {
    node_df$x[i] <- lncrna_x[nm]
  }
}

na_idx <- which(is.na(node_df$x) | is.na(node_df$y))
if (length(na_idx) > 0) {
  node_df$x[na_idx] <- runif(length(na_idx), 0, 1)
  node_df$y[na_idx] <- 1
}

cat("\n========== Step 5: 绘制分层网络图 ==========\n")

g <- graph_from_data_frame(edges, vertices = node_df, directed = FALSE)

NODE_COLORS <- c(mRNA = "#E64B35", miRNA = "#4DBBD5", lncRNA = "#00A087")

layout_coords <- as.matrix(node_df[, c("x", "y")])
rownames(layout_coords) <- node_df$name

tg <- as_tbl_graph(g)
tg <- tg %>%
  activate(nodes) %>%
  mutate(
    x    = node_df$x[match(name, node_df$name)],
    y    = node_df$y[match(name, node_df$name)],
    type = node_df$type[match(name, node_df$name)]
  )

LABEL_SIZE <- 5

p_net <- ggraph(tg, layout = "manual",
                x = node_df$x[match(igraph::V(g)$name, node_df$name)],
                y = node_df$y[match(igraph::V(g)$name, node_df$name)]) +
  geom_edge_link(
    aes(color = factor(type)),
    edge_width = 0.5,
    edge_alpha = 0.6
  ) +
  scale_edge_color_manual(
    values = c("lncRNA-miRNA" = "#888888", "miRNA-mRNA" = "#555555"),
    name = "Edge type"
  ) +
  geom_node_point(
    aes(color = type, size = type),
    show.legend = TRUE
  ) +
  scale_color_manual(values = NODE_COLORS, name = "Node type") +
  scale_size_manual(values = c(mRNA = 11, miRNA = 6, lncRNA = 9), guide = "none") +

  geom_node_text(
    aes(label = ifelse(type == "mRNA", name, "")),
    nudge_y = -0.22,
    size = LABEL_SIZE,
    fontface = "bold",
    color = "black",
    show.legend = FALSE
  ) +

  geom_node_text(
    aes(label = ifelse(type == "lncRNA", name, "")),
    repel        = TRUE,
    nudge_y      = 0.25,
    direction    = "x",
    size         = LABEL_SIZE,
    fontface     = "bold",
    color        = "black",
    box.padding  = 0.4,
    point.padding = 0.3,
    force        = 2,
    segment.size = 0.2,
    segment.color = "gray60",
    segment.alpha = 0.6,
    min.segment.length = 0.1,
    max.overlaps = Inf,
    show.legend  = FALSE
  ) +

  geom_node_text(
    aes(label = ifelse(type == "miRNA",
                       gsub("^hsa-", "", name),
                       "")),
    repel        = TRUE,
    size         = LABEL_SIZE,
    fontface     = "bold",
    color        = "black",
    max.overlaps = Inf,
    box.padding  = 0.7,
    point.padding = 0.6,
    force        = 5,
    force_pull   = 0.3,
    segment.size = 0.25,
    segment.color = "gray50",
    segment.alpha = 0.7,
    min.segment.length = 0,
    show.legend  = FALSE
  ) +

  annotate("text", x = -0.18, y = LAYER_Y["lncRNA"], label = "lncRNA",
           color = "#00A087", size = 6, fontface = "bold", hjust = 0) +
  annotate("text", x = -0.18, y = LAYER_Y["miRNA"], label = "miRNA",
           color = "#4DBBD5", size = 6, fontface = "bold", hjust = 0) +
  annotate("text", x = -0.18, y = LAYER_Y["mRNA"], label = "mRNA",
           color = "#E64B35", size = 6, fontface = "bold", hjust = 0) +

  coord_cartesian(xlim = c(-0.25, 1.20), ylim = c(-0.6, 3.0), clip = "off") +
  labs(title = "lncRNA-miRNA-mRNA Regulatory Network",
       subtitle = "Three-layer ceRNA network of 5 stroke biomarkers (miRNA names abbreviated, hsa- prefix omitted)") +
  theme_void(base_size = 13) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key       = element_rect(fill = "white", color = NA),

    plot.title       = element_text(hjust = 0.5, size = 16, face = "bold", color = "black"),
    plot.subtitle    = element_text(hjust = 0.5, size = 11, color = "gray30"),
    legend.position  = "bottom",
    legend.box       = "horizontal",
    legend.title     = element_text(size = 12, face = "bold", color = "black"),
    legend.text      = element_text(size = 11, color = "black"),
    plot.margin      = margin(20, 30, 20, 30)
  )

ggsave(file.path(fig.path, "LncRNA_miRNA_mRNA_Network.pdf"),
       p_net, width = 18, height = 12)
ggsave(file.path(fig.path, "LncRNA_miRNA_mRNA_Network.png"),
       p_net, width = 18, height = 12, dpi = 300)
cat("✅ 分层网络图已保存\n")

cat("\n========== Step 6: 度数统计 ==========\n")

deg <- igraph::degree(g)
deg_df <- data.frame(
  Node   = names(deg),
  Degree = as.numeric(deg),
  Type   = node_df$type[match(names(deg), node_df$name)],
  stringsAsFactors = FALSE
)
deg_df <- deg_df[order(-deg_df$Degree), ]
write.csv(deg_df, file.path(res.path, "network_degree.csv"), row.names = FALSE)

cat("Top 10 hub 节点:\n")
print(head(deg_df, 10))

cat("\n========== Step 7: 度分布图 ==========\n")

p_deg <- ggplot(deg_df, aes(x = Degree, fill = Type)) +
  geom_histogram(binwidth = 1, alpha = 0.8, position = "dodge", color = "white") +
  scale_fill_manual(values = NODE_COLORS) +
  labs(title = "Node Degree Distribution",
       x = "Degree", y = "Count") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "top")

ggsave(file.path(fig.path, "Degree_Distribution.pdf"), p_deg, width = 8, height = 5)
ggsave(file.path(fig.path, "Degree_Distribution.png"), p_deg, width = 8, height = 5, dpi = 300)
cat("✅ 度分布图已保存\n")

cat("\n========== Step 8: Cytoscape 兼容文件 ==========\n")

igraph::write_graph(g, file.path(res.path, "network.graphml"), format = "graphml")
write.csv(edges, file.path(res.path, "network_for_cytoscape_edges.csv"), row.names = FALSE)
write.csv(node_df, file.path(res.path, "network_for_cytoscape_nodes.csv"), row.names = FALSE)
cat("✅ 已导出 GraphML + CSV\n")

source_info <- data.frame(
  Component = c("miRNA prediction", "lncRNA prediction"),
  Method = c(
    if (multiMir_available)
      "multiMiR miRTarBase (experimental, primary) + TargetScan ∪ miRDB (predicted, supplement) + literature fallback if needed; per scheme document alternative"
    else
      "Built-in literature database (fallback)",
    if (encori_available)
      "ENCORI/Starbase API (real-time query)"
    else
      "Built-in ENCORI/Starbase literature database (fallback)"
  ),
  Reference = c(
    "Ru Y et al. Bioinformatics 2014 (multiMiR); Huang HY et al. Nucleic Acids Res 2022 (miRTarBase); Agarwal V et al. eLife 2015 (TargetScan); Chen Y et al. Nucleic Acids Res 2020 (miRDB)",
    "Li JH et al. Nucleic Acids Res 2014 (ENCORI/Starbase)"
  ),
  stringsAsFactors = FALSE
)
write.csv(source_info, file.path(res.path, "data_sources.csv"), row.names = FALSE)

mrna_mirna_table <- data.frame()
for (gene in names(mirna_results)) {
  for (m in mirna_results[[gene]]$mirnas) {
    mrna_mirna_table <- rbind(mrna_mirna_table, data.frame(
      mRNA   = gene,
      miRNA  = m,
      Source = mirna_results[[gene]]$source,
      stringsAsFactors = FALSE
    ))
  }
}
write.csv(mrna_mirna_table, file.path(res.path, "mRNA_miRNA_pairs.csv"), row.names = FALSE)

mirna_lncrna_table <- data.frame()
for (m in names(lncrna_results)) {
  for (lnc in lncrna_results[[m]]$lncrnas) {
    mirna_lncrna_table <- rbind(mirna_lncrna_table, data.frame(
      miRNA  = m,
      lncRNA = lnc,
      Source = lncrna_results[[m]]$source,
      stringsAsFactors = FALSE
    ))
  }
}
write.csv(mirna_lncrna_table, file.path(res.path, "miRNA_lncRNA_pairs.csv"), row.names = FALSE)

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("    第十一部分 lncRNA-miRNA-mRNA 网络构建 完成\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("网络规模:\n")
cat(sprintf("  mRNA   : %d\n", sum(node_df$type == "mRNA")))
cat(sprintf("  miRNA  : %d\n", sum(node_df$type == "miRNA")))
cat(sprintf("  lncRNA : %d\n", sum(node_df$type == "lncRNA")))
cat(sprintf("  总节点 : %d\n", nrow(node_df)))
cat(sprintf("  总边数 : %d (mRNA-miRNA: %d, miRNA-lncRNA: %d)\n\n",
            nrow(edges),
            sum(edges$type == "miRNA-mRNA"),
            sum(edges$type == "lncRNA-miRNA")))

cat("数据来源:\n")
cat(sprintf("  miRNA  : %s\n",
            ifelse(multiMir_available,
                   "multiMiR miRTarBase (primary) + TargetScan ∪ miRDB / literature (supplement)",
                   "Built-in literature database (fallback)")))
cat(sprintf("  lncRNA : %s\n\n",
            ifelse(encori_available,
                   "ENCORI/Starbase API (real-time)",
                   "Built-in ENCORI/Starbase literature (fallback)")))

cat("输出文件:\n")
cat("  Figures/LncRNA_miRNA_mRNA_Network.pdf/png  (分层网络主图)\n")
cat("  Figures/Degree_Distribution.pdf/png         (度分布)\n")
cat("  Results/network_edges.csv                   (边列表)\n")
cat("  Results/network_nodes.csv                   (节点列表)\n")
cat("  Results/network_degree.csv                  (度统计)\n")
cat("  Results/mRNA_miRNA_pairs.csv                (mRNA-miRNA 对)\n")
cat("  Results/miRNA_lncRNA_pairs.csv              (miRNA-lncRNA 对)\n")
cat("  Results/network.graphml                     (Cytoscape 导入文件)\n")
cat("  Results/data_sources.csv                    (数据来源记录)\n\n")

cat("✅ 第十一部分分析完成\n")

sink(file.path(res.path, "session_info.txt"))
sessionInfo()
sink()
