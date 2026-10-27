
cat("\n========================================\n")
cat("15. 分子对接结果分析 (CB-Dock 2 后处理)\n")
cat("========================================\n\n")

rm(list = ls())

options(bitmapType = "cairo")
options(device = "pdf")

setwd("")

output_dir <- "./15_Docking"
res.path <- file.path(output_dir, "Results")
fig.path <- file.path(output_dir, "Figures")

for (path in c(output_dir, res.path, fig.path)) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

cat("\n========== 加载 R 包 ==========\n")

packages <- c("ggplot2", "dplyr", "tidyr", "stringr", "forcats",
              "pheatmap", "RColorBrewer", "ggrepel")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org", quiet = TRUE)
    library(pkg, character.only = TRUE)
  }
  cat("✅", pkg, "\n")
}

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate

cat("\n========== 7 个对接组合定义 ==========\n")

dock_pairs <- tibble::tribble(
  ~Protein,  ~Ligand,        ~PDB_ID,       ~Compound_NSC, ~Spearman_r, ~Mechanism,
  "CYP1B1",  "AZD6738",      "3PM0",        "780250",      0.486,       "ATR kinase inhibitor",
  "CYP1B1",  "Ouabain",      "3PM0",        "757425",      0.482,       "Na/K-ATPase inhibitor",
  "CYP1B1",  "Hypothemycin", "3PM0",        "354462",     -0.464,       "Covalent kinase inhibitor",
  "ADM",     "Ciclosporin",  "AlphaFold",   "290193",     -0.462,       "Calcineurin inhibitor",
  "CYP1B1",  "Rapamycin",    "3PM0",        "758664",      0.431,       "mTOR inhibitor",
  "ABCA1",   "Cordycepin",   "AlphaFold",   "63984",      -0.424,       "Adenosine analog",
  "HSPA1A",  "Streptozocin", "1HJO",        "85998",      -0.416,       "DNA alkylating agent"
)

print(dock_pairs)

cat("\n========== Step 1: 解析 CB-Dock 2 CurPockets ==========\n")

parse_curpockets <- function(txt_path) {
  if (!file.exists(txt_path)) {
    return(list(error = paste("File not found:", txt_path)))
  }

  lines <- readLines(txt_path, warn = FALSE)
  lines <- lines[nchar(trimws(lines)) > 0]

  if (length(lines) < 2) {
    return(list(error = "File too short"))
  }

  hdr <- strsplit(trimws(lines[1]), "\\s+")[[1]]
  score_col <- which(tolower(hdr) == "score")
  if (length(score_col) == 0) score_col <- 9
  vol_col <- which(tolower(hdr) %in% c("cavities_volume", "cavity_volume"))
  if (length(vol_col) == 0) vol_col <- 2

  scores <- c()
  volumes <- c()
  cavity_ids <- c()
  contact_residues <- c()
  for (i in 2:length(lines)) {
    parts <- strsplit(trimws(lines[i]), "\\s+")[[1]]
    if (length(parts) < score_col) next
    s <- suppressWarnings(as.numeric(parts[score_col]))
    v <- suppressWarnings(as.numeric(parts[vol_col]))
    cid <- parts[1]
    contacts <- if (length(parts) > score_col) {
      paste(parts[(score_col + 1):length(parts)], collapse = " ")
    } else { NA }
    if (!is.na(s)) {
      scores <- c(scores, s)
      volumes <- c(volumes, v)
      cavity_ids <- c(cavity_ids, cid)
      contact_residues <- c(contact_residues, contacts)
    }
  }

  if (length(scores) == 0) {
    return(list(error = "No valid scores parsed"))
  }

  best_idx <- which.min(scores)

  list(
    n_cavities = length(scores),
    all_scores = scores,
    all_volumes = volumes,
    all_cavity_ids = cavity_ids,
    best_score = scores[best_idx],
    best_cavity_id = cavity_ids[best_idx],
    best_cavity_volume = volumes[best_idx],
    best_contact_residues = contact_residues[best_idx],
    error = NA
  )
}

docking_base <- file.path(output_dir, "docking_results")
if (!dir.exists(docking_base)) {
  stop("❌ docking_results/ 目录不存在: ", docking_base)
}
cat("docking_results/ 目录存在, 文件数:",
    length(list.files(docking_base)), "\n")

dock_pairs$txt_path           <- NA_character_
dock_pairs$complex_path       <- NA_character_
dock_pairs$ligand_path        <- NA_character_
dock_pairs$best_vina          <- NA_real_
dock_pairs$best_cavity_id     <- NA_character_
dock_pairs$best_cavity_volume <- NA_real_
dock_pairs$n_cavities         <- NA_integer_
dock_pairs$all_scores_str     <- NA_character_
dock_pairs$contact_residues   <- NA_character_

for (i in seq_len(nrow(dock_pairs))) {
  pair_name <- paste0(dock_pairs$Protein[i], "_", dock_pairs$Ligand[i])

  txt_files <- list.files(docking_base,
    pattern = paste0("^", pair_name, ".*CurPockets.*\\.txt$"),
    full.names = TRUE)
  cmplx_files <- list.files(docking_base,
    pattern = paste0("^", pair_name, ".*complex.*\\.pdb$"),
    full.names = TRUE)
  lig_files <- list.files(docking_base,
    pattern = paste0("^", pair_name, ".*ligand.*\\.pdb$"),
    full.names = TRUE)

  if (length(txt_files) > 0)   dock_pairs$txt_path[i]     <- txt_files[1]
  if (length(cmplx_files) > 0) dock_pairs$complex_path[i] <- cmplx_files[1]
  if (length(lig_files) > 0)   dock_pairs$ligand_path[i]  <- lig_files[1]

  if (is.na(dock_pairs$txt_path[i])) {
    cat(sprintf("  [%d] %-25s: ❌ 未找到 CurPockets.txt\n", i, pair_name))
    next
  }

  result <- parse_curpockets(dock_pairs$txt_path[i])
  if (!is.null(result$error) && !is.na(result$error)) {
    cat(sprintf("  [%d] %-25s: ❌ 解析失败: %s\n", i, pair_name, result$error))
    next
  }

  dock_pairs$best_vina[i]          <- result$best_score
  dock_pairs$best_cavity_id[i]     <- result$best_cavity_id
  dock_pairs$best_cavity_volume[i] <- result$best_cavity_volume
  dock_pairs$n_cavities[i]         <- result$n_cavities
  dock_pairs$all_scores_str[i]     <- paste(round(result$all_scores, 2), collapse = "; ")
  dock_pairs$contact_residues[i]   <- result$best_contact_residues

  cat(sprintf("  [%d] %-25s: Best Vina = %5.2f kcal/mol  (cavity %s, %d Å³)\n",
              i, pair_name, result$best_score, result$best_cavity_id,
              result$best_cavity_volume))
}

cat("\n========== Step 2: 结合能汇总表 (Table 3) ==========\n")

n_complete <- sum(!is.na(dock_pairs$best_vina))
cat(sprintf("成功解析: %d / %d\n\n", n_complete, nrow(dock_pairs)))

classify_strength <- function(v) {
  if (is.na(v)) return("N/A")
  if (v <= -8.5) return("Very Strong")
  if (v <= -7.5) return("Strong")
  if (v <= -6.0) return("Moderate")
  return("Weak")
}
dock_pairs$Binding_Strength <- sapply(dock_pairs$best_vina, classify_strength)

table3 <- dock_pairs %>%
  select(Protein, Ligand, PDB_ID, Compound_NSC, Mechanism,
         Spearman_r, best_vina, best_cavity_id, best_cavity_volume,
         Binding_Strength, all_scores_str)
colnames(table3) <- c("Target_Protein", "Compound", "PDB_ID", "NSC",
                      "Mechanism", "CellMiner_r", "Vina_Score_kcal_mol",
                      "Best_Cavity", "Cavity_Volume_A3",
                      "Binding_Strength", "All_5_Pose_Scores")

print(table3)
write.csv(table3, file.path(res.path, "binding_energy_summary.csv"),
          row.names = FALSE)

table3_simple <- table3 %>%
  select(Target_Protein, Compound, PDB_ID, Vina_Score_kcal_mol,
         CellMiner_r, Binding_Strength)
write.csv(table3_simple, file.path(res.path, "Table3_simplified.csv"),
          row.names = FALSE)

cat("\n========== Step 3: 结合能柱状图 ==========\n")

valid_df <- dock_pairs %>% filter(!is.na(best_vina))

if (nrow(valid_df) > 0) {
  valid_df$pair_label <- paste0(valid_df$Protein, " - ", valid_df$Ligand)
  valid_df <- valid_df %>% arrange(best_vina)
  valid_df$pair_label <- factor(valid_df$pair_label, levels = valid_df$pair_label)

  protein_colors <- c(
    CYP1B1 = "#E64B35",
    HSPA1A = "#3C5488",
    ABCA1  = "#00A087",
    ADM    = "#4DBBD5",
    LRG1   = "#F39B7F"
  )

  p_bar <- ggplot(valid_df,
                  aes(x = pair_label, y = best_vina, fill = Protein)) +
    geom_bar(stat = "identity", width = 0.65,
             color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f", best_vina)),
              vjust = 1.5, color = "white",
              size = 5, fontface = "bold") +
    geom_hline(yintercept = -7, linetype = "dashed",
               color = "gray40", linewidth = 0.4) +
    geom_hline(yintercept = -8.5, linetype = "dashed",
               color = "darkred", linewidth = 0.4) +
    scale_fill_manual(values = protein_colors, name = "Target protein") +
    scale_y_continuous(
      breaks = seq(-12, 0, by = 2),
      expand = expansion(mult = c(0.08, 0.02))
    ) +
    labs(title = "Molecular Docking: Binding Affinities",
         subtitle = sprintf("AutoDock Vina (CB-Dock 2) blind docking, %d protein-ligand pairs",
                            nrow(valid_df)),
         x = NULL,
         y = "Binding energy (kcal/mol)") +
    theme_bw(base_size = 13) +
    theme(plot.title    = element_text(hjust = 0.5, size = 15, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray40"),
          axis.text.x   = element_text(angle = 30, hjust = 1, size = 11, face = "bold"),
          axis.text.y   = element_text(size = 11),
          legend.position = "right",
          legend.title  = element_text(size = 11, face = "bold"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor   = element_blank())

  ggsave(file.path(fig.path, "Binding_Energy_Barplot.pdf"),
         p_bar, width = 11, height = 6)
  ggsave(file.path(fig.path, "Binding_Energy_Barplot.png"),
         p_bar, width = 11, height = 6, dpi = 300)
  cat("✅ 结合能柱状图已保存\n")
}

cat("\n========== Step 4: 结合能热图 ==========\n")

if (nrow(valid_df) > 0) {
  proteins <- unique(dock_pairs$Protein)
  ligands  <- unique(dock_pairs$Ligand)

  be_mat <- matrix(NA_real_,
                   nrow = length(proteins),
                   ncol = length(ligands),
                   dimnames = list(proteins, ligands))
  for (i in seq_len(nrow(valid_df))) {
    be_mat[valid_df$Protein[i], valid_df$Ligand[i]] <- valid_df$best_vina[i]
  }

  breaks <- seq(min(be_mat, na.rm = TRUE) - 0.5,
                max(be_mat, na.rm = TRUE) + 0.5,
                length.out = 101)

  draw_hm <- function() {
    pheatmap(be_mat,
             color = colorRampPalette(c("#A50F15", "#FB6A4A", "#FDD49E", "#FFFFFF"))(100),
             breaks = breaks,
             main = "Docking Binding Energy (kcal/mol)\n(more negative = stronger binding)",
             display_numbers = TRUE,
             number_format = "%.1f",
             fontsize_number = 12,
             fontsize_row = 12,
             fontsize_col = 11,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             border_color = "gray80",
             angle_col = 45,
             cellwidth = 70,
             cellheight = 50,
             na_col = "gray95")
  }

  pdf(file.path(fig.path, "Binding_Energy_Heatmap.pdf"), width = 12, height = 5)
  draw_hm()
  dev.off()
  png(file.path(fig.path, "Binding_Energy_Heatmap.png"),
      width = 12, height = 5, units = "in", res = 300)
  draw_hm()
  dev.off()
  cat("✅ 结合能热图已保存\n")
}

cat("\n========== Step 5: Vina score vs CellMiner |r| 散点图 ==========\n")

if (nrow(valid_df) > 0) {
  valid_df$abs_r <- abs(valid_df$Spearman_r)

  p_scatter <- ggplot(valid_df,
                      aes(x = abs_r, y = best_vina, color = Protein, label = Ligand)) +
    geom_point(size = 5, alpha = 0.85) +
    ggrepel::geom_text_repel(size = 4, fontface = "bold",
                             box.padding = 0.5, point.padding = 0.4,
                             segment.color = "gray60",
                             show.legend = FALSE) +
    scale_color_manual(values = c(
      CYP1B1 = "#E64B35", HSPA1A = "#3C5488",
      ABCA1  = "#00A087", ADM    = "#4DBBD5"
    )) +
    labs(title = "CellMiner Correlation vs Docking Binding Energy",
         x = "CellMiner |Spearman r|",
         y = "Vina binding energy (kcal/mol)") +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "right")

  ggsave(file.path(fig.path, "Vina_vs_CellMiner.pdf"),
         p_scatter, width = 9, height = 6)
  ggsave(file.path(fig.path, "Vina_vs_CellMiner.png"),
         p_scatter, width = 9, height = 6, dpi = 300)
  cat("✅ 散点图已保存\n")
}

cat("\n========== 生成报告 ==========\n")

report_lines <- c(
  "# 第十五部分 (分子对接) 报告",
  paste("生成时间:", Sys.time()),
  "",
  "## 对接策略",
  "每个 ADMET 筛选出的 Tier 2 候选化合物与其 CellMiner 最强相关基因的蛋白结构进行 blind docking。",
  "",
  "## 蛋白结构来源",
  "- CYP1B1: PDB 3PM0 (X-ray, 2.7 Å, 含 HEM 血红素辅基)",
  "- HSPA1A: PDB 1HJO (X-ray, 1.84 Å, 含 ADP, NBD ATPase 结构域)",
  "- ABCA1 : AlphaFold AF-O95477 (全长 2261 aa, 膜蛋白)",
  "- ADM   : AlphaFold AF-P35318 (185 aa, 信号肽 + 成熟肽)",
  "",
  "## 对接工具",
  "- 工具: CB-Dock 2 (http://clab.labshare.cn/cb-dock2/)",
  "- 算法: AutoDock Vina + CurPocket cavity detection",
  "- 模式: Auto Blind Docking",
  "",
  "## 结合能结果"
)

if (n_complete > 0) {
  report_lines <- c(report_lines, "",
    "| Protein | Ligand | Vina (kcal/mol) | Cavity | Volume (A3) | |r| | Strength |",
    "|---|---|---|---|---|---|---|")
  ordered <- valid_df %>% arrange(best_vina)
  for (i in seq_len(nrow(ordered))) {
    report_lines <- c(report_lines,
      sprintf("| %s | %s | **%.2f** | %s | %d | %.3f | %s |",
              ordered$Protein[i], ordered$Ligand[i],
              ordered$best_vina[i], ordered$best_cavity_id[i],
              ordered$best_cavity_volume[i],
              abs(ordered$Spearman_r[i]),
              ordered$Binding_Strength[i]))
  }
}

report_lines <- c(report_lines,
  "",
  "## 输出文件",
  "",
  "### Results/",
  "- binding_energy_summary.csv  (全部信息, 11 列)",
  "- Table3_simplified.csv",
  "",
  "### Figures/",
  "- Binding_Energy_Barplot.pdf/png  (柱状图, 按 Vina 排序)",
  "- Binding_Energy_Heatmap.pdf/png  (蛋白×配体热图)",
  "- Vina_vs_CellMiner.pdf/png       (Vina vs |r| 散点图)",
  "",
  "## 下一步",
  "1. PyMoL 可视化 Top 3 强结合对接复合物 (论文 Fig.13a)"
)

writeLines(report_lines, file.path(res.path, "Docking_Report.md"))

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("    第十五部分 (分子对接分析) 完成\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat(sprintf("对接组合: %d 个 (成功解析: %d 个)\n\n", nrow(dock_pairs), n_complete))

if (n_complete > 0) {
  cat("结合能汇总 (按 Vina 从最强到最弱):\n")
  ordered_print <- valid_df %>%
    arrange(best_vina) %>%
    select(Protein, Ligand, best_vina, best_cavity_id,
           best_cavity_volume, Binding_Strength)
  print(ordered_print)
}

cat("\n输出文件:\n")
cat("  Results/binding_energy_summary.csv\n")
cat("  Results/Table3_simplified.csv\n")
cat("  Figures/Binding_Energy_Barplot.pdf/png\n")
cat("  Figures/Binding_Energy_Heatmap.pdf/png\n")
cat("  Figures/Vina_vs_CellMiner.pdf/png\n")

cat("\n✅ 第十五部分 (分子对接) 数据分析完成\n")

sink(file.path(res.path, "session_info.txt"))
sessionInfo()
sink()
