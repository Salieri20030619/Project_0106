
cat("\n========================================\n")
cat("14b. ADMET 筛选 v4 (ADMETlab 3.0, 策略 B 放宽)\n")
cat("========================================\n\n")

rm(list = ls())

options(bitmapType = "cairo")
options(device = "pdf")

setwd("")

output_dir <- "./14_Drug_Prediction"
res.path <- file.path(output_dir, "Results")
fig.path <- file.path(output_dir, "Figures")

cat("========== 加载 R 包 ==========\n")

packages <- c("ggplot2", "dplyr", "tidyr", "stringr", "forcats", "pheatmap", "RColorBrewer")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org", quiet = TRUE)
    library(pkg, character.only = TRUE)
  }
  cat("OK", pkg, "\n")
}

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate

cat("\n========== 读取输入文件 ==========\n")

admet_candidates <- c(
  file.path(res.path, "admet_result.csv"),
  list.files(res.path, pattern = "^tmp.*\\.csv$", full.names = TRUE),
  list.files(res.path, pattern = "admet.*\\.csv", full.names = TRUE, ignore.case = TRUE)
)

admet_file <- NULL
for (f in unique(admet_candidates)) {
  if (!file.exists(f)) next
  first_line <- readLines(f, n = 1)
  if (grepl("raw_smiles|^smiles", first_line, ignore.case = TRUE)) {
    admet_file <- f
    break
  }
}

if (is.null(admet_file)) {
  stop("ADMETlab 结果文件未找到. 请把 ADMETlab 下载的 CSV 放到 Results/admet_result.csv")
}
cat("OK 使用 ADMET 文件:", admet_file, "\n")

target_path <- file.path(res.path, "admet_result.csv")
if (admet_file != target_path) {
  file.copy(admet_file, target_path, overwrite = TRUE)
  cat("已归档到:", target_path, "\n")
}

admet_raw <- read.csv(admet_file, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("ADMET 数据: %d 行 x %d 列\n", nrow(admet_raw), ncol(admet_raw)))

candidates_file <- file.path(res.path, "candidate_drugs_with_SMILES.csv")
if (!file.exists(candidates_file)) {
  stop("candidate_drugs_with_SMILES.csv 不存在")
}
candidates <- read.csv(candidates_file, stringsAsFactors = FALSE)
cat(sprintf("CellMiner 候选: %d 行\n", nrow(candidates)))

sig_pairs_file <- file.path(res.path, "significant_drug_gene_pairs.csv")
sig_pairs <- read.csv(sig_pairs_file, stringsAsFactors = FALSE)
cat(sprintf("显著药物-基因对: %d 行\n", nrow(sig_pairs)))

cat("\n========== ADMET 列名自动检测 ==========\n")

detect_col <- function(patterns, cols) {
  for (p in patterns) {
    hits <- grep(p, cols, ignore.case = FALSE, value = TRUE)
    if (length(hits) > 0) return(hits[1])
  }
  for (p in patterns) {
    hits <- grep(p, cols, ignore.case = TRUE, value = TRUE)
    if (length(hits) > 0) return(hits[1])
  }
  return(NA_character_)
}

all_cols <- colnames(admet_raw)

col_smiles   <- detect_col(c("^raw_smiles$", "^smiles$", "^SMILES$"), all_cols)
col_hia      <- detect_col(c("^hia$", "^HIA$"), all_cols)
col_bbb      <- detect_col(c("^BBB$", "^bbb$"), all_cols)
col_1a2      <- detect_col(c("^CYP1A2-inh$", "^CYP1A2_inh$"), all_cols)
col_2c9      <- detect_col(c("^CYP2C9-inh$", "^CYP2C9_inh$"), all_cols)
col_2c19     <- detect_col(c("^CYP2C19-inh$", "^CYP2C19_inh$"), all_cols)
col_2d6      <- detect_col(c("^CYP2D6-inh$", "^CYP2D6_inh$"), all_cols)
col_3a4      <- detect_col(c("^CYP3A4-inh$", "^CYP3A4_inh$"), all_cols)
col_herg     <- detect_col(c("^hERG$", "^HERG$"), all_cols)
col_ames     <- detect_col(c("^Ames$", "^AMES$"), all_cols)
col_lipinski <- detect_col(c("^Lipinski$", "^lipinski$"), all_cols)

cat("检测到的列:\n")
cat(sprintf("  SMILES   = %s\n", col_smiles))
cat(sprintf("  hia      = %s\n", col_hia))
cat(sprintf("  BBB      = %s\n", col_bbb))
cat(sprintf("  CYP1A2   = %s\n", col_1a2))
cat(sprintf("  CYP2C9   = %s\n", col_2c9))
cat(sprintf("  CYP2C19  = %s\n", col_2c19))
cat(sprintf("  CYP2D6   = %s\n", col_2d6))
cat(sprintf("  CYP3A4   = %s\n", col_3a4))
cat(sprintf("  hERG     = %s\n", col_herg))
cat(sprintf("  Ames     = %s\n", col_ames))
cat(sprintf("  Lipinski = %s\n", col_lipinski))

required <- list(SMILES = col_smiles, hia = col_hia, BBB = col_bbb,
                  hERG = col_herg, Ames = col_ames, Lipinski = col_lipinski)
missing_keys <- names(required)[sapply(required, is.na)]
if (length(missing_keys) > 0) {
  cat(sprintf("WARN 关键列缺失: %s\n", paste(missing_keys, collapse = ", ")))
  cat("可用列 (前 40 个):\n")
  print(head(all_cols, 40))
  stop("请检查 CSV 格式")
}

cat("\n========== 提取关键 ADMET 端点 ==========\n")

admet_df <- data.frame(
  SMILES       = as.character(admet_raw[[col_smiles]]),
  HIA          = as.numeric(admet_raw[[col_hia]]),
  BBB          = as.numeric(admet_raw[[col_bbb]]),
  CYP1A2_inh   = as.numeric(admet_raw[[col_1a2]]),
  CYP2C9_inh   = as.numeric(admet_raw[[col_2c9]]),
  CYP2C19_inh  = as.numeric(admet_raw[[col_2c19]]),
  CYP2D6_inh   = as.numeric(admet_raw[[col_2d6]]),
  CYP3A4_inh   = as.numeric(admet_raw[[col_3a4]]),
  hERG         = as.numeric(admet_raw[[col_herg]]),
  Ames         = as.numeric(admet_raw[[col_ames]]),
  Lipinski_raw = admet_raw[[col_lipinski]],
  stringsAsFactors = FALSE
)

cat("ADMET 端点数值摘要:\n")
print(summary(admet_df[, c("HIA", "BBB", "hERG", "Ames")]))
cat(sprintf("\nLipinski 原始取值:\n"))
print(table(admet_df$Lipinski_raw, useNA = "ifany"))

lipinski_violations <- rep(NA_real_, nrow(admet_df))
lip_vals <- admet_df$Lipinski_raw

if (all(is.na(lip_vals) | lip_vals %in% c("0", "1", "2", "3", "4", "5",
                                           0, 1, 2, 3, 4, 5))) {
  lipinski_violations <- suppressWarnings(as.numeric(lip_vals))
  cat("\nOK Lipinski 解析为数字 (违反规则数)\n")
} else if (any(grepl("accept|reject", lip_vals, ignore.case = TRUE))) {
  lipinski_violations <- ifelse(
    tolower(trimws(as.character(lip_vals))) == "accepted", 0,
    ifelse(tolower(trimws(as.character(lip_vals))) == "rejected", 1, NA)
  )
  cat("\nOK Lipinski 解析为字符串 (Accepted=0, Rejected=1)\n")
} else {
  lipinski_violations <- suppressWarnings(as.numeric(lip_vals))
  cat("\nWARN Lipinski 格式未知, 尝试数字解析\n")
}

admet_df$Lipinski_violations <- lipinski_violations
admet_df$pass_Lipinski       <- admet_df$Lipinski_violations == 0
admet_df$pass_Lipinski_relax <- admet_df$Lipinski_violations <= 1

cat(sprintf("Lipinski = 0 (通过): %d / %d\n",
            sum(admet_df$pass_Lipinski, na.rm = TRUE), nrow(admet_df)))

admet_df$pass_HIA       <- admet_df$HIA  < 0.3
admet_df$pass_BBB       <- admet_df$BBB  < 0.3
admet_df$pass_hERG      <- admet_df$hERG < 0.3
admet_df$pass_AMES_strict  <- admet_df$Ames < 0.5
admet_df$pass_AMES_relax   <- admet_df$Ames < 0.7

admet_df$n_CYP_clean <- (admet_df$CYP1A2_inh  < 0.5) +
                        (admet_df$CYP2C9_inh  < 0.5) +
                        (admet_df$CYP2C19_inh < 0.5) +
                        (admet_df$CYP2D6_inh  < 0.5) +
                        (admet_df$CYP3A4_inh  < 0.5)
admet_df$pass_CYP_all     <- admet_df$n_CYP_clean == 5
admet_df$pass_CYP_relaxed <- admet_df$n_CYP_clean >= 3

admet_df$Tier1 <- with(admet_df,
  pass_HIA & pass_BBB & pass_CYP_all & pass_hERG & pass_AMES_strict & pass_Lipinski)

n_nonsafety_pass <- with(admet_df, pass_HIA + pass_BBB + pass_CYP_relaxed)
admet_df$Tier2 <- admet_df$pass_hERG &
                  admet_df$pass_AMES_relax &
                  (n_nonsafety_pass >= 2) &
                  admet_df$pass_Lipinski_relax

admet_df$Tier3 <- admet_df$pass_hERG & admet_df$pass_AMES_relax & admet_df$pass_Lipinski_relax

cat("\n========== 合并 CellMiner 信息 ==========\n")

admet_with_drug <- admet_df %>%
  left_join(candidates %>% select(SMILES, Drug_Name, Drug_NSC, max_abs_cor),
            by = "SMILES")

n_unmatched <- sum(is.na(admet_with_drug$Drug_Name))
cat(sprintf("SMILES 直接匹配: %d / %d\n",
            nrow(admet_with_drug) - n_unmatched, nrow(admet_with_drug)))

if (n_unmatched > nrow(admet_with_drug) * 0.3) {
  cat(sprintf("WARN %d 行未匹配 (>30%%), 尝试按顺序对齐\n", n_unmatched))
  if (nrow(admet_df) == nrow(candidates)) {
    admet_with_drug$Drug_Name   <- candidates$Drug_Name[seq_len(nrow(admet_df))]
    admet_with_drug$Drug_NSC    <- candidates$Drug_NSC[seq_len(nrow(admet_df))]
    admet_with_drug$max_abs_cor <- candidates$max_abs_cor[seq_len(nrow(admet_df))]
    cat("  OK 按顺序对齐成功\n")
  } else {
    cat(sprintf("  WARN 行数不匹配 (ADMET: %d, candidates: %d)\n",
                nrow(admet_df), nrow(candidates)))
  }
}

still_unmatched <- is.na(admet_with_drug$Drug_Name)
if (any(still_unmatched)) {
  admet_with_drug$Drug_Name[still_unmatched] <-
    paste0("Unknown_", which(still_unmatched))
}

write.csv(admet_with_drug,
          file.path(res.path, "admet_full_results.csv"),
          row.names = FALSE)

cat("\n========== 筛选漏斗 ==========\n")

funnel <- data.frame(
  Step = c("All compounds",
           "hia good (<0.3)",
           "BBB- (<0.3)",
           "CYP1A2 clean (<0.5)",
           "CYP2C9 clean (<0.5)",
           "CYP2C19 clean (<0.5)",
           "CYP2D6 clean (<0.5)",
           "CYP3A4 clean (<0.5)",
           "CYP >=3/5 clean (relaxed)",
           "hERG low (<0.3)*",
           "AMES <0.5 (strict)",
           "AMES <0.7 (relaxed)*",
           "Lipinski = 0 (strict)",
           "Lipinski <= 1 (relaxed)",
           "Tier 1 PASS (all strict scheme)",
           "Tier 2 PASS (safety-core + relaxed)",
           "Tier 3 PASS (hERG+AMES+Lipinski)"),
  Count = c(nrow(admet_df),
            sum(admet_df$pass_HIA, na.rm = TRUE),
            sum(admet_df$pass_BBB, na.rm = TRUE),
            sum(admet_df$CYP1A2_inh  < 0.5, na.rm = TRUE),
            sum(admet_df$CYP2C9_inh  < 0.5, na.rm = TRUE),
            sum(admet_df$CYP2C19_inh < 0.5, na.rm = TRUE),
            sum(admet_df$CYP2D6_inh  < 0.5, na.rm = TRUE),
            sum(admet_df$CYP3A4_inh  < 0.5, na.rm = TRUE),
            sum(admet_df$pass_CYP_relaxed, na.rm = TRUE),
            sum(admet_df$pass_hERG, na.rm = TRUE),
            sum(admet_df$pass_AMES_strict, na.rm = TRUE),
            sum(admet_df$pass_AMES_relax, na.rm = TRUE),
            sum(admet_df$pass_Lipinski, na.rm = TRUE),
            sum(admet_df$pass_Lipinski_relax, na.rm = TRUE),
            sum(admet_df$Tier1, na.rm = TRUE),
            sum(admet_df$Tier2, na.rm = TRUE),
            sum(admet_df$Tier3, na.rm = TRUE)),
  stringsAsFactors = FALSE
)

print(funnel)
write.csv(funnel, file.path(res.path, "admet_filter_funnel.csv"), row.names = FALSE)

cat(sprintf("\nTier 1 (方案严格):   %d\n", sum(admet_df$Tier1, na.rm = TRUE)))
cat(sprintf("Tier 2 (底线+放宽):  %d  <- 论文最终候选\n", sum(admet_df$Tier2, na.rm = TRUE)))
cat(sprintf("Tier 3 (核心安全):   %d\n", sum(admet_df$Tier3, na.rm = TRUE)))

cat("\n========== Tier 2 通过化合物 ==========\n")

tier2_pass <- admet_with_drug %>%
  filter(Tier2) %>%
  arrange(desc(max_abs_cor))

cat(sprintf("Tier 2 通过化合物数: %d\n", nrow(tier2_pass)))

if (nrow(tier2_pass) > 0) {
  print_cols <- intersect(c("Drug_Name", "Drug_NSC", "HIA", "BBB", "n_CYP_clean",
                            "hERG", "Ames", "Lipinski_violations", "max_abs_cor"),
                          colnames(tier2_pass))
  print(tier2_pass[, print_cols])
  write.csv(tier2_pass, file.path(res.path, "ADMET_filtered_drugs_Tier2.csv"),
            row.names = FALSE)

  if ("Drug_NSC" %in% colnames(tier2_pass) && "Drug_NSC" %in% colnames(sig_pairs)) {
    tier2_with_genes <- sig_pairs %>%
      inner_join(tier2_pass %>% select(Drug_NSC) %>% distinct(),
                 by = "Drug_NSC") %>%
      arrange(Gene, desc(abs(Cor)))
    write.csv(tier2_with_genes,
              file.path(res.path, "ADMET_filtered_drug_gene_pairs.csv"),
              row.names = FALSE)
    cat(sprintf("Tier 2 对应的药物-基因对: %d\n", nrow(tier2_with_genes)))
  }
}

cat("\n========== Fig.14a: 筛选漏斗图 ==========\n")

funnel_plot_df <- funnel
funnel_plot_df$Step <- factor(funnel_plot_df$Step, levels = funnel_plot_df$Step)
funnel_plot_df$Pct <- funnel_plot_df$Count / funnel_plot_df$Count[1] * 100
funnel_plot_df$Color <- ifelse(funnel_plot_df$Step == "Tier 2 PASS (safety-core + relaxed)", "#E64B35",
                         ifelse(funnel_plot_df$Step == "All compounds", "#3C5488",
                         ifelse(grepl("Tier", funnel_plot_df$Step), "#F39B7F",
                         ifelse(grepl("\\*", funnel_plot_df$Step), "#DC0000", "#4DBBD5"))))

p_funnel <- ggplot(funnel_plot_df,
                   aes(x = Count, y = fct_rev(Step), fill = Color)) +
  geom_bar(stat = "identity", width = 0.72) +
  geom_text(aes(label = sprintf("%d (%.0f%%)", Count, Pct)),
            hjust = -0.08, size = 3.6, fontface = "bold") +
  scale_fill_identity() +
  scale_x_continuous(expand = expansion(mult = c(0, 0.28))) +
  labs(title = "ADMET Screening Funnel",
       subtitle = sprintf("Filtering %d candidates; * = mandatory safety baseline",
                          nrow(admet_df)),
       x = "Number of compounds", y = NULL) +
  theme_bw(base_size = 12) +
  theme(plot.title    = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
        axis.text.y   = element_text(size = 10, face = "bold", color = "black"),
        axis.text.x   = element_text(size = 9),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank())

ggsave(file.path(fig.path, "ADMET_Filter_Funnel.pdf"),
       p_funnel, width = 11, height = 8)
ggsave(file.path(fig.path, "ADMET_Filter_Funnel.png"),
       p_funnel, width = 11, height = 8, dpi = 300)
cat("OK ADMET 筛选漏斗图已保存\n")

cat("\n========== Fig.14b: ADMET 端点热图 ==========\n")

if (nrow(tier2_pass) > 0) {
  hm_data <- tier2_pass %>%
    select(Drug_Name, HIA, BBB, CYP1A2_inh, CYP2C9_inh, CYP2C19_inh,
           CYP2D6_inh, CYP3A4_inh, hERG, Ames) %>%
    distinct()

  rownames(hm_data) <- make.unique(as.character(hm_data$Drug_Name))
  hm_mat <- as.matrix(hm_data[, -1])

  draw_admet_hm <- function() {
    pheatmap(hm_mat,
             color = colorRampPalette(c("#00A087", "white", "#E64B35"))(100),
             breaks = seq(0, 1, length.out = 101),
             main = "ADMET Profile of Tier-2 Filtered Drugs\n(green = safe, red = risk)",
             display_numbers = round(hm_mat, 2),
             number_format = "%.2f",
             fontsize_number = 10,
             fontsize_row = 11,
             fontsize_col = 11,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             border_color = "gray80",
             angle_col = 45,
             cellwidth = 40,
             cellheight = 24)
  }

  h_inch <- max(4, nrow(hm_mat) * 0.55 + 3)
  pdf(file.path(fig.path, "ADMET_Profile_Heatmap.pdf"),
      width = 11, height = h_inch)
  draw_admet_hm()
  dev.off()

  png(file.path(fig.path, "ADMET_Profile_Heatmap.png"),
      width = 11, height = h_inch, units = "in", res = 300)
  draw_admet_hm()
  dev.off()
  cat("OK ADMET 端点热图已保存\n")
} else {
  cat("WARN 无 Tier 2 化合物, 跳过热图\n")
}

report_lines <- c(
  "# 第十四部分 报告",
  paste("生成时间:", Sys.time()),
  paste("ADMET 数据源:", admet_file),
  "",
  "## 方法学要点 (策略 B 放宽版)",
  "",
  "方案原文要求 hERG + Ames 为'必须满足的安全性底线'。本研究保留 hERG 严格阈值",
  "(< 0.3, QT 延长为致命风险不可放宽), 同时对 Ames 阈值从 < 0.5 放宽至 < 0.7。",
  "放宽理由:",
  "  (a) ADMETlab 3.0 对 Ames 概率 0.5-0.7 区间定义为'中等风险', 不等于明确阳性;",
  "  (b) NCI-60 候选药物以细胞毒抗癌化合物为主, 其抗癌机制本身涉及 DNA 损伤,",
  "      Ames 阳性在此类化合物中难以完全避免;",
  "  (c) 心脏安全 (hERG) 是 FDA 药物下架的主要原因之一, 优先级高于 Ames。",
  "",
  "## 筛选阈值",
  "- 吸收 (hia < 0.3): P(HIA+ 吸收差) < 0.3 = 大概率高吸收",
  "- 分布 (BBB < 0.3): P(穿透) < 0.3 = 不穿透",
  "- 代谢 (CYP >= 3/5): 5 种 CYP 中至少 3 项抑制概率 < 0.5",
  "- **毒性 hERG < 0.3 (严格底线)**",
  "- **毒性 Ames < 0.7 (放宽底线)**",
  "- 类药性 Lipinski <= 1 (允许 1 条轻微违反)",
  "",
  "## 三层筛选策略",
  "- Tier 1 (方案严格版): 5 项 ADMET 全部通过 + Lipinski = 0",
  "- **Tier 2 (方案底线版, 放宽实用, 论文用)**: hERG 严格 + Ames 放宽 + 非底线 >= 2/3",
  "- Tier 3 (核心安全): hERG + Ames(relax) + Lipinski(relax)",
  "",
  "## 筛选结果",
  paste("- 输入化合物:", nrow(admet_df)),
  paste("- Tier 1 (方案严格):", sum(admet_df$Tier1, na.rm = TRUE)),
  paste("- **Tier 2 (方案底线+放宽, 推荐):", sum(admet_df$Tier2, na.rm = TRUE), "**  <- 论文用"),
  paste("- Tier 3 (核心安全):", sum(admet_df$Tier3, na.rm = TRUE)),
  ""
)

if (nrow(tier2_pass) > 0) {
  report_lines <- c(report_lines, "## Tier 2 通过的最终候选化合物")
  for (i in seq_len(nrow(tier2_pass))) {
    d <- tier2_pass[i, ]
    report_lines <- c(report_lines,
      sprintf("- **%s** (NSC %s): hia=%.2f, BBB=%.2f, CYP_clean=%d/5, hERG=%.2f, Ames=%.2f, Lipinski=%s, max|cor|=%.3f",
              d$Drug_Name,
              ifelse(is.na(d$Drug_NSC), "-", as.character(d$Drug_NSC)),
              d$HIA, d$BBB, d$n_CYP_clean, d$hERG, d$Ames,
              as.character(d$Lipinski_violations),
              ifelse(is.na(d$max_abs_cor), 0, d$max_abs_cor)))
  }
}

report_lines <- c(report_lines,
  "",
  "## 输出文件",
  "### Results/",
  "- admet_full_results.csv",
  "- admet_filter_funnel.csv",
  "- ADMET_filtered_drugs_Tier2.csv",
  "- ADMET_filtered_drug_gene_pairs.csv",
  "### Figures/",
  "- ADMET_Filter_Funnel.pdf/png",
  "- ADMET_Profile_Heatmap.pdf/png",
  "",
  "## 下一步",
  "基于 Tier 2 化合物进入第十五部分 (CB-Dock / AutoDock Vina 分子对接)"
)

writeLines(report_lines, file.path(res.path, "ADMET_Filter_Report.md"))

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("    第十四部分 完成\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("筛选漏斗:\n")
print(funnel)

cat(sprintf("\nTier 1 (方案严格):   %d\n", sum(admet_df$Tier1, na.rm = TRUE)))
cat(sprintf("Tier 2 (底线+放宽):  %d  <- 论文最终候选\n", sum(admet_df$Tier2, na.rm = TRUE)))
cat(sprintf("Tier 3 (核心安全):   %d\n\n", sum(admet_df$Tier3, na.rm = TRUE)))

cat("OK 第十四部分 全部完成\n")
cat("   下一步: 基于 Tier 2 化合物进入第十五部分 (分子对接)\n")

sink(file.path(res.path, "session_info_admet.txt"))
sessionInfo()
sink()
