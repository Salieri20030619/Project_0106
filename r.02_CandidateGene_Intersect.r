
cat("\n========================================\n")
cat("02. Candidate Gene Identification (Intersection of DEGs)\n")
cat("========================================\n\n")

rm(list = ls())

options(bitmapType = "cairo")

setwd("")

input_dir <- "./01_DEGs/Results"

output_dir <- "./02_CandidateGene_Intersect"
res.path <- file.path(output_dir, "Results")
fig.path <- file.path(output_dir, "Figures")

for (path in c(output_dir, res.path, fig.path)) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

cat("Input directory (DEGs):", input_dir, "\n")
cat("Output directory (Candidate genes):", output_dir, "\n")
cat("Results directory:", res.path, "\n")
cat("Figures directory:", fig.path, "\n\n")

cat("\n========== Load R packages ==========\n")

if (!require("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram", repos = "http://cran.us.r-project.org")
  library(VennDiagram)
}
cat("Loaded: VennDiagram\n")

packages <- c("ggplot2", "gridExtra", "grid")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
    library(pkg, character.only = TRUE)
  } else {
    cat("Loaded:", pkg, "\n")
  }
}

cat("\n========== 1. Read all DEGs ==========\n")

deg1_file <- file.path(input_dir, "DEG1_GSE16561_significant.csv")
deg1 <- read.csv(deg1_file, stringsAsFactors = FALSE)
deg1_up <- deg1$Gene[deg1$logFC > 0]
deg1_down <- deg1$Gene[deg1$logFC < 0]
deg1_all <- deg1$Gene
cat("DEG1 (Blood): Up =", length(deg1_up), ", Down =", length(deg1_down), ", Total =", length(deg1_all), "\n")

deg2_file <- file.path(input_dir, "DEG2_Brain_C57_significant.csv")
deg2 <- read.csv(deg2_file, stringsAsFactors = FALSE)
deg2_up <- deg2$Gene[deg2$log2FoldChange > 0]
deg2_down <- deg2$Gene[deg2$log2FoldChange < 0]
deg2_all <- deg2$Gene
cat("DEG2 (Brain_C57): Up =", length(deg2_up), ", Down =", length(deg2_down), ", Total =", length(deg2_all), "\n")

deg3_file <- file.path(input_dir, "DEG3_Brain_IL21_significant.csv")
deg3 <- read.csv(deg3_file, stringsAsFactors = FALSE)
deg3_up <- deg3$Gene[deg3$log2FoldChange > 0]
deg3_down <- deg3$Gene[deg3$log2FoldChange < 0]
deg3_all <- deg3$Gene
cat("DEG3 (Brain_IL21): Up =", length(deg3_up), ", Down =", length(deg3_down), ", Total =", length(deg3_all), "\n")

deg4_file <- file.path(input_dir, "DEG4_Spleen_C57_significant.csv")
deg4 <- read.csv(deg4_file, stringsAsFactors = FALSE)
deg4_up <- deg4$Gene[deg4$log2FoldChange > 0]
deg4_down <- deg4$Gene[deg4$log2FoldChange < 0]
deg4_all <- deg4$Gene
cat("DEG4 (Spleen_C57): Up =", length(deg4_up), ", Down =", length(deg4_down), ", Total =", length(deg4_all), "\n")

deg5_file <- file.path(input_dir, "DEG5_Spleen_IL21_significant.csv")
deg5 <- read.csv(deg5_file, stringsAsFactors = FALSE)
deg5_up <- deg5$Gene[deg5$log2FoldChange > 0]
deg5_down <- deg5$Gene[deg5$log2FoldChange < 0]
deg5_all <- deg5$Gene
cat("DEG5 (Spleen_IL21): Up =", length(deg5_up), ", Down =", length(deg5_down), ", Total =", length(deg5_all), "\n")

intersect_multiple <- function(...) {
  args <- list(...)
  if (length(args) == 0) return(character(0))
  result <- args[[1]]
  for (i in 2:length(args)) {
    result <- intersect(result, args[[i]])
  }
  return(result)
}

cat("\n========== 2. Compute Combination 124 intersection ==========\n")
cat("Combination 124: DEG1 (Blood) + DEG2 (Brain_C57) + DEG4 (Spleen_C57)\n")

up_124 <- intersect_multiple(deg1_up, deg2_up, deg4_up)
down_124 <- intersect_multiple(deg1_down, deg2_down, deg4_down)

cat("   Up-regulated genes:", length(up_124), "\n")
cat("   Down-regulated genes:", length(down_124), "\n")
if (length(up_124) > 0) cat("   Up-regulated:", paste(up_124, collapse = ", "), "\n")
if (length(down_124) > 0) cat("   Down-regulated:", paste(down_124, collapse = ", "), "\n")

candidate_up <- up_124
candidate_down <- down_124
source_type <- "combination_124"

candidate_genes <- c(candidate_up, candidate_down)

cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("Final candidate gene results:\n")
cat(paste(rep("=", 50), collapse = ""), "\n")
cat("Up-regulated genes:", length(candidate_up), "\n")
cat("Down-regulated genes:", length(candidate_down), "\n")
cat("Total candidate genes:", length(candidate_genes), "\n")

if (length(candidate_up) > 0) {
  cat("\nUp-regulated gene list:\n")
  print(candidate_up)
}
if (length(candidate_down) > 0) {
  cat("\nDown-regulated gene list:\n")
  print(candidate_down)
}

cat("\n========== 3. Draw Combination 124 Venn diagram ==========\n")

combo124_colors <- c("#E64B35", "#4DBBD5", "#00A087")
combo124_names  <- c("DEG1_Blood", "DEG2_Brain_C57", "DEG4_Spleen_C57")

venn124_up <- venn.diagram(
  x = list("DEG1_Blood"      = deg1_up,
           "DEG2_Brain_C57"  = deg2_up,
           "DEG4_Spleen_C57" = deg4_up),
  category.names = combo124_names,
  filename = NULL, output = TRUE,
  scaled = FALSE, euler.d = FALSE,
  col = combo124_colors, fill = combo124_colors, alpha = 0.5, lwd = 2,
  main = paste0("Up-regulated (n=", length(up_124), ")"),
  main.cex = 1.8, main.fontface = "bold",
  cat.cex = 1.5, cat.fontface = "bold", cat.col = combo124_colors,
  cat.dist = c(0.07, 0.07, 0.05),
  cat.pos = c(-20, 20, 180),
  cex = 1.6, fontface = "bold", margin = 0.12,
  disable.logging = TRUE
)

venn124_down <- venn.diagram(
  x = list("DEG1_Blood"      = deg1_down,
           "DEG2_Brain_C57"  = deg2_down,
           "DEG4_Spleen_C57" = deg4_down),
  category.names = combo124_names,
  filename = NULL, output = TRUE,
  scaled = FALSE, euler.d = FALSE,
  col = combo124_colors, fill = combo124_colors, alpha = 0.5, lwd = 2,
  main = paste0("Down-regulated (n=", length(down_124), ")"),
  main.cex = 1.8, main.fontface = "bold",
  cat.cex = 1.5, cat.fontface = "bold", cat.col = combo124_colors,
  cat.dist = c(0.07, 0.07, 0.05),
  cat.pos = c(-20, 20, 180),
  cex = 1.6, fontface = "bold", margin = 0.12,
  disable.logging = TRUE
)

pdf(file.path(fig.path, "Venn_Diagram_Combo124.pdf"), width = 20, height = 10)
grid.newpage()

grid.text("Combination 124 (Blood + Brain_C57 + Spleen_C57)",
          x = 0.5, y = 0.96, gp = gpar(fontsize = 20, fontface = "bold"))

pushViewport(viewport(x = 0.25, y = 0.46, width = 0.5, height = 0.88))
grid.draw(venn124_up)
popViewport()

pushViewport(viewport(x = 0.75, y = 0.46, width = 0.5, height = 0.88))
grid.draw(venn124_down)
popViewport()
dev.off()

png(file.path(fig.path, "Venn_Diagram_Combo124.png"),
    width = 20, height = 10, units = "in", res = 300)
grid.newpage()
grid.text("Combination 124 (Blood + Brain_C57 + Spleen_C57)",
          x = 0.5, y = 0.96, gp = gpar(fontsize = 20, fontface = "bold"))
pushViewport(viewport(x = 0.25, y = 0.46, width = 0.5, height = 0.88))
grid.draw(venn124_up)
popViewport()
pushViewport(viewport(x = 0.75, y = 0.46, width = 0.5, height = 0.88))
grid.draw(venn124_down)
popViewport()
dev.off()

cat("Combination 124 Up/Down side-by-side Venn diagram saved\n")

cat("\n========== 4. Save results to 02_CandidateGene_Intersect ==========\n")

if (length(candidate_genes) > 0) {

  writeLines(candidate_genes, file.path(res.path, "final_candidate_genes.txt"))

  candidate_df <- data.frame(
    Gene = candidate_genes,
    Regulation = ifelse(candidate_genes %in% candidate_up, "Up", "Down"),
    stringsAsFactors = FALSE
  )
  write.csv(candidate_df, file.path(res.path, "final_candidate_genes.csv"), row.names = FALSE)

  if (length(candidate_up) > 0) {
    writeLines(candidate_up, file.path(res.path, "final_candidate_up.txt"))
  }
  if (length(candidate_down) > 0) {
    writeLines(candidate_down, file.path(res.path, "final_candidate_down.txt"))
  }

  cat("Candidate gene list saved to:", res.path, "\n")
}

cat("\n========== 5. Save combination 124 results ==========\n")

if (length(up_124) > 0) {
  writeLines(up_124, file.path(res.path, "combo_124_up_genes.txt"))
}
if (length(down_124) > 0) {
  writeLines(down_124, file.path(res.path, "combo_124_down_genes.txt"))
}

cat("Combination 124 results saved\n")

cat("\n========== 6. Results summary ==========\n")

summary_df <- data.frame(
  Dataset = c("DEG1_Blood", "DEG2_Brain_C57", "DEG3_Brain_IL21", "DEG4_Spleen_C57", "DEG5_Spleen_IL21",
              "Combo124_Up", "Combo124_Down", "Final_Candidate_Genes"),
  Total_Genes = c(length(deg1_all), length(deg2_all), length(deg3_all),
                  length(deg4_all), length(deg5_all),
                  length(up_124), length(down_124),
                  length(candidate_genes)),
  Upregulated = c(length(deg1_up), length(deg2_up), length(deg3_up),
                  length(deg4_up), length(deg5_up),
                  length(up_124), NA,
                  length(candidate_up)),
  Downregulated = c(length(deg1_down), length(deg2_down), length(deg3_down),
                    length(deg4_down), length(deg5_down),
                    NA, length(down_124),
                    length(candidate_down))
)

write.csv(summary_df, file.path(res.path, "DEGs_summary_with_candidates.csv"), row.names = FALSE)

cat("\nDEGs and candidate genes summary:\n")
print(summary_df)

cat("\n========== 7. Candidate gene detailed information ==========\n")

candidate_detail <- data.frame(
  Gene = candidate_genes,
  Regulation = ifelse(candidate_genes %in% candidate_up, "Up", "Down"),
  In_DEG1 = candidate_genes %in% deg1_all,
  In_DEG2 = candidate_genes %in% deg2_all,
  In_DEG3 = candidate_genes %in% deg3_all,
  In_DEG4 = candidate_genes %in% deg4_all,
  In_DEG5 = candidate_genes %in% deg5_all,
  stringsAsFactors = FALSE
)

write.csv(candidate_detail, file.path(res.path, "candidate_genes_detail.csv"), row.names = FALSE)

cat("Candidate gene details saved\n")

cat("\n========== 8. Save RData file ==========\n")

save(list = c("candidate_genes", "candidate_up", "candidate_down",
              "candidate_detail", "source_type", "summary_df",
              "up_124", "down_124",
              "deg1_all", "deg2_all", "deg3_all", "deg4_all", "deg5_all",
              "deg1_up", "deg1_down", "deg2_up", "deg2_down",
              "deg3_up", "deg3_down", "deg4_up", "deg4_down",
              "deg5_up", "deg5_down"),
     file = file.path(res.path, "candidate_genes_results.RData"))

cat("RData file saved\n")

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("        Candidate gene identification complete!\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")
cat("Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\nInput directory (DEGs):", input_dir, "\n")
cat("Output directory (Candidate genes):", output_dir, "\n")
cat("\nCandidate gene source:\n")
cat("   - Combination 124 (DEG1 + DEG2 + DEG4): Blood + Brain_C57 + Spleen_C57\n")
cat("\nFinal candidate gene statistics:\n")
cat("   - Up-regulated genes:", length(candidate_up), "\n")
cat("   - Down-regulated genes:", length(candidate_down), "\n")
cat("   - Total:", length(candidate_genes), "\n")
cat("\nGenerated files:\n")
cat("   - Results/final_candidate_genes.txt\n")
cat("   - Results/final_candidate_genes.csv\n")
cat("   - Results/candidate_genes_detail.csv\n")
cat("   - Results/DEGs_summary_with_candidates.csv\n")
cat("   - Results/candidate_genes_results.RData\n")
cat("   - Results/combo_124_up_genes.txt\n")
cat("   - Results/combo_124_down_genes.txt\n")
cat("   - Figures/Venn_Diagram_Combo124.pdf/png\n")
cat("\nAll processing complete!\n")
