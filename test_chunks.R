# Test script to verify the new chunks work correctly
library(ggplot2)
library(dplyr)

# Load the data
de_p12_p1 <- readRDS("results/de_results_P12_vs_P1.rds")
hallmark_p28_p1 <- readRDS("results/hallmark_enrichment_P28_vs_P1.rds")
cfg_analysis <- list(padj_threshold = 0.05, log2fc_threshold = 1)

cat("=== Testing Volcano Plot Chunk ===\n")
if (is.null(de_p12_p1)) {
  cat("DE results are NULL\n")
} else if (!all(c("log2FoldChange", "padj") %in% names(de_p12_p1))) {
  cat("DE results missing required columns\n")
} else {
  padj_cut <- cfg_analysis$padj_threshold
  lfc_cut <- cfg_analysis$log2fc_threshold
  
  plot_df <- de_p12_p1 |>
    filter(!is.na(log2FoldChange), !is.na(padj), padj > 0) |>
    mutate(
      sig = case_when(
        padj <= padj_cut & log2FoldChange >= lfc_cut ~ "Up",
        padj <= padj_cut & log2FoldChange <= -lfc_cut ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  cat("Volcano plot data prepared successfully:\n")
  cat("- Total points:", nrow(plot_df), "\n")
  cat("- Significance levels:", paste(unique(plot_df$sig), collapse = ", "), "\n")
  cat("- Upregulated genes:", sum(plot_df$sig == "Up", na.rm = TRUE), "\n")
  cat("- Downregulated genes:", sum(plot_df$sig == "Down", na.rm = TRUE), "\n")
  cat("- Non-significant genes:", sum(plot_df$sig == "NS", na.rm = TRUE), "\n\n")
}

cat("=== Testing Hallmark Barplot Chunk ===\n")
if (is.null(hallmark_p28_p1)) {
  cat("Hallmark results are NULL - will show placeholder\n")
} else {
  df <- hallmark_p28_p1
  
  if ("p.adjust" %in% names(df) && !"padj" %in% names(df)) {
    df$padj <- df$p.adjust
  }
  
  df_clean <- df |>
    filter(!is.na(padj), padj > 0, padj <= cfg_analysis$padj_threshold) |>
    arrange(padj) |>
    slice_head(n = 15) |>
    mutate(Pathway = factor(Description, levels = rev(Description)))
  
  if (nrow(df_clean) == 0) {
    cat("No pathways pass significance threshold - will show placeholder\n")
  } else {
    cat("Hallmark barplot data prepared successfully:\n")
