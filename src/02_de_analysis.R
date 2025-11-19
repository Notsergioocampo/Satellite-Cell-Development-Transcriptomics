#!/usr/bin/env Rscript
# Satellite Cell Development Transcriptomics Pipeline
# Script 02: Differential Expression Analysis using limma
# MIT-Level Research Pipeline

suppressPackageStartupMessages({
  library(Biobase)
  library(limma)
  library(readr)
  library(dplyr)
})

# ============================
# Load data
# ============================

# Load raw ExpressionSet
raw_data_file <- "data/raw/GSE65927_raw.RData"
if (!file.exists(raw_data_file)) {
  stop("Raw data file not found: ", raw_data_file)
}

load(raw_data_file)

# Check that gse object exists
if (!exists("gse")) {
  stop("ExpressionSet 'gse' not found in loaded data")
}

message("Loaded ExpressionSet with ", nrow(exprs(gse)), " genes and ", ncol(exprs(gse)), " samples")

# Extract expression matrix
expr_mat <- Biobase::exprs(gse)
message("Expression matrix dimensions: ", nrow(expr_mat), " x ", ncol(expr_mat))

# Load sample metadata from phenotype file (contains stage info)
phenotype_file <- "data/metadata/GSE65927_phenotype.csv"
if (!file.exists(phenotype_file)) {
  stop("Phenotype file not found: ", phenotype_file)
}

phenotype_data <- read_csv(phenotype_file, show_col_types = FALSE)

# Extract stage from title column - be more specific to avoid confusion
extract_stage <- function(title) {
  if (grepl("P28", title)) return("P28")
  if (grepl("P12", title)) return("P12")
  if (grepl("P1", title)) return("P1")
  return("unknown")
}

# Create sample_info with proper stage information
sample_info <- data.frame(
  sample_id = phenotype_data$geo_accession,
  stage = sapply(phenotype_data$title, extract_stage),
  replicate = phenotype_data$title,
  stringsAsFactors = FALSE
)

message("Created sample info with ", nrow(sample_info), " samples")
message("Sample stages: ", paste(table(sample_info$stage), collapse = ", "))

# ============================
# Validate data
# ============================

# Check required columns
required_cols <- c("sample_id", "stage")
missing_cols <- setdiff(required_cols, names(sample_info))
if (length(missing_cols) > 0) {
  stop("Missing required columns in sample_info: ", paste(missing_cols, collapse = ", "))
}

# Check that all sample IDs match
expr_samples <- colnames(expr_mat)
meta_samples <- sample_info$sample_id

missing_in_meta <- setdiff(expr_samples, meta_samples)
missing_in_expr <- setdiff(meta_samples, expr_samples)

if (length(missing_in_meta) > 0) {
  stop("Samples in expression matrix not found in metadata: ", paste(missing_in_meta, collapse = ", "))
}

if (length(missing_in_expr) > 0) {
  warning("Samples in metadata not found in expression matrix: ", paste(missing_in_expr, collapse = ", "))
  # Remove samples not in expression matrix
  sample_info <- sample_info[sample_info$sample_id %in% expr_samples, ]
}

# Reorder sample_info to match expression matrix
sample_info <- sample_info[match(expr_samples, sample_info$sample_id), ]

# Check for NA stages
if (any(is.na(sample_info$stage))) {
  stop("Some samples have missing stage information")
}

# Create stage factor with proper levels
sample_info$stage <- factor(sample_info$stage, levels = c("P1", "P12", "P28"))
message("Sample stages: ", paste(table(sample_info$stage), collapse = ", "))

# ============================
# Create design matrix
# ============================

design <- model.matrix(~ 0 + stage, data = sample_info)
colnames(design) <- c("P1", "P12", "P28")
message("Design matrix:")
print(design)

# ============================
# Fit limma model
# ============================

message("Fitting limma model...")
fit <- limma::lmFit(expr_mat, design)

# Define contrasts
cont <- limma::makeContrasts(
  P12_vs_P1  = P12 - P1,
  P28_vs_P1  = P28 - P1,
  P28_vs_P12 = P28 - P12,
  levels = design
)

message("Contrasts:")
print(cont)

# Apply contrasts and empirical Bayes
fit2 <- limma::contrasts.fit(fit, cont)
fit2 <- limma::eBayes(fit2)

# ============================
# Extract results for each contrast
# ============================

# Create results directory if it doesn't exist
if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}

contrasts_names <- c("P12_vs_P1", "P28_vs_P1", "P28_vs_P12")

for (contrast in contrasts_names) {
  message("\nProcessing contrast: ", contrast)
  
  # Get results for all genes
  res <- limma::topTable(fit2, coef = contrast, number = Inf, sort.by = "P")
  
  # Standardize column names
  res <- res %>%
    rename(
      log2FoldChange = logFC,
      padj = adj.P.Val,
      pvalue = P.Value
    )
  
  # Ensure required columns exist
  required_cols <- c("log2FoldChange", "padj", "pvalue")
  missing_cols <- setdiff(required_cols, names(res))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in results: ", paste(missing_cols, collapse = ", "))
  }
  
  message("  Genes in results: ", nrow(res))
  message("  Significant genes (FDR < 0.05): ", sum(res$padj < 0.05, na.rm = TRUE))
  message("  Upregulated genes (FDR < 0.05, log2FC > 1): ", 
          sum(res$padj < 0.05 & res$log2FoldChange > 1, na.rm = TRUE))
  message("  Downregulated genes (FDR < 0.05, log2FC < -1): ", 
          sum(res$padj < 0.05 & res$log2FoldChange < -1, na.rm = TRUE))
  
  # Save results
  output_file <- paste0("results/de_results_", contrast, ".rds")
  saveRDS(res, file = output_file)
  message("  Saved to: ", output_file)
}

message("\nDifferential expression analysis completed successfully!")
