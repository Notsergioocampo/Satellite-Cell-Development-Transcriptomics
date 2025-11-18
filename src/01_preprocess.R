#!/usr/bin/env Rscript
# Satellite Cell Development Transcriptomics Pipeline
# Script 01: Data Preprocessing and Quality Control
# MIT-Level Research Pipeline

suppressPackageStartupMessages({
  library(limma)
  library(affy)
  library(preprocessCore)
  library(Biobase)
  library(yaml)
  library(logging)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(gridExtra)
})

# -----------------------------
# Helper functions (define FIRST)
# -----------------------------

perform_qc_metrics <- function(expression_matrix, sample_info, stage) {
  list(
    stage              = stage,
    n_samples          = ncol(expression_matrix),
    n_features         = nrow(expression_matrix),
    missing_values     = sum(is.na(expression_matrix)),
    mean_expression    = mean(expression_matrix, na.rm = TRUE),
    median_expression  = median(expression_matrix, na.rm = TRUE),
    expression_range   = range(expression_matrix, na.rm = TRUE),
    sample_correlations = cor(expression_matrix, method = "spearman"),
    variance_distribution = apply(expression_matrix, 1, var)
  )
}

map_probes_to_genes <- function(feature_data, expression_matrix) {
  # Extract gene symbols from feature data
  if ("Gene.symbol" %in% colnames(feature_data)) {
    gene_symbols <- feature_data$Gene.symbol
  } else if ("gene_assignment" %in% colnames(feature_data)) {
    gene_symbols <- sapply(feature_data$gene_assignment, function(x) {
      # best effort: take text before first colon
      sub(":.*", "", x)
    })
  } else {
    # Use probe IDs as gene names if no mapping available
    gene_symbols <- rownames(feature_data)
  }

  valid_mapping <- !is.na(gene_symbols) & gene_symbols != ""

  data.frame(
    probe_id    = rownames(feature_data)[valid_mapping],
    gene_symbol = gene_symbols[valid_mapping],
    stringsAsFactors = FALSE
  )
}

aggregate_probes <- function(expression_matrix, gene_mapping) {
  # Match probes in expression matrix with mapping
  common_probes <- intersect(rownames(expression_matrix), gene_mapping$probe_id)

  if (length(common_probes) == 0) {
    warning("No matching probes found between expression matrix and mapping; returning original matrix.")
    return(expression_matrix)
  }

  # Subset expression matrix and mapping, **preserving order**
  expr_subset    <- expression_matrix[common_probes, , drop = FALSE]
  mapping_subset <- gene_mapping[match(common_probes, gene_mapping$probe_id), , drop = FALSE]

  # Safety check
  stopifnot(all(rownames(expr_subset) == mapping_subset$probe_id))

  # Replace probe IDs with gene symbols
  rownames(expr_subset) <- mapping_subset$gene_symbol

  # Aggregate multiple probes per gene using median
  gene_list <- tapply(seq_len(nrow(expr_subset)),
                      rownames(expr_subset),
                      function(i) apply(expr_subset[i, , drop = FALSE], 2, median))

  gene_expression <- do.call(rbind, gene_list)
  return(gene_expression)
}

perform_pca_analysis <- function(expression_matrix, sample_info) {
  pca_result <- prcomp(t(expression_matrix), center = TRUE, scale. = TRUE)

  pca_df <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    PC3 = if (ncol(pca_result$x) >= 3) pca_result$x[, 3] else NA_real_,
    stage     = sample_info$stage,
    replicate = sample_info$replicate
  )

  var_imp <- summary(pca_result)$importance[2, ]
  variance_explained <- var_imp[seq_len(min(3, length(var_imp)))]

  list(
    pca_result        = pca_result,
    pca_df            = pca_df,
    variance_explained = variance_explained
  )
}

perform_hierarchical_clustering <- function(expression_matrix, sample_info) {
  distance_matrix <- dist(t(expression_matrix), method = "euclidean")
  hc_result       <- hclust(distance_matrix, method = "ward.D2")
  dendrogram      <- as.dendrogram(hc_result)

  list(
    hc_result       = hc_result,
    dendrogram      = dendrogram,
    distance_matrix = distance_matrix
  )
}

detect_outliers <- function(expression_matrix, sample_info) {
  # Mahalanobis on samples (columns)
  center     <- rowMeans(expression_matrix)
  cov_matrix <- cov(t(expression_matrix))

  mahal_dist <- mahalanobis(t(expression_matrix), center, cov_matrix)
  df         <- ncol(expression_matrix)
  threshold  <- qchisq(0.95, df = df)

  list(
    mahal_distances = mahal_dist,
    outliers        = which(mahal_dist > threshold),
    threshold       = threshold
  )
}

create_qc_plots <- function(qc_before, qc_after, pca_results, clustering_results, outlier_results) {
  # PCA plot
  ve <- pca_results$variance_explained
  pc1_lab <- if (length(ve) >= 1) paste0("PC1 (", round(ve[1] * 100, 1), "%)") else "PC1"
  pc2_lab <- if (length(ve) >= 2) paste0("PC2 (", round(ve[2] * 100, 1), "%)") else "PC2"

  pca_plot <- ggplot(pca_results$pca_df, aes(x = PC1, y = PC2, color = stage)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_text(aes(label = replicate), vjust = -1, hjust = 0.5, size = 3) +
    labs(
      title = "PCA Analysis of Satellite Cell Development",
      x     = pc1_lab,
      y     = pc2_lab
    ) +
    theme_minimal() +
    scale_color_brewer(palette = "Set2")

  # Variance distribution plot
  variance_plot <- ggplot(
    data.frame(variance = log10(qc_after$variance_distribution + 1e-8)),
    aes(x = variance)
  ) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    labs(
      title = "Gene Variance Distribution",
      x     = "Log10(Variance)",
      y     = "Frequency"
    ) +
    theme_minimal()

  # Correlation heatmap object (not combined in ggsave, but available)
  cor_matrix <- qc_after$sample_correlations
  cor_heatmap <- Heatmap(
    cor_matrix,
    name  = "Correlation",
    col   = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
    show_row_names    = TRUE,
    show_column_names = TRUE,
    row_names_rot     = 45,
    column_names_rot  = 45
  )

  combined <- gridExtra::grid.arrange(pca_plot, variance_plot, ncol = 2)

  list(
    pca_plot           = pca_plot,
    correlation_heatmap = cor_heatmap,
    variance_plot      = variance_plot,
    combined_plot      = combined
  )
}

generate_preprocessing_report <- function(qc_before, qc_after, final_expression, sample_info) {
  list(
    timestamp            = as.character(Sys.time()),
    initial_samples      = qc_before$n_samples,
    initial_features     = qc_before$n_features,
    final_samples        = qc_after$n_samples,
    final_features       = qc_after$n_features,
    samples_removed      = qc_before$n_samples - qc_after$n_samples,
    features_removed     = qc_before$n_features - qc_after$n_features,
    missing_values_before = qc_before$missing_values,
    missing_values_after  = qc_after$missing_values,
    developmental_stages = unique(sample_info$stage),
    replicates_per_stage = as.list(table(sample_info$stage)),
    expression_range     = range(final_expression, na.rm = TRUE)
  )
}

# -----------------------------
# Logging configuration
# -----------------------------

logfile <- file.path("results", "logs", "01_preprocess.log")
if (!dir.exists(dirname(logfile))) dir.create(dirname(logfile), recursive = TRUE)

basicConfig()  # console
addHandler(writeToFile, file = logfile, level = "INFO")

loginfo("Starting satellite cell transcriptomics preprocessing pipeline")

# -----------------------------
# Load configuration
# -----------------------------

if (!file.exists("config.yaml")) {
  logerror("config.yaml not found in project root.")
  stop("config.yaml not found. Run from project root.")
}

config <- yaml::read_yaml("config.yaml")
geo_id <- config$dataset$geo_id

# -----------------------------
# Load raw data
# -----------------------------

loginfo("Loading raw ExpressionSet from data/raw ...")
raw_file <- file.path("data", "raw", paste0(geo_id, "_raw.RData"))
if (!file.exists(raw_file)) {
  stop("Raw data file not found: ", raw_file, "\nRun 00_download.R first.")
}
load(raw_file)  # expects object 'gse'

if (!exists("gse")) {
  stop("Object 'gse' not found in raw data file. Check 00_download.R output.")
}

# Load sample information
sample_file <- file.path("data", "metadata", "sample_info.csv")
if (!file.exists(sample_file)) {
  stop("Sample info file not found: ", sample_file)
}
sample_info <- read.csv(sample_file, stringsAsFactors = FALSE)

# Ensure sample_info is aligned with expression columns
expr_matrix <- exprs(gse)
if (!"sample_id" %in% colnames(sample_info)) {
  stop("sample_info must contain a 'sample_id' column.")
}
idx <- match(colnames(expr_matrix), sample_info$sample_id)
if (any(is.na(idx))) {
  logwarn("Some samples in ExpressionSet are missing in sample_info; they will be dropped from QC annotations.")
}
sample_info <- sample_info[idx, ]
rownames(sample_info) <- sample_info$sample_id

# -----------------------------
# Step 1: Background & normalization
# -----------------------------

loginfo("Performing background correction and normalization...")

# NOTE: GEO matrices are usually preprocessed already.
# To avoid limma::backgroundCorrect errors on plain matrices, we skip explicit bg correction.
loginfo("Skipping explicit limma::backgroundCorrect (GEO matrix is usually background-corrected).")
bg_corrected <- expr_matrix

loginfo("Applying quantile normalization...")
normalized <- preprocessCore::normalize.quantiles(bg_corrected)

loginfo("Applying log2 transformation...")
log_transformed <- log2(normalized + 1)

rownames(log_transformed) <- rownames(expr_matrix)
colnames(log_transformed) <- colnames(expr_matrix)

# QC before preprocessing (using original expression)
qc_before <- perform_qc_metrics(expr_matrix, sample_info, "before")

# -----------------------------
# Step 2: Probe-to-Gene Mapping
# -----------------------------

loginfo("Performing probe-to-gene mapping...")

feature_file <- file.path("data", "metadata", paste0(geo_id, "_features.csv"))
if (!file.exists(feature_file)) {
  stop("Feature annotation file not found: ", feature_file)
}
feature_data <- read.csv(feature_file, stringsAsFactors = FALSE, row.names = 1)

gene_mapping   <- map_probes_to_genes(feature_data, log_transformed)
gene_expression <- aggregate_probes(log_transformed, gene_mapping)

# -----------------------------
# Step 3: Quality Control
# -----------------------------

loginfo("Performing comprehensive quality control...")

qc_after        <- perform_qc_metrics(gene_expression, sample_info, "after")
pca_results     <- perform_pca_analysis(gene_expression, sample_info)
clustering_results <- perform_hierarchical_clustering(gene_expression, sample_info)
outlier_results <- detect_outliers(gene_expression, sample_info)

# -----------------------------
# Step 4: Filter low-quality features
# -----------------------------

loginfo("Filtering low-quality features by variance and mean expression...")

variance_filter <- apply(gene_expression, 1, var) > config$analysis$variance_threshold
filtered_expression <- gene_expression[variance_filter, , drop = FALSE]

expression_filter <- rowMeans(filtered_expression) > config$analysis$min_counts
final_expression  <- filtered_expression[expression_filter, , drop = FALSE]

loginfo(paste("Filtered from", nrow(gene_expression), "to", nrow(final_expression), "genes."))

# -----------------------------
# Step 5: QC plots & save data
# -----------------------------

loginfo("Generating quality control plots...")

qc_plots <- create_qc_plots(qc_before, qc_after, pca_results, clustering_results, outlier_results)

if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)
qc_plot_file <- file.path("figures", "quality_control.pdf")
ggsave(qc_plot_file, qc_plots$combined_plot,
       width = 16, height = 12, dpi = config$visualization$figure_dpi)
loginfo(paste("QC plots saved to:", qc_plot_file))

# Aggregate QC results into a single object for saving
qc_results <- list(
  before      = qc_before,
  after       = qc_after,
  pca         = pca_results,
  clustering  = clustering_results,
  outliers    = outlier_results
)

processed_file <- file.path("data", "processed", "preprocessed_expression.RData")
save(final_expression, sample_info, gene_mapping, qc_results, file = processed_file)
loginfo(paste("Preprocessed data saved to:", processed_file))

# -----------------------------
# Optional: preprocessing report
# -----------------------------

loginfo("Attempting to generate preprocessing report (if template exists)...")
preprocessing_report <- generate_preprocessing_report(qc_before, qc_after, final_expression, sample_info)

report_tpl  <- "reports/preprocessing_report.Rmd"
report_file <- file.path("results", "preprocessing_report.html")

if (file.exists(report_tpl)) {
  dir.create("results", showWarnings = FALSE, recursive = TRUE)
  tryCatch({
    rmarkdown::render(
      input       = report_tpl,
      output_file = report_file,
      params      = list(report_data = preprocessing_report),
      quiet       = TRUE
    )
    loginfo(paste("Preprocessing report generated:", report_file))
  }, error = function(e) {
    logwarn(paste("Could not render preprocessing report:", conditionMessage(e)))
  })
} else {
  logwarn(paste("Report template not found at", report_tpl, "- skipping report generation."))
}

loginfo("Preprocessing pipeline completed successfully")

