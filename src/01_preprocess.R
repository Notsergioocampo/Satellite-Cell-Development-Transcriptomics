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
})

# Configure logging
logfile <- file.path("results", "logs", "01_preprocess.log")
if (!dir.exists(dirname(logfile))) dir.create(dirname(logfile), recursive = TRUE)
basicConfig(level = "INFO", file = logfile)

loginfo("Starting satellite cell transcriptomics preprocessing pipeline")

# Load configuration
config <- yaml::read_yaml("config.yaml")
geo_id <- config$dataset$geo_id

# Load raw data
loginfo("Loading raw data...")
raw_file <- file.path("data/raw", paste0(geo_id, "_raw.RData"))
if (!file.exists(raw_file)) {
  stop("Raw data file not found. Run 00_download.R first.")
}
load(raw_file)

# Load sample information
sample_file <- file.path("data/metadata", "sample_info.csv")
sample_info <- read.csv(sample_file, stringsAsFactors = FALSE)

# Step 1: Background Correction and Normalization
loginfo("Performing background correction and normalization...")

# Extract expression matrix
expr_matrix <- exprs(gse)

# Quality control before preprocessing
qc_before <- perform_qc_metrics(expr_matrix, sample_info, "before")

# Background correction using normexp method
loginfo("Applying normexp background correction...")
bg_corrected <- backgroundCorrect(expr_matrix, method = "normexp")

# Quantile normalization
loginfo("Applying quantile normalization...")
normalized <- normalize.quantiles(bg_corrected)

# Log2 transformation
loginfo("Applying log2 transformation...")
log_transformed <- log2(normalized + 1)

# Step 2: Probe-to-Gene Mapping
loginfo("Performing probe-to-gene mapping...")

# Load feature annotation
feature_file <- file.path("data/metadata", paste0(geo_id, "_features.csv"))
feature_data <- read.csv(feature_file, stringsAsFactors = FALSE)

# Map probes to gene symbols
gene_mapping <- map_probes_to_genes(feature_data, expr_matrix)

# Aggregate multiple probes per gene using median
gene_expression <- aggregate_probes(log_transformed, gene_mapping)

# Step 3: Quality Control Analysis
loginfo("Performing comprehensive quality control...")

# Quality control after preprocessing
qc_after <- perform_qc_metrics(gene_expression, sample_info, "after")

# PCA analysis
pca_results <- perform_pca_analysis(gene_expression, sample_info)

# Hierarchical clustering
clustering_results <- perform_hierarchical_clustering(gene_expression, sample_info)

# Outlier detection
outlier_results <- detect_outliers(gene_expression, sample_info)

# Step 4: Filter Low-Quality Features
loginfo("Filtering low-quality features...")

# Remove features with low variance
variance_filter <- apply(gene_expression, 1, var) > config$analysis$variance_threshold
filtered_expression <- gene_expression[variance_filter, ]

# Remove features with low expression
expression_filter <- rowMeans(filtered_expression) > config$analysis$min_counts
final_expression <- filtered_expression[expression_filter, ]

# Step 5: Generate Quality Control Reports
loginfo("Generating quality control reports...")

# Create QC plots
qc_plots <- create_qc_plots(qc_before, qc_after, pca_results, clustering_results, outlier_results)

# Save QC plots
qc_plot_file <- file.path("figures", "quality_control.pdf")
if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)
ggsave(qc_plot_file, qc_plots$combined_plot, width = 16, height = 12, dpi = config$visualization$figure_dpi)

# Save processed data
processed_file <- file.path("data/processed", "preprocessed_expression.RData")
save(final_expression, sample_info, gene_mapping, qc_results, file = processed_file)
loginfo(paste("Preprocessed data saved to:", processed_file))

# Generate preprocessing report
preprocessing_report <- generate_preprocessing_report(qc_before, qc_after, final_expression, sample_info)
report_file <- file.path("results", "preprocessing_report.html")
rmarkdown::render(
  input = system.file("rmarkdown", "templates", "preprocessing_report", package = "rmarkdown"),
  output_file = report_file,
  params = list(report_data = preprocessing_report)
)

loginfo("Preprocessing pipeline completed successfully")

# Helper functions
perform_qc_metrics <- function(expression_matrix, sample_info, stage) {
  metrics <- list(
    stage = stage,
    n_samples = ncol(expression_matrix),
    n_features = nrow(expression_matrix),
    missing_values = sum(is.na(expression_matrix)),
    mean_expression = mean(expression_matrix, na.rm = TRUE),
    median_expression = median(expression_matrix, na.rm = TRUE),
    expression_range = range(expression_matrix, na.rm = TRUE),
    sample_correlations = cor(expression_matrix, method = "spearman"),
    variance_distribution = apply(expression_matrix, 1, var)
  )
  return(metrics)
}

map_probes_to_genes <- function(feature_data, expression_matrix) {
  # Extract gene symbols from feature data
  if ("Gene.symbol" %in% colnames(feature_data)) {
    gene_symbols <- feature_data$Gene.symbol
  } else if ("gene_assignment" %in% colnames(feature_data)) {
    gene_symbols <- sapply(feature_data$gene_assignment, function(x) {
      str_extract(x, "([^:]*):.*", 1)
    })
  } else {
    # Use probe IDs as gene names if no mapping available
    gene_symbols <- rownames(feature_data)
  }
  
  # Remove NA gene symbols
  valid_mapping <- !is.na(gene_symbols) & gene_symbols != ""
  
  mapping <- data.frame(
    probe_id = rownames(feature_data)[valid_mapping],
    gene_symbol = gene_symbols[valid_mapping],
    stringsAsFactors = FALSE
  )
  
  return(mapping)
}

aggregate_probes <- function(expression_matrix, gene_mapping) {
  # Match probes in expression matrix with mapping
  common_probes <- intersect(rownames(expression_matrix), gene_mapping$probe_id)
  
  if (length(common_probes) == 0) {
    warning("No matching probes found between expression matrix and mapping")
    return(expression_matrix)
  }
  
  # Subset expression matrix and mapping
  expr_subset <- expression_matrix[common_probes, ]
  mapping_subset <- gene_mapping[gene_mapping$probe_id %in% common_probes, ]
  
  # Create a matrix with gene symbols as rownames
  rownames(expr_subset) <- mapping_subset$gene_symbol
  
  # Aggregate multiple probes per gene using median
  gene_expression <- tapply(seq_len(nrow(expr_subset)), 
                           rownames(expr_subset), 
                           function(i) apply(expr_subset[i, , drop = FALSE], 2, median))
  
  # Convert to matrix
  gene_expression <- do.call(rbind, gene_expression)
  
  return(gene_expression)
}

perform_pca_analysis <- function(expression_matrix, sample_info) {
  # Perform PCA
  pca_result <- prcomp(t(expression_matrix), center = TRUE, scale. = TRUE)
  
  # Create PCA data frame
  pca_df <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    PC3 = pca_result$x[, 3],
    stage = sample_info$stage,
    replicate = sample_info$replicate
  )
  
  # Calculate variance explained
  variance_explained <- summary(pca_result)$importance[2, 1:3]
  
  return(list(
    pca_result = pca_result,
    pca_df = pca_df,
    variance_explained = variance_explained
  ))
}

perform_hierarchical_clustering <- function(expression_matrix, sample_info) {
  # Calculate distance matrix
  distance_matrix <- dist(t(expression_matrix), method = "euclidean")
  
  # Perform hierarchical clustering
  hc_result <- hclust(distance_matrix, method = "ward.D2")
  
  # Create clustering dendrogram
  dendrogram <- as.dendrogram(hc_result)
  
  return(list(
    hc_result = hc_result,
    dendrogram = dendrogram,
    distance_matrix = distance_matrix
  ))
}

detect_outliers <- function(expression_matrix, sample_info) {
  # Use Mahalanobis distance for outlier detection
  center <- colMeans(expression_matrix)
  cov_matrix <- cov(expression_matrix)
  
  # Calculate Mahalanobis distances
  mahal_dist <- mahalanobis(t(expression_matrix), center, cov_matrix)
  
  # Determine outliers (using chi-square threshold)
  threshold <- qchisq(0.95, df = nrow(expression_matrix))
  outliers <- which(mahal_dist > threshold)
  
  return(list(
    mahal_distances = mahal_dist,
    outliers = outliers,
    threshold = threshold
  ))
}

create_qc_plots <- function(qc_before, qc_after, pca_results, clustering_results, outlier_results) {
  # PCA plot
  pca_plot <- ggplot(pca_results$pca_df, aes(x = PC1, y = PC2, color = stage)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_text(aes(label = replicate), vjust = -1, hjust = 0.5, size = 3) +
    labs(title = "PCA Analysis of Satellite Cell Development",
         x = paste0("PC1 (", round(pca_results$variance_explained[1] * 100, 1), "%)"),
         y = paste0("PC2 (", round(pca_results$variance_explained[2] * 100, 1), "%)")) +
    theme_minimal() +
    scale_color_brewer(palette = "Set2")
  
  # Sample correlation heatmap
  cor_matrix <- qc_after$sample_correlations
  cor_heatmap <- Heatmap(cor_matrix,
                        name = "Correlation",
                        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
                        show_row_names = TRUE,
                        show_column_names = TRUE,
                        row_names_rot = 45,
                        column_names_rot = 45)
  
  # Variance distribution plot
  variance_plot <- ggplot(data.frame(variance = log10(qc_after$variance_distribution)), 
                         aes(x = variance)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    labs(title = "Gene Variance Distribution",
         x = "Log10(Variance)",
         y = "Frequency") +
    theme_minimal()
  
  return(list(
    pca_plot = pca_plot,
    correlation_heatmap = cor_heatmap,
    variance_plot = variance_plot,
    combined_plot = gridExtra::grid.arrange(pca_plot, variance_plot, ncol = 2)
  ))
}

generate_preprocessing_report <- function(qc_before, qc_after, final_expression, sample_info) {
  report <- list(
    timestamp = Sys.time(),
    initial_samples = qc_before$n_samples,
    initial_features = qc_before$n_features,
    final_samples = qc_after$n_samples,
    final_features = qc_after$n_features,
    samples_removed = qc_before$n_samples - qc_after$n_samples,
    features_removed = qc_before$n_features - qc_after$n_features,
    missing_values_before = qc_before$missing_values,
    missing_values_after = qc_after$missing_values,
    developmental_stages = unique(sample_info$stage),
    replicates_per_stage = table(sample_info$stage),
    expression_range = range(final_expression, na.rm = TRUE)
  )
  
  return(report)
}
