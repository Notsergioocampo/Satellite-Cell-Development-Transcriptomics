#!/usr/bin/env Rscript
# Satellite Cell Development Transcriptomics Pipeline
# Script 02: Differential Expression Analysis
# MIT-Level Research Pipeline

suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(limma)
  library(yaml)
  library(logging)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(EnhancedVolcano)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
})

# Configure logging
logfile <- file.path("results", "logs", "02_de_analysis.log")
if (!dir.exists(dirname(logfile))) dir.create(dirname(logfile), recursive = TRUE)
basicConfig(level = "INFO", file = logfile)

loginfo("Starting differential expression analysis pipeline")

# Load configuration
config <- yaml::read_yaml("config.yaml")

# Load preprocessed data
loginfo("Loading preprocessed data...")
processed_file <- file.path("data/processed", "preprocessed_expression.RData")
if (!file.exists(processed_file)) {
  stop("Preprocessed data not found. Run 01_preprocess.R first.")
}
load(processed_file)

# Create DESeq2 dataset
loginfo("Creating DESeq2 dataset...")
dds <- create_deseq_dataset(final_expression, sample_info)

# Perform differential expression analysis
loginfo("Performing differential expression analysis...")

# Define contrasts for developmental stages
contrasts <- list(
  P12_vs_P1 = c("stage", "P12", "P1"),
  P28_vs_P12 = c("stage", "P28", "P12"),
  P28_vs_P1 = c("stage", "P28", "P1")
)

# Run differential expression for each contrast
de_results <- list()
for (contrast_name in names(contrasts)) {
  loginfo(paste("Analyzing contrast:", contrast_name))
  
  contrast <- contrasts[[contrast_name]]
  res <- results(dds, contrast = contrast, alpha = config$analysis$padj_threshold)
  
  # Filter significant results
  sig_genes <- res[which(res$padj < config$analysis$padj_threshold & 
                        abs(res$log2FoldChange) > config$analysis$log2fc_threshold), ]
  
  de_results[[contrast_name]] <- list(
    full_results = res,
    significant_genes = sig_genes,
    n_upregulated = sum(res$log2FoldChange > config$analysis$log2fc_threshold & res$padj < config$analysis$padj_threshold),
    n_downregulated = sum(res$log2FoldChange < -config$analysis$log2fc_threshold & res$padj < config$analysis$padj_threshold),
    n_total_significant = nrow(sig_genes)
  )
  
  loginfo(paste("Found", nrow(sig_genes), "significant genes for", contrast_name))
}

# Perform time-course analysis
loginfo("Performing time-course differential expression analysis...")
timecourse_results <- perform_timecourse_analysis(final_expression, sample_info)

# Identify developmental stage-specific genes
loginfo("Identifying stage-specific gene expression patterns...")
stage_specific_genes <- identify_stage_specific_genes(final_expression, sample_info, de_results)

# Create volcano plots
loginfo("Creating volcano plots...")
volcano_plots <- create_volcano_plots(de_results, config$analysis$padj_threshold, config$analysis$log2fc_threshold)

# Create heatmaps
loginfo("Creating differential expression heatmaps...")
de_heatmaps <- create_de_heatmaps(final_expression, sample_info, de_results)

# Perform gene set enrichment analysis
loginfo("Performing gene set enrichment analysis...")
gsea_results <- perform_gsea_analysis(de_results, config$analysis$padj_threshold)

# Create trajectory analysis
loginfo("Creating gene expression trajectories...")
trajectory_plots <- create_expression_trajectories(final_expression, sample_info, de_results)

# Save results
loginfo("Saving differential expression results...")

# Save DE results
de_file <- file.path("results", "de_analysis_results.RData")
save(de_results, timecourse_results, stage_specific_genes, gsea_results, file = de_file)

# Save plots
volcano_file <- file.path("figures", "volcano_plots.pdf")
ggsave(volcano_file, volcano_plots$combined_plot, width = 16, height = 12, dpi = config$visualization$figure_dpi)

heatmap_file <- file.path("figures", "de_heatmaps.pdf")
ggsave(heatmap_file, de_heatmaps$combined_plot, width = 14, height = 10, dpi = config$visualization$figure_dpi)

trajectory_file <- file.path("figures", "expression_trajectories.pdf")
ggsave(trajectory_file, trajectory_plots$combined_plot, width = 14, height = 10, dpi = config$visualization$figure_dpi)

# Generate comprehensive DE report
de_report <- generate_de_report(de_results, timecourse_results, stage_specific_genes, gsea_results)
report_file <- file.path("results", "de_analysis_report.html")
rmarkdown::render(
  input = system.file("rmarkdown", "templates", "de_report", package = "rmarkdown"),
  output_file = report_file,
  params = list(report_data = de_report)
)

loginfo("Differential expression analysis pipeline completed successfully")

# Helper functions
create_deseq_dataset <- function(expression_matrix, sample_info) {
  # Create count matrix (convert from log2 back to counts approximation)
  count_matrix <- round(2^expression_matrix - 1)
  count_matrix[count_matrix < 0] <- 0
  
  # Create DESeqDataSet
  sample_conditions <- factor(sample_info$stage, levels = c("P1", "P12", "P28"))
  colData <- DataFrame(condition = sample_conditions)
  
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = colData,
    design = ~ condition
  )
  
  # Filter low count genes
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  # Run DESeq
  dds <- DESeq(dds)
  
  return(dds)
}

perform_timecourse_analysis <- function(expression_matrix, sample_info) {
  # Create time points
  time_points <- ifelse(sample_info$stage == "P1", 1, 
                       ifelse(sample_info$stage == "P12", 12, 28))
  
  # Fit linear models for each gene
  timecourse_results <- list()
  
  for (gene in rownames(expression_matrix)) {
    gene_expr <- as.numeric(expression_matrix[gene, ])
    
    # Fit polynomial model (quadratic)
    model <- lm(gene_expr ~ poly(time_points, 2))
    
    # Calculate significance of time trend
    anova_result <- anova(model)
    
    timecourse_results[[gene]] <- list(
      gene = gene,
      linear_coef = coef(model)[2],
      quadratic_coef = coef(model)[3],
      linear_pval = anova_result$"Pr(>F)"[2],
      quadratic_pval = anova_result$"Pr(>F)"[3],
      r_squared = summary(model)$r.squared,
      trend = determine_trend(coef(model)[2], coef(model)[3])
    )
  }
  
  # Convert to data frame
  timecourse_df <- do.call(rbind, lapply(timecourse_results, function(x) {
    data.frame(
      gene = x$gene,
      linear_coef = x$linear_coef,
      quadratic_coef = x$quadratic_coef,
      linear_pval = x$linear_pval,
      quadratic_pval = x$quadratic_pval,
      r_squared = x$r_squared,
      trend = x$trend,
      stringsAsFactors = FALSE
    )
  }))
  
  # Adjust p-values
  timecourse_df$linear_padj <- p.adjust(timecourse_df$linear_pval, method = "BH")
  timecourse_df$quadratic_padj <- p.adjust(timecourse_df$quadratic_pval, method = "BH")
  
  return(timecourse_df)
}

identify_stage_specific_genes <- function(expression_matrix, sample_info, de_results) {
  stage_specific <- list()
  
  # For each developmental stage
  for (stage in unique(sample_info$stage)) {
    # Get genes highly expressed in this stage
    stage_samples <- sample_info$stage == stage
    other_samples <- sample_info$stage != stage
    
    stage_means <- rowMeans(expression_matrix[, stage_samples])
    other_means <- rowMeans(expression_matrix[, other_samples])
    
    # Calculate fold change
    fold_changes <- log2(stage_means + 1) - log2(other_means + 1)
    
    # Get significant genes from DE analysis
    significant_genes <- unique(unlist(lapply(de_results, function(x) {
      rownames(x$significant_genes)
    })))
    
    # Filter for stage-specific genes
    stage_genes <- names(fold_changes)[which(fold_changes > 1 & 
                                           names(fold_changes) %in% significant_genes)]
    
    stage_specific[[stage]] <- list(
      genes = stage_genes,
      n_genes = length(stage_genes),
      top_genes = head(stage_genes[order(fold_changes[stage_genes], decreasing = TRUE)], 20)
    )
  }
  
  return(stage_specific)
}

create_volcano_plots <- function(de_results, padj_threshold, log2fc_threshold) {
  volcano_plots <- list()
  
  for (contrast_name in names(de_results)) {
    res <- de_results[[contrast_name]]$full_results
    
    # Create volcano plot
    volcano_plot <- EnhancedVolcano(res,
                                   lab = rownames(res),
                                   x = 'log2FoldChange',
                                   y = 'padj',
                                   title = paste("Volcano Plot:", contrast_name),
                                   pCutoff = padj_threshold,
                                   FCcutoff = log2fc_threshold,
                                   pointSize = 2.0,
                                   labSize = 4.0,
                                   legendPosition = "right",
                                   legendLabSize = 12,
                                   legendIconSize = 4.0,
                                   col = c("black", "red", "green3", "blue"),
                                   colAlpha = 0.8)
    
    volcano_plots[[contrast_name]] <- volcano_plot
  }
  
  # Combine plots
  combined_plot <- gridExtra::grid.arrange(grobs = volcano_plots, ncol = 2)
  
  return(list(
    individual_plots = volcano_plots,
    combined_plot = combined_plot
  ))
}

create_de_heatmaps <- function(expression_matrix, sample_info, de_results) {
  # Get all significant genes
  all_sig_genes <- unique(unlist(lapply(de_results, function(x) rownames(x$significant_genes))))
  
  if (length(all_sig_genes) == 0) {
    warning("No significant genes found for heatmap")
    return(NULL)
  }
  
  # Subset expression matrix
  sig_expression <- expression_matrix[all_sig_genes, ]
  
  # Order samples by developmental stage
  stage_order <- order(factor(sample_info$stage, levels = c("P1", "P12", "P28")))
  sig_expression_ordered <- sig_expression[, stage_order]
  sample_info_ordered <- sample_info[stage_order, ]
  
  # Z-score normalization
  sig_expression_zscore <- t(scale(t(sig_expression_ordered)))
  
  # Create heatmap
  heatmap_plot <- pheatmap(sig_expression_zscore,
                           annotation_col = data.frame(stage = sample_info_ordered$stage),
                           show_rownames = FALSE,
                           cluster_rows = TRUE,
                           cluster_cols = FALSE,
                           main = "Differential Expression Heatmap",
                           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
  
  return(list(
    heatmap_plot = heatmap_plot,
    sig_genes = all_sig_genes,
    sig_expression = sig_expression_zscore
  ))
}

perform_gsea_analysis <- function(de_results, padj_threshold) {
  gsea_results <- list()
  
  # Load gene sets (using clusterProfiler and org.Mm.eg.db)
  tryCatch({
    library(clusterProfiler)
    library(org.Mm.eg.db)
    
    for (contrast_name in names(de_results)) {
      res <- de_results[[contrast_name]]$full_results
      
      # Prepare gene list for GSEA
      gene_list <- res$log2FoldChange
      names(gene_list) <- rownames(res)
      gene_list <- sort(gene_list, decreasing = TRUE)
      
      # Perform GO enrichment
      go_enrich <- enrichGO(gene = names(gene_list)[gene_list > 0],
                           OrgDb = org.Mm.eg.db,
                           ont = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = padj_threshold,
                           qvalueCutoff = padj_threshold)
      
      # Perform KEGG enrichment
      kegg_enrich <- enrichKEGG(gene = names(gene_list)[gene_list > 0],
                               organism = "mmu",
                               pvalueCutoff = padj_threshold,
                               qvalueCutoff = padj_threshold)
      
      gsea_results[[contrast_name]] <- list(
        go_enrichment = go_enrich,
        kegg_enrichment = kegg_enrich,
        gene_list = gene_list
      )
    }
    
  }, error = function(e) {
    logwarn(paste("GSEA analysis failed:", e$message))
    gsea_results <- list(error = e$message)
  })
  
  return(gsea_results)
}

create_expression_trajectories <- function(expression_matrix, sample_info, de_results) {
  # Get top differentially expressed genes
  top_genes <- unique(unlist(lapply(de_results, function(x) {
    head(rownames(x$significant_genes), 50)
  })))
  
  # Create time points
  time_points <- ifelse(sample_info$stage == "P1", 1, 
                       ifelse(sample_info$stage == "P12", 12, 28))
  
  # Calculate mean expression per stage
  trajectory_data <- list()
  for (gene in top_genes) {
    if (gene %in% rownames(expression_matrix)) {
      gene_expr <- expression_matrix[gene, ]
      
      stage_means <- aggregate(gene_expr, list(stage = sample_info$stage), mean)
      colnames(stage_means) <- c("stage", "expression")
      stage_means$time_point <- c(1, 12, 28)
      stage_means$gene <- gene
      
      trajectory_data[[gene]] <- stage_means
    }
  }
  
  # Combine all trajectories
  all_trajectories <- do.call(rbind, trajectory_data)
  
  # Create trajectory plot
  trajectory_plot <- ggplot(all_trajectories, aes(x = time_point, y = expression, color = gene)) +
    geom_line(size = 1.2, alpha = 0.8) +
    geom_point(size = 3) +
    labs(title = "Gene Expression Trajectories During Satellite Cell Development",
         x = "Postnatal Day",
         y = "Normalized Expression") +
    theme_minimal() +
    scale_x_continuous(breaks = c(1, 12, 28), labels = c("P1", "P12", "P28")) +
    theme(legend.position = "none")
  
  return(list(
    trajectory_plot = trajectory_plot,
    trajectory_data = all_trajectories
  ))
}

determine_trend <- function(linear_coef, quadratic_coef) {
  if (abs(quadratic_coef) > abs(linear_coef)) {
    if (quadratic_coef > 0) return("quadratic_up") else return("quadratic_down")
  } else {
    if (linear_coef > 0) return("linear_up") else return("linear_down")
  }
}

generate_de_report <- function(de_results, timecourse_results, stage_specific_genes, gsea_results) {
  report <- list(
    timestamp = Sys.time(),
    contrasts_analyzed = names(de_results),
    total_significant_genes = sum(sapply(de_results, function(x) x$n_total_significant)),
    timecourse_genes = nrow(timecourse_results),
    stage_specific_analysis = lapply(stage_specific_genes, function(x) x$n_genes),
    gsea_performed = !is.null(gsea_results$error),
    top_go_terms = if (!is.null(gsea_results$P12_vs_P1$go_enrichment)) {
      head(gsea_results$P12_vs_P1$go_enrichment@result$Description, 10)
    } else {
      "GSEA analysis failed"
    }
  )
  
  return(report)
}
