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
  library(gridExtra)
})

# ============================
# Helper functions (define FIRST)
# ============================

create_deseq_dataset <- function(expression_matrix, sample_info) {
  # Create count matrix (approximate counts from log2 expr)
  loginfo("Reconstructing approximate counts from log2-normalized expression...")
  count_matrix <- round(2^expression_matrix - 1)
  count_matrix[count_matrix < 0] <- 0

  # DESeq2 design
  sample_conditions <- factor(sample_info$stage, levels = c("P1", "P12", "P28"))
  colData <- S4Vectors::DataFrame(stage = sample_conditions)

  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData   = colData,
    design    = ~ stage
  )

  # Filter low-count genes
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]

  loginfo(paste("Kept", nrow(dds), "genes after count filtering."))

  dds <- DESeq(dds)
  return(dds)
}

determine_trend <- function(linear_coef, quadratic_coef) {
  if (is.na(linear_coef) || is.na(quadratic_coef)) return("unknown")
  if (abs(quadratic_coef) > abs(linear_coef)) {
    if (quadratic_coef > 0) "quadratic_up" else "quadratic_down"
  } else {
    if (linear_coef > 0) "linear_up" else "linear_down"
  }
}

perform_timecourse_analysis <- function(expression_matrix, sample_info) {
  loginfo("Running polynomial time-course models for each gene (this may take a while)...")

  time_points <- ifelse(
    sample_info$stage == "P1", 1,
    ifelse(sample_info$stage == "P12", 12, 28)
  )

  timecourse_results <- lapply(rownames(expression_matrix), function(gene) {
    gene_expr <- as.numeric(expression_matrix[gene, ])

    model <- tryCatch(
      lm(gene_expr ~ poly(time_points, 2)),
      error = function(e) NULL
    )

    if (is.null(model)) {
      return(list(
        gene           = gene,
        linear_coef    = NA_real_,
        quadratic_coef = NA_real_,
        linear_pval    = NA_real_,
        quadratic_pval = NA_real_,
        r_squared      = NA_real_,
        trend          = "failed"
      ))
    }

    anova_result <- anova(model)

    list(
      gene           = gene,
      linear_coef    = coef(model)[2],
      quadratic_coef = coef(model)[3],
      linear_pval    = anova_result$"Pr(>F)"[2],
      quadratic_pval = anova_result$"Pr(>F)"[3],
      r_squared      = summary(model)$r.squared,
      trend          = determine_trend(coef(model)[2], coef(model)[3])
    )
  })

  timecourse_df <- do.call(rbind, lapply(timecourse_results, function(x) {
    data.frame(
      gene           = x$gene,
      linear_coef    = x$linear_coef,
      quadratic_coef = x$quadratic_coef,
      linear_pval    = x$linear_pval,
      quadratic_pval = x$quadratic_pval,
      r_squared      = x$r_squared,
      trend          = x$trend,
      stringsAsFactors = FALSE
    )
  }))

  timecourse_df$linear_padj    <- p.adjust(timecourse_df$linear_pval,    method = "BH")
  timecourse_df$quadratic_padj <- p.adjust(timecourse_df$quadratic_pval, method = "BH")

  timecourse_df
}

identify_stage_specific_genes <- function(expression_matrix, sample_info, de_results) {
  stage_specific <- list()

  sig_genes_all <- unique(unlist(lapply(de_results, function(x) {
    rownames(x$significant_genes)
  })))

  for (stage in unique(sample_info$stage)) {
    stage_samples <- sample_info$stage == stage
    other_samples <- sample_info$stage != stage

    stage_means <- rowMeans(expression_matrix[, stage_samples, drop = FALSE])
    other_means <- rowMeans(expression_matrix[, other_samples, drop = FALSE])

    fold_changes <- log2(stage_means + 1) - log2(other_means + 1)

    stage_genes <- names(fold_changes)[
      which(
        fold_changes > 1 &
          names(fold_changes) %in% sig_genes_all
      )
    ]

    stage_specific[[stage]] <- list(
      genes     = stage_genes,
      n_genes   = length(stage_genes),
      top_genes = head(
        stage_genes[order(fold_changes[stage_genes], decreasing = TRUE)],
        20
      )
    )
  }

  stage_specific
}

create_volcano_plots <- function(de_results, padj_threshold, log2fc_threshold) {
  volcano_plots <- list()

  for (contrast_name in names(de_results)) {
    res <- de_results[[contrast_name]]$full_results

    volcano_plot <- EnhancedVolcano(
      as.data.frame(res),
      lab      = rownames(res),
      x        = "log2FoldChange",
      y        = "padj",
      title    = paste("Volcano Plot:", contrast_name),
      pCutoff  = padj_threshold,
      FCcutoff = log2fc_threshold,
      pointSize      = 2.0,
      labSize        = 3.5,
      legendPosition = "right",
      legendLabSize  = 10,
      legendIconSize = 3.0,
      col            = c("black", "red", "green3", "blue"),
      colAlpha       = 0.8
    )

    volcano_plots[[contrast_name]] <- volcano_plot
  }

  combined_plot <- gridExtra::grid.arrange(grobs = volcano_plots, ncol = 2)

  list(
    individual_plots = volcano_plots,
    combined_plot    = combined_plot
  )
}

create_de_heatmaps <- function(expression_matrix, sample_info, de_results) {
  all_sig_genes <- unique(unlist(lapply(de_results, function(x) {
    rownames(x$significant_genes)
  })))

  if (length(all_sig_genes) == 0) {
    warning("No significant genes found for heatmap; skipping heatmap generation.")
    return(NULL)
  }

  all_sig_genes <- intersect(all_sig_genes, rownames(expression_matrix))
  if (length(all_sig_genes) == 0) {
    warning("Significant genes not found in expression matrix; skipping heatmap.")
    return(NULL)
  }

  sig_expression <- expression_matrix[all_sig_genes, , drop = FALSE]

  stage_order <- order(factor(sample_info$stage, levels = c("P1", "P12", "P28")))
  sig_expression_ordered <- sig_expression[, stage_order, drop = FALSE]
  sample_info_ordered    <- sample_info[stage_order, , drop = FALSE]

  sig_expression_zscore <- t(scale(t(sig_expression_ordered)))

  annotation_col <- data.frame(stage = sample_info_ordered$stage)
  rownames(annotation_col) <- rownames(sample_info_ordered)

  heatmap_plot <- pheatmap(
    sig_expression_zscore,
    annotation_col = annotation_col,
    show_rownames  = FALSE,
    cluster_rows   = TRUE,
    cluster_cols   = FALSE,
    main           = "Differential Expression Heatmap",
    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
  )

  list(
    heatmap_plot       = heatmap_plot,
    sig_genes          = all_sig_genes,
    sig_expression_z   = sig_expression_zscore
  )
}

perform_gsea_analysis <- function(de_results, padj_threshold) {
  gsea_results <- list()

  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Mm.eg.db)
  })

  for (contrast_name in names(de_results)) {
    res <- de_results[[contrast_name]]$full_results

    gene_list <- res$log2FoldChange
    names(gene_list) <- rownames(res)
    gene_list <- sort(gene_list, decreasing = TRUE)

    up_genes <- names(gene_list)[which(gene_list > 0)]
    up_genes <- bitr(up_genes, fromType = "SYMBOL",
                     toType = "ENTREZID", OrgDb = org.Mm.eg.db)
    up_entrez <- unique(up_genes$ENTREZID)

    go_enrich   <- NULL
    kegg_enrich <- NULL

    if (length(up_entrez) > 10) {
      go_enrich <- tryCatch(
        enrichGO(
          gene          = up_entrez,
          OrgDb         = org.Mm.eg.db,
          ont           = "BP",
          pAdjustMethod = "BH",
          pvalueCutoff  = padj_threshold,
          qvalueCutoff  = padj_threshold,
          readable      = TRUE
        ),
        error = function(e) NULL
      )

      kegg_enrich <- tryCatch(
        enrichKEGG(
          gene          = up_entrez,
          organism      = "mmu",
          pvalueCutoff  = padj_threshold,
          qvalueCutoff  = padj_threshold
        ),
        error = function(e) NULL
      )
    }

    gsea_results[[contrast_name]] <- list(
      go_enrichment   = go_enrich,
      kegg_enrichment = kegg_enrich,
      gene_list       = gene_list
    )
  }

  gsea_results
}

create_expression_trajectories <- function(expression_matrix, sample_info, de_results) {
  top_genes <- unique(unlist(lapply(de_results, function(x) {
    head(rownames(x$significant_genes), 50)
  })))

  top_genes <- intersect(top_genes, rownames(expression_matrix))
  if (length(top_genes) == 0) {
    warning("No top genes for trajectories; skipping trajectory plot.")
    return(NULL)
  }

  time_points <- ifelse(
    sample_info$stage == "P1", 1,
    ifelse(sample_info$stage == "P12", 12, 28)
  )

  trajectory_data <- lapply(top_genes, function(gene) {
    gene_expr <- expression_matrix[gene, ]
    stage_means <- aggregate(
      as.numeric(gene_expr),
      list(stage = sample_info$stage),
      mean
    )
    colnames(stage_means) <- c("stage", "expression")
    stage_means$time_point <- c(1, 12, 28)
    stage_means$gene       <- gene
    stage_means
  })

  all_trajectories <- do.call(rbind, trajectory_data)

  trajectory_plot <- ggplot(
    all_trajectories,
    aes(x = time_point, y = expression, group = gene, color = gene)
  ) +
    geom_line(size = 1.0, alpha = 0.6) +
    geom_point(size = 2) +
    labs(
      title = "Gene Expression Trajectories During Satellite Cell Development",
      x     = "Postnatal Day",
      y     = "Normalized Expression"
    ) +
    theme_minimal() +
    scale_x_continuous(
      breaks = c(1, 12, 28),
      labels = c("P1", "P12", "P28")
    ) +
    theme(legend.position = "none")

  list(
    trajectory_plot = trajectory_plot,
    trajectory_data = all_trajectories
  )
}

generate_de_report <- function(de_results, timecourse_results, stage_specific_genes, gsea_results) {
  # Figure out if GSEA actually ran
  gsea_ok <- length(gsea_results) > 0 && !("error" %in% names(gsea_results))

  top_go_terms <- "GSEA not available"
  if (gsea_ok && "P12_vs_P1" %in% names(gsea_results)) {
    go_obj <- gsea_results$P12_vs_P1$go_enrichment
    if (!is.null(go_obj) && nrow(go_obj@result) > 0) {
      top_go_terms <- head(go_obj@result$Description, 10)
    }
  }

  list(
    timestamp               = Sys.time(),
    contrasts_analyzed      = names(de_results),
    total_significant_genes = sum(sapply(de_results, function(x) x$n_total_significant)),
    timecourse_genes        = nrow(timecourse_results),
    stage_specific_analysis = lapply(stage_specific_genes, function(x) x$n_genes),
    gsea_performed          = gsea_ok,
    top_go_terms            = top_go_terms
  )
}

# ============================
# Logging configuration
# ============================

logfile <- file.path("results", "logs", "02_de_analysis.log")
if (!dir.exists(dirname(logfile))) dir.create(dirname(logfile), recursive = TRUE)

basicConfig()  # console
addHandler(writeToFile, file = logfile, level = "INFO")

loginfo("Starting differential expression analysis pipeline")

# ============================
# Load configuration & data
# ============================

if (!file.exists("config.yaml")) {
  logerror("config.yaml not found in project root.")
  stop("config.yaml not found. Run from project root.")
}

config <- yaml::read_yaml("config.yaml")

loginfo("Loading preprocessed data...")
processed_file <- file.path("data", "processed", "preprocessed_expression.RData")
if (!file.exists(processed_file)) {
  stop("Preprocessed data not found: ", processed_file, "\nRun 01_preprocess.R first.")
}
load(processed_file)  # loads final_expression, sample_info, gene_mapping, qc_results

if (!exists("final_expression") || !exists("sample_info")) {
  stop("final_expression or sample_info missing from preprocessed data.")
}

# ============================
# Create DESeq2 dataset
# ============================

loginfo("Creating DESeq2 dataset from preprocessed expression...")
dds <- create_deseq_dataset(final_expression, sample_info)

# ============================
# Differential expression per contrast
# ============================

loginfo("Performing differential expression analysis...")

contrasts <- list(
  P12_vs_P1   = c("stage", "P12", "P1"),
  P28_vs_P12  = c("stage", "P28", "P12"),
  P28_vs_P1   = c("stage", "P28", "P1")
)

de_results <- list()

for (contrast_name in names(contrasts)) {
  loginfo(paste("Analyzing contrast:", contrast_name))

  contrast <- contrasts[[contrast_name]]

  res <- results(
    dds,
    contrast = contrast,
    alpha    = config$analysis$padj_threshold
  )

  res_df <- as.data.frame(res)
  res_df <- res_df[complete.cases(res_df$log2FoldChange) & !is.na(res_df$padj), ]

  sig_genes <- res_df[
    which(
      res_df$padj < config$analysis$padj_threshold &
        abs(res_df$log2FoldChange) > config$analysis$log2fc_threshold
    ),
  ]

  de_results[[contrast_name]] <- list(
    full_results      = res_df,
    significant_genes = sig_genes,
    n_upregulated     = sum(
      res_df$log2FoldChange > config$analysis$log2fc_threshold &
        res_df$padj < config$analysis$padj_threshold
    ),
    n_downregulated   = sum(
      res_df$log2FoldChange < -config$analysis$log2fc_threshold &
        res_df$padj < config$analysis$padj_threshold
    ),
    n_total_significant = nrow(sig_genes)
  )

  loginfo(paste("Found", nrow(sig_genes), "significant genes for", contrast_name))

  # Save per-contrast RDS for Rmd report usage
  rds_name <- switch(
    contrast_name,
    P12_vs_P1  = "results/de_results_P12_vs_P1.rds",
    P28_vs_P12 = "results/de_results_P28_vs_P12.rds",
    P28_vs_P1  = "results/de_results_P28_vs_P1.rds"
  )
  if (!is.null(rds_name)) {
    saveRDS(res_df, file = rds_name)
    loginfo(paste("Saved DE results for", contrast_name, "to", rds_name))
  }
}

# ============================
# Time-course & stage-specific
# ============================

loginfo("Performing time-course differential expression analysis...")
timecourse_results <- perform_timecourse_analysis(final_expression, sample_info)

loginfo("Identifying stage-specific gene expression patterns...")
stage_specific_genes <- identify_stage_specific_genes(final_expression, sample_info, de_results)

# ============================
# Plots
# ============================

loginfo("Creating volcano plots...")
volcano_plots <- create_volcano_plots(
  de_results,
  config$analysis$padj_threshold,
  config$analysis$log2fc_threshold
)

loginfo("Creating differential expression heatmaps...")
de_heatmaps <- create_de_heatmaps(final_expression, sample_info, de_results)

loginfo("Creating gene expression trajectories...")
trajectory_plots <- create_expression_trajectories(final_expression, sample_info, de_results)

# ============================
# GSEA
# ============================

loginfo("Performing gene set enrichment analysis...")
gsea_results <- tryCatch(
  perform_gsea_analysis(de_results, config$analysis$padj_threshold),
  error = function(e) {
    logwarn(paste("GSEA analysis failed:", e$message))
    list(error = e$message)
  }
)

# ============================
# Save results & figures
# ============================

loginfo("Saving differential expression results...")

if (!dir.exists("results"))  dir.create("results", recursive = TRUE)
if (!dir.exists("figures"))  dir.create("figures", recursive = TRUE)

de_file <- file.path("results", "de_analysis_results.RData")
save(de_results, timecourse_results, stage_specific_genes, gsea_results, file = de_file)
loginfo(paste("DE analysis results saved to:", de_file))

# Volcano plots
volcano_file <- file.path("figures", "volcano_plots.pdf")
ggsave(
  filename = volcano_file,
  plot     = volcano_plots$combined_plot,
  width    = 16,
  height   = 12,
  dpi      = config$visualization$figure_dpi
)
loginfo(paste("Volcano plots saved to:", volcano_file))

# Heatmaps (if available)
heatmap_file <- file.path("figures", "de_heatmaps.pdf")
if (!is.null(de_heatmaps)) {
  pdf(heatmap_file, width = 10, height = 8)
  print(de_heatmaps$heatmap_plot)
  dev.off()
  loginfo(paste("DE heatmaps saved to:", heatmap_file))
} else {
  logwarn("No DE heatmaps generated; skipping heatmap PDF.")
}

# Trajectories (if available)
trajectory_file <- file.path("figures", "expression_trajectories.pdf")
if (!is.null(trajectory_plots)) {
  ggsave(
    filename = trajectory_file,
    plot     = trajectory_plots$trajectory_plot,
    width    = 14,
    height   = 10,
    dpi      = config$visualization$figure_dpi
  )
  loginfo(paste("Expression trajectories saved to:", trajectory_file))
} else {
  logwarn("No trajectory plots generated; skipping trajectory PDF.")
}

# ============================
# DE report (optional template)
# ============================

loginfo("Attempting to generate DE analysis report (if template exists)...")

de_report <- generate_de_report(
  de_results,
  timecourse_results,
  stage_specific_genes,
  gsea_results
)

report_tpl  <- "reports/de_analysis_report.Rmd"
report_file <- file.path("results", "de_analysis_report.html")

if (file.exists(report_tpl)) {
  tryCatch({
    rmarkdown::render(
      input       = report_tpl,
      output_file = report_file,
      params      = list(report_data = de_report),
      quiet       = TRUE
    )
    loginfo(paste("DE analysis report generated:", report_file))
  }, error = function(e) {
    logwarn(paste("Could not render DE report:", conditionMessage(e)))
  })
} else {
  logwarn(paste("DE report template not found at", report_tpl, "- skipping report generation."))
}

loginfo("Differential expression analysis pipeline completed successfully")
