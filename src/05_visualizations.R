#!/usr/bin/env Rscript
# Satellite Cell Development Transcriptomics Pipeline
# Script 05: Advanced Visualizations and Trajectory Analysis (Simplified & Robust)
# MIT-Level Research Pipeline

suppressPackageStartupMessages({
  library(ggplot2)
  library(viridis)
  library(umap)
  library(Rtsne)
  library(dbscan)
  library(yaml)
  library(logging)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
  library(gridExtra)
  library(cowplot)
  library(rmarkdown)
})

# --------------------------------------------------
# Logging
# --------------------------------------------------
logfile <- file.path("results", "logs", "05_visualizations.log")
if (!dir.exists(dirname(logfile))) dir.create(dirname(logfile), recursive = TRUE)

logging::basicConfig()  # set up default logger
addHandler(writeToFile, file = logfile, level = "INFO")

loginfo("Starting advanced visualizations and trajectory analysis pipeline (simplified)")

# --------------------------------------------------
# Load configuration & data
# --------------------------------------------------
config <- yaml::read_yaml("config.yaml")

loginfo("Loading preprocessed data...")
processed_file <- file.path("data/processed", "preprocessed_expression.RData")
if (!file.exists(processed_file)) {
  stop("Preprocessed data not found. Run 01_preprocess.R first.")
}
load(processed_file)   # expects: final_expression, sample_info

loginfo("Loading DE results...")
de_file <- file.path("results", "de_analysis_results.RData")
if (!file.exists(de_file)) {
  stop("DE results not found. Run 02_de_analysis.R first.")
}
load(de_file)          # expects: de_results, timecourse_results, stage_specific_genes, gsea_results

loginfo("Loading enrichment results...")
enrichment_file <- file.path("results", "enrichment_analysis_results.RData")
if (file.exists(enrichment_file)) {
  load(enrichment_file)  # expects: enrichment_results
} else {
  enrichment_results <- list()
  logwarn("Enrichment results file not found. Skipping enrichment-based plots.")
}

loginfo("Loading network results...")
network_file <- file.path("results", "network_inference_results.RData")
if (file.exists(network_file)) {
  load(network_file)      # expects: network_results
} else {
  network_results <- list()
  logwarn("Network results file not found. Skipping network-based plots.")
}

# --------------------------------------------------
# Step 1: Manifold learning (UMAP + t-SNE only)
# --------------------------------------------------
loginfo("Performing manifold learning (UMAP, t-SNE)...")

umap_results <- perform_umap_analysis(final_expression, sample_info)
tsne_results <- perform_tsne_analysis(final_expression, sample_info)

# --------------------------------------------------
# Step 2: Simple pseudotime (stage order based)
# --------------------------------------------------
loginfo("Computing simple pseudotime based on stage ordering...")
pseudotime_results <- simple_pseudotime(sample_info)

# --------------------------------------------------
# Step 3: Developmental trajectory mapping
# --------------------------------------------------
loginfo("Creating developmental trajectory maps...")
trajectory_maps <- create_developmental_trajectory_maps(final_expression, sample_info, pseudotime_results)

# --------------------------------------------------
# Step 4: Cell state discovery (kmeans, hclust, dbscan)
# --------------------------------------------------
loginfo("Discovering cell states...")
cell_states <- discover_cell_states(final_expression, sample_info)

# --------------------------------------------------
# Step 5: Plot creation
# --------------------------------------------------
loginfo("Creating manifold plots...")
manifold_plots <- create_manifold_plots(umap_results, tsne_results, sample_info)

loginfo("Creating trajectory plots...")
trajectory_plots <- create_trajectory_plots(pseudotime_results, trajectory_maps)

loginfo("Creating cell state plots...")
cell_state_plots <- create_cell_state_plots(cell_states)

loginfo("Creating integrated multi-panel figures...")
integrated_plots <- create_integrated_figures(manifold_plots, trajectory_plots, cell_state_plots)

loginfo("Creating publication-style figures...")
publication_figures <- create_publication_figures(integrated_plots, de_results, enrichment_results, network_results)

# --------------------------------------------------
# Step 6: Save all plots & interactive data
# --------------------------------------------------
if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)

loginfo("Saving manifold, trajectory, and cell state plots...")
save_individual_plots(manifold_plots, "manifold_learning")
save_individual_plots(trajectory_plots, "trajectory_analysis")
save_individual_plots(cell_state_plots, "cell_state_analysis")
save_individual_plots(integrated_plots, "integrated_analysis")
save_publication_figures(publication_figures)

loginfo("Preparing interactive data object...")
interactive_data <- prepare_interactive_data(umap_results, tsne_results, pseudotime_results, cell_states, sample_info)
interactive_file <- file.path("results", "interactive_visualization_data.RData")
save(interactive_data, file = interactive_file)

# --------------------------------------------------
# Step 7: Visualization report
# --------------------------------------------------
loginfo("Generating visualization report...")
visualization_report <- generate_visualization_report(manifold_plots, pseudotime_results, cell_states, publication_figures)
report_file <- file.path("results", "visualization_report.html")

rmarkdown::render(
  input = system.file("rmarkdown", "templates", "html_vignette", "resources", "vignette.Rmd", package = "rmarkdown"),
  output_file = report_file,
  params = list(report_data = visualization_report),
  quiet = TRUE
)

loginfo("Advanced visualizations and trajectory analysis pipeline completed successfully")

# ==================================================
# Helper functions
# ==================================================

perform_umap_analysis <- function(expression_matrix, sample_info) {
  loginfo("Running UMAP...")
  expr_t <- t(expression_matrix)

  set.seed(123)
  umap_result <- umap(expr_t, n_neighbors = 7, min_dist = 0.3, n_components = 2)

  umap_df <- data.frame(
    UMAP1 = umap_result$layout[, 1],
    UMAP2 = umap_result$layout[, 2],
    stage = sample_info$stage,
    replicate = sample_info$replicate
  )

  list(
    embedding = umap_result,
    coordinates = umap_df,
    parameters = list(n_neighbors = 7, min_dist = 0.3)
  )
}

perform_tsne_analysis <- function(expression_matrix, sample_info) {
  loginfo("Running t-SNE...")
  expr_t <- t(expression_matrix)

  set.seed(123)
  tsne_result <- Rtsne(expr_t, dims = 2, perplexity = 3, theta = 0.5, verbose = FALSE)

  tsne_df <- data.frame(
    tSNE1 = tsne_result$Y[, 1],
    tSNE2 = tsne_result$Y[, 2],
    stage = sample_info$stage,
    replicate = sample_info$replicate
  )

  list(
    embedding = tsne_result,
    coordinates = tsne_df,
    parameters = list(perplexity = 3, theta = 0.5)
  )
}

simple_pseudotime <- function(sample_info) {
  # Assign pseudotime based on stage: P1 = 1, P12 = 2, P28 = 3
  stage_levels <- c("P1", "P12", "P28")
  pt_values <- as.numeric(factor(sample_info$stage, levels = stage_levels))

  data.frame(
    cell_id = sample_info$sample_id %||% seq_len(nrow(sample_info)),
    pseudotime = pt_values,
    stage = sample_info$stage,
    replicate = sample_info$replicate
  )
}

create_developmental_trajectory_maps <- function(expression_matrix, sample_info, pseudotime_results) {
  loginfo("Building developmental trajectory data...")

  stage_order <- c("P1", "P12", "P28")
  key_genes <- head(rownames(expression_matrix), 50)

  stage_means <- lapply(stage_order, function(stage) {
    idx <- sample_info$stage == stage
    rowMeans(expression_matrix[key_genes, idx, drop = FALSE])
  })
  names(stage_means) <- stage_order

  trajectory_data <- do.call(rbind, lapply(seq_along(stage_order), function(i) {
    stage <- stage_order[i]
    data.frame(
      stage = stage,
      time_point = i,
      mean_expression = stage_means[[stage]],
      gene = names(stage_means[[stage]]),
      stringsAsFactors = FALSE
    )
  }))

  list(
    stage_trajectory = trajectory_data,
    pseudotime_mapping = aggregate(pseudotime ~ stage, data = pseudotime_results, mean)
  )
}

discover_cell_states <- function(expression_matrix, sample_info) {
  loginfo("Running clustering for cell state discovery...")
  expr_t <- t(expression_matrix)

  set.seed(123)
  kmeans_result <- kmeans(expr_t, centers = 3, nstart = 25)

  hc_result <- hclust(dist(expr_t), method = "ward.D2")
  hc_clusters <- cutree(hc_result, k = 3)

  set.seed(123)
  dbscan_result <- dbscan(expr_t, eps = 0.5, minPts = 2)

  kmeans_df <- data.frame(
    cell_id = colnames(expression_matrix),
    cluster = kmeans_result$cluster,
    stage = sample_info$stage
  )

  hierarchical_df <- data.frame(
    cell_id = colnames(expression_matrix),
    cluster = hc_clusters,
    stage = sample_info$stage
  )

  dbscan_df <- data.frame(
    cell_id = colnames(expression_matrix),
    cluster = dbscan_result$cluster,
    stage = sample_info$stage
  )

  consensus_clusters <- calculate_consensus_clusters(list(
    kmeans = kmeans_df,
    hierarchical = hierarchical_df,
    dbscan = dbscan_df
  ))

  consensus_df <- data.frame(
    cell_id = colnames(expression_matrix),
    cluster = consensus_clusters,
    stage = sample_info$stage
  )

  list(
    kmeans = kmeans_df,
    hierarchical = hierarchical_df,
    dbscan = dbscan_df,
    consensus = consensus_df
  )
}

calculate_consensus_clusters <- function(cell_states) {
  cluster_matrix <- cbind(
    kmeans = cell_states$kmeans$cluster,
    hierarchical = cell_states$hierarchical$cluster,
    dbscan = cell_states$dbscan$cluster
  )

  apply(cluster_matrix, 1, function(x) {
    tab <- table(x)
    as.numeric(names(tab)[which.max(tab)])
  })
}

create_manifold_plots <- function(umap_results, tsne_results, sample_info) {
  plots <- list()

  if (!is.null(umap_results)) {
    p <- ggplot(umap_results$coordinates, aes(x = UMAP1, y = UMAP2, color = stage)) +
      geom_point(size = 5, alpha = 0.85) +
      geom_text(aes(label = replicate), vjust = -1, size = 3) +
      scale_color_viridis_d() +
      labs(
        title = "UMAP: Satellite Cell Development",
        subtitle = "P1 → P12 → P28",
        x = "UMAP1", y = "UMAP2"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
    plots$umap <- p
  }

  if (!is.null(tsne_results)) {
    p <- ggplot(tsne_results$coordinates, aes(x = tSNE1, y = tSNE2, color = stage)) +
      geom_point(size = 5, alpha = 0.85) +
      geom_text(aes(label = replicate), vjust = -1, size = 3) +
      scale_color_viridis_d() +
      labs(
        title = "t-SNE: Satellite Cell Development",
        x = "tSNE1", y = "tSNE2"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
    plots$tsne <- p
  }

  plots
}

create_trajectory_plots <- function(pseudotime_results, trajectory_maps) {
  plots <- list()

  if (!is.null(pseudotime_results)) {
    p <- ggplot(pseudotime_results, aes(x = stage, y = pseudotime, fill = stage)) +
      geom_boxplot(alpha = 0.7) +
      geom_jitter(width = 0.1, alpha = 0.7) +
      scale_fill_viridis_d() +
      labs(
        title = "Pseudotime by Developmental Stage",
        x = "Stage", y = "Pseudotime"
      ) +
      theme_minimal() +
      theme(legend.position = "none")
    plots$pseudotime <- p
  }

  if (!is.null(trajectory_maps$stage_trajectory)) {
    p <- ggplot(trajectory_maps$stage_trajectory,
                aes(x = time_point, y = mean_expression, group = gene, color = gene)) +
      geom_line(alpha = 0.5) +
      scale_x_continuous(breaks = c(1, 2, 3), labels = c("P1", "P12", "P28")) +
      labs(
        title = "Gene Expression Trajectories (subset of genes)",
        x = "Stage", y = "Mean expression"
      ) +
      theme_minimal() +
      theme(legend.position = "none")
    plots$stage_trajectory <- p
  }

  plots
}

create_cell_state_plots <- function(cell_states) {
  plots <- list()

  if (!is.null(cell_states$consensus)) {
    p <- ggplot(cell_states$consensus,
                aes(x = stage, fill = factor(cluster))) +
      geom_bar(position = "fill") +
      scale_fill_viridis_d(name = "Cell state") +
      labs(
        title = "Consensus Cell States by Stage",
        x = "Stage", y = "Proportion"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
    plots$consensus <- p
  }

  if (!is.null(cell_states$kmeans)) {
    p <- ggplot(cell_states$kmeans,
                aes(x = stage, fill = factor(cluster))) +
      geom_bar(position = "dodge") +
      scale_fill_viridis_d(name = "Cell state") +
      labs(
        title = "K-means Cell States by Stage",
        x = "Stage", y = "Count"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
    plots$kmeans <- p
  }

  plots
}

create_integrated_figures <- function(manifold_plots, trajectory_plots, cell_state_plots) {
  integrated <- list()

  if (length(manifold_plots) > 0) {
    integrated$manifold_panel <- gridExtra::grid.arrange(
      grobs = manifold_plots,
      ncol = 2,
      top = "Manifold Learning (UMAP, t-SNE)"
    )
  }

  if (length(trajectory_plots) > 0) {
    integrated$trajectory_panel <- gridExtra::grid.arrange(
      grobs = trajectory_plots,
      ncol = 1,
      top = "Developmental Trajectory Analysis"
    )
  }

  if (length(cell_state_plots) > 0) {
    integrated$cell_state_panel <- gridExtra::grid.arrange(
      grobs = cell_state_plots,
      ncol = 1,
      top = "Cell State Discovery"
    )
  }

  integrated
}

create_publication_figures <- function(integrated_plots, de_results, enrichment_results, network_results) {
  figs <- list()

  # Figure 1 – overview
  figs$figure_1_overview <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 1),
              fill = "grey90", color = NA) +
    annotate("text", x = 0.5, y = 0.75, label = "Satellite Cell Development", size = 7, fontface = "bold") +
    annotate("text", x = 0.5, y = 0.55, label = "P1 → P12 → P28", size = 6) +
    annotate("text", x = 0.5, y = 0.35, label = "Transcriptomic Systems Pipeline", size = 5) +
    theme_void() +
    labs(title = "Figure 1: Experimental Overview")

  # Figure 2 – DE summary
  de_summary <- data.frame(
    contrast = names(de_results),
    upregulated = sapply(de_results, function(x) x$n_upregulated),
    downregulated = sapply(de_results, function(x) x$n_downregulated),
    total_significant = sapply(de_results, function(x) x$n_total_significant)
  )

  figs$figure_2_differential_expression <- ggplot(de_summary,
                                                  aes(x = contrast, y = total_significant)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_text(aes(label = total_significant), vjust = -0.5, size = 4) +
    theme_minimal() +
    labs(
      title = "Differential Expression Summary",
      x = "Contrast", y = "Number of significant genes"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Figure 3 – manifold & trajectories if available
  if (!is.null(integrated_plots$manifold_panel)) {
    figs$figure_3_manifold <- integrated_plots$manifold_panel
  }
  if (!is.null(integrated_plots$trajectory_panel)) {
    figs$figure_4_trajectory <- integrated_plots$trajectory_panel
  }

  figs
}

save_individual_plots <- function(plots, prefix) {
  if (!length(plots)) return(invisible(NULL))
  if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)

  for (nm in names(plots)) {
    f <- file.path("figures", paste0(prefix, "_", nm, ".pdf"))
    ggsave(f, plots[[nm]], width = 8, height = 6, dpi = 300)
  }
}

save_publication_figures <- function(publication_figures) {
  if (!length(publication_figures)) return(invisible(NULL))
  if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)

  for (nm in names(publication_figures)) {
    f <- file.path("figures", paste0(nm, ".pdf"))
    ggsave(f, publication_figures[[nm]], width = 10, height = 7, dpi = 300)
  }
}

prepare_interactive_data <- function(umap_results, tsne_results, pseudotime_results, cell_states, sample_info) {
  list(
    umap = if (!is.null(umap_results)) umap_results$coordinates else NULL,
    tsne = if (!is.null(tsne_results)) tsne_results$coordinates else NULL,
    pseudotime = pseudotime_results,
    cell_states = cell_states$consensus,
    metadata = sample_info
  )
}

generate_visualization_report <- function(manifold_results, pseudotime_results, cell_states, publication_figures) {
  list(
    timestamp = Sys.time(),
    manifold_methods_used = length(manifold_results),
    pseudotime_analysis_performed = !is.null(pseudotime_results),
    cell_states_identified = if (!is.null(cell_states$consensus)) {
      length(unique(cell_states$consensus$cluster))
    } else {
      0
    },
    publication_figures_generated = length(publication_figures),
    visualization_quality = "Publication-ready (simplified)",
    interactive_data_prepared = TRUE
  )
}
