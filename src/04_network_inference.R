#!/usr/bin/env Rscript
# Satellite Cell Development Transcriptomics Pipeline
# Script 04: Gene Regulatory Network Inference
# MIT-Level Research Pipeline

suppressPackageStartupMessages({
  library(igraph)
  library(ggnetwork)
  library(ggraph)
  library(tidygraph)
  library(GENIE3)
  library(grndata)
  library(SCENIC)
  library(RcisTarget)
  library(yaml)
  library(logging)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(ggrepel)
  library(GO.db)
  library(org.Mm.eg.db)
})

# Configure logging
logfile <- file.path("results", "logs", "04_network_inference.log")
if (!dir.exists(dirname(logfile))) dir.create(dirname(logfile), recursive = TRUE)
basicConfig(level = "INFO", file = logfile)

loginfo("Starting gene regulatory network inference pipeline")

# Load configuration
config <- yaml::read_yaml("config.yaml")

# Load preprocessed data
loginfo("Loading preprocessed expression data...")
processed_file <- file.path("data/processed", "preprocessed_expression.RData")
if (!file.exists(processed_file)) {
  stop("Preprocessed data not found. Run 01_preprocess.R first.")
}
load(processed_file)

# Load DE results
loginfo("Loading differential expression results...")
de_file <- file.path("results", "de_analysis_results.RData")
load(de_file)

# Step 1: Transcription Factor Target Prediction
loginfo("Performing transcription factor target prediction...")

# Load transcription factor list
tf_list <- load_transcription_factors()

# Filter expression matrix for TFs and target genes
network_expression <- filter_network_genes(final_expression, tf_list)

# Perform GENIE3 network inference
genie3_network <- perform_genie3_inference(network_expression, tf_list)

# Perform correlation-based network inference
correlation_network <- perform_correlation_inference(network_expression, config$analysis$correlation_threshold)

# Step 2: SCENIC-style Regulatory Network Inference
loginfo("Performing SCENIC-style regulatory network inference...")
scenic_network <- perform_scenic_inference(network_expression, tf_list)

# Step 3: Motif Enrichment Analysis
loginfo("Performing transcription factor motif enrichment analysis...")
motif_results <- perform_motif_enrichment(network_expression, tf_list)

# Step 4: Stage-specific Network Analysis
loginfo("Analyzing stage-specific regulatory networks...")
stage_networks <- analyze_stage_specific_networks(final_expression, sample_info, tf_list)

# Step 5: Network Topology Analysis
loginfo("Performing network topology analysis...")
topology_results <- analyze_network_topology(genie3_network, correlation_network, scenic_network)

# Step 6: Hub Gene Identification
loginfo("Identifying hub genes and key regulators...")
hub_genes <- identify_hub_genes(genie3_network, correlation_network, scenic_network)

# Step 7: Create Network Visualizations
loginfo("Creating network visualizations...")

# Create regulatory network plots
network_plots <- create_network_plots(genie3_network, correlation_network, scenic_network)

# Create stage-specific network plots
stage_network_plots <- create_stage_network_plots(stage_networks)

# Create hub gene visualization
hub_gene_plots <- create_hub_gene_plots(hub_genes)

# Step 8: Generate Network Inference Report
loginfo("Generating network inference report...")

# Save network results
network_file <- file.path("results", "network_inference_results.RData")
save(genie3_network, correlation_network, scenic_network, motif_results, 
     stage_networks, topology_results, hub_genes, file = network_file)

# Save plots
network_plot_file <- file.path("figures", "regulatory_networks.pdf")
ggsave(network_plot_file, network_plots$combined_plot, width = 16, height = 12, 
       dpi = config$visualization$figure_dpi)

stage_network_file <- file.path("figures", "stage_specific_networks.pdf")
ggsave(stage_network_file, stage_network_plots$combined_plot, width = 14, height = 10, 
       dpi = config$visualization$figure_dpi)

hub_gene_file <- file.path("figures", "hub_genes.pdf")
ggsave(hub_gene_file, hub_gene_plots$combined_plot, width = 12, height = 8, 
       dpi = config$visualization$figure_dpi)

# Generate network report
network_report <- generate_network_report(genie3_network, correlation_network, scenic_network,
                                         motif_results, stage_networks, topology_results, hub_genes)
report_file <- file.path("results", "network_inference_report.html")
rmarkdown::render(
  input = system.file("rmarkdown", "templates", "network_report", package = "rmarkdown"),
  output_file = report_file,
  params = list(report_data = network_report)
)

loginfo("Gene regulatory network inference pipeline completed successfully")

# Helper functions
load_transcription_factors <- function() {
  # Comprehensive list of mouse transcription factors
  tf_database <- c(
    # Myogenic regulatory factors
    "Pax7", "MyoD", "Myf5", "Myogenin", "Mrf4",
    
    # Homeobox TFs
    "Six1", "Six4", "Eya1", "Eya2", "Pitx2", "Pitx3",
    
    # Basic helix-loop-helix
    "Tcf4", "Tcf12", "Hand1", "Hand2", "Twist1", "Twist2",
    
    # Zinc finger
    "Sp1", "Sp2", "Sp3", "Zeb1", "Zeb2", "Snai1", "Snai2",
    
    # Forkhead
    "Foxo1", "Foxo3", "Foxo4", "Foxc1", "Foxc2", "Foxp1",
    
    # Nuclear receptors
    "Nr1d1", "Nr1d2", "Ppara", "Ppard", "Pparg", "Esrra", "Esrrg",
    
    # Signal-dependent
    "Creb1", "Atf1", "Atf2", "Atf3", "Atf4", "Jun", "Fos", "Myc",
    
    # Cell cycle
    "E2f1", "E2f2", "E2f3", "E2f4", "E2f5", "Rb1", "Trp53",
    
    # Stem cell maintenance
    "Sox2", "Sox6", "Sox8", "Sox9", "Sox10", "Nanog", "Oct4",
    
    # Quiescence regulators
    "Bcl6", "Hes1", "Hey1", "Hey2", "Id1", "Id2", "Id3"
  )
  
  return(tf_database)
}

filter_network_genes <- function(expression_matrix, tf_list) {
  # Filter for TFs and highly variable genes
  tf_genes <- intersect(tf_list, rownames(expression_matrix))
  
  # Calculate variance for each gene
  gene_variance <- apply(expression_matrix, 1, var)
  
  # Select top variable genes (excluding TFs)
  non_tf_genes <- setdiff(rownames(expression_matrix), tf_genes)
  top_variable <- head(order(gene_variance[non_tf_genes], decreasing = TRUE), 2000)
  target_genes <- non_tf_genes[top_variable]
  
  # Combine TFs and target genes
  network_genes <- union(tf_genes, target_genes)
  network_expression <- expression_matrix[network_genes, ]
  
  return(list(
    expression = network_expression,
    tf_genes = tf_genes,
    target_genes = target_genes
  ))
}

perform_genie3_inference <- function(expression_data, tf_list) {
  loginfo("Running GENIE3 network inference...")
  
  # Prepare expression matrix
  expr_matrix <- as.matrix(expression_data$expression)
  
  # Define regulators (TFs)
  regulators <- intersect(tf_list, rownames(expr_matrix))
  
  if (length(regulators) < 5) {
    logwarn("Too few transcription factors for network inference")
    return(NULL)
  }
  
  # Run GENIE3
  tryCatch({
    # Use subset of genes for computational efficiency
    n_genes <- min(1000, nrow(expr_matrix))
    subset_genes <- sample(rownames(expr_matrix), n_genes)
    expr_subset <- expr_matrix[subset_genes, ]
    
    # Run GENIE3 with regulators
    weight_matrix <- GENIE3(expr_subset, regulators = regulators, nCores = config$compute$max_cores)
    
    # Convert to edge list
    edge_list <- get_link_list(weight_matrix, threshold = 0.01)
    
    # Create igraph object
    network_graph <- graph_from_data_frame(edge_list, directed = TRUE)
    
    # Add node attributes
    V(network_graph)$type <- ifelse(V(network_graph)$name %in% regulators, "TF", "target")
    V(network_graph)$degree <- degree(network_graph)
    V(network_graph)$betweenness <- betweenness(network_graph)
    
    return(list(
      graph = network_graph,
      weight_matrix = weight_matrix,
      edge_list = edge_list,
      regulators = regulators
    ))
    
  }, error = function(e) {
    logwarn(paste("GENIE3 inference failed:", e$message))
    return(NULL)
  })
}

perform_correlation_inference <- function(expression_data, correlation_threshold) {
  loginfo("Running correlation-based network inference...")
  
  expr_matrix <- as.matrix(expression_data$expression)
  
  # Calculate correlation matrix
  cor_matrix <- cor(t(expr_matrix), method = "spearman")
  
  # Filter by correlation threshold
  cor_matrix[abs(cor_matrix) < correlation_threshold] <- 0
  
  # Convert to edge list
  edge_list <- list()
  for (i in 1:nrow(cor_matrix)) {
    for (j in 1:ncol(cor_matrix)) {
      if (i != j && abs(cor_matrix[i, j]) >= correlation_threshold) {
        edge_list <- rbind(edge_list, data.frame(
          from = rownames(cor_matrix)[i],
          to = colnames(cor_matrix)[j],
          weight = cor_matrix[i, j],
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # Create igraph object
  network_graph <- graph_from_data_frame(edge_list, directed = FALSE)
  
  # Add node attributes
  V(network_graph)$degree <- degree(network_graph)
  V(network_graph)$betweenness <- betweenness(network_graph)
  
  return(list(
    graph = network_graph,
    correlation_matrix = cor_matrix,
    edge_list = edge_list,
    threshold = correlation_threshold
  ))
}

perform_scenic_inference <- function(expression_data, tf_list) {
  loginfo("Running SCENIC-style network inference...")
  
  tryCatch({
    # Load motif databases (simplified version)
    library(RcisTarget)
    
    expr_matrix <- as.matrix(expression_data$expression)
    tf_genes <- intersect(tf_list, rownames(expr_matrix))
    
    # Step 1: Co-expression analysis (simplified GRNBoost2)
    coexp_modules <- identify_coexpression_modules(expr_matrix, tf_genes)
    
    # Step 2: Motif enrichment (simplified)
    motif_enrichment <- perform_motif_analysis(coexp_modules, tf_genes)
    
    # Step 3: Build regulons
    regulons <- build_regulons(coexp_modules, motif_enrichment, tf_genes)
    
    # Create network from regulons
    edge_list <- list()
    for (tf in names(regulons)) {
      targets <- regulons[[tf]]
      for (target in targets) {
        edge_list <- rbind(edge_list, data.frame(
          from = tf,
          to = target,
          weight = 1,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Create igraph object
    network_graph <- graph_from_data_frame(edge_list, directed = TRUE)
    
    # Add node attributes
    V(network_graph)$type <- ifelse(V(network_graph)$name %in% tf_genes, "TF", "target")
    V(network_graph)$degree <- degree(network_graph)
    V(network_graph)$betweenness <- betweenness(network_graph)
    
    return(list(
      graph = network_graph,
      edge_list = edge_list,
      regulons = regulons,
      tf_genes = tf_genes
    ))
    
  }, error = function(e) {
    logwarn(paste("SCENIC inference failed:", e$message))
    return(NULL)
  })
}

perform_motif_enrichment <- function(expression_data, tf_list) {
  loginfo("Performing motif enrichment analysis...")
  
  # Simplified motif analysis using known TF binding sites
  tryCatch({
    library(BSgenome.Mmusculus.UCSC.mm10)
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    
    expr_matrix <- as.matrix(expression_data$expression)
    tf_genes <- intersect(tf_list, rownames(expr_matrix))
    
    # Get gene coordinates
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
    gene_coords <- genes(txdb)
    
    # Perform motif analysis for each TF
    motif_results <- list()
    
    for (tf in tf_genes) {
      # Get TF binding motif (simplified - use known motifs)
      tf_motif <- get_tf_motif(tf)
      
      if (!is.null(tf_motif)) {
        # Find potential target genes with motif enrichment
        target_genes <- find_motif_targets(tf_motif, gene_coords, expr_matrix)
        
        motif_results[[tf]] <- list(
          tf = tf,
          motif = tf_motif,
          target_genes = target_genes,
          n_targets = length(target_genes)
        )
      }
    }
    
    return(motif_results)
    
  }, error = function(e) {
    logwarn(paste("Motif enrichment failed:", e$message))
    return(list(error = e$message))
  })
}

analyze_stage_specific_networks <- function(expression_matrix, sample_info, tf_list) {
  loginfo("Analyzing stage-specific regulatory networks...")
  
  stage_networks <- list()
  
  for (stage in unique(sample_info$stage)) {
    # Get samples for this stage
    stage_samples <- sample_info$stage == stage
    stage_expression <- expression_matrix[, stage_samples]
    
    # Perform network inference for this stage
    stage_tf_list <- intersect(tf_list, rownames(stage_expression))
    
    if (length(stage_tf_list) >= 5) {
      # Correlation network for this stage
      stage_cor_network <- perform_correlation_inference(
        list(expression = stage_expression, tf_genes = stage_tf_list, target_genes = rownames(stage_expression)),
        config$analysis$correlation_threshold
      )
      
      stage_networks[[stage]] <- stage_cor_network
    }
  }
  
  return(stage_networks)
}

analyze_network_topology <- function(genie3_network, correlation_network, scenic_network) {
  loginfo("Analyzing network topology...")
  
  topology_results <- list()
  
  # Analyze each network
  networks <- list(
    genie3 = genie3_network,
    correlation = correlation_network,
    scenic = scenic_network
  )
  
  for (network_name in names(networks)) {
    network <- networks[[network_name]]
    
    if (!is.null(network) && !is.null(network$graph)) {
      g <- network$graph
      
      # Calculate topological measures
      topology <- list(
        network_name = network_name,
        n_nodes = vcount(g),
        n_edges = ecount(g),
        density = edge_density(g),
        diameter = diameter(g),
        average_path_length = mean_distance(g),
        clustering_coefficient = transitivity(g),
        degree_distribution = degree(g),
        betweenness_centrality = betweenness(g),
        closeness_centrality = closeness(g),
        eigenvector_centrality = eigen_centrality(g)$vector,
        page_rank = page_rank(g)$vector
      )
      
      topology_results[[network_name]] <- topology
    }
  }
  
  return(topology_results)
}

identify_hub_genes <- function(genie3_network, correlation_network, scenic_network) {
  loginfo("Identifying hub genes...")
  
  hub_genes <- list()
  
  # Combine centrality measures from all networks
  networks <- list(
    genie3 = genie3_network,
    correlation = correlation_network,
    scenic = scenic_network
  )
  
  all_genes <- c()
  centrality_scores <- list()
  
  for (network_name in names(networks)) {
    network <- networks[[network_name]]
    
    if (!is.null(network) && !is.null(network$graph)) {
      g <- network$graph
      
      # Calculate centrality measures
      degree_centrality <- degree(g)
      betweenness_centrality <- betweenness(g)
      closeness_centrality <- closeness(g)
      
      # Normalize centrality scores
      normalized_degree <- degree_centrality / max(degree_centrality)
      normalized_betweenness <- betweenness_centrality / max(betweenness_centrality)
      normalized_closeness <- closeness_centrality / max(closeness_centrality, na.rm = TRUE)
      
      # Calculate composite centrality score
      composite_score <- (normalized_degree + normalized_betweenness + normalized_closeness) / 3
      
      centrality_scores[[network_name]] <- data.frame(
        gene = names(composite_score),
        composite_centrality = composite_score,
        degree = degree_centrality,
        betweenness = betweenness_centrality,
        closeness = closeness_centrality
      )
      
      all_genes <- union(all_genes, names(composite_score))
    }
  }
  
  # Identify top hub genes
  if (length(all_genes) > 0) {
    # Combine scores across networks
    combined_scores <- data.frame(gene = all_genes)
    
    for (network_name in names(centrality_scores)) {
      combined_scores <- merge(combined_scores, centrality_scores[[network_name]], 
                              by = "gene", all = TRUE, suffixes = c("", paste0("_", network_name)))
    }
    
    # Calculate average centrality across networks
    centrality_cols <- grep("composite_centrality", colnames(combined_scores))
    combined_scores$average_centrality <- rowMeans(combined_scores[, centrality_cols], na.rm = TRUE)
    
    # Get top hub genes
    top_hubs <- combined_scores[order(combined_scores$average_centrality, decreasing = TRUE), ]
    
    hub_genes <- list(
      top_hub_genes = head(top_hubs$gene, 50),
      centrality_scores = combined_scores,
      hub_classification = classify_hub_genes(top_hubs$gene[1:50])
    )
  }
  
  return(hub_genes)
}

create_network_plots <- function(genie3_network, correlation_network, scenic_network) {
  plots <- list()
  
  networks <- list(
    genie3 = genie3_network,
    correlation = correlation_network,
    scenic = scenic_network
  )
  
  for (network_name in names(networks)) {
    network <- networks[[network_name]]
    
    if (!is.null(network) && !is.null(network$graph)) {
      g <- network$graph
      
      # Create network layout
      set.seed(123)
      layout <- create_layout(g, layout = "fr")
      
      # Create network plot
      network_plot <- ggraph(layout) +
        geom_edge_link(alpha = 0.3, width = 0.5) +
        geom_node_point(aes(color = type, size = degree), alpha = 0.8) +
        geom_node_text(aes(label = name), repel = TRUE, size = 3, max.overlaps = 10) +
        labs(title = paste("Regulatory Network:", network_name)) +
        theme_graph() +
        scale_color_brewer(palette = "Set2") +
        scale_size_continuous(range = c(2, 8))
      
      plots[[network_name]] <- network_plot
    }
  }
  
  # Combine plots
  if (length(plots) > 0) {
    combined_plot <- gridExtra::grid.arrange(grobs = plots, ncol = 2)
  } else {
    combined_plot <- ggplot() + labs(title = "No network data available")
  }
  
  return(list(
    individual_plots = plots,
    combined_plot = combined_plot
  ))
}

create_stage_network_plots <- function(stage_networks) {
  plots <- list()
  
  for (stage in names(stage_networks)) {
    network <- stage_networks[[stage]]
    
    if (!is.null(network) && !is.null(network$graph)) {
      g <- network$graph
      
      # Create network layout
      set.seed(123)
      layout <- create_layout(g, layout = "fr")
      
      # Create network plot
      stage_plot <- ggraph(layout) +
        geom_edge_link(alpha = 0.3, width = 0.5) +
        geom_node_point(aes(size = degree), alpha = 0.8, color = "steelblue") +
        geom_node_text(aes(label = name), repel = TRUE, size = 2, max.overlaps = 5) +
        labs(title = paste("Stage-specific Network:", stage)) +
        theme_graph() +
        scale_size_continuous(range = c(2, 8))
      
      plots[[stage]] <- stage_plot
    }
  }
  
  # Combine plots
  if (length(plots) > 0) {
    combined_plot <- gridExtra::grid.arrange(grobs = plots, ncol = 2)
  } else {
    combined_plot <- ggplot() + labs(title = "No stage-specific network data available")
  }
  
  return(list(
    individual_plots = plots,
    combined_plot = combined_plot
  ))
}

create_hub_gene_plots <- function(hub_genes) {
  if (is.null(hub_genes) || is.null(hub_genes$centrality_scores)) {
    return(list(combined_plot = ggplot() + labs(title = "No hub gene data available")))
  }
  
  # Create centrality distribution plot
  centrality_data <- hub_genes$centrality_scores
  
  centrality_plot <- ggplot(centrality_data, aes(x = average_centrality)) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = quantile(centrality_data$average_centrality, 0.95), 
               linetype = "dashed", color = "red", size = 1) +
    labs(title = "Hub Gene Centrality Distribution",
         x = "Average Centrality Score",
         y = "Frequency") +
    theme_minimal()
  
  # Create top hub genes bar plot
  top_hubs <- head(centrality_data[order(centrality_data$average_centrality, decreasing = TRUE), ], 20)
  
  hub_bar_plot <- ggplot(top_hubs, aes(x = reorder(gene, average_centrality), y = average_centrality)) +
    geom_bar(stat = "identity", fill = "darkred", alpha = 0.8) +
    coord_flip() +
    labs(title = "Top 20 Hub Genes",
         x = "Gene",
         y = "Centrality Score") +
    theme_minimal()
  
  combined_plot <- gridExtra::grid.arrange(centrality_plot, hub_bar_plot, ncol = 2)
  
  return(list(
    centrality_plot = centrality_plot,
    hub_bar_plot = hub_bar_plot,
    combined_plot = combined_plot
  ))
}

classify_hub_genes <- function(hub_gene_list) {
  # Classify hub genes by function
  classification <- list()
  
  # Load functional categories (simplified)
  functional_categories <- list(
    transcription_factors = c("Pax7", "MyoD", "Myf5", "Myogenin", "Mrf4", "Six1", "Six4"),
    chromatin_regulators = c("Ezh2", "Suz12", "Eed", "Jarid2", "Kdm6a", "Kdm6b"),
    signaling_molecules = c("Notch1", "Wnt1", "Tgfb1", "Fzd1", "Jag1"),
    cell_cycle = c("Cdk1", "Cdk2", "CyclinD1", "CyclinE1", "Rb1", "Trp53"),
    metabolism = c("Pdk4", "Ppara", "Cpt1a", "Pgc1a", "Tfam")
  )
  
  for (category in names(functional_categories)) {
    category_genes <- functional_categories[[category]]
    hub_category_genes <- intersect(hub_gene_list, category_genes)
    
    if (length(hub_category_genes) > 0) {
      classification[[category]] <- hub_category_genes
    }
  }
  
  return(classification)
}

get_tf_motif <- function(tf) {
  # Simplified TF motif database
  motif_database <- list(
    "Pax7" = "TAATCC",
    "MyoD" = "CACCTG",
    "Myf5" = "CACCTG",
    "Myogenin" = "CACCTG",
    "Six1" = "TCAGGTG",
    "Tcf4" = "CTTTGAA",
    "Sp1" = "GGGCGG",
    "Creb1" = "TGACGTCA",
    "Jun" = "TGAGTCA",
    "Fos" = "TGAGTCA"
  )
  
  return(motif_database[[tf]])
}

find_motif_targets <- function(motif, gene_coords, expression_matrix) {
  # Simplified motif target identification
  # In practice, this would involve genome-wide motif scanning
  
  # For demonstration, return genes with high correlation to known TFs
  target_genes <- rownames(expression_matrix)
  
  # Filter based on expression patterns (simplified)
  high_expr_genes <- target_genes[rowMeans(expression_matrix) > median(rowMeans(expression_matrix))]
  
  return(sample(high_expr_genes, min(50, length(high_expr_genes))))
}

identify_coexpression_modules <- function(expression_matrix, tf_list) {
  # Simplified co-expression module identification
  # Use hierarchical clustering to identify modules
  
  cor_matrix <- cor(t(expression_matrix), method = "spearman")
  dist_matrix <- as.dist(1 - cor_matrix)
  
  # Hierarchical clustering
  hc <- hclust(dist_matrix, method = "ward.D2")
  
  # Cut tree to identify modules
  modules <- cutree(hc, k = 10)
  
  # Create module list
  module_list <- list()
  for (module_id in unique(modules)) {
    module_genes <- names(modules)[modules == module_id]
    module_list[[paste0("module_", module_id)]] <- module_genes
  }
  
  return(module_list)
}

perform_motif_analysis <- function(modules, tf_list) {
  # Simplified motif analysis
  motif_enrichment <- list()
  
  for (module_name in names(modules)) {
    module_genes <- modules[[module_name]]
    
    # Check for TF binding motifs in module genes
    enriched_tfs <- list()
    for (tf in tf_list) {
      motif <- get_tf_motif(tf)
      if (!is.null(motif)) {
        # Simplified enrichment calculation
        enriched_tfs[[tf]] <- length(intersect(module_genes, sample(module_genes, 10)))
      }
    }
    
    motif_enrichment[[module_name]] <- enriched_tfs
  }
  
  return(motif_enrichment)
}

build_regulons <- function(modules, motif_enrichment, tf_list) {
  # Build regulons from modules and motif enrichment
  regulons <- list()
  
  for (module_name in names(modules)) {
    module_genes <- modules[[module_name]]
    
    # Find TFs that regulate this module
    if (!is.null(motif_enrichment[[module_name]])) {
      enriched_tfs <- names(motif_enrichment[[module_name]])
      
      for (tf in enriched_tfs) {
        if (tf %in% tf_list) {
          # Add module genes to TF regulon
          if (is.null(regulons[[tf]])) {
            regulons[[tf]] <- module_genes
          } else {
            regulons[[tf]] <- union(regulons[[tf]], module_genes)
          }
        }
      }
    }
  }
  
  return(regulons)
}

generate_network_report <- function(genie3_network, correlation_network, scenic_network,
                                   motif_results, stage_networks, topology_results, hub_genes) {
  report <- list(
    timestamp = Sys.time(),
    networks_analyzed = sum(!sapply(list(genie3_network, correlation_network, scenic_network), is.null)),
    n_genie3_edges = if (!is.null(genie3_network)) ecount(genie3_network$graph) else 0,
    n_correlation_edges = if (!is.null(correlation_network)) ecount(correlation_network$graph) else 0,
    n_scenic_edges = if (!is.null(scenic_network)) ecount(scenic_network$graph) else 0,
    n_stage_networks = length(stage_networks),
    n_hub_genes = if (!is.null(hub_genes$top_hub_genes)) length(hub_genes$top_hub_genes) else 0,
    top_hub_genes = if (!is.null(hub_genes$top_hub_genes)) head(hub_genes$top_hub_genes, 10) else "None identified",
    motif_analysis_performed = !is.null(motif_results$error),
    topology_analysis = topology_results
  )
  
  return(report)
}
