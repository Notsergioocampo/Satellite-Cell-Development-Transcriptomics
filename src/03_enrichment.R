#!/usr/bin/env Rscript
# Satellite Cell Development Transcriptomics Pipeline
# Script 03: Pathway Enrichment and Hallmark Analysis
# MIT-Level Research Pipeline

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(DOSE)
  library(enrichplot)
  library(pathview)
  library(GSVA)
  library(GSEABase)
  library(yaml)
  library(logging)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(ggforce)
  library(ggraph)
  library(igraph)
})

# Configure logging
logfile <- file.path("results", "logs", "03_enrichment.log")
if (!dir.exists(dirname(logfile))) dir.create(dirname(logfile), recursive = TRUE)
basicConfig(level = "INFO", file = logfile)

loginfo("Starting pathway enrichment analysis pipeline")

# Load configuration
config <- yaml::read_yaml("config.yaml")

# Load DE results
loginfo("Loading differential expression results...")
de_file <- file.path("results", "de_analysis_results.RData")
if (!file.exists(de_file)) {
  stop("DE analysis results not found. Run 02_de_analysis.R first.")
}
load(de_file)

# Load preprocessed data
loginfo("Loading preprocessed expression data...")
processed_file <- file.path("data/processed", "preprocessed_expression.RData")
load(processed_file)

# Step 1: Comprehensive Gene Set Enrichment Analysis
loginfo("Performing comprehensive gene set enrichment analysis...")

# Perform GO enrichment for each contrast
go_enrichment_results <- perform_go_enrichment(de_results, config$analysis$padj_threshold)

# Perform KEGG pathway enrichment
kegg_enrichment_results <- perform_kegg_enrichment(de_results, config$analysis$padj_threshold)

# Perform Reactome pathway enrichment
reactome_enrichment_results <- perform_reactome_enrichment(de_results, config$analysis$padj_threshold)

# Perform MSigDB hallmark enrichment
hallmark_enrichment_results <- perform_msigdb_hallmark_enrichment(de_results, config$analysis$padj_threshold)

# Step 2: Hallmark Pathway Velocity Analysis
loginfo("Performing hallmark pathway velocity analysis...")
pathway_velocity_results <- perform_pathway_velocity_analysis(final_expression, sample_info, config$analysis$padj_threshold)

# Step 3: Gene Set Variation Analysis (GSVA)
loginfo("Performing gene set variation analysis...")
gsva_results <- perform_gsva_analysis(final_expression, sample_info)

# Step 4: Pathway Dynamics Over Developmental Time
loginfo("Analyzing pathway dynamics over developmental time...")
pathway_dynamics <- analyze_pathway_dynamics(final_expression, sample_info, hallmark_enrichment_results)

# Step 5: Quiescence Signature Discovery
loginfo("Discovering quiescence-specific gene signatures...")
quiescence_signatures <- discover_quiescence_signatures(final_expression, sample_info, de_results)

# Step 6: Create Enrichment Visualizations
loginfo("Creating enrichment visualizations...")

# Create enrichment dot plots
enrichment_plots <- create_enrichment_plots(go_enrichment_results, kegg_enrichment_results, 
                                           reactome_enrichment_results, hallmark_enrichment_results)

# Create pathway velocity plots
velocity_plots <- create_pathway_velocity_plots(pathway_velocity_results)

# Create GSVA heatmaps
gsva_heatmaps <- create_gsva_heatmaps(gsva_results)

# Create pathway network
pathway_network <- create_pathway_network(hallmark_enrichment_results)

# Step 7: Generate Comprehensive Enrichment Report
loginfo("Generating comprehensive enrichment report...")

# Save enrichment results
enrichment_file <- file.path("results", "enrichment_analysis_results.RData")
save(go_enrichment_results, kegg_enrichment_results, reactome_enrichment_results,
     hallmark_enrichment_results, pathway_velocity_results, gsva_results, 
     pathway_dynamics, quiescence_signatures, file = enrichment_file)

# Save plots
enrichment_plot_file <- file.path("figures", "enrichment_analysis.pdf")
ggsave(enrichment_plot_file, enrichment_plots$combined_plot, width = 16, height = 12, 
       dpi = config$visualization$figure_dpi)

velocity_plot_file <- file.path("figures", "pathway_velocity.pdf")
ggsave(velocity_plot_file, velocity_plots$combined_plot, width = 14, height = 10, 
       dpi = config$visualization$figure_dpi)

gsva_heatmap_file <- file.path("figures", "gsva_heatmaps.pdf")
ggsave(gsva_heatmap_file, gsva_heatmaps$combined_plot, width = 12, height = 10, 
       dpi = config$visualization$figure_dpi)

# Generate enrichment report
enrichment_report <- generate_enrichment_report(go_enrichment_results, kegg_enrichment_results,
                                               reactome_enrichment_results, hallmark_enrichment_results,
                                               pathway_velocity_results, quiescence_signatures)
report_file <- file.path("results", "enrichment_analysis_report.html")
rmarkdown::render(
  input = system.file("rmarkdown", "templates", "enrichment_report", package = "rmarkdown"),
  output_file = report_file,
  params = list(report_data = enrichment_report)
)

loginfo("Pathway enrichment analysis pipeline completed successfully")

# Helper functions
perform_go_enrichment <- function(de_results, padj_threshold) {
  go_results <- list()
  
  for (contrast_name in names(de_results)) {
    res <- de_results[[contrast_name]]$full_results
    
    # Get significant genes
    sig_genes <- rownames(res)[which(res$padj < padj_threshold & abs(res$log2FoldChange) > 1)]
    
    if (length(sig_genes) > 0) {
      # Convert to Entrez IDs
      entrez_ids <- map_gene_symbols_to_entrez(sig_genes)
      
      # Perform GO enrichment
      ego <- enrichGO(
        gene = entrez_ids,
        OrgDb = org.Mm.eg.db,
        ont = "ALL",
        pAdjustMethod = "BH",
        pvalueCutoff = padj_threshold,
        qvalueCutoff = padj_threshold,
        minGSSize = config$analysis$min_geneset_size,
        maxGSSize = config$analysis$max_geneset_size
      )
      
      go_results[[contrast_name]] <- ego
    }
  }
  
  return(go_results)
}

perform_kegg_enrichment <- function(de_results, padj_threshold) {
  kegg_results <- list()
  
  for (contrast_name in names(de_results)) {
    res <- de_results[[contrast_name]]$full_results
    
    # Get significant genes
    sig_genes <- rownames(res)[which(res$padj < padj_threshold & abs(res$log2FoldChange) > 1)]
    
    if (length(sig_genes) > 0) {
      # Convert to Entrez IDs
      entrez_ids <- map_gene_symbols_to_entrez(sig_genes)
      
      # Perform KEGG enrichment
      kke <- enrichKEGG(
        gene = entrez_ids,
        organism = "mmu",
        pvalueCutoff = padj_threshold,
        qvalueCutoff = padj_threshold,
        minGSSize = config$analysis$min_geneset_size,
        maxGSSize = config$analysis$max_geneset_size
      )
      
      kegg_results[[contrast_name]] <- kke
    }
  }
  
  return(kegg_results)
}

perform_reactome_enrichment <- function(de_results, padj_threshold) {
  reactome_results <- list()
  
  # Load Reactome pathways (if available)
  tryCatch({
    library(reactome.db)
    
    for (contrast_name in names(de_results)) {
      res <- de_results[[contrast_name]]$full_results
      
      # Get significant genes
      sig_genes <- rownames(res)[which(res$padj < padj_threshold & abs(res$log2FoldChange) > 1)]
      
      if (length(sig_genes) > 0) {
        # Convert to Entrez IDs
        entrez_ids <- map_gene_symbols_to_entrez(sig_genes)
        
        # Perform Reactome enrichment
        reactome_enrich <- enrichPathway(
          gene = entrez_ids,
          pvalueCutoff = padj_threshold,
          qvalueCutoff = padj_threshold
        )
        
        reactome_results[[contrast_name]] <- reactome_enrich
      }
    }
    
  }, error = function(e) {
    logwarn(paste("Reactome enrichment failed:", e$message))
    reactome_results <- list(error = "Reactome database not available")
  })
  
  return(reactome_results)
}

perform_msigdb_hallmark_enrichment <- function(de_results, padj_threshold) {
  hallmark_results <- list()
  
  # Load MSigDB hallmark gene sets
  tryCatch({
    library(msigdbr)
    
    # Get mouse hallmark gene sets
    m_df <- msigdbr(species = "Mus musculus", category = "H")
    
    for (contrast_name in names(de_results)) {
      res <- de_results[[contrast_name]]$full_results
      
      # Get significant genes
      sig_genes <- rownames(res)[which(res$padj < padj_threshold & abs(res$log2FoldChange) > 1)]
      
      if (length(sig_genes) > 0) {
        # Convert to Entrez IDs
        entrez_ids <- map_gene_symbols_to_entrez(sig_genes)
        
        # Perform hallmark enrichment
        hallmark_enrich <- enricher(
          gene = entrez_ids,
          TERM2GENE = m_df[, c("gs_name", "gene_symbol")],
          pAdjustMethod = "BH",
          pvalueCutoff = padj_threshold,
          qvalueCutoff = padj_threshold
        )
        
        hallmark_results[[contrast_name]] <- hallmark_enrich
      }
    }
    
  }, error = function(e) {
    logwarn(paste("MSigDB hallmark enrichment failed:", e$message))
    hallmark_results <- list(error = "MSigDB not available")
  })
  
  return(hallmark_results)
}

perform_pathway_velocity_analysis <- function(expression_matrix, sample_info, padj_threshold) {
  # Create time points
  time_points <- ifelse(sample_info$stage == "P1", 1, 
                       ifelse(sample_info$stage == "P12", 12, 28))
  
  # Load hallmark gene sets
  tryCatch({
    library(msigdbr)
    m_df <- msigdbr(species = "Mus musculus", category = "H")
    
    # Calculate pathway activity scores using ssGSEA
    pathway_scores <- list()
    
    for (pathway in unique(m_df$gs_name)) {
      pathway_genes <- m_df$gene_symbol[m_df$gs_name == pathway]
      common_genes <- intersect(pathway_genes, rownames(expression_matrix))
      
      if (length(common_genes) > 10) {
        # Calculate single-sample GSEA score
        pathway_expr <- expression_matrix[common_genes, ]
        
        # Simple ssGSEA implementation
        ranks <- apply(pathway_expr, 2, function(x) rank(x))
        pathway_score <- colMeans(ranks[common_genes, , drop = FALSE])
        
        pathway_scores[[pathway]] <- pathway_score
      }
    }
    
    # Convert to matrix
    pathway_matrix <- do.call(rbind, pathway_scores)
    
    # Calculate velocity (rate of change) for each pathway
    velocity_results <- list()
    
    for (pathway in rownames(pathway_matrix)) {
      scores <- pathway_matrix[pathway, ]
      
      # Fit polynomial model to capture dynamics
      model <- lm(scores ~ poly(time_points, 2))
      
      # Calculate first derivative (velocity) at each time point
      velocity <- predict(model, newdata = data.frame(time_points = time_points), 
                         type = "terms")[, "poly(time_points, 2)1"]
      
      # Calculate acceleration (second derivative)
      acceleration <- predict(model, newdata = data.frame(time_points = time_points), 
                             type = "terms")[, "poly(time_points, 2)2"]
      
      velocity_results[[pathway]] <- list(
        pathway = pathway,
        scores = scores,
        velocity = velocity,
        acceleration = acceleration,
        time_points = time_points,
        model = model
      )
    }
    
    return(velocity_results)
    
  }, error = function(e) {
    logwarn(paste("Pathway velocity analysis failed:", e$message))
    return(list(error = e$message))
  })
}

perform_gsva_analysis <- function(expression_matrix, sample_info) {
  # Load gene sets
  tryCatch({
    library(msigdbr)
    m_df <- msigdbr(species = "Mus musculus", category = "H")
    
    # Create gene set collection
    gene_sets <- list()
    for (pathway in unique(m_df$gs_name)) {
      pathway_genes <- m_df$gene_symbol[m_df$gs_name == pathway]
      common_genes <- intersect(pathway_genes, rownames(expression_matrix))
      if (length(common_genes) > 10) {
        gene_sets[[pathway]] <- common_genes
      }
    }
    
    # Perform GSVA
    gsva_result <- gsva(as.matrix(expression_matrix), 
                       gene_sets, 
                       method = "gsva",
                       kcdf = "Gaussian")
    
    return(list(
      gsva_scores = gsva_result,
      gene_sets = gene_sets,
      sample_info = sample_info
    ))
    
  }, error = function(e) {
    logwarn(paste("GSVA analysis failed:", e$message))
    return(list(error = e$message))
  })
}

analyze_pathway_dynamics <- function(expression_matrix, sample_info, hallmark_results) {
  # Analyze how pathways change over developmental time
  dynamics_results <- list()
  
  # For each contrast, analyze pathway changes
  for (contrast_name in names(hallmark_results)) {
    if (!is.null(hallmark_results[[contrast_name]]) && 
        !is.character(hallmark_results[[contrast_name]])) {
      
      enrich_result <- hallmark_results[[contrast_name]]
      
      # Get significant pathways
      sig_pathways <- enrich_result@result[enrich_result@result$p.adjust < 0.05, ]
      
      if (nrow(sig_pathways) > 0) {
        # Analyze direction of change
        up_pathways <- sig_pathways[sig_pathways$GeneRatio > 0, ]
        down_pathways <- sig_pathways[sig_pathways$GeneRatio < 0, ]
        
        dynamics_results[[contrast_name]] <- list(
          significant_pathways = sig_pathways,
          upregulated_pathways = up_pathways,
          downregulated_pathways = down_pathways,
          n_total = nrow(sig_pathways),
          n_up = nrow(up_pathways),
          n_down = nrow(down_pathways)
        )
      }
    }
  }
  
  return(dynamics_results)
}

discover_quiescence_signatures <- function(expression_matrix, sample_info, de_results) {
  # Identify genes associated with quiescence (P28 vs P12 comparison)
  quiescence_results <- list()
  
  if ("P28_vs_P12" %in% names(de_results)) {
    p28_vs_p12 <- de_results[["P28_vs_P12"]]$full_results
    
    # Get genes upregulated in P28 (quiescent state)
    quiescence_genes <- rownames(p28_vs_p12)[which(p28_vs_p12$log2FoldChange > 1 & 
                                                   p28_vs_p12$padj < 0.05)]
    
    # Identify transcription factors
    tfs <- identify_transcription_factors(quiescence_genes)
    
    # Identify chromatin regulators
    chromatin_regs <- identify_chromatin_regulators(quiescence_genes)
    
    # Identify metabolic genes
    metabolic_genes <- identify_metabolic_genes(quiescence_genes)
    
    # Identify signaling pathway genes
    signaling_genes <- identify_signaling_genes(quiescence_genes)
    
    quiescence_results <- list(
      quiescence_genes = quiescence_genes,
      n_quiescence_genes = length(quiescence_genes),
      transcription_factors = tfs,
      chromatin_regulators = chromatin_regs,
      metabolic_genes = metabolic_genes,
      signaling_genes = signaling_genes,
      top_candidates = head(quiescence_genes, 20)
    )
  }
  
  return(quiescence_results)
}

create_enrichment_plots <- function(go_results, kegg_results, reactome_results, hallmark_results) {
  plots <- list()
  
  # Create dot plots for each enrichment type
  for (contrast_name in names(go_results)) {
    if (!is.null(go_results[[contrast_name]])) {
      go_plot <- dotplot(go_results[[contrast_name]], showCategory = 10, title = paste("GO Enrichment:", contrast_name))
      plots[[paste0("go_", contrast_name)]] <- go_plot
    }
  }
  
  for (contrast_name in names(kegg_results)) {
    if (!is.null(kegg_results[[contrast_name]])) {
      kegg_plot <- dotplot(kegg_results[[contrast_name]], showCategory = 10, title = paste("KEGG Enrichment:", contrast_name))
      plots[[paste0("kegg_", contrast_name)]] <- kegg_plot
    }
  }
  
  for (contrast_name in names(hallmark_results)) {
    if (!is.null(hallmark_results[[contrast_name]]) && !is.character(hallmark_results[[contrast_name]])) {
      hallmark_plot <- dotplot(hallmark_results[[contrast_name]], showCategory = 10, title = paste("Hallmark Enrichment:", contrast_name))
      plots[[paste0("hallmark_", contrast_name)]] <- hallmark_plot
    }
  }
  
  # Combine plots
  if (length(plots) > 0) {
    combined_plot <- gridExtra::grid.arrange(grobs = plots, ncol = 2)
  } else {
    combined_plot <- ggplot() + labs(title = "No enrichment results available")
  }
  
  return(list(
    individual_plots = plots,
    combined_plot = combined_plot
  ))
}

create_pathway_velocity_plots <- function(velocity_results) {
  if (is.null(velocity_results) || is.character(velocity_results)) {
    return(list(combined_plot = ggplot() + labs(title = "Pathway velocity analysis failed")))
  }
  
  plots <- list()
  
  # Select top pathways with highest velocity
  top_pathways <- head(names(velocity_results), 12)
  
  for (pathway in top_pathways) {
    if (!is.null(velocity_results[[pathway]])) {
      result <- velocity_results[[pathway]]
      
      # Create velocity plot
      velocity_df <- data.frame(
        time_point = result$time_points,
        pathway_score = result$scores,
        velocity = result$velocity,
        acceleration = result$acceleration
      )
      
      velocity_plot <- ggplot(velocity_df, aes(x = time_point)) +
        geom_line(aes(y = pathway_score, color = "Pathway Score"), size = 1.2) +
        geom_line(aes(y = velocity * 10, color = "Velocity (×10)"), size = 1, linetype = "dashed") +
        labs(title = pathway,
             x = "Postnatal Day",
             y = "Score / Velocity") +
        theme_minimal() +
        scale_x_continuous(breaks = c(1, 12, 28), labels = c("P1", "P12", "P28")) +
        scale_color_brewer(palette = "Set1")
      
      plots[[pathway]] <- velocity_plot
    }
  }
  
  # Combine plots
  if (length(plots) > 0) {
    combined_plot <- gridExtra::grid.arrange(grobs = plots, ncol = 3)
  } else {
    combined_plot <- ggplot() + labs(title = "No velocity results available")
  }
  
  return(list(
    individual_plots = plots,
    combined_plot = combined_plot
  ))
}

create_gsva_heatmaps <- function(gsva_results) {
  if (is.null(gsva_results$gsva_scores) || is.character(gsva_results)) {
    return(list(combined_plot = ggplot() + labs(title = "GSVA analysis failed")))
  }
  
  # Get GSVA scores
  gsva_scores <- gsva_results$gsva_scores
  sample_info <- gsva_results$sample_info
  
  # Order samples by developmental stage
  stage_order <- order(factor(sample_info$stage, levels = c("P1", "P12", "P28")))
  gsva_ordered <- gsva_scores[, stage_order]
  sample_info_ordered <- sample_info[stage_order, ]
  
  # Create heatmap
  heatmap_plot <- Heatmap(gsva_ordered,
                         name = "GSVA Score",
                         col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
                         show_row_names = TRUE,
                         show_column_names = FALSE,
                         cluster_rows = TRUE,
                         cluster_cols = FALSE,
                         top_annotation = HeatmapAnnotation(
                           stage = sample_info_ordered$stage,
                           col = list(stage = c("P1" = "#E41A1C", "P12" = "#377EB8", "P28" = "#4DAF4A"))
                         ))
  
  return(list(
    heatmap_plot = heatmap_plot,
    gsva_scores = gsva_scores
  ))
}

create_pathway_network <- function(hallmark_results) {
  # Create network of enriched pathways
  network_data <- list()
  
  for (contrast_name in names(hallmark_results)) {
    if (!is.null(hallmark_results[[contrast_name]]) && !is.character(hallmark_results[[contrast_name]])) {
      enrich_result <- hallmark_results[[contrast_name]]
      
      # Get significant pathways
      sig_pathways <- enrich_result@result[enrich_result@result$p.adjust < 0.05, ]
      
      if (nrow(sig_pathways) > 0) {
        network_data[[contrast_name]] <- sig_pathways
      }
    }
  }
  
  # Create simple network visualization
  if (length(network_data) > 0) {
    # Extract pathway names and create edges based on shared genes
    all_pathways <- unique(unlist(lapply(network_data, function(x) x$Description)))
    
    # Create network plot
    network_plot <- ggplot(data.frame(pathway = all_pathways), aes(x = 1, y = seq_along(pathway))) +
      geom_point(size = 5, color = "steelblue") +
      geom_text(aes(label = pathway), hjust = -0.1, size = 3) +
      labs(title = "Hallmark Pathway Network",
           x = "", y = "") +
      theme_minimal() +
      theme(axis.text = element_blank(), axis.ticks = element_blank())
    
    return(list(
      network_plot = network_plot,
      network_data = network_data
    ))
  } else {
    return(list(
      network_plot = ggplot() + labs(title = "No pathway network data available"),
      network_data = NULL
    ))
  }
}

map_gene_symbols_to_entrez <- function(gene_symbols) {
  # Map gene symbols to Entrez IDs
  tryCatch({
    entrez_ids <- mapIds(org.Mm.eg.db, keys = gene_symbols, keytype = "SYMBOL", 
                        column = "ENTREZID")
    return(na.omit(entrez_ids))
  }, error = function(e) {
    logwarn(paste("Gene mapping failed:", e$message))
    return(gene_symbols)  # Return original symbols if mapping fails
  })
}

identify_transcription_factors <- function(gene_list) {
  # Load transcription factor database (simplified)
  tf_database <- c("Pax7", "MyoD", "Myf5", "Myogenin", "Mrf4", "Six1", "Six4", 
                   "Eya1", "Eya2", "Tcf4", "Sp1", "Creb1", "Atf1", "Jun", "Fos")
  
  tfs <- intersect(gene_list, tf_database)
  return(tfs)
}

identify_chromatin_regulators <- function(gene_list) {
  # Load chromatin regulator database (simplified)
  chromatin_database <- c("Ezh2", "Suz12", "Eed", "Jarid2", "Mtf2", "Rnf2", "Rybp",
                         "Kdm6a", "Kdm6b", "Kdm4a", "Kdm4b", "Hdac1", "Hdac2", "Hdac3")
  
  chromatin_regs <- intersect(gene_list, chromatin_database)
  return(chromatin_regs)
}

identify_metabolic_genes <- function(gene_list) {
  # Load metabolic gene database (simplified)
  metabolic_database <- c("Pdk4", "Ppara", "Cpt1a", "Cpt1b", "Acadm", "Acadl", "Pgc1a",
                         "Tfam", "Nrf1", "Nrf2", "Sod1", "Sod2", "Cat", "Gpx1")
  
  metabolic_genes <- intersect(gene_list, metabolic_database)
  return(metabolic_genes)
}

identify_signaling_genes <- function(gene_list) {
  # Load signaling pathway gene database (simplified)
  signaling_database <- c("Notch1", "Notch2", "Notch3", "Jag1", "Jag2", "Dll1", "Dll4",
                         "Wnt1", "Wnt3a", "Wnt5a", "Wnt7a", "Fzd1", "Fzd2", "Fzd3",
                         "Tgfb1", "Tgfb2", "Tgfb3", "Tgfbr1", "Tgfbr2", "Smad2", "Smad3")
  
  signaling_genes <- intersect(gene_list, signaling_database)
  return(signaling_genes)
}

generate_enrichment_report <- function(go_results, kegg_results, reactome_results, 
                                     hallmark_results, velocity_results, quiescence_results) {
  report <- list(
    timestamp = Sys.time(),
    go_enrichments = length(go_results),
    kegg_enrichments = length(kegg_results),
    reactome_enrichments = length(reactome_results),
    hallmark_enrichments = length(hallmark_results),
    pathway_velocity_performed = !is.null(velocity_results$error),
    quiescence_signatures_found = if (!is.null(quiescence_results$n_quiescence_genes)) {
      quiescence_results$n_quiescence_genes
    } else {
      0
    },
    top_quiescence_tfs = if (!is.null(quiescence_results$transcription_factors)) {
      quiescence_results$transcription_factors
    } else {
      "None identified"
    }
  )
  
  return(report)
}
