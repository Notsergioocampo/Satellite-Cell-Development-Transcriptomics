#!/usr/bin/env Rscript
# Satellite Cell Development Transcriptomics Pipeline
# Interactive Shiny Dashboard
# MIT-Level Research Pipeline

library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(plotly)
library(DT)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(ComplexHeatmap)
library(circlize)

# Load analysis results
load_results <- function() {
  results <- list()
  
  # Load preprocessed data
  if (file.exists("data/processed/preprocessed_expression.RData")) {
    load("data/processed/preprocessed_expression.RData")
    results$expression <- final_expression
    results$sample_info <- sample_info
  }
  
  # Load DE results
  if (file.exists("results/de_analysis_results.RData")) {
    load("results/de_analysis_results.RData")
    results$de_results <- de_results
  }
  
  # Load enrichment results
  if (file.exists("results/enrichment_analysis_results.RData")) {
    load("results/enrichment_analysis_results.RData")
    results$enrichment <- list(
      go = go_enrichment_results,
      kegg = kegg_enrichment_results,
      hallmark = hallmark_enrichment_results
    )
  }
  
  # Load network results
  if (file.exists("results/network_inference_results.RData")) {
    load("results/network_inference_results.RData")
    results$network <- list(
      genie3 = genie3_network,
      correlation = correlation_network,
      scenic = scenic_network
    )
  }
  
  # Load visualization results
  if (file.exists("results/interactive_visualization_data.RData")) {
    load("results/interactive_visualization_data.RData")
    results$interactive <- interactive_data
  }
  
  return(results)
}

# UI Definition
ui <- dashboardPage(
  dashboardHeader(
    title = "Satellite Cell Development Explorer",
    titleWidth = 350
  ),
  
  dashboardSidebar(
    width = 350,
    
    sidebarMenu(
      menuItem("Overview", tabName = "overview", icon = icon("home")),
      menuItem("Expression Explorer", tabName = "expression", icon = icon("chart-line")),
      menuItem("Differential Expression", tabName = "de", icon = icon("not-equal")),
      menuItem("Pathway Enrichment", tabName = "enrichment", icon = icon("project-diagram")),
      menuItem("Regulatory Networks", tabName = "networks", icon = icon("network-wired")),
      menuItem("Trajectory Analysis", tabName = "trajectory", icon = icon("route")),
      menuItem("Cell State Discovery", tabName = "cellstates", icon = icon("microscope")),
      menuItem("Download Results", tabName = "download", icon = icon("download"))
    ),
    
    # Global controls
    hr(),
    h4("Global Settings", style = "padding-left: 15px;"),
    
    fluidRow(
      column(12,
        selectInput("color_palette", "Color Palette:",
                    choices = c("Viridis" = "viridis", "Set2" = "set2", "Dark2" = "dark2", "Paired" = "paired"),
                    selected = "viridis"),
        
        sliderInput("point_size", "Point Size:",
                    min = 1, max = 10, value = 3, step = 0.5),
        
        sliderInput("alpha_value", "Transparency:",
                    min = 0.1, max = 1, value = 0.8, step = 0.1),
        
        checkboxInput("show_labels", "Show Sample Labels", value = FALSE),
        
        actionButton("refresh_data", "Refresh Data", icon = icon("refresh"), 
                    class = "btn-primary btn-block")
      )
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper, .right-side {
          background-color: #f8f9fa;
        }
        .box {
          border-radius: 10px;
          box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }
        .info-box {
          border-radius: 10px;
        }
        .small-box {
          border-radius: 10px;
        }
      "))
    ),
    
    tabItems(
      # Overview Tab
      tabItem(tabName = "overview",
        fluidRow(
          valueBoxOutput("n_genes"),
          valueBoxOutput("n_samples"),
          valueBoxOutput("n_de_genes")
        ),
        
        fluidRow(
          box(
            width = 12,
            title = "Pipeline Overview",
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            collapsed = FALSE,
            tags$h4("Satellite Cell Development Transcriptomics Pipeline"),
            tags$p("This interactive dashboard explores the molecular mechanisms governing satellite cell development across postnatal stages (P1 → P12 → P28)."),
            tags$br(),
            tags$h5("Key Features:"),
            tags$ul(
              tags$li("Differential expression analysis across developmental stages"),
              tags$li("Pathway enrichment and hallmark analysis"),
              tags$li("Gene regulatory network inference"),
              tags$li("Manifold learning and trajectory analysis"),
              tags$li("Quiescence signature discovery")
            ),
            tags$br(),
            tags$h5("Dataset Information:"),
            tags$ul(
              tags$li("GEO Accession: GSE65927"),
              tags$li("Platform: Illumina MouseRef-8 v2.0"),
              tags$li("Species: Mus musculus"),
              tags$li("Stages: P1, P12, P28 (3 replicates each)")
            )
          )
        ),
        
        fluidRow(
          box(
            width = 6,
            title = "Quick Navigation",
            status = "info",
            solidHeader = TRUE,
            tags$h5("Analysis Modules:"),
            actionButton("go_expression", "Expression Explorer", icon = icon("chart-line"), 
                        width = "100%", class = "btn-info"),
            br(), br(),
            actionButton("go_de", "Differential Expression", icon = icon("not-equal"), 
                        width = "100%", class = "btn-success"),
            br(), br(),
            actionButton("go_enrichment", "Pathway Enrichment", icon = icon("project-diagram"), 
                        width = "100%", class = "btn-warning")
          ),
          
          box(
            width = 6,
            title = "Advanced Analysis",
            status = "danger",
            solidHeader = TRUE,
            tags$h5("Systems Biology:"),
            actionButton("go_networks", "Regulatory Networks", icon = icon("network-wired"), 
                        width = "100%", class = "btn-danger"),
            br(), br(),
            actionButton("go_trajectory", "Trajectory Analysis", icon = icon("route"), 
                        width = "100%", class = "btn-primary"),
            br(), br(),
            actionButton("go_cellstates", "Cell States", icon = icon("microscope"), 
                        width = "100%", class = "btn-dark")
          )
        )
      ),
      
      # Expression Explorer Tab
      tabItem(tabName = "expression",
        fluidRow(
          box(
            width = 4,
            title = "Gene Selection",
            status = "primary",
            solidHeader = TRUE,
            textInput("gene_input", "Enter Gene Symbol:", 
                     placeholder = "e.g., Pax7, MyoD, Myf5"),
            actionButton("search_gene", "Search", icon = icon("search")),
            br(), br(),
            selectInput("stage_filter", "Filter by Stage:",
                       choices = c("All", "P1", "P12", "P28"),
                       selected = "All"),
            br(),
            actionButton("show_top_genes", "Show Top Variable Genes", 
                        icon = icon("star"), class = "btn-success")
          ),
          
          box(
            width = 8,
            title = "Gene Expression Visualization",
            status = "primary",
            solidHeader = TRUE,
            plotlyOutput("gene_expression_plot"),
            br(),
            downloadButton("download_expression_plot", "Download Plot")
          )
        ),
        
        fluidRow(
          box(
            width = 12,
            title = "Expression Heatmap",
            status = "info",
            solidHeader = TRUE,
            plotOutput("expression_heatmap", height = "600px"),
            br(),
            downloadButton("download_heatmap", "Download Heatmap")
          )
        )
      ),
      
      # Differential Expression Tab
      tabItem(tabName = "de",
        fluidRow(
          box(
            width = 4,
            title = "Contrast Selection",
            status = "success",
            solidHeader = TRUE,
            selectInput("de_contrast", "Select Contrast:",
                       choices = c("P12 vs P1", "P28 vs P12", "P28 vs P1"),
                       selected = "P12 vs P1"),
            br(),
            sliderInput("fc_threshold", "Fold Change Threshold:",
                       min = 0.5, max = 3, value = 1, step = 0.1),
            br(),
            sliderInput("p_threshold", "P-value Threshold:",
                       min = 0.001, max = 0.1, value = 0.05, step = 0.001),
            br(),
            actionButton("update_de", "Update Analysis", icon = icon("refresh"), 
                        class = "btn-success")
          ),
          
          box(
            width = 8,
            title = "Volcano Plot",
            status = "success",
            solidHeader = TRUE,
            plotlyOutput("volcano_plot"),
            br(),
            downloadButton("download_volcano", "Download Volcano Plot")
          )
        ),
        
        fluidRow(
          box(
            width = 6,
            title = "Top Upregulated Genes",
            status = "info",
            solidHeader = TRUE,
            DTOutput("up_genes_table")
          ),
          
          box(
            width = 6,
            title = "Top Downregulated Genes",
            status = "warning",
            solidHeader = TRUE,
            DTOutput("down_genes_table")
          )
        )
      ),
      
      # Pathway Enrichment Tab
      tabItem(tabName = "enrichment",
        fluidRow(
          box(
            width = 4,
            title = "Enrichment Settings",
            status = "warning",
            solidHeader = TRUE,
            selectInput("enrichment_type", "Enrichment Type:",
                       choices = c("GO Biological Process", "KEGG Pathways", "Hallmark Gene Sets"),
                       selected = "Hallmark Gene Sets"),
            br(),
            selectInput("enrichment_contrast", "Contrast:",
                       choices = c("P12 vs P1", "P28 vs P12", "P28 vs P1"),
                       selected = "P28 vs P1"),
            br(),
            sliderInput("top_pathways", "Top Pathways:",
                       min = 5, max = 50, value = 20, step = 5),
            br(),
            actionButton("run_enrichment", "Run Enrichment", icon = icon("play"), 
                        class = "btn-warning")
          ),
          
          box(
            width = 8,
            title = "Enrichment Results",
            status = "warning",
            solidHeader = TRUE,
            plotlyOutput("enrichment_plot"),
            br(),
            downloadButton("download_enrichment", "Download Enrichment Plot")
          )
        ),
        
        fluidRow(
          box(
            width = 12,
            title = "Pathway Details",
            status = "info",
            solidHeader = TRUE,
            DTOutput("pathway_details_table")
          )
        )
      ),
      
      # Regulatory Networks Tab
      tabItem(tabName = "networks",
        fluidRow(
          box(
            width = 4,
            title = "Network Settings",
            status = "danger",
            solidHeader = TRUE,
            selectInput("network_type", "Network Type:",
                       choices = c("GENIE3", "Correlation", "SCENIC"),
                       selected = "GENIE3"),
            br(),
            selectInput("network_stage", "Developmental Stage:",
                       choices = c("All Stages", "P1", "P12", "P28"),
                       selected = "All Stages"),
            br(),
            sliderInput("edge_threshold", "Edge Weight Threshold:",
                       min = 0.1, max = 1, value = 0.5, step = 0.1),
            br(),
            sliderInput("max_nodes", "Max Nodes to Display:",
                       min = 10, max = 200, value = 50, step = 10),
            br(),
            actionButton("generate_network", "Generate Network", icon = icon("project-diagram"), 
                        class = "btn-danger")
          ),
          
          box(
            width = 8,
            title = "Gene Regulatory Network",
            status = "danger",
            solidHeader = TRUE,
            visNetworkOutput("network_plot"),
            br(),
            downloadButton("download_network", "Download Network")
          )
        ),
        
        fluidRow(
          box(
            width = 6,
            title = "Hub Genes",
            status = "info",
            solidHeader = TRUE,
            DTOutput("hub_genes_table")
          ),
          
          box(
            width = 6,
            title = "Network Statistics",
            status = "success",
            solidHeader = TRUE,
            tableOutput("network_stats")
          )
        )
      ),
      
      # Trajectory Analysis Tab
      tabItem(tabName = "trajectory",
        fluidRow(
          box(
            width = 4,
            title = "Trajectory Settings",
            status = "primary",
            solidHeader = TRUE,
            selectInput("trajectory_method", "Manifold Learning Method:",
                       choices = c("UMAP", "t-SNE", "PHATE", "Diffusion Map"),
                       selected = "UMAP"),
            br(),
            selectInput("color_by", "Color By:",
                       choices = c("Developmental Stage", "Pseudotime", "Cell State"),
                       selected = "Developmental Stage"),
            br(),
            actionButton("generate_trajectory", "Generate Trajectory", icon = icon("route"), 
                        class = "btn-primary")
          ),
          
          box(
            width = 8,
            title = "Developmental Trajectory",
            status = "primary",
            solidHeader = TRUE,
            plotlyOutput("trajectory_plot"),
            br(),
            downloadButton("download_trajectory", "Download Trajectory Plot")
          )
        ),
        
        fluidRow(
          box(
            width = 12,
            title = "Pseudotime Analysis",
            status = "info",
            solidHeader = TRUE,
            plotlyOutput("pseudotime_plot"),
            br(),
            downloadButton("download_pseudotime", "Download Pseudotime Plot")
          )
        )
      ),
      
      # Cell State Discovery Tab
      tabItem(tabName = "cellstates",
        fluidRow(
          box(
            width = 4,
            title = "Clustering Settings",
            status = "dark",
            solidHeader = TRUE,
            selectInput("clustering_method", "Clustering Method:",
                       choices = c("K-means", "Hierarchical", "DBSCAN", "Consensus"),
                       selected = "Consensus"),
            br(),
            sliderInput("n_clusters", "Number of Clusters:",
                       min = 2, max = 10, value = 5, step = 1),
            br(),
            actionButton("run_clustering", "Run Clustering", icon = icon("microscope"), 
                        class = "btn-dark")
          ),
          
          box(
            width = 8,
            title = "Cell State Visualization",
            status = "dark",
            solidHeader = TRUE,
            plotlyOutput("cell_state_plot"),
            br(),
            downloadButton("download_cell_states", "Download Cell States")
          )
        ),
        
        fluidRow(
          box(
            width = 6,
            title = "Cluster Composition",
            status = "info",
            solidHeader = TRUE,
            tableOutput("cluster_composition")
          ),
          
          box(
            width = 6,
            title = "Marker Genes",
            status = "warning",
            solidHeader = TRUE,
            DTOutput("marker_genes_table")
          )
        )
      ),
      
      # Download Results Tab
      tabItem(tabName = "download",
        fluidRow(
          box(
            width = 12,
            title = "Download Analysis Results",
            status = "success",
            solidHeader = TRUE,
            tags$h4("Available Downloads:"),
            br(),
            
            fluidRow(
              column(6,
                downloadButton("download_de_results", "Differential Expression Results", 
                              class = "btn-block btn-lg"),
                br(), br(),
                downloadButton("download_enrichment_results", "Pathway Enrichment Results", 
                              class = "btn-block btn-lg"),
                br(), br(),
                downloadButton("download_network_results", "Network Inference Results", 
                              class = "btn-block btn-lg")
              ),
              
              column(6,
                downloadButton("download_trajectory_results", "Trajectory Analysis Results", 
                              class = "btn-block btn-lg"),
                br(), br(),
                downloadButton("download_cell_state_results", "Cell State Results", 
                              class = "btn-block btn-lg"),
                br(), br(),
                downloadButton("download_all_results", "All Results (ZIP)", 
                              class = "btn-block btn-lg btn-primary")
              )
            ),
            
            br(),
            tags$h4("Session Information:"),
            verbatimTextOutput("session_info")
          )
        )
      )
    )
  )
)

# Server Logic
server <- function(input, output, session) {
  
  # Load results
  results <- reactive({
    load_results()
  })
  
  # Navigation buttons
  observeEvent(input$go_expression, {
    updateTabItems(session, "sidebarMenu", "expression")
  })
  
  observeEvent(input$go_de, {
    updateTabItems(session, "sidebarMenu", "de")
  })
  
  observeEvent(input$go_enrichment, {
    updateTabItems(session, "sidebarMenu", "enrichment")
  })
  
  observeEvent(input$go_networks, {
    updateTabItems(session, "sidebarMenu", "networks")
  })
  
  observeEvent(input$go_trajectory, {
    updateTabItems(session, "sidebarMenu", "trajectory")
  })
  
  observeEvent(input$go_cellstates, {
    updateTabItems(session, "sidebarMenu", "cellstates")
  })
  
  # Overview boxes
  output$n_genes <- renderValueBox({
    res <- results()
    if (!is.null(res$expression)) {
      valueBox(
        value = nrow(res$expression),
        subtitle = "Genes Analyzed",
        icon = icon("dna"),
        color = "blue"
      )
    }
  })
  
  output$n_samples <- renderValueBox({
    res <- results()
    if (!is.null(res$sample_info)) {
      valueBox(
        value = nrow(res$sample_info),
        subtitle = "Samples Analyzed",
        icon = icon("flask"),
        color = "green"
      )
    }
  })
  
  output$n_de_genes <- renderValueBox({
    res <- results()
    if (!is.null(res$de_results)) {
      total_de <- sum(sapply(res$de_results, function(x) x$n_total_significant))
      valueBox(
        value = total_de,
        subtitle = "Differentially Expressed Genes",
        icon = icon("not-equal"),
        color = "yellow"
      )
    }
  })
  
  # Gene expression plot
  output$gene_expression_plot <- renderPlotly({
    res <- results()
    if (!is.null(res$expression) && !is.null(input$gene_input) && input$gene_input != "") {
      gene <- input$gene_input
      
      if (gene %in% rownames(res$expression)) {
        expr_data <- data.frame(
          sample = colnames(res$expression),
          expression = as.numeric(res$expression[gene, ]),
          stage = res$sample_info$stage
        )
        
        if (input$stage_filter != "All") {
          expr_data <- expr_data[expr_data$stage == input$stage_filter, ]
        }
        
        p <- ggplot(expr_data, aes(x = stage, y = expression, color = stage)) +
          geom_boxplot(alpha = input$alpha_value) +
          geom_jitter(width = 0.2, size = input$point_size, alpha = input$alpha_value) +
          labs(title = paste("Expression of", gene),
               x = "Developmental Stage",
               y = "Normalized Expression") +
          theme_minimal() +
          scale_color_viridis_d()
        
        if (input$show_labels) {
          p <- p + geom_text(aes(label = sample), vjust = -1, size = 3)
        }
        
        ggplotly(p)
      }
    }
  })
  
  # Expression heatmap
  output$expression_heatmap <- renderPlot({
    res <- results()
    if (!is.null(res$expression)) {
      # Select top variable genes
      gene_var <- apply(res$expression, 1, var)
      top_genes <- head(order(gene_var, decreasing = TRUE), 50)
      
      heatmap_data <- res$expression[top_genes, ]
      
      # Z-score normalization
      heatmap_data <- t(scale(t(heatmap_data)))
      
      # Create annotation
      annotation <- data.frame(
        stage = res$sample_info$stage,
        row.names = colnames(heatmap_data)
      )
      
      Heatmap(heatmap_data,
               name = "Z-score",
               col = colorRampPalette(c("blue", "white", "red"))(100),
               show_row_names = TRUE,
               show_column_names = FALSE,
               cluster_rows = TRUE,
               cluster_columns = FALSE,
               top_annotation = HeatmapAnnotation(stage = annotation$stage,
                                                 col = list(stage = c("P1" = "#E41A1C", "P12" = "#377EB8", "P28" = "#4DAF4A"))))
    }
  })
  
  # Volcano plot
  output$volcano_plot <- renderPlotly({
    res <- results()
    if (!is.null(res$de_results) && input$update_de > 0) {
      contrast <- gsub(" vs ", "_vs_", input$de_contrast)
      
      if (contrast %in% names(res$de_results)) {
        de_data <- res$de_results[[contrast]]$full_results
        
        # Add significance categories
        de_data$significance <- ifelse(
          de_data$padj < input$p_threshold & abs(de_data$log2FoldChange) > input$fc_threshold,
          ifelse(de_data$log2FoldChange > 0, "Upregulated", "Downregulated"),
          "Not significant"
        )
        
        p <- ggplot(de_data, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
          geom_point(alpha = input$alpha_value, size = input$point_size) +
          geom_vline(xintercept = c(-input$fc_threshold, input$fc_threshold), linetype = "dashed", alpha = 0.5) +
          geom_hline(yintercept = -log10(input$p_threshold), linetype = "dashed", alpha = 0.5) +
          labs(title = paste("Volcano Plot:", input$de_contrast),
               x = "Log2 Fold Change",
               y = "-Log10 Adjusted P-value") +
          theme_minimal() +
          scale_color_manual(values = c("Not significant" = "gray", "Upregulated" = "red", "Downregulated" = "blue"))
        
        if (input$show_labels) {
          # Label significant genes
          sig_genes <- de_data[de_data$significance != "Not significant", ]
          top_genes <- head(sig_genes[order(abs(sig_genes$log2FoldChange), decreasing = TRUE), ], 10)
          
          p <- p + geom_text_repel(data = top_genes, aes(label = rownames(top_genes)), 
                                  size = 3, max.overlaps = 5)
        }
        
        ggplotly(p)
      }
    }
  })
  
  # Top upregulated genes table
  output$up_genes_table <- renderDT({
    res <- results()
    if (!is.null(res$de_results) && input$update_de > 0) {
      contrast <- gsub(" vs ", "_vs_", input$de_contrast)
      
      if (contrast %in% names(res$de_results)) {
        up_genes <- res$de_results[[contrast]]$full_results
        up_genes <- up_genes[up_genes$log2FoldChange > input$fc_threshold & 
                           up_genes$padj < input$p_threshold, ]
        up_genes <- up_genes[order(up_genes$log2FoldChange, decreasing = TRUE), ]
        up_genes <- head(up_genes, 20)
        
        datatable(
          data.frame(
            Gene = rownames(up_genes),
            Log2FC = round(up_genes$log2FoldChange, 2),
            P_value = signif(up_genes$padj, 3),
            stringsAsFactors = FALSE
          ),
          options = list(pageLength = 10, scrollX = TRUE),
          rownames = FALSE
        )
      }
    }
  })
  
  # Top downregulated genes table
  output$down_genes_table <- renderDT({
    res <- results()
    if (!is.null(res$de_results) && input$update_de > 0) {
      contrast <- gsub(" vs ", "_vs_", input$de_contrast)
      
      if (contrast %in% names(res$de_results)) {
        down_genes <- res$de_results[[contrast]]$full_results
        down_genes <- down_genes[down_genes$log2FoldChange < -input$fc_threshold & 
                               down_genes$padj < input$p_threshold, ]
        down_genes <- down_genes[order(down_genes$log2FoldChange), ]
        down_genes <- head(down_genes, 20)
        
        datatable(
          data.frame(
            Gene = rownames(down_genes),
            Log2FC = round(down_genes$log2FoldChange, 2),
            P_value = signif(down_genes$padj, 3),
            stringsAsFactors = FALSE
          ),
          options = list(pageLength = 10, scrollX = TRUE),
          rownames = FALSE
        )
      }
    }
  })
  
  # Enrichment plot
  output$enrichment_plot <- renderPlotly({
    res <- results()
    if (!is.null(res$enrichment) && input$run_enrichment > 0) {
      enrichment_data <- NULL
      
      if (input$enrichment_type == "Hallmark Gene Sets") {
        contrast <- gsub(" vs ", "_vs_", input$enrichment_contrast)
        if (!is.null(res$enrichment$hallmark) && contrast %in% names(res$enrichment$hallmark)) {
          hallmark_res <- res$enrichment$hallmark[[contrast]]
          if (!is.null(hallmark_res) && !is.character(hallmark_res)) {
            enrichment_data <- hallmark_res@result
            enrichment_data <- head(enrichment_data[order(enrichment_data$p.adjust), ], input$top_pathways)
          }
        }
      }
      
      if (!is.null(enrichment_data) && nrow(enrichment_data) > 0) {
        p <- ggplot(enrichment_data, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                                        size = Count, color = p.adjust)) +
          geom_point(alpha = input$alpha_value) +
          labs(title = paste(input$enrichment_type, "-", input$enrichment_contrast),
               x = "Gene Ratio",
               y = "Pathway",
               size = "Gene Count",
               color = "Adjusted P-value") +
          theme_minimal() +
          scale_color_viridis_c()
        
        ggplotly(p)
      }
    }
  })
  
  # Pathway details table
  output$pathway_details_table <- renderDT({
    res <- results()
    if (!is.null(res$enrichment) && input$run_enrichment > 0) {
      enrichment_data <- NULL
      
      if (input$enrichment_type == "Hallmark Gene Sets") {
        contrast <- gsub(" vs ", "_vs_", input$enrichment_contrast)
        if (!is.null(res$enrichment$hallmark) && contrast %in% names(res$enrichment$hallmark)) {
          hallmark_res <- res$enrichment$hallmark[[contrast]]
          if (!is.null(hallmark_res) && !is.character(hallmark_res)) {
            enrichment_data <- hallmark_res@result
            enrichment_data <- head(enrichment_data[order(enrichment_data$p.adjust), ], input$top_pathways)
          }
        }
      }
      
      if (!is.null(enrichment_data) && nrow(enrichment_data) > 0) {
        display_data <- data.frame(
          Pathway = enrichment_data$Description,
          GeneRatio = enrichment_data$GeneRatio,
          Pvalue = signif(enrichment_data$p.adjust, 3),
          Count = enrichment_data$Count,
          stringsAsFactors = FALSE
        )
        
        datatable(display_data, options = list(pageLength = 10, scrollX = TRUE))
      }
    }
  })
  
  # Network plot
  output$network_plot <- renderVisNetwork({
    res <- results()
    if (!is.null(res$network) && input$generate_network > 0) {
      network_data <- NULL
      
      if (input$network_type == "GENIE3" && !is.null(res$network$genie3)) {
        network_data <- res$network$genie3
      } else if (input$network_type == "Correlation" && !is.null(res$network$correlation)) {
        network_data <- res$network$correlation
      } else if (input$network_type == "SCENIC" && !is.null(res$network$scenic)) {
        network_data <- res$network$scenic
      }
      
      if (!is.null(network_data) && !is.null(network_data$graph)) {
        g <- network_data$graph
        
        # Filter by threshold
        edge_list <- as_data_frame(g, what = "edges")
        if (nrow(edge_list) > 0 && "weight" %in% colnames(edge_list)) {
          edge_list <- edge_list[edge_list$weight >= input$edge_threshold, ]
        }
        
        # Limit nodes
        if (nrow(edge_list) > 0) {
          nodes_to_keep <- unique(c(edge_list$from, edge_list$to))
          if (length(nodes_to_keep) > input$max_nodes) {
            nodes_to_keep <- sample(nodes_to_keep, input$max_nodes)
            edge_list <- edge_list[edge_list$from %in% nodes_to_keep & 
                                 edge_list$to %in% nodes_to_keep, ]
          }
          
          # Create visNetwork
          nodes <- data.frame(
            id = nodes_to_keep,
            label = nodes_to_keep,
            group = ifelse(nodes_to_keep %in% load_transcription_factors(), "TF", "Target"),
            size = 10
          )
          
          edges <- data.frame(
            from = edge_list$from,
            to = edge_list$to,
            width = ifelse("weight" %in% colnames(edge_list), edge_list$weight * 5, 1),
            arrows = "to"
          )
          
          visNetwork(nodes, edges) %>%
            visNodes(shape = "dot") %>%
            visEdges(smooth = list(type = "continuous")) %>%
            visLayout(randomSeed = 123) %>%
            visOptions(highlightNearest = list(enabled = TRUE, degree = 1))
        }
      }
    }
  })
  
  # Hub genes table
  output$hub_genes_table <- renderDT({
    res <- results()
    if (!is.null(res$network) && input$generate_network > 0) {
      # Simple hub gene identification based on degree
      network_data <- NULL
      
      if (input$network_type == "GENIE3" && !is.null(res$network$genie3)) {
        network_data <- res$network$genie3
      } else if (input$network_type == "Correlation" && !is.null(res$network$correlation)) {
        network_data <- res$network$correlation
      } else if (input$network_type == "SCENIC" && !is.null(res$network$scenic)) {
        network_data <- res$network$scenic
      }
      
      if (!is.null(network_data) && !is.null(network_data$graph)) {
        g <- network_data$graph
        degree_centrality <- degree(g)
        
        hub_genes <- data.frame(
          Gene = names(degree_centrality),
          Degree = as.numeric(degree_centrality),
          stringsAsFactors = FALSE
        )
        
        hub_genes <- hub_genes[order(hub_genes$Degree, decreasing = TRUE), ]
        hub_genes <- head(hub_genes, 20)
        
        datatable(hub_genes, options = list(pageLength = 10), rownames = FALSE)
      }
    }
  })
  
  # Network statistics
  output$network_stats <- renderTable({
    res <- results()
    if (!is.null(res$network) && input$generate_network > 0) {
      network_data <- NULL
      
      if (input$network_type == "GENIE3" && !is.null(res$network$genie3)) {
        network_data <- res$network$genie3
      } else if (input$network_type == "Correlation" && !is.null(res$network$correlation)) {
        network_data <- res$network$correlation
      } else if (input$network_type == "SCENIC" && !is.null(res$network$scenic)) {
        network_data <- res$network$scenic
      }
      
      if (!is.null(network_data) && !is.null(network_data$graph)) {
        g <- network_data$graph
        
        data.frame(
          Metric = c("Nodes", "Edges", "Density", "Average Degree"),
          Value = c(
            vcount(g),
            ecount(g),
            round(edge_density(g), 3),
            round(mean(degree(g)), 2)
          )
        )
      }
    }
  })
  
  # Trajectory plot
  output$trajectory_plot <- renderPlotly({
    res <- results()
    if (!is.null(res$interactive) && input$generate_trajectory > 0) {
      traj_data <- NULL
      
      if (input$trajectory_method == "UMAP" && !is.null(res$interactive$umap)) {
        traj_data <- res$interactive$umap
      } else if (input$trajectory_method == "t-SNE" && !is.null(res$interactive$tsne)) {
        traj_data <- res$interactive$tsne
      }
      
      if (!is.null(traj_data)) {
        color_var <- "stage"
        if (input$color_by == "Pseudotime" && !is.null(res$interactive$pseudotime)) {
          traj_data$pseudotime <- res$interactive$pseudotime$pseudotime
          color_var <- "pseudotime"
        } else if (input$color_by == "Cell State" && !is.null(res$interactive$cell_states)) {
          traj_data$cell_state <- res$interactive$cell_states$cluster
          color_var <- "cell_state"
        }
        
        coords <- colnames(traj_data)[1:2]
        
        p <- ggplot(traj_data, aes_string(x = coords[1], y = coords[2], color = color_var)) +
          geom_point(size = input$point_size, alpha = input$alpha_value) +
          labs(title = paste(input$trajectory_method, "Trajectory"),
               x = coords[1],
               y = coords[2],
               color = input$color_by) +
          theme_minimal()
        
        if (input$show_labels) {
          p <- p + geom_text(aes(label = replicate), vjust = -1, size = 3)
        }
        
        ggplotly(p)
      }
    }
  })
  
  # Pseudotime plot
  output$pseudotime_plot <- renderPlotly({
    res <- results()
    if (!is.null(res$interactive) && !is.null(res$interactive$pseudotime) && input$generate_trajectory > 0) {
      pseudo_data <- res$interactive$pseudotime
      
      p <- ggplot(pseudo_data, aes(x = stage, y = pseudotime, color = stage)) +
        geom_boxplot(alpha = input$alpha_value) +
        geom_jitter(width = 0.2, size = input$point_size, alpha = input$alpha_value) +
        labs(title = "Pseudotime Distribution",
             x = "Developmental Stage",
             y = "Pseudotime") +
        theme_minimal() +
        scale_color_viridis_d()
      
      ggplotly(p)
    }
  })
  
  # Cell state plot
  output$cell_state_plot <- renderPlotly({
    res <- results()
    if (!is.null(res$interactive) && !is.null(res$interactive$cell_states) && input$run_clustering > 0) {
      cell_data <- res$interactive$cell_states
      
      if (input$clustering_method == "Consensus") {
        method_data <- cell_data$consensus
      } else if (input$clustering_method == "K-means") {
        method_data <- cell_data$kmeans
      } else {
        method_data <- cell_data$consensus  # Default
      }
      
      if (!is.null(res$interactive$umap)) {
        plot_data <- merge(res$interactive$umap, method_data, by.x = "sample_id", by.y = "cell_id")
        
        p <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = factor(cluster))) +
          geom_point(size = input$point_size, alpha = input$alpha_value) +
          labs(title = paste(input$clustering_method, "Cell States"),
               x = "UMAP1",
               y = "UMAP2",
               color = "Cell State") +
          theme_minimal() +
          scale_color_viridis_d()
        
        ggplotly(p)
      }
    }
  })
  
  # Cluster composition
  output$cluster_composition <- renderTable({
    res <- results()
    if (!is.null(res$interactive) && !is.null(res$interactive$cell_states) && input$run_clustering > 0) {
      cell_data <- res$interactive$cell_states
      
      if (input$clustering_method == "Consensus") {
        method_data <- cell_data$consensus
      } else if (input$clustering_method == "K-means") {
        method_data <- cell_data$kmeans
      } else {
        method_data <- cell_data$consensus  # Default
      }
      
      composition <- table(method_data$cluster, method_data$stage)
      as.data.frame.matrix(composition)
    }
  })
  
  # Download handlers
  output$download_expression_plot <- downloadHandler(
    filename = function() {
      paste("gene_expression_", input$gene_input, ".pdf", sep = "")
    },
    content = function(file) {
      # Implementation for downloading expression plot
      ggsave(file, width = 10, height = 6, dpi = 300)
    }
  )
  
  output$download_volcano <- downloadHandler(
    filename = function() {
      paste("volcano_plot_", gsub(" vs ", "_vs_", input$de_contrast), ".pdf", sep = "")
    },
    content = function(file) {
      # Implementation for downloading volcano plot
      ggsave(file, width = 10, height = 8, dpi = 300)
    }
  )
  
  output$download_enrichment <- downloadHandler(
    filename = function() {
      paste("enrichment_", gsub(" ", "_", tolower(input$enrichment_type)), ".pdf", sep = "")
    },
    content = function(file) {
      # Implementation for downloading enrichment plot
      ggsave(file, width = 10, height = 8, dpi = 300)
    }
  )
  
  output$download_de_results <- downloadHandler(
    filename = function() {
      "differential_expression_results.zip"
    },
    content = function(file) {
      # Create zip file with DE results
      tmpdir <- tempdir()
      # Implementation for creating zip file
      zip(file, files = list.files(tmpdir, full.names = TRUE))
    }
  )
  
  output$download_all_results <- downloadHandler(
    filename = function() {
      paste("satellite_cell_analysis_results_", Sys.Date(), ".zip", sep = "")
    },
    content = function(file) {
      # Create comprehensive zip file
      # Implementation for creating complete results zip
    }
  )
  
  # Session info
  output$session_info <- renderPrint({
    sessionInfo()
  })
}

# Load transcription factors helper function
load_transcription_factors <- function() {
  c("Pax7", "MyoD", "Myf5", "Myogenin", "Mrf4", "Six1", "Six4", 
    "Eya1", "Eya2", "Tcf4", "Tcf12", "Hand1", "Hand2", "Twist1", "Twist2",
    "Sp1", "Sp2", "Sp3", "Zeb1", "Zeb2", "Snai1", "Snai2",
    "Foxo1", "Foxo3", "Foxo4", "Foxc1", "Foxc2", "Foxp1",
    "Nr1d1", "Nr1d2", "Ppara", "Ppard", "Pparg", "Esrra", "Esrrg",
    "Creb1", "Atf1", "Atf2", "Atf3", "Atf4", "Jun", "Fos", "Myc",
    "E2f1", "E2f2", "E2f3", "E2f4", "E2f5", "Rb1", "Trp53",
    "Sox2", "Sox6", "Sox8", "Sox9", "Sox10", "Nanog", "Oct4",
    "Bcl6", "Hes1", "Hey1", "Hey2", "Id1", "Id2", "Id3")
}

# Run the application
shinyApp(ui = ui, server = server)
