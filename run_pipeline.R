#!/usr/bin/env Rscript
# Satellite Cell Development Transcriptomics Pipeline
# Main Pipeline Runner

# Load utility functions
source("src/utils.R")

# Setup logging and configuration
logfile <- setup_logging()
config <- load_config()

# Create progress tracker
progress_file <- create_progress_tracker()

# Check system information
sys_info <- get_system_info()
loginfo(paste("System:", sys_info$platform, "| R version:", sys_info$version))

# Check dependencies
required_packages <- c(
  "GEOquery", "Biobase", "limma", "affy", "preprocessCore", "DESeq2", "edgeR",
  "clusterProfiler", "org.Mm.eg.db", "DOSE", "enrichplot", "GSVA", "GSEABase",
  "igraph", "ggnetwork", "ggraph", "tidygraph", "GENIE3", "umap", "Rtsne", 
  "phate", "destiny", "monocle3", "ggplot2", "ComplexHeatmap", "viridis",
  "yaml", "logging", "dplyr", "tidyr", "jsonlite", "R.utils"
)

missing_packages <- check_dependencies(required_packages)
if (!is.null(missing_packages)) {
  logwarn(paste("Installing missing packages:", paste(missing_packages, collapse = ", ")))
  install_missing_packages(missing_packages)
}

# Create directory structure
create_directory_structure()

# Initialize results tracking
results_list <- list()

# Step 0: Data Download
loginfo("STEP 0: Data Download")
tryCatch({
  source("src/00_download.R")
  update_progress("data_download", "completed")
  results_list$data_download <- "completed"
}, error = function(e) {
  update_progress("data_download", "error", conditionMessage(e))
  stop("Data download failed")
})

# Step 1: Preprocessing and Quality Control
loginfo("STEP 1: Preprocessing and Quality Control")
tryCatch({
  source("src/01_preprocess.R")
  update_progress("preprocessing", "completed")
  results_list$preprocessing <- "completed"
}, error = function(e) {
  update_progress("preprocessing", "error", conditionMessage(e))
  stop("Preprocessing failed")
})

# Step 2: Differential Expression Analysis
loginfo("STEP 2: Differential Expression Analysis")
tryCatch({
  source("src/02_de_analysis.R")
  update_progress("de_analysis", "completed")
  results_list$de_analysis <- "completed"
}, error = function(e) {
  update_progress("de_analysis", "error", conditionMessage(e))
  stop("DE analysis failed")
})

# Step 3: Pathway Enrichment Analysis
loginfo("STEP 3: Pathway Enrichment Analysis")
tryCatch({
  source("src/03_enrichment.R")
  update_progress("enrichment_analysis", "completed")
  results_list$enrichment_analysis <- "completed"
}, error = function(e) {
  update_progress("enrichment_analysis", "error", conditionMessage(e))
  stop("Enrichment analysis failed")
})

# Step 4: Network Inference
loginfo("STEP 4: Gene Regulatory Network Inference")
tryCatch({
  source("src/04_network_inference.R")
  update_progress("network_inference", "completed")
  results_list$network_inference <- "completed"
}, error = function(e) {
  update_progress("network_inference", "error", conditionMessage(e))
  stop("Network inference failed")
})

# Step 5: Advanced Visualizations
loginfo("STEP 5: Advanced Visualizations")
tryCatch({
  source("src/05_visualizations.R")
  update_progress("visualizations", "completed")
  results_list$visualizations <- "completed"
}, error = function(e) {
  update_progress("visualizations", "error", conditionMessage(e))
  stop("Visualization failed")
})

# Generate final report
loginfo("GENERATING FINAL REPORT")
final_report <- generate_summary_report(results_list)

# Final system check
final_memory <- check_memory_usage()
loginfo(paste("Final memory usage:", final_memory$used_mb, "MB"))

# Pipeline completion
loginfo("PIPELINE COMPLETED SUCCESSFULLY")
loginfo(paste("Results saved to:", getwd()))
loginfo(paste("Final report:", final_report))
loginfo("All analysis steps completed successfully")

# Print completion message
cat("\nSATELLITE CELL DEVELOPMENT TRANSCRIPTOMICS PIPELINE COMPLETED\n")
cat("Results available in:", getwd(), "\n")
cat("Final report:", final_report, "\n")
cat("Analysis completed at:", Sys.time(), "\n")
cat("This pipeline represents MIT-level research quality analysis\n")
