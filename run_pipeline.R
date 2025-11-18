#!/usr/bin/env Rscript
# Satellite Cell Development Transcriptomics Pipeline
# Main Pipeline Runner (MIT-level, but not MIT-level fragile 😅)

# -------------------------
# 0. Basic safety checks
# -------------------------
if (!file.exists("src/utils.R")) {
  stop("utils.R not found. Make sure you are running this script from the project root (where src/ lives).")
}

# Load utility functions (logging, config, helpers)
source("src/utils.R")

# -------------------------
# 1. Setup logging & config
# -------------------------
logfile <- setup_logging()
loginfo(paste("Logging initialized. Log file:", logfile))

config <- load_config()
loginfo("Configuration loaded from config.yaml")

# Create progress tracker
progress_file <- create_progress_tracker()
loginfo(paste("Progress tracker initialized:", progress_file))

# System info
sys_info <- get_system_info()
loginfo(paste("System:", sys_info$platform, "| R version:", sys_info$version))

# -------------------------
# 2. Dependency management
# -------------------------
required_packages <- c(
  "GEOquery", "Biobase", "limma", "affy", "preprocessCore",
  "DESeq2", "edgeR",
  "clusterProfiler", "org.Mm.eg.db", "DOSE", "enrichplot", "GSVA", "GSEABase",
  "igraph", "ggnetwork", "ggraph", "tidygraph", "GENIE3",
  "umap", "Rtsne", "phate", "destiny", "monocle3",
  "ggplot2", "ComplexHeatmap", "viridis",
  "yaml", "logging", "dplyr", "tidyr", "jsonlite", "R.utils"
)

# check_dependencies() should be defined in utils.R
missing_packages <- check_dependencies(required_packages)

if (!is.null(missing_packages) && length(missing_packages) > 0) {
  logwarn(paste("Missing packages detected:", paste(missing_packages, collapse = ", ")))
  loginfo("Attempting to install missing packages via install_missing_packages()")

  install_missing_packages(missing_packages)

  # Re-check after installation
  missing_after <- check_dependencies(required_packages)
  if (!is.null(missing_after) && length(missing_after) > 0) {
    logerror(paste(
      "Some packages are still missing after installation attempt:",
      paste(missing_after, collapse = ", ")
    ))
    stop("Unresolved package dependencies. Please install these manually and re-run the pipeline.")
  } else {
    loginfo("All required packages are now installed.")
  }
} else {
  loginfo("All required packages are already installed.")
}

# Optionally load them (if your code expects them to be attached)
suppressPackageStartupMessages(
  invisible(lapply(required_packages, function(p) {
    # Some Bioc packages might not be on CRAN, but library() works once installed
    tryCatch(library(p, character.only = TRUE), error = function(e) {
      logwarn(paste("Could not attach package:", p, "->", e$message))
    })
  }))
)

# -------------------------
# 3. Directory structure
# -------------------------
create_directory_structure()
loginfo("Directory structure verified/created.")

# Initialize results tracking
results_list <- list()

# Helper to safely source a step
run_step <- function(step_name, script_path, key) {
  loginfo(paste("STEP", step_name, ":", key))
  if (!file.exists(script_path)) {
    msg <- paste("Script not found:", script_path)
    logerror(msg)
    update_progress(key, "error", msg)
    stop(msg)
  }

  tryCatch({
    source(script_path)
    update_progress(key, "completed")
    results_list[[key]] <<- "completed"
    loginfo(paste("Step", key, "completed successfully."))
  }, error = function(e) {
    err_msg <- conditionMessage(e)
    logerror(paste("Step", key, "failed:", err_msg))
    update_progress(key, "error", err_msg)
    stop(paste("Pipeline halted at step", key, "->", err_msg))
  })
}

# -------------------------
# 4. Run pipeline steps
# -------------------------

# Step 0: Data Download
run_step("0: Data Download", "src/00_download.R", "data_download")

# Step 1: Preprocessing and Quality Control
run_step("1: Preprocessing and Quality Control", "src/01_preprocess.R", "preprocessing")

# Step 2: Differential Expression Analysis
run_step("2: Differential Expression Analysis", "src/02_de_analysis.R", "de_analysis")

# Step 3: Pathway Enrichment Analysis
run_step("3: Pathway Enrichment Analysis", "src/03_enrichment.R", "enrichment_analysis")

# Step 4: Network Inference
run_step("4: Gene Regulatory Network Inference", "src/04_network_inference.R", "network_inference")

# Step 5: Advanced Visualizations
run_step("5: Advanced Visualizations", "src/05_visualizations.R", "visualizations")

# -------------------------
# 5. Final report & checks
# -------------------------
loginfo("GENERATING FINAL SUMMARY REPORT")
final_report <- NULL
tryCatch({
  final_report <- generate_summary_report(results_list)
  loginfo(paste("Summary report generated:", final_report))
}, error = function(e) {
  logwarn(paste("Could not generate summary report:", conditionMessage(e)))
})

final_memory <- check_memory_usage()
loginfo(paste("Final memory usage:", final_memory$used_mb, "MB"))

loginfo("PIPELINE COMPLETED SUCCESSFULLY")
loginfo(paste("Results saved to:", getwd()))

if (!is.null(final_report)) {
  loginfo(paste("Final report:", final_report))
}

# -------------------------
# 6. Console completion message
# -------------------------
cat("\n=============================================\n")
cat("SATELLITE CELL DEVELOPMENT PIPELINE COMPLETED\n")
cat("=============================================\n\n")
cat("Results directory: ", getwd(), "\n", sep = "")
if (!is.null(final_report)) {
  cat("Final report:      ", final_report, "\n", sep = "")
}
cat("Completion time:   ", as.character(Sys.time()), "\n", sep = "")
cat("Note: This pipeline is designed for MIT-level research quality analysis.\n\n")
