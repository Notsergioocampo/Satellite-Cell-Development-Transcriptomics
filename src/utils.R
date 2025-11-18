#!/usr/bin/env Rscript
# Satellite Cell Development Transcriptomics Pipeline
# Utility Functions
# MIT-Level Research Pipeline

suppressPackageStartupMessages({
  library(logging)
  library(yaml)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(jsonlite)
  library(R.utils)
  library(tools)
  library(parallel)
  # For palettes
  library(viridis)
  library(RColorBrewer)
})

# --------------------------------------------------
# Logging configuration
# --------------------------------------------------
setup_logging <- function(log_dir = "results/logs", log_level = "INFO") {
  if (!dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE)
  }

  logfile <- file.path(log_dir, paste0(format(Sys.time(), "%Y%m%d"), "_pipeline.log"))

  # Basic console logger
  basicConfig()
  # Attach file handler
  addHandler(writeToFile, file = logfile, level = log_level)

  loginfo("Logging initialized")
  return(logfile)
}

# --------------------------------------------------
# Configuration management
# --------------------------------------------------
load_config <- function(config_file = "config.yaml") {
  if (!file.exists(config_file)) {
    stop(paste("Configuration file not found:", config_file))
  }

  config <- yaml::read_yaml(config_file)
  loginfo(paste("Configuration loaded from:", config_file))

  return(config)
}

# --------------------------------------------------
# Data validation
# --------------------------------------------------
validate_expression_matrix <- function(expression_matrix, sample_info) {
  loginfo("Validating expression matrix...")

  validation_results <- list(
    timestamp           = Sys.time(),
    n_genes             = nrow(expression_matrix),
    n_samples           = ncol(expression_matrix),
    missing_values      = sum(is.na(expression_matrix)),
    negative_values     = sum(expression_matrix < 0, na.rm = TRUE),
    infinite_values     = sum(!is.finite(expression_matrix)),
    sample_matching     = all(colnames(expression_matrix) == sample_info$sample_id),
    gene_names_unique   = length(unique(rownames(expression_matrix))) == nrow(expression_matrix),
    sample_names_unique = length(unique(colnames(expression_matrix))) == ncol(expression_matrix)
  )

  issues <- c()
  if (validation_results$missing_values > 0)
    issues <- c(issues, paste("Found", validation_results$missing_values, "missing values"))
  if (validation_results$negative_values > 0)
    issues <- c(issues, paste("Found", validation_results$negative_values, "negative values"))
  if (validation_results$infinite_values > 0)
    issues <- c(issues, paste("Found", validation_results$infinite_values, "infinite values"))
  if (!validation_results$sample_matching)
    issues <- c(issues, "Sample names don't match between expression matrix and sample info")
  if (!validation_results$gene_names_unique)
    issues <- c(issues, "Gene names are not unique")
  if (!validation_results$sample_names_unique)
    issues <- c(issues, "Sample names are not unique")

  validation_results$issues <- issues
  validation_results$status <- if (length(issues) == 0) "PASS" else "WARN"

  loginfo(paste("Validation completed. Status:", validation_results$status))
  return(validation_results)
}

# --------------------------------------------------
# Directory structure
# --------------------------------------------------
create_directory_structure <- function(base_dir = ".") {
  loginfo("Creating directory structure...")

  directories <- c(
    "data/raw",
    "data/processed",
    "data/metadata",
    "src",
    "notebooks",
    "figures",
    "results",
    "results/logs",
    "shiny_app",
    "tests",
    "docs"
  )

  created_dirs <- c()
  for (dir in directories) {
    full_path <- file.path(base_dir, dir)
    if (!dir.exists(full_path)) {
      dir.create(full_path, recursive = TRUE)
      created_dirs <- c(created_dirs, full_path)
    }
  }

  if (length(created_dirs) > 0) {
    loginfo(paste("Created directories:", paste(created_dirs, collapse = ", ")))
  } else {
    loginfo("All directories already exist")
  }

  return(created_dirs)
}

# --------------------------------------------------
# Progress tracking
# --------------------------------------------------
create_progress_tracker <- function() {
  if (!dir.exists("results")) dir.create("results", recursive = TRUE)
  progress_file <- file.path("results", "pipeline_progress.json")

  initial_progress <- list(
    timestamp      = Sys.time(),
    status         = "initialized",
    completed_steps = list(),
    current_step   = NULL,
    errors         = list(),
    warnings       = list()
  )

  write_json(initial_progress, progress_file, auto_unbox = TRUE, pretty = TRUE)
  return(progress_file)
}

update_progress <- function(step_name, status = "completed", error = NULL, warning = NULL) {
  progress_file <- file.path("results", "pipeline_progress.json")

  if (file.exists(progress_file)) {
    progress <- read_json(progress_file)
  } else {
    create_progress_tracker()
    progress <- read_json(progress_file)
  }

  progress$timestamp    <- as.character(Sys.time())
  progress$current_step <- step_name

  if (status == "completed") {
    progress$completed_steps <- unique(c(progress$completed_steps, step_name))
    loginfo(paste("Step completed:", step_name))
  } else if (status == "error") {
    progress$errors <- c(progress$errors, list(
      step      = step_name,
      message   = error,
      timestamp = as.character(Sys.time())
    ))
    logerror(paste("Error in step:", step_name, "-", error))
  } else if (status == "warning") {
    progress$warnings <- c(progress$warnings, list(
      step      = step_name,
      message   = warning,
      timestamp = as.character(Sys.time())
    ))
    logwarn(paste("Warning in step:", step_name, "-", warning))
  } else if (status == "running") {
    loginfo(paste("Step started:", step_name))
  }

  write_json(progress, progress_file, auto_unbox = TRUE, pretty = TRUE)
  return(progress)
}

# --------------------------------------------------
# Memory management
# --------------------------------------------------
check_memory_usage <- function() {
  memory_info <- gc()
  memory_stats <- list(
    timestamp      = Sys.time(),
    used_mb        = round(sum(memory_info[, 2]) / 1024, 2),
    gc_trigger_mb  = round(sum(memory_info[, 3]) / 1024, 2),
    max_used_mb    = round(sum(memory_info[, 4]) / 1024, 2)
  )
  loginfo(paste("Memory usage:", memory_stats$used_mb, "MB"))
  return(memory_stats)
}

optimize_memory <- function() {
  loginfo("Optimizing memory usage...")
  gc()
  loginfo("Memory optimization completed")
}

# --------------------------------------------------
# Parallel processing
# --------------------------------------------------
setup_parallel_processing <- function(n_cores = NULL) {
  if (is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }

  loginfo(paste("Setting up parallel processing with", n_cores, "cores"))

  if (requireNamespace("doParallel", quietly = TRUE)) {
    library(doParallel)
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    loginfo("Parallel processing initialized with doParallel")
    return(cl)
  } else {
    logwarn("doParallel not available, using sequential processing")
    return(NULL)
  }
}

cleanup_parallel_processing <- function(cl) {
  if (!is.null(cl)) {
    stopCluster(cl)
    loginfo("Parallel processing cleaned up")
  }
}

# --------------------------------------------------
# Data export
# --------------------------------------------------
export_results <- function(results, output_dir = "results",
                           format = c("RData", "csv", "json")) {
  loginfo(paste("Exporting results to", output_dir))
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  format <- match.arg(format, several.ok = TRUE)
  exported_files <- list(rdata = character(), csv = character(), json = character())

  for (result_name in names(results)) {
    result_data <- results[[result_name]]

    if ("RData" %in% format) {
      rdata_file <- file.path(output_dir, paste0(result_name, ".RData"))
      save(result_data, file = rdata_file)
      exported_files$rdata <- c(exported_files$rdata, rdata_file)
    }

    if ("csv" %in% format && is.data.frame(result_data)) {
      csv_file <- file.path(output_dir, paste0(result_name, ".csv"))
      write.csv(result_data, file = csv_file, row.names = TRUE)
      exported_files$csv <- c(exported_files$csv, csv_file)
    }

    if ("json" %in% format) {
      json_file <- file.path(output_dir, paste0(result_name, ".json"))
      write_json(result_data, json_file, auto_unbox = TRUE, pretty = TRUE)
      exported_files$json <- c(exported_files$json, json_file)
    }
  }

  loginfo(paste("Exported", length(unlist(exported_files)), "files"))
  return(exported_files)
}

# --------------------------------------------------
# Stats helpers
# --------------------------------------------------
calculate_effect_size <- function(group1, group2, method = "cohen_d") {
  if (method == "cohen_d") {
    pooled_sd <- sqrt(
      ((length(group1) - 1) * var(group1) +
         (length(group2) - 1) * var(group2)) /
        (length(group1) + length(group2) - 2)
    )
    effect_size <- (mean(group1) - mean(group2)) / pooled_sd
  } else if (method == "log2_fold_change") {
    effect_size <- log2(mean(group1) + 1) - log2(mean(group2) + 1)
  } else {
    stop(paste("Unknown effect size method:", method))
  }

  return(effect_size)
}

perform_multiple_testing_correction <- function(p_values, method = "BH") {
  if (!method %in% c("BH", "bonferroni")) {
    stop(paste("Unknown multiple testing correction method:", method))
  }
  adjusted_p <- p.adjust(p_values, method = method)
  return(adjusted_p)
}

# --------------------------------------------------
# Gene annotation
# --------------------------------------------------
annotate_genes <- function(gene_list, species = "Mus musculus") {
  loginfo(paste("Annotating", length(gene_list), "genes"))

  tryCatch({
    if (species == "Mus musculus") {
      suppressPackageStartupMessages(library(org.Mm.eg.db))

      annotations <- list(
        entrez_id     = mapIds(org.Mm.eg.db, keys = gene_list,
                               keytype = "SYMBOL", column = "ENTREZID"),
        gene_name     = mapIds(org.Mm.eg.db, keys = gene_list,
                               keytype = "SYMBOL", column = "GENENAME"),
        go_terms      = mapIds(org.Mm.eg.db, keys = gene_list,
                               keytype = "SYMBOL", column = "GO"),
        kegg_pathways = mapIds(org.Mm.eg.db, keys = gene_list,
                               keytype = "SYMBOL", column = "PATH")
      )

      annotation_df <- data.frame(
        gene_symbol   = gene_list,
        entrez_id     = annotations$entrez_id,
        gene_name     = annotations$gene_name,
        go_terms      = annotations$go_terms,
        kegg_pathways = annotations$kegg_pathways,
        stringsAsFactors = FALSE
      )

      return(annotation_df)
    } else {
      stop(paste("Species not supported:", species))
    }
  }, error = function(e) {
    logwarn(paste("Gene annotation failed:", e$message))
    return(data.frame(gene_symbol = gene_list, stringsAsFactors = FALSE))
  })
}

# --------------------------------------------------
# QC utilities
# --------------------------------------------------
calculate_qc_metrics <- function(expression_matrix) {
  metrics <- list(
    n_genes                  = nrow(expression_matrix),
    n_samples                = ncol(expression_matrix),
    total_counts             = sum(expression_matrix, na.rm = TRUE),
    mean_counts_per_gene     = mean(rowMeans(expression_matrix, na.rm = TRUE)),
    median_counts_per_gene   = median(rowMeans(expression_matrix, na.rm = TRUE)),
    mean_counts_per_sample   = mean(colMeans(expression_matrix, na.rm = TRUE)),
    median_counts_per_sample = median(colMeans(expression_matrix, na.rm = TRUE)),
    max_expression           = max(expression_matrix, na.rm = TRUE),
    min_expression           = min(expression_matrix, na.rm = TRUE),
    zero_counts              = sum(expression_matrix == 0, na.rm = TRUE),
    zero_proportion          = mean(expression_matrix == 0, na.rm = TRUE)
  )
  return(metrics)
}

detect_outliers <- function(expression_matrix, method = "mahalanobis", threshold = 0.95) {
  if (method == "mahalanobis") {
    center     <- colMeans(expression_matrix)
    cov_matrix <- cov(expression_matrix)

    distances <- mahalanobis(t(expression_matrix), center, cov_matrix)
    threshold_value <- qchisq(threshold, df = nrow(expression_matrix))
    outliers <- which(distances > threshold_value)

    return(list(
      distances = distances,
      outliers  = outliers,
      threshold = threshold_value,
      method    = method
    ))
  } else {
    stop(paste("Unknown outlier detection method:", method))
  }
}

# --------------------------------------------------
# Visualization helpers
# --------------------------------------------------
create_color_palette <- function(n_colors, palette_type = "viridis", alpha = 1.0) {
  if (palette_type == "viridis") {
    colors <- viridis::viridis(n_colors, alpha = alpha)
  } else if (palette_type == "brewer") {
    colors <- RColorBrewer::brewer.pal(min(n_colors, 8), "Set2")
    if (n_colors > length(colors)) {
      colors <- colorRampPalette(colors)(n_colors)
    }
  } else if (palette_type == "custom") {
    colors <- c(
      "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
      "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"
    )
    if (n_colors > length(colors)) {
      colors <- grDevices::colorRampPalette(colors)(n_colors)
    } else {
      colors <- colors[1:n_colors]
    }
  } else {
    colors <- grDevices::rainbow(n_colors, alpha = alpha)
  }

  return(colors)
}

save_plot_with_metadata <- function(plot, filename,
                                    width = 10, height = 8, dpi = 300,
                                    metadata = NULL) {
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  ggsave(filename, plot, width = width, height = height, dpi = dpi)

  if (!is.null(metadata)) {
    metadata_file <- sub("\\.pdf$", "_metadata.json", filename)
    metadata$filename  <- basename(filename)
    metadata$timestamp <- as.character(Sys.time())
    metadata$file_size <- file.size(filename)

    write_json(metadata, metadata_file, auto_unbox = TRUE, pretty = TRUE)
  }

  loginfo(paste("Plot saved:", filename))
}

# --------------------------------------------------
# Pipeline step wrapper
# --------------------------------------------------
run_pipeline_step <- function(step_function, step_name, ...) {
  loginfo(paste("Starting step:", step_name))
  update_progress(step_name, status = "running")

  tryCatch({
    result <- step_function(...)
    update_progress(step_name, status = "completed")
    loginfo(paste("Step completed successfully:", step_name))
    return(result)
  }, error = function(e) {
    error_msg <- conditionMessage(e)
    update_progress(step_name, status = "error", error = error_msg)
    logerror(paste("Step failed:", step_name, "-", error_msg))
    stop(e)
  })
}

# --------------------------------------------------
# Summary report
# --------------------------------------------------
generate_summary_report <- function(results_list,
                                    output_file = "results/pipeline_summary.html") {
  loginfo("Generating pipeline summary report...")

  summary_stats <- list(
    timestamp      = Sys.time(),
    total_steps    = length(results_list),
    completed_steps = names(results_list),
    memory_usage   = check_memory_usage(),
    pipeline_status = "completed"
  )

  html_content <- create_html_report(summary_stats, results_list)
  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
  writeLines(html_content, output_file)

  loginfo(paste("Summary report generated:", output_file))
  return(output_file)
}

create_html_report <- function(summary_stats, results_list) {
  html <- paste0(
    "<!DOCTYPE html>\n<html>\n<head>\n",
    "<title>Satellite Cell Development Transcriptomics Pipeline Report</title>\n",
    "<style>\n",
    "body { font-family: Arial, sans-serif; margin: 40px; }\n",
    "h1 { color: #2E86AB; }\n",
    "h2 { color: #A23B72; }\n",
    ".summary { background-color: #f0f0f0; padding: 20px; border-radius: 5px; }\n",
    ".step { margin-bottom: 20px; padding: 15px; border-left: 4px solid #2E86AB; }\n",
    "</style>\n</head>\n<body>\n",
    "<h1>Satellite Cell Development Transcriptomics Pipeline Report</h1>\n",
    "<div class='summary'>\n",
    "<h2>Pipeline Summary</h2>\n",
    "<p><strong>Timestamp:</strong> ", summary_stats$timestamp, "</p>\n",
    "<p><strong>Total Steps:</strong> ", summary_stats$total_steps, "</p>\n",
    "<p><strong>Memory Usage:</strong> ", summary_stats$memory_usage$used_mb, " MB</p>\n",
    "<p><strong>Status:</strong> ", summary_stats$pipeline_status, "</p>\n",
    "</div>\n",
    "<h2>Analysis Steps</h2>\n"
  )

  for (step_name in names(results_list)) {
    html <- paste0(
      html,
      "<div class='step'>\n",
      "<h3>", step_name, "</h3>\n",
      "<p>Status: Completed</p>\n",
      "</div>\n"
    )
  }

  html <- paste0(html, "</body>\n</html>\n")
  return(html)
}

# --------------------------------------------------
# Error handling
# --------------------------------------------------
handle_error <- function(error, context = "") {
  error_timestamp <- Sys.time()
  error_message   <- paste0("[", error_timestamp, "] ERROR in ", context, ": ",
                            conditionMessage(error))

  logerror(error_message)

  tb <- capture.output(traceback())

  error_report <- list(
    timestamp = as.character(error_timestamp),
    context   = context,
    message   = conditionMessage(error),
    call      = tryCatch(deparse(sys.call(-1)), error = function(e) NA_character_),
    traceback = tb
  )

  dir.create("results/logs", showWarnings = FALSE, recursive = TRUE)
  error_file <- file.path(
    "results", "logs",
    paste0("error_", format(error_timestamp, "%Y%m%d_%H%M%S"), ".json")
  )
  write_json(error_report, error_file, auto_unbox = TRUE, pretty = TRUE)

  return(error_report)
}

safe_execute <- function(expr, context = "") {
  tryCatch(expr, error = function(e) {
    handle_error(e, context)
    return(NULL)
  })
}

# --------------------------------------------------
# Dependency management
# --------------------------------------------------
check_dependencies <- function(required_packages) {
  missing_packages <- vapply(required_packages, function(pkg) {
    !requireNamespace(pkg, quietly = TRUE)
  }, logical(1L))

  missing <- required_packages[missing_packages]
  if (length(missing) > 0) {
    logwarn(paste("Missing packages:", paste(missing, collapse = ", ")))
    return(missing)
  }

  loginfo("All required packages are available")
  return(NULL)
}

install_missing_packages <- function(packages) {
  if (length(packages) > 0) {
    loginfo(paste("Installing missing packages (CRAN only):", paste(packages, collapse = ", ")))
    for (pkg in packages) {
      tryCatch({
        install.packages(pkg, dependencies = TRUE)
        loginfo(paste("Installed package:", pkg))
      }, error = function(e) {
        logerror(paste("Failed to install package:", pkg, "-", conditionMessage(e)))
      })
    }
  }
}

# --------------------------------------------------
# System info
# --------------------------------------------------
get_system_info <- function() {
  mem_limit <- tryCatch({
    if ("memory.limit" %in% ls()) memory.limit() else NA
  }, error = function(e) NA)

  sys_info <- list(
    timestamp         = Sys.time(),
    platform          = R.version$platform,
    version           = R.version.string,
    user              = Sys.info()[["user"]],
    working_directory = getwd(),
    available_cores   = parallel::detectCores(),
    memory_limit      = mem_limit,
    session_info      = capture.output(sessionInfo())
  )

  return(sys_info)
}
