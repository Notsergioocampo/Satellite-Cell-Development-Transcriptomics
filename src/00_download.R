#!/usr/bin/env Rscript
# Satellite Cell Development Transcriptomics Pipeline
# Script 00: Data Download and Import
# MIT-Level Research Pipeline

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(yaml)
  library(logging)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# -----------------------------
# Helper functions (define FIRST)
# -----------------------------

extract_developmental_stage <- function(sample_names) {
  # Extract developmental stage (P1, P12, P28) from sample names
  stages <- str_extract(sample_names, "P[0-9]+")
  ifelse(is.na(stages), "unknown", stages)
}

extract_replicate_info <- function(sample_names) {
  # Extract replicate information from sample names (e.g. "rep1" or trailing digit)
  replicates <- str_extract(sample_names, "rep[0-9]+|[0-9]$")
  ifelse(is.na(replicates), "1", replicates)
}

validate_geo_data <- function(gse, sample_info, geo_id) {
  # Comprehensive data validation
  validation <- list(
    timestamp       = as.character(Sys.time()),
    dataset_id      = geo_id,
    platform        = annotation(gse),
    n_samples       = ncol(gse),
    n_features      = nrow(gse),
    n_stages        = length(unique(sample_info$stage)),
    stages_present  = unique(sample_info$stage),
    missing_values  = sum(is.na(exprs(gse))),
    expression_range = range(exprs(gse), na.rm = TRUE),
    qc_pass         = TRUE,
    warnings        = character(0)
  )

  # Check for required developmental stages
  required_stages <- c("P1", "P12", "P28")
  missing_stages  <- setdiff(required_stages, validation$stages_present)

  if (length(missing_stages) > 0) {
    validation$warnings <- c(
      validation$warnings,
      paste("Missing developmental stages:",
            paste(missing_stages, collapse = ", "))
    )
    validation$qc_pass <- FALSE
  }

  # Check data quality metrics
  if (validation$missing_values > 0) {
    validation$warnings <- c(
      validation$warnings,
      "Dataset contains missing values"
    )
  }

  validation
}

# -----------------------------
# Logging configuration
# -----------------------------

log_dir <- file.path("results", "logs")
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)

logfile <- file.path(log_dir, "00_download.log")

# Initialize logging: basicConfig to console + file handler
basicConfig()  # sets root logger at INFO to console
addHandler(writeToFile, file = logfile, level = "INFO")

loginfo("Starting satellite cell transcriptomics data download pipeline")

# -----------------------------
# Load configuration
# -----------------------------

if (!file.exists("config.yaml")) {
  logerror("config.yaml not found in project root.")
  stop("config.yaml not found. Make sure you run this script from the project root.")
}

config   <- yaml::read_yaml("config.yaml")
geo_id   <- config$dataset$geo_id
platform <- config$dataset$platform

loginfo(paste("Configured GEO dataset:", geo_id))
loginfo(paste("Expected platform:", platform))

# -----------------------------
# Create data directories
# -----------------------------

dirs <- c("data/raw", "data/processed", "data/metadata")
for (dir in dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    loginfo(paste("Created directory:", dir))
  }
}

# -----------------------------
# Download GEO dataset
# -----------------------------

tryCatch({
  loginfo("Connecting to GEO database via GEOquery::getGEO ...")
  gse <- getGEO(geo_id, GSEMatrix = TRUE, AnnotGPL = TRUE)

  if (length(gse) == 0) {
    stop("No data retrieved from GEO for ID: ", geo_id)
  }

  # Many GSEs return a list; use the first ExpressionSet
  gse <- gse[[1]]

  loginfo(paste("Successfully downloaded dataset:", geo_id))
  loginfo(paste("Platform (from ExpressionSet annotation):", annotation(gse)))
  loginfo(paste("Number of samples:", ncol(gse)))
  loginfo(paste("Number of features:", nrow(gse)))

  # Save raw data
  raw_file <- file.path("data", "raw", paste0(geo_id, "_raw.RData"))
  save(gse, file = raw_file)
  loginfo(paste("Raw ExpressionSet saved to:", raw_file))

  # Extract and save phenotype data
  pheno_data <- pData(gse)
  pheno_file <- file.path("data", "metadata", paste0(geo_id, "_phenotype.csv"))
  write.csv(pheno_data, file = pheno_file, row.names = TRUE)
  loginfo(paste("Phenotype data saved to:", pheno_file))

  # Extract and save feature data
  feature_data <- fData(gse)
  feature_file <- file.path("data", "metadata", paste0(geo_id, "_features.csv"))
  write.csv(feature_data, file = feature_file, row.names = TRUE)
  loginfo(paste("Feature data saved to:", feature_file))

  # Process sample information
  sample_ids <- colnames(gse)
  sample_info <- data.frame(
    sample_id = sample_ids,
    stage     = extract_developmental_stage(sample_ids),
    replicate = extract_replicate_info(sample_ids),
    stringsAsFactors = FALSE
  )

  sample_file <- file.path("data", "metadata", "sample_info.csv")
  write.csv(sample_info, file = sample_file, row.names = FALSE)
  loginfo(paste("Sample information saved to:", sample_file))

  # Validate data integrity
  validation_results <- validate_geo_data(gse, sample_info, geo_id)
  validation_file    <- file.path("results", "data_validation.json")
  if (!dir.exists("results")) dir.create("results", recursive = TRUE)
  jsonlite::write_json(validation_results, validation_file,
                       auto_unbox = TRUE, pretty = TRUE)
  loginfo(paste("Data validation completed. Report saved to:", validation_file))

}, error = function(e) {
  msg <- paste("Error downloading GEO data:", conditionMessage(e))
  logerror(msg)
  stop(msg)
})

loginfo("Data download pipeline completed successfully")
