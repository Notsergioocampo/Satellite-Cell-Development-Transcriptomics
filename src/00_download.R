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

# Configure logging
logfile <- file.path("results", "logs", "00_download.log")
if (!dir.exists(dirname(logfile))) dir.create(dirname(logfile), recursive = TRUE)
basicConfig(level = "INFO", file = logfile)

loginfo("Starting satellite cell transcriptomics data download pipeline")

# Load configuration
config <- yaml::read_yaml("config.yaml")
geo_id <- config$dataset$geo_id
platform <- config$dataset$platform

loginfo(paste("Downloading GEO dataset:", geo_id))

# Create data directories
dirs <- c("data/raw", "data/processed", "data/metadata")
for (dir in dirs) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# Download GEO dataset with error handling
tryCatch({
  loginfo("Connecting to GEO database...")
  gse <- getGEO(geo_id, GSEMatrix = TRUE, AnnotGPL = TRUE)
  
  if (length(gse) > 0) {
    gse <- gse[[1]]  # Extract first dataset
    
    loginfo(paste("Successfully downloaded dataset:", geo_id))
    loginfo(paste("Platform:", annotation(gse)$platform))
    loginfo(paste("Samples:", ncol(gse)))
    loginfo(paste("Features:", nrow(gse)))
    
    # Save raw data
    raw_file <- file.path("data/raw", paste0(geo_id, "_raw.RData"))
    save(gse, file = raw_file)
    loginfo(paste("Raw data saved to:", raw_file))
    
    # Extract and save phenotype data
    pheno_data <- pData(gse)
    pheno_file <- file.path("data/metadata", paste0(geo_id, "_phenotype.csv"))
    write.csv(pheno_data, file = pheno_file, row.names = TRUE)
    loginfo(paste("Phenotype data saved to:", pheno_file))
    
    # Extract and save feature data
    feature_data <- fData(gse)
    feature_file <- file.path("data/metadata", paste0(geo_id, "_features.csv"))
    write.csv(feature_data, file = feature_file, row.names = TRUE)
    loginfo(paste("Feature data saved to:", feature_file))
    
    # Process sample information
    sample_info <- data.frame(
      sample_id = colnames(gse),
      stage = extract_developmental_stage(colnames(gse)),
      replicate = extract_replicate_info(colnames(gse)),
      stringsAsFactors = FALSE
    )
    
    sample_file <- file.path("data/metadata", "sample_info.csv")
    write.csv(sample_info, file = sample_file, row.names = FALSE)
    loginfo(paste("Sample information saved to:", sample_file))
    
    # Validate data integrity
    validation_results <- validate_geo_data(gse, sample_info)
    validation_file <- file.path("results", "data_validation.json")
    jsonlite::write_json(validation_results, validation_file, auto_unbox = TRUE, pretty = TRUE)
    loginfo("Data validation completed")
    
  } else {
    stop("No data retrieved from GEO")
  }
  
}, error = function(e) {
  logerror(paste("Error downloading GEO data:", e$message))
  stop(e)
})

loginfo("Data download pipeline completed successfully")

# Helper functions
extract_developmental_stage <- function(sample_names) {
  # Extract developmental stage (P1, P12, P28) from sample names
  stages <- str_extract(sample_names, "P[0-9]+")
  return(ifelse(is.na(stages), "unknown", stages))
}

extract_replicate_info <- function(sample_names) {
  # Extract replicate information from sample names
  replicates <- str_extract(sample_names, "rep[0-9]+|[0-9]$")
  return(ifelse(is.na(replicates), "1", replicates))
}

validate_geo_data <- function(gse, sample_info) {
  # Comprehensive data validation
  validation <- list(
    timestamp = Sys.time(),
    dataset_id = geo_id,
    platform = annotation(gse)$platform,
    n_samples = ncol(gse),
    n_features = nrow(gse),
    n_stages = length(unique(sample_info$stage)),
    stages_present = unique(sample_info$stage),
    missing_values = sum(is.na(exprs(gse))),
    expression_range = range(exprs(gse), na.rm = TRUE),
    qc_pass = TRUE
  )
  
  # Check for required developmental stages
  required_stages <- c("P1", "P12", "P28")
  missing_stages <- setdiff(required_stages, validation$stages_present)
  
  if (length(missing_stages) > 0) {
    validation$warnings <- paste("Missing developmental stages:", paste(missing_stages, collapse = ", "))
    validation$qc_pass <- FALSE
  }
  
  # Check data quality metrics
  if (validation$missing_values > 0) {
    validation$warnings <- c(validation$warnings, "Dataset contains missing values")
  }
  
  return(validation)
}
