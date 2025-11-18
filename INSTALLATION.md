# Installation Guide - Satellite Cell Development Transcriptomics Pipeline

## Quick Start (Recommended)

```bash
# Clone the repository
git clone https://github.com/yourusername/satellite-cell-development-transcriptomics.git
cd satellite-cell-development-transcriptomics

# Run the automated installation
Rscript install_dependencies.R

# Run the complete pipeline
Rscript run_pipeline.R
```

## Detailed Installation

### Prerequisites

- **R ≥ 4.0.0** (recommended: 4.2.0+)
- **Bioconductor ≥ 3.14**
- **8GB+ RAM** (16GB recommended for large datasets)
- **Multi-core processor** (4+ cores recommended)
- **Internet connection** for package downloads

### Step 1: Install R and RStudio

**Ubuntu/Debian:**
```bash
sudo apt update
sudo apt install r-base r-base-dev
```

**macOS:**
```bash
brew install r
```

**Windows:** Download from [CRAN](https://cran.r-project.org/)

### Step 2: Install Bioconductor

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.16")
```

### Step 3: Install Required Packages

Run this R script to install all dependencies:

```r
# install_dependencies.R
cat("Installing required packages for Satellite Cell Development Pipeline...\n")

# CRAN packages
cran_packages <- c(
  "devtools", "remotes", "yaml", "logging", "jsonlite", "R.utils",
  "dplyr", "tidyr", "ggplot2", "ggrepel", "ggpubr", "viridis",
  "umap", "Rtsne", "phate", "destiny", "monocle3",
  "igraph", "ggraph", "tidygraph", "ggnetwork",
  "plotly", "DT", "shiny", "shinydashboard", "shinyWidgets",
  "ComplexHeatmap", "circlize", "RColorBrewer",
  "ggforce", "ggthemes", "gridExtra", "cowplot",
  "ggalluvial", "ggbeeswarm", "ggridges", "ggExtra"
)

# Install CRAN packages
for (pkg in cran_packages) {
  tryCatch({
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      cat(paste("Installed", pkg, "\n"))
    } else {
      cat(paste(pkg, "already available\n"))
    }
  }, error = function(e) {
    cat(paste("Failed to install", pkg, ":", e$message, "\n"))
    cat(paste("Installed", pkg, "\n"))
    } else {
      cat(paste(pkg, "already available\n"))
    }
  }, error = function(e) {
    cat(paste("Failed to install", pkg, ":", e$message, "\n"))
  })
}

# Bioconductor packages
bioc_packages <- c(
  "GEOquery", "Biobase", "limma", "affy", "preprocessCore",
  "DESeq2", "edgeR", "clusterProfiler", "org.Mm.eg.db",
  "DOSE", "enrichplot", "GSVA", "GSEABase", "pathview",
  "reactome.db", "msigdbr", "RcisTarget", "BSgenome.Mmusculus.UCSC.mm10",
  "TxDb.Mmusculus.UCSC.mm10.knownGene"
)

# Install Bioconductor packages
for (pkg in bioc_packages) {
  tryCatch({
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg)
      cat(paste("✓ Installed", pkg, "\n"))
    } else {
      cat(paste("✓", pkg, "already available\n"))
    }
  }, error = function(e) {
    cat(paste("✗ Failed to install", pkg, ":", e$message, "\n"))
  })
}

# Optional packages (may fail on some systems)
optional_packages <- c("GENIE3", "grndata", "SCENIC")

for (pkg in optional_packages) {
  tryCatch({
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (pkg == "GENIE3") {
        install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/GENIE3_1.20.0.tar.gz", repos = NULL, type = "source")
      } else {
        BiocManager::install(pkg)
      }
      cat(paste("✓ Installed", pkg, "\n"))
    } else {
      cat(paste("✓", pkg, "already available\n"))
    }
  }, error = function(e) {
    cat(paste("Optional package", pkg, "failed:", e$message, "\n"))
    cat(paste("  Pipeline will continue without", pkg, "functionality\n"))
  })
}

cat("\nInstallation complete\n")
cat("You can now run: Rscript run_pipeline.R\n")
```

### Step 4: Verify Installation

```r
# verification.R
cat("Verifying installation...\n\n")

# Check R version
cat("R Version:", R.version.string, "\n")

# Check available cores
cat("Available cores:", parallel::detectCores(), "\n")

# Check memory
cat("Memory limit:", if (exists("memory.limit")) memory.limit() else "Not available", "\n\n")

# Test package loading
test_packages <- c("ggplot2", "DESeq2", "clusterProfiler", "umap", "monocle3")

for (pkg in test_packages) {
  tryCatch({
    library(pkg, character.only = TRUE)
    cat(paste(pkg, "loaded successfully\n"))
  }, error = function(e) {
    cat(paste(pkg, "failed to load:", e$message, "\n"))
  })
}

cat("\n✅ Verification complete!\n")
```

## Docker Installation (Alternative)

```dockerfile
# Dockerfile
FROM rocker/rstudio:4.2.0

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages(c('devtools', 'BiocManager'))"
RUN R -e "BiocManager::install(c('GEOquery', 'Biobase', 'limma', 'DESeq2', 'clusterProfiler', 'org.Mm.eg.db', 'monocle3', 'umap', 'ComplexHeatmap'))"

# Copy pipeline
COPY . /home/rstudio/satellite-cell-development-transcriptomics
WORKDIR /home/rstudio/satellite-cell-development-transcriptomics

# Expose port
EXPOSE 8787

# Run pipeline
CMD ["Rscript", "run_pipeline.R"]
```

Build and run:
```bash
docker build -t satellite-cell-pipeline .
docker run -p 8787:8787 satellite-cell-pipeline
```

## Troubleshooting

### Common Issues

**1. Package Installation Fails**
```r
# Try installing from source
install.packages("package_name", type = "source")

# Or use BiocManager
BiocManager::install("package_name")
```

**2. Memory Issues**
```r
# Increase memory limit (Windows)
memory.limit(size = 16000)

# Monitor memory usage
gc()
```

**3. Parallel Processing Issues**
```r
# Set number of cores
options(mc.cores = 4)
```

**4. Bioconductor Package Issues**
```r
# Update Bioconductor
BiocManager::install(version = "3.16")

# Install specific version
BiocManager::install("package_name", version = "3.16")
```

### System-Specific Notes

**macOS:**
- Install Xcode command line tools: `xcode-select --install`
- May need to install gfortran for some packages

**Ubuntu/Debian:**
- Install build-essential: `sudo apt install build-essential`
- Install libxml2-dev, libcurl4-openssl-dev, libssl-dev

**Windows:**
- Install Rtools for package compilation
- Some packages may require manual installation

## Performance Optimization

### Memory Management
```r
# Configure for large datasets
options(java.parameters = "-Xmx8g")
Sys.setenv("R_MAX_VSIZE" = "16Gb")
```

### Parallel Processing
```r
# Enable parallel processing
library(doParallel)
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)
```

### Package Caching
```r
# Enable package caching
options(repos = c(CRAN = "https://cloud.r-project.org/"))
```

## Getting Help

### Documentation
- Full documentation: `docs/` directory
- Function help: `?function_name`
- Package vignettes: `browseVignettes("package_name")`

### Community Support
- GitHub Issues: Report bugs and feature requests
- Bioconductor Support: bioconductor.org/support/
- Stack Overflow: Tag with R and bioconductor

### Contact
- Email: research@institution.edu
- Issues: GitHub repository issues page

## Next Steps

1. **Run the pipeline**: `Rscript run_pipeline.R`
2. **Explore results**: Open `shiny_app/app.R` in RStudio
3. **Generate report**: Knit `research_report.Rmd`
4. **Customize analysis**: Edit `config.yaml`

## Citation

If you use this pipeline in your research, please cite:

```
Satellite Cell Development Transcriptomics Pipeline (2024)
MIT-Level Research Pipeline
https://github.com/yourusername/satellite-cell-development-transcriptomics
```

---

**Happy analyzing!**
