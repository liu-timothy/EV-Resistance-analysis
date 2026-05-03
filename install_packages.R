# Set CRAN mirror - use HTTPS with proper URL
options(repos = c(CRAN = "https://cran.r-project.org"))

# Install BiocManager from CRAN first
install.packages("BiocManager")

# Load BiocManager
library(BiocManager)

# List of packages
pkg <- c("DESeq2", "msigdbr", "stringr", "dplyr", "fgsea", "ggplot2", "ggrepel", 
         "gg.gap", "ggbreak", "openxlsx", "tibble", "tidyr", "pheatmap",
         "ComplexHeatmap", "circlize", "RColorBrewer", "grid", "gridExtra", "Gviz", 
         "GenomicFeatures", "reshape2", "rtracklayer", "M3C", "BiocManager",
         "ggpubr", "glue", "gprofiler2", "AnnotationDbi", "org.Hs.eg.db", "ggvenn",
         "clusterProfiler", "forcats", "patchwork", "ggpattern", "purrr", "ComplexUpset",
         "ReactomePA", "Sushi")

# Install all packages using BiocManager
BiocManager::install(pkg, update = FALSE, ask = FALSE)

# Verify installations
installed <- rownames(installed.packages())
missing <- pkg[!pkg %in% installed]
if (length(missing) > 0) {
    cat("The following packages failed to install:\n")
    print(missing)
} else {
    cat("All packages installed successfully!\n")
}