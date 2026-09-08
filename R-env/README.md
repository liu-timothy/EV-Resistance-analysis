# Analysis environments and R pipeline

## Files

- `RNAseq Individual Analysis - QC, DESeq2, GSEA.r`: bulk RNA-seq workflow covering count QC, PCA, DESeq2, MA/volcano plots, Hallmark GSEA, GO/KEGG/Reactome enrichment, and intersection analyses.
- `install_packages.R`: installs the R/Bioconductor packages required by the broader workflow and reports packages that failed to install.
- `environment.yml`: Conda environment for the Python/Jupyter notebooks, despite this directory's historical `R-env` name.

## Running the R workflow

```bash
Rscript R-env/install_packages.R
Rscript "R-env/RNAseq Individual Analysis - QC, DESeq2, GSEA.r"
```

Before running, edit the paths marked `(EDIT)`. The script currently expects external count/TPM matrices, metadata, GENCODE 28 annotation, and `functions/plot_functions.R`; these files are not committed here. It also uses a Windows-specific OneDrive base path.

The script's sample-code convention is `<cell line><phenotype>_<replicate>`: `v`, `r`, or `u` for 647V, RT112, or UMUC1; followed by `p` or `r` for parental or resistant.

## Environment caveat

`environment.yml` is a fully pinned Linux export with a machine-specific `prefix`. It is useful as an exact provenance snapshot. For a portable environment, consider maintaining a smaller `environment-portable.yml` containing only direct dependencies and no build strings or prefix.

