# EV Resistance Analysis

This project uses bulk RNA-seq and qPCR data to identify and validate genes associated with Enfortumab Vedotin (EV) sensitivity in bladder cancer cell lines. Candidate genes are selected using 3 approaches: pathway enrichment, overlap across cell lines, and overlap with a CRISPR screen.

## Project structure

``` text
EV-Resistance-analysis/
├── identifying-hits/       RNA-seq data and notebooks for selecting candidate genes
├── qPCR-validation/        qPCR validation data, notebooks, and output figures
├── R-env/                  R script for the original RNA-seq analysis
├── environment.yml         Conda environment for the Python notebooks
└── README.md               Project overview and file guide
```

### `identifying-hits/`

Contains the three approaches used to identify candidate EV-sensitivity genes.

#### `bulkRNAseq_gene-overlaps/`

| Notebook | Purpose |
|------------------------------------|------------------------------------|
| `hits_by_pathway_enrichment.ipynb` | Prioritizes differentially expressed genes found in relevant enriched pathways |
| `hits_by_cross-cell-line_overlap.ipynb` | Finds genes shared across the RT112, UMUC1, and 647V resistant cell-line comparisons |
| `hits_by_crispr_rnaseq_overlap.ipynb` | Compares RNA-seq hits with genes identified in the EV-vs-PBS CRISPR screen |

#### `data/`

Contains processed DESeq2 results and figures used by the hit-identification notebooks.

| Path | Description |
|------------------------------------|------------------------------------|
| `rp_vs_rr/` | RT112 parental vs resistant (2F2) |
| `up_vs_ur/` | UMUC1 parental vs resistant (HT8) |
| `vp_vs_vr/` | 647V parental vs resistant (1C2) |
| `combined-lines_summary.*` | Combined differential-expression results from all three cell lines |
| `EV_vs_PBS_DESeq2_results (Screen DESEQ2).txt` | Differential-expression results from the CRISPR screen |

Each cell-line folder includes full and fold-change-shrunken DESeq2 tables, upregulated/downregulated gene lists, summary statistics, and MA, volcano, and GSEA plots.

### `qPCR-validation/`

Contains the experimental data and analyses used to validate selected candidate genes.

| Folder | Description |
|------------------------------------|------------------------------------|
| `qpcr_data/` | Raw qPCR instrument-export CSV files |
| `bulkRNAseq_data/` | Reduced RNA-seq tables matched to the qPCR comparisons |
| `scripts/` | Notebooks comparing qPCR Ct values and expression changes with RNA-seq results |
| `outs/` | Generated qPCR-versus-RNA-seq plots, grouped by hit source and expression direction |

The two main validation notebooks are:

- `compare_qpcr_rnaseq_log2fc.ipynb`: compares qPCR relative expression with RNA-seq log2 fold change.
- `compare_qpcr_ct_rnaseq_basemean.ipynb`: compares qPCR Ct values with RNA-seq base mean as a QC check.

### `R-env/`

Contains `RNAseq Individual Analysis - QC, DESeq2, GSEA.r`, the original R workflow for count QC, PCA, differential-expression analysis, fold-change shrinkage, enrichment analysis, and result visualization. The script references external raw counts, metadata, annotations, and helper functions that are not included in this repository.
