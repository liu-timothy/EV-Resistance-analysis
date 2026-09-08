# EV Resistance Analysis

Analysis workspace for comparing extracellular-vesicle (EV) resistant bladder cancer cell lines with their parental lines, identifying shared RNA-seq signals, comparing those signals with a CRISPR screen, and validating selected genes by qPCR.

The repository contains analysis code, processed differential-expression results, and generated figures. It does **not** currently contain the raw sequencing inputs referenced by the R pipeline.

## Project map

```text
EV-Resistance-analysis/
├── R-env/                       R environment setup and bulk RNA-seq pipeline
├── data/                        DESeq2/CRISPR tables and RNA-seq figures
│   ├── rp_vs_rr/                RT112 parental vs resistant
│   ├── up_vs_ur/                UMUC1 parental vs resistant
│   └── vp_vs_vr/                647V parental vs resistant
├── bulkRNAseq_gene-overlaps/    Notebooks for gene-list overlap and prioritization
├── qPCR-validation/             qPCR inputs, matched RNA-seq tables, notebooks, figures
├── output.png                   Legacy/generated top-level figure (provenance unspecified)
└── README.md                    This guide
```

More detailed inventories are available in [`data/README.md`](data/README.md), [`bulkRNAseq_gene-overlaps/README.md`](bulkRNAseq_gene-overlaps/README.md), [`qPCR-validation/README.md`](qPCR-validation/README.md), and [`R-env/README.md`](R-env/README.md).

## Naming conventions

Comparison names encode the cell line and phenotype:

| Code | Meaning |
|---|---|
| `v` | 647V |
| `r` | RT112 |
| `u` | UMUC1 |
| `p` | parental |
| `r` (second position) | resistant |

Thus, `vp_vs_vr` means **647V parental vs 647V resistant**, `rp_vs_rr` means **RT112 parental vs RT112 resistant**, and `up_vs_ur` means **UMUC1 parental vs UMUC1 resistant**. In DESeq2 result names, `p_vs_r` indicates that fold changes describe the second condition relative to the first; confirm the contrast in the producing analysis before interpreting the sign.

Known resistant derivatives represented in the qPCR files are HT8 (UMUC1), 2F2 early/late (RT112), and 1C2 (647V).

## Analysis workflow

1. The R script performs count QC, PCA, DESeq2, fold-change shrinkage, GSEA/ORA, and overlap analyses.
2. Processed results and plots are stored under `data/<comparison>/`.
3. Python notebooks compare genes across cell lines and against the EV-vs-PBS CRISPR screen.
4. qPCR notebooks calculate/plot expression changes and compare them with bulk RNA-seq results.

## Getting started

### Python notebooks

The checked-in Conda environment contains Python 3.10, Jupyter, pandas, NumPy, matplotlib, seaborn, scikit-learn, scipy, `matplotlib-venn`, and `nbqol`.

```bash
conda env create -f R-env/environment.yml
conda activate python-env
jupyter lab
```

Open notebooks from `bulkRNAseq_gene-overlaps/` or `qPCR-validation/scripts/`. They use `nbqol.path_to_git_root` to locate this repository.

### R analysis

Install the required CRAN/Bioconductor packages in an R environment:

```bash
Rscript R-env/install_packages.R
```

Then edit every location marked `(EDIT)` in `R-env/RNAseq Individual Analysis - QC, DESeq2, GSEA.r`. The current script expects external feature-count, TPM, metadata, annotation, and plotting-function files at a Windows/OneDrive path; those inputs are not included here.

## File types

| Extension | Role |
|---|---|
| `.R` / `.r` | Package installation and RNA-seq analysis code |
| `.ipynb` | Exploratory Python analyses and plot generation |
| `.csv` | Comma-separated qPCR inputs or combined/trimmed RNA-seq tables |
| `.txt` | Tab-separated DESeq2/CRISPR results, gene lists, or run statistics |
| `.xlsx` | Spreadsheet copies of DESeq2 result tables |
| `.pdf` | RNA-seq MA, volcano, and GSEA plots |
| `.png` | qPCR/RNA-seq comparison figures and exploratory outputs |
| `.yml` | Reproducible Conda environment specification |

## Reproducibility notes

- `environment.yml` was exported with Linux build pins and an absolute `prefix`; recreating it on macOS or Windows may require removing build strings and the final `prefix` entry.
- `Ct_basemean_QC_check.ipynb` and `qpcr-plots.ipynb` reference a CSV in `/Users/timothy/Downloads/`. Place that source file in `qPCR-validation/qpcr_data/` and update the notebook path before a clean rerun.
- Notebooks contain saved outputs. Restart the kernel and run all cells when generating results for reporting.
- Raw RNA-seq data are not versioned in this repository. Record their source, genome/annotation version, and checksums separately when publishing or transferring the project.

## Contribution and organization guidelines

- Put reusable analysis code in `R-env/` or the appropriate notebook folder; put generated results under the matching output folder.
- Use lowercase, underscore-separated names for new files and avoid spaces where practical.
- Keep raw inputs separate from derived outputs. For new analyses, prefer `data/raw/`, `data/processed/`, and `results/` rather than mixing them.
- Add a short Markdown cell at the top of each new notebook stating its purpose, inputs, outputs, and filtering thresholds.
- Do not commit local session files (`.RData`, `.Rhistory`), OS metadata, notebook checkpoints, or virtual environments.

