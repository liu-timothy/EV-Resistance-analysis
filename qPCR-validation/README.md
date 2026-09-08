# qPCR validation

Inputs, notebooks, and figures used to compare qPCR measurements with bulk RNA-seq fold changes.

## Folders

| Folder | Contents |
|---|---|
| `qpcr_data/` | Instrument-export CSVs containing sample, target, Ct, and quality-control fields |
| `bulkRNAseq_data/` | Reduced, CSV-formatted DESeq2 tables for 647V, RT112, and UMUC1 comparisons |
| `scripts/` | qPCR calculation/QC notebooks; also contains some figures generated during notebook development |
| `outs/` | Generated qPCR-versus-RNA-seq plots |

`outs/` groups many candidate-gene figures by source and expected direction:

- `crispr_upregulated/` and `crispr_downregulated/`: candidates originating from the CRISPR screen.
- `tim_upregulated/` and `tim_downregulated/`: candidate lists analyzed by Tim.
- `will_upregulated/` and `will_downregulated/`: candidate lists analyzed by Will.
- Files directly under `outs/`: cell-line comparisons and selected individual genes.

## Notebooks

- `scripts/qpcr-plots.ipynb`: calculates qPCR relative expression and generates qPCR-versus-RNA-seq plots for parental/resistant and CRISPR-related comparisons.
- `scripts/Ct_basemean_QC_check.ipynb`: compares qPCR Ct behavior with RNA-seq base mean as a quality-control check.

Both notebooks currently reference an additional redo dataset by an absolute path in `/Users/timothy/Downloads/`; that file is not in this repository. Move a versioned copy into `qpcr_data/` and update the path for portable reruns.

## Input conventions

qPCR CSVs use fields such as `Sample Name`, `Target Name`, `CT`, `Ct Mean`, and instrument QC flags (`HIGHSD`, `NOAMP`, `EXPFAIL`). The exact set of columns varies between exports. Bulk RNA-seq CSVs contain `baseMean`, `log2FoldChange`, `lfcSE`, `pvalue`, `padj`, `GeneSymbol`, and `gene_type`.

