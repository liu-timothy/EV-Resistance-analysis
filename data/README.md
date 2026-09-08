# Data and bulk RNA-seq results

This directory holds processed differential-expression outputs for three parental-versus-resistant cell-line comparisons plus results from an EV-vs-PBS CRISPR screen.

## Comparison folders

| Folder | Comparison |
|---|---|
| `rp_vs_rr/` | RT112 parental vs RT112 resistant (2F2) |
| `up_vs_ur/` | UMUC1 parental vs UMUC1 resistant (HT8) |
| `vp_vs_vr/` | 647V parental vs 647V resistant (1C2) |

Each comparison folder follows the same pattern, where `<comparison>` is the folder name:

| Pattern | Contents |
|---|---|
| `<comparison>_summary.txt` / `.xlsx` | Full DESeq2 results using raw log2 fold changes |
| `<comparison>_lfcs_summary.txt` / `.xlsx` | Results using shrunken log2 fold changes (`lfcs`) |
| `<comparison>_up_genes.txt` | Significant upregulated gene list from raw results |
| `<comparison>_down_genes.txt` | Significant downregulated gene list from raw results |
| `<comparison>_lfcs_up_genes.txt` | Upregulated genes selected using shrunken fold changes |
| `<comparison>_lfcs_down_genes.txt` | Downregulated genes selected using shrunken fold changes |
| `<comparison>_stats.txt` | DESeq2 filtering/significance summary |
| `* MA Plot *.pdf` | MA plots for raw or shrunken fold changes |
| `* Volcano Plot *.pdf` | Volcano plots for raw or shrunken fold changes |
| `* GSEA Plot.pdf` | Gene-set enrichment visualization |

Typical result columns are `baseMean`, `log2FoldChange`, `lfcSE`, `stat` (raw summaries only), `pvalue`, `padj`, `GeneSymbol`, and `gene_type`. Rows in the shrunken text tables retain Ensembl IDs as row names.

## Top-level files

- `combined-lines_summary.csv`: raw DESeq2 summaries combined across RT112, UMUC1, and 647V, with a `line` column identifying the cell line.
- `EV_vs_PBS_DESeq2_results (Screen DESEQ2).txt`: differential-expression results from the CRISPR-screen EV-vs-PBS comparison, used by `EV_CRISPR_overlap_analysis.ipynb`.

These are processed analysis products, not raw sequencing reads or feature-count matrices.

