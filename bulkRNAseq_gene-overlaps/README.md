# Bulk RNA-seq gene-overlap analyses

Python/Jupyter notebooks for finding shared differentially expressed genes and prioritizing candidates.

| Notebook | Purpose | Primary inputs |
|---|---|---|
| `tim_analysis_overlapping-genes.ipynb` | Combines the three cell-line DESeq2 tables, selects top up/down genes, and explores pairwise overlaps | `data/{rp_vs_rr,up_vs_ur,vp_vs_vr}/*_summary.txt` |
| `final_T10-genes.ipynb` | Checks shortlisted/top genes for fold-change magnitude and statistical significance | `data/combined-lines_summary.csv` |
| `EV_CRISPR_overlap_analysis.ipynb` | Compares bulk RNA-seq hits with the EV-vs-PBS CRISPR-screen results and explores gene overlap | CRISPR result table plus `data/rp_vs_rr/rp_vs_rr_summary.txt` |

The notebooks use `nbqol.path_to_git_root` to construct repository-relative paths. They are exploratory records and include saved cell outputs; rerun all cells in order before relying on regenerated results.

