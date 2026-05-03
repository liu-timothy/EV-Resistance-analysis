# Adapted by William  - Chou Lab, 2026
# To edit this script, CTRL+F "(EDIT)" and replace with appropriate code.

# Script Functions:
# 1. Gene Count Quality Control
# 2. PCA
# 3. DESeq2 Analysis (Volcano, MA Plots)
# 4. GSEA - Hallmarks (Barplots)
# 5. GSEA and ORA - GO, KEGG, Reactome (Master Table)
# 6. Intersection Analysis (Upset Plots, Dot Plots, BP GO Plots)

# Note: GSEA and ORA often give similar results, but GSEA can be particularly good when changes are weak (either as a result of poor sample quality,
# weak biological phenomenon, etc.). ORA is primarily dictated by how you select significant genes (established in your p-adjusted and log2FC threshold).

# Color Scheme (colorblind-friendly):
#D81B60
#1E88E5
#FFC107
#004D40

# packages
pkg <- c("DESeq2", "msigdbr", "stringr", "dplyr", "fgsea", "ggplot2", "ggrepel", 
          "gg.gap", "ggbreak", "openxlsx", "tibble", "tidyr", "pheatmap",
          "ComplexHeatmap", "circlize", "RColorBrewer", "grid", "gridExtra", "Gviz", 
          "GenomicFeatures", "reshape2", "rtracklayer", "M3C", "BiocManager",
          "ggpubr", "glue", "gprofiler2", "AnnotationDbi", "org.Hs.eg.db", "ggvenn",
          "clusterProfiler", "forcats", "patchwork", "ggpattern", "purrr", "ComplexUpset",
          "ReactomePA", "Sushi")

# excluded packages:         

pkg <- c("DESeq2", "ggplot2", "ComplexHeatmap")

for (x in pkg) {
  library(x, character.only = TRUE)
}

# directories
dir_base <- "C:/Users/Will/OneDrive - UCSF/" # (EDIT)
input_base_dir <- paste0(dir_base, "RNAseq 2026_04_23_WD_EV_Resistance_Analysis/") # (EDIT)
output_base_dir <- paste0(dir_base, "RNAseq 2026_04_23_WD_EV_Resistance_Analysis/output_files/") # (EDIT): leave the /output_files/ at the end
de_dir <- paste0(output_base_dir, "de/")
dir.create(de_dir)
fc_dir <- paste0(input_base_dir, "featurecounts_gencode28/")
metadata_dir <- paste0(input_base_dir,"metadata/")
fn_gene_info <- paste0(input_base_dir, "info/gencode28_ensembl2sym.txt")
fn_samples_647v <- paste0(metadata_dir, "metadata_647v.txt")
fn_gene_counts_647v <- paste(fc_dir, "counts_gene_id_647v.tsv", sep = "")
fn_gene_tpm_647v <- paste(fc_dir, "counts_tpm_gene_id_647v.tsv", sep = "")
fn_samples_rt112 <- paste0(metadata_dir, "metadata_rt112.txt")
fn_gene_counts_rt112 <- paste(fc_dir, "counts_gene_id_rt112.tsv", sep = "")
fn_gene_tpm_rt112 <- paste(fc_dir, "counts_tpm_gene_id_rt112.tsv", sep = "")
fn_samples_umuc1 <- paste0(metadata_dir, "metadata_umuc1.txt")
fn_gene_counts_umuc1 <- paste(fc_dir, "counts_gene_id_umuc1.tsv", sep = "")
fn_gene_tpm_umuc1 <- paste(fc_dir, "counts_tpm_gene_id_umuc1.tsv", sep = "")
fn_plot_functions <- paste0(input_base_dir, "functions/plot_functions.R")

# gene info and pathways
gene_info <- read.delim(fn_gene_info, header = T, stringsAsFactors = F, strip.white = T)
pathways.hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% split(x = .$gene_symbol, f = .$gs_name)

# remove mitochondrial RNA, rRNA, and rRNA pseudogenes
genes_to_keep <- gene_info$ensembl_ID[which(!gene_info$chr %in% c("chrM") &
                                              !gene_info$type %in% c("rRNA", "rRNA_pseudogene", "Mt_rRNA"))]

# sample naming scheme:
# character 1:
  # v = 647v
  # r = rt112
  # u = umuc1
# character 2: 
  # p = parental
  # r = resistant
# character 3, 4: 
  # _1, _2, _3 = replicate 1, 2, or 3

# dataset cleaning function
process_dataset <- function(count_file,
                            tpm_file,
                            metadata_file,
                            genes_to_keep,
                            sample_filter_fn = NULL,
                            rename_fn,
                            batch_label) {
  
  # load tpm
  tpm <- read.delim(tpm_file, sep = "\t", header = TRUE,
                    stringsAsFactors = FALSE, check.names = FALSE)
  rownames(tpm) <- tpm$feature_id
  tpm$feature_id <- NULL
  
  # load counts
  counts <- read.delim(count_file, sep = "\t", header = TRUE,
                       stringsAsFactors = FALSE, check.names = FALSE)
  rownames(counts) <- counts$feature_id
  counts$feature_id <- NULL
  
  # gene filtering (removing chrM, rRNA, rRNA pseudogenes, Mt rRNA)
  genes_present <- intersect(genes_to_keep, rownames(counts))
  counts <- counts[genes_present, ]
  tpm    <- tpm[genes_present, ]
  
  # metadata
  samples <- read.delim(metadata_file, header = TRUE,
                        stringsAsFactors = FALSE)
  
  # dataset-specific filtering (optional)
  if (!is.null(sample_filter_fn)) {
    samples <- sample_filter_fn(samples)
  }
  
  # keep original IDs
  samples$sample_id_old <- samples$sample_id
  
  # generate new sample names
  samples$sample_id_new <- rename_fn(samples$sample_id_old)
  
  # align samples
  common_samples <- intersect(samples$sample_id_old, colnames(counts))
  
  counts <- counts[, common_samples, drop = FALSE]
  tpm    <- tpm[, common_samples, drop = FALSE]
  
  samples <- samples[samples$sample_id_old %in% common_samples, ]
  
  # rename matrices
  name_map <- setNames(samples$sample_id_new, samples$sample_id_old)
  
  colnames(counts) <- name_map[colnames(counts)]
  colnames(tpm)    <- name_map[colnames(tpm)]
  
  # finalize metadata
  samples$sample_id <- samples$sample_id_new
  
  # add treatment and batch to metadata
  samples$treatment <- factor(
    ifelse(substr(samples$sample_id, 2, 2) == "p", "Parental", "Resistant"),
    levels = c("Parental", "Resistant")
  )

  samples$batch <- factor(batch_label)

  # enforce ordering
  samples <- samples[match(colnames(counts), samples$sample_id), ]
  
  # sanity checks
  stopifnot(all(colnames(counts) == samples$sample_id))
  stopifnot(sum(is.na(counts)) == 0)
  
  return(list(
    counts = counts,
    tpm = tpm,
    samples = samples
  ))
}

rename_647v <- function(ids) {
  paste0(
    "v",
    ifelse(substr(ids, 5, 5) == "1", "p", "r"),
    "_",
    substr(ids, 7, 7)
  )
}

filter_rt112 <- function(df) {
  df[grepl("^RT-[12]-", df$sample_id), ]
}

rename_rt112 <- function(ids) {
  parts <- strsplit(ids, "-")
  sapply(parts, function(x) {
    paste0(
      "r",
      ifelse(x[2] == "1", "p", "r"),
      "_",
      x[3]
    )
  })
}

filter_umuc1 <- function(df) {
  df[grepl("^(UM1|HT8)_", df$sample_id), ]
}

rename_umuc1 <- function(ids) {
  parts <- strsplit(ids, "_")
  sapply(parts, function(x) {
    paste0(
      "u",
      ifelse(x[1] == "UM1", "p", "r"),
      "_",
      x[2]
    )
  })
}

res_647v <- process_dataset(
  fn_gene_counts_647v,
  fn_gene_tpm_647v,
  fn_samples_647v,
  genes_to_keep,
  sample_filter_fn = NULL,
  rename_fn = rename_647v,
  batch_label = "1"
)

res_rt112 <- process_dataset(
  fn_gene_counts_rt112,
  fn_gene_tpm_rt112,
  fn_samples_rt112,
  genes_to_keep,
  sample_filter_fn = filter_rt112,
  rename_fn = rename_rt112,
  batch_label = "2"
)

res_umuc1 <- process_dataset(
  fn_gene_counts_umuc1,
  fn_gene_tpm_umuc1,
  fn_samples_umuc1,
  genes_to_keep,
  sample_filter_fn = filter_umuc1,
  rename_fn = rename_umuc1,
  batch_label = "3"
)

master_counts <- list(
  res_647v$counts,
  res_rt112$counts,
  res_umuc1$counts
)

master_samples <- list(
  res_647v$samples,
  res_rt112$samples,
  res_umuc1$samples
)

# union + zero fill
all_genes <- Reduce(union, lapply(master_counts, rownames))

master_counts <- lapply(master_counts, function(mat) {
  missing <- setdiff(all_genes, rownames(mat))
  if (length(missing) > 0) {
    zero_mat <- matrix(0, nrow = length(missing), ncol = ncol(mat),
                       dimnames = list(missing, colnames(mat)))
    mat <- rbind(mat, zero_mat)
  }
  mat[all_genes, ]
})

combined_counts <- do.call(cbind, master_counts)
combined_samples <- do.call(rbind, master_samples)

# enforce alignment
combined_samples <- combined_samples[
  match(colnames(combined_counts), combined_samples$sample_id), ]

samples_info <- combined_samples
View(samples_info)
count_matrix <- as.matrix(combined_counts)
View(count_matrix)

######################
### Gene Counts QC ###
######################
fn_hist_gene_counts_mean <- paste0(de_dir, "Gene Counts Distribution (All Samples).png")
png(fn_hist_gene_counts_mean, width = 600, height = 400)
hist(log2(rowMeans(count_matrix) + 1), breaks = 10000, xlab = "log2(gene counts)", main = "")
dev.off()

######################
###### PCA Plot ######
######################
run_pca_plot <- function(tpm_matrix, samples_info, dataset_name, percentile = 10, outdir) {
  
  samples_to_use <- colnames(tpm_matrix)
  
  # log transform
  mat <- log2(tpm_matrix + 1)
  mat <- mat[rowSums(mat > 0) > 0, samples_to_use, drop = FALSE]
  
  # align metadata
  cluster_samples <- samples_info[
    samples_info$sample_id %in% samples_to_use, ]
  cluster_samples <- cluster_samples[
    match(samples_to_use, cluster_samples$sample_id), ]
  
  # feature selection
  mat_for_pca <- M3C::featurefilter(
    mydata = mat,
    percentile = percentile,
    method = "MAD"
  )
  
  mat_for_cluster <- t(mat_for_pca[[1]])
  
  # PCA
  pca <- prcomp(mat_for_cluster)
  
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)
  percentVar <- percentVar[1:3]
  
  d_pca <- data.frame(
    PC1 = pca$x[,1],
    PC2 = pca$x[,2],
    PC3 = pca$x[,3]
  )
  
  d_pca <- cbind(d_pca, cluster_samples)
  
  # plot
  p <- ggplot(d_pca, aes(x = PC1, y = PC2)) +
    geom_point(aes(fill = treatment, color = treatment), size = 3) +
    geom_text_repel(aes(label = sample_id), size = 3, max.overlaps = 50) +
    theme(
      text = element_text(size = 15),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      strip.background = element_rect(fill = "#FFFFFF", color = "#FFFFFF"),
      panel.background = element_rect(fill = "#FFFFFF", color = "#000000")
    ) +
    ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
    ggtitle(paste0("PCA - ", dataset_name))
  
  print(p)
  
  # save
  fn <- paste0(outdir,
               sprintf("PCA_%s_top%spercent.png", dataset_name, percentile))
  
  png(fn, width = 600, height = 600)
  print(p)
  dev.off()
  
  return(p)
}

run_pca_plot(res_647v$tpm,  res_647v$samples,  "647v",  percentile = 10, outdir = de_dir)
run_pca_plot(res_rt112$tpm, res_rt112$samples, "rt112", percentile = 10, outdir = de_dir)
run_pca_plot(res_umuc1$tpm, res_umuc1$samples, "umuc1", percentile = 10, outdir = de_dir)

######################
### DESeq2 Analysis ##
######################
# settings
padj_threshold <- 0.05 # (EDIT): adjusted p-value threshold for significance, particularly important for ORA
log2fc_threshold <- 0.5 # (EDIT): log2 fold change threshold for significance, particularly important for ORA
gsea_padj_threshold <- 0.05 # (EDIT): adjusted p-value threshold for GSEA significance
nTop <- 50 # number of top pathways to display in GSEA barplot

unique(samples_info$treatment) # checks for unique treatment values and puts them into a list

# helper function: map gene symbols to entrez ids
map_symbols_to_entrez <- function(symbols) {
  bitr(
    unique(na.omit(symbols)),
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  ) %>%
    distinct(SYMBOL, ENTREZID) %>%
    rename(GeneSymbol = SYMBOL)
}

# helper function: tidy ora results
tidy_ora <- function(x, method, database, ontology, direction) {
  if (is.null(x) || nrow(as.data.frame(x)) == 0) return(tibble())
  
  as.data.frame(x) %>%
    as_tibble() %>%
    mutate(
      method = method,
      database = database,
      ontology = ontology,
      direction = direction,
      raw_score = -log10(p.adjust),
      rank_score = raw_score
    )
}

# helper function: tidy gsea results
tidy_gsea <- function(x, method, database, ontology, direction) {
  if (is.null(x) || nrow(as.data.frame(x)) == 0) return(tibble())
  
  as.data.frame(x) %>%
    as_tibble() %>%
    mutate(
      method = method,
      database = database,
      ontology = ontology,
      direction = direction,
      raw_score = NES,
      rank_score = abs(NES)
    )
}

comparisons <- c("vp_vs_vr", "rp_vs_rr", "up_vs_ur") # (EDIT)

all_deseq2_results_list <- list() # empty lists to be filled with DESeq2 results
all_gsea_results_list <- list() # empty lists to be filled with GSEA results
all_vst_results_list <- list() # empty lists to be filled with variance-stabilizing transformed normalized counts
all_pathways_tidy_list <- list() # empty lists to be filled with tidy pathway analysis results

for (comparison in comparisons) {
    
    # establish treatment and control conditions (splitting string and using each element as a condition)
    exp <- str_split_fixed(comparison, "_vs_", 2)[, 2]
    ctrl <- str_split_fixed(comparison, "_vs_", 2)[, 1]

    # validate extracted conditions
    print(comparison)
    print(exp)
    print(ctrl)

    # create output directory
    dir_out <- paste0(de_dir, comparison, "/")
    dir.create(dir_out)

    # output directories
    fn_stats <- paste0(dir_out, sprintf("%s_stats.txt", comparison))
    fn_summary_txt <- paste0(dir_out, sprintf("%s_summary.txt", comparison))
    fn_summary_xlsx <- paste0(dir_out, sprintf("%s_summary.xlsx", comparison))
    fn_de_up <- paste0(dir_out, sprintf("%s_up_genes.txt", comparison))
    fn_de_down <- paste0(dir_out, sprintf("%s_down_genes.txt", comparison))
    fn_summary_txt_lfc <- paste0(dir_out, sprintf("%s_lfcs_summary.txt", comparison))
    fn_summary_xlsx_lfc <- paste0(dir_out, sprintf("%s_lfcs_summary.xlsx", comparison))
    fn_de_up_lfc <- paste0(dir_out, sprintf("%s_lfcs_up_genes.txt", comparison))
    fn_de_down_lfc <- paste0(dir_out, sprintf("%s_lfcs_down_genes.txt", comparison))

    # specify the conditions to perform DE analysis on
    sampleCondition <- samples_info[,c("sample_id", "treatment")] # using sample information to generate a new metadata table
    sampleCondition$condition <- substr(sampleCondition$sample_id, 1, 2) # creating a new column for the treatment
    sampleCondition <- sampleCondition[which(sampleCondition$condition %in% c(ctrl, exp)),] # filtering for the treatments and controls of interest
    rownames(sampleCondition) <- sampleCondition$sample_id # setting the row names to the sample ids

    cts <- count_matrix[, rownames(sampleCondition)] # extracting the count matrix for the samples of interest using the metadata table

    # ensuring that the rownames of the sample table and the colnames of the count matrix match
    stopifnot(all(rownames(sampleCondition) == colnames(cts)))

    # construct dds for DE analysis from count matrix (required for DESeq2); design can be modified for batch effects
    dds <- DESeqDataSetFromMatrix(countData = cts, colData = sampleCondition, design = ~ condition) # (EDIT)

    keep <- rowSums(counts(dds)) >= 10 # establishing a filter that removes genes with less than 10 reads total
    dds <- dds[keep,] # performing pre-filtering on the dds object
    dds$condition <- factor(dds$condition, levels = c(ctrl, exp)) # orders the dds object by condition


    dds <- DESeq(dds) # performing DESeq2 analysis
    res <- results(dds, alpha = padj_threshold) # extracting the results from the DESeq2 analysis
    summary(res) # printing the summary of the results
    
    # saving the summary to a text file
    sink(fn_stats)
    print(summary(res))
    sink()

    resultsNames(dds) # establishes the model coefficients (beta values in the generalized linear model) and design (batch, condition, etc.) that DESeq2 used 
    
    # for design(dds) <- ~ condition,           log2(mean expression) = beta0 + beta1 * condition
    # for design(dds) <- ~ condition + batch,   log2(mean expression) = beta0 + beta1 * condition + beta2 * batch
  
    # log fold change shrinkage (makes high variant, low expressed genes less prominent in datasets; this ultimately improves estimates of true effect size and visualization)
    # note that shrinkage does not require alteration of padj since the model is not altered; it is simply a post-estimation stabilization step
    resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm") # performs log fold change shrinkage

    # save DE results (logFC)
    summary1 <- as.data.frame(res)
    summary1 <- summary1[with(summary1, order(padj)),]
    summary1$GeneSymbol <- gene_info$name[match(rownames(summary1), gene_info$ensembl_ID)]
    summary1$gene_type <- gene_info$type[match(rownames(summary1), gene_info$ensembl_ID)]

    write.table(summary1, fn_summary_txt, sep = "\t", row.names = F, col.names = T, quote = F)
    write.xlsx(summary1, fn_summary_xlsx, overwrite = T)

    genes_up <- summary1$GeneSymbol[which(summary1$padj < padj_threshold & summary1$log2FoldChange > log2fc_threshold)]
    genes_down <- summary1$GeneSymbol[which(summary1$padj < padj_threshold & summary1$log2FoldChange < -log2fc_threshold)]

    write.table(genes_up, fn_de_up, row.names = TRUE, col.names = TRUE, quote = FALSE)
    write.table(genes_down, fn_de_down, row.names = TRUE, col.names = TRUE, quote = FALSE)

    # save DE results (logFC - shrunk)
    summary2 <- as.data.frame(resLFC)
    summary2 <- summary2[with(summary2, order(padj)),]
    summary2$GeneSymbol <- gene_info$name[match(rownames(summary2), gene_info$ensembl_ID)]
    summary2$gene_type <- gene_info$type[match(rownames(summary2), gene_info$ensembl_ID)]

    write.table(summary2, fn_summary_txt_lfc, sep = "\t", row.names = TRUE, col.names = T, quote = F)
    write.xlsx(summary2, fn_summary_xlsx_lfc, overwrite = T)

    genes_up <- summary2$GeneSymbol[which(summary2$padj < padj_threshold & summary2$log2FoldChange > log2fc_threshold)]
    genes_down <- summary2$GeneSymbol[which(summary2$padj < padj_threshold & summary2$log2FoldChange < -log2fc_threshold)]

    write.table(genes_up, fn_de_up_lfc, row.names = TRUE, col.names = F, quote = F)
    write.table(genes_down, fn_de_down_lfc, row.names = TRUE, col.names = F, quote = F)

    # generate combined summary table with both logFC and shrunk logFC
    # this table is used for all downstream analyses
    # using raw versus shrunken LFC is context-dependent; in this set of code, it is only used for the volcano plot
    colnames(resLFC)[colnames(resLFC) == "log2FoldChange"] <- "log2FoldChange_shrunk"
    colnames(resLFC)[colnames(resLFC) == "lfcSE"] <- "lfcSE_shrunk"
    
    summary <- summary1
    summary$log2FoldChange_shrunk <- resLFC[rownames(summary), "log2FoldChange_shrunk"]
    summary$lfcSE_shrunk <- resLFC[rownames(summary), "lfcSE_shrunk"]

    # combine DESeq2 results for the downstream analysis
    deseq2_results_df <- as.data.frame(summary)
    deseq2_results_df$comparison <- comparison
    all_deseq2_results_list[[comparison]] <- deseq2_results_df

    # extract the variance-stabilized transformation of counts (i.e., normalizing counts) for downstream use (e.g., gsea heatmaps, etc.)
    # NOT Z SCORED
      vsd <- vst(dds, blind = FALSE)

      mat_vst <- assay(vsd)

      all_vst_results_list[[comparison]] <- list(
        vst = mat_vst
      )

    ######################
    #### Volcano Plot ####
    ######################
    fn_volcano <- paste0(dir_out, sprintf("%s Volcano Plot (Raw LFC).pdf", comparison))
    fn_volcano_shrunk <- paste0(dir_out, sprintf("%s Volcano Plot (Shrunk LFC).pdf", comparison))

    # filter by padj and LFC that meet thresholds
    summary_filtered <- summary[which(!is.na(summary$padj)),]
    summary_filtered <- mutate(summary_filtered, 
                              sig = case_when(
                                summary_filtered$padj < padj_threshold & summary_filtered$log2FoldChange > log2fc_threshold ~ "Up",
                                summary_filtered$padj < padj_threshold & summary_filtered$log2FoldChange < -log2fc_threshold ~ "Down",
                                TRUE ~ "NS"
                              ))

    summary_filtered_ranked <- summary_filtered[order(summary_filtered$stat),] # order by Wald statistic

    # identify top genes for labeling 
    sig_volc <- summary_filtered %>% filter(sig != "NS")
    label_volc <- bind_rows(
      sig_volc %>% arrange(desc(stat)) %>% head(10),
      sig_volc %>% arrange(stat) %>% head(10),
      summary_filtered %>% filter(GeneSymbol == "FOLH1")
    ) %>%
      distinct(GeneSymbol, .keep_all = TRUE)

    # volcano plot (raw LFC)
    # visualization
    volc_raw <- ggplot(summary_filtered, aes(log2FoldChange, -log10(padj))) +
      geom_point(aes(col = sig), alpha = 0.7) +
      geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed", color = "black", alpha = 0.25) +
      geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "black", alpha = 0.25) +
      scale_color_manual(values = c("NS" = "gray", "Up" = "#D81B60", "Down" = "#1E88E5")) +
      geom_label_repel(
        data = label_volc,
        aes(label = GeneSymbol),
        size = 2,
        label.size = 0.25,      # border thickness
        label.padding = 0.15,   # space inside box
        box.padding = 0.5,      # space around label
        point.padding = 0.3,    # space around point
        fill = "white",         # box fill
        color = "black",
        segment.color = "black", # line from label to point
        segment.size = 0.3,      # line thickness
        min.segment.length = 0,  # always draw line
        max.overlaps = Inf
      ) +
      theme_bw() +
      theme(text = element_text(size = 12), 
            legend.title = element_blank(),
            legend = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      xlab("log2FoldChange (shrunk)") +
      ylab("-log10(padj)")

    pdf(fn_volcano, width = 8, height = 5)
    print(volc_raw)
    dev.off()

    # volcano plot (shrunk LFC)
    # visualization
    volc_shrunk <- ggplot(summary_filtered, aes(log2FoldChange_shrunk, -log10(padj))) +
      geom_point(aes(col = sig), alpha = 0.7) +
      geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed", color = "black", alpha = 0.25) +
      geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "black", alpha = 0.25) +
      scale_color_manual(values = c("NS" = "gray", "Up" = "#D81B60", "Down" = "#1E88E5")) +
      geom_label_repel(
        data = label_volc,
        aes(label = GeneSymbol),
        size = 2,
        label.size = 0.25,      # border thickness
        label.padding = 0.15,   # space inside box
        box.padding = 0.5,      # space around label
        point.padding = 0.3,    # space around point
        fill = "white",         # box fill
        color = "black",
        segment.color = "black", # line from label to point
        segment.size = 0.3,      # line thickness
        min.segment.length = 0,  # always draw line
        max.overlaps = Inf
      ) +
      theme_bw() +
      theme(text = element_text(size = 12), 
            legend.title = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      xlab("log2FoldChange (shrunk)") +
      ylab("-log10(padj)")

    pdf(fn_volcano_shrunk, width = 8, height = 5)
    print(volc_shrunk)
    dev.off()

    ######################
    ###### MA Plot #######
    ######################
    fn_ma_raw <- paste0(dir_out, sprintf("%s MA Plot (Raw LFC).pdf", comparison))
    fn_ma_shrunk <- paste0(dir_out, sprintf("%s MA Plot (Shrunk LFC).pdf", comparison))

    # filters out genes with no baseMean
    ma_df <- summary_filtered
    ma_df <- ma_df[!is.na(ma_df$baseMean) & ma_df$baseMean > 0, ]

    # identify top genes for labeling
    label_ma <- ma_df %>%
      filter(!is.na(padj),
            padj < padj_threshold,
            abs(log2FoldChange_shrunk) > log2fc_threshold,
            baseMean > 50,
            sig %in% c("Up","Down")) %>%
      mutate(score = -log10(padj) * abs(log2FoldChange_shrunk)) %>%
      group_by(sig) %>%
      arrange(desc(score)) %>%
      slice_head(n = 30) %>%
      ungroup() %>%
      distinct(GeneSymbol, .keep_all = TRUE)

    # MA plot (raw LFC)
    p_ma_raw <- ggplot(ma_df, aes(x = log10(baseMean), y = log2FoldChange)) +
      geom_point(aes(color = sig), alpha = 0.6, size = 1) +
      geom_hline(yintercept = c(-log2fc_threshold, log2fc_threshold),
                linetype = "dashed", color = "black", alpha = 0.25) +
      scale_color_manual(values = c("NS" = "gray", "Up" = "#D81B60", "Down" = "#1E88E5")) +
      geom_text_repel(
        data = label_ma,
        aes(label = GeneSymbol),
        size = 2,
        max.overlaps = Inf,
        box.padding = 0.4,
        point.padding = 0.2,
        segment.size = 0.3,
        segment.color = "black"
      ) +
      theme_bw() +
      theme(
        text = element_text(size = 12),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      xlab("log10(baseMean)") +
      ylab("log2FoldChange")

    pdf(fn_ma_raw, width = 8, height = 5)
    print(p_ma_raw)
    dev.off()

    # MA plot (shrunk LFC)
    p_ma_shrunk <- ggplot(ma_df, aes(x = log10(baseMean), y = log2FoldChange_shrunk)) +
      geom_point(aes(color = sig), alpha = 0.6, size = 1) +
      geom_hline(yintercept = c(-log2fc_threshold, log2fc_threshold),
                linetype = "dashed", color = "black", alpha = 0.25) +
      scale_color_manual(values = c("NS" = "gray", "Up" = "#D81B60", "Down" = "#1E88E5")) +
      geom_text_repel(
        data = label_ma,
        aes(label = GeneSymbol),
        size = 2,
        max.overlaps = Inf,
        box.padding = 0.4,
        point.padding = 0.2,
        segment.size = 0.3,
        segment.color = "black"
      ) +
      theme_bw() +
      theme(
        text = element_text(size = 12),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      xlab("log10(baseMean)") +
      ylab("log2FoldChange (shrunk)")

    pdf(fn_ma_shrunk, width = 8, height = 5)
    print(p_ma_shrunk)
    dev.off()

    ######################
    ######## GSEA ########
    ######################
    ## QUESTION: Gene set enrichment analysis
    ## The changes that I see in my gene sets are most similar to what known list of changes?
    # prepare gene lists for gsea
    gsea_map <- map_symbols_to_entrez(summary_filtered$GeneSymbol)

    # map gene symbols to entrez ids and rank by statistic, maintaining only the highest stat value for entrez ids with multiple gene symbols
    gsea_rank_df <- summary_filtered %>%
      select(GeneSymbol, stat) %>%
      inner_join(gsea_map, by = "GeneSymbol") %>%
      group_by(ENTREZID) %>%
      summarise(stat = stat[which.max(abs(stat))]) %>%
      arrange(desc(stat))

    geneList_entrez <- setNames(gsea_rank_df$stat, gsea_rank_df$ENTREZID)
    geneList_entrez <- sort(geneList_entrez, decreasing = TRUE)

    # gsea - fgsea - hallmark
    genes_rank <- summary_filtered$stat # establish genes_rank
    names(genes_rank) <- summary_filtered$GeneSymbol # establish names for genes_rank

    genes_rank_averaged <- summary_filtered %>% # averaging repeated GeneSymbols
      group_by(GeneSymbol) %>%
      summarise(stat = mean(stat, na.rm = TRUE)) %>%
      ungroup()

    genes_rank <- setNames(genes_rank_averaged$stat, genes_rank_averaged$GeneSymbol) # convert to a named numeric vector for fgsea

    set.seed(10) # since fgsea estimates enrichment p-values via random permutation, setting a seed ensures reproducibility
    fgseaRes <- fgsea(pathways = pathways.hallmark, stats = genes_rank, nperm = 100000, maxSize = 500, nproc = 6)# run fgsea with the filtered genes_rank

    fgsea_df <- as.data.frame(fgseaRes) # combine gsea results for the hallmark dot plot
    fgsea_df$comparison <- comparison
    all_gsea_results_list[[comparison]] <- fgsea_df

    theDF <- fgseaRes

    theDF$theLabel <- NA
    theDF$theLabel[theDF$padj < gsea_padj_threshold & !is.na(theDF$padj)] <- theDF$pathway[theDF$padj < gsea_padj_threshold & !is.na(theDF$padj)]

    
    theDF_top <- theDF[order(abs(theDF$NES), decreasing = T)[1:nTop],] # visualization
    theDF_top$pathway_clean <- str_remove(theDF_top$pathway, "^HALLMARK_")

    p_top <- ggplot(theDF_top, aes(reorder(pathway_clean,-NES),NES)) +
      geom_col(aes(fill = theDF_top$padj < gsea_padj_threshold)) +
      coord_flip() +
      labs(x="", y="Normalized Enrichment Score",
           title="",
           fill = sprintf("padj < %s", gsea_padj_threshold)) +
      theme_classic() +
      theme(legend.position = "right",
            axis.text = element_text(size = 7),
            axis.title = element_text(size = 9),
            axis.line.x = element_line(),
            axis.line.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "transparent",colour = NA),
            plot.background = element_rect(fill = "transparent",colour = NA)) +
      scale_fill_manual(values = list("TRUE" = "#D81B60", "FALSE" = "gray"))

    print(p_top)

    fn_gsea_top <- paste0(dir_out, sprintf("%s GSEA Plot.pdf", comparison))
    pdf(fn_gsea_top, width = 8, height = 5)
    print(p_top)
    dev.off()
    
    # gsea - go - bp
    gsea_go_bp <- gseGO(
      geneList = geneList_entrez,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      verbose = FALSE
      )
    # gsea - go - mf
    gsea_go_mf <- gseGO(
      geneList = geneList_entrez,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "MF",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      verbose = FALSE
      )
    # gsea - go - cc
    gsea_go_cc <- gseGO(
      geneList = geneList_entrez,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "CC",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      verbose = FALSE
      )
  # gsea - kegg
    gsea_kegg <- gseKEGG(
      geneList = geneList_entrez,
      organism = "hsa",
      keyType = "ncbi-geneid",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      verbose = FALSE
      )
    # gsea - reactome
    gsea_reactome <- ReactomePA::gsePathway(
      geneList = geneList_entrez,
      organism = "human",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05
      )
    
    #######################
    ######### ORA #########
    #######################
    ## QUESTION: Over-representation analysis
    ## WHICH PATHWAYS ARE OVERREPRESENTED IN MY LIST OF SIGNIFICANT GENES (p.adjusted < 0.05, log2FoldChange > 0.5 or < -0.5)?

    # prepare gene lists for ora
    deg_tbl <- summary2 %>%
      filter(!is.na(padj)) %>% # remove genes with no p-value
      mutate(
        GeneSymbol = ifelse(is.na(GeneSymbol), NA, GeneSymbol) # replace empty gene symbols with NA
      )
    # filter for significant genes for ora
    genes_up_sym <- deg_tbl$GeneSymbol[ 
      deg_tbl$padj < padj_threshold & deg_tbl$log2FoldChange > log2fc_threshold
    ]

    genes_dn_sym <- deg_tbl$GeneSymbol[
      deg_tbl$padj < padj_threshold & deg_tbl$log2FoldChange < -log2fc_threshold
    ]

    # gene universe for ora: all tested genes that map to ENTREZ
    universe_map <- map_symbols_to_entrez(deg_tbl$GeneSymbol)
    universe_entrez <- unique(universe_map$ENTREZID)

    length(universe_entrez)
    
    # gene up/down for ora: up/downregulated genes that map to ENTREZ
    up_entrez <- map_symbols_to_entrez(genes_up_sym)$ENTREZID %>% unique()
    dn_entrez <- map_symbols_to_entrez(genes_dn_sym)$ENTREZID %>% unique()

    # ora - go - bp - up
    ora_go_bp_up <- enrichGO(
      gene = up_entrez,
      universe = universe_entrez,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = FALSE
      )
    # ora - go - mf - up
    ora_go_mf_up <- enrichGO(
      gene = up_entrez,
      universe = universe_entrez,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "MF",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = FALSE
      )
    # ora - go - cc - up
    ora_go_cc_up <- enrichGO(
      gene = up_entrez,
      universe = universe_entrez,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "CC",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = FALSE
      )
    # ora - kegg - up
    ora_kegg_up <- enrichKEGG(
      gene = up_entrez,
      organism = "hsa",
      keyType = "ncbi-geneid",
      universe = universe_entrez,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05
      )
    # ora - reactome - up
    ora_reactome_up <- ReactomePA::enrichPathway(
      gene = up_entrez,
      universe = universe_entrez,
      organism = "human",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = FALSE
      )
    # ora - go - bp - dn
    ora_go_bp_dn <- enrichGO(
      gene = dn_entrez,
      universe = universe_entrez,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = FALSE
      )
    # ora - go - mf - dn
    ora_go_mf_dn <- enrichGO(
      gene = dn_entrez,
      universe = universe_entrez,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "MF",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = FALSE
      )
    # ora - go - cc - dn
    ora_go_cc_dn <- enrichGO(
      gene = dn_entrez,
      universe = universe_entrez,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "CC",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = FALSE
      )
    # ora - kegg - dn
    ora_kegg_dn <- enrichKEGG(
      gene = dn_entrez,
      organism = "hsa",
      keyType = "ncbi-geneid",
      universe = universe_entrez,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05
      )
    # ora - reactome - dn
    ora_reactome_dn <- ReactomePA::enrichPathway(
      gene = dn_entrez,
      universe = universe_entrez,
      organism = "human",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = FALSE
      )

    # combine results into tidy dataframe
    pathway_table <- bind_rows(
      tidy_ora(ora_go_bp_up, "ORA", "GO", "BP", "Up"),
      tidy_ora(ora_go_mf_up, "ORA", "GO", "MF", "Up"),
      tidy_ora(ora_go_cc_up, "ORA", "GO", "CC", "Up"),
      tidy_ora(ora_kegg_up, "ORA", "KEGG", NA, "Up"),
      tidy_ora(ora_reactome_up, "ORA", "Reactome", NA, "Up"),

      tidy_ora(ora_go_bp_dn, "ORA", "GO", "BP", "Down"),
      tidy_ora(ora_go_mf_dn, "ORA", "GO", "MF", "Down"),
      tidy_ora(ora_go_cc_dn, "ORA", "GO", "CC", "Down"),
      tidy_ora(ora_kegg_dn, "ORA", "KEGG", NA, "Down"),
      tidy_ora(ora_reactome_dn, "ORA", "Reactome", NA, "Down"),

      tidy_gsea(gsea_go_bp, "GSEA", "GO", "BP", "Ranked"),
      tidy_gsea(gsea_go_mf, "GSEA", "GO", "MF", "Ranked"),
      tidy_gsea(gsea_go_cc, "GSEA", "GO", "CC", "Ranked"),
      tidy_gsea(gsea_kegg, "GSEA", "KEGG", NA, "Ranked"),
      tidy_gsea(gsea_reactome, "GSEA", "Reactome", NA, "Ranked")) %>%
      group_by(method, database, ontology, direction) %>%
      mutate(rank_within_group = dense_rank(desc(rank_score))) %>%
      ungroup() %>%
      arrange(method, database, ontology, rank_within_group)

    pathway_table$comparison <- comparison
    all_pathways_tidy_list[[comparison]] <- pathway_table
  }


######################
## Master Dataframe ##
######################
# combine and save deseq2 results
master_df_deseq2 <- bind_rows(all_deseq2_results_list)
write.xlsx(master_df_deseq2, paste0(de_dir, "master_df_deseq2.xlsx"), overwrite = T)

# combine and save gsea results
master_df_gsea <- bind_rows(all_gsea_results_list)
write.xlsx(master_df_gsea, paste0(de_dir, "master_df_gsea.xlsx"), overwrite = T)

# combine and save pathway results
master_df_pathway <- bind_rows(all_pathways_tidy_list)
write.xlsx(master_df_pathway, paste0(de_dir, "master_df_pathway.xlsx"), overwrite = T)

# combine and save vst counts
wb <- createWorkbook()
for (comparison in names(all_vst_results_list)) {
  
  mat <- all_vst_results_list[[comparison]]$vst
  
  df <- as.data.frame(mat)
  df$GeneSymbol <- gene_info$name[match(rownames(df), gene_info$ensembl_ID)]
  df <- df[, c("GeneSymbol", setdiff(colnames(df), "GeneSymbol"))]
  
  addWorksheet(wb, comparison)
  writeData(wb, comparison, df)
}
saveWorkbook(wb, paste0(de_dir, "master_vst.xlsx"), overwrite = TRUE)

######################
#### Common Genes ####
######################
# intersection approach using individual cell line deseq2 results
# pros:
# strict/conservative: ensures that genes are significant in each cell line
# ease of interpretation: specific genes are differentially expressed across all cell lines
# does not assume effect size is shared
# cons: 
# low statistical power: data is split across individual deseq2 analyses
# false negative odds increase: genes that are missing significance may be dropped
# sensitivity to noise: a noisy dataset can drown out real signals
# ignores shared information: data is weaker when split into individual cell lines

master_df_deseq2_1 <- read.xlsx(paste0(de_dir, "master_df_deseq2.xlsx")) # (EDIT): change this to the name of your desired master dataframe (individual comparisons versus intersection comparisons)
View(master_df_deseq2_1)
master_df_deseq2_2 <- read.table(paste0(de_dir, "EV_vs_PBS_DESeq2_results.txt"), header = TRUE, sep = "\t") # (EDIT): change this to the name of your desired master dataframe (individual comparisons versus intersection comparisons)
View(master_df_deseq2_2)

cols <- c("GeneSymbol", "log2FoldChange", "padj", "comparison")

master_df_deseq2_1 <- master_df_deseq2_1 %>%
  select(all_of(cols))
View(master_df_deseq2_1)
master_df_deseq2_2 <- master_df_deseq2_2 %>%
  rename(GeneSymbol = `gene`) %>%
  mutate(comparison = "evrs")

master_df_deseq2_2 <- master_df_deseq2_2 %>%
  select(all_of(cols))
View(master_df_deseq2_2)

master_df_deseq2 <- bind_rows(master_df_deseq2_1, master_df_deseq2_2)
View(master_df_deseq2)
# visualization settings
comparison_labels <- c( # (EDIT)
  "vp_vs_vr" = "EVR 647V - Chou",
  "rp_vs_rr" = "EVR RT112 - Chou",
  "up_vs_ur" = "EVR UMUC1 - Chou",
  "evrs" = "EVR RT112 - Chu"
)

comparison_colors <- c( # (EDIT)
  "vp_vs_vr" = "#D81B60",
  "rp_vs_rr" = "#1E88E5",
  "up_vs_ur" = "#FFC107",
  "evrs" = "#004D40"
)

intersection_comparisons <- c("vp_vs_vr", "rp_vs_rr", "up_vs_ur", "evrs") # (EDIT)

# filter data by padj < 0.05, abs(log2FC) > 0.5, and protein coding genes only
filtered_df <- master_df_deseq2 %>%
  filter(
    comparison %in% intersection_comparisons,
    padj < 0.05,
    abs(log2FoldChange) > 0.5,
    !is.na(GeneSymbol)
  )

# generate list of up/down genes
up_list <- filtered_df %>%
  filter(log2FoldChange > 0) %>%
  group_by(comparison) %>%
  summarise(genes = list(unique(GeneSymbol))) %>%
  deframe()

down_list <- filtered_df %>%
  filter(log2FoldChange < 0) %>%
  group_by(comparison) %>%
  summarise(genes = list(unique(GeneSymbol))) %>%
  deframe()

# convert list of up/down genes for upset plots (i.e., venn diagrams)
list_to_binary_df <- function(gene_list) {
  all_genes <- unique(unlist(gene_list))
  df <- data.frame(gene = all_genes)
  for (comp in names(gene_list)) {
    df[[comp]] <- all_genes %in% gene_list[[comp]]
  }
  df
}

up_df <- list_to_binary_df(up_list)
dn_df <- list_to_binary_df(down_list)

# upset plots
p_up <- upset(
  up_df,
  intersect = intersection_comparisons,
  sort_intersections_by = "degree"
) +
  plot_annotation(
    title = "Common Upregulated Genes (padj < 0.05, log2FC > 0.5)"
  ) &
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  )

p_dn <- upset(
  dn_df,
  intersect = intersection_comparisons,
  sort_intersections_by = "degree"
) +
  plot_annotation(
    title = "Common Downregulated Genes (padj < 0.05, log2FC < -0.5)"
  ) &
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  )

# identify common genes between the cell lines
common_up <- Reduce(intersect, up_list)
common_dn <- Reduce(intersect, down_list)

common_up_long <- filtered_df %>%
  filter(GeneSymbol %in% common_up, log2FoldChange > 0)

common_dn_long <- filtered_df %>%
  filter(GeneSymbol %in% common_dn, log2FoldChange < 0)

# rank genes using a score (score = mean - spread of data)
# up common genes
top_up_genes <- common_up_long %>%
  group_by(GeneSymbol) %>%
  summarise(
    mean_log2FC = mean(log2FoldChange),
    spread = max(log2FoldChange) - min(log2FoldChange),
    score = mean_log2FC - spread,
    .groups = "drop"
  ) %>%
  arrange(desc(score)) %>%
  slice_head(n = 50) %>%
  pull(GeneSymbol)

# down common genes
top_dn_genes <- common_dn_long %>%
  group_by(GeneSymbol) %>%
  summarise(
    mean_log2FC = mean(log2FoldChange),
    spread = max(log2FoldChange) - min(log2FoldChange),
    score = -mean_log2FC - spread,
    .groups = "drop"
  ) %>%
  arrange(desc(score)) %>%
  slice_head(n = 50) %>%
  pull(GeneSymbol)

# filtering for plotting
plot_up_df <- common_up_long %>%
  filter(GeneSymbol %in% top_up_genes)

plot_dn_df <- common_dn_long %>%
  filter(GeneSymbol %in% top_dn_genes)

# dot plots
p_up_dot <- ggplot(plot_up_df,
                   aes(x = log2FoldChange,
                       y = reorder(GeneSymbol, log2FoldChange),
                       color = comparison)) +
  theme_classic() +
  theme(
    panel.grid.major.y = element_line(color = "grey90", size = 0.2),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 8)
  ) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 8)) +
  geom_point(size = 2, alpha = 0.9) +
  scale_color_manual(values = comparison_colors,
                     labels = comparison_labels) +
  labs(title = sprintf("Top %s Common Upregulated Genes", nTop),
       x = "log2FoldChange", y = "", color = "Cell Line")

p_dn_dot <- ggplot(plot_dn_df,
                   aes(x = log2FoldChange,
                       y = reorder(GeneSymbol, log2FoldChange),
                       color = comparison)) +
  theme_classic() +
  theme(
    panel.grid.major.y = element_line(color = "grey90", size = 0.2),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 8)
  ) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 8)) +
  geom_point(size = 2, alpha = 0.9) +
  scale_color_manual(values = comparison_colors,
                     labels = comparison_labels) +
  labs(title = sprintf("Top %s Common Downregulated Genes", nTop),
       x = "log2FoldChange", y = "", color = "Cell Line")

# gene ontology
up_entrez <- bitr(common_up, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
dn_entrez <- bitr(common_dn, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

ego_up <- enrichGO(
  gene = up_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

ego_dn <- enrichGO(
  gene = dn_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

# helper function for prettier GO plots
pretty_go_plot <- function(ego, title) {
  
  # shorten / wrap GO terms
  ego@result$Description <- str_wrap(ego@result$Description, width = 40)
  
  p <- dotplot(ego, showCategory = 12) +   # fewer categories = cleaner
    scale_color_gradient(low = "#1E88E5", high = "#D81B60") +
    ggtitle(title) +
    
    # cleaner theme
    theme_classic(base_size = 8) +
    theme(
      axis.text.y = element_text(size = 7),   # smaller GO labels
      axis.text.x = element_text(size = 7),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 7)
    ) +
    
    # nicer p-value labels
    scale_color_continuous(labels = scales::scientific_format(digits = 2))
  
  return(p)
}

# save to pdf
pdf(paste0(de_dir, "Intersection Analysis (Cell Line DESeq2 Matrix Comparisons).pdf"), width = 8, height = 5)

# upset plots
print(p_up)
print(p_dn)

# dot plots
print(p_up_dot)
print(p_dn_dot)

# gene ontology plots
print(pretty_go_plot(ego_up, "GO Enrichment (Upregulated)"))

print(pretty_go_plot(ego_dn, "GO Enrichment (Downregulated)"))

dev.off()




evrs <- read.table(paste0(de_dir, "EV_vs_PBS_DESeq2_results.txt"), header = TRUE, sep = "\t") # (EDIT): change this to the name of your desired master dataframe (individual comparisons versus intersection comparisons)

# gsea on ev resistance/sensitivity screen
evrs_filtered <- evrs[which(!is.na(evrs$padj)),]
evrs_filtered <- mutate(evrs_filtered, 
                              sig = case_when(
                                evrs_filtered$padj < padj_threshold & evrs_filtered$log2FoldChange > log2fc_threshold ~ "Up",
                                evrs_filtered$padj < padj_threshold & evrs_filtered$log2FoldChange < -log2fc_threshold ~ "Down",
                                TRUE ~ "NS"
                              ))

evrs_mapped <- map_symbols_to_entrez(evrs_filtered$gene)

# map gene symbols to entrez ids and rank by statistic, maintaining only the highest stat value for entrez ids with multiple gene symbols
gsea_rank_evrs <- evrs_filtered %>%
  select(gene, stat) %>%
  inner_join(evrs_mapped, by = c("gene" = "GeneSymbol")) %>%
  group_by(ENTREZID) %>%
  summarise(stat = stat[which.max(abs(stat))]) %>%
  arrange(desc(stat))

geneList_entrez_evrs <- setNames(gsea_rank_evrs$stat, gsea_rank_evrs$ENTREZID)
geneList_entrez_evrs <- sort(geneList_entrez_evrs, decreasing = TRUE)

# gsea - fgsea - hallmark
genes_rank <- evrs_filtered$stat # establish genes_rank
names(genes_rank) <- evrs_filtered$gene # establish names for genes_rank

genes_rank_averaged <- evrs_filtered %>% # averaging repeated GeneSymbols
  group_by(gene) %>%
  summarise(stat = mean(stat, na.rm = TRUE)) %>%
  ungroup()

genes_rank <- setNames(genes_rank_averaged$stat, genes_rank_averaged$gene) # convert to a named numeric vector for fgsea

set.seed(10) # since fgsea estimates enrichment p-values via random permutation, setting a seed ensures reproducibility
fgseaRes <- fgsea(pathways = pathways.hallmark, stats = genes_rank, nperm = 100000, maxSize = 500, nproc = 6)# run fgsea with the filtered genes_rank

fgsea_df <- as.data.frame(fgseaRes) # combine gsea results for the hallmark dot plot

theDF <- fgseaRes

theDF$theLabel <- NA
theDF$theLabel[theDF$padj < gsea_padj_threshold & !is.na(theDF$padj)] <- theDF$pathway[theDF$padj < gsea_padj_threshold & !is.na(theDF$padj)]

write.xlsx(theDF, paste0(de_dir, "EV Resistance Screen GSEA Results.xlsx"), overwrite = T)

theDF_top <- theDF[order(abs(theDF$NES), decreasing = T)[1:nTop],] # visualization
theDF_top$pathway_clean <- str_remove(theDF_top$pathway, "^HALLMARK_")

p_top <- ggplot(theDF_top, aes(reorder(pathway_clean,-NES),NES)) +
      geom_col(aes(fill = theDF_top$padj < gsea_padj_threshold)) +
      coord_flip() +
      labs(x="", y="Normalized Enrichment Score",
           title="",
           fill = sprintf("padj < %s", gsea_padj_threshold)) +
      theme_classic() +
      theme(legend.position = "right",
            axis.text = element_text(size = 7),
            axis.title = element_text(size = 9),
            axis.line.x = element_line(),
            axis.line.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "transparent",colour = NA),
            plot.background = element_rect(fill = "transparent",colour = NA)) +
      scale_fill_manual(values = list("TRUE" = "#D81B60", "FALSE" = "gray"))

print(p_top)

fn_gsea <- paste0(de_dir, "EV Resistance Screen GSEA Plot.pdf")
pdf(fn_gsea, width = 8, height = 5)
print(p_top)
dev.off()

# filter by padj and LFC that meet thresholds
evrs_filtered <- evrs_filtered[order(evrs_filtered$stat),] # order by Wald statistic

sig_volc <- evrs_filtered %>% filter(sig != "NS")
label_volc <- bind_rows(
      sig_volc %>% arrange(desc(stat)) %>% head(15),
      sig_volc %>% arrange(stat) %>% head(15)
    ) %>%
      distinct(gene, .keep_all = TRUE)

# volcano plot (raw LFC)
# visualization
volc_raw <- ggplot(evrs_filtered, aes(log2FoldChange, -log10(padj))) +
      geom_point(aes(col = sig), alpha = 0.7) +
      geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed", color = "black", alpha = 0.25) +
      geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "black", alpha = 0.25) +
      scale_color_manual(values = c("NS" = "gray", "Up" = "#D81B60", "Down" = "#1E88E5")) +
      geom_label_repel(
        data = label_volc,
        aes(label = gene),
        size = 2,
        label.size = 0.25,      # border thickness
        label.padding = 0.15,   # space inside box
        box.padding = 0.5,      # space around label
        point.padding = 0.3,    # space around point
        fill = "white",         # box fill
        color = "black",
        segment.color = "black", # line from label to point
        segment.size = 0.3,      # line thickness
        min.segment.length = 0,  # always draw line
        max.overlaps = Inf
      ) +
      theme_bw() +
      theme(text = element_text(size = 12), 
            legend.title = element_blank(),
            legend = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      xlab("log2FoldChange") +
      ylab("-log10(padj)")

fn_volcano <- paste0(de_dir, "EV Resistance Screen Volcano Plot.pdf")
pdf(fn_volcano, width = 8, height = 5)
print(volc_raw)
dev.off()

# filters out genes with no baseMean
ma_df <- evrs_filtered
ma_df <- ma_df[!is.na(ma_df$baseMean) & ma_df$baseMean > 0, ]

# identify top genes for labeling
label_ma <- ma_df %>%
      filter(!is.na(padj),
            padj < padj_threshold,
            abs(log2FoldChange) > log2fc_threshold,
            baseMean > 50,
            sig %in% c("Up","Down")) %>%
      mutate(score = -log10(padj) * abs(log2FoldChange)) %>%
      group_by(sig) %>%
      arrange(desc(score)) %>%
      slice_head(n = 30) %>%
      ungroup() %>%
      distinct(gene, .keep_all = TRUE)
   
# MA plot (raw LFC)
p_ma_raw <- ggplot(ma_df, aes(x = log10(baseMean), y = log2FoldChange)) +
      geom_point(aes(color = sig), alpha = 0.6, size = 1) +
      geom_hline(yintercept = c(-log2fc_threshold, log2fc_threshold),
                linetype = "dashed", color = "black", alpha = 0.25) +
      scale_color_manual(values = c("NS" = "gray", "Up" = "#D81B60", "Down" = "#1E88E5")) +
      geom_text_repel(
        data = label_ma,
        aes(label = gene),
        size = 2,
        max.overlaps = Inf,
        box.padding = 0.4,
        point.padding = 0.2,
        segment.size = 0.3,
        segment.color = "black"
      ) +
      theme_bw() +
      theme(
        text = element_text(size = 12),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      xlab("log10(baseMean)") +
      ylab("log2FoldChange")
      
fn_ma_raw <- paste0(de_dir, "EV Resistance Screen MA Plot.pdf")
pdf(fn_ma_raw, width = 8, height = 5)
print(p_ma_raw)
dev.off()