#!/usr/bin/env Rscript
# 00a_deseq2_full_results.R
# Regenerates the full per-gene DESeq2 results (all genes) for each study,
# reproducing the DEG pipeline that produced data/degs (apeglm-shrunken LFC,
# ref = Control, >=10 total counts, padj<0.05 & |log2FC|>1).
# Writes outputs/differential_expression/Study_N/deseq2_full_results.csv
# (consumed by 04_deg_volcano, 05_overlap_analysis, 06_enrichment) and
# verifies the significant subset matches data/degs.

local({
    args <- commandArgs(trailingOnly = FALSE)
    f <- grep("--file=", args, value = TRUE)
    d <- if (length(f)) dirname(normalizePath(sub("--file=", "", f))) else getwd()
    source(file.path(d, "config.R"))
})

suppressPackageStartupMessages({
    library(DESeq2)
    library(apeglm)
    library(readr)
    library(readxl)
    library(dplyr)
})

study_table <- read_csv(file.path(METADATA_DIR, "study_table.csv"), show_col_types = FALSE)

verify <- list()

for (sid in 1:8) {
    counts_file <- file.path(COUNTS_DIR, paste0(sid, ".txt"))
    raw <- read_tsv(counts_file, comment = "#", show_col_types = FALSE)
    srr_cols <- grep("SRR", colnames(raw), value = TRUE)
    mat <- as.matrix(raw[, srr_cols])
    rownames(mat) <- raw$Geneid
    colnames(mat) <- sub(".*?(SRR[0-9]+).*", "\\1", colnames(mat))   # extract SRR id from any path/suffix
    mat[is.na(mat) | mat < 0] <- 0

    st <- study_table %>% filter(`Study No.` == sid)
    control_srr <- st$`SRA Accession`[st$`Control/Sample` == "Control"]
    sample_srr  <- st$`SRA Accession`[st$`Control/Sample` == "Sample"]
    samples <- c(control_srr, sample_srr)
    stopifnot(all(samples %in% colnames(mat)))
    mat <- mat[, samples, drop = FALSE]

    colData <- data.frame(
        row.names = samples,
        condition = factor(c(rep("Control", length(control_srr)),
                             rep("Sample",  length(sample_srr))),
                           levels = c("Control", "Sample"))
    )

    dds <- DESeqDataSetFromMatrix(countData = mat, colData = colData, design = ~ condition)
    dds <- dds[rowSums(counts(dds)) >= 10, ]
    dds <- DESeq(dds, quiet = TRUE)
    resLFC <- lfcShrink(dds, coef = "condition_Sample_vs_Control", type = "apeglm")

    res_df <- as.data.frame(resLFC)
    res_df$gene <- rownames(res_df)
    res_df <- res_df[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]

    out_dir <- file.path(DE_DIR, paste0("Study_", sid))
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    write_csv(res_df, file.path(out_dir, "deseq2_full_results.csv"))

    # Verify significant subset vs data/degs
    up   <- res_df %>% filter(!is.na(padj), padj < PADJ_THRESHOLD, log2FoldChange >  LOG2FC_THRESHOLD)
    down <- res_df %>% filter(!is.na(padj), padj < PADJ_THRESHOLD, log2FoldChange < -LOG2FC_THRESHOLD)
    ref_up   <- read_excel(file.path(DEGS_DIR, paste0(sid, "_Upregulated_Genes.xlsx")))$gene
    ref_down <- read_excel(file.path(DEGS_DIR, paste0(sid, "_Downregulated_Genes.xlsx")))$gene
    verify[[sid]] <- data.frame(
        Study = sid,
        rerun_up = nrow(up), ref_up = length(ref_up),
        rerun_down = nrow(down), ref_down = length(ref_down),
        up_match = length(intersect(up$gene, ref_up)),
        down_match = length(intersect(down$gene, ref_down)),
        total_genes = nrow(res_df)
    )
    cat(sprintf("Study %d: rerun %d up / %d down | ref %d / %d | overlap %d / %d | genes %d\n",
                sid, nrow(up), length(ref_up), nrow(down), length(ref_down),
                length(intersect(up$gene, ref_up)), length(intersect(down$gene, ref_down)), nrow(res_df)))
}

vdf <- bind_rows(verify)
write_csv(vdf, file.path(DE_DIR, "deg_reproduction_check.csv"))
cat("\n=== DEG reproduction vs data/degs ===\n")
print(vdf)
