#!/usr/bin/env Rscript
# 08_marker_heatmap.R
# Marker-gene log2FC heatmap (mycotoxin / autophagy / secondary-metabolism markers)
# across all 8 studies. Marker set from tables/supplementary/s2_markergenes.xlsx.
# Non-FGSG studies use the KO-based ortholog of each FGSG marker; markers with no
# ortholog or no test result are shown as NA (grey).

local({
    args <- commandArgs(trailingOnly = FALSE)
    f <- grep("--file=", args, value = TRUE)
    d <- if (length(f)) dirname(normalizePath(sub("--file=", "", f))) else getwd()
    source(file.path(d, "config.R"))
})

suppressPackageStartupMessages({
    library(dplyr); library(readr); library(readxl)
    library(pheatmap); library(svglite); library(RColorBrewer)
})

# Markers (FGSG reference IDs)
markers <- read_excel(file.path(BASE_DIR, "tables", "supplementary", "s2_markergenes.xlsx"))
markers <- markers %>% select(gene_id = Gene_ID, gene_name = Gene_Name, category = Functional_Category)

# FGSG -> native ortholog inverters (one FGSG may have several native orthologs)
fox <- read_csv(file.path(RESOURCE_DIR, "fox_to_fgr_orthology.csv"), show_col_types = FALSE)  # FOXG_ID, FGSG_ID
fpu <- read_csv(file.path(RESOURCE_DIR, "fpu_to_fgr_orthology.csv"), show_col_types = FALSE)  # FPSE_ID, FGSG_ID
fgsg_to_fox <- split(fox$FOXG_ID, fox$FGSG_ID)
fgsg_to_fpu <- split(fpu$FPSE_ID, fpu$FGSG_ID)

org_code_v <- setNames(STUDY_INFO$org_code, as.character(STUDY_INFO$study_id))
bca_v      <- setNames(STUDY_INFO$bca,      as.character(STUDY_INFO$study_id))

# log2FC lookup per study
lfc_of <- function(study_id, fgsg_id) {
    res <- read_csv(file.path(DE_DIR, paste0("Study_", study_id), "deseq2_full_results.csv"),
                    show_col_types = FALSE)
    v <- setNames(res$log2FoldChange, res$gene)
    oc <- org_code_v[[as.character(study_id)]]
    cand <- if (oc == "fgr") fgsg_id
            else if (oc == "fox") fgsg_to_fox[[fgsg_id]]
            else if (oc == "fpu") fgsg_to_fpu[[fgsg_id]]
            else fgsg_id
    cand <- cand[cand %in% names(v)]
    if (length(cand) == 0) return(NA_real_)
    # if several orthologs are tested, take the most strongly regulated
    vals <- v[cand]
    vals[which.max(abs(vals))]
}

mat <- matrix(NA_real_, nrow = nrow(markers), ncol = 8,
              dimnames = list(markers$gene_name, paste0("S", 1:8)))
for (j in 1:8) for (i in seq_len(nrow(markers))) mat[i, j] <- lfc_of(j, markers$gene_id[i])
colnames(mat) <- unname(bca_v[as.character(1:8)])

row_anno <- data.frame(Cluster = markers$category, row.names = markers$gene_name)
col_anno <- data.frame(BCA_Type = STUDY_INFO$bca_type, row.names = colnames(mat))
cat_levels <- unique(markers$category)
cluster_cols <- setNames(brewer.pal(max(3, length(cat_levels)), "Set2")[seq_along(cat_levels)], cat_levels)
anno_colors <- list(BCA_Type = COLORS_BCA, Cluster = cluster_cols)

rng <- max(abs(mat), na.rm = TRUE)
breaks <- seq(-rng, rng, length.out = 101)
pal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

for (ext in c("png", "svg")) {
    fp <- file.path(FIG_DIR, paste0("marker_gene_heatmap.", ext))
    if (ext == "png") png(fp, width = 11, height = 7, units = "in", res = 300) else svglite(fp, width = 11, height = 7)
    pheatmap(mat,
             color = pal, breaks = breaks, na_col = "grey85",
             cluster_rows = FALSE, cluster_cols = FALSE,
             display_numbers = matrix(ifelse(is.na(mat), "NA", sprintf("%.1f", mat)), nrow = nrow(mat)),
             number_color = "black", fontsize_number = 8,
             annotation_row = row_anno, annotation_col = col_anno,
             annotation_colors = anno_colors, gaps_row = cumsum(rle(markers$category)$lengths),
             gaps_col = 4, main = "", fontsize = 11)
    dev.off()
}

write_csv(as.data.frame(mat) %>% mutate(gene = rownames(mat), .before = 1),
          file.path(TABLE_DIR, "marker_gene_log2fc.csv"))
cat("Marker heatmap written. log2FC matrix:\n"); print(round(mat, 2))
