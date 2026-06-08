#!/usr/bin/env Rscript
# 11_supplementary_tables.R — assemble Supplementary Tables S1-S9 from pipeline outputs.
# Run after 04-09. Writes xlsx files to outputs/supplementary/.
# S1 study metadata | S2 marker genes | S3 DEG lists | S4 DEG KEGG | S5 DEG GO |
# S6 rewired proteins | S7 rewired KEGG | S8 rewired GO | S9 Louvain modules.

local({
    args <- commandArgs(trailingOnly = FALSE)
    f <- grep("--file=", args, value = TRUE)
    d <- if (length(f)) dirname(normalizePath(sub("--file=", "", f))) else getwd()
    source(file.path(d, "config.R"))
})
suppressPackageStartupMessages({ library(readr); library(readxl); library(openxlsx); library(dplyr) })

SUPP_DIR <- file.path(OUTPUT_DIR, "supplementary"); if (!dir.exists(SUPP_DIR)) dir.create(SUPP_DIR, recursive = TRUE)
sheets_from_csvs <- function(pattern, strip) {
    files <- list.files(ENRICH_DIR, pattern = pattern, full.names = TRUE)
    setNames(lapply(files, read_csv, show_col_types = FALSE),
             gsub(strip, "", gsub("\\.csv$", "", basename(files))))
}
save_book <- function(sheets, path) {
    wb <- createWorkbook()
    for (nm in names(sheets)) { s <- substr(nm, 1, 31); addWorksheet(wb, s); writeData(wb, s, sheets[[nm]]) }
    saveWorkbook(wb, path, overwrite = TRUE)
}

# S1 — study metadata / accessions
write.xlsx(read_csv(file.path(METADATA_DIR, "study_table.csv"), show_col_types = FALSE),
           file.path(SUPP_DIR, "S1_study_metadata.xlsx"))

# S2 — marker gene list (canonical F. graminearum IDs)
markers <- data.frame(
  Gene_ID = c("FGSG_03535","FGSG_03537","FGSG_02395","FGSG_02396","FGSG_02398",
              "FGSG_08900","FGSG_10740","FGSG_08491","FGSG_02320","FGSG_02324","FGSG_02329"),
  Gene_Name = c("TRI4","TRI5","ZEA1","ZEA2","ZEB2","ATG3","ATG8","ATG13","aurR1","PKS12","aurS"),
  Functional_Category = c(rep("Trichothecene Biosynthesis",2), rep("Zearalenone Biosynthesis",3),
                          rep("Autophagy",3), rep("Aurofusarin Biosynthesis",3)),
  stringsAsFactors = FALSE)
write.xlsx(markers, file.path(SUPP_DIR, "S2_marker_genes.xlsx"))

# S3 — DEG lists per study (data/degs is the authoritative base)
s3 <- list()
for (s in 1:8) {
    up <- read_excel(file.path(DEGS_DIR, paste0(s, "_Upregulated_Genes.xlsx")));   up$direction <- "up"
    dn <- read_excel(file.path(DEGS_DIR, paste0(s, "_Downregulated_Genes.xlsx"))); dn$direction <- "down"
    s3[[paste0("Study_", s)]] <- bind_rows(up, dn)
}
save_book(s3, file.path(SUPP_DIR, "S3_DEG_lists.xlsx"))

# S4 / S5 — DEG KEGG / GO enrichment
save_book(sheets_from_csvs("^deg_kegg_enrichment_.*csv$", "deg_kegg_enrichment_"), file.path(SUPP_DIR, "S4_DEG_KEGG_enrichment.xlsx"))
save_book(sheets_from_csvs("^deg_go_enrichment_.*csv$",   "deg_go_enrichment_"),   file.path(SUPP_DIR, "S5_DEG_GO_enrichment.xlsx"))

# S6 — rewired proteins
write.xlsx(read_csv(file.path(REWIRED_DIR, "all_rewired_proteins.csv"), show_col_types = FALSE),
           file.path(SUPP_DIR, "S6_rewired_proteins.xlsx"))

# S7 / S8 — rewired KEGG / GO enrichment
save_book(sheets_from_csvs("^kegg_enrichment_.*csv$", "kegg_enrichment_"), file.path(SUPP_DIR, "S7_rewired_KEGG_enrichment.xlsx"))
save_book(sheets_from_csvs("^go_enrichment_.*csv$",   "go_enrichment_"),   file.path(SUPP_DIR, "S8_rewired_GO_enrichment.xlsx"))

# S9 — Louvain modules
write.xlsx(read_csv(file.path(MODULE_DIR, "full_module_membership.csv"), show_col_types = FALSE),
           file.path(SUPP_DIR, "S9_modules.xlsx"))

cat("Supplementary tables S1-S9 written to", SUPP_DIR, "\n")
