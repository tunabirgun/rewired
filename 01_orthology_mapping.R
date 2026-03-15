#!/usr/bin/env Rscript
# 01_orthology_mapping.R
# Creates orthology mappings using KEGG KO identifiers
# Maps F. oxysporum (FOXG) and F. pseudograminearum (FPSE) to F. graminearum (FGSG)

# Source shared configuration
local({
    args <- commandArgs(trailingOnly = FALSE)
    f <- grep("--file=", args, value = TRUE)
    d <- if (length(f)) dirname(normalizePath(sub("--file=", "", f))) else getwd()
    source(file.path(d, "00_setup.R"))
})

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
})

# --- Load KEGG KO Files ---
load_kegg_ko <- function(filepath, prefix) {
    df <- read_tsv(filepath, col_names = c("gene", "ko"), col_types = "cc", show_col_types = FALSE)
    df$gene <- gsub(paste0("^", prefix, ":"), "", df$gene)
    df$ko <- gsub("^ko:", "", df$ko)
    return(df)
}

fgr_ko <- load_kegg_ko(file.path(STRING_DIR, "fgr_ko.txt"), "fgr")
fox_ko <- load_kegg_ko(file.path(STRING_DIR, "fox_ko.txt"), "fox")
fpu_ko <- load_kegg_ko(file.path(STRING_DIR, "fpu_ko.txt"), "fpu")

# --- Create Orthology Mappings ---

# FOX to FGR mapping via KO
fox_to_fgr <- fox_ko %>%
    inner_join(fgr_ko, by = "ko", suffix = c("_fox", "_fgr")) %>%
    select(FOXG_ID = gene_fox, FGSG_ID = gene_fgr, KO = ko) %>%
    distinct()

# Handle one-to-many mappings: pick one FGR ortholog per FOXG.
# Strategy: select the lowest FGSG numeric ID for deterministic, reproducible
# arbitration. This avoids non-determinism from tibble row order and ensures
# consistent results across runs. The choice of "lowest" is arbitrary but
# stable; biological equivalence is assumed among KO-linked paralogs.
fox_to_fgr_unique <- fox_to_fgr %>%
    mutate(fgsg_num = as.numeric(gsub("FGSG_", "", FGSG_ID))) %>%
    group_by(FOXG_ID) %>%
    slice_min(fgsg_num, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(FOXG_ID, FGSG_ID, KO)

# FPU to FGR mapping via KO
fpu_to_fgr <- fpu_ko %>%
    inner_join(fgr_ko, by = "ko", suffix = c("_fpu", "_fgr")) %>%
    select(FPSE_ID = gene_fpu, FGSG_ID = gene_fgr, KO = ko) %>%
    distinct()

# Handle one-to-many mappings (same deterministic strategy as above)
fpu_to_fgr_unique <- fpu_to_fgr %>%
    mutate(fgsg_num = as.numeric(gsub("FGSG_", "", FGSG_ID))) %>%
    group_by(FPSE_ID) %>%
    slice_min(fgsg_num, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(FPSE_ID, FGSG_ID, KO)

# --- Load STRING Alias Files for ID Conversion ---
load_string_aliases <- function(filepath, taxid) {
    df <- read_tsv(filepath, col_names = c("string_id", "alias", "source"),
                   col_types = "ccc", comment = "#", show_col_types = FALSE)
    # Extract gene IDs (FGSG, FOXG, FPSE)
    kegg_short <- df %>%
        filter(source == "KEGG_KEGGID_SHORT") %>%
        select(string_id, gene_id = alias)
    return(kegg_short)
}

fgr_aliases <- load_string_aliases(file.path(STRING_DIR, "229533.protein.aliases.v12.0.txt"), "229533")
fox_aliases <- load_string_aliases(file.path(STRING_DIR, "426428.protein.aliases.v12.0.txt"), "426428")

# --- Create Complete Mapping Tables ---

# FGR: STRING ID to FGSG
fgr_complete <- fgr_aliases %>%
    rename(FGSG_ID = gene_id, FGR_STRING_ID = string_id)

# FOX: STRING ID to FGSG (via KO orthology)
fox_complete <- fox_aliases %>%
    rename(FOXG_ID = gene_id, FOX_STRING_ID = string_id) %>%
    left_join(fox_to_fgr_unique %>% select(FOXG_ID, FGSG_ID), by = "FOXG_ID")

# --- Save Mapping Files ---
write_csv(fgr_complete, file.path(RESOURCE_DIR, "fgr_string_to_gene.csv"))
write_csv(fox_complete, file.path(RESOURCE_DIR, "fox_string_to_fgr.csv"))
write_csv(fpu_to_fgr_unique, file.path(RESOURCE_DIR, "fpu_to_fgr_orthology.csv"))
write_csv(fox_to_fgr_unique, file.path(RESOURCE_DIR, "fox_to_fgr_orthology.csv"))

# Create a unified mapping table for all studies
# Study 1-3: F. graminearum (direct)
# Study 4,5,7,8: F. oxysporum (use fox_to_fgr)
# Study 6: F. pseudograminearum (use fpu_to_fgr)
study_organism_map <- data.frame(
    study_id = STUDY_INFO$study_id,
    organism = STUDY_INFO$org_code,
    taxid = c("229533", "229533", "229533", "426428", "426428", "1028729", "426428", "426428"),
    stringsAsFactors = FALSE
)
write_csv(study_organism_map, file.path(RESOURCE_DIR, "study_organism_mapping.csv"))

# --- Summary Statistics ---
summary_stats <- data.frame(
    Organism = c("F. graminearum", "F. oxysporum", "F. pseudograminearum"),
    Total_Genes_KO = c(n_distinct(fgr_ko$gene), n_distinct(fox_ko$gene), n_distinct(fpu_ko$gene)),
    Genes_Mapped_to_FGR = c(n_distinct(fgr_ko$gene),
                            n_distinct(fox_to_fgr_unique$FOXG_ID),
                            n_distinct(fpu_to_fgr_unique$FPSE_ID)),
    stringsAsFactors = FALSE
)
summary_stats$Mapping_Rate <- round(100 * summary_stats$Genes_Mapped_to_FGR / summary_stats$Total_Genes_KO, 1)
write_csv(summary_stats, file.path(RESOURCE_DIR, "orthology_summary_stats.csv"))

mapped_fox <- sum(!is.na(fox_complete$FGSG_ID))
