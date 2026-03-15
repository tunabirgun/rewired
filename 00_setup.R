#!/usr/bin/env Rscript
# 00_setup.R
# Shared configuration sourced by all pipeline scripts
# Defines paths, packages, STUDY_INFO, and utility functions
#
# Reference: Differential interactome methodology adapted from
# Gulfidan et al. (2020), PLoS ONE 15(3):e0228187
# STRING database: v12.0 (https://string-db.org)

# --- Path Resolution ---
get_script_dir <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
        return(dirname(normalizePath(sub("--file=", "", file_arg))))
    }
    if (exists(".script_dir_override")) return(.script_dir_override)
    return(getwd())
}

SCRIPT_DIR <- get_script_dir()
BASE_DIR   <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)

# Validate project root — require counts/ or degs/ subdirectory
if (!dir.exists(file.path(BASE_DIR, "counts")) && !dir.exists(file.path(BASE_DIR, "degs"))) {
    stop("Cannot locate project root. Expected 'counts/' or 'degs/' under: ", BASE_DIR,
         "\nRun scripts from the publication_scripts/ directory or set .script_dir_override.")
}

# Directory tree
COUNTS_DIR    <- file.path(BASE_DIR, "counts")
DEGS_DIR      <- file.path(BASE_DIR, "degs")
METADATA_DIR  <- file.path(BASE_DIR, "metadata")
STRING_DIR    <- file.path(BASE_DIR, "string_data")
RESOURCE_DIR  <- file.path(BASE_DIR, "resources")
OUTPUT_DIR    <- file.path(BASE_DIR, "outputs")
DI_DIR        <- file.path(OUTPUT_DIR, "differential_interactome")
DE_DIR        <- file.path(OUTPUT_DIR, "differential_expression")
REWIRED_DIR   <- file.path(OUTPUT_DIR, "rewired_analysis")
ENRICH_DIR    <- file.path(OUTPUT_DIR, "enrichment")
FIG_DIR       <- file.path(OUTPUT_DIR, "figures")
TABLE_DIR     <- file.path(OUTPUT_DIR, "tables")
NETWORK_DIR   <- file.path(OUTPUT_DIR, "networks")
MODULE_DIR    <- file.path(OUTPUT_DIR, "modules")

# Create output directories
for (d in c(OUTPUT_DIR, DI_DIR, DE_DIR, REWIRED_DIR, ENRICH_DIR,
            FIG_DIR, TABLE_DIR, NETWORK_DIR, MODULE_DIR)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# --- Study Metadata ---
STUDY_INFO <- data.frame(
    study_id   = 1:8,
    organism   = c("F. graminearum", "F. graminearum", "F. graminearum",
                   "F. oxysporum", "F. oxysporum", "F. pseudograminearum",
                   "F. oxysporum", "F. oxysporum"),
    bca        = c("T. gamsii", "T. atroviride", "T. hamatum",
                   "A. pullulans", "P. aeruginosa", "B. velezensis",
                   "B. subtilis", "B. amyloliquefaciens"),
    bca_type   = c("Fungal", "Fungal", "Fungal", "Fungal",
                   "Bacterial", "Bacterial", "Bacterial", "Bacterial"),
    org_code   = c("fgr", "fgr", "fgr", "fox", "fox", "fpu", "fox", "fox"),
    gene_prefix = c("FGSG", "FGSG", "FGSG", "FOXG", "FOXG", "FPSE", "FOXG", "FOXG"),
    stringsAsFactors = FALSE
)

# --- Thresholds ---
# Q-value and DELTA thresholds from Gulfidan et al. (2020), PLoS ONE 15(3):e0228187
LOG2FC_THRESHOLD <- 1.0
PADJ_THRESHOLD   <- 0.05
QCRIT1           <- 0.10   # Repression threshold: Q <= 0.10 indicates control-dominant state
QCRIT2           <- 0.90   # Activation threshold: Q >= 0.90 indicates treatment-dominant state
DELTA            <- 0.20   # Minimum frequency: pmax(NNC, NND) >= DELTA filters noise
DEGREE_PERCENTILE <- 0.75  # Top 25% differential degree

# Global seed for reproducibility (affects Louvain clustering, network layouts, etc.)
set.seed(42)

# --- Color Palettes ---
COLORS_FUNGAL    <- c("#4DAF4A", "#377EB8", "#984EA3", "#FF7F00")
COLORS_BACTERIAL <- c("#E41A1C", "#F781BF", "#A65628", "#999999")
COLORS_BCA       <- c(Fungal = "#2E7D32", Bacterial = "#C62828")
COLOR_ACTIVATED  <- "#E41A1C"
COLOR_REPRESSED  <- "#377EB8"
COLOR_HUB        <- "#FF7F00"
COLOR_REGULAR    <- "#999999"
SIG_COLORS       <- c(NS = "grey70", Up = "#D55E00", Down = "#0072B2")

# --- Utility Functions ---

# Save a figure as PNG + SVG + description file
save_fig <- function(plot_obj, filename, width = 10, height = 8,
                     title = "", description = "", is_grid = FALSE) {
    png_path <- file.path(FIG_DIR, paste0(filename, ".png"))
    svg_path <- file.path(FIG_DIR, paste0(filename, ".svg"))

    png(png_path, width = width, height = height, units = "in", res = 300)
    if (is_grid) grid::grid.draw(plot_obj) else print(plot_obj)
    dev.off()

    svglite::svglite(svg_path, width = width, height = height)
    if (is_grid) grid::grid.draw(plot_obj) else print(plot_obj)
    dev.off()

    if (nzchar(title) || nzchar(description)) {
        writeLines(c(
            paste("Title:", title), "",
            paste("Description:", description)
        ), file.path(FIG_DIR, paste0(filename, "_description.txt")))
    }
}

# Map gene IDs to FGR orthologs
map_to_fgr <- function(gene_ids, organism,
                       fox_mapping = NULL, fpu_mapping = NULL) {
    if (organism == "fgr") {
        return(data.frame(original_id = gene_ids, fgsg_id = gene_ids,
                          stringsAsFactors = FALSE))
    }
    if (organism == "fox" && !is.null(fox_mapping)) {
        return(fox_mapping %>%
                   dplyr::select(original_id = FOXG_ID, fgsg_id = FGSG_ID) %>%
                   dplyr::filter(original_id %in% gene_ids))
    }
    if (organism == "fpu" && !is.null(fpu_mapping)) {
        return(fpu_mapping %>%
                   dplyr::select(original_id = FPSE_ID, fgsg_id = FGSG_ID) %>%
                   dplyr::filter(original_id %in% gene_ids))
    }
    data.frame(original_id = character(), fgsg_id = character(),
               stringsAsFactors = FALSE)
}

# Load DEG gene lists from degs/ xlsx files
load_deg_genes <- function(study_id) {
    up_file   <- file.path(DEGS_DIR, paste0(study_id, "_Upregulated_Genes.xlsx"))
    down_file <- file.path(DEGS_DIR, paste0(study_id, "_Downregulated_Genes.xlsx"))
    genes <- character()
    for (f in c(up_file, down_file)) {
        if (file.exists(f)) {
            df <- readxl::read_excel(f)
            names(df) <- tolower(names(df))
            if ("gene" %in% names(df)) genes <- c(genes, df$gene)
        }
    }
    unique(genes)
}
