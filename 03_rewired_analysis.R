#!/usr/bin/env Rscript
# 03_rewired_analysis.R
# Identifies "Rewired" proteins: high connectivity change, stable expression
# Rewired = proteins with significant PPI changes but NOT DEGs
#
# Definition: A protein is "rewired" if it meets both criteria:
#   1. High differential degree — in the top 25% (DEGREE_PERCENTILE) of proteins
#      by total number of repressed + activated interactions (Ntotal).
#   2. NOT a DEG — the gene is not differentially expressed (|log2FC| > 1, padj < 0.05).
# Rationale: These proteins change their interaction partners without transcriptional
# change, suggesting post-transcriptional or structural rewiring of the interactome.

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
    library(readxl)
    library(openxlsx)
    library(ggplot2)
    library(tidyr)
})

# --- Load Orthology Mappings ---
study_org_map <- read_csv(file.path(RESOURCE_DIR, "study_organism_mapping.csv"), show_col_types = FALSE)
fox_mapping <- read_csv(file.path(RESOURCE_DIR, "fox_to_fgr_orthology.csv"), show_col_types = FALSE)
fpu_mapping <- read_csv(file.path(RESOURCE_DIR, "fpu_to_fgr_orthology.csv"), show_col_types = FALSE)

# --- Helper: Load DEG data for a study and return gene IDs with direction ---
load_deg_data <- function(study_id) {
    up_file <- file.path(DEGS_DIR, paste0(study_id, "_Upregulated_Genes.xlsx"))
    down_file <- file.path(DEGS_DIR, paste0(study_id, "_Downregulated_Genes.xlsx"))

    deg_list <- list()

    if (file.exists(up_file)) {
        df <- read_excel(up_file)
        names(df) <- tolower(names(df))
        if ("gene" %in% names(df) && "log2foldchange" %in% names(df)) {
            deg_list[[1]] <- df %>%
                select(gene_id = gene, logFC = log2foldchange) %>%
                mutate(direction = "up")
        }
    }

    if (file.exists(down_file)) {
        df <- read_excel(down_file)
        names(df) <- tolower(names(df))
        if ("gene" %in% names(df) && "log2foldchange" %in% names(df)) {
            deg_list[[2]] <- df %>%
                select(gene_id = gene, logFC = log2foldchange) %>%
                mutate(direction = "down")
        }
    }

    if (length(deg_list) > 0) {
        return(bind_rows(deg_list))
    }
    return(data.frame(gene_id = character(), logFC = numeric(), direction = character()))
}

# --- Helper: Map gene IDs to FGR (local version using loaded mappings) ---
map_genes_to_fgr <- function(gene_ids, organism) {
    if (organism == "fgr") {
        return(data.frame(original_id = gene_ids, fgsg_id = gene_ids, stringsAsFactors = FALSE))
    } else if (organism == "fox") {
        mapping <- fox_mapping %>%
            filter(FOXG_ID %in% gene_ids) %>%
            select(original_id = FOXG_ID, fgsg_id = FGSG_ID)
        return(mapping)
    } else if (organism == "fpu") {
        mapping <- fpu_mapping %>%
            filter(FPSE_ID %in% gene_ids) %>%
            select(original_id = FPSE_ID, fgsg_id = FGSG_ID)
        return(mapping)
    }
    return(data.frame(original_id = character(), fgsg_id = character()))
}

# --- Main Analysis ---
all_rewired <- list()
all_deg_sets <- list()

for (study_id in 1:8) {
    # Get study info
    organism <- study_org_map$organism[study_org_map$study_id == study_id]

    # Load differential interactome results
    di_file <- file.path(DI_DIR, paste0("Study_", study_id), "protein_stats.csv")
    if (!file.exists(di_file)) {
        next
    }

    prot_stats <- read_csv(di_file, show_col_types = FALSE)

    # Load DEG data
    deg_data <- load_deg_data(study_id)

    if (nrow(deg_data) == 0) {
        deg_fgsg <- character()
    } else {
        # Map DEGs to FGR
        deg_mapping <- map_genes_to_fgr(deg_data$gene_id, organism)
        deg_fgsg <- unique(deg_mapping$fgsg_id)
    }

    all_deg_sets[[study_id]] <- deg_fgsg

    # Calculate differential degree threshold
    if (nrow(prot_stats) == 0) next

    degree_threshold <- quantile(prot_stats$Ntotal, DEGREE_PERCENTILE, na.rm = TRUE)

    # Identify rewired proteins:
    # 1. High differential degree (top percentile)
    # 2. NOT in DEG list (stable expression)
    high_degree_prots <- prot_stats %>%
        filter(Ntotal >= degree_threshold)

    rewired <- high_degree_prots %>%
        filter(!(Protein %in% deg_fgsg))

    # Also get non-rewired high-degree proteins (DEGs with high connectivity)
    deg_high_degree <- high_degree_prots %>%
        filter(Protein %in% deg_fgsg)

    if (nrow(rewired) > 0) {
        rewired$Study <- study_id

        # Load interaction data to classify rewired proteins
        int_file <- file.path(DI_DIR, paste0("Study_", study_id), "all_interactions.csv")
        if (file.exists(int_file)) {
            interactions <- read_csv(int_file, show_col_types = FALSE)

            # Classify rewired proteins by dominant interaction type
            rewired <- rewired %>%
                rowwise() %>%
                mutate(
                    Classification = case_when(
                        Nrep > Nact ~ "Mostly_Repressed",
                        Nact > Nrep ~ "Mostly_Activated",
                        TRUE ~ "Balanced"
                    )
                ) %>%
                ungroup()
        }

        all_rewired[[study_id]] <- rewired
    }

    # Save per-study results
    study_output <- file.path(REWIRED_DIR, paste0("Study_", study_id))
    if (!dir.exists(study_output)) dir.create(study_output)

    if (nrow(rewired) > 0) {
        write_csv(rewired, file.path(study_output, "rewired_proteins.csv"))
    }
    write_csv(high_degree_prots, file.path(study_output, "high_degree_proteins.csv"))
}

# --- Aggregate Results ---
if (length(all_rewired) > 0) {
    combined_rewired <- bind_rows(all_rewired)

    # Add BCA type from STUDY_INFO
    combined_rewired <- combined_rewired %>%
        mutate(BCA_Type = STUDY_INFO$bca_type[match(Study, STUDY_INFO$study_id)])

    write_csv(combined_rewired, file.path(REWIRED_DIR, "all_rewired_proteins.csv"))

    # Cross-study analysis: count occurrences of each protein across studies
    protein_frequency <- combined_rewired %>%
        group_by(Protein, BCA_Type) %>%
        summarise(
            n_studies = n(),
            studies = paste(unique(Study), collapse = ","),
            mean_degree = mean(Ntotal),
            mean_repressed = mean(Nrep),
            mean_activated = mean(Nact),
            .groups = "drop"
        ) %>%
        arrange(desc(n_studies), desc(mean_degree))

    write_csv(protein_frequency, file.path(REWIRED_DIR, "rewired_protein_frequency.csv"))

    # Separate by BCA type
    fungal_rewired <- combined_rewired %>% filter(BCA_Type == "Fungal")
    bacterial_rewired <- combined_rewired %>% filter(BCA_Type == "Bacterial")

    write_csv(fungal_rewired, file.path(REWIRED_DIR, "fungal_rewired_proteins.csv"))
    write_csv(bacterial_rewired, file.path(REWIRED_DIR, "bacterial_rewired_proteins.csv"))

    # Unique proteins per BCA type
    fungal_unique <- unique(fungal_rewired$Protein)
    bacterial_unique <- unique(bacterial_rewired$Protein)

    # Overlap analysis
    overlap <- intersect(fungal_unique, bacterial_unique)
    fungal_specific <- setdiff(fungal_unique, bacterial_unique)
    bacterial_specific <- setdiff(bacterial_unique, fungal_unique)

    # Save overlap results
    overlap_df <- data.frame(
        Category = c("Shared", "Fungal_Specific", "Bacterial_Specific"),
        Count = c(length(overlap), length(fungal_specific), length(bacterial_specific)),
        Proteins = c(
            paste(overlap, collapse = ";"),
            paste(fungal_specific, collapse = ";"),
            paste(bacterial_specific, collapse = ";")
        ),
        stringsAsFactors = FALSE
    )
    write_csv(overlap_df, file.path(REWIRED_DIR, "rewired_overlap_summary.csv"))

    # Detailed lists
    write_csv(data.frame(Protein = overlap), file.path(REWIRED_DIR, "rewired_shared.csv"))
    write_csv(data.frame(Protein = fungal_specific), file.path(REWIRED_DIR, "rewired_fungal_specific.csv"))
    write_csv(data.frame(Protein = bacterial_specific), file.path(REWIRED_DIR, "rewired_bacterial_specific.csv"))

    # Per-study rewired protein summary
    per_study_summary <- combined_rewired %>%
        group_by(Study, BCA_Type) %>%
        summarise(
            n_rewired = n(),
            n_mostly_repressed = sum(Classification == "Mostly_Repressed", na.rm = TRUE),
            n_mostly_activated = sum(Classification == "Mostly_Activated", na.rm = TRUE),
            n_balanced = sum(Classification == "Balanced", na.rm = TRUE),
            mean_degree = mean(Ntotal),
            .groups = "drop"
        )

    write_csv(per_study_summary, file.path(REWIRED_DIR, "per_study_rewired_summary.csv"))

    # --- Cross-Study Behavior Matrix ---
    all_proteins <- unique(combined_rewired$Protein)
    behavior_matrix <- data.frame(Protein = all_proteins)

    for (sid in 1:8) {
        if (sid %in% names(all_rewired)) {
            study_rewired <- all_rewired[[sid]]
            behavior_matrix[[paste0("Study_", sid, "_Degree")]] <-
                study_rewired$Ntotal[match(all_proteins, study_rewired$Protein)]
            behavior_matrix[[paste0("Study_", sid, "_Nrep")]] <-
                study_rewired$Nrep[match(all_proteins, study_rewired$Protein)]
            behavior_matrix[[paste0("Study_", sid, "_Nact")]] <-
                study_rewired$Nact[match(all_proteins, study_rewired$Protein)]
        } else {
            behavior_matrix[[paste0("Study_", sid, "_Degree")]] <- NA
            behavior_matrix[[paste0("Study_", sid, "_Nrep")]] <- NA
            behavior_matrix[[paste0("Study_", sid, "_Nact")]] <- NA
        }
    }

    # Add summary columns
    behavior_matrix <- behavior_matrix %>%
        rowwise() %>%
        mutate(
            n_studies_rewired = sum(!is.na(c_across(ends_with("_Degree")))),
            fungal_count = sum(!is.na(c_across(matches("Study_[1-4]_Degree")))),
            bacterial_count = sum(!is.na(c_across(matches("Study_[5-8]_Degree"))))
        ) %>%
        ungroup() %>%
        arrange(desc(n_studies_rewired))

    write_csv(behavior_matrix, file.path(REWIRED_DIR, "cross_study_behavior_matrix.csv"))

    # Create supplementary table for publication
    supp_table <- behavior_matrix %>%
        select(Protein, n_studies_rewired, fungal_count, bacterial_count) %>%
        mutate(
            Specificity = case_when(
                fungal_count > 0 & bacterial_count == 0 ~ "Fungal_Specific",
                bacterial_count > 0 & fungal_count == 0 ~ "Bacterial_Specific",
                fungal_count > 0 & bacterial_count > 0 ~ "Shared",
                TRUE ~ "Unknown"
            )
        )

    write_csv(supp_table, file.path(REWIRED_DIR, "supplementary_table_rewired_summary.csv"))

    print(per_study_summary)

}
