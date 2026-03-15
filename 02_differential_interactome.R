#!/usr/bin/env Rscript
# 02_differential_interactome.R
# Differential Interactome Analysis using F. graminearum as reference
# Maps all studies to FGR orthology before running the algorithm
#
# Algorithm reference: Gulfidan et al. (2020), PLoS ONE 15(3):e0228187
#
# Steps:
#   1. Binarize expression: for each sample, genes with expression >= mean + 1 SD
#      are set to +1 (high), <= mean - 1 SD to -1 (low), else 0 (normal).
#   2. For each PPI edge and each of 9 co-expression states (s1, s2) in {-1,0,+1}^2,
#      count frequencies in control (NNC) and treatment (NND) groups.
#   3. Q-value = NND / (NNC + NND); Q=0.5 if both zero (no change).
#   4. Filter: Q <= 0.10 (repressed) or Q >= 0.90 (activated), AND
#      pmax(NNC, NND) >= DELTA (minimum frequency threshold to exclude noise).
# PPI source: STRING v12.0, F. graminearum (taxid 229533), combined_score >= 400

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
    library(tictoc)
    library(tibble)
})

# --- Load Resources ---
METADATA_FILE <- file.path(METADATA_DIR, "study_metadata.csv")
study_meta <- read_csv(METADATA_FILE, show_col_types = FALSE)

# Study-organism mapping
study_org_map <- read_csv(file.path(RESOURCE_DIR, "study_organism_mapping.csv"), show_col_types = FALSE)

# Orthology mappings
fgr_mapping <- read_csv(file.path(RESOURCE_DIR, "fgr_string_to_gene.csv"), show_col_types = FALSE)
fox_mapping <- read_csv(file.path(RESOURCE_DIR, "fox_string_to_fgr.csv"), show_col_types = FALSE)
fpu_mapping <- read_csv(file.path(RESOURCE_DIR, "fpu_to_fgr_orthology.csv"), show_col_types = FALSE)

# F. graminearum PPI network (reference)
ppi_raw <- read_delim(file.path(STRING_DIR, "229533.protein.links.v12.0.txt"),
                      delim = " ", show_col_types = FALSE)
# Filter for high-confidence interactions (score >= 400)
ppi_filtered <- ppi_raw %>%
    filter(combined_score >= 400) %>%
    mutate(
        protein1 = gsub("^229533\\.", "", protein1),
        protein2 = gsub("^229533\\.", "", protein2)
    )

# Map STRING IDs to FGSG IDs
string_to_fgsg <- fgr_mapping %>%
    mutate(string_id_clean = gsub("^229533\\.", "", FGR_STRING_ID)) %>%
    select(string_id_clean, FGSG_ID) %>%
    distinct()

ppi_network <- ppi_filtered %>%
    inner_join(string_to_fgsg, by = c("protein1" = "string_id_clean")) %>%
    rename(ProteinA = FGSG_ID) %>%
    inner_join(string_to_fgsg, by = c("protein2" = "string_id_clean")) %>%
    rename(ProteinB = FGSG_ID) %>%
    select(ProteinA, ProteinB, combined_score) %>%
    filter(ProteinA != ProteinB) %>%
    distinct()

# --- Helper: Map gene IDs to FGR orthologs (local version using loaded mappings) ---
map_to_fgr_local <- function(gene_ids, organism) {
    if (organism == "fgr") {
        return(data.frame(original_id = gene_ids, fgsg_id = gene_ids, stringsAsFactors = FALSE))
    } else if (organism == "fox") {
        mapping <- fox_mapping %>%
            select(original_id = FOXG_ID, fgsg_id = FGSG_ID) %>%
            filter(original_id %in% gene_ids)
        return(mapping)
    } else if (organism == "fpu") {
        mapping <- fpu_mapping %>%
            select(original_id = FPSE_ID, fgsg_id = FGSG_ID) %>%
            filter(original_id %in% gene_ids)
        return(mapping)
    }
    return(data.frame(original_id = character(), fgsg_id = character(), stringsAsFactors = FALSE))
}

# --- Core Differential Interactome Algorithm ---
run_diff_interactome <- function(expdata, ppi_df, X0, X1) {
    # expdata: matrix with genes as rows, samples as columns
    # ppi_df: data frame with ProteinA, ProteinB columns (FGSG IDs)
    # X0: number of control samples
    # X1: number of case samples

    if (nrow(expdata) == 0 || nrow(ppi_df) == 0) {
        return(NULL)
    }

    genes <- rownames(expdata)

    # Map PPI to row indices
    INT <- ppi_df %>%
        mutate(
            idx_a = match(ProteinA, genes),
            idx_b = match(ProteinB, genes)
        ) %>%
        filter(!is.na(idx_a) & !is.na(idx_b)) %>%
        select(ProteinA, ProteinB, idx_a, idx_b)

    if (nrow(INT) == 0) {
        return(NULL)
    }

    # Binarization: +1 (high), 0 (normal), -1 (low)
    Data <- as.matrix(expdata)
    AvgExp <- colMeans(Data)
    StdDev <- apply(Data, 2, sd)

    NData <- matrix(0, nrow = nrow(Data), ncol = ncol(Data))
    for (i in 1:ncol(Data)) {
        # Guard: if a sample has zero variance (SD=0), all genes are at the mean,
        # so binarization yields all zeros (normal state)
        if (StdDev[i] == 0) next
        NData[, i] <- ifelse(Data[, i] >= AvgExp[i] + StdDev[i], 1,
                             ifelse(Data[, i] <= AvgExp[i] - StdDev[i], -1, 0))
    }

    NCData <- NData[, 1:X0, drop = FALSE]
    NDData <- NData[, (X0 + 1):(X0 + X1), drop = FALSE]

    # Calculate interaction state frequencies
    states <- expand.grid(s1 = c(-1, 0, 1), s2 = c(-1, 0, 1))
    state_names <- paste0("[", states$s1, " ", states$s2, "]")

    calc_freq <- function(idx_a, idx_b, data_mat, s1, s2) {
        sum(data_mat[idx_a, ] == s1 & data_mat[idx_b, ] == s2)
    }

    results_list <- list()

    for (i in 1:nrow(states)) {
        s1 <- states$s1[i]
        s2 <- states$s2[i]
        state_name <- state_names[i]

        NC <- sapply(1:nrow(INT), function(j) calc_freq(INT$idx_a[j], INT$idx_b[j], NCData, s1, s2))
        ND <- sapply(1:nrow(INT), function(j) calc_freq(INT$idx_a[j], INT$idx_b[j], NDData, s1, s2))

        NNC <- NC / X0
        NND <- ND / X1

        # Q-value calculation
        q_value <- ifelse((NNC + NND) != 0, NND / (NNC + NND), 0.5)

        res_df <- data.frame(
            ProteinA = INT$ProteinA,
            ProteinB = INT$ProteinB,
            NN = NC,
            NT = ND,
            q_value = q_value,
            State = state_name,
            NNC = NNC,
            NND = NND,
            stringsAsFactors = FALSE
        )

        # Filter significant interactions
        sig_df <- res_df %>%
            filter((q_value <= QCRIT1 | q_value >= QCRIT2) & pmax(NNC, NND) >= DELTA)

        results_list[[i]] <- sig_df
    }

    RESS <- bind_rows(results_list)
    RESS$q_value <- signif(RESS$q_value, digits = 3)
    RESS <- RESS %>% arrange(q_value)

    return(RESS)
}

# --- Main Processing Loop ---
all_results <- list()

for (study_id in 1:8) {
    # Get study info
    study_info <- study_meta %>% filter(Study_No == study_id)
    organism <- study_org_map$organism[study_org_map$study_id == study_id]

    # Load counts
    counts_file <- file.path(COUNTS_DIR, paste0(study_id, ".txt"))
    if (!file.exists(counts_file)) {
        next
    }

    counts_raw <- read_tsv(counts_file, comment = "#", show_col_types = FALSE)

    # Extract gene IDs and count columns
    gene_col <- which(colnames(counts_raw) == "Geneid")
    if (length(gene_col) == 0) gene_col <- 1
    count_cols <- (gene_col + 5 + 1):ncol(counts_raw)  # Skip Geneid, Chr, Start, End, Strand, Length

    gene_ids <- counts_raw[[gene_col]]
    count_matrix <- as.matrix(counts_raw[, count_cols])
    rownames(count_matrix) <- gene_ids

    # Map to FGR orthologs
    gene_mapping <- map_to_fgr_local(gene_ids, organism)

    if (nrow(gene_mapping) == 0) {
        next
    }

    # Create mapped expression matrix
    mapped_idx <- match(gene_mapping$original_id, gene_ids)
    mapped_matrix <- count_matrix[mapped_idx, , drop = FALSE]
    rownames(mapped_matrix) <- gene_mapping$fgsg_id

    # Handle duplicate mappings (multiple genes -> same ortholog): sum counts
    if (anyDuplicated(rownames(mapped_matrix))) {
        mapped_df <- as.data.frame(mapped_matrix)
        mapped_df$gene <- rownames(mapped_matrix)
        mapped_df <- mapped_df %>%
            group_by(gene) %>%
            summarise(across(everything(), sum)) %>%
            ungroup()
        mapped_matrix <- as.matrix(mapped_df[, -1])
        rownames(mapped_matrix) <- mapped_df$gene
    }

    # Determine control and case sample indices
    col_names <- colnames(mapped_matrix)

    # Parse SRR IDs from column names
    srr_ids <- gsub(".*?(SRR[0-9]+).*", "\\1", col_names)

    # Load detailed metadata to identify controls
    detailed_meta <- read_csv(file.path(METADATA_DIR, "study_table.csv"), show_col_types = FALSE)
    study_detailed <- detailed_meta %>%
        filter(`Study No.` == study_id)

    control_srrs <- study_detailed$`SRA Accession`[study_detailed$`Control/Sample` == "Control"]
    sample_srrs <- study_detailed$`SRA Accession`[study_detailed$`Control/Sample` == "Sample"]

    control_idx <- which(srr_ids %in% control_srrs)
    sample_idx <- which(srr_ids %in% sample_srrs)

    if (length(control_idx) == 0 || length(sample_idx) == 0) {
        next
    }

    # Reorder: controls first, then samples
    reordered_matrix <- mapped_matrix[, c(control_idx, sample_idx), drop = FALSE]
    X0 <- length(control_idx)
    X1 <- length(sample_idx)

    # Run differential interactome
    tic()
    results <- run_diff_interactome(reordered_matrix, ppi_network, X0, X1)
    toc()

    if (is.null(results) || nrow(results) == 0) {
        next
    }

    # Classify interactions
    repressed <- results %>% filter(q_value <= QCRIT1)
    activated <- results %>% filter(q_value >= QCRIT2)

    # Calculate protein-level statistics
    all_prots <- unique(c(results$ProteinA, results$ProteinB))
    prot_stats <- data.frame(Protein = all_prots, stringsAsFactors = FALSE) %>%
        rowwise() %>%
        mutate(
            Nrep = sum(repressed$ProteinA == Protein) + sum(repressed$ProteinB == Protein),
            Nact = sum(activated$ProteinA == Protein) + sum(activated$ProteinB == Protein),
            Ntotal = Nrep + Nact
        ) %>%
        ungroup() %>%
        arrange(desc(Ntotal))

    # Save results
    study_output_dir <- file.path(DI_DIR, paste0("Study_", study_id))
    if (!dir.exists(study_output_dir)) dir.create(study_output_dir)

    wb <- createWorkbook()

    addWorksheet(wb, "All_Interactions")
    writeData(wb, "All_Interactions", results)

    addWorksheet(wb, "Repressed")
    writeData(wb, "Repressed", repressed)

    addWorksheet(wb, "Activated")
    writeData(wb, "Activated", activated)

    addWorksheet(wb, "Protein_Stats")
    writeData(wb, "Protein_Stats", prot_stats)

    saveWorkbook(wb, file.path(study_output_dir, "DiffInt_Results.xlsx"), overwrite = TRUE)

    # Save CSVs for downstream analysis
    write_csv(results, file.path(study_output_dir, "all_interactions.csv"))
    write_csv(repressed, file.path(study_output_dir, "repressed_interactions.csv"))
    write_csv(activated, file.path(study_output_dir, "activated_interactions.csv"))
    write_csv(prot_stats, file.path(study_output_dir, "protein_stats.csv"))

    # Store for aggregation
    all_results[[as.character(study_id)]] <- list(
        study_id = study_id,
        bca_type = study_info$BCA_Type[1],
        all = results,
        repressed = repressed,
        activated = activated,
        prot_stats = prot_stats
    )
}

# --- Aggregate Results ---

# Combine all protein stats
all_prot_stats <- lapply(names(all_results), function(id) {
    df <- all_results[[id]]$prot_stats
    df$Study <- id
    df$BCA_Type <- all_results[[id]]$bca_type
    return(df)
}) %>% bind_rows()

# Combine all interactions
all_interactions <- lapply(names(all_results), function(id) {
    df <- all_results[[id]]$all
    df$Study <- id
    df$BCA_Type <- all_results[[id]]$bca_type
    return(df)
}) %>% bind_rows()

# Save aggregated results
write_csv(all_prot_stats, file.path(DI_DIR, "all_studies_protein_stats.csv"))
write_csv(all_interactions, file.path(DI_DIR, "all_studies_interactions.csv"))

# Summary table
summary_df <- data.frame(
    Study = names(all_results),
    BCA_Type = sapply(all_results, function(x) x$bca_type),
    Total_Interactions = sapply(all_results, function(x) nrow(x$all)),
    Repressed = sapply(all_results, function(x) nrow(x$repressed)),
    Activated = sapply(all_results, function(x) nrow(x$activated)),
    Unique_Proteins = sapply(all_results, function(x) nrow(x$prot_stats)),
    stringsAsFactors = FALSE
)

write_csv(summary_df, file.path(DI_DIR, "summary_statistics.csv"))

print(summary_df)
