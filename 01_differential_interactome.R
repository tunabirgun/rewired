#!/usr/bin/env Rscript
# 01_differential_interactome.R
# Calculates the "Differential Interactome" based on PCC differences.
# Methodology Reference: Gulfidan et al. (2020). Sci Rep 10, 3272.

library(dplyr)
library(readr)
library(igraph)

# --- USER CONFIGURATION ---
EXPRESSION_FILE <- "path/to/normalized_expression_matrix.csv"  # Genes x Samples (e.g., TPM)
INTERACTOME_FILE <- "path/to/reference_interactome.csv"        # Columns: Protein_A, Protein_B
METADATA_FILE <- "path/to/sample_metadata.csv"                 # Columns: Sample_ID, Condition

# Output
OUTPUT_FILE <- "differential_interactome_scores.csv"

# --- 1. Load Data ---
cat("Loading data...\n")
expr_data <- read.csv(EXPRESSION_FILE, row.names = 1, check.names = FALSE)
interactome <- read.csv(INTERACTOME_FILE)
metadata <- read.csv(METADATA_FILE)

# Define Control and Treatment samples
control_samples <- metadata$Sample_ID[metadata$Condition == "Control"]
case_samples <- metadata$Sample_ID[metadata$Condition == "Treatment"]

if (length(control_samples) < 3 || length(case_samples) < 3) {
    stop("Error: Need at least 3 replicates per condition for reliable correlation.")
}

# --- 2. Define Correlation Function ---
# Calculates Pearson Correlation Coefficient (PCC) for a pair of genes
calc_correlation <- function(geneA, geneB, expression_matrix) {
    if (!geneA %in% rownames(expression_matrix) || !geneB %in% rownames(expression_matrix)) {
        return(NA)
    }
    vecA <- as.numeric(expression_matrix[geneA, ])
    vecB <- as.numeric(expression_matrix[geneB, ])
    cor(vecA, vecB, method = "pearson")
}

# --- 3. Differential Interactome Calculation (RelA Method) ---
cat("Calculating Differential Correlations...\n")

results <- data.frame(
    Protein_A = interactome$Protein_A,
    Protein_B = interactome$Protein_B,
    PCC_Control = NA,
    PCC_Case = NA,
    Delta_PCC = NA
)

# Loop through reference interactions
# (Note: In production, this step is parallelized)
for (i in 1:nrow(results)) {
    pA <- results$Protein_A[i]
    pB <- results$Protein_B[i]
    
    # Calculate PCC in Control
    pcc_ctrl <- calc_correlation(pA, pB, expr_data[, control_samples])
    
    # Calculate PCC in Case (Treatment)
    pcc_case <- calc_correlation(pA, pB, expr_data[, case_samples])
    
    if (!is.na(pcc_ctrl) && !is.na(pcc_case)) {
        results$PCC_Control[i] <- pcc_ctrl
        results$PCC_Case[i] <- pcc_case
        # The Metric: Absolute change in correlation
        results$Delta_PCC[i] <- abs(pcc_case - pcc_ctrl)
    }
}

# --- 4. Identify Significant Rewiring ---
# Threshold: Delta PCC distribution typically filters for > 0.5 or Top X percentile
delta_cutoff <- 0.5
rewired_edges <- results %>%
    filter(Delta_PCC > delta_cutoff)

cat(sprintf("Identified %d rewired interactions (Delta PCC > %.1f)\n", 
            nrow(rewired_edges), delta_cutoff))

# --- 5. Calculate Differential Degree (k_diff) ---
# k_diff = Sum of rewired edges for a protein
all_proteins <- unique(c(interactome$Protein_A, interactome$Protein_B))
degree_stats <- data.frame(Protein = all_proteins, k_diff = 0)

for (prot in all_proteins) {
    degree_stats$k_diff[degree_stats$Protein == prot] <- sum(
        (rewired_edges$Protein_A == prot) | (rewired_edges$Protein_B == prot)
    )
}

# Normalize (Optional but recommended by Gulfidan et al.)
# This script outputs raw k_diff for simplicity
degree_stats <- degree_stats %>% arrange(desc(k_diff))

write.csv(degree_stats, OUTPUT_FILE, row.names = FALSE)
cat("Done. Saved differential degree scores.\n")
