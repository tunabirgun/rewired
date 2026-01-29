#!/usr/bin/env Rscript
# 02_rewired_classification.R
# distinguishes "Rewired Proteins" from "Transcriptional Targets" (DEGs).

library(dplyr)

# --- USER CONFIGURATION ---
DIFF_DEGREE_FILE <- "differential_interactome_scores.csv" # Output from Step 01
DEG_FILE <- "path/to/deseq2_results.csv"                  # Must have columns: gene, log2FoldChange, padj

# Thresholds
DEGREE_PERCENTILE <- 0.75  # Top 25% of differential connectivity are candidate hubs
LOGFC_CUTOFF <- 1.0        # Stability threshold
PADJ_CUTOFF <- 0.05        # Significance threshold

# --- 1. Load Data ---
k_stats <- read.csv(DIFF_DEGREE_FILE)
deg_data <- read.csv(DEG_FILE)

# --- 2. Define DEGs (Transcriptional Targets) ---
# Proteins that change in ABUNDANCE
deg_list <- deg_data %>%
    filter(padj < PADJ_CUTOFF & abs(log2FoldChange) > LOGFC_CUTOFF) %>%
    pull(gene)

cat(sprintf("Loaded %d Differentially Expressed Genes (DEGs)\n", length(deg_list)))

# --- 3. Define Rewired Proteins (Conformational/Functional Targets) ---
# Logic: High network change (Top 25% k_diff) BUT Stable expression (|logFC| < 1)

# Determine connectivity threshold
k_threshold <- quantile(k_stats$k_diff, DEGREE_PERCENTILE, na.rm = TRUE)

# Step A: Find high-connectivity variance proteins
high_variance_proteins <- k_stats %>%
    filter(k_diff >= k_threshold)

# Step B: REMOVE DEGs from this list
# We want "hidden" targets that do NOT change expression, only interactions
rewired_proteins <- high_variance_proteins %>%
    filter(!Protein %in% deg_list)

cat(sprintf("\nClassification Results:\n"))
cat(sprintf("  High Variance Network Nodes: %d\n", nrow(high_variance_proteins)))
cat(sprintf("  Filtered (Removed DEGs): -%d\n", nrow(high_variance_proteins) - nrow(rewired_proteins)))
cat(sprintf("  Final Rewired Proteins: %d\n", nrow(rewired_proteins)))

# --- 4. Export ---
write.csv(rewired_proteins, "candidates_rewired_proteins.csv", row.names = FALSE)
cat("Saved candidate rewired list.\n")
