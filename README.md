# Systems-Level Host Defense Analysis Pipeline

This repository contains the core analytical code used to identify "Rewired" proteins and functional defense modules in *Fusarium* species. The scripts provided here are simplified, non-looping versions intended to demonstrate the underlying methodology logic. They are designed to be adaptable to any user-provided RNA-seq and Interactome dataset.

## Methodology Reference

The **Differential Interactome** methodology used to calculate protein-protein interaction rewiring is based on the "RelA" approach described in:

> **Gulfidan, G., Turanli, B., Beklen, H. et al. (2020).**
> "Pan-cancer mapping of differential protein-protein interactions".
> *Scientific Reports* 10, 3272.
> [https://doi.org/10.1038/s41598-020-60127-x](https://doi.org/10.1038/s41598-020-60127-x)

## Repository Structure

The workflow is divided into three logical phases:

1.  **`01_differential_interactome.R`**
    *   **Input:** Normalized Gene Expression Matrix (TPM/FPKM) and a Reference Interactome (PPI).
    *   **Process:** Calculates Pearson Correlation Coefficients (PCC) for all interactions in Control vs. Treatment conditions. Defines "Rewiring" as the absolute difference in correlation ($\Delta PCC$).
    *   **Output:** A list of proteins ranked by their "Differential Degree" ($k_{diff}$ — the sum of rewired edges).

2.  **`02_rewired_classification.R`**
    *   **Input:** Output from Step 01 and a list of Differentially Expressed Genes (DEGs) from standard tools (e.g., DESeq2).
    *   **Classification:**
        *   **Rewired Proteins:** Top 25% highest Differential Degree AND Stable Expression ($|logFC| < 1$).
        *   **Transcriptional Targets:** Significant DEGs ($P_{adj} < 0.05$).
    *   **Goal:** To distinguish proteins undergoing *conformational/functional* changes (Rewiring) from those undergoing *abundance* changes (Expression).

3.  **`03_module_discovery.R`**
    *   **Input:** The sub-network of identified Rewired Proteins.
    *   **Process:** Applies **Louvain Modularity** detection to identify dense functional clusters (modules) within the rewired network.
    *   **Output:** List of functional modules representative of the system's coordinated defense response.

## Usage

These scripts use placeholders (e.g., `"path/to/expression_matrix.csv"`) and should be updated with your specific file paths.

```R
# Example: Run the differential interactome calculation
source("01_differential_interactome.R")
```

## Requirements
*   R (v4.0+)
*   `igraph`
*   `dplyr`
*   `WGCNA` (optional, for fast correlation) or base R
