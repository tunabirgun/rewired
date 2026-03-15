# Differential Interactome Meta-Analysis Pipeline

The pipeline identifies **rewired proteins** — proteins that undergo significant changes in their interaction landscape without proportional changes in gene expression — across eight independent RNA-seq studies comparing *Fusarium* response to fungal and bacterial biological control agents (BCAs).

---

## Methodology Reference

The **Differential Interactome** methodology (Q-value framework) used to quantify protein-protein interaction rewiring is adapted from:

> **Gulfidan, G., Turanli, B., Beklen, H. et al. (2020).**
> "Pan-cancer mapping of differential protein-protein interactions."
> *PLoS ONE* 15(3): e0228187.
> [https://doi.org/10.1371/journal.pone.0228187](https://doi.org/10.1371/journal.pone.0228187)

Protein-protein interaction data: **STRING v12.0** (combined score ≥ 400, taxid 229533 for *F. graminearum*).

---

## Repository Structure

```
scripts/
├── 00_setup.R                  # Shared configuration (sourced by all scripts)
├── 01_orthology_mapping.R      # Cross-species gene ID mapping via KEGG KO
├── 02_differential_interactome.R  # Core Q-value differential PPI algorithm
├── 03_rewired_analysis.R       # Rewired protein identification and classification
├── 04_deg_volcano.R            # Volcano plots with rewired protein overlay
├── 05_overlap_analysis.R       # Venn/UpSet cross-study overlap analysis
├── 06_enrichment.R             # KEGG + GO hypergeometric enrichment
└── 07_network_module.R         # Network construction, Louvain modules, hub genes
```

---

## Script Descriptions

### `00_setup.R` — Shared Configuration

Sourced at the top of every pipeline script. Defines all file paths, analysis thresholds, study metadata (8 studies × organism/BCA/condition), shared color palettes, and utility functions (`save_fig()`, `ensure_dir()`).

**No inputs required.** Sets `BASE_DIR` relative to the script location automatically.

---

### `01_orthology_mapping.R` — Cross-Species Orthology

- **Input:** KEGG KO assignment files (`fgr_ko.txt`, `fox_ko.txt`, `fpu_ko.txt`)
- **Process:** Maps *F. oxysporum* (FOXG) and *F. pseudograminearum* (FPSE) gene identifiers to *F. graminearum* (FGSG) IDs via shared KEGG Orthology groups. One-to-many mappings are resolved deterministically by selecting the lowest FGSG number.
- **Output:** `resources/fox_to_fgr_orthology.csv`, `resources/fpu_to_fgr_orthology.csv`

---

### `02_differential_interactome.R` — Core Q-value Algorithm

- **Input:** featureCounts expression matrices (`counts/1.txt` … `counts/8.txt`), STRING PPI network, orthology mappings from Step 01
- **Process:** For each STRING edge in each study:
  1. Binarizes gene expression per sample (1 if expression > mean − 1 SD, else 0)
  2. Counts co-expression in Control (*N*<sub>NC</sub>) and Treatment (*N*<sub>ND</sub>) conditions
  3. Computes Q-value = *N*<sub>ND</sub> / (*N*<sub>NC</sub> + *N*<sub>ND</sub>)
  4. Classifies edges: **Activated** (Q ≥ 0.90), **Repressed** (Q ≤ 0.10)
  5. Applies minimum frequency filter: max(*N*<sub>NC</sub>, *N*<sub>ND</sub>) ≥ Δ = 0.20
- **Output:** Per-study `protein_stats.csv` and `interaction_stats.csv` under `outputs/differential_interactome/Study_N/`

---

### `03_rewired_analysis.R` — Rewired Protein Identification

- **Input:** `protein_stats.csv` from Step 02, DEG lists (`degs/`)
- **Process:** A protein is classified as **rewired** if it meets both criteria:
  1. Differential degree (*N*<sub>total</sub>) in the **top 25%** across the study
  2. **Not a DEG**: |log₂FC| ≤ 1.0 or P<sub>adj</sub> ≥ 0.05
- **Output:** `rewired_proteins.csv` per study under `outputs/rewired_analysis/Study_N/`; aggregated `outputs/tables/rewired_summary.csv`

---

### `04_deg_volcano.R` — Volcano Plots *(optional)*

- **Input:** DEG lists, rewired protein lists from Step 03
- **Process:** Generates volcano plots (log₂FC vs −log₁₀ P<sub>adj</sub>) for each of the 8 studies. Rewired proteins are overlaid as highlighted points; top-ranked rewired proteins are labeled.
- **Output:** `outputs/figures/volcano_*.png/.svg`; panel figure `Fig1_volcano_panel_all_studies`

---

### `05_overlap_analysis.R` — Cross-Study Overlap *(optional)*

- **Input:** Rewired protein lists from Step 03, DEG lists
- **Process:** Computes pairwise and multi-way overlaps. Generates four-way Venn diagrams (fungal BCAs; bacterial BCAs separately) and an UpSet plot for all 8 studies combined.
- **Output:** `outputs/figures/FigS2_venn_fungal_rewired`, `FigS3_venn_bacterial_rewired`, `Fig3_upset_rewired_proteins`

---

### `06_enrichment.R` — Functional Enrichment

- **Input:** Rewired protein lists from Step 03; KEGG KO file (`fgr_ko.txt`); UniProt GO annotations (`fgr_uniprot_go.tsv`); STRING alias table for ID mapping
- **Process:** Hypergeometric enrichment test for KEGG pathways and GO terms (MF/CC/BP). P-value computed as P(X ≥ k) = `phyper(k − 1, m, n, N, lower.tail = FALSE)`. FDR correction by Benjamini–Hochberg. Pathway names fetched live from the KEGG REST API (rest.kegg.jp).
- **Output:** Enrichment tables and dot plots under `outputs/enrichment/` and `outputs/figures/Fig6_kegg_*`, `Fig7_go_*`

---

### `07_network_module.R` — Network Construction and Modules

- **Input:** Rewired protein lists from Step 03, STRING PPI network
- **Process:**
  1. Builds aggregated PPI sub-networks for fungal and bacterial rewired proteins
  2. Applies **Louvain community detection** (`igraph::cluster_louvain()`, seed = 42, default resolution = 1)
  3. Annotates modules with functional labels and hub genes (highest intra-module degree)
  4. Saves node/edge/hub CSVs for downstream network visualization
- **Output:** Network CSVs under `outputs/networks/`; module table `outputs/tables/module_summary.csv`; figures `Fig4a/b`, `Fig5_module_traits_heatmap`

---

## How to Run

Scripts are designed to be run sequentially. Each script sources `00_setup.R` automatically.

```bash
# Run full pipeline in order
Rscript 01_orthology_mapping.R
Rscript 02_differential_interactome.R
Rscript 03_rewired_analysis.R
Rscript 04_deg_volcano.R       # optional
Rscript 05_overlap_analysis.R  # optional
Rscript 06_enrichment.R
Rscript 07_network_module.R
```

> **Note:** Scripts 04 and 05 are optional visualization scripts and do not affect downstream analysis.

---

## Required Input Data

```
project_root/
├── counts/
│   └── 1.txt … 8.txt              # featureCounts output (tab-separated, with header)
├── degs/
│   └── {study_id}_{Up,Down}regulated_Genes.xlsx
├── metadata/
│   ├── study_metadata.csv          # Sample-level metadata (study_id, condition, replicate)
│   └── study_table.csv             # Study-level metadata (organism, BCA, type)
├── string_data/
│   ├── 229533.protein.links.v12.0.txt   # F. graminearum STRING PPI
│   ├── 229533.protein.aliases.v12.0.txt # F. graminearum protein aliases
│   ├── 426428.protein.aliases.v12.0.txt # F. oxysporum protein aliases
│   ├── fgr_ko.txt                  # F. graminearum KEGG KO assignments
│   ├── fox_ko.txt                  # F. oxysporum KEGG KO assignments
│   ├── fpu_ko.txt                  # F. pseudograminearum KEGG KO assignments
│   └── fgr_uniprot_go.tsv          # F. graminearum UniProt GO annotations
├── resources/                      # Auto-generated by 01_orthology_mapping.R
└── scripts/                        # This pipeline
```

---

## Required R Packages

```r
install.packages(c(
  "dplyr", "readr", "readxl", "openxlsx", "tidyr", "tibble",
  "ggplot2", "ggrepel", "ggraph", "igraph",
  "pheatmap", "VennDiagram", "UpSetR", "RColorBrewer",
  "gridExtra", "svglite", "httr", "stringr", "futile.logger"
))
```

Optional: `RCy3` for Cytoscape integration.

---

## Key Thresholds

| Parameter                           | Value                          | Source                 |
| ----------------------------------- | ------------------------------ | ---------------------- |
| Q-value — repression                | ≤ 0.10                         | Gulfidan et al. (2020) |
| Q-value — activation                | ≥ 0.90                         | Gulfidan et al. (2020) |
| Δ (minimum co-expression frequency) | ≥ 0.20                         | Gulfidan et al. (2020) |
| DEG log₂FC threshold                | > 1.0                          | Standard               |
| DEG adjusted p-value                | < 0.05                         | Standard               |
| Rewired protein percentile          | Top 25% of differential degree | This study             |
| STRING PPI score                    | ≥ 400 (medium confidence)      | STRING v12.0           |

---

## Output Structure

All outputs are written to `outputs/` relative to the project root:

```
outputs/
├── differential_interactome/Study_N/   # Q-value tables per study
├── differential_expression/            # DEG summaries
├── rewired_analysis/Study_N/           # Rewired protein lists per study
├── enrichment/                         # KEGG + GO enrichment results
├── networks/                           # Node/edge/hub CSVs for network figures
├── modules/                            # Louvain module assignments
├── figures/                            # All generated figures (PNG + SVG)
└── tables/                             # Summary tables
```
