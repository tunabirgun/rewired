End-to-end pipeline, from raw SRA reads to the differential-interactome / rewired-protein
results and every manuscript figure and supplementary table. Scripts run in numeric order.
The backbone (01–07) takes SRA → rewired proteins; the downstream scripts (08–13) produce
enrichment, modules, figures, and supplementary tables. Figure and supplementary generation are
kept in their own scripts, separate from the analytical backbone.

The differential-interactome (Q-value) methodology is adapted from **Gulfidan et al. (2020),
"Pan-cancer mapping of differential protein–protein interactions," *Sci Rep* 10:3272,
doi:10.1038/s41598-020-60127-x.**

## Layout

```text
codes/
  config.R                      shared configuration (paths, thresholds, study metadata, seed=42)
  environment.yml               conda environment (CLI tools + R/Bioconductor packages)
  inputs/                       bundled small inputs (sample sheet, metadata, KEGG KO, UniProt GO)
  .gitignore                    excludes the large STRING files from commits

  # ---- backbone: SRA -> rewired proteins ----
  01_download_sra.sh            prefetch + fasterq-dump        -> data/raw/<study>/*.fastq.gz
  02_preprocess.sh              FastQC, fastp, SortMeRNA, BBMap, MultiQC -> data/clean/<study>/clean/
  03_align_quantify.sh          STAR + featureCounts           -> data/counts/<study>.txt
  04_differential_expression.R  DESeq2 (apeglm shrinkage)       -> data/degs/, outputs/differential_expression/
  05_orthology_mapping.R        KO-based FOXG/FPSE -> FGSG       -> resources/
  06_differential_interactome.R Q-value differential interactome-> outputs/differential_interactome/
  07_rewired_analysis.R         rewired-protein classification  -> outputs/rewired_analysis/

  # ---- downstream: enrichment, modules, figures, supplementary ----
  08_enrichment.R               KEGG/GO hypergeometric enrichment (DEGs + rewired)
  09_network_modules.R          aggregated networks + Louvain modules
  10_figures_volcano.R          volcano panel with rewired overlay (Fig 4)
  11_figures_overlaps.R         DEG/rewired Venn + UpSet + classification (Figs 2, 3, 5, 6)
  12_figures_marker.R           marker-gene heatmap (Fig 8)
  13_supplementary_tables.R     Supplementary Tables S1-S9
```

## Inputs

The small inputs are bundled in `codes/inputs/`; the large STRING network and the genome references
are not committed (too large for git) and are downloaded once. `config.R` reads the small inputs from
`codes/inputs/` and the large STRING files from `../string_data/`.

**Bundled in `codes/inputs/` (provided):**

- `study_table.csv`, `study_metadata.csv` — sample sheet / study metadata; the sample sheet drives the
  download loop and the control/treatment assignment (one row per run).
- `fgr_ko.txt`, `fox_ko.txt`, `fpu_ko.txt` — KEGG Orthology assignments. Re-downloadable from the KEGG
  REST API, e.g. `https://rest.kegg.jp/link/ko/fgr` (and `/fox`, `/fpu`).
- `fgr_uniprot_go.tsv` — *F. graminearum* GO annotations (UniProt, <https://www.uniprot.org>).

**Download once (not in git):**

- **STRING v12.0** → place in `string_data/` at the project root. From <https://string-db.org>
  (version 12.0): `229533.protein.links.v12.0.txt` + `229533.protein.aliases.v12.0.txt`
  (*F. graminearum*, taxid 229533) and `426428.protein.aliases.v12.0.txt` (*F. oxysporum*, taxid
  426428). The `.protein.links` file is ~100 MB, which is why it is kept out of the repository.
- **Genome references** → `data/ref/<asm>/` (FASTA `*_genomic.fna` + GTF `*_genomic.gtf`). NCBI RefSeq:
  `fgr` = GCF_000240135.3, `fox` = GCF_000149955.1, `fpu` = GCF_000303195.2
  (<https://www.ncbi.nlm.nih.gov/datasets/>).
- **SortMeRNA default DB** `smr_v4.3_default_db.fasta` — from the SortMeRNA 4.3.4 release
  (<https://github.com/sortmerna/sortmerna>); export `$SORTMERNA_DB`.

`data/degs/` (the DESeq2 output of step 04) is the authoritative DEG base for all downstream steps.

## Run order

```bash
# backbone (run the heavy CLI steps on Linux/HPC)
bash codes/01_download_sra.sh
THREADS=8  SORTMERNA_DB=~/sortmerna/database/smr_v4.3_default_db.fasta  bash codes/02_preprocess.sh
THREADS=12 bash codes/03_align_quantify.sh
Rscript codes/04_differential_expression.R
Rscript codes/05_orthology_mapping.R
Rscript codes/06_differential_interactome.R     # heavy: STRING co-expression over 8 studies
Rscript codes/07_rewired_analysis.R
# downstream
Rscript codes/08_enrichment.R
Rscript codes/09_network_modules.R
Rscript codes/10_figures_volcano.R
Rscript codes/11_figures_overlaps.R
Rscript codes/12_figures_marker.R
Rscript codes/13_supplementary_tables.R
```

## Key parameters to check (all centralised in `config.R`)

| Stage | Parameter | Value | Where |
|---|---|---|---|
| featureCounts | mode | `-p -C -t exon -g gene_id` | `03_align_quantify.sh` |
| DESeq2 | low-count filter | total counts ≥ 10 | `04` |
| DESeq2 | reference level | `Control` (vs `Sample`) | `04` |
| DESeq2 | LFC shrinkage | `apeglm` | `04` |
| DEG call | significance | `padj < 0.05` and `abs(log2FC) > 1` | `config.R` (`PADJ_THRESHOLD`, `LOG2FC_THRESHOLD`) |
| Orthology | one-to-many tie-break | lowest FGSG numeric ID | `05` |
| Interactome | STRING confidence | combined score ≥ 400 | `06` |
| Interactome | binarisation | +1/0/−1 at per-sample mean ± 1 SD | `06` |
| Interactome | Q cut-offs | activated Q ≥ 0.90, repressed Q ≤ 0.10 | `config.R` (`QCRIT1`, `QCRIT2`) |
| Interactome | min frequency | `max(N_ctrl, N_trt) ≥ 0.20` | `config.R` (`DELTA`) |
| Rewired | differential-degree cut | top 25th percentile | `config.R` (`DEGREE_PERCENTILE`) |
| Rewired | exclusion | not a DEG | `07` |
| Enrichment | test / FDR | hypergeometric (`phyper`), Benjamini–Hochberg < 0.05 | `08` |
| Modules | Louvain | seed = 42, default resolution; keep size ≥ 3 + sig. term | `config.R`, `09` |

Other things to watch: STAR needs ≈ 32 GB RAM for index/alignment; step 06 is the slowest (STRING
co-expression across all studies); set `THREADS` per machine. **Study 3** is dual RNA-seq — only
*F. graminearum*-mapped reads are kept, which the F. graminearum-only STAR index enforces. The
global seed (42) in `config.R` fixes Louvain clustering and layout. To re-run only the downstream
analysis without touching the (authoritative) DEGs, start from step 05.

## Methods summary (matches the manuscript)

- **Trimming** fastp (default). **rRNA** SortMeRNA default DB. **Alignment** STAR to the per-host
  RefSeq genome. **Counts** featureCounts (paired, exon, gene_id).
- **DEGs** DESeq2 with apeglm LFC shrinkage; ref = Control; genes with < 10 total counts removed;
  `padj < 0.05` and `|log2FC| > 1`.
- **Orthology** FOXG / FPSE mapped to FGSG by shared KEGG Orthology; lowest FGSG ID on ties.
- **Differential interactome** STRING v12.0 (≥ 400). Per sample, expression discretised to +1/0/−1
  at mean ± 1 SD; co-expression counted over the nine pairwise states; Q = N_trt / (N_ctrl + N_trt);
  activated Q ≥ 0.90, repressed Q ≤ 0.10, `max(N_ctrl, N_trt) ≥ 0.20`.
- **Rewired proteins** top 25th percentile of differential degree and not a DEG.
- **Enrichment** hypergeometric (`phyper`), BH FDR < 0.05; KEGG names via the KEGG REST API; GO from
  UniProt. **Modules** Louvain (`igraph`, seed = 42); modules with ≥ 3 members and a significant
  KEGG/GO term retained.

## Tool / package versions

Validated environment (the versions this pipeline was run and checked with):

- **CLI:** sra-tools (bioconda current), FastQC 0.12.1, MultiQC 1.30, fastp 1.0.1, SortMeRNA 4.3.4,
  BBMap 38.96, STAR 2.7.11b, samtools 1.22.1, Subread/featureCounts 2.1.1.
- **R 4.5.3:** DESeq2 1.50.2, apeglm 1.32.0, dplyr 1.2.1, readr 2.2.0, readxl 1.4.5, openxlsx 4.2.8.1,
  tidyr 1.3.2, tibble 3.3.1, ggplot2 4.0.2, ggrepel 0.9.8, ggraph 2.2.2, igraph 2.2.3, pheatmap 1.0.13,
  VennDiagram 1.8.2, UpSetR 1.4.0, RColorBrewer 1.1.3, svglite 2.2.2, httr 1.4.8, tictoc 1.2.1,
  stringr 1.6.0, gridExtra 2.3.

The manuscript Methods cite the original analysis versions (e.g. R 4.5.2, DESeq2 1.44.0, igraph 2.0.3,
httr 1.4.7); the results reproduce across these minor version differences.

## Figure / table provenance

- Fig 4 (volcano) — `10_figures_volcano.R`
- Fig 2 (DEG Venn), Fig 3 (DEG UpSet), Fig 5 (rewired UpSet), Fig 6 (rewired classification) — `11_figures_overlaps.R`
- Fig 8 (marker heatmap) — `12_figures_marker.R`
- Fig 1 (DEG KEGG) — `08_enrichment.R`; Fig 7 (rewired PPI networks) — `09_network_modules.R`
- Tables 2–10 — `05–09`; Supplementary S1–S9 — `13_supplementary_tables.R`
