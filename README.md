# Rewired

R and shell pipeline for differential-interactome analysis in *Fusarium*: public RNA-seq reads to differential expression, orthology mapping, rewired-protein classification, enrichment, network modules, figures and supplementary tables.

## Setup

Run the sequencing stages (01 to 03) on Linux or an HPC node. Create the environment, then run every script from the repository root.

```sh
conda env create -f environment.yml
conda activate fusarium-biocontrol
```

`config.R` is sourced by all R stages. It defines the directory layout, study metadata, analysis thresholds and random seed; check it against the intended study before running. Thread counts for the shell stages are set with the `THREADS` environment variable.

## Inputs

| Input | Source | Location |
|---|---|---|
| Sample sheet and study metadata: `study_table.csv`, `study_metadata.csv` | Bundled | `inputs/` |
| KEGG Orthology assignments: `fgr_ko.txt`, `fox_ko.txt`, `fpu_ko.txt` | Bundled (KEGG REST, `link/ko/<org>`) | `inputs/` |
| *F. graminearum* GO annotations: `fgr_uniprot_go.tsv` | Bundled (UniProt) | `inputs/` |
| STRING v12.0 network and aliases: `229533.protein.links.v12.0.txt`, `229533.protein.aliases.v12.0.txt` (*F. graminearum*), `426428.protein.aliases.v12.0.txt` (*F. oxysporum*) | [STRING](https://string-db.org), not bundled (links file ~100 MB) | `string_data/` |
| Genome FASTA (`*_genomic.fna`) and GTF (`*_genomic.gtf`): *F. graminearum* GCF_000240135.3, *F. oxysporum* GCF_000149955.1, *F. pseudograminearum* GCF_000303195.2 | [NCBI RefSeq](https://www.ncbi.nlm.nih.gov/datasets/), not bundled | `data/ref/fgr/`, `data/ref/fox/`, `data/ref/fpu/` |
| SortMeRNA database `smr_v4.3_default_db.fasta` | [SortMeRNA 4.3.4 release](https://github.com/sortmerna/sortmerna/releases/tag/v4.3.4), not bundled | path exported as `SORTMERNA_DB` |
| Reference DEG lists: `<study>_Upregulated_Genes.xlsx`, `<study>_Downregulated_Genes.xlsx` | Not bundled and not written by the scripts | `data/degs/` |

The sample sheet (`inputs/study_table.csv`, one row per run) drives the SRA download and the control and treatment assignment. `data/degs/` is the authoritative DEG base for stages 04, 07 and 13. `09_network_modules.R` reads `229533.protein.aliases.v12.0.txt` from `inputs/`, so place a copy there before running that stage.

## Run order

Run each stage separately and inspect its outputs before starting the next. Invoke `.sh` scripts with `bash` and `.R` scripts with `Rscript`.

| Stage | Purpose |
|---|---|
| `01_download_sra.sh` to `03_align_quantify.sh` | Download, preprocess, align and count reads (`data/counts/`) |
| `04_differential_expression.R` | DESeq2 with apeglm shrinkage; writes per-gene results to `outputs/differential_expression/` and checks the significant genes against `data/degs/` |
| `05_orthology_mapping.R` to `07_rewired_analysis.R` | Orthology to *F. graminearum*, differential interactomes, rewired-protein classification |
| `08_enrichment.R`, `09_network_modules.R` | Functional enrichment, network modules |
| `10_figures_volcano.R` to `13_supplementary_tables.R` | Figures and supplementary tables |

`data/`, `string_data/`, `resources/` and `outputs/` are generated or downloaded locally and excluded by `.gitignore`. `config.R` also accepts a layout in which the scripts sit in a `codes/` directory beside `data/`.

## Interpretation

Before interpreting results, review the study-specific assumptions, orthology tie-breaking, Q-value thresholds, enrichment backgrounds and comparison direction in the source. `environment.yml` defines the dependency set; reproduction under other package versions has not been verified.

## Methodology

The Q-value approach is adapted from Gulfidan et al. (2020), "Pan-cancer mapping of differential protein-protein interactions," *Scientific Reports* 10:3272, [doi:10.1038/s41598-020-60127-x](https://doi.org/10.1038/s41598-020-60127-x).
