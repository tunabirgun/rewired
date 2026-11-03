# Rewired

R and shell pipeline for differential-interactome analysis in *Fusarium*, from public sequencing reads through differential expression, orthology mapping, rewired-protein classification and downstream figures/tables.

The Q-value methodology is adapted from Gulfidan et al. (2020), “Pan-cancer mapping of differential protein–protein interactions,” *Scientific Reports* 10:3272, [doi:10.1038/s41598-020-60127-x](https://doi.org/10.1038/s41598-020-60127-x).

## Setup

Use Linux/HPC for the sequencing stages. Create the environment from `environment.yml`, activate `fusarium-biocontrol`, and run from the repository root. `config.R` defines paths, study metadata, thresholds and the seed. Check these against the intended study before execution; size thread counts to available resources.

## Inputs

Small inputs are bundled in `inputs/`; download the large STRING network and genome references separately. `config.R` reads the small inputs from `inputs/` and the large STRING files from `string_data/` at the project root.

**Bundled in `inputs/` (provided):**

- `study_table.csv`, `study_metadata.csv` — sample sheet / study metadata; the sample sheet drives the
  download loop and the control/treatment assignment (one row per run).
- `fgr_ko.txt`, `fox_ko.txt`, `fpu_ko.txt` — KEGG Orthology assignments. Re-downloadable from the KEGG
  REST API, e.g. `https://rest.kegg.jp/link/ko/fgr` (and `/fox`, `/fpu`).
- `fgr_uniprot_go.tsv` — *F. graminearum* GO annotations (UniProt, <https://www.uniprot.org>).

**Download once (not in git):**

- **STRING v12.0** → place in `string_data/` at the project root (the repository root). From <https://string-db.org>
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

Run each numbered stage separately and inspect its outputs before continuing. Invoke `.sh` scripts with `bash` and `.R` scripts with `Rscript`.

| Stages | Purpose |
|---|---|
| `01_download_sra.sh`–`03_align_quantify.sh` | Download, preprocess, align and count reads |
| `04_differential_expression.R` | Differential expression; `data/degs/` becomes the downstream DEG source |
| `05_orthology_mapping.R`–`07_rewired_analysis.R` | Orthology, differential interactomes and rewired-protein classification |
| `08_enrichment.R`–`09_network_modules.R` | Functional enrichment and network modules |
| `10_figures_volcano.R`–`13_supplementary_tables.R` | Figures and supplementary tables |

The full script filenames and executable parameters are in the repository. Preserve `data/`, `resources/`, `string_data/` and `outputs/` outside commits according to `.gitignore`. The scripts also support the documented `codes/` project layout detected by `config.R`.

Review study-specific assumptions, orthology tie-breaking, Q-value thresholds, enrichment backgrounds and comparison direction in the source before interpreting results. The environment manifests describe the dependency setup; this README does not establish reproduction across different package versions.
