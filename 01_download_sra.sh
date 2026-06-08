#!/usr/bin/env bash
# 01_download_sra.sh — download raw reads from NCBI SRA, convert to gzipped paired FASTQ.
# Tool: sra-tools (prefetch, fasterq-dump).
# Driven by the sample sheet metadata/study_table.csv
#   columns: Study No., Study, Organism, Biocontrol Agent, Experimental Condition,
#            SRA Accession, Control/Sample, BCA Type
# Output: data/raw/<study>/<SRR>_1.fastq.gz , <SRR>_2.fastq.gz
set -euo pipefail

SAMPLES="${1:-metadata/study_table.csv}"
RAW="data/raw"
THREADS="${THREADS:-8}"
mkdir -p "$RAW"

# (study, SRR) pairs from the sample sheet (col 1 = study, col 6 = SRA accession)
tail -n +2 "$SAMPLES" | awk -F',' 'NF>=6 {gsub(/[ \r]/,"",$1); gsub(/[ \r]/,"",$6); print $1"\t"$6}' \
| while IFS=$'\t' read -r study srr; do
    [ -z "${srr:-}" ] && continue
    out="$RAW/$study"; mkdir -p "$out"
    if [ -s "$out/${srr}_1.fastq.gz" ]; then
        echo "[$study] $srr already present, skipping"; continue
    fi
    echo "[$study] prefetch $srr"
    prefetch "$srr" -O "$out"
    echo "[$study] fasterq-dump $srr"
    fasterq-dump "$out/$srr/$srr.sra" -O "$out" --split-files -e "$THREADS"
    gzip -f "$out/${srr}_1.fastq" "$out/${srr}_2.fastq"
done

echo "Done. Raw paired FASTQ written to $RAW/<study>/"
