#!/usr/bin/env bash
# 02_preprocess.sh — read QC, adapter/quality trimming, rRNA removal, paired-end repair.
# Tools: FastQC 0.12.1, MultiQC 1.30, fastp 1.0.1, SortMeRNA 4.3.4, BBMap 38.96 (repair.sh).
# Input : data/raw/<study>/<SRR>_{1,2}.fastq.gz   (from 01_download_sra.sh)
# Output: data/clean/<study>/clean/<SRR>_clean_R{1,2}.fq.gz  (+ QC reports)
set -euo pipefail

RAW="data/raw"
WORK="data/clean"
SMRNA_DB="${SORTMERNA_DB:-$HOME/sortmerna/database/smr_v4.3_default_db.fasta}"  # smr_v4.3_default_db.fasta
THREADS="${THREADS:-8}"

for studydir in "$RAW"/*/; do
  study=$(basename "$studydir")
  out="$WORK/$study"; mkdir -p "$out"/{pre_qc,fastp,rrna,clean,post_qc}
  for r1 in "$studydir"*_1.fastq.gz; do
    [ -e "$r1" ] || continue
    srr=$(basename "$r1" _1.fastq.gz); r2="${studydir}${srr}_2.fastq.gz"

    # 1) pre-trim QC
    fastqc -t "$THREADS" -o "$out/pre_qc" "$r1" "$r2"

    # 2) adapter/quality trimming (fastp, default parameters)
    fastp -i "$r1" -I "$r2" \
          -o "$out/fastp/${srr}_1.trim.fastq.gz" -O "$out/fastp/${srr}_2.trim.fastq.gz" \
          -h "$out/fastp/${srr}.html" -j "$out/fastp/${srr}.json" --thread 4

    # 3) rRNA removal (SortMeRNA, default database, keep non-rRNA reads)
    #    --paired_in treats a pair as rRNA only if both mates align; --fastx writes FASTQ.
    rm -rf "$out/rrna/${srr}_wd"
    sortmerna --ref "$SMRNA_DB" \
      --reads "$out/fastp/${srr}_1.trim.fastq.gz" \
      --reads "$out/fastp/${srr}_2.trim.fastq.gz" \
      --paired_in --fastx --threads "$THREADS" \
      --other "$out/rrna/${srr}_clean" --aligned "$out/rrna/${srr}_rrna" \
      --workdir "$out/rrna/${srr}_wd"

    # 4) re-pair the interleaved non-rRNA output into R1/R2 (BBMap repair.sh)
    #    (SortMeRNA emits one interleaved file when --out2 is not given; repair.sh splits it.)
    repair.sh in="$out/rrna/${srr}_clean.fq.gz" \
              out1="$out/clean/${srr}_clean_R1.fq.gz" \
              out2="$out/clean/${srr}_clean_R2.fq.gz" overwrite=t

    # 5) post-clean QC
    fastqc -t "$THREADS" -o "$out/post_qc" \
           "$out/clean/${srr}_clean_R1.fq.gz" "$out/clean/${srr}_clean_R2.fq.gz"
  done
  multiqc -f -o "$out" "$out/pre_qc" "$out/post_qc" "$out/fastp" || true
done

echo "Done. Cleaned paired reads in $WORK/<study>/clean/ ; MultiQC reports per study."
