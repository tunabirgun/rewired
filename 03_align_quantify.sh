#!/usr/bin/env bash
# 03_align_quantify.sh — align cleaned reads (STAR) and produce gene-level counts (featureCounts).
# Tools: STAR 2.7.11b, samtools 1.22.1, featureCounts (subread) 2.1.1.
# Per-study host -> RefSeq assembly (manuscript Methods):
#   F. graminearum       GCF_000240135.3 (ASM24013v3)   studies 1,2,3
#   F. oxysporum         GCF_000149955.1 (ASM14995v2)   studies 4,5,7,8
#   F. pseudograminearum GCF_000303195.2 (FP7)          study 6
# Reference FASTA (*_genomic.fna) and annotation (*_genomic.gtf) must be placed in
#   data/ref/<asm>/   (download once from NCBI RefSeq for each assembly).
# For Study 3 (dual RNA-seq with T. hamatum) only reads mapping to F. graminearum are kept,
# which the F. graminearum-only reference index already enforces.
# Output: data/counts/<study>.txt  (featureCounts table; the DEG pipeline reads these)
set -euo pipefail

CLEAN="data/clean"
REFROOT="data/ref"
COUNTS="data/counts"
THREADS="${THREADS:-12}"
mkdir -p "$COUNTS"

declare -A ASM=(   [1]=fgr [2]=fgr [3]=fgr [4]=fox [5]=fox [6]=fpu [7]=fox [8]=fox )

build_index () {       # build a STAR index once per assembly
  local asm=$1; local dir="$REFROOT/$asm"
  [ -f "$dir/star_index/SAindex" ] && return
  mkdir -p "$dir/star_index"
  STAR --runMode genomeGenerate --genomeDir "$dir/star_index" \
       --genomeFastaFiles "$dir"/*_genomic.fna --sjdbGTFfile "$dir"/*_genomic.gtf \
       --genomeSAindexNbases 11 --runThreadN "$THREADS"
}

for studydir in "$CLEAN"/*/; do
  study=$(basename "$studydir"); asm="${ASM[$study]}"; dir="$REFROOT/$asm"
  build_index "$asm"
  mapped="${studydir}mapped"; mkdir -p "$mapped"
  for r1 in "${studydir}clean/"*_clean_R1.fq.gz; do
    [ -e "$r1" ] || continue
    srr=$(basename "$r1" _clean_R1.fq.gz); r2="${studydir}clean/${srr}_clean_R2.fq.gz"
    STAR --runMode alignReads --genomeDir "$dir/star_index" \
         --readFilesIn "$r1" "$r2" --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate --runThreadN "$THREADS" \
         --outFileNamePrefix "$mapped/${srr}_"
  done
  # gene-level counts (paired-end, exon, gene_id) — matches manuscript Methods
  featureCounts -a "$dir"/*_genomic.gtf -o "$COUNTS/${study}.txt" \
                -p -C -t exon -g gene_id -T 6 \
                "$mapped"/*_Aligned.sortedByCoord.out.bam
done

echo "Done. Count matrices in $COUNTS/<study>.txt"
