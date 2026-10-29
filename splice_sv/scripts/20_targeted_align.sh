#!/usr/bin/env bash
# Targeted alignment of the proband reads to hg38, keeping only BUB1B-anchored
# read pairs. Detects structural variants and mobile-element insertions that a
# small-variant VCF cannot (findings/BUB1B_STRUCTURAL.md, BUB1B_KMER_CN.md, which
# leave insertions and balanced events untested).
#
# Each lane is streamed from HuggingFace, aligned by bwa mem, and filtered on the
# fly to read pairs where either end maps within BUB1B +/- FLANK. The reads are
# never stored; only the tiny anchored BAM is kept. Mates that map elsewhere (a
# mobile element, another chromosome) are retained because their partner anchors
# in BUB1B -- that is the insertion / translocation signal.
#
# Full hg38 is required: a chr15-only reference would let paralog and mobile-element
# reads mismap onto BUB1B and manufacture false discordant signal.
set -euo pipefail

REF="${MVA_HG38:-/Volumes/IanSSD/mva-hackathon/ref/hg38/hg38.fna}"
OUT="${MVA_RESULTS:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/results}/align"
TOKEN="$(cat ~/.cache/huggingface/token)"
BASE="https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data/resolve/main"
REGION="chr15"; GSTART=40160984; GEND=40221137; FLANK=5000
LO=$((GSTART-FLANK)); HI=$((GEND+FLANK))
THREADS="${MVA_THREADS:-12}"
mkdir -p "$OUT"

[ -f "$REF.bwt" ] || { echo "MISSING bwa index for $REF (run: bwa index $REF)"; exit 1; }

# Download each lane pair to the SSD, align, then delete: robust to a mid-stream
# drop (a truncated bwa input would silently under-call), resumable, and 21GB per
# lane fits the free space. The fastqs are removed immediately after alignment so
# the 85GB release never accumulates.
dl() { curl -sL --retry 5 -H "Authorization: Bearer $TOKEN" -o "$2" "$BASE/$1"; }
TMP="$OUT/fastq_tmp"; mkdir -p "$TMP"

for L in 1 2 3 4; do
  R1="WGS_EX2312012_HGWCNDSX7_S16_L00${L}_R1_001.fastq.gz"
  R2="WGS_EX2312012_HGWCNDSX7_S16_L00${L}_R2_001.fastq.gz"
  bam="$OUT/bub1b_L00${L}.bam"
  [ -f "$bam" ] && { echo "skip L00${L} (exists)"; continue; }
  echo "[$(date +%H:%M:%S)] downloading lane $L ..."
  dl "$R1" "$TMP/$R1"; dl "$R2" "$TMP/$R2"
  echo "[$(date +%H:%M:%S)] aligning lane $L ..."
  # Keep a SAM record if the read OR its mate lands in BUB1B+/-FLANK. Header kept.
  bwa mem -t "$THREADS" -R "@RG\tID:L00${L}\tSM:PROBAND01\tPL:ILLUMINA\tLB:WGS" \
      "$REF" "$TMP/$R1" "$TMP/$R2" 2>"$OUT/bwa_L00${L}.log" \
    | awk -v c="$REGION" -v lo="$LO" -v hi="$HI" '
        /^@/ {print; next}
        {
          rn=$3; rp=$4; mn=$7; mp=$8;
          if (mn=="=") mn=rn;
          if ((rn==c && rp>=lo && rp<=hi) || (mn==c && mp>=lo && mp<=hi)) print
        }' \
    | samtools sort -@ 4 -o "$bam" -
  samtools index "$bam"
  rm -f "$TMP/$R1" "$TMP/$R2"
  echo "[$(date +%H:%M:%S)] lane $L done: $(samtools view -c "$bam") records"
done
rmdir "$TMP" 2>/dev/null || true

echo "merging lanes ..."
samtools merge -f "$OUT/bub1b_merged.bam" "$OUT"/bub1b_L00*.bam
samtools index "$OUT/bub1b_merged.bam"
echo "MERGED: $(samtools view -c "$OUT/bub1b_merged.bam") records in $OUT/bub1b_merged.bam"
echo "ALIGN_DONE" > "$OUT/align.status"
