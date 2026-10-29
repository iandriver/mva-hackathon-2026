#!/usr/bin/env bash
# Delly SV calls on the targeted BUB1B BAM. Delly detects deletions, insertions,
# duplications, inversions and translocations from split reads and discordant
# pairs. Run under Rosetta (bioconda osx-64 build) via the conda 'delly' env.
set -euo pipefail
REF="${MVA_HG38:-/Volumes/IanSSD/mva-hackathon/ref/hg38/hg38.fna}"
OUT="${MVA_RESULTS:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/results}/align"
DELLY="${MVA_DELLY:-/opt/homebrew/Caskroom/miniforge/base/envs/delly/bin/delly}"
BAM="$OUT/bub1b_merged.bam"
[ -f "$BAM" ] || { echo "MISSING $BAM (run 20_targeted_align.sh)"; exit 1; }

"$DELLY" call -g "$REF" -o "$OUT/delly_bub1b.bcf" "$BAM" 2> "$OUT/delly.log" || true
if command -v bcftools >/dev/null; then
  bcftools view "$OUT/delly_bub1b.bcf" 2>/dev/null > "$OUT/delly_bub1b.vcf" || true
fi
echo "=== Delly calls touching BUB1B (chr15:40,160,984-40,221,137) ==="
if [ -f "$OUT/delly_bub1b.vcf" ]; then
  awk -F'\t' '/^#/{next} $1=="chr15" && $2>=40150000 && $2<=40231000 {
     match($8,/SVTYPE=[^;]+/); st=substr($8,RSTART+7,RLENGTH-7);
     match($8,/END=[0-9]+/); en=substr($8,RSTART+4,RLENGTH-4);
     printf "  %s:%s  %s  END=%s  FILTER=%s\n",$1,$2,st,en,$7 }' "$OUT/delly_bub1b.vcf"
  n=$(awk -F'\t' '/^#/{next} $1=="chr15" && $2>=40150000 && $2<=40231000' "$OUT/delly_bub1b.vcf" | wc -l)
  echo "  ($n calls in/near BUB1B)"
else
  echo "  (no delly VCF produced; see $OUT/delly.log)"
fi
