#!/usr/bin/env bash
# Reproduces the VCF-level facts: record count, reference build, FILTER summary.
# Establishes reference build, sample, contig naming and record count.
#
# The header is written to a temp FILE, not piped: under `set -o pipefail`, any downstream
# reader that exits early (head, grep -m1) raises SIGPIPE upstream and aborts the script.
source "$(dirname "$0")/env.sh"
need_cmd bcftools "brew install bcftools"
need_vcf
out="$MVA_RESULTS/01_vcf_facts.txt"
hdr=$(mktemp); trap 'rm -f "$hdr"' EXIT
bcftools view -h "$VCF" > "$hdr"
{
  echo "# VCF facts  ($(date -u +%Y-%m-%dT%H:%M:%SZ))"
  echo "# source: $VCF"
  echo
  echo "## reference build"
  grep -E '^##reference' "$hdr" || echo "(none declared)"
  echo
  echo "## sample"
  bcftools query -l "$VCF"
  echo
  echo "## contig naming"
  echo "contigs declared: $(grep -c '^##contig' "$hdr")"
  awk '/^##contig/ && ++n<=3' "$hdr"
  first_contig=$(awk -F'ID=|,' '/^##contig/{print $2; exit}' "$hdr")
  case "$first_contig" in
    chr*) echo "PREFIX: present ('$first_contig')" ;;
    *)    echo "PREFIX: ABSENT ('$first_contig') -- the Track 1 submission format requires" \
               "chrN, so the prefix must be added when writing the CSV. Silent scoring failure if missed." ;;
  esac
  echo
  echo "## caller"
  grep -oE 'ID=VariantFiltration|Version="[0-9.]+"' "$hdr" | sort -u
  echo
  echo "## FILTER values present"
  bcftools query -f '%FILTER\n' "$VCF" | sort | uniq -c | sort -rn
  echo "NOTE: LowQual calls are RETAINED in this file. A pipeline that filters to PASS by"
  echo "default discards candidates silently."
  echo
  echo "## record count"
  echo "records: $(bcftools view -H "$VCF" | wc -l | tr -d ' ')"
} > "$out"
cat "$out"
echo >&2 "wrote $out"
