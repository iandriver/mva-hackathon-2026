#!/usr/bin/env bash
# Shared locations. Override by exporting before running any script.
# Gated data never lives in this repository (see DATA_HANDLING.md).
set -euo pipefail
MVA_DATA="${MVA_DATA:-/Volumes/IanSSD/mva-hackathon/data}"
MVA_REF="${MVA_REF:-/Volumes/IanSSD/mva-hackathon/ref}"
MVA_REPO="${MVA_REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
MVA_RESULTS="${MVA_RESULTS:-$MVA_REPO/results}"
VCF="$MVA_DATA/WGS_EX2312012_HGWCNDSX7.vcf.gz"
mkdir -p "$MVA_RESULTS"

# A missing input names the command that fetches it, rather than producing a traceback.
need_file() {
  [ -f "$1" ] && return 0
  echo "MISSING: $1" >&2
  echo "  fetch: $2" >&2
  exit 1
}
need_cmd() {
  command -v "$1" >/dev/null 2>&1 && return 0
  echo "MISSING COMMAND: $1" >&2
  echo "  install: $2" >&2
  exit 1
}
need_vcf() {
  need_file "$VCF" "hf download SageBio/mva-hackathon-2026-data --repo-type dataset \\
      --local-dir \"\$MVA_DATA/..\" WGS_EX2312012_HGWCNDSX7.vcf.gz WGS_EX2312012_HGWCNDSX7.vcf.gz.tbi
      (requires an accepted Data Transfer Agreement on the dataset page)"
}
