#!/usr/bin/env bash
# Unbiased ClinVar P/LP intersection over the whole VCF, with ACMG SF v3.2 triage.
# Genotype-bearing output stays gitignored and is never committed.
# Unbiased: intersects the whole proband VCF against the ENTIRE ClinVar P/LP set,
# then judges reportability against ACMG SF v3.2 afterwards.
source "$(dirname "$0")/env.sh"
need_cmd bcftools "brew install bcftools"
need_vcf
mkdir -p "$MVA_REF"
CV="$MVA_REF/clinvar.vcf.gz"
if [ ! -f "$CV" ]; then
  echo "fetching ClinVar GRCh38 VCF (~190 MB)..." >&2
  curl -sL -o "$CV"     https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
  curl -sL -o "$CV.tbi" https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
fi
echo "ClinVar release: $(bcftools view -h "$CV" | grep -m1 '^##fileDate')" >&2

PLP="$MVA_REF/clinvar_plp.vcf.gz"
if [ ! -f "$PLP" ]; then
  bcftools view -i 'INFO/CLNSIG~"Pathogenic" || INFO/CLNSIG~"Likely_pathogenic"' "$CV" \
    | bcftools view -e 'INFO/CLNSIG~"Conflicting"' -Oz -o "$PLP"
  bcftools index -t -f "$PLP"
fi
echo "ClinVar P/LP records: $(bcftools view -H "$PLP" | wc -l | tr -d ' ')" >&2

HITS="$MVA_RESULTS/03_proband_clinvar_plp.vcf.gz"
bcftools isec -n=2 -w1 -Oz -o "$HITS" "$VCF" "$PLP"
bcftools index -t -f "$HITS"

out="$MVA_RESULTS/03_secondary_findings.tsv"
{
  printf "chrom\tpos\tref\talt\tfilter\tgt\tad\tdp\tmq\tmqrs\tqd\tgene\tclnsig\tclnrevstat\tmc\tclndn\n"
  bcftools annotate -a "$PLP" \
    -c INFO/CLNSIG,INFO/CLNDN,INFO/GENEINFO,INFO/CLNREVSTAT,INFO/MC "$HITS" \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%GT\t%AD\t%DP]\t%MQ\t%MQRankSum\t%QD\t%GENEINFO\t%CLNSIG\t%CLNREVSTAT\t%MC\t%CLNDN\n'
} > "$out"
echo >&2 "exact REF/ALT matches: $(( $(wc -l < "$out") - 1 ))"
echo >&2 "wrote $out"
echo >&2 ""
echo >&2 "Next: scripts/04_acmg_triage.py decides reportability against ACMG SF v3.2."
