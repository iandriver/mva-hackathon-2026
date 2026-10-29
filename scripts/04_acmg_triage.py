#!/usr/bin/env python3
"""Decides which ClinVar P/LP hits are reportable ACMG SF v3.2 secondary findings.

The gene list is fetched from NCBI rather than recalled, then every hit is scored on
population frequency, allele balance and mapping-quality bias so that a ClinVar label
alone cannot carry a finding into the report.
"""
import csv, json, os, re, sys, time, urllib.request
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import gnomad_variant, af_of

OUT = os.environ.get("MVA_RESULTS", os.path.join(os.path.dirname(__file__), "..", "results"))
HITS = os.path.join(OUT, "03_secondary_findings.tsv")

# Two distinct things are being separated here, and conflating them is a mistake:
#   - a CALL that is not trustworthy (bad allele balance, mapping bias)  -> REJECT_CALL
#   - a trustworthy call that is simply too common to be a penetrant
#     dominant cause, i.e. a carrier or benign variant                   -> TOO_COMMON
# Only the first is an artifact. A recessive carrier allele is a real finding that is
# merely not reportable under ACMG SF.
MAX_AF_DOMINANT = 1e-3          # above this, no penetrant dominant role
MIN_VAF, MAX_VAF = 0.30, 0.70   # germline heterozygote
MIN_MQ = 50.0                   # mapping quality
MIN_MQRS = -3.0                 # ref/alt mapping-quality bias

def acmg_sf_genes():
    """ACMG SF v3.2 gene symbols, scraped from the NCBI ClinVar page."""
    html = urllib.request.urlopen(
        "https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/", timeout=90).read().decode("utf8", "replace")
    genes = set(re.findall(r'/clinvar/\?term=([A-Z0-9]+)\[gene\]', html))
    genes |= set(re.findall(r'gene/\d+"[^>]*>([A-Z][A-Z0-9]{1,9})<', html))
    if len(genes) < 60:
        sys.exit(f"ACMG SF scrape returned only {len(genes)} genes; refusing to triage "
                 f"against an incomplete list. Check https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/")
    return genes

def main():
    if not os.path.exists(HITS):
        sys.exit(f"MISSING: {HITS}\n  run: scripts/03_secondary_findings.sh")
    sf = acmg_sf_genes()
    print(f"ACMG SF v3.2: {len(sf)} genes", file=sys.stderr)

    verdicts = []
    for r in csv.DictReader(open(HITS), delimiter="\t"):
        gene = (r["gene"] or "").split(":")[0]
        ad = r["ad"].split(",")
        vaf = int(ad[1]) / max(1, int(ad[0]) + int(ad[1])) if len(ad) > 1 else None
        status, payload = gnomad_variant(f"{r['chrom']}-{r['pos']}-{r['ref']}-{r['alt']}")
        time.sleep(0.4)
        af = af_of(status, payload)

        def num(field):
            try: return float(r.get(field, "") or "nan")
            except ValueError: return float("nan")
        mq, mqrs = num("mq"), num("mqrs")

        call_problems, freq_problems = [], []
        if vaf is not None and not (MIN_VAF <= vaf <= MAX_VAF):
            call_problems.append(f"VAF {vaf:.3f} - not a germline heterozygote")
        if mq == mq and mq < MIN_MQ:
            call_problems.append(f"MQ {mq:.2f} - low mapping quality")
        if mqrs == mqrs and mqrs < MIN_MQRS:
            call_problems.append(f"MQRankSum {mqrs:.3f} - ref/alt mapping-quality bias")
        if af is not None and af > MAX_AF_DOMINANT:
            freq_problems.append(f"gnomAD AF {af:.3g} - too common for a penetrant dominant allele")
        reasons = call_problems + freq_problems
        if status == "failed":
            # Frequency is one of the two rejection channels. Without it the triage is
            # incomplete, so the row is held rather than passed or rejected.
            verdict = "INCONCLUSIVE"
            reasons.append(f"gnomAD query failed ({payload}); frequency filter not applied")
        elif call_problems:
            verdict = "REJECT_CALL"          # the call itself is not trustworthy
        elif freq_problems:
            verdict = "TOO_COMMON"           # real call, but no penetrant dominant role
        else:
            verdict = "ACMG_SF_REPORTABLE" if gene in sf else "NOT_SF_GENE"
        verdicts.append({"locus": f"{r['chrom']}:{r['pos']} {r['ref']}>{r['alt']}",
                         "gene": gene, "clnsig": r["clnsig"], "consequence": r["mc"],
                         "vaf": round(vaf, 3) if vaf else None, "dp": r["dp"],
                         "gnomad_af": af, "gnomad_status": status, "in_acmg_sf": gene in sf,
                         "mq": None if mq != mq else mq, "mqrs": None if mqrs != mqrs else mqrs,
                         "verdict": verdict, "reasons": reasons})

    json.dump(verdicts, open(os.path.join(OUT, "04_acmg_triage.json"), "w"), indent=2)
    print(f"\n{'locus':<26}{'gene':<10}{'VAF':>7}{'gnomAD':>11}  verdict")
    for v in verdicts:
        af = v["gnomad_af"]
        shown = ("FAILED" if v["gnomad_status"] == "failed"
                 else f"{af:.2e}" if af else "absent")
        print(f"{v['locus']:<26}{v['gene']:<10}{(v['vaf'] or 0):>7.3f}"
              f"{shown:>13}  {v['verdict']}")
        for why in v["reasons"]:
            print(f"{'':<54}  - {why}")
    n_sf = sum(1 for v in verdicts if v["verdict"] == "ACMG_SF_REPORTABLE")
    n_bad = sum(1 for v in verdicts if v["verdict"] == "INCONCLUSIVE")
    print(f"\nACMG SF v3.2 reportable secondary findings: {n_sf}", file=sys.stderr)
    if n_bad:
        print(f"WARNING: {n_bad} row(s) INCONCLUSIVE because gnomAD was unreachable. "
              f"The triage is incomplete; re-run before relying on it.", file=sys.stderr)
        sys.exit(2)

if __name__ == "__main__":
    main()
