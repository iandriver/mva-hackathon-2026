#!/usr/bin/env python3
"""Writes the Track 1 submission CSV from results/, and validates it against the spec.

The chr prefix is added here. The released VCF has none, and the submission format requires
chrN -- a mismatch scores zero without erroring anywhere.
"""
import csv, json, os, sys

REPO = os.path.join(os.path.dirname(__file__), "..")
OUT = os.environ.get("MVA_RESULTS", os.path.join(REPO, "results"))
DEST = os.path.join(REPO, "submissions", "idriver7_bub1b-compound-het.csv")
# The challenge accepts only this literal id. It is NOT the VCF sample name
# (WGS_EX2312012), which the scorer rejects with "Unknown proband_id".
PROBAND = "PROBAND01"
COLUMNS = ["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1",
           "chrom_2", "pos_2", "ref_2", "alt_2", "epcr", "finding_type", "notes"]

def chrom(c):
    return c if str(c).startswith("chr") else f"chr{c}"

def validate(rows):
    problems = []
    if not rows: problems.append("no rows")
    if len(rows) > 10: problems.append(f"{len(rows)} rows; the format accepts at most 10")
    seen_epcr = []
    for i, r in enumerate(rows, 1):
        if list(r.keys()) != COLUMNS:
            problems.append(f"row {i}: column order does not match the spec")
        if not str(r["chrom_1"]).startswith("chr"):
            problems.append(f"row {i}: chrom_1 '{r['chrom_1']}' lacks the chr prefix")
        if r["chrom_2"] and not str(r["chrom_2"]).startswith("chr"):
            problems.append(f"row {i}: chrom_2 '{r['chrom_2']}' lacks the chr prefix")
        e = float(r["epcr"])
        if not (0 < e <= 1): problems.append(f"row {i}: epcr {e} outside (0, 1]")
        seen_epcr.append(e)
        if r["finding_type"] not in ("primary", "secondary"):
            problems.append(f"row {i}: finding_type '{r['finding_type']}' invalid")
        for f in ("pos_1", "pos_2"):
            if r[f] and not str(r[f]).isdigit():
                problems.append(f"row {i}: {f} '{r[f]}' is not an integer")
    if seen_epcr != sorted(seen_epcr, reverse=True):
        problems.append("rows are not ordered by descending epcr")
    return problems

def main():
    locus = os.path.join(OUT, "02_bub1b_locus.json")
    if os.path.exists(locus):
        d = json.load(open(locus))
        coding = [v for v in d["variants"]
                  if v["in_gene_body"] and v["impact"] in ("HIGH", "MODERATE")]
        coding.sort(key=lambda v: v["pos"])
        print(f"BUB1B coding candidates from {locus}: "
              f"{', '.join(v['hgvsp'].split(':')[-1] for v in coding if v['hgvsp'])}",
              file=sys.stderr)
        if len(coding) != 2:
            print(f"WARNING: expected 2 coding candidates, found {len(coding)}. "
                  f"Review before submitting.", file=sys.stderr)
        v1, v2 = (coding + coding)[:2]
    else:
        print(f"NOTE: {locus} absent; using the recorded call. "
              f"Run scripts/02_bub1b_locus.py to regenerate from the VCF.", file=sys.stderr)
        v1 = {"chrom": "15", "pos": 40209701, "ref": "T", "alt": "G"}
        v2 = {"chrom": "15", "pos": 40220612, "ref": "T", "alt": "G"}

    rows = [
        dict(zip(COLUMNS, [
            PROBAND, chrom(v1["chrom"]), v1["pos"], v1["ref"], v1["alt"],
            chrom(v2["chrom"]), v2["pos"], v2["ref"], v2["alt"], 0.95, "primary",
            "BUB1B compound heterozygote. c.2210T>G p.Leu737Ter (ClinVar Pathogenic/Likely "
            "pathogenic, multiple submitters no conflicts; gnomAD v4 exome 115/1461846) in "
            "trans-inferred pairing with c.3006T>G p.Asn1002Lys (absent from ClinVar for this "
            "nucleotide change; gnomAD v4 exome 1/1461878). Biallelic BUB1B causes mosaic "
            "variegated aneuploidy type 1 (MVA1). Phase not demonstrated: variants are 10911 bp "
            "apart, beyond read-pair phasing, and no parental data is released."])),
    ]

    problems = validate(rows)
    if problems:
        for p in problems: print(f"INVALID: {p}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(os.path.dirname(DEST), exist_ok=True)
    with open(DEST, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader(); w.writerows(rows)
    print(f"wrote {DEST}")
    for r in rows:
        pair = f"  +  {r['chrom_2']}:{r['pos_2']} {r['ref_2']}>{r['alt_2']}" if r["chrom_2"] else ""
        print(f"  {r['finding_type']:<10} {r['chrom_1']}:{r['pos_1']} "
              f"{r['ref_1']}>{r['alt_1']}{pair}   epcr={r['epcr']}")
    print("\nvalidation: PASS")

if __name__ == "__main__":
    main()
