#!/usr/bin/env python3
"""Reproduces findings/BUB1B_LOCUS.md.

Pulls every variant in the BUB1B gene span from the proband VCF, annotates via the
Ensembl VEP REST API, and attaches gnomAD v4 frequencies and ClinVar classifications.
Gene bounds are fetched from Ensembl rather than hardcoded.
"""
import json, os, subprocess, sys, urllib.request, urllib.parse, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import gnomad_variant, af_of, fmt_af

DATA = os.environ.get("MVA_DATA", "/Volumes/IanSSD/mva-hackathon/data")
VCF = os.path.join(DATA, "WGS_EX2312012_HGWCNDSX7.vcf.gz")
OUT = os.environ.get("MVA_RESULTS", os.path.join(os.path.dirname(__file__), "..", "results"))
os.makedirs(OUT, exist_ok=True)
FLANK = 5000

def post(url, payload, timeout=180):
    req = urllib.request.Request(url, data=json.dumps(payload).encode(),
                                 headers={"Content-Type": "application/json",
                                          "Accept": "application/json"})
    return json.load(urllib.request.urlopen(req, timeout=timeout))

def gene_bounds(symbol="BUB1B"):
    u = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{symbol}?content-type=application/json"
    g = json.load(urllib.request.urlopen(u, timeout=60))
    assert g["assembly_name"] == "GRCh38", f"expected GRCh38, got {g['assembly_name']}"
    return g["seq_region_name"], g["start"], g["end"], g["id"], g["canonical_transcript"]

def main():
    if not os.path.exists(VCF):
        sys.exit(f"MISSING: {VCF}\n  see scripts/env.sh for the fetch command")
    chrom, start, end, gid, tx = gene_bounds()
    print(f"BUB1B  {chrom}:{start:,}-{end:,}  {gid}  canonical {tx}", file=sys.stderr)

    region = f"{chrom}:{start-FLANK}-{end+FLANK}"
    q = subprocess.run(["bcftools", "query", "-r", region,
                        "-f", "%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t[%GT\t%AD\t%DP]\n", VCF],
                       capture_output=True, text=True, check=True).stdout.strip().splitlines()
    rows = [l.split("\t") for l in q if l]
    print(f"{len(rows)} variants in gene span +/- {FLANK}bp", file=sys.stderr)

    vep = post("https://rest.ensembl.org/vep/human/region"
               "?canonical=1&hgvs=1&numbers=1&check_existing=1&variant_class=1",
               {"variants": [f"{r[0]} {r[1]} . {r[2]} {r[3]} . . ." for r in rows]})
    by_pos = {}
    for v in vep:
        by_pos.setdefault(v["start"], []).append(v)

    out = []
    for r in rows:
        c, pos, ref, alt, qual, filt, gt, ad, dp = r
        pos_i = int(pos)
        cand = by_pos.get(pos_i) or by_pos.get(pos_i + 1) or []
        tc = []
        for v in cand:
            tc += [t for t in v.get("transcript_consequences", []) if t.get("gene_symbol") == "BUB1B"]
        can = [t for t in tc if t.get("canonical")] or tc
        t = can[0] if can else {}
        in_gene = start <= pos_i <= end
        status, payload = gnomad_variant(f"{c}-{pos}-{ref}-{alt}"); time.sleep(0.4)
        if status == "failed":
            print(f"  gnomAD query FAILED for {c}:{pos} {ref}>{alt} -- {payload}", file=sys.stderr)
        # payload is a dict only when status == "ok"
        g = payload if status == "ok" else None
        preds = {x["id"]: x["value"] for x in (g or {}).get("in_silico_predictors") or []}
        out.append({
            "chrom": c, "pos": pos_i, "ref": ref, "alt": alt, "qual": qual, "filter": filt,
            "gt": gt, "ad": ad, "dp": dp, "in_gene_body": in_gene,
            "consequence": (t.get("consequence_terms") or [None])[0],
            "impact": t.get("impact"), "exon": t.get("exon"), "intron": t.get("intron"),
            "hgvsc": t.get("hgvsc"), "hgvsp": t.get("hgvsp"),
            "gnomad_status": status,
            "gnomad_exome_af": ((g or {}).get("exome") or {}).get("af"),
            "gnomad_genome_af": ((g or {}).get("genome") or {}).get("af"),
            "rsids": (g or {}).get("rsids"), "predictors": preds,
        })

    path = os.path.join(OUT, "02_bub1b_locus.json")
    json.dump({"gene": {"id": gid, "chrom": chrom, "start": start, "end": end,
                        "canonical_transcript": tx}, "variants": out},
              open(path, "w"), indent=2)

    # Report exome and genome separately: deep-intronic variants are not exome-captured,
    # so an exome-only column shows common variants as "absent".
    print(f"\n{'pos':>10} {'change':<10} {'consequence':<22} {'impact':<9} "
          f"{'exome':>10} {'genome':>10}  hgvsp")
    for v in sorted(out, key=lambda x: x["pos"]):
        if not v["in_gene_body"]:
            continue
        def cell(af):
            if v["gnomad_status"] == "failed": return "FAILED"
            if v["gnomad_status"] == "absent": return "absent"
            return f"{af:.2e}" if af else "n/c"      # n/c = not covered by that callset
        print(f"{v['pos']:>10} {v['ref']+'>'+v['alt']:<10} {str(v['consequence']):<22} "
              f"{str(v['impact']):<9} {cell(v['gnomad_exome_af']):>10} "
              f"{cell(v['gnomad_genome_af']):>10}  {v['hgvsp'] or ''}")
    coding = [v for v in out if v["in_gene_body"] and v["impact"] in ("HIGH", "MODERATE")]
    nfail = sum(1 for v in out if v["gnomad_status"] == "failed")
    print(f"\ncoding candidates: {len(coding)}", file=sys.stderr)
    if nfail:
        print(f"WARNING: {nfail}/{len(out)} gnomAD queries failed. Those rows show "
              f"'QUERY FAILED', not 'absent' -- do not read them as evidence of rarity. "
              f"Re-run when the service recovers.", file=sys.stderr)
    print(f"wrote {path}", file=sys.stderr)

if __name__ == "__main__":
    main()
