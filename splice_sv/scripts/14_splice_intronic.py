#!/usr/bin/env python3
"""Scores the proband's BUB1B variants for splice disruption with SpliceAI.

Closes caveat 4 of findings/BUB1B_LOCUS.md: the twelve intronic variants were
never evaluated for splice impact, and MECHANISM.md 2 records that published
MVA1 second alleles are repeatedly non-coding and splice-disrupting.

SpliceAI runs entirely locally. The weights ship in the package and no network
call is made, so no proband data leaves this machine. Ensembl is queried only
for the public gene model.

Scores are interpreted against calibration/bub1b_spliceai_calibration.csv,
a 582-variant panel of public gnomAD and literature variants in this same gene,
built and validated in ../alphagenome_control (public data). On that panel
SpliceAI recovered 3/3 published pathogenic MVA1 splice alleles at 6.9%
background using its own recommended DS >= 0.5 cutoff.

Outputs carrying genotypes go to MVA_RESULTS, which is gitignored.

SpliceAI: Jaganathan et al. 2019, Cell, doi:10.1016/j.cell.2018.12.015.
Source under PolyForm Strict 1.0.0, models under CC BY-NC 4.0, non-commercial
use. Run unmodified and not redistributed.
"""
import json, os, re, subprocess, sys, urllib.request

DATA = os.environ.get("MVA_DATA", "/Volumes/IanSSD/mva-hackathon/data")
REF = os.environ.get("MVA_REF", "/Volumes/IanSSD/mva-hackathon/ref")
VCF = os.path.join(DATA, "WGS_EX2312012_HGWCNDSX7.vcf.gz")
FASTA = os.path.join(REF, "chr15.fa.bgz")
OUT = os.environ.get("MVA_RESULTS",
                     os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                  "..", "..", "results"))  # repo-root results/, gitignored
CALIB = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     "calibration", "bub1b_spliceai_calibration.csv")
SPLICEAI = os.environ.get(
    "MVA_SPLICEAI",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                 "alphagenome_control", ".venv-spliceai", "bin", "spliceai"))
DIST = int(os.environ.get("MVA_SPLICEAI_DIST", "50"))
os.makedirs(OUT, exist_ok=True)


def need(path, how):
    if not os.path.exists(path):
        sys.exit(f"MISSING: {path}\n  {how}")


def gene_model(symbol="BUB1B"):
    """Canonical transcript exon bounds, from Ensembl. Never hardcoded."""
    u = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{symbol}?content-type=application/json"
    g = json.load(urllib.request.urlopen(u, timeout=60))
    assert g["assembly_name"] == "GRCh38", f"expected GRCh38, got {g['assembly_name']}"
    tx = g["canonical_transcript"].split(".")[0]
    u = f"https://rest.ensembl.org/lookup/id/{tx}?expand=1;content-type=application/json"
    t = json.load(urllib.request.urlopen(u, timeout=60))
    exons = sorted(({"start": e["start"], "end": e["end"]} for e in t["Exon"]),
                   key=lambda e: e["start"])
    return g["seq_region_name"], g["start"], g["end"], t["id"], exons


def offset(pos, exons):
    """0 if exonic, else signed distance into the intron: +n past a donor,
    -n before an acceptor, whichever is nearer."""
    for e in exons:
        if e["start"] <= pos <= e["end"]:
            return ("exonic", 0)
    best = None
    for e in exons:
        if pos > e["end"] and (best is None or pos - e["end"] < abs(best[1])):
            best = ("donor", pos - e["end"])
        if pos < e["start"] and (best is None or e["start"] - pos < abs(best[1])):
            best = ("acceptor", -(e["start"] - pos))
    return best


def main():
    need(VCF, "see scripts/env.sh for the fetch command")
    need(FASTA, "curl -sL https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr15.fa.gz"
                f" | gunzip -c | bgzip > {FASTA} && samtools faidx {FASTA}")
    need(SPLICEAI, "uv venv --python 3.11 ENV && uv pip install 'tensorflow==2.15.*' "
                   "'keras==2.15.*' spliceai pyfaidx biopython 'setuptools<81'; "
                   "then export MVA_SPLICEAI=ENV/bin/spliceai")
    need(CALIB, "regenerate from ../alphagenome_control")

    chrom, gstart, gend, tx, exons = gene_model()
    print(f"BUB1B {chrom}:{gstart:,}-{gend:,}  {tx}  {len(exons)} exons", file=sys.stderr)

    q = subprocess.run(
        ["bcftools", "query", "-r", f"{chrom}:{gstart}-{gend}",
         "-f", "%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t[%GT\t%AD\t%DP]\n", VCF],
        capture_output=True, text=True, check=True).stdout.strip().splitlines()

    variants = []
    for line in q:
        c, pos, ref, alt, qual, filt, gt, ad, dp = line.split("\t")
        kind, off = offset(int(pos), exons)
        vclass = "SNV" if (len(ref) == 1 and len(alt) == 1) else "indel"
        variants.append(dict(chrom=c if c.startswith("chr") else "chr" + c,
                             pos=int(pos), ref=ref, alt=alt, gt=gt, ad=ad, dp=dp,
                             qual=qual, filt=filt, kind=kind, off=off,
                             vclass=vclass))
    n_ex = sum(v["kind"] == "exonic" for v in variants)
    n_snv = sum(v["vclass"] == "SNV" for v in variants)
    print(f"{len(variants)} variants in the gene body: {n_ex} exonic, "
          f"{len(variants)-n_ex} intronic; {n_snv} SNVs, "
          f"{len(variants)-n_snv} indels", file=sys.stderr)
    assert len(variants) == 14, f"expected the 14 gene-body variants, got {len(variants)}"

    in_vcf = os.path.join(OUT, "bub1b_proband_snvs.vcf")
    with open(in_vcf, "w") as f:
        f.write("##fileformat=VCFv4.2\n##contig=<ID=chr15,length=101991189>\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for v in sorted(variants, key=lambda v: v["pos"]):
            f.write(f"{v['chrom']}\t{v['pos']}\t.\t{v['ref']}\t{v['alt']}\t.\t.\t.\n")
    print(f"scoring {len(variants)} proband variants", file=sys.stderr)

    out_vcf = os.path.join(OUT, f"bub1b_proband_spliceai_D{DIST}.vcf")
    subprocess.run([SPLICEAI, "-I", in_vcf, "-O", out_vcf, "-R", FASTA,
                    "-A", "grch38", "-D", str(DIST)], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    scored = {}
    for line in open(out_vcf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        m = re.search(r"SpliceAI=([^;\s]+)", f[7])
        if not m:
            continue
        p = m.group(1).split("|")
        ds = [float(x) for x in p[2:6]]
        scored[(int(f[1]), f[3], f[4])] = dict(ds=max(ds), AG=ds[0], AL=ds[1], DG=ds[2], DL=ds[3])

    import csv
    bg = []
    with open(CALIB) as f:
        for row in csv.DictReader(f):
            if row["panel_set"] == "background":
                bg.append(float(row["spliceai"]))
    bg.sort()
    n_bg = len(bg)

    rows = []
    for v in variants:
        s = scored.get((v["pos"], v["ref"], v["alt"]))
        if s is None:
            continue
        rank = sum(1 for b in bg if b >= s["ds"]) + 1
        rows.append({**v, **s, "rank": rank, "n_bg": n_bg,
                     "pct": 100 * rank / (n_bg + 1)})
    rows.sort(key=lambda r: -r["ds"])

    with open(os.path.join(OUT, "bub1b_proband_splice_scores.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    thr = 0.5                                            # SpliceAI's own recommended cutoff
    print(f"\nSpliceAI D={DIST}, ranked against {n_bg} public gnomAD intronic "
          f"variants in this gene\n")
    hdr = (f"{'position':<14}{'ref>alt':<9}{'zone':<16}{'GT':<6}{'DS':>6}"
           f"{'AG':>6}{'AL':>6}{'DG':>6}{'DL':>6}{'rank':>10}{'':>4}")
    print(hdr); print("-" * len(hdr))
    for r in rows:
        zone = "exonic" if r["kind"] == "exonic" else f"{r['kind']} {r['off']:+d}"
        flag = "  <-- DS>=0.5" if r["ds"] >= thr else ""
        ra = f"{r['ref']}>{r['alt']}"
        print(f"{r['chrom']}:{r['pos']:<9}{ra:<9}{zone:<16}{r['gt']:<6}"
              f"{r['ds']:>6.2f}{r['AG']:>6.2f}{r['AL']:>6.2f}{r['DG']:>6.2f}{r['DL']:>6.2f}"
              f"{r['rank']:>6}/{n_bg+1}{flag}")

    n_hit = sum(1 for r in rows if r["ds"] >= thr)
    print(f"\n{n_hit} of {len(rows)} at or above SpliceAI's recommended DS >= {thr}.")
    print(f"On the calibration panel that cutoff recovered 3/3 published pathogenic")
    print(f"MVA1 splice alleles and passed 6.9% of background variants.")
    print(f"\nwrote {OUT}/bub1b_proband_splice_scores.csv (gitignored: carries genotypes)")


if __name__ == "__main__":
    main()
