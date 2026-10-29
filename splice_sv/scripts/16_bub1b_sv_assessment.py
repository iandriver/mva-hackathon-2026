#!/usr/bin/env python3
"""What the release can and cannot say about structural variants in BUB1B.

MECHANISM.md 2 records that the Kawakami MVA1 second allele was a 3.0 kb Alu/SVA
insertion in intron 2, invisible to short-read genome sequencing and found only
by RNA-seq plus adaptive-sampling long reads. SPLICE_INTRONIC.md left structural
variants as the untested class. This asks what this VCF supports.

Three parts:
  1. Whether structural variants were called at all.
  2. The largest event the caller actually emitted anywhere, as an empirical
     size ceiling.
  3. Whether a large heterozygous deletion is excluded. It would show as halved
     read depth and a loss of heterozygous calls, both visible in a small-variant
     VCF. Where there are no variant positions there is no depth information, so
     the uncovered windows are mapped and the exons inside them named.

Absence is stored with a reason and a power number, never inferred.
"""
import json, os, subprocess, sys, urllib.request
from statistics import median

DATA = os.environ.get("MVA_DATA", "/Volumes/IanSSD/mva-hackathon/data")
VCF = os.path.join(DATA, "WGS_EX2312012_HGWCNDSX7.vcf.gz")
OUT = os.environ.get("MVA_RESULTS",
                     os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                  "..", "..", "results"))  # repo-root results/, gitignored
os.makedirs(OUT, exist_ok=True)
SV_INFO = {"SVTYPE", "SVLEN", "END", "CIPOS", "CIEND", "IMPRECISE", "MATEID", "CN"}


def q(args):
    return subprocess.run(args, capture_output=True, text=True, check=True).stdout


def gene_model(symbol="BUB1B"):
    u = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{symbol}?content-type=application/json"
    g = json.load(urllib.request.urlopen(u, timeout=60))
    assert g["assembly_name"] == "GRCh38", f"expected GRCh38, got {g['assembly_name']}"
    tx = g["canonical_transcript"].split(".")[0]
    t = json.load(urllib.request.urlopen(
        f"https://rest.ensembl.org/lookup/id/{tx}?expand=1;content-type=application/json",
        timeout=60))
    exons = sorted(({"start": e["start"], "end": e["end"]} for e in t["Exon"]),
                   key=lambda e: e["start"])
    return g["seq_region_name"], g["start"], g["end"], exons


def main():
    if not os.path.exists(VCF):
        sys.exit(f"MISSING: {VCF}\n  see scripts/env.sh")

    print("=== 1. Were structural variants called? ===\n")
    hdr = q(["bcftools", "view", "-h", VCF])
    declared = {l.split("ID=")[1].split(",")[0] for l in hdr.splitlines()
                if l.startswith("##INFO=")}
    sv_fields = sorted(declared & SV_INFO)
    alt_lines = [l for l in hdr.splitlines() if l.startswith("##ALT=")]
    src = [l for l in hdr.splitlines() if l.startswith(("##source=", "##GATKCommandLine"))]
    print(f"  pipeline: {'GATK' if any('GATK' in s for s in src) else 'unknown'}"
          f" ({len(src)} tool records)")
    print(f"  SV INFO fields declared: {sv_fields or 'none'}")
    print(f"  symbolic ALT definitions: {len(alt_lines) or 'none'}")
    n_sym = int(q(["bash", "-c",
                   f"bcftools view -H {VCF} | awk -F'\\t' '$5 ~ /[<\\[\\]]/' | wc -l"]).strip())
    print(f"  records with a symbolic or breakend ALT: {n_sym}")
    called = bool(sv_fields) or bool(alt_lines) or n_sym > 0
    print(f"\n  VERDICT: structural variants were {'called' if called else 'NEVER CALLED'}."
          f" This is a small-variant callset.")

    print("\n=== 2. Largest event the caller emitted, genome-wide ===\n")
    sizes = q(["bash", "-c",
               f"bcftools query -f '%REF\\t%ALT\\n' {VCF} | "
               "awk -F'\\t' '{n=length($2)-length($1); if(n!=0) print (n<0?-n:n)}'"]).split()
    sizes = sorted(int(s) for s in sizes)
    over = [s for s in sizes if s > 100]
    print(f"  {len(sizes):,} indels   largest {sizes[-1]} bp   "
          f"{len(over):,} over 100 bp ({100*len(over)/len(sizes):.2f}%)")
    print(f"  An insertion the size of the Kawakami allele (3.0 kb) is "
          f"{'within' if sizes[-1] >= 3000 else 'ABOVE'} this caller's demonstrated range.")

    print("\n=== 3. Is a large heterozygous deletion excluded in BUB1B? ===\n")
    chrom, gstart, gend, exons = gene_model()
    bg = [int(x) for x in q(["bcftools", "query", "-r", chrom, "-f", "[%DP]\n", VCF]).split()
          if x.isdigit()]
    bg.sort()
    bg_med = median(bg)
    print(f"  chr15 background depth: median {bg_med:.0f}  "
          f"p10 {bg[int(.10*len(bg))]}  p90 {bg[int(.90*len(bg))]}  (n={len(bg):,})")
    print(f"  a heterozygous deletion would put depth near {bg_med/2:.0f}\n")

    rows = []
    for line in q(["bcftools", "query", "-r", f"{chrom}:{gstart}-{gend}",
                   "-f", "%POS\t%REF\t%ALT\t[%GT]\t[%DP]\t%FILTER\n", VCF]).strip().splitlines():
        pos, ref, alt, gt, dp, filt = line.split("\t")
        rows.append(dict(pos=int(pos), gt=gt, dp=int(dp) if dp.isdigit() else None,
                         filt=filt, het=gt in ("0/1", "1/0", "0|1", "1|0")))
    ok = [r for r in rows if r["dp"] and r["filt"] == "PASS"]
    print(f"  {len(rows)} variants in the gene body, {len(ok)} PASS with depth")
    print(f"  their depth: median {median(r['dp'] for r in ok):.0f}, "
          f"range {min(r['dp'] for r in ok)}-{max(r['dp'] for r in ok)}")
    n_het = sum(r["het"] for r in ok)
    print(f"  heterozygous calls: {n_het} of {len(ok)}. A heterozygous call cannot sit")
    print(f"  inside a heterozygous deletion, so each one excludes a deletion there")
    print(f"  regardless of depth.")

    # Only non-heterozygous positions can be hemizygous. Judge them against the
    # background spread rather than a fixed fraction of the median: low depth is
    # common genome-wide and a bare threshold manufactures hits.
    pct = lambda d: 100 * sum(1 for b in bg if b <= d) / len(bg)
    hom = [r for r in ok if not r["het"]]
    print(f"\n  Of the {len(hom)} homozygous PASS positions, depth vs the chr15 "
          f"distribution:")
    susp = []
    for r in sorted(hom, key=lambda r: r["dp"]):
        p_here = pct(r["dp"])
        tag = ""
        if r["dp"] < bg_med / 2 * 1.3 and p_here < 5:
            tag = "  <-- consistent with hemizygosity"
            susp.append(r)
        print(f"    {r['pos']}  depth {r['dp']:>3}  ({p_here:4.1f}th pct of chr15){tag}")
    if not susp:
        print(f"\n  None. The lowest is at the {pct(min(r['dp'] for r in hom)):.0f}th "
              f"percentile of chr15, which is ordinary coverage, not half-coverage.")

    print("\n  Windows with NO informative position, where a deletion could hide:\n")
    marks = sorted({gstart - 1} | {r["pos"] for r in ok} | {gend + 1})
    gaps = []
    for a, b in zip(marks, marks[1:]):
        span = b - a - 1
        if span < 1000:
            continue
        inside = [i for i, e in enumerate(exons, 1) if a < e["start"] and e["end"] < b]
        gaps.append((a, b, span, inside))
    gaps.sort(key=lambda g: -g[2])
    print(f"    {'window':<26}{'span':>8}   exons fully inside")
    for a, b, span, inside in gaps:
        ex = f"{len(inside)}: {inside[0]}-{inside[-1]}" if inside else "none"
        print(f"    {a+1}-{b-1:<12} {span:>7,}   {ex}")
    covered = {i for i, e in enumerate(exons, 1)} - {i for _, _, _, ins in gaps for i in ins}
    print(f"\n    {len(exons)-len(covered)} of {len(exons)} exons sit inside an "
          f"uninformative window.")

    print("\n=== Detection floor ===\n")
    print(f"  Excluded: any heterozygous deletion overlapping one of the {len(ok)} PASS")
    print(f"  variant positions, since all show depth near the chr15 median and "
          f"{n_het} are heterozygous.")
    print(f"  Not excluded: a deletion confined to one of the {len(gaps)} windows above,")
    print(f"  the largest being {gaps[0][2]:,} bp. Also not excluded: any insertion, any")
    print(f"  balanced event, and any deletion the caller never attempted to find.")
    print("\n  Resolving this needs the alignments. The release ships FASTQs only")
    print("  (8 files, ~85 GB), so it requires aligning from reads.")


if __name__ == "__main__":
    main()
