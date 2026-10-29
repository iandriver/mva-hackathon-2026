#!/usr/bin/env python3
"""Detection floor for scripts/14_splice_intronic.py.

14 returned no splice signal for any of the proband's fourteen BUB1B gene-body
variants. That is only a finding if SpliceAI could have seen a signal there.
Twelve of the fourteen are deep intronic, 905 to 3653 bp from the nearest exon,
and the positive control validated near-splice detection, not deep intronic
detection. So the calibration panel does not license the negative on its own.

This scans every possible SNV in a window around each proband variant and asks
what the highest attainable score in that neighbourhood is. If nothing there can
score, the model is blind at these positions and the negative means little. If
neighbouring changes do score, the model is reading sequence there and the
proband's quiet scores are informative.

Public reference sequence only. Proband positions define the window centres, so
outputs go to MVA_RESULTS, which is gitignored.
"""
import csv, os, re, subprocess, sys

DATA = os.environ.get("MVA_DATA", "/Volumes/IanSSD/mva-hackathon/data")
REF = os.environ.get("MVA_REF", "/Volumes/IanSSD/mva-hackathon/ref")
FASTA = os.path.join(REF, "chr15.fa.bgz")
OUT = os.environ.get("MVA_RESULTS",
                     os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                  "..", "..", "results"))  # repo-root results/, gitignored
SPLICEAI = os.environ.get(
    "MVA_SPLICEAI",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                 "alphagenome_control", ".venv-spliceai", "bin", "spliceai"))
WINDOW = int(os.environ.get("MVA_FLOOR_WINDOW", "50"))
BASES = "ACGT"


def main():
    scores_csv = os.path.join(OUT, "bub1b_proband_splice_scores.csv")
    if not os.path.exists(scores_csv):
        sys.exit(f"MISSING: {scores_csv}\n  run scripts/14_splice_intronic.py first")
    import pysam
    fa = pysam.FastaFile(FASTA)

    obs = list(csv.DictReader(open(scores_csv)))
    centres = [(int(r["pos"]), r["kind"], int(r["off"]), float(r["ds"])) for r in obs]

    seen, rows = set(), []
    for pos, kind, off, _ in centres:
        for p in range(pos - WINDOW, pos + WINDOW + 1):
            ref = fa.fetch("chr15", p - 1, p).upper()
            if ref not in BASES:
                continue
            for alt in BASES:
                if alt == ref or (p, ref, alt) in seen:
                    continue
                seen.add((p, ref, alt))
                rows.append((p, ref, alt, pos))
    print(f"scanning {len(rows):,} possible SNVs, +/-{WINDOW} bp around "
          f"{len(centres)} proband variants", file=sys.stderr)

    in_vcf = os.path.join(OUT, "bub1b_floor_scan.vcf")
    with open(in_vcf, "w") as f:
        f.write("##fileformat=VCFv4.2\n##contig=<ID=chr15,length=101991189>\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for p, ref, alt, _ in sorted(rows):
            f.write(f"chr15\t{p}\t.\t{ref}\t{alt}\t.\t.\t.\n")

    out_vcf = os.path.join(OUT, "bub1b_floor_scan_spliceai.vcf")
    subprocess.run([SPLICEAI, "-I", in_vcf, "-O", out_vcf, "-R", FASTA,
                    "-A", "grch38", "-D", "50"], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    scored = {}
    for line in open(out_vcf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        m = re.search(r"SpliceAI=([^;\s]+)", f[7])
        if not m:
            continue
        scored[(int(f[1]), f[3], f[4])] = max(float(x) for x in m.group(1).split("|")[2:6])

    by_centre = {}
    for p, ref, alt, centre in rows:
        s = scored.get((p, ref, alt))
        if s is not None:
            by_centre.setdefault(centre, []).append(s)

    print(f"\nWhat SpliceAI could have detected within +/-{WINDOW} bp "
          f"({len(scored):,} SNVs scored)\n")
    hdr = (f"{'centre':<12}{'zone':<17}{'observed':>9}{'nbhd max':>10}"
           f"{'nbhd p99':>10}{'n>=0.5':>8}{'n>=0.2':>8}  verdict")
    print(hdr); print("-" * len(hdr))
    blind = []
    for pos, kind, off, ds in sorted(centres, key=lambda c: c[0]):
        v = sorted(by_centre.get(pos, []))
        if not v:
            continue
        n = len(v)
        mx = v[-1]
        p99 = v[int(0.99 * (n - 1))]
        n5 = sum(1 for x in v if x >= 0.5)
        n2 = sum(1 for x in v if x >= 0.2)
        zone = "exonic" if kind == "exonic" else f"{kind} {off:+d}"
        verdict = "informative" if mx >= 0.5 else ("weak" if mx >= 0.2 else "BLIND")
        if verdict == "BLIND":
            blind.append(pos)
        print(f"chr15:{pos:<6}{zone:<17}{ds:>9.2f}{mx:>10.2f}{p99:>10.2f}"
              f"{n5:>8}{n2:>8}  {verdict}")

    allv = [s for v in by_centre.values() for s in v]
    print(f"\nAcross the whole scan: max {max(allv):.2f}, "
          f"{sum(1 for x in allv if x >= 0.5)} of {len(allv):,} at DS >= 0.5 "
          f"({100*sum(1 for x in allv if x >= 0.5)/len(allv):.2f}%).")
    if blind:
        print(f"\n{len(blind)} of {len(centres)} neighbourhoods are BLIND: no possible SNV")
        print("within the window reaches DS 0.2. At those positions a negative result")
        print("carries no information, because nothing there could have scored.")
    else:
        print("\nEvery neighbourhood contains changes SpliceAI scores highly, so the")
        print("model is reading sequence at all of these positions and the proband's")
        print("quiet scores are informative rather than an absence of sensitivity.")
    print(f"\nwrote {OUT}/bub1b_floor_scan_spliceai.vcf (gitignored)")


if __name__ == "__main__":
    main()
