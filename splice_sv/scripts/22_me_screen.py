#!/usr/bin/env python3
"""Breakpoint and mobile-element screen on the targeted BUB1B BAM.

Targets the class BUB1B_STRUCTURAL.md and BUB1B_KMER_CN.md leave open: an
insertion (e.g. a Kawakami-type Alu/SVA) or a balanced event. A mobile-element
insertion leaves three signatures at its breakpoint:
  - a pileup of soft-clipped reads at one position,
  - one-end-anchored / discordant pairs whose mate maps to a mobile element or is
    unmapped,
  - clipped or mate sequence with homology to an Alu / L1 / SVA consensus.

This measures all three, then aligns the breakpoint-associated sequences to the
mobile-element consensus (bwa) to type any hit. Coverage per window is also
reported as a cross-check against the k-mer copy-number profile.

Reads the targeted BAM (proband-derived, gitignored). Writes a public summary
table (breakpoint positions and support counts, no read sequences).
"""
import os, subprocess, sys, tempfile
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.environ.get("MVA_RESULTS", os.path.join(HERE, "..", "..", "results"))
ALN = os.path.join(RES, "align")
BAM = os.path.join(ALN, "bub1b_merged.bam")
ME = os.environ.get("MVA_ME", "/Volumes/IanSSD/mva-hackathon/ref/me/me_consensus.fa")
PUBOUT = os.path.join(HERE, "..", "alphagenome_control", "results")
CHROM, GSTART, GEND, FLANK = "chr15", 40160984, 40221137, 5000
MIN_CLIP_LEN = 20          # soft-clip length to count
MIN_CLIP_SUPPORT = 5       # clipped reads at a position to flag a breakpoint
MIN_MAPQ = 20

def main():
    import pysam
    if not os.path.exists(BAM):
        sys.exit(f"MISSING {BAM}; run 20_targeted_align.sh")
    bam = pysam.AlignmentFile(BAM, "rb")

    lo, hi = GSTART - FLANK, GEND + FLANK
    clip_at = defaultdict(int)                 # breakpoint pos -> clipped read count
    clip_seqs = []                             # (pos, seq) clipped segments
    disc = []                                  # discordant/one-end-anchored records
    disc_seqs = []                             # their sequences (candidate ME/other)
    win = defaultdict(int); WIN = 500
    n_reads = 0

    for r in bam.fetch(CHROM, lo, hi):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        n_reads += 1
        if r.reference_start is not None:
            win[((r.reference_start - GSTART)//WIN)] += 1
        # soft-clips
        if r.cigartuples and r.mapping_quality >= MIN_MAPQ:
            # leading clip
            if r.cigartuples[0][0] == 4 and r.cigartuples[0][1] >= MIN_CLIP_LEN:
                clip_at[r.reference_start] += 1
                clip_seqs.append((r.reference_start, r.query_sequence[:r.cigartuples[0][1]]))
            if r.cigartuples[-1][0] == 4 and r.cigartuples[-1][1] >= MIN_CLIP_LEN:
                clip_at[r.reference_end] += 1
                clip_seqs.append((r.reference_end, r.query_sequence[-r.cigartuples[-1][1]:]))
        # discordant / one-end-anchored: mate on other chrom, far, or unmapped
        anchored_here = (r.reference_name == CHROM and lo <= r.reference_start <= hi)
        if anchored_here:
            far = (r.next_reference_name != CHROM) or \
                  (r.next_reference_start is not None and abs(r.next_reference_start - r.reference_start) > 1000) or \
                  r.mate_is_unmapped
            if far and r.mapping_quality >= MIN_MAPQ:
                disc.append((r.reference_start, r.next_reference_name, r.next_reference_start,
                             r.mate_is_unmapped))
                if r.query_sequence:
                    disc_seqs.append((r.reference_start, r.query_sequence))

    print(f"targeted BAM: {n_reads:,} primary reads in {CHROM}:{lo}-{hi}\n")

    # Coverage cross-check vs k-mer CN.
    print("=== coverage per 500bp window (cross-check k-mer CN) ===")
    covs = sorted(win.values())
    if covs:
        med = covs[len(covs)//2]
        print(f"  windows {len(win)}  median reads/window {med}  "
              f"min {covs[0]}  max {covs[-1]}")
        lowc = [ (GSTART+w*WIN) for w,c in sorted(win.items()) if c < med*0.35 ]
        print(f"  windows below 35% of median (possible deletion): "
              f"{len(lowc)} {lowc[:8] if lowc else ''}")
    print()

    # Breakpoint clip clusters.
    print("=== soft-clip breakpoint clusters ===")
    bp = sorted([(p, n) for p, n in clip_at.items() if n >= MIN_CLIP_SUPPORT],
                key=lambda x: -x[1])
    if not bp:
        print(f"  none with >= {MIN_CLIP_SUPPORT} clipped reads.\n")
    else:
        for p, n in bp[:20]:
            print(f"  chr15:{p}  {n} clipped reads")
        print()

    # Discordant summary.
    print("=== discordant / one-end-anchored pairs (anchored in BUB1B) ===")
    print(f"  {len(disc)} reads with a mate on another chrom, >1kb away, or unmapped")
    # cluster by 500bp bin
    dclust = defaultdict(int)
    for pos, mn, mp, un in disc:
        dclust[(pos//500)*500] += 1
    hot = sorted([(k, v) for k, v in dclust.items() if v >= MIN_CLIP_SUPPORT], key=lambda x: -x[1])
    for k, v in hot[:15]:
        print(f"  ~chr15:{k}-{k+499}  {v} discordant reads")
    if not hot:
        print("  no cluster >= threshold")
    print()

    # Mobile-element homology: align clipped + discordant seqs to ME consensus.
    print("=== mobile-element homology of breakpoint reads ===")
    seqs = [(f"clip_{i}_{p}", s) for i, (p, s) in enumerate(clip_seqs) if len(s) >= 20] + \
           [(f"disc_{i}_{p}", s) for i, (p, s) in enumerate(disc_seqs)]
    if not os.path.exists(ME + ".bwt"):
        subprocess.run(["bwa", "index", ME], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    me_hits = defaultdict(list)
    if seqs:
        with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fq:
            for name, s in seqs:
                fq.write(f">{name}\n{s}\n")
            fqp = fq.name
        sam = subprocess.run(["bwa", "mem", "-k", "15", "-T", "20", ME, fqp],
                             capture_output=True, text=True).stdout
        os.unlink(fqp)
        for line in sam.splitlines():
            if line.startswith("@"): continue
            f = line.split("\t")
            if len(f) < 6 or f[2] == "*": continue
            flag = int(f[1])
            if flag & 0x4: continue
            qpos = f[0].rsplit("_", 1)[-1]
            me_hits[f[2]].append((f[0], qpos, f[5]))
    if me_hits:
        for me, hits in sorted(me_hits.items(), key=lambda x: -len(x[1])):
            bp_positions = defaultdict(int)
            for name, qpos, cig in hits:
                bp_positions[qpos] += 1
            print(f"  {me}: {len(hits)} breakpoint reads match; "
                  f"positions {dict(list(bp_positions.items())[:6])}")
    else:
        print("  no breakpoint reads match Alu / L1 / SVA consensus.")
    print()

    # Verdict + floor.
    os.makedirs(PUBOUT, exist_ok=True)
    out = os.path.join(PUBOUT, "bub1b_sv_screen.csv")
    with open(out, "w") as f:
        f.write("type,position,support\n")
        for p, n in bp:
            f.write(f"softclip_breakpoint,chr15:{p},{n}\n")
        for k, v in hot:
            f.write(f"discordant_cluster,chr15:{k}-{k+499},{v}\n")
        for me, hits in me_hits.items():
            f.write(f"me_homology_{me},multiple,{len(hits)}\n")
    print(f"wrote {out}")

    print("\n=== Detection floor ===")
    print(f"  A mobile-element insertion >~300bp (Alu) leaves both a soft-clip pileup")
    print(f"  and a discordant/ME-homologous cluster; either would appear above at")
    print(f"  >= {MIN_CLIP_SUPPORT} reads given ~40x depth. Not detectable: events fully")
    print(f"  contained within one read, and small indels already in the VCF. A")
    print(f"  balanced translocation shows as a discordant cluster to another chrom.")
    if not bp and not hot and not me_hits:
        print("\n  VERDICT: no breakpoint, discordant cluster, or ME homology in BUB1B.")
    else:
        print("\n  VERDICT: signal(s) above -- inspect the clusters and ME hits.")

if __name__ == "__main__":
    main()
