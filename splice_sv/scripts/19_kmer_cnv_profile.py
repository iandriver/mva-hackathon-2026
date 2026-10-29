#!/usr/bin/env python3
"""Copy-number profile of BUB1B from the streamed k-mer counts (scripts 17, 18).

Detects a heterozygous exon-level deletion in the BUB1B exons that carry no depth
information in the VCF (findings/BUB1B_STRUCTURAL.md). A het deletion shows as a
window at ~half the gene's diploid per-k-mer depth.

Normalization is SELF-referential, to the gene's own diploid level, for two
reasons found empirically:
  - The external control regions (script 17) are repeat-contaminated: their
    "unique-within-region" k-mers recur genome-wide, giving a baseline ~25x too
    high. They are reported here only as the diagnostic that revealed this, not
    used for normalization.
  - Many intronic windows are themselves repeat-dense (per-k-mer depth orders of
    magnitude above the diploid mode). These are flagged UNINTERPRETABLE rather
    than reported as high copy number.

The diploid baseline is the median per-k-mer depth of the exon windows, which
resists the minority of repeat-inflated windows. It is CN 2 because the VCF shows
normal-depth heterozygous SNVs across the gene, excluding a uniform loss.

Reads the gitignored per-file counts; writes a public CN profile (no genotypes).
"""
import glob, json, math, os, sys
from statistics import median

HERE = os.path.dirname(os.path.abspath(__file__))
PANEL = os.path.join(HERE, "..", "alphagenome_control", "data", "kmer_panel.json")
RES = os.environ.get("MVA_RESULTS", os.path.join(HERE, "..", "..", "results"))
COUNTD = os.path.join(RES, "kmer")
PUBOUT = os.path.join(HERE, "..", "alphagenome_control", "results")
REPEAT_FACTOR = 3.0        # window per-k-mer depth above this * baseline = repeat-dense

def main():
    panel = json.load(open(PANEL))
    win_meta = {wm["start"]: wm for wm in panel["windows"]}

    files = sorted(glob.glob(os.path.join(COUNTD, "count_*.json")))
    if len(files) < 8:
        print(f"WARNING: only {len(files)}/8 per-file counts present", file=sys.stderr)
    if not files:
        sys.exit(f"no counts in {COUNTD}; run script 18 first")

    win_tot, ctrl_tot, n_reads, nck = {}, 0, 0, None
    for f in files:
        d = json.load(open(f))
        n_reads += d["n_reads"]; ctrl_tot += d["control_total"]; nck = d["n_control_kmers"]
        for w, c in d["window_counts"].items():
            win_tot[int(w)] = win_tot.get(int(w), 0) + c
    print(f"merged {len(files)} files, {n_reads:,} reads")

    pk = {w: win_tot.get(w, 0) / win_meta[w]["n_kmers"] for w in win_meta}

    # Diploid baseline: median per-k-mer depth of exon windows (robust to a
    # minority of repeat-inflated windows).
    exon_w = sorted(w for w, wm in win_meta.items() if wm["exons"])
    baseline = median(pk[w] for w in exon_w)
    repeat_cut = REPEAT_FACTOR * baseline
    print(f"diploid baseline (median exon per-k-mer depth): {baseline:.1f} = CN 2.0")
    print(f"repeat-dense cutoff: per-k-mer depth > {repeat_cut:.0f}\n")

    # External controls: diagnostic only.
    print(f"external control k-mers: {ctrl_tot/nck:.0f} hits/k-mer = {2*(ctrl_tot/nck)/baseline:.0f}x "
          f"the diploid level -> repeat-contaminated, NOT used for normalization\n")

    rows = []
    for w in sorted(win_meta):
        wm = win_meta[w]
        cn = 2 * pk[w] / baseline
        repeat = pk[w] > repeat_cut
        rows.append(dict(start=w, end=wm["end"], exons=wm["exons"], n_kmers=wm["n_kmers"],
                         per_kmer=pk[w], cn=cn, repeat=repeat))

    # The question: any assessable exon window at ~half the diploid level?
    print("=== Exon windows (the question: any at CN ~1.0?) ===\n")
    print(f"  {'window':<24}{'exon':<8}{'per_kmer':>9}{'CN':>7}  status")
    exon_assessable, exon_repeat, low_exons = [], [], []
    exons_seen = {}
    for r in rows:
        if not r["exons"]:
            continue
        for e in r["exons"]:
            exons_seen.setdefault(e, []).append(r)
        tag = "REPEAT-DENSE (uninterpretable)" if r["repeat"] else ""
        if not r["repeat"]:
            exon_assessable.append(r)
            if r["cn"] < 1.4:
                low_exons.append(r); tag = "<-- below diploid"
        else:
            exon_repeat.append(r)
        ex = ",".join(map(str, r["exons"]))
        print(f"  chr15:{r['start']}-{r['end']:<7}{ex:<8}{r['per_kmer']:>9.0f}{r['cn']:>7.2f}  {tag}")

    assessable_cn = sorted(r["cn"] for r in exon_assessable)
    print(f"\nassessable exon windows: {len(exon_assessable)}  "
          f"CN min {assessable_cn[0]:.2f}  median {median(assessable_cn):.2f}  max {assessable_cn[-1]:.2f}")

    # Which whole exons have NO assessable window (every window repeat-dense)?
    exon_ids = sorted(exons_seen)
    unassessable_exons = [e for e in exon_ids
                          if all(r["repeat"] for r in exons_seen[e])]
    print(f"exons with at least one assessable window: "
          f"{len(exon_ids)-len(unassessable_exons)}/{len(exon_ids)}")
    if unassessable_exons:
        print(f"exons uninterpretable (all windows repeat-dense): {unassessable_exons}")

    # Detection floor for the assessable exon windows.
    floors = []
    for r in exon_assessable:
        lam = r["n_kmers"] * baseline          # expected hits at CN2
        floors.append(3 * 2.0 / math.sqrt(lam) if lam > 0 else float("inf"))
    print(f"\n=== Detection floor ===")
    print(f"  Median 3-sigma CN resolution across assessable exon windows: {median(floors):.3f}")
    print(f"  A heterozygous deletion (CN 2->1) is far above this, so any assessable")
    print(f"  exon carrying one would read near CN 1.0.")

    # Write public CN table.
    os.makedirs(PUBOUT, exist_ok=True)
    out = os.path.join(PUBOUT, "bub1b_kmer_cn_profile.csv")
    with open(out, "w") as f:
        f.write("start,end,exons,n_kmers,per_kmer_depth,copy_number,repeat_dense\n")
        for r in rows:
            f.write(f"{r['start']},{r['end']},{'|'.join(map(str,r['exons']))},{r['n_kmers']},"
                    f"{r['per_kmer']:.2f},{r['cn']:.3f},{int(r['repeat'])}\n")
    print(f"\nwrote {out}")

    if low_exons:
        print(f"\nVERDICT: {len(low_exons)} assessable exon window(s) below CN 1.4 - INSPECT:")
        for r in low_exons:
            print(f"  exons {r['exons']} at CN {r['cn']:.2f}")
    else:
        print(f"\nVERDICT: no assessable exon window is at hemizygous depth. No evidence")
        print(f"  of a heterozygous exon-level deletion in the {len(exon_ids)-len(unassessable_exons)} "
              f"assessable exons.")
        if unassessable_exons:
            print(f"  Exons {unassessable_exons} remain uninterpretable (repeat-dense k-mers).")

if __name__ == "__main__":
    main()
