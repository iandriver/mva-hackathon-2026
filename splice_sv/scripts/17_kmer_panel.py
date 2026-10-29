#!/usr/bin/env python3
"""Build the k-mer panel for an alignment-free copy-number profile of BUB1B.

Goal: detect a heterozygous exon-level deletion in the 21 BUB1B exons that carry
no depth information in the VCF (see findings/BUB1B_STRUCTURAL.md). A het deletion
shows as a contiguous run of windows at ~0.5x the gene's own coverage.

Method is self-normalizing and needs no whole-genome download:
  target   31-mers tiling the BUB1B gene body, binned into windows
  control  31-mers from several unique chr15 regions of known diploid copy number,
           to anchor the 1.0 baseline and measure count noise

Uniqueness is enforced two ways, because a repeat/paralog k-mer inflates counts:
  1. drop any 31-mer occurring more than once within chr15 (self-collision)
  2. drop low-complexity 31-mers (a single base >=70%, or a 2-mer period)
The residual genome-wide-repeat risk is handled at count time by the k-mer
coverage histogram (repeat k-mers form a high-count tail), and stated as a caveat.

Reference only. No reads, no proband data. Writes the panel as JSON.
"""
import json, os, sys
from collections import Counter

REF = os.environ.get("MVA_REF", "/Volumes/IanSSD/mva-hackathon/ref")
FASTA = os.path.join(REF, "chr15.fa.bgz")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "alphagenome_control", "data")
K = 31
# BUB1B gene body, GRCh38 (Ensembl); +/- 2 kb flank so the profile has context.
GENE_START, GENE_END = 40160984, 40221137
FLANK = 2000
WIN = 500                              # window size for the depth profile
# Control: unique single-copy chr15 windows well away from the gene and centromere.
CONTROL_REGIONS = [(35_000_000, 35_050_000), (55_000_000, 55_050_000),
                   (75_000_000, 75_050_000), (90_000_000, 90_050_000)]

COMP = str.maketrans("ACGT", "TGCA")
def canon(km):
    rc = km.translate(COMP)[::-1]
    return km if km <= rc else rc

def low_complexity(km):
    c = Counter(km)
    if max(c.values()) / len(km) >= 0.70:
        return True
    if len(set(km[i:i+2] for i in range(0, len(km)-1, 2))) <= 2:  # dinucleotide period
        return True
    return False

def main():
    try:
        import pysam
    except ImportError:
        sys.exit("need pysam: use the spliceai venv or pip install pysam")
    if not os.path.exists(FASTA):
        sys.exit(f"MISSING: {FASTA}")
    fa = pysam.FastaFile(FASTA)

    # Pass 1: count every canonical 31-mer across the whole gene+flank+controls so
    # we can drop the ones that collide within chr15.
    spans = [("gene", GENE_START - FLANK, GENE_END + FLANK)] + \
            [(f"ctrl{i}", a, b) for i, (a, b) in enumerate(CONTROL_REGIONS)]
    global_counts = Counter()
    seqs = {}
    for name, a, b in spans:
        s = fa.fetch("chr15", a - 1, b).upper()
        seqs[name] = (a, s)
        for i in range(len(s) - K + 1):
            km = s[i:i+K]
            if "N" in km:
                continue
            global_counts[canon(km)] += 1

    def usable(km):
        return global_counts[canon(km)] == 1 and not low_complexity(km)

    # Target: assign each usable gene k-mer to a window by its start position.
    a0, gseq = seqs["gene"]
    target = {}                        # window_start -> list of canonical kmers
    kept = dropped_dup = dropped_lc = 0
    for i in range(len(gseq) - K + 1):
        pos = a0 + i                   # 1-based genomic position of k-mer start
        if not (GENE_START <= pos <= GENE_END):
            continue
        km = gseq[i:i+K]
        if "N" in km:
            continue
        if global_counts[canon(km)] != 1:
            dropped_dup += 1; continue
        if low_complexity(km):
            dropped_lc += 1; continue
        w = ((pos - GENE_START) // WIN) * WIN + GENE_START
        target.setdefault(w, []).append(canon(km))
        kept += 1

    control = []
    for name in seqs:
        if name == "gene":
            continue
        a, s = seqs[name]
        for i in range(len(s) - K + 1):
            km = s[i:i+K]
            if "N" in km or global_counts[canon(km)] != 1 or low_complexity(km):
                continue
            control.append(canon(km))

    # Exon annotation for the windows, from Ensembl (so exon labels are not hardcoded),
    # with retries and a fallback to the exon table already on disk.
    import urllib.request, time as _t
    def _get(url):
        for a in range(4):
            try:
                return json.load(urllib.request.urlopen(url, timeout=60))
            except Exception as e:
                if a == 3:
                    raise
                _t.sleep(3 * (a + 1))
    tx, exons = None, None
    try:
        tx = _get("https://rest.ensembl.org/lookup/symbol/homo_sapiens/BUB1B"
                  "?content-type=application/json")["canonical_transcript"].split(".")[0]
        t = _get(f"https://rest.ensembl.org/lookup/id/{tx}"
                 "?expand=1;content-type=application/json")
        exons = sorted(((e["start"], e["end"]) for e in t["Exon"]))
    except Exception as e:
        cached = os.path.join(OUT, "exons.json")
        if not os.path.exists(cached):
            raise
        ex = json.load(open(cached))
        exons = sorted((e["start"], e["end"]) for e in ex)
        tx = "ENST00000287598 (from cached exons.json; Ensembl unreachable)"
        print(f"Ensembl unreachable ({e}); used cached {cached}", file=sys.stderr)

    win_meta = []
    for w in sorted(target):
        we = min(w + WIN - 1, GENE_END)
        ex = [i for i, (s, e) in enumerate(exons, 1) if not (e < w or s > we)]
        win_meta.append(dict(start=w, end=we, n_kmers=len(target[w]), exons=ex))

    os.makedirs(OUT, exist_ok=True)
    panel = dict(k=K, win=WIN, gene=("chr15", GENE_START, GENE_END),
                 transcript=tx, n_exons=len(exons),
                 target={str(w): kms for w, kms in target.items()},
                 windows=win_meta, control=control)
    with open(os.path.join(OUT, "kmer_panel.json"), "w") as f:
        json.dump(panel, f)

    n_win_kmers = sum(len(v) for v in target.values())
    print(f"target windows: {len(target)}  ({n_win_kmers:,} unique k-mers)")
    print(f"  kept {kept:,}   dropped self-collision {dropped_dup:,}   low-complexity {dropped_lc:,}")
    print(f"control k-mers: {len(control):,} across {len(CONTROL_REGIONS)} chr15 regions")
    empty = [wm for wm in win_meta if wm["n_kmers"] == 0]
    exon_win = [wm for wm in win_meta if wm["exons"]]
    print(f"windows overlapping an exon: {len(exon_win)}/{len(win_meta)}")
    if empty:
        print(f"WARNING: {len(empty)} windows have no usable k-mer (repeat-dense); "
              f"those cannot be profiled")
    print(f"\nwrote {os.path.join(OUT, 'kmer_panel.json')}")

if __name__ == "__main__":
    main()
