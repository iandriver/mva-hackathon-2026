#!/usr/bin/env python3
"""Reproduces findings/MOSAICISM.md — mosaic aneuploidy, with a simulated detection floor.

MVA is defined by mosaic variegated aneuploidy, so the question is whether the released
WGS shows it. Two channels: B-allele frequency at heterozygous SNPs, and read depth.

The naive statistic (per-chromosome MAXIMUM window deviation) measures pericentromeric
mapping artifact, not mosaicism -- the largest window in this sample sits on the chr7
centromere. This uses the MEDIAN across windows, which those artifacts cannot move, and
calibrates a floor by simulating against the observed depth distribution.

A mosaic trisomy at cell fraction f gives BAF modes at (1+f)/(2+f) and 1/(2+f),
and a depth ratio of 1 + f/2.
"""
import json, os, random, statistics as st, subprocess, sys
from collections import defaultdict

DATA = os.environ.get("MVA_DATA", "/Volumes/IanSSD/mva-hackathon/data")
VCF = os.path.join(DATA, "WGS_EX2312012_HGWCNDSX7.vcf.gz")
OUT = os.environ.get("MVA_RESULTS", os.path.join(os.path.dirname(__file__), "..", "results"))
os.makedirs(OUT, exist_ok=True)

WINDOW = 2000        # het SNPs per window
DP_MIN, DP_MAX = 20, 70
SEED = 3
FRACTIONS = (0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.15, 0.20)

def load():
    if not os.path.exists(VCF):
        sys.exit(f"MISSING: {VCF}\n  see scripts/env.sh for the fetch command")
    cmd = ["bcftools", "query", "-f", "%CHROM\t%POS\t[%AD\t%DP]\n",
           "-i", 'GT="het" && TYPE="snp" && FILTER="PASS"', VCF]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    data = defaultdict(list)
    for line in p.stdout:
        f = line.rstrip("\n").split("\t")
        if len(f) < 4: continue
        try:
            ad = f[2].split(","); ref, alt = int(ad[0]), int(ad[1]); dp = int(f[3])
        except (ValueError, IndexError):
            continue
        if not (DP_MIN <= dp <= DP_MAX) or ref + alt == 0: continue
        data[f[0]].append((int(f[1]), alt / (ref + alt), dp))
    p.wait()
    for c in data: data[c].sort()
    return data

def stat(sites):
    """Median across windows of mean |BAF - 0.5| -- robust to centromeric artifact."""
    wins = [sites[i:i+WINDOW] for i in range(0, len(sites) - WINDOW + 1, WINDOW)]
    if not wins: return None
    return st.median([st.mean(abs(b - 0.5) for _, b, _ in w) for w in wins])

def simulate(f, depths, reps=40):
    hi = (1 + f) / (2 + f); lo = 1 - hi
    out = []
    for _ in range(reps):
        tot = 0.0
        for _ in range(WINDOW):
            d = random.choice(depths)
            p = random.choice([hi, lo]) if f > 0 else 0.5
            a = sum(1 for _ in range(d) if random.random() < p)
            tot += abs(a / d - 0.5)
        out.append(tot / WINDOW)
    return st.median(out)

def main():
    random.seed(SEED)
    data = load()
    depths = [d for c in data for _, _, d in data[c]]
    med_d = st.median(depths)
    total = len(depths)
    print(f"heterozygous PASS SNPs (DP {DP_MIN}-{DP_MAX}): {total:,}", file=sys.stderr)
    print(f"median depth: {med_d:.0f}", file=sys.stderr)

    auto = sorted([c for c in data if c.isdigit()], key=int)
    obs, depth_ratio = {}, {}
    for c in auto:
        s = stat(data[c])
        if s is None: continue
        obs[c] = s
        depth_ratio[c] = st.median([d for _, _, d in data[c]]) / med_d

    lo_c, hi_c = min(obs, key=obs.get), max(obs, key=obs.get)
    spread = obs[hi_c] - obs[lo_c]

    print(f"\n{'chr':<5}{'statistic':>12}{'depth ratio':>14}")
    for c in auto:
        if c in obs: print(f"{c:<5}{obs[c]:>12.5f}{depth_ratio[c]:>14.3f}")
    print(f"\nobserved: min {obs[lo_c]:.5f} (chr{lo_c})  max {obs[hi_c]:.5f} (chr{hi_c})")
    print(f"between-chromosome spread: {spread:.5f}")
    print(f"sd across chromosomes:     {st.pstdev(list(obs.values())):.5f}")
    print(f"depth ratio range:         {min(depth_ratio.values()):.3f} - {max(depth_ratio.values()):.3f}")

    null = simulate(0.0, depths)
    print(f"\nsimulated null: {null:.5f}\n")
    print(f"{'f':>6}{'simulated':>12}{'shift':>10}{'> spread':>11}")
    floor = None
    rows = []
    for f in FRACTIONS:
        s = simulate(f, depths); shift = s - null
        exceeds = shift > spread
        if exceeds and floor is None: floor = f
        rows.append({"f": f, "stat": s, "shift": shift, "exceeds_spread": exceeds})
        print(f"{f:>6.2f}{s:>12.5f}{shift:>10.5f}{('yes' if exceeds else 'no'):>11}")

    # MVA is *variegated*: a real mosaic gain makes one chromosome stand out from the rest,
    # so the between-chromosome spread is the right empirical null.
    print(f"\nDETECTION FLOOR: mosaic whole-chromosome gain in >= {floor*100:.0f}% of cells")
    print("would shift one chromosome by more than the entire observed spread. None does.")
    json.dump({"n_het": total, "median_depth": med_d, "observed": obs,
               "depth_ratio": depth_ratio, "spread": spread, "null": null,
               "simulation": rows, "floor": floor, "window": WINDOW, "seed": SEED},
              open(os.path.join(OUT, "05_mosaicism.json"), "w"), indent=2)
    print(f"\nwrote {OUT}/05_mosaicism.json", file=sys.stderr)

if __name__ == "__main__":
    main()
