"""Does AlphaGenome add information beyond 'how close is it to a splice site'?

Two comparisons:
  1. A proximity baseline, score = -|offset|, ranked over the same background.
  2. A within-stratum test: rank each positive only against background variants
     at the SAME offset, which removes position entirely.
"""
import pathlib
import numpy as np, pandas as pd

w = pd.read_csv(pathlib.Path(__file__).parent/'results'/'per_variant_scores.csv')
bg = w[w.panel_set=='background']; pos = w[w.panel_set=='positive']
M = 'SPLICE_SITES'

print("=== 1. Proximity baseline vs AlphaGenome, same background (n=463) ===\n")
print(f"{'allele':<22} {'offset':>7} {'AlphaGenome':>14} {'proximity':>14}")
for _, p in pos.iterrows():
    ag  = int((bg[M] >= p[M]).sum()) + 1
    # proximity: smaller |offset| ranks higher; ties counted as beating it (conservative)
    prox = int((bg.panel_off.abs() <= abs(p.panel_off)).sum()) + 1
    print(f"{p.panel_name:<22} {int(p.panel_off):>+7} {ag:>8}/464    {prox:>8}/464")
print()
print("    Proximity ties are counted in the baseline's favour, so these are")
print("    upper bounds on the baseline's rank (i.e. best case for it).\n")

print("=== 2. Within-offset test: same distance from the splice site ===\n")
for _, p in pos.iterrows():
    same = bg[bg.panel_off == p.panel_off]
    if len(same) == 0:
        print(f"{p.panel_name}: no background at offset {int(p.panel_off):+d}"); continue
    rank = int((same[M] >= p[M]).sum()) + 1
    print(f"{p.panel_name:<22} offset {int(p.panel_off):+3d}   "
          f"rank {rank:>3}/{len(same)+1:<4} among same-offset variants   "
          f"(their median {same[M].median():.4f}, this {p[M]:.4f})")
print()

print("=== 3. Is offset alone predictive of AlphaGenome score in the background? ===\n")
from scipy.stats import spearmanr
r, pv = spearmanr(bg.panel_off.abs(), bg[M])
print(f"    Spearman(|offset|, {M}) over background: rho = {r:.3f}, p = {pv:.2g}")
print(f"    median {M} by offset band:")
for lo, hi in [(3,5),(6,10),(11,15),(16,20)]:
    s = bg[(bg.panel_off.abs()>=lo) & (bg.panel_off.abs()<=hi)]
    print(f"      |off| {lo:>2}-{hi:<2}  n={len(s):<4} median {s[M].median():.4f}  "
          f"p95 {s[M].quantile(.95):.4f}  max {s[M].max():.4f}")
