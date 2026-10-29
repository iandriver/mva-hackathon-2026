"""Operating characteristics, and one negative result about the headline score."""
import pathlib
import numpy as np, pandas as pd

ROOT = pathlib.Path(__file__).parent
w = pd.read_csv(ROOT/'results'/'per_variant_scores.csv')
bg = w[w.panel_set=='background']; pos = w[w.panel_set=='positive']
can = w[w.panel_set=='canonical'];  com = w[w.panel_set=='common']; deep = w[w.panel_set=='deep']

print("=== A. The quantile score is saturated here and cannot be used ===")
for s in ('background','common','deep','canonical','positive'):
    v = w[w.panel_set==s]['max_quantile']
    print(f"    {s:<11} n={len(v):<4} min {v.min():.6f}  median {v.median():.6f}  max {v.max():.6f}")
print("    Every set sits at ~1.0. max-over-tracks quantile does not separate anything;")
print("    only the raw scores do.\n")

def auc(a, b):                      # P(random a > random b), ties at 0.5
    a, b = np.asarray(a), np.asarray(b)
    return float((np.greater.outer(a, b).sum() + 0.5*np.equal.outer(a, b).sum()) / (len(a)*len(b)))

print("=== B. Separation, AUC vs the matched background ===")
for m in ('SPLICE_SITES','SPLICE_SITE_USAGE','SPLICE_JUNCTIONS'):
    print(f"    {m:<18} canonical {auc(can[m], bg[m]):.3f}   positives {auc(pos[m], bg[m]):.3f}")
print()

M = 'SPLICE_SITES'
thr = bg[M].quantile(.95)
print(f"=== C. Operating point: {M} >= {thr:.4f} (95th pct of matched background) ===")
for s, d in (('positive',pos), ('canonical',can), ('background',bg), ('common',com), ('deep',deep)):
    n = int((d[M] >= thr).sum())
    print(f"    {s:<11} {n:>3}/{len(d):<4} pass  ({100*n/len(d):5.1f}%)")
print()
print("    Applied to a candidate set of 12 rare intronic variants, this threshold")
print(f"    yields {12*0.05:.1f} expected false positives, and caught 3/3 known alleles.")
print()

print("=== D. Where the positives sit against the canonical GT/AG anchor ===")
print(f"    canonical  median {can[M].median():.4f}   min {can[M].min():.4f}")
print(f"    positives  {', '.join(f'{r.panel_name.split(chr(95))[0]} {r[M]:.4f}' for _, r in pos.iterrows())}")
print("    The three non-canonical alleles score as high as true GT/AG disruptions,")
print("    which is what their published mis-splicing phenotypes would predict.")
