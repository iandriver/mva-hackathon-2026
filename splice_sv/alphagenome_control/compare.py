"""AlphaGenome vs SpliceAI on the same panel, same background, same statistic."""
import pathlib, re
import numpy as np, pandas as pd

ROOT = pathlib.Path(__file__).parent
ag = pd.read_csv(ROOT/'results'/'per_variant_scores.csv')

rows = []
for line in open(ROOT/'results'/'spliceai_D50.vcf'):
    if line.startswith('#'): continue
    f = line.rstrip('\n').split('\t')
    m = re.search(r'SpliceAI=([^;\s]+)', f[7])
    if not m: continue
    p = m.group(1).split('|')
    ds = [float(x) for x in p[2:6]]                    # DS_AG, DS_AL, DS_DG, DS_DL
    rows.append(dict(panel_name=f[2], spliceai=max(ds),
                     DS_AG=ds[0], DS_AL=ds[1], DS_DG=ds[2], DS_DL=ds[3]))
sa = pd.DataFrame(rows)
print(f"SpliceAI scored {len(sa)} / {len(ag)} panel variants")

w = ag.merge(sa, on='panel_name', how='inner')
w.to_csv(ROOT/'results'/'compare_ag_spliceai.csv', index=False)
bg  = w[w.panel_set=='background']; pos = w[w.panel_set=='positive']
can = w[w.panel_set=='canonical'];  com = w[w.panel_set=='common']; deep = w[w.panel_set=='deep']

def auc(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return float((np.greater.outer(a,b).sum() + 0.5*np.equal.outer(a,b).sum())/(len(a)*len(b)))

AG, SA = 'SPLICE_SITES', 'spliceai'

print(f"\n=== Ranks of the three positives in the matched background (n={len(bg)}) ===\n")
print(f"{'allele':<22} {'offset':>6}  {'AlphaGenome':>13}  {'SpliceAI':>13}  {'proximity':>11}")
for _, p in pos.iterrows():
    r_ag = int((bg[AG] >= p[AG]).sum())+1
    r_sa = int((bg[SA] >= p[SA]).sum())+1
    r_px = int((bg.panel_off.abs() <= abs(p.panel_off)).sum())+1
    print(f"{p.panel_name:<22} {int(p.panel_off):>+6}  {r_ag:>5}/{len(bg)+1}   "
          f"{r_sa:>5}/{len(bg)+1}   {r_px:>4}/{len(bg)+1}")
print()
print(f"{'':<22} {'':>6}  {'score':>13}  {'score':>13}")
for _, p in pos.iterrows():
    print(f"{p.panel_name:<22} {'':>6}  {p[AG]:>13.4f}  {p[SA]:>13.4f}")

print("\n=== Separation from the matched background (AUC) ===\n")
print(f"{'':<12} {'AlphaGenome':>12} {'SpliceAI':>10}")
for nm, d in (('positives',pos), ('canonical',can)):
    print(f"{nm:<12} {auc(d[AG],bg[AG]):>12.3f} {auc(d[SA],bg[SA]):>10.3f}")
print(f"{'common':<12} {auc(com[AG],bg[AG]):>12.3f} {auc(com[SA],bg[SA]):>10.3f}   (want low)")
print(f"{'deep':<12} {auc(deep[AG],bg[AG]):>12.3f} {auc(deep[SA],bg[SA]):>10.3f}   (want low)")

print("\n=== Operating point at each tool's own background 95th percentile ===\n")
for nm, col in (('AlphaGenome',AG), ('SpliceAI',SA)):
    thr = bg[col].quantile(.95)
    print(f"  {nm}  threshold {thr:.4f}")
    for s, d in (('positive',pos),('canonical',can),('background',bg),('common',com),('deep',deep)):
        n = int((d[col] >= thr).sum())
        print(f"      {s:<11} {n:>3}/{len(d):<4} ({100*n/len(d):5.1f}%)")
    print()

print("=== Agreement between the two models across the whole panel ===\n")
from scipy.stats import spearmanr, pearsonr
r_s,_ = spearmanr(w[AG], w[SA]); r_p,_ = pearsonr(w[AG], w[SA])
print(f"  Spearman {r_s:.3f}   Pearson {r_p:.3f}   (n={len(w)})")
thr_ag, thr_sa = bg[AG].quantile(.95), bg[SA].quantile(.95)
a = w[AG] >= thr_ag; b = w[SA] >= thr_sa
print(f"  Both flag {int((a&b).sum())}   AG only {int((a&~b).sum())}   "
      f"SpliceAI only {int((b&~a).sum())}   neither {int((~a&~b).sum())}")

print("\n=== SpliceAI's called direction on the positives ===\n")
for _, p in pos.iterrows():
    print(f"  {p.panel_name:<22} AG_gain {p.DS_AG:.2f}  AG_loss {p.DS_AL:.2f}  "
          f"DG_gain {p.DS_DG:.2f}  DG_loss {p.DS_DL:.2f}")
