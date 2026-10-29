"""Does AlphaGenome separate known pathogenic MVA1 splice alleles from their
positionally matched neighbours in BUB1B?

Per variant, each splice head is reduced to max |raw_score| over its tracks.
Positives are ranked inside the matched background (intronic, |offset| 3-20),
which is the comparison a user would actually make: this candidate against its
positional peers. p is the exact one-sided rank p-value, (rank)/(n_bg+1).
"""
import pathlib
import numpy as np, pandas as pd

ROOT = pathlib.Path(__file__).parent
df = pd.concat([pd.read_parquet(p) for p in sorted((ROOT/'results'/'shards').glob('*.parquet'))],
               ignore_index=True)
print(f"rows {len(df):,}   variants {df.panel_name.nunique()}   "
      f"outputs {sorted(df.output_type.unique())}\n")

df['abs_raw'] = df.raw_score.abs()
# The literature positives carry no gnomAD AF in the panel; pandas groupby drops
# NaN keys by default and would silently delete exactly the rows under test.
df['panel_af'] = df['panel_af'].fillna(-1.0)
per = (df.groupby(['panel_name','panel_set','panel_off','panel_kind','panel_af','output_type'])
         .agg(abs_raw=('abs_raw','max'), quantile=('quantile_score','max'))
         .reset_index())
wide = per.pivot_table(index=['panel_name','panel_set','panel_off','panel_kind','panel_af'],
                       columns='output_type', values='abs_raw').reset_index()
q = per.groupby(['panel_name'])['quantile'].max().rename('max_quantile')
wide = wide.merge(q, on='panel_name')
METRICS = ['SPLICE_SITES','SPLICE_SITE_USAGE','SPLICE_JUNCTIONS','max_quantile']
wide.to_csv(ROOT/'results'/'per_variant_scores.csv', index=False)

print("=== 1. Validity anchor: does the model read splice sites here at all? ===")
print("    median max|raw| by set\n")
print(wide.groupby('panel_set')[METRICS].median().round(4).to_string(), "\n")
can = wide[wide.panel_set=='canonical']; bg = wide[wide.panel_set=='background']
for m in METRICS[:3]:
    hi = (can[m].median() > bg[m].median())
    print(f"    {m:<18} canonical {can[m].median():.4f}  vs background {bg[m].median():.4f}   "
          f"{'PASS' if hi else 'FAIL'}")
print()

print("=== 2. The test: positives ranked within the matched background ===")
print("    background = intronic SNVs at |offset| 3-20, n =", len(bg), "\n")
pos = wide[wide.panel_set=='positive']
rows = []
for _, p in pos.iterrows():
    for m in METRICS:
        n_ge = int((bg[m] >= p[m]).sum())
        rank = n_ge + 1
        rows.append(dict(variant=p.panel_name, metric=m, score=p[m],
                         rank=rank, n_bg=len(bg), pct=100*rank/(len(bg)+1),
                         p=rank/(len(bg)+1),
                         bg_p95=bg[m].quantile(.95), bg_median=bg[m].median()))
res = pd.DataFrame(rows)
res.to_csv(ROOT/'results'/'positive_ranks.csv', index=False)
for v in res.variant.unique():
    s = res[res.variant==v]
    off = int(pos[pos.panel_name==v].panel_off.iloc[0])
    print(f"  {v}  (offset {off:+d})")
    for _, r in s.iterrows():
        star = '  <-- top 5%' if r.pct <= 5 else ''
        print(f"      {r.metric:<18} score {r.score:8.4f}   rank {int(r['rank']):>4}/{r.n_bg+1}"
              f"   top {r.pct:5.2f}%   p={r.p:.4f}{star}")
    print()

print("=== 3. Detection floor: what a variant must reach to clear the background ===")
print("    (95th and 99th percentile of the matched background)\n")
for m in METRICS:
    print(f"    {m:<18} p95 {bg[m].quantile(.95):8.4f}    p99 {bg[m].quantile(.99):8.4f}"
          f"    max {bg[m].max():8.4f}")
print()

print("=== 4. Joint test: all three positives, one-sided ===")
for m in METRICS:
    ranks = res[res.metric==m]['rank'].values
    # P(all three land at rank <= max observed) under the null of random ranks
    worst = ranks.max(); N = len(bg)+1
    print(f"    {m:<18} ranks {sorted(ranks)}   worst top {100*worst/N:.2f}%   "
          f"P(all 3 <= that by chance) = {(worst/N)**3:.2e}")
