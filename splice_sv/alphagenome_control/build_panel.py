"""Build the variant panel for the AlphaGenome splice positive control.

Every variant here is public: three published MVA1 splice alleles from the
literature, and gnomAD v4 SNVs in BUB1B. No patient data is read or written.

Sets
  positive    3 literature alleles with experimentally reported mis-splicing
  canonical   gnomAD SNVs at intron offset +/-1,2 (the GT..AG dinucleotides).
              Expected-high anchor: if these do not score high the model is
              not reading splice sites at this locus and the test is void.
  background  gnomAD intronic SNVs at |offset| 3-20, positives excluded.
              The matched null -- same positional class as the positives.
  common      gnomAD intronic SNVs at AF >= 1e-3, any offset. Expected-low.
  deep        gnomAD deep intronic SNVs, |offset| > 100. Floor check.
"""
import json, pathlib

D = pathlib.Path(__file__).parent / 'data'
snvs = json.load(open(D / 'bub1b_snvs.json'))

# GRCh38, ENST00000287598.11. Coordinates resolved via Ensembl variant_recoder
# and independently confirmed against the Ensembl exon table (see verify below).
POSITIVES = [
    dict(name='RioFrio_c.2386-11A>G', pos=40212488, ref='A', alt='G',
         cite='Rio Frio 2010 NEJM 10.1056/NEJMoa1006565',
         effect='de novo acceptor preferred over the authentic site'),
    dict(name='Kawakami_c.1402-5A>G', pos=40200239, ref='A', alt='G',
         cite='Kawakami 2023 J Hum Genet 10.1038/s10038-023-01211-8',
         effect='mis-splicing with 4-bp insertion; BUBR1 severely diminished'),
    dict(name='Laberko_c.1288+5G>A', pos=40196779, ref='G', alt='A',
         cite='Laberko 2019 Orphanet J Rare Dis 10.1186/s13023-019-1073-x',
         effect='donor +5 in a compound het with c.498_505del'),
]

exons = json.load(open(D / 'exons.json'))

def offset(pos):
    for e in exons:
        if e['start'] <= pos <= e['end']:
            return ('exonic', 0)
    best = None
    for e in exons:
        if pos > e['end']:
            n = pos - e['end']
            if best is None or n < abs(best[1]):
                best = ('donor', n)
        if pos < e['start']:
            n = e['start'] - pos
            if best is None or n < abs(best[1]):
                best = ('acceptor', -n)
    return best

# Confirm the recoder coordinates land where the papers say they do.
for p in POSITIVES:
    p['kind'], p['off'] = offset(p['pos'])
    p['set'] = 'positive'
expect = {'RioFrio_c.2386-11A>G': ('acceptor', -11),
          'Kawakami_c.1402-5A>G': ('acceptor', -5),
          'Laberko_c.1288+5G>A': ('donor', 5)}
for p in POSITIVES:
    got = (p['kind'], p['off'])
    assert got == expect[p['name']], f"{p['name']}: exon table says {got}, HGVS says {expect[p['name']]}"

pos_keys = {(p['pos'], p['ref'], p['alt']) for p in POSITIVES}
intronic = [r for r in snvs if r['kind'] != 'exonic']

def tag(r, s):
    r = dict(r); r['set'] = s; r['name'] = r['vid']; return r

panel = list(POSITIVES)
panel += [tag(r, 'canonical')  for r in intronic if 1 <= abs(r['off']) <= 2]
panel += [tag(r, 'background') for r in intronic if 3 <= abs(r['off']) <= 20
          and (r['pos'], r['ref'], r['alt']) not in pos_keys]
panel += [tag(r, 'common')     for r in intronic if r['af'] >= 1e-3
          and not (3 <= abs(r['off']) <= 20) and abs(r['off']) > 2]
panel += [tag(r, 'deep')       for r in intronic if abs(r['off']) > 100][:60]

seen, dedup = set(), []
for r in panel:
    k = (r['pos'], r['ref'], r['alt'])
    if k in seen:
        continue
    seen.add(k); dedup.append(r)

json.dump(dedup, open(D / 'panel.json', 'w'), indent=1)

from collections import Counter
print("panel:", len(dedup), dict(Counter(r['set'] for r in dedup)))
print()
print("positives, and whether gnomAD v4 holds them:")
by_key = {(r['pos'], r['ref'], r['alt']): r for r in snvs}
for p in POSITIVES:
    g = by_key.get((p['pos'], p['ref'], p['alt']))
    print(f"  {p['name']:<22} chr15:{p['pos']} {p['ref']}>{p['alt']}  {p['kind']}{p['off']:+d}  "
          f"gnomAD AF={g['af']:.3g}" if g else
          f"  {p['name']:<22} chr15:{p['pos']} {p['ref']}>{p['alt']}  {p['kind']}{p['off']:+d}  absent from gnomAD v4")
