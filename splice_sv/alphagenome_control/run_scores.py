"""Score the panel with AlphaGenome's three splice heads.

Writes one parquet shard per chunk so an interrupted run resumes rather than
restarts. Public variants only.
"""
import json, os, pathlib, sys, time
import pandas as pd
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

ROOT = pathlib.Path(__file__).parent
OUT = ROOT / 'results' / 'shards'
OUT.mkdir(parents=True, exist_ok=True)

SEQLEN = dna_client.SEQUENCE_LENGTH_100KB
SCORERS = [variant_scorers.RECOMMENDED_VARIANT_SCORERS[k]
           for k in ('SPLICE_SITES', 'SPLICE_SITE_USAGE', 'SPLICE_JUNCTIONS')]
CHUNK = 10

panel = json.load(open(ROOT / 'data' / 'panel.json'))
model = dna_client.create(os.environ['ALPHAGENOME_API_KEY'])

chunks = [panel[i:i + CHUNK] for i in range(0, len(panel), CHUNK)]
print(f"{len(panel)} variants in {len(chunks)} chunks of {CHUNK}", flush=True)

for ci, chunk in enumerate(chunks):
    shard = OUT / f'shard_{ci:04d}.parquet'
    if shard.exists():
        continue
    variants, intervals = [], []
    for r in chunk:
        v = genome.Variant(chromosome='chr15', position=r['pos'],
                           reference_bases=r['ref'], alternate_bases=r['alt'])
        variants.append(v)
        intervals.append(v.reference_interval.resize(SEQLEN))

    for attempt in range(5):
        try:
            scores = model.score_variants(intervals=intervals, variants=variants,
                                          variant_scorers=SCORERS,
                                          progress_bar=False, max_workers=4)
            break
        except Exception as e:                       # rate limit / transient RPC
            wait = 20 * (attempt + 1)
            print(f"  chunk {ci} attempt {attempt} failed: {type(e).__name__}: "
                  f"{str(e)[:180]} -- retry in {wait}s", flush=True)
            time.sleep(wait)
    else:
        print(f"  chunk {ci} ABANDONED", flush=True)
        continue

    frames = []
    for r, per_variant in zip(chunk, scores):
        df = variant_scorers.tidy_scores(per_variant)
        df = df[df['gene_name'] == 'BUB1B']          # ignore neighbouring genes
        df['panel_name'] = r['name']
        df['panel_set'] = r['set']
        df['panel_off'] = r['off']
        df['panel_kind'] = r['kind']
        df['panel_af'] = r.get('af')
        frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    for c in ('variant_id', 'scored_interval', 'variant_scorer'):
        out[c] = out[c].astype(str)          # these hold objects, not strings
    out.to_parquet(shard)
    print(f"  shard {ci+1}/{len(chunks)} written", flush=True)

print("DONE", flush=True)
