#!/usr/bin/env python3
"""Reproduces findings/INTERFERON_HYPOTHESIS.md checks 1 and 2.

The BUB1B-KO consensus signature's top genes are an interferon module. Two questions
decide whether that means anything:

  CHECK 1 (specificity)  Is it specific to BUB1B, or does any CRISPR knockout do this?
                         Scored against 300 randomly sampled consensus KO signatures.
  CHECK 2 (mechanism)    Do other spindle-checkpoint knockouts behave the same way?

Metric: mean percentile-rank of Hallmark interferon genes within each signature.
0.5 means no shift, so the background mean is itself a check that the metric is unbiased.
"""
import bisect, gzip, io, json, os, random, statistics as st, sys, urllib.request
from concurrent.futures import ThreadPoolExecutor

META = "https://maayanlab.cloud/sigcom-lincs/metadata-api/"
CONSENSUS = ("https://s3.amazonaws.com/lincs-dcic/LINCS-sigs-2021/means/xpr_mean_sigs/"
             "LINCS_L1000_CRISPRKO_{g}_mean_sig.tsv.gz")
LIB_MEAN_XPR = "98c6a1ef-e3b7-5a96-9913-480acf78a577"
ENRICHR = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020"
OUT = os.environ.get("MVA_RESULTS", os.path.join(os.path.dirname(__file__), "..", "results"))
CACHE = os.environ.get("MVA_CACHE", "/Volumes/IanSSD/mva-hackathon/consensus")
os.makedirs(OUT, exist_ok=True); os.makedirs(CACHE, exist_ok=True)

TARGET = "BUB1B"
SAC = ["BUB1", "TTK", "CENPE", "AURKB", "PLK1", "KIF11", "NDC80", "SPC25", "NUF2", "ZWINT", "INCENP"]
N_BACKGROUND = 300
SEED = 11
MIN_GENES, MIN_ISG_HITS = 5000, 50

def hallmark_ifn():
    sets = {}
    for l in urllib.request.urlopen(ENRICHR, timeout=120).read().decode().splitlines():
        p = l.rstrip("\n").split("\t")
        sets[p[0]] = [g for g in p[2:] if g]
    return set(sets["Interferon Alpha Response"]) | set(sets["Interferon Gamma Response"])

def all_consensus_genes():
    genes, skip = [], 0
    while True:
        req = urllib.request.Request(META + "signatures/find",
            data=json.dumps({"filter": {"where": {"library": LIB_MEAN_XPR},
                                        "fields": ["meta.pert_name"],
                                        "limit": 1000, "skip": skip}}).encode(),
            headers={"Content-Type": "application/json"})
        b = json.load(urllib.request.urlopen(req, timeout=180))
        if not b: break
        genes += [x["meta"].get("pert_name") for x in b]
        skip += 1000
        if skip > 12000: break
    return sorted({g for g in genes if g})

def fetch(gene):
    p = os.path.join(CACHE, f"LINCS_L1000_CRISPRKO_{gene}_mean_sig.tsv.gz")
    if not os.path.exists(p):
        try:
            urllib.request.urlretrieve(CONSENSUS.format(g=gene), p)
        except Exception:
            return None
    return p

def score(gene, isg):
    p = fetch(gene)
    if not p: return None
    try:
        with gzip.open(p, "rt") as f:
            f.readline()
            vals = {}
            for l in f:
                a = l.rstrip("\n").split("\t")
                if len(a) > 1 and a[1]: vals[a[0]] = float(a[1])
    except Exception:
        return None
    if len(vals) < MIN_GENES: return None
    srt = sorted(vals.items(), key=lambda kv: kv[1])
    n = len(srt)
    pr = {g: i / (n - 1) for i, (g, _) in enumerate(srt)}
    hit = [pr[g] for g in isg if g in pr]
    return st.mean(hit) if len(hit) >= MIN_ISG_HITS else None

def main():
    isg = hallmark_ifn()
    print(f"Hallmark interferon alpha+gamma union: {len(isg)} genes", file=sys.stderr)
    genes = all_consensus_genes()
    print(f"consensus KO signatures available: {len(genes):,}", file=sys.stderr)
    random.seed(SEED)
    pool = [g for g in genes if g != TARGET and g not in SAC]
    background = random.sample(pool, N_BACKGROUND)

    want = [TARGET] + [g for g in SAC if g in genes] + background
    with ThreadPoolExecutor(max_workers=12) as ex:
        scores = dict(zip(want, ex.map(lambda g: score(g, isg), want)))

    bg = sorted(v for g, v in scores.items() if g in background and v is not None)
    mean, sd = st.mean(bg), st.pstdev(bg)
    print(f"\nbackground: n={len(bg)}  mean={mean:.4f}  sd={sd:.4f}   (0.5 = no shift)")

    def report(g):
        v = scores.get(g)
        if v is None: return None
        return {"gene": g, "score": v,
                "percentile": bisect.bisect_left(bg, v) / len(bg) * 100,
                "z": (v - mean) / sd}

    t = report(TARGET)
    print(f"\nCHECK 1 -- specificity")
    print(f"  {TARGET}: {t['score']:.4f}  {t['percentile']:.1f}th percentile  z={t['z']:+.2f}")
    print(f"  {'PASS' if t['z'] > 1.5 else 'FAIL'}: interferon activation is "
          f"{'specific to this knockout' if t['z'] > 1.5 else 'not distinguishable from generic knockout stress'}")

    print(f"\nCHECK 2 -- do other spindle-checkpoint knockouts agree?")
    print(f"  {'gene':<9}{'score':>9}{'pct':>8}{'z':>8}")
    sac_rows = []
    for g in [TARGET] + SAC:
        r = report(g)
        if r:
            sac_rows.append(r)
            print(f"  {g:<9}{r['score']:>9.4f}{r['percentile']:>8.1f}{r['z']:>+8.2f}")

    json.dump({"background": {"n": len(bg), "mean": mean, "sd": sd},
               "target": t, "sac_panel": sac_rows, "seed": SEED},
              open(os.path.join(OUT, "08_interferon_specificity.json"), "w"), indent=2)
    print(f"\nwrote {OUT}/08_interferon_specificity.json", file=sys.stderr)
    print("\nNOTE: any partition of the SAC panel (checkpoint signalling vs attachment)\n"
          "read off this table is POST-HOC. It was not tested prospectively.", file=sys.stderr)

if __name__ == "__main__":
    main()
