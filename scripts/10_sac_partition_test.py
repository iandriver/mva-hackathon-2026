#!/usr/bin/env python3
"""Prospective test of the post-hoc split in findings/INTERFERON_HYPOTHESIS.md check 2.

THE HYPOTHESIS, generated post-hoc and therefore not evidence for itself:
  interferon activation after knockout tracks with SPINDLE-CHECKPOINT SIGNALLING genes
  (loss abrogates the checkpoint, cells proceed through mitosis missegregating, producing
  viable aneuploid daughters that form micronuclei) and NOT with OUTER-KINETOCHORE
  ATTACHMENT genes (loss prevents stable attachment, holds the checkpoint on, and arrests
  cells before they can produce aneuploid daughters).

THE TEST:
  1. The two groups below are defined from prior annotation of protein function.
     They are fixed in this file BEFORE any score is computed.
  2. Every gene whose score was seen while generating the hypothesis is EXCLUDED.
     The script aborts if any leaks in.
  3. Prediction registered in advance: GROUP A mean > GROUP B mean.
  4. Significance by label permutation, which makes no distributional assumption.

A result that fails this test retires the split. That is the point of running it.
"""
import bisect, gzip, json, os, random, statistics as st, sys, urllib.request
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# ---- pre-registered, from annotation, not from data -------------------------------
GROUP_A_SAC_SIGNALLING = [       # loss abrogates the checkpoint
    "MAD1L1", "MAD2L1", "BUB3", "TRIP13", "ZW10", "ZWILCH", "KNTC1", "CASC5", "KNL1",
    "CENPF", "CCDC99", "SPDL1",
]
GROUP_B_ATTACHMENT = [           # loss prevents attachment and arrests cells
    "SPC24", "DSN1", "NSL1", "MIS12", "PMF1", "SKA1", "SKA2", "SKA3", "NDC80L",
]
PREDICTION = "GROUP_A > GROUP_B"

# every gene whose ISG score was inspected while forming the hypothesis
GENERATING_SET = {"BUB1B", "BUB1", "TTK", "CENPE", "AURKB", "PLK1", "KIF11",
                  "NDC80", "SPC25", "NUF2", "ZWINT", "INCENP"}
# -----------------------------------------------------------------------------------

META = "https://maayanlab.cloud/sigcom-lincs/metadata-api/"
CONSENSUS = ("https://s3.amazonaws.com/lincs-dcic/LINCS-sigs-2021/means/xpr_mean_sigs/"
             "LINCS_L1000_CRISPRKO_{g}_mean_sig.tsv.gz")
LIB = "98c6a1ef-e3b7-5a96-9913-480acf78a577"
ENRICHR = ("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&"
           "libraryName=MSigDB_Hallmark_2020")
OUT = os.environ.get("MVA_RESULTS", os.path.join(os.path.dirname(__file__), "..", "results"))
CACHE = os.environ.get("MVA_CACHE", "/Volumes/IanSSD/mva-hackathon/consensus")
os.makedirs(OUT, exist_ok=True); os.makedirs(CACHE, exist_ok=True)
N_BACKGROUND, SEED, N_PERM = 300, 11, 20000

def hallmark_ifn():
    sets = {}
    for l in urllib.request.urlopen(ENRICHR, timeout=120).read().decode().splitlines():
        p = l.rstrip("\n").split("\t"); sets[p[0]] = [g for g in p[2:] if g]
    return set(sets["Interferon Alpha Response"]) | set(sets["Interferon Gamma Response"])

def available_genes():
    out, skip = [], 0
    while True:
        req = urllib.request.Request(META + "signatures/find",
            data=json.dumps({"filter": {"where": {"library": LIB},
                                        "fields": ["meta.pert_name"],
                                        "limit": 1000, "skip": skip}}).encode(),
            headers={"Content-Type": "application/json"})
        b = json.load(urllib.request.urlopen(req, timeout=180))
        if not b: break
        out += [x["meta"].get("pert_name") for x in b]; skip += 1000
        if skip > 12000: break
    return {g for g in out if g}

def score(gene, isg):
    p = os.path.join(CACHE, f"LINCS_L1000_CRISPRKO_{gene}_mean_sig.tsv.gz")
    if not os.path.exists(p):
        try: urllib.request.urlretrieve(CONSENSUS.format(g=gene), p)
        except Exception: return None
    try:
        with gzip.open(p, "rt") as f:
            f.readline(); vals = {}
            for l in f:
                a = l.rstrip("\n").split("\t")
                if len(a) > 1 and a[1]: vals[a[0]] = float(a[1])
    except Exception: return None
    if len(vals) < 5000: return None
    srt = sorted(vals.items(), key=lambda kv: kv[1]); n = len(srt)
    pr = {g: i / (n - 1) for i, (g, _) in enumerate(srt)}
    hit = [pr[g] for g in isg if g in pr]
    return st.mean(hit) if len(hit) >= 50 else None

def main():
    print(__doc__.split("A result that fails")[0].strip())
    print("=" * 74)
    print(f"PRE-REGISTERED PREDICTION: {PREDICTION}")
    print(f"GROUP A (checkpoint signalling): {', '.join(GROUP_A_SAC_SIGNALLING)}")
    print(f"GROUP B (attachment structural): {', '.join(GROUP_B_ATTACHMENT)}")
    print(f"EXCLUDED (used to form the hypothesis): {', '.join(sorted(GENERATING_SET))}")
    print("=" * 74)

    leak = (set(GROUP_A_SAC_SIGNALLING) | set(GROUP_B_ATTACHMENT)) & GENERATING_SET
    if leak:
        sys.exit(f"ABORT: {leak} appear in both the test set and the generating set. "
                 f"The test would not be prospective.")

    isg = hallmark_ifn()
    have = available_genes()
    A = [g for g in GROUP_A_SAC_SIGNALLING if g in have]
    B = [g for g in GROUP_B_ATTACHMENT if g in have]
    print(f"\nwith a consensus signature:  A {len(A)}/{len(GROUP_A_SAC_SIGNALLING)} "
          f"({', '.join(A)})   B {len(B)}/{len(GROUP_B_ATTACHMENT)} ({', '.join(B)})")
    if len(A) < 3 or len(B) < 3:
        json.dump({"prediction": PREDICTION, "verdict": "UNTESTABLE",
                   "reason": "LINCS consensus CRISPR-KO coverage",
                   "group_a_available": A, "group_b_available": B,
                   "group_a_requested": GROUP_A_SAC_SIGNALLING,
                   "group_b_requested": GROUP_B_ATTACHMENT,
                   "library_size": len(have)},
                  open(os.path.join(OUT, "10_sac_partition_test.json"), "w"), indent=2)
        print(f"\nUNTESTABLE: the consensus CRISPR-KO library covers {len(have):,} genes but "
              f"contains {len(A)} of {len(GROUP_A_SAC_SIGNALLING)} group-A and {len(B)} of "
              f"{len(GROUP_B_ATTACHMENT)} group-B genes.\n"
              f"The twelve genes that generated the hypothesis are almost the only "
              f"checkpoint genes the library holds, so the split cannot be tested on genes "
              f"held out from it. It stays post-hoc.\n"
              f"wrote {OUT}/10_sac_partition_test.json")
        sys.exit(3)

    random.seed(SEED)
    bg_pool = sorted(have - GENERATING_SET - set(A) - set(B))
    background = random.sample(bg_pool, N_BACKGROUND)
    targets = A + B + background
    with ThreadPoolExecutor(max_workers=12) as ex:
        sc = dict(zip(targets, ex.map(lambda g: score(g, isg), targets)))

    bg = sorted(v for g in background if (v := sc.get(g)) is not None)
    mean, sd = st.mean(bg), st.pstdev(bg)
    pct = lambda v: bisect.bisect_left(bg, v) / len(bg) * 100
    print(f"\nbackground n={len(bg)}  mean={mean:.4f}  sd={sd:.4f}")

    rows = []
    print(f"\n{'gene':<9}{'group':<6}{'score':>9}{'pct':>8}{'z':>8}")
    for grp, genes in (("A", A), ("B", B)):
        for g in genes:
            v = sc.get(g)
            if v is None:
                print(f"{g:<9}{grp:<6}{'no signature':>25}"); continue
            rows.append({"gene": g, "group": grp, "score": v,
                         "percentile": pct(v), "z": (v - mean) / sd})
            print(f"{g:<9}{grp:<6}{v:>9.4f}{pct(v):>8.1f}{(v-mean)/sd:>+8.2f}")

    a = [r["score"] for r in rows if r["group"] == "A"]
    b = [r["score"] for r in rows if r["group"] == "B"]
    if len(a) < 3 or len(b) < 3:
        sys.exit(f"ABORT: only {len(a)} A and {len(b)} B genes scored.")
    diff = st.mean(a) - st.mean(b)

    pool = a + b; na = len(a)
    hits = 0
    for _ in range(N_PERM):
        random.shuffle(pool)
        if st.mean(pool[:na]) - st.mean(pool[na:]) >= diff: hits += 1
    p = (hits + 1) / (N_PERM + 1)

    print(f"\n{'':<20}{'n':>4}{'mean':>10}{'median':>10}")
    print(f"{'A signalling':<20}{len(a):>4}{st.mean(a):>10.4f}{st.median(a):>10.4f}")
    print(f"{'B attachment':<20}{len(b):>4}{st.mean(b):>10.4f}{st.median(b):>10.4f}")
    print(f"\ndifference A - B = {diff:+.4f}")
    print(f"one-sided permutation p = {p:.4f}  ({N_PERM:,} permutations)")

    held = diff > 0 and p < 0.05
    print("\n" + "=" * 74)
    if held:
        print("PREDICTION HELD on genes not used to form it. The split survives a "
              "prospective test.")
    elif diff > 0:
        print(f"PREDICTION DIRECTIONALLY CORRECT BUT NOT SIGNIFICANT (p = {p:.3f}). "
              f"With n={len(a)} vs {len(b)} this is underpowered; the split is not "
              f"established.")
    else:
        print("PREDICTION FAILED. The direction reverses on held-out genes. The "
              "signalling-vs-attachment split should be retired, not reported.")
    print("=" * 74)

    json.dump({"prediction": PREDICTION, "group_a": A, "group_b": B,
               "excluded": sorted(GENERATING_SET), "background_n": len(bg),
               "background_mean": mean, "background_sd": sd, "rows": rows,
               "mean_a": st.mean(a), "mean_b": st.mean(b), "difference": diff,
               "permutation_p": p, "n_permutations": N_PERM, "seed": SEED,
               "verdict": "held" if held else ("directional" if diff > 0 else "failed")},
              open(os.path.join(OUT, "10_sac_partition_test.json"), "w"), indent=2)
    print(f"wrote {OUT}/10_sac_partition_test.json")

if __name__ == "__main__":
    main()
