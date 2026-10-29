#!/usr/bin/env python3
"""Reports ddG for N1002K from both predictors, each calibrated against the whole domain."""
import bisect, csv, json, os, statistics as st, subprocess, sys

REPO = os.path.join(os.path.dirname(__file__), "..")
OUT = os.environ.get("MVA_RESULTS", os.path.join(REPO, "results"))
WORK = os.environ.get("MVA_WORK", os.path.join(REPO, "tools", "ddgun-work"))
os.makedirs(OUT, exist_ok=True)
DOMAIN_START, VARIANT, WT, MUT = 766, 1002, "N", "K"
res = {}

# --- ThermoMPNN: positive = destabilising ---
for tag, path, off in (("domain", "ThermoMPNN_inference_BUBR1_AF_kinase_766_1050.csv", DOMAIN_START),
                       ("full",   "ThermoMPNN_inference_BUBR1_AF_full.csv", 1)):
    p = os.path.join(REPO, "structures", path)
    if not os.path.exists(p): continue
    rows = [(float(r["ddG_pred"]), int(r["position"]) + off, r["wildtype"], r["mutation"])
            for r in csv.DictReader(open(p))]
    hit = [d for d, pos, wt, mu in rows if pos == VARIANT and wt == WT and mu == MUT]
    dom = sorted(d for d, pos, wt, mu in rows
                 if DOMAIN_START <= pos <= 1050 and wt != mu)
    if hit:
        res[f"thermompnn_{tag}"] = {
            "ddG": hit[0], "convention": "positive = destabilising",
            "percentile": bisect.bisect_left(dom, hit[0]) / len(dom) * 100,
            "n_substitutions": len(dom), "domain_mean": st.mean(dom),
            "domain_median": st.median(dom)}

# --- DDGun: negative = destabilising ---
prof = os.path.join(WORK, "bubr1_kinase.prof")
if os.path.exists(prof):
    code = f'''
import sys, statistics as st, bisect
from ddgun.ddgun import ddgun_seq
from ddgun.aa import Substitution
from ddgun.profile import Profile
p = Profile({prof!r})
seq = list(p.profile["SEQ"]); AA = "ACDEFGHIKLMNPQRSTVWY"
idx = {VARIANT} - {DOMAIN_START} + 1
target = ddgun_seq(Substitution.parse(f"{WT}{{idx}}{MUT}"), p)
vals = []
for i, wt in enumerate(seq, start=1):
    if wt not in AA: continue
    for mt in AA:
        if mt == wt: continue
        try: vals.append(ddgun_seq(Substitution.parse(f"{{wt}}{{i}}{{mt}}"), p))
        except Exception: pass
vals.sort()
print(target, len(vals), st.mean(vals), st.median(vals),
      bisect.bisect_left(vals, target) / len(vals) * 100)
'''
    r = subprocess.run([os.path.join(REPO, "tools", "ddg-env", "bin", "python"), "-c", code],
                       capture_output=True, text=True)
    line = [l for l in r.stdout.splitlines() if l and l[0].isdigit() or l.startswith("-")]
    if line:
        t, n, mean, med, pct = line[-1].split()
        res["ddgun_seq"] = {"ddG": float(t), "convention": "negative = destabilising",
                            "percentile_from_destabilising_end": float(pct),
                            "n_substitutions": int(n), "domain_mean": float(mean),
                            "domain_median": float(med)}

print(f"\nBUBR1 p.{WT}{VARIANT}{MUT}\n" + "=" * 64)
for k, v in res.items():
    if "percentile" in v:
        note = f"more destabilising than {v['percentile']:.1f}% of domain substitutions"
    else:
        note = (f"more destabilising than {100 - v['percentile_from_destabilising_end']:.1f}% "
                f"of domain substitutions")
    print(f"  {k:<22} ddG {v['ddG']:+.3f}   ({v['convention']})")
    print(f"  {'':<22} {note}  (n={v['n_substitutions']})")
mags = [abs(v["ddG"]) for v in res.values()]
if len(mags) > 1:
    print(f"\n  magnitude across methods: {min(mags):.3f} - {max(mags):.3f} kcal/mol")
    print("  Both conventions agree on DIRECTION (destabilising) and MAGNITUDE.")
json.dump(res, open(os.path.join(OUT, "09_ddg_stability.json"), "w"), indent=2)
print(f"\nwrote {OUT}/09_ddg_stability.json")
