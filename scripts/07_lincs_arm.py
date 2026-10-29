#!/usr/bin/env python3
"""Reproduces findings/TWO_ARMS.md "Arm A" -- LINCS signature reversal, and why it fails.

Three queries against SigCom-LINCS, in the order that matters:
  1. POSITIVE CONTROL  -- does the BUB1B-KO signature recover BUB1B knockout? (it does, rank 2)
  2. DIRECTION CONTROL -- does BUB1B OVEREXPRESSION reverse it? (it does not; splits both ways)
  3. compound reversers, with the cell-line confound quantified against library background

Control 2 is the one that decides whether the arm means anything. Run it before believing
any candidate list.

NOTE: the enrich endpoint answers a 307 to POST; urllib drops the body on redirect, so
requests must follow it manually.
"""
import collections, gzip, io, json, os, sys, urllib.error, urllib.parse, urllib.request

META = "https://maayanlab.cloud/sigcom-lincs/metadata-api/"
DATA_API = "https://maayanlab.cloud/sigcom-lincs/data-api/api/v1/"
CONSENSUS = ("https://s3.amazonaws.com/lincs-dcic/LINCS-sigs-2021/means/xpr_mean_sigs/"
             "LINCS_L1000_CRISPRKO_{gene}_mean_sig.tsv.gz")
OUT = os.environ.get("MVA_RESULTS", os.path.join(os.path.dirname(__file__), "..", "results"))
os.makedirs(OUT, exist_ok=True)
N_SET = 250

def post(url, payload, timeout=600):
    body = json.dumps(payload).encode()
    req = urllib.request.Request(url, data=body,
                                 headers={"Content-Type": "application/json"})
    try:
        return json.load(urllib.request.urlopen(req, timeout=timeout))
    except urllib.error.HTTPError as e:
        if e.code in (301, 302, 303, 307, 308):
            new = e.headers["Location"]
            req2 = urllib.request.Request(new, data=body,
                                          headers={"Content-Type": "application/json"})
            return json.load(urllib.request.urlopen(req2, timeout=timeout))
        raise

def consensus_signature(gene="BUB1B"):
    raw = urllib.request.urlopen(CONSENSUS.format(gene=gene), timeout=300).read()
    vals = {}
    with gzip.open(io.BytesIO(raw), "rt") as f:
        f.readline()
        for l in f:
            p = l.rstrip("\n").split("\t")
            if len(p) > 1 and p[1]:
                vals[p[0]] = float(p[1])
    return vals

def to_entities(up, down):
    ents = post(META + "entities/find",
                {"filter": {"where": {"meta.symbol": {"inq": up + down}},
                            "fields": ["id", "meta.symbol"], "limit": 2000}})
    U, D = set(up), set(down)
    return {"up_entities": [e["id"] for e in ents if e["meta"]["symbol"] in U],
            "down_entities": [e["id"] for e in ents if e["meta"]["symbol"] in D]}

def resolve(uuids):
    out = {}
    for i in range(0, len(uuids), 100):
        for s in post(META + "signatures/find",
                      {"filter": {"where": {"id": {"inq": uuids[i:i+100]}}, "limit": 200}}):
            out[s["id"]] = s["meta"]
    return out

def query(fe, database, limit):
    res = post(DATA_API + "enrich/ranktwosided", {**fe, "limit": limit, "database": database})
    rev = sorted([x for x in res["results"] if x["type"] == "reversers"], key=lambda x: x["z-sum"])
    mim = sorted([x for x in res["results"] if x["type"] == "mimickers"], key=lambda x: -x["z-sum"])
    return rev, mim

def rank_of(lst, meta, name):
    return [(i + 1, x["z-sum"], meta.get(x["uuid"], {}).get("cell_line"),
             meta.get(x["uuid"], {}).get("pert_time"))
            for i, x in enumerate(lst) if meta.get(x["uuid"], {}).get("pert_name") == name]

def main():
    sig = consensus_signature("BUB1B")
    srt = sorted(sig.items(), key=lambda kv: -kv[1])
    up = [g for g, _ in srt[:N_SET]]; down = [g for g, _ in srt[-N_SET:]]
    print(f"consensus signature: {len(sig):,} genes")
    print(f"  BUB1B in its own DOWN set: {'BUB1B' in down}   <- internal check")
    print(f"  top up genes: {', '.join(up[:12])}")
    fe = to_entities(up, down)
    print(f"  mapped {len(fe['up_entities'])} up / {len(fe['down_entities'])} down\n")

    print("=" * 72)
    print("1. POSITIVE CONTROL -- CRISPR library (l1000_xpr)")
    rev, mim = query(fe, "l1000_xpr", 200)
    m = resolve([x["uuid"] for x in mim])
    hits = rank_of(mim, m, "BUB1B")
    for r, z, cl, t in hits:
        print(f"   BUB1B KO recovered as mimicker: rank {r}, z-sum {z:.2f}, {cl}")
    print(f"   PASS: {bool(hits)}\n")

    print("=" * 72)
    print("2. DIRECTION CONTROL -- overexpression library (l1000_oe)")
    print("   Premise under test: BUB1B OE should REVERSE BUB1B KO.")
    rev_oe, mim_oe = query(fe, "l1000_oe", 2000)
    mr = resolve([x["uuid"] for x in rev_oe]); mm = resolve([x["uuid"] for x in mim_oe])
    as_rev = rank_of(rev_oe, mr, "BUB1B"); as_mim = rank_of(mim_oe, mm, "BUB1B")
    for r, z, cl, t in as_rev: print(f"   reverser  rank {r:>5}  z {z:+.2f}  {cl} {t}")
    for r, z, cl, t in as_mim: print(f"   MIMICKER  rank {r:>5}  z {z:+.2f}  {cl} {t}")
    consistent = len(as_rev) / max(1, len(as_rev) + len(as_mim))
    print(f"\n   {len(as_rev)} reverse, {len(as_mim)} mimic -> {consistent*100:.0f}% consistent")
    print(f"   {'PASS' if consistent >= 0.75 else 'FAIL'}: loss and gain are "
          f"{'mirror images' if consistent >= 0.75 else 'NOT mirror images'}; "
          f"{'reversal is meaningful' if consistent >= 0.75 else 'reversing the loss signature has no validated meaning'}\n")

    print("=" * 72)
    print("3. compound reversers (l1000_cp) -- with the confound quantified")
    rev_cp, _ = query(fe, "l1000_cp", 200)
    mc = resolve([x["uuid"] for x in rev_cp])
    rows = [((mc.get(x["uuid"], {}).get("pert_name") or "?"), x["z-sum"],
             mc.get(x["uuid"], {}).get("cell_line", "")) for x in rev_cp]
    for i, (n, z, cl) in enumerate(rows[:15], 1):
        print(f"   {i:>2}. {n:<26} z {z:7.2f}  {cl}")
    counts = collections.Counter(cl for _, _, cl in rows)
    top_cl, top_n = counts.most_common(1)[0]

    def count(where):
        u = META + "signatures/count?where=" + urllib.parse.quote(json.dumps(where))
        return json.load(urllib.request.urlopen(u, timeout=120))["count"]
    bg = count({"meta.cell_line": top_cl, "meta.pert_type": "Chemical"})
    allc = count({"meta.pert_type": "Chemical"})
    print(f"\n   top cell line {top_cl}: {top_n}/{len(rows)} = {top_n/len(rows)*100:.1f}% of hits")
    print(f"   background rate in library: {bg:,}/{allc:,} = {bg/allc*100:.2f}%")
    print(f"   ENRICHMENT: {(top_n/len(rows))/(bg/allc):.1f}x  <- context confound")

    json.dump({"positive_control_rank": hits, "direction_control":
               {"as_reverser": as_rev, "as_mimicker": as_mim, "consistency": consistent},
               "top_reversers": rows[:25],
               "confound": {"cell_line": top_cl, "in_hits": top_n,
                            "background_frac": bg/allc,
                            "enrichment": (top_n/len(rows))/(bg/allc)}},
              open(os.path.join(OUT, "07_lincs_arm.json"), "w"), indent=2)
    print(f"\nwrote {OUT}/07_lincs_arm.json", file=sys.stderr)

if __name__ == "__main__":
    main()
