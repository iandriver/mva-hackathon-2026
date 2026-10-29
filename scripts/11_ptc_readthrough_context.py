#!/usr/bin/env python3
"""Sequence context of the BUB1B p.Leu737Ter PTC, for assessing nonsense readthrough.

Readthrough of a premature termination codon is strongly context-dependent. This extracts
the identity of the stop codon and the flanking positions that Tarraso et al. 2020
(doi:10.1242/dmm.043281) identify as influential (-9, -8, -3, -2, +4, +13, +14; the first
base of the stop is +1), and compares them with the PYGM p.R50X context that study tested.

It also reports NMD competence, because readthrough needs a transcript to act on.
"""
import json, os, sys, urllib.request

OUT = os.environ.get("MVA_RESULTS", os.path.join(os.path.dirname(__file__), "..", "results"))
os.makedirs(OUT, exist_ok=True)
TX, CODON, CDS_POS = "ENST00000287598", 737, 2210
PTC_EXON, TOTAL_EXONS = 17, 23

def cds(tx):
    u = f"https://rest.ensembl.org/sequence/id/{tx}?type=cds;content-type=application/json"
    return json.load(urllib.request.urlopen(u, timeout=90))["seq"]

def main():
    s = cds(TX)
    start = (CODON - 1) * 3
    wt = s[start:start+3]
    off = (CDS_POS - 1) - start
    mut = list(wt); mut[off] = "G"; mut = "".join(mut)
    stops = {"TAA": "ochre", "TAG": "amber", "TGA": "opal"}
    print(f"transcript {TX}   CDS {len(s)} nt, {len(s)//3} codons")
    print(f"codon {CODON}: {wt} -> {mut}   ({stops.get(mut, 'NOT A STOP')})")
    if mut not in stops:
        sys.exit("the substitution does not create a stop codon; check the coordinates")

    ctx = {}
    for p in range(-9, 0):
        ctx[p] = s[start + p]
    for i, b in enumerate(mut, start=1):
        ctx[i] = b
    for p in range(4, 15):
        ctx[p] = s[start + 3 + (p - 4)]

    print(f"\ncontext (first base of the stop = +1):")
    print("  " + "".join(f"{p:>4}" for p in range(-9, 0)) + " |" +
          "".join(f"{p:>4}" for p in range(1, 15)))
    print("  " + "".join(f"{ctx[p]:>4}" for p in range(-9, 0)) + " |" +
          "".join(f"{ctx[p]:>4}" for p in range(1, 15)))

    # PYGM p.R50X, the context Tarraso et al. tested and found NO readthrough in
    pygm = {"stop": "TGA", -9: "G", -1: "G", -3: "C", -2: "C", 4: "not C"}
    print(f"\ncomparison with PYGM p.R50X (Tarraso et al. 2020), where four readthrough")
    print(f"agents produced no detectable protein across three cellular models:")
    rows = [
        ("stop codon", mut, pygm["stop"], "TGA is the most readthrough-permissive stop"),
        ("-9", ctx[-9], pygm[-9], "G reported favourable"),
        ("-3", ctx[-3], pygm[-3], "C reported favourable"),
        ("-2", ctx[-2], pygm[-2], "C reported to COUNTERACT readthrough"),
        ("-1", ctx[-1], pygm[-1], "G reported favourable"),
        ("+4", ctx[4], "not C", "C reported favourable; anything else is not"),
    ]
    print(f"\n  {'position':<12}{'this PTC':<11}{'PYGM R50X':<12}note")
    for pos, ours, theirs, note in rows:
        print(f"  {pos:<12}{ours:<11}{theirs:<12}{note}")

    fav = [r[0] for r in rows if r[0] == "stop codon" and r[1] == "TGA"]
    fav += [r[0] for r in rows if r[0] in ("-9", "-3", "-1") and r[1] == r[2]]
    unfav = ["+4 is not C"] if ctx[4] != "C" else []
    avoided = ["-2 is not C (avoids the counteracting base)"] if ctx[-2] != "C" else []

    codons_to_native_stop = len(s) // 3 - CODON
    print(f"\nNMD competence:")
    print(f"  PTC in exon {PTC_EXON} of {TOTAL_EXONS}; {codons_to_native_stop} codons before the "
          f"native stop.")
    print(f"  A PTC upstream of the last exon-exon junction triggers nonsense-mediated decay,")
    print(f"  so the transcript is degraded and there is less template for readthrough to act")
    print(f"  on. This is the principal obstacle, independent of sequence context.")

    verdict = ("Context is partially favourable (TGA, G at -1, and it avoids the counteracting "
               "C at -2) but not strongly so (+4 is not C, -3 and -9 differ from the "
               "favourable bases). The comparable PYGM TGA context yielded no detectable "
               "readthrough with amlexanox, ataluren, RTC13 or G418. Combined with NMD "
               "competence, readthrough monotherapy is not a well-supported candidate here.")
    print(f"\nVERDICT: {verdict}")
    json.dump({"transcript": TX, "codon": CODON, "wt_codon": wt, "ptc": mut,
               "stop_type": stops[mut], "context": {str(k): v for k, v in ctx.items()},
               "favourable": fav, "unfavourable": unfav, "avoided": avoided,
               "ptc_exon": PTC_EXON, "total_exons": TOTAL_EXONS,
               "codons_to_native_stop": codons_to_native_stop,
               "nmd_competent": True, "verdict": verdict},
              open(os.path.join(OUT, "11_ptc_readthrough_context.json"), "w"), indent=2)
    print(f"\nwrote {OUT}/11_ptc_readthrough_context.json")

if __name__ == "__main__":
    main()
