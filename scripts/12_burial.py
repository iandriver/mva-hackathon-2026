#!/usr/bin/env python3
"""Is BUBR1 residue 1002 buried or surface-exposed?

This gates the chemical-chaperone hypothesis. A destabilising substitution at a BURIED
position disrupts the hydrophobic core, and the misfolded species is what quality control
degrades -- the situation chaperones address. A substitution at an EXPOSED position is
usually tolerated, the protein largely keeps its fold, and there is less for a chaperone
to rescue even when a ddG predictor calls it destabilising.

Relative solvent accessibility via Shrake-Rupley on the AlphaFold kinase domain, with
Tien et al. 2013 theoretical maximum ASA values. Convention: RSA < 0.20 buried,
0.20-0.50 intermediate, > 0.50 exposed.
"""
import json, os, sys
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

REPO = os.path.join(os.path.dirname(__file__), "..")
OUT = os.environ.get("MVA_RESULTS", os.path.join(REPO, "results"))
PDB = os.path.join(REPO, "structures", "BUBR1_AF_kinase_766_1050.pdb")
TARGET = 1002
# Tien et al. 2013, theoretical maxima (doi:10.1371/journal.pone.0080635)
MAXASA = {"ALA":129,"ARG":274,"ASN":195,"ASP":193,"CYS":167,"GLN":225,"GLU":223,"GLY":104,
          "HIS":224,"ILE":197,"LEU":201,"LYS":236,"MET":224,"PHE":240,"PRO":159,"SER":155,
          "THR":172,"TRP":285,"TYR":263,"VAL":174}

def main():
    if not os.path.exists(PDB):
        sys.exit(f"MISSING: {PDB}\n  run scripts/09_ddg_stability.sh first")
    st = PDBParser(QUIET=True).get_structure("m", PDB)
    ShrakeRupley().compute(st[0], level="R")
    rsa = {}
    for res in st[0].get_residues():
        n = res.get_resname()
        if n in MAXASA:
            rsa[res.id[1]] = res.sasa / MAXASA[n]
    vals = sorted(rsa.values())
    t = rsa.get(TARGET)
    if t is None: sys.exit(f"residue {TARGET} not in the model")

    def band(v):
        return "buried" if v < 0.20 else ("intermediate" if v <= 0.50 else "exposed")

    import bisect
    pct = bisect.bisect_left(vals, t) / len(vals) * 100
    print(f"kinase domain 766-1050: {len(rsa)} residues scored")
    print(f"\nresidue {TARGET} (Asn):")
    print(f"  SASA          {st[0]['A'][TARGET].sasa:.1f} A^2")
    print(f"  RSA           {t:.3f}")
    print(f"  classified    {band(t).upper()}")
    print(f"  percentile    {pct:.1f}% of domain residues are less exposed")

    buried = sum(1 for v in vals if v < 0.20)
    print(f"\n  domain composition: {buried} buried, "
          f"{sum(1 for v in vals if 0.20 <= v <= 0.50)} intermediate, "
          f"{sum(1 for v in vals if v > 0.50)} exposed")

    # neighbours within 8 A of the CB (or CA for Gly) as a packing cross-check
    import math
    tgt = st[0]["A"][TARGET]
    ref = tgt["CB"] if "CB" in tgt else tgt["CA"]
    nb = 0
    for res in st[0].get_residues():
        if res.id[1] == TARGET: continue
        a = res["CB"] if "CB" in res else (res["CA"] if "CA" in res else None)
        if a is not None and (a - ref) < 8.0: nb += 1
    print(f"  packing: {nb} residues with CB within 8 A")

    if t < 0.20:
        verdict = ("BURIED. A charge-introducing substitution here disrupts the core, the "
                   "misfolded species is a degradation substrate, and a chaperone hypothesis "
                   "has a target.")
    elif t <= 0.50:
        verdict = ("INTERMEDIATE. Partially packed. Destabilisation is plausible but the "
                   "protein is unlikely to be grossly misfolded; chaperone rescue is a weaker "
                   "proposition than for a buried substitution.")
    else:
        verdict = ("EXPOSED. Asn to Lys at a solvent-exposed position is usually well "
                   "tolerated. The ddG prediction is likely dominated by conservation rather "
                   "than packing, and the chemical-chaperone hypothesis loses its target: "
                   "there is little misfolding to rescue.")
    print(f"\nVERDICT: {verdict}")
    json.dump({"residue": TARGET, "sasa": st[0]["A"][TARGET].sasa, "rsa": t,
               "class": band(t), "percentile_less_exposed": pct,
               "neighbours_8A": nb, "n_residues": len(rsa), "verdict": verdict},
              open(os.path.join(OUT, "12_burial.json"), "w"), indent=2)
    print(f"wrote {OUT}/12_burial.json")

if __name__ == "__main__":
    main()
