#!/usr/bin/env python3
"""Which genome-editing modalities could correct this child's two BUB1B alleles?

Base editors perform TRANSITIONS only:
  cytosine base editors (CBE)  C -> T   (and G -> A on the opposite strand)
  adenine base editors  (ABE)  A -> G   (and T -> C on the opposite strand)
A correction requiring a TRANSVERSION (purine <-> pyrimidine) is outside base-editor
chemistry and needs prime editing or nuclease-plus-template repair.

This checks each variant on both strands.
"""
import json, os

OUT = os.environ.get("MVA_RESULTS", os.path.join(os.path.dirname(__file__), "..", "results"))
os.makedirs(OUT, exist_ok=True)
COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}
PURINES, PYRIMIDINES = {"A", "G"}, {"C", "T"}

VARIANTS = [
    {"name": "p.Leu737Ter", "hgvs": "c.2210T>G", "ref": "T", "alt": "G",
     "note": "nonsense, null allele"},
    {"name": "p.Asn1002Lys", "hgvs": "c.3006T>G", "ref": "T", "alt": "G",
     "note": "missense, hypomorphic allele"},
]

def kind(frm, to):
    both_pur = frm in PURINES and to in PURINES
    both_pyr = frm in PYRIMIDINES and to in PYRIMIDINES
    return "transition" if (both_pur or both_pyr) else "transversion"

def editable(frm, to):
    """Can a base editor make this change, on either strand?"""
    hits = []
    if (frm, to) == ("C", "T"): hits.append("CBE, sense strand")
    if (frm, to) == ("A", "G"): hits.append("ABE, sense strand")
    f, t = COMP[frm], COMP[to]
    if (f, t) == ("C", "T"): hits.append("CBE, antisense strand")
    if (f, t) == ("A", "G"): hits.append("ABE, antisense strand")
    return hits

rows = []
print("Correcting each variant means reverting ALT back to REF.\n")
for v in VARIANTS:
    frm, to = v["alt"], v["ref"]              # correction direction
    k = kind(frm, to)
    hits = editable(frm, to)
    rows.append({**v, "correction": f"{frm}>{to}", "class": k,
                 "base_editable": bool(hits), "editors": hits})
    print(f"{v['name']}  ({v['hgvs']}, {v['note']})")
    print(f"  correction needed : {frm} -> {to}  on the sense strand")
    print(f"                      {COMP[frm]} -> {COMP[to]}  on the antisense strand")
    print(f"  class             : {k.upper()}")
    print(f"  base-editable     : {'YES via ' + ', '.join(hits) if hits else 'NO'}")
    print()

n_edit = sum(1 for r in rows if r["base_editable"])
print("=" * 70)
if n_edit == 0:
    print("NEITHER allele is base-editable. Both pathogenic changes are T>G, so correcting")
    print("either requires a G>T transversion (C>A on the antisense strand). Base editors do")
    print("transitions only. Correction would require prime editing or nuclease plus a repair")
    print("template, neither of which is an approved medicine.")
else:
    print(f"{n_edit} of {len(rows)} alleles fall within base-editor chemistry.")
print("=" * 70)
json.dump({"variants": rows, "n_base_editable": n_edit,
           "conclusion": ("neither allele is base-editable; both corrections are "
                          "transversions" if n_edit == 0 else "partially editable")},
          open(os.path.join(OUT, "13_editability.json"), "w"), indent=2)
print(f"wrote {OUT}/13_editability.json")
