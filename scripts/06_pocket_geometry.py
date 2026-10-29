#!/usr/bin/env python3
"""Reproduces findings/TWO_ARMS.md "Arm B" -- the check run BEFORE docking anything.

Asks two questions of the BUBR1 pseudokinase ATP site on the AlphaFold model:
is there a pocket, and is the variant near it. Both answers were no, which is why
no ligands were docked.
"""
import json, math, os, sys, urllib.request

OUT = os.environ.get("MVA_RESULTS", os.path.join(os.path.dirname(__file__), "..", "results"))
STRUCT = os.path.join(os.path.dirname(__file__), "..", "structures")
os.makedirs(OUT, exist_ok=True); os.makedirs(STRUCT, exist_ok=True)

UNIPROT = "O60566"
ATP_SITE = list(range(772, 781)) + [795]   # UniProt binding-site annotations
VARIANT = 1002
PROBE = 3.3        # a water-sized probe; points further than this from any atom are open
RADIUS = 8.0
GRID = 1.0

def alphafold_model():
    path = os.path.join(STRUCT, "BUBR1_AF_full.pdb")
    if not os.path.exists(path):
        meta = json.load(urllib.request.urlopen(
            f"https://alphafold.ebi.ac.uk/api/prediction/{UNIPROT}", timeout=90))[0]
        urllib.request.urlretrieve(meta["pdbUrl"], path)
    return path

def atoms(path, lo=None, hi=None):
    out = []
    for l in open(path):
        if l.startswith("ATOM"):
            r = int(l[22:26])
            if lo and not (lo <= r <= hi): continue
            out.append((r, l[12:16].strip(), float(l[30:38]), float(l[38:46]), float(l[46:54]),
                        float(l[60:66])))
    return out

def main():
    path = alphafold_model()
    dom = atoms(path, 766, 1050)
    if not dom: sys.exit("no kinase-domain atoms parsed")

    ca = {r: b for r, n, *_ , b in [(a[0], a[1], a[2], a[3], a[4], a[5]) for a in dom] if n == "CA"}
    print(f"pLDDT at {VARIANT}: {ca.get(VARIANT)}")
    print(f"pLDDT mean, kinase domain 766-1050: "
          f"{sum(ca.values())/len(ca):.1f}  (n={len(ca)})")

    site = [a for a in dom if a[0] in ATP_SITE]
    cx = sum(a[2] for a in site)/len(site)
    cy = sum(a[3] for a in site)/len(site)
    cz = sum(a[4] for a in site)/len(site)
    print(f"\nATP-site centroid: ({cx:.2f}, {cy:.2f}, {cz:.2f})  from residues {sorted({a[0] for a in site})}")

    free = tot = 0
    steps = [i*GRID for i in range(-int(RADIUS/GRID), int(RADIUS/GRID)+1)]
    for dx in steps:
        for dy in steps:
            for dz in steps:
                if math.dist((dx, dy, dz), (0, 0, 0)) > RADIUS: continue
                tot += 1
                p = (cx+dx, cy+dy, cz+dz)
                if min(math.dist(p, (a[2], a[3], a[4])) for a in dom) > PROBE:
                    free += 1
    print(f"\ncavity probe ({RADIUS:.0f} A sphere, {GRID:.0f} A grid, {PROBE} A probe):")
    print(f"  grid points sampled     : {tot}")
    print(f"  probe-accessible        : {free}  ({free/tot*100:.1f}%)")
    print(f"  approx accessible volume: {free} A^3")
    print(f"  (a druggable kinase ATP site is typically 300-600 A^3)")

    v = [a for a in dom if a[0] == VARIANT]
    vx = sum(a[2] for a in v)/len(v); vy = sum(a[3] for a in v)/len(v); vz = sum(a[4] for a in v)/len(v)
    dist = math.dist((cx, cy, cz), (vx, vy, vz))
    print(f"\ndistance ATP-site centroid -> residue {VARIANT}: {dist:.1f} A")

    verdict = ("NOT a well-posed docking target: pocket is occluded and the variant is remote"
               if free < 150 and dist > 15 else "pocket and variant warrant docking")
    print(f"\nVERDICT: {verdict}")
    json.dump({"plddt_variant": ca.get(VARIANT), "accessible_volume_A3": free,
               "grid_points": tot, "atp_centroid": [cx, cy, cz],
               "dist_atp_to_variant_A": dist, "verdict": verdict},
              open(os.path.join(OUT, "06_pocket_geometry.json"), "w"), indent=2)

if __name__ == "__main__":
    main()
