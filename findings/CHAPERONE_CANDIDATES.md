# Chemical and pharmacological chaperones for p.Asn1002Lys: refused

Run 2026-09-07. Reproduce the structural gate with `scripts/12_burial.py`.

Following PTC_READTHROUGH.md, the hypomorphic allele is the better target: it makes
full-length protein destabilised by about 1.5 kcal/mol. Chaperone approaches are the natural
candidate class, and several agents are approved.

## The gate the hypothesis passes

A ddG prediction is not evidence that a protein is degraded. A destabilising substitution at a
**buried** position disrupts the core and creates a quality-control substrate; the same
substitution at an **exposed** position is usually tolerated and leaves little to rescue.

Shrake-Rupley on the AlphaFold kinase domain, normalised by Tien et al. 2013 maximum ASA
(doi:10.1371/journal.pone.0080635):

| | |
|---|---|
| residue 1002 SASA | 31.5 A^2 |
| **RSA** | **0.162** |
| classification | **buried** (< 0.20) |
| neighbours with CB within 8 A | 8 |
| domain composition | 143 buried, 97 intermediate, 45 exposed |

Buried, though only just, and moderately rather than densely packed. That is consistent with
the modest ddG: a deeply buried substitution introducing charge would score far above 1.5
kcal/mol. The picture is a hypomorph, not a dead protein, and it matches a milder clinical
presentation.

So the hypothesis has a target. It fails on the next gate.

## Why it fails

Every pharmacological chaperone that works in the literature shares two properties, and BUBR1
has neither.

**They act on secretory-pathway proteins.** CFTR correctors relieve an ER folding defect
(Baatallah et al. 2021, doi:10.1007/s00018-021-03994-5). NB-DGJ rescues lysosomal
beta-galactosidase trafficking in GM1-gangliosidosis (Mohamed et al. 2020,
doi:10.1007/s00439-020-02153-3). MC4R chaperones rescue a receptor retained intracellularly
(Huang et al. 2017, doi:10.1016/j.bbadis.2017.03.001). BUBR1 is cytosolic and nuclear. There is
no ER quality-control bottleneck to relieve, which is the mechanism these exploit.

**They bind a specific site on the target.** VX-809 and VX-445 occupy defined pockets on CFTR
MSD1 and NBD1; NB-DGJ is an active-site-directed iminosugar. **BUBR1 offers no such site.**
Our own pocket analysis (TWO_ARMS.md, Arm B) found the pseudokinase ATP site occluded at
roughly 63 A^3 against the 300-600 A^3 typical of a druggable kinase pocket, with residue 1002
a further 33 A away. The geometry that closed the docking arm closes this one: a
pharmacological chaperone needs somewhere to bind.

Generic chemical chaperones avoid the binding-site problem by acting on solvent properties
rather than the target, and pay for it. Huang et al. describe their therapeutic application as
"very limited due to their non-specific mechanism of action and, for chemical chaperone, high
dosage needed to be effective". 4-phenylbutyrate, TMAO and TUDCA fall here.

## The one branch that survives the pocket argument, and where it stops

Network-level proteostasis modulators, such as HSP co-inducers (arimoclomol), do not bind the
target at all; they amplify the endogenous heat shock response. The pocket argument does not
touch them.

They require BUBR1 N1002K to be a chaperone client. A targeted search returned no report of
BUBR1 as an HSP70 or HSP90 client. Without that, the approach has no mechanistic connection to
this variant, and proposing it would be reasoning from drug class rather than from mechanism,
which is what the direction gate exists to prevent.

## Consequence

**BUBR1 has no small-molecule binding site, established two independent ways** — an occluded
pocket found geometrically, and the absence of any working chaperone precedent for a cytosolic
protein without a binding site. That claim is scoped to small molecules; MODALITY_MAP.md works
through antibodies, oligonucleotides, editing, gene addition and peptides separately, and
reaches the same place by a longer route. Any Track 2 candidate must act downstream of the
protein rather than on it.

The downstream consequence we quantified is the type-I interferon module
(INTERFERON_HYPOTHESIS.md), which passed specificity and drug-direction controls and failed on
safety direction in a child with an existing cancer predisposition.

Five candidate classes have now been assessed and refused: aneuploidy-selective cytotoxics,
LINCS signature reversal, ATP-pocket docking, nonsense readthrough, and chaperones.
