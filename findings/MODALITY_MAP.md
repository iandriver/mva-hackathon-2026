# Which modalities are actually ruled out, and why

Run 2026-09-07. Earlier documents said "BUBR1 is not directly druggable". That claim was
scoped to small molecules and stated more broadly than the evidence supported. This is the
modality-by-modality version.

Three distinct kinds of exclusion, kept separate because they have different remedies:

- **Biology** — the modality cannot work on this target, whatever the delivery or approval status.
- **Delivery** — the mechanism is sound but nothing reaches the affected tissues.
- **Brief** — scientifically live, but not an *existing approved medication*, which is what
  Track 2 asks for.

## Ruled out on biology

| Modality | Why |
|---|---|
| **Small molecule, direct** | No binding site. The pseudokinase ATP pocket is occluded at ~63 A^3 against 300-600 typical for a druggable kinase site, and residue 1002 sits 33 A away (TWO_ARMS.md). Independently corroborated by the absence of any chaperone precedent for a cytosolic protein lacking a pocket (CHAPERONE_CANDIDATES.md). |
| **Antibodies and biologics** | BUBR1 is cytosolic and nuclear. Conventional antibodies do not cross the plasma membrane, so the target is unreachable. Intrabodies exist as research reagents; none is approved, and none has a systemic delivery route. |
| **siRNA, shRNA, PROTACs, molecular glue degraders** | All reduce protein. This child needs *more* BUBR1 function, not less. Wrong direction, and the direction gate that governs everything else here applies to modality choice too. |
| **Base editing** | Both pathogenic alleles are T>G, so correcting either requires a **G>T transversion** (C>A on the antisense strand). Cytosine and adenine base editors perform transitions only. Neither allele is base-editable. Computed in `scripts/13_editability.py`. |
| ~~**AAV gene addition**~~ | **Moved. See the correction below: this argument was against one vector, not against gene delivery.** |

## Ruled out on delivery

| Modality | Why |
|---|---|
| **Antisense oligonucleotides** | The only sound sub-modes here are NMD evasion on the PTC allele and translational upregulation of the hypomorph; knockdown is the wrong direction and splice modulation has nothing to correct, since both variants are exonic rather than splice-disrupting. ASO biodistribution is essentially liver, kidney and CNS-by-intrathecal. MVA affects growth, brain, kidney, muscle and every proliferating compartment at once. No route covers that tissue set. |
| **Prime editing, nuclease plus repair template** | The only editing chemistry that *can* make a G>T transversion, so not excluded on biology. Excluded on reaching the affected tissues, and on the brief below. |
| **mRNA therapy** | Transient expression requiring indefinite repeat dosing to every dividing tissue. No approved product does chronic intracellular protein replacement. |
| **Peptides and PPI stabilisers** | The functional output to restore is the BubR1-KARD-PP2A-B56 interaction (Gama Braga et al. 2020, doi:10.1016/j.celrep.2020.108397). Stabilising a protein-protein interaction is substantially harder than blocking one, and peptides do not cross the plasma membrane. |
| **Ex vivo lentiviral HSC gene therapy** | Addresses the haematopoietic compartment only, leaving growth, brain and solid-tumour risk untouched. HSCT in MVA1 has been tried: engraftment, then graft rejection at day +47 and death at day +156 (Laberko et al. 2019, doi:10.1186/s13023-019-1073-x). |

## Live on biology, excluded by the brief

| Modality | Status |
|---|---|
| **Suppressor tRNA therapy** | Engineered tRNAs that decode the PTC. Mechanistically the cleanest answer to a TGA nonsense allele and does not depend on the poor sequence context that weakens small-molecule readthrough (PTC_READTHROUGH.md). Entirely preclinical; delivery unsolved. |
| **Prime editing** | See above. Real, not approved. |
| **CRISPRa** | Transcriptional upregulation would raise the hypomorph, but also the PTC allele whose message is degraded anyway. No approved product; delivery unsolved. |

## What this leaves

Track 2 asks for **existing, approved medications**. Every modality that survives the biology
and delivery gates fails the brief, and every modality that satisfies the brief is a small
molecule, which fails on the absent binding site.

**The answer, if there is one, must be a small molecule acting downstream of BUBR1 rather than
on it.** That is the same conclusion the chaperone analysis reached, arrived at independently
and with the alternatives now closed by name rather than by omission.

The downstream space that has been searched is the type-I interferon module, which passed
specificity and drug-direction controls and failed the safety-direction gate in a child with
an existing cancer predisposition (INTERFERON_HYPOTHESIS.md). Aneuploidy-tolerance biology and
cancer surveillance remain unsearched.

## Correction: gene delivery, revisited

The AAV argument above was correct and too narrow. Episomal persistence is a property of
standard AAV, not of gene delivery, and I wrote it as though it closed the modality.

**Persistence is solvable.** Lentiviral vectors integrate. Transposon systems (Sleeping
Beauty, piggyBac) integrate. Genome editing is permanent by construction. Any of these
survives cell division, so "the vector is diluted through mitosis" argues against AAV
specifically and not against correcting a gene needed in dividing cells.

**LNP is real and is not the answer here, for two separate reasons.** Conventional LNPs
"primarily accumulate in the liver following systemic administration, limiting their success
to liver-centric conditions" (Geisler et al. 2024, doi:10.1016/j.jconrel.2024.05.036), driven
by ApoE-mediated tropism (Isaac et al. 2025, doi:10.1002/smll.202503034). Retargeting works:
antibody-conjugated LNPs roughly double placental delivery, and corticosteroid-derived LNPs
reach 97-99% spleen selectivity. But these redirect to *one* organ at a time, in mice, and
preclinically. Separately, the usual LNP cargo is mRNA, which is transient, so it does not
address a chronic requirement in dividing cells whatever the tropism.

**So the binding constraint is not persistence. It is tissue breadth, and timing.**

*Breadth.* Aneuploidy is generated at every mitosis, so correction is needed in neural
progenitors, growth plate, kidney, muscle and the haematopoietic compartment at once. No
platform, viral or lipid, reaches that set. Single-organ redirection does not compose into
whole-body coverage.

*Timing.* Growth restriction in MVA begins before birth, and this child's presentation
includes it, so a substantial part of the developmental phenotype was established in utero.
Even flawless delivery after diagnosis addresses ongoing aneuploidy generation and future
cancer risk; it does not reverse growth restriction already sustained. Any gene-therapy
proposal should say which of those two it is aiming at.

### A correction I owe to ex vivo HSC gene therapy

The delivery table dismisses ex vivo lentiviral HSC gene therapy partly by citing the MVA1
transplant that ended in graft rejection and death (Laberko et al. 2019). That reasoning is
wrong. **That failure was allogeneic graft rejection from a matched unrelated donor.
Autologous, gene-corrected cells carry no allogeneic mismatch and would not face it.**

Autologous HSC correction remains narrow: it addresses the MDS and leukaemia risk, which is
real in MVA1 and was the indication in that case, and leaves growth, brain and solid-tumour
risk untouched. It is also not an approved medicine for this indication. But it is a more
defensible option than the earlier table implied, and it was dismissed with the wrong
argument.

## Correction to the earlier claim

CHAPERONE_CANDIDATES.md ends "BUBR1 is not directly druggable, established two independent
ways". That is correct **for small molecules** and was written without the modality scope
stated. The accurate claim is: BUBR1 has no small-molecule binding site, is unreachable by
antibodies, is not base-editable, and is a poor fit for AAV. Modalities that could work in
principle are excluded by delivery or by not existing as medicines yet.
