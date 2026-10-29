# Track 2 groundwork: BUB1B / BUBR1 mechanism and the repurposing direction problem

Status: data-independent groundwork. Written 2026-09-03, before dataset access.
Everything here is from public sources; nothing uses the gated hackathon data.

Sources searched via PubMed (NCBI), ChEMBL (EMBL-EBI v34), and ClinicalTrials.gov (NIH/NLM).
Open Targets was globally rate-limited at the time of writing and is listed under Pending.

## 0. Target identity

| Field | Value |
|---|---|
| Gene | BUB1B (HGNC:1149) |
| Protein | BUBR1, mitotic checkpoint serine/threonine-protein kinase BUB1 beta |
| UniProt | O60566 |
| ChEMBL target | CHEMBL4295998 (SINGLE PROTEIN, Homo sapiens) |
| Domains | Kinase-like (IPR011009), Mad3/Bub1_I (IPR013212, PF08311) |
| Experimental structures | 2WVI, 3SI5, 4GGD, 5JJA, 5K6S, 5KHU, 5LCW, 5SWF, 6TLJ |
| Predicted structure | AlphaFoldDB O60566 |
| Localisation (GO) | kinetochore, outer kinetochore, mitotic checkpoint complex, APC/C, centrosome |
| Process (GO) | mitotic SAC signalling, metaphase/anaphase transition, centromeric sister chromatid cohesion |

Nine PDB entries exist, but **none resolves the pseudokinase domain (766-1050)**, the two that
map to UniProt 1-1050 model nothing past residue 345 and are 3.8/4.8 A cryo-EM. Structural work
on a variant in that domain rests on the AlphaFold model. See findings/N1002K_STABILITY.md
("Correction to an earlier claim").

## 1. Established mechanism chain

Biallelic BUB1B loss of function
  -> reduced BUBR1 protein
  -> weakened spindle assembly checkpoint
  -> premature chromatid separation (PCS)
  -> mosaic aneuploidy, variegated across chromosomes and tissues
  -> (a) apoptosis of proliferating progenitors -> microcephaly, IUGR, growth failure
     (b) cancer predisposition -> embryonal rhabdomyosarcoma, Wilms tumour, ALL, MDS

Anchors:
- Hanks et al. 2004, Nat Genet, biallelic BUB1B mutations cause MVA. doi:10.1038/ng1449
- Sterling et al. 2023, Front Cell Dev Biol. BubR1-deficient cortical progenitors show
  chromosomal segregation defects, DNA damage, P53 activation and apoptosis; *co-deleting P53
  gives only minimal rescue*, so a powerful P53-independent death pathway also operates.
  doi:10.3389/fcell.2023.1282182
- Pun et al. 2025, J Am Heart Assoc, cardiac BubR1 deletion stalls heart development after
  looping; scRNA-seq shows suppressed cardiogenesis genes and Wnt hyperactivation.
  doi:10.1161/JAHA.124.038286
- Miyamoto & Matsuura 2015, Oncotarget, a ciliopathy component in PCS/MVA (primary cilium,
  mitotic kinesin, PLK1). doi:10.18632/oncotarget.5244
- Schmid et al. 2014, Chromosome Res. BubR1-/- mouse embryos are >50% aneuploid with elevated
  PCS; stochastic massive apoptosis, growth retardation, death E7.5–13.5.
  doi:10.1007/s10577-014-9432-x

Clinical subtypes (García-Castillo et al. 2008, Am J Med Genet A, doi:10.1002/ajmg.a.32315):
patients in whom only one BUB1B mutation was detected were the *most* severely affected
(Dandy-Walker 7/8, cataracts 6/6, Wilms 7/8); biallelic-mutation patients were moderate;
BUB1B-negative patients were mild with no cancer. Read with care ,  "monoallelic" there means
one allele found by exon screening, and the second allele is likely a non-coding variant the
screen missed (see §2).

## 2. Allele architecture is the decisive unknown

A recurring feature of published MVA1 genotypes: the second allele is frequently a
**non-coding splice-disrupting variant** producing reduced but non-zero BUBR1.

- c.2386-11A>G, homozygous intronic, creates a de novo splice site *preferred over* the
  authentic one; adult-onset GI neoplasia phenotype. Rio Frio et al. 2010, NEJM.
  doi:10.1056/NEJMoa1006565
- c.1402-5A>G (intron 10) causing mis-splicing with 4-bp insertion, *in trans* with a paternal
  3.0 kb Alu/SVA transposon insertion in intron 2 found only by RNA-seq + adaptive-sampling
  long-read sequencing. Immunoblot: BUBR1 severely diminished in patient cells.
  Kawakami et al. 2023, J Hum Genet. doi:10.1038/s10038-023-01211-8
- c.498_505delAAACTTTA + c.1288+5G>A; MDS with monosomy 7 at age 9.
  Laberko et al. 2019, Orphanet J Rare Dis. doi:10.1186/s13023-019-1073-x
- c.2363_2364del (p.S788Cfs*29) + a variant in an upstream regulatory element.
  Qu et al. 2025, Zhonghua Yi Xue Yi Chuan Xue Za Zhi. doi:10.3760/cma.j.cn511374-20240716-00393
- c.2679A>T (p.Arg893Ser), 12% mosaic aneuploidy on karyotype.
  Pavone et al. 2022, Neurol Sci. doi:10.1007/s10072-022-06247-w

Consequence for repurposing: **the therapeutic variable is BUBR1 protein abundance**, and in
several published cases at least one allele is leaky. That makes dose restoration a coherent
goal rather than a fantasy. Which strategy applies depends on which two alleles this child
carries:

| Allele architecture | Repurposing goal |
|---|---|
| (i) both null / frameshift | no protein to rescue; target downstream consequences only |
| (ii) leaky splice + null | splice modulation / boost correctly-spliced fraction |
| (iii) missense with residual protein | stabilise the folded variant (chaperone / degradation) |

This is the first thing to determine from the VCF, and it reorganises the whole report.

Note that (ii) is the case where the leaderboard's "compound-het pair" framing matters most,
and note the Kawakami case, where the second allele was invisible to genome sequencing and
required RNA-seq. Worth checking whether this dataset's confirmed answer has the same shape.

## 3. The direction problem: the central scientific point

The literature's ready-made "drugs for aneuploidy" set is aneuploidy-*selective killing*:

- Tang, Williams, Siegel & Amon 2011, Cell. AICAR (energy stress / AMPK), 17-AAG (HSP90,
  protein-folding inhibition) and chloroquine (autophagy inhibition) are lethal to aneuploid
  but not euploid cells. doi:10.1016/j.cell.2011.01.017
  Commentary: Manchado & Malumbres 2011, Cell. doi:10.1016/j.cell.2011.01.037
- Extended: energy-stress AMPK activation sensitises MPS1/TTK inhibition in TNBC.
  Lim et al. 2024, Oncol Rep. doi:10.3892/or.2024.8760
- Caveat on generality: trisomic embryonic stem cells were *not* uniformly sensitive to
  proteotoxic or energy stress, though they were sensitised to genotoxins.
  Zhang et al. 2019, Sci China Life Sci. doi:10.1007/s11427-019-9666-y

These are the first hits for "aneuploidy drug" and they read like a repurposing answer.
They are the wrong direction for constitutional MVA.

In this child, a large fraction of *normal somatic tissue* is aneuploid, and the MVA
phenotype is already driven by excess apoptosis of aneuploid progenitors (Sterling 2023).
An aneuploidy-selective cytotoxic agent given systemically pushes in the same direction as
the disease mechanism.

Two goals must therefore be kept separate, and any candidate must declare which one it serves:

**Goal A, disease modification.** Raise BUBR1 function, reduce ongoing missegregation, or
reduce progenitor apoptosis. Aneuploidy-selective killing is contraindicated here.
Note Sterling 2023: blocking P53 alone does not rescue, so an anti-apoptotic strategy aimed
at the P53 axis has a published reason to fail.

**Goal B, tumour surveillance and treatment.** Aneuploidy-selective agents are on-direction
against an established tumour. But MVA patients tolerate cytotoxic therapy poorly:
- Embryonal RMS in a PCS/MVA infant required *reduced-intensity* chemotherapy.
  Nishitani-Isa et al. 2019, Pediatr Int. doi:10.1111/ped.13849
- HSCT for MVA1-associated MDS: engraftment achieved, then graft rejection at day +47,
  pancytopenia, sepsis, death at day +156. Laberko et al. 2019. doi:10.1186/s13023-019-1073-x

A Goal-B proposal that ignores this hypersensitivity is proposing a dose the patient
population has already demonstrated it cannot take.

## 4. Bounded negative: no therapeutic trial exists for this space

ClinicalTrials.gov query over "mosaic variegated aneuploidy OR premature chromatid separation
OR chromosomal instability OR aneuploidy": **159 studies, zero interventional trials of a
therapy for constitutional aneuploidy or MVA.** The 159 are prenatal aneuploidy screening,
PGT-A / IVF embryo selection, and CIN as a cancer biomarker.

Detection floor: the same query returns single-site feasibility studies with n=30 and
terminated PHASE4 studies with n=34. Any registered interventional MVA therapy trial, however
small, would have been visible. It is absent, not undetected.

Consequences: no comparator arm, no precedent to cite for safety, and no competitor programme
to position against. This is a genuinely empty therapeutic space rather than a crowded one.

## 5. Pending

- **Open Targets**: tractability buckets, known drugs, target-prioritisation factors,
  associated diseases. Blocked on a global rate limit; retry.
- **Boltz structural work**: blocked on knowing the actual variants. Once the alleles are
  known, if either is missense, predict the structural consequence (kinase domain vs TPR vs
  KEN-box / BUB3 interface) and use that to discriminate strategy (iii) from (i).
- **Chemical matter**: no approved drug is known to target BUB1B directly; confirm via Open
  Targets known-drugs before asserting it.

## 6. What changes when the data arrives

One fact reorganises this document: whether the child's two BUB1B alleles leave any residual
BUBR1 protein. Until that is known, Goal A candidate selection cannot be narrowed, and the
structural work has no input.

Nothing above requires the gated data, and nothing above should be presented as a
recommendation for this child's care.
