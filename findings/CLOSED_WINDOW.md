# Aneuploidy tolerance and cancer surveillance: the window is closed in both directions

Run 2026-09-07. This completes the downstream search and explains the shape of every
negative that preceded it.

## Aneuploidy tolerance

Tolerance is the cellular property that lets a cell survive and propagate an aneuploid genome.
The mechanisms are known and specific:

- **BCL9L loss** reduces basal caspase-2, preventing cleavage of MDM2 and BID, and permits
  aneuploidy tolerance in both TP53 wild-type and mutant cells (Lopez-Garcia et al. 2017,
  Cancer Cell, doi:10.1016/j.ccell.2016.11.001).
- **Cyclin D1 overexpression** overrides an intact p53/p21 checkpoint to permit tolerance of
  genome doubling (Crockford et al. 2017, Ann Oncol, doi:10.1093/annonc/mdw612).

Both papers describe tolerance as a **cancer** mechanism. Lopez-Garcia frame it as what
"permits the propagation of aneuploid genomes", driving "cancer evolution, intratumor
heterogeneity, and drug resistance", and propose exploiting the BCL9L/caspase-2/BID axis to
*limit* cancer diversity.

That framing is the problem. In this child, tolerance is simultaneously:

| | |
|---|---|
| what his **somatic tissue needs more of** | aneuploid progenitors are dying. That apoptosis is what produces the growth restriction, and in published models the microcephaly (Sterling et al. 2023, doi:10.3389/fcell.2023.1282182) |
| what his **premalignant cells need less of** | tolerance is precisely what allows an aneuploid cell to persist and become a rhabdomyosarcoma or a Wilms tumour |

**Raising tolerance treats the developmental disease and accelerates the cancer risk.
Lowering it does the reverse.** These are the same cells and the same mechanism, so the two
effects cannot be separated by tissue, by timing, or by dose.

The therapeutic window on aneuploidy tolerance in MVA is closed in both directions, and it is
closed for a mechanistic reason rather than for want of a drug.

This is the general form of the direction problem first noted in MECHANISM.md §3. There the
observation was that the aneuploidy-selective compounds (AICAR, 17-AAG, chloroquine) point the
wrong way for constitutional MVA. The stronger statement is that **no setting of the tolerance
dial is safe here**, so the whole axis is unavailable rather than merely mis-signed.

## Cancer surveillance and chemoprevention

The precedent worth testing is aspirin in Lynch syndrome, the clearest case of an approved
drug reducing cancer incidence in a germline predisposition syndrome. It does not transfer:

- **Different mechanism.** Lynch is mismatch-repair deficiency; the proposed aspirin effect
  runs through COX/inflammation and immune surveillance of a hypermutated, neoantigen-rich
  tumour. MVA is chromosomal instability, which is a different mutational process entirely.
- **Different tumour spectrum and age.** Lynch is adult-onset carcinoma. MVA is childhood
  embryonal tumours, mostly in the first years of life.
- **Contraindicated in the patient.** Aspirin in young children carries Reye syndrome risk.
- **The evidence is not definitive even in Lynch.** "No definitive data show efficacy of
  chemopreventive drugs, although aspirin is a promising drug" (Guillen-Ponce et al. 2012,
  doi:10.1586/era.12.114; see also Ricciardiello et al. 2016, doi:10.1038/nrgastro.2016.56).

What exists for childhood cancer predisposition is **surveillance imaging**, not medication.
That is real clinical value and it is outside a drug-repurposing brief.

## Why every arm failed, in one picture

The three places an intervention could act, and what closes each:

| Level | Intervention | Status |
|---|---|---|
| **Upstream** | reduce the missegregation rate | Would help the developmental phenotype and the cancer risk **at once**, because it prevents new aneuploidy rather than adjudicating between existing aneuploid cells. Requires restoring BUBR1 function. No small-molecule site, and every other modality fails on biology, delivery, or the approved-medicine brief (MODALITY_MAP.md). |
| **Downstream, tolerance** | raise or lower it | Closed in both directions, above. |
| **Downstream, consequence** | the type-I interferon module | Passed specificity and drug-direction controls, failed the safety-direction gate in a cancer-predisposed child (INTERFERON_HYPOTHESIS.md). |

The single place where the interests of the somatic tissue and the tumour risk **align** is
upstream, at the missegregation rate. That is exactly the level at which the target is
undruggable. Every level where a drug could act is a level where the two interests oppose.

That is the structural reason MVA has no treatment beyond symptom management, and it is a
sharper statement than "no candidate was found".

## Scope

I did not find this argument stated in the MVA or aneuploidy literature, but the searches were
targeted rather than exhaustive and absence from them is not proof. The components are all
published; the combination is what is being asserted here.
