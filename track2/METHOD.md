# Track 2 method: BUB1B signature reversal, with an internal direction control

Status: design fixed, data-independent. The LINCS content it depends on has been verified to
exist (§1). Written 2026-09-03.

## Why this design

The Track 2 report needs a mechanism-to-drug link that is reproducible and does not depend on
the child's transcriptome, the hackathon releases FASTQ + VCF, not RNA. Perturbation
signature reversal supplies that link from public data: the disease direction is defined by
BUB1B loss-of-function signatures, and candidate drugs are those whose own signature runs
against it.

The design's distinguishing feature is that the direction is controlled internally rather than
assumed. See §3.

## 1. The resource exists and covers the relevant contexts

SigCom-LINCS (`https://maayanlab.cloud/sigcom-lincs/metadata-api/`, public, no key) holds
**139 BUB1B perturbation signatures**:

| Perturbation | n |
|---|---|
| shRNA knockdown | 89 |
| CRISPR knockout | 32 |
| **Overexpression** | **18** |

Including a consensus signature `LINCS_L1000_CRISPRKO_BUB1B_mean_sig`, and coverage across
27 contexts. Most are cancer lines (VCAP 27, PC3 19, A375 13, MCF7 11, A549 10, HT29 9 …),
but two non-cancer contexts are present:

- **NPC, neural progenitor cells, 3 signatures.** This is the cell type in which the MVA
  microcephaly phenotype originates: BubR1-deficient cortical progenitors undergo massive
  apoptosis (Sterling et al. 2023, doi:10.3389/fcell.2023.1282182).
- ASC, adipose stromal cells, 3 signatures.

## 2. Pipeline

1. **Define the loss direction.** Consensus CRISPR-KO signature plus the shRNA set; carry the
   NPC subset separately as the disease-relevant context rather than pooling it away.
2. **Define the restoration direction.** The 18 overexpression signatures are what *restored*
   BUBR1 looks like transcriptionally.
3. **Query for reversers** of the loss signature (negative z-sum in SigCom-LINCS) across the
   compound libraries.
4. **Filter to approved drugs** via ChEMBL (`max_phase=4`), carrying black-box-warning and
   withdrawal flags forward, this is a paediatric patient.
5. **Structural arm**, only if the child carries a missense allele (see MECHANISM.md §2).
6. **Mechanism plausibility** and a literature dossier per candidate gene/drug entering the
   record. "No evidence found" is a valid dossier answer; blank is not.

## 3. The internal direction control

Two queries that must agree:

- drugs that **reverse** the BUB1B knockdown/knockout signature
- drugs that **mimic** the BUB1B overexpression signature

These are independent queries against different signature sets answering the same question.
Candidates appearing in both are supported; candidates appearing in only one are a defect or a
context effect, and get reported as such rather than dropped silently.

Most submissions will run one direction. Running both is what makes the direction claim
checkable instead of asserted, and given that the naive reading of the aneuploidy literature
points the *wrong way* (MECHANISM.md §3), a checkable direction is the point.

## 4. The failure mode to screen for

Reversing a BUB1B-knockdown signature will surface antimitotics and general cytotoxics, because
knocking down a mitotic checkpoint kinase perturbs proliferation genes and anything that
arrests the cell cycle anticorrelates with that. Those are not Goal A candidates, they are the
direction problem wearing a signature. Screen by asking whether the candidate's proposed
mechanism touches BUBR1 abundance, kinetochore function, or progenitor survival, rather than
whether its correlation coefficient is large.

## 5. Stated limits

- L1000 measures 978 landmark genes with the rest inferred; this is not a full transcriptome.
- CRISPR knockout is complete loss. The child is very likely partial loss (MECHANISM.md §2).
  The signature is therefore a dose mismatch, in the direction of overstating the effect.
- NPC n=3 is thin. Any NPC-specific claim needs its own noise estimate before it is stated.
- Signature reversal has a weak record of predicting clinical efficacy. The output is a ranked
  hypothesis list. The hackathon asks for exactly that ,  "hypotheses for follow-up, not
  evidence that a medicine works", and the report should not exceed it.

## 6. Tooling

**Ma'ayan Lab tools skills**, cloned to `tools/maayan-lab-tools-skills/`, 14 skills in the
Agent Skills format (one `SKILL.md` each), enabled by symlinking into `~/.claude/skills/` or a
project `.claude/skills/`. Relevant here:

| Skill | Use |
|---|---|
| SigCom-LINCS | the reversal/mimicry engine above |
| Enrichr, ChEA3, KEA3 | pathway / TF / kinase enrichment on the derived signature |
| Rummagene | 793,703 gene sets mined from PMC, finds published BUB1B/aneuploidy sets |
| ARCHS4, GeneRanger | expression context for candidate targets |
| Geneshot, Harmonizome | literature and integrated-resource gene associations |
| DrugEnrichr, Perturb-Seqr | drug-set enrichment; perturbations with similar effects |

**PandaDock** (PyPI `pandadock` 4.1.1, Python ≥3.8, needs RDKit; torch/torch-geometric for the
GNN), narrow role, and only in the missense branch. Reported benchmarks: top-1 pose ≤2 Å in
30–34% of an 815-complex redocking set (best-of-9: 50–57%); affinity Pearson 0.531 pooled on
PDBbind v2020 but **within-target median R = 0.31**. The authors state the GNN should not be
used for pose selection.

So its value is as a second, physics-based opinion cross-checking Boltz's co-folding, not as
evidence on its own. Agreement between two structurally different methods is worth more than
either; a within-target R of 0.31 will not rank-order a candidate list on its own.

Neither Boltz nor PandaDock predicts the quantity that matters in the missense branch, the
**destabilisation of the variant protein itself** (ΔΔG of folding). **Gap closed 2026-09-03**
with ThermoMPNN + DDGun; see findings/N1002K_STABILITY.md.

BUB1B's nine PDB entries do **not** cover the pseudokinase domain; the structural arm rests on
the AlphaFold model (pLDDT 91.1 at residue 1002). Correction recorded in
findings/N1002K_STABILITY.md.
