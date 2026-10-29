# Track 2: Why MVA has no repurposing candidate, and what that rules in

Submitted by @idriver7. Proband PROBAND01, BUB1B compound heterozygote.

---

## Summary

This report does not propose a medication. It argues that no repurposed approved drug can
currently serve this child, gives a structural reason rather than an inventory of failures,
and states what a candidate would have to look like to survive.

Five candidate classes were assessed and refused, each against a control fixed before the
result was seen. The refusals are not independent. They share a cause:

> In mosaic variegated aneuploidy, the one level at which the interests of the child's somatic
> tissue and the interests of his cancer risk **align** is upstream, at the rate of chromosome
> missegregation. That is precisely the level at which the target is undruggable. Every level a
> drug can reach is a level at which those two interests **oppose**.

That is a sharper claim than "we found nothing", it is falsifiable, and it predicts which
future work can succeed.

---

## 1. What is disrupted

The proband is a compound heterozygote for `BUB1B`, which encodes BUBR1, a core component of
the spindle assembly checkpoint. Biallelic BUB1B causes mosaic variegated aneuploidy type 1
(Hanks et al. 2004, doi:10.1038/ng1449).

| Allele | Variant | Consequence |
|---|---|---|
| 1 | `c.2210T>G` p.Leu737Ter | nonsense, exon 17/23, ClinVar Pathogenic/Likely pathogenic |
| 2 | `c.3006T>G` p.Asn1002Lys | missense in the kinase domain, gnomAD 1 / 1,461,878 |

The two alleles are broken differently, and that determines everything downstream.

p.Leu737Ter truncates 29 residues before the kinase domain begins, in exon 17 of 23, 314
codons before the native stop. The transcript is NMD-competent, so this allele is
functionally null and its message is largely degraded.

p.Asn1002Lys sits inside the kinase domain at a buried position (relative solvent
accessibility 0.162, 8 neighbours within 8 A). Two ddG predictors with opposite sign
conventions and unrelated inductive biases agree it is destabilising by about 1.5 kcal/mol:
ThermoMPNN gives +1.507 on the domain and +1.506 on the full-length model, DDGun -1.533. Those
are the 75th and 63rd percentiles of all substitutions in this domain. Above average, and
outside the extreme tail. In the same 661-sequence alignment, asparagine is 88.2% conserved at
1002 and lysine never appears.

So the druggable allele is the hypomorph, not the null. One allele makes full-length
protein that is mildly destabilised; the other makes a message that is mostly degraded.

The damage is conformational rather than catalytic. Human BUBR1's kinase domain is a
pseudokinase: catalysis is dispensable for error-free chromosome segregation, while the
residues that contact ATP in conventional kinases are required for conformational stability
(Suijkerbuijk et al. 2012, doi:10.1016/j.devcel.2012.03.009). The domain is needed to promote
KARD phosphorylation and PP2A-B56 recruitment to the outer kinetochore (Gama Braga et al. 2020,
doi:10.1016/j.celrep.2020.108397). That argument predicted a stability effect before the ddG
was computed, and a stability effect is what both methods found.

---

## 2. The method: a gate that can refuse

Every candidate class was put through the same three questions, fixed in advance:

1. **Does the mechanism support it?** Not "is there a drug for this pathway".
2. **Which direction does it push, and is that the direction this patient needs?**
3. **What control would show me I am wrong, and does it pass?**

Question 2 is the one that does the work, and it is the one most easily skipped. A gene up in
disease is often compensatory rather than causative, and the two imply opposite drugs. In a
constitutional disorder whose normal tissue carries the disease phenotype, a drug validated in
cancer can be exactly inverted.

The gate is implemented, not just described. `scripts/` contains thirteen scripts that
reproduce every number here; three exist mainly to refuse conclusions.

---

## 3. What was refused, and why

### 3.1 Aneuploidy-selective compounds

The literature offers a ready-made answer. AICAR, 17-AAG and chloroquine are selectively lethal
to aneuploid cells (Tang, Williams, Siegel & Amon 2011, doi:10.1016/j.cell.2011.01.017). They
are the first hits for "aneuploidy drug" and they read like a repurposing result.

They are inverted here. A large fraction of this child's **normal somatic tissue** is
aneuploid, and the phenotype is already driven by excess apoptosis of aneuploid progenitors.
BubR1-deficient cortical progenitors die en masse, and deleting P53 gives only minimal rescue,
so a P53-independent death pathway also operates (Sterling et al. 2023,
doi:10.3389/fcell.2023.1282182). Giving an aneuploidy-selective cytotoxic systemically pushes
in the same direction as the disease.

### 3.2 LINCS signature reversal

Positive control passed. The consensus BUB1B knockout signature recovers **BUB1B knockout at
rank 2 of roughly 70,000** CRISPR signatures, in a cell line other than those driving the
consensus, and BUB1B appears in its own bottom-15 down genes.

Direction control failed. If reversal is meaningful, BUB1B **overexpression** should reverse
BUB1B knockout. Across 17,326 overexpression signatures, two reverse it and two mimic it, with
A375 on both sides at different timepoints. Loss and gain are not mirror images for this gene,
so "reverse the loss signature" has no validated meaning.

Independently, the compound hit list is confounded: THP1 cells are 13.0% of the top 200
reversers against a 2.06% library background, a 6.3-fold enrichment, and the top ranks are
uncharacterised BRD-* and STK33 tool compounds.

### 3.3 Structural docking of the pseudokinase pocket

Checked before docking anything. Probing an 8 A sphere at the annotated ATP site on a 1 A grid,
**63 of 2,109 points are more than 3.3 A from any protein atom**. That is roughly 63 A^3 of
accessible volume, against the 300 to 600 A^3 typical of a druggable kinase site. Residue 1002
sits **33.3 A** away. No ligands were docked. Scores against an occluded site 33 A from the variant would have
been a ranked list with no meaning.

### 3.4 Nonsense readthrough of the null allele

`c.2210T>G` turns codon 737 from TTA into **TGA**, the most readthrough-permissive stop, with G
at -1 and no C at -2. But +4 is A rather than C, and -3 and -9 are unfavourable. Tarraso et al.
tested a comparable TGA context in PYGM p.R50X with amlexanox, ataluren, RTC13 and G418 across
three cellular models and found no detectable readthrough (doi:10.1242/dmm.043281).

The larger obstacle is NMD. The transcript is degraded, so there is little template. The
standard fix is readthrough plus NMD inhibition, which leaves a repurposing brief entirely
since no NMD inhibitor is approved.

### 3.5 Chemical and pharmacological chaperones

Passes the first gate. Residue 1002 is buried, so the variant should make a quality-control
substrate rather than a tolerated surface change.

Fails the second. Every working pharmacological chaperone acts on a **secretory-pathway**
protein and binds a **defined site** on it. CFTR correctors bind MSD1 and NBD1
(doi:10.1007/s00018-021-03994-5); NB-DGJ occupies the beta-galactosidase active site
(doi:10.1007/s00439-020-02153-3); MC4R chaperones bind the receptor
(doi:10.1016/j.bbadis.2017.03.001). BUBR1 is cytosolic, so there is no ER bottleneck to
relieve, and section 3.3 already established there is nowhere to bind. Generic chemical
chaperones avoid the binding-site problem and pay for it; Huang et al. describe their
therapeutic application as very limited, non-specific and high-dose.

---

## 4. The structural result

### 4.1 Aneuploidy tolerance is closed in both directions

Tolerance is the property that lets a cell survive an aneuploid genome. Its mechanisms are
known: BCL9L loss reduces caspase-2 and prevents MDM2 and BID cleavage (Lopez-Garcia et al.
2017, doi:10.1016/j.ccell.2016.11.001); cyclin D1 overexpression overrides an intact p53/p21
checkpoint (Crockford et al. 2017, doi:10.1093/annonc/mdw612). Both describe tolerance as how
**cancers** survive aneuploidy.

In this child, tolerance is at once what his dying somatic progenitors need more of and what
his premalignant cells need less of. Same cells, same mechanism. Raising it treats the
developmental disease and accelerates the cancer risk; lowering it does the reverse. The two
effects cannot be separated by tissue, timing or dose.

### 4.2 The one aligned level is the undruggable one

| Level | Would it help? | Available? |
|---|---|---|
| **Upstream**: reduce the missegregation rate | Helps the phenotype **and** the cancer risk, because it prevents new aneuploidy rather than adjudicating between existing aneuploid cells | Requires BUBR1 function. No small-molecule site (3.3). All other modalities fail on biology, delivery or the approved-medicine brief (section 5) |
| **Downstream**: aneuploidy tolerance | Helps one and harms the other, whichever way it is set | Closed (4.1) |
| **Downstream**: the type-I interferon consequence | Would reduce inflammatory burden | Failed the safety gate (4.3) |

This is the structural reason MVA has no treatment beyond symptom management.

### 4.3 The interferon module, and why it still fails

The consensus BUB1B knockout signature's strongest up-genes are a type-I interferon module
(IFI6, ISG15, IRF7, IFI44L, MX1) alongside CDKN1A. That is specific rather than generic
knockout stress. Against 300 randomly sampled consensus knockouts, with a background mean of
0.4939 and sd 0.1045 against a null of 0.5, **BUB1B sits at the 98.3rd percentile, z = +2.03**.
Its paralog and checkpoint partner **BUB1 is the highest of twelve at 98.7**. JAK inhibitors
reverse the ISG-restricted signature consistently: 8 signatures across 4 drugs and 6 cell
lines on the reverser side against 3 on the mimicker side.

It fails the direction gate. This child **has already had a rhabdomyosarcoma** and carries a
constitutional cancer-predisposition syndrome. Type-I interferon and JAK-STAT signalling are
central to tumour immunosurveillance, and the module measured here is a response to the
instability that causes his cancer risk. Nothing establishes it as pathology rather than
defence. In ORAL Surveillance, tofacitinib versus a TNF inhibitor in 4,362 patients gave a malignancy
hazard ratio of **1.48 (95% CI 1.04 to 2.09)**, with noninferiority unmet
(doi:10.1056/NEJMoa2109927). That is in adults **without** a predisposition syndrome.

An earlier mechanistic reading of this result, that interferon tracks with checkpoint
signalling rather than kinetochore attachment, was generated post-hoc and **could not be
tested**: the consensus library holds none of the 21 held-out genes needed. It remains
labelled post-hoc rather than promoted.

---

## 5. Modalities, and what "undruggable" means here

"No small-molecule pocket" does not mean "undruggable". Worked through separately:

Excluded on biology. Small molecules (no site). Antibodies (BUBR1 is cytosolic; no approved
intracellular antibody modality). siRNA and degraders (reduce protein; wrong direction).
**Base editing**: both alleles are T>G, so correcting either needs a **G>T transversion**, and
base editors perform transitions only. Neither allele is base-editable.

Excluded on delivery, not persistence. Persistence is solvable, since lentivirus and
transposons integrate and editing is permanent. The binding constraints are breadth and timing.
Aneuploidy arises at every mitosis, so correction is needed across neural progenitors, growth
plate, kidney, muscle and marrow at once. LNPs are liver-dominant through ApoE tropism
(doi:10.1016/j.jconrel.2024.05.036). Retargeting is real but redirects to one organ at a time
in mice, reaching 97 to 99% spleen selectivity in the best case (doi:10.1002/smll.202503034).
Single-organ redirection does not compose into whole-body coverage. Separately, growth
restriction in MVA begins before birth, and this child's presentation includes it, so much of the
developmental phenotype predates any possible intervention. A gene-therapy proposal must state
whether it targets ongoing aneuploidy and cancer risk or claims to reverse what is done.

Live on biology, outside the brief. Suppressor tRNA therapy is the cleanest mechanistic
answer to a TGA nonsense allele and bypasses the poor readthrough context entirely. Prime
editing is the only chemistry that can make the required transversion. Neither is an approved
medicine.

Narrower than a cure, and better than we first argued. Autologous ex vivo HSC gene
correction addresses the MDS and leukaemia risk specifically. The MVA1 transplant that ended in
graft rejection and death used a matched unrelated donor (doi:10.1186/s13023-019-1073-x);
autologous corrected cells carry no allogeneic mismatch and would not face that failure mode.

---

## 6. What would change the answer

The thesis is falsifiable. Any of these would break it:

1. **A cryptic pocket on BUBR1.** The pocket analysis used a single AlphaFold conformation, and
   pseudokinases are often modelled closed. Molecular dynamics or a fragment screen finding an
   openable site would reopen the upstream level, which is the aligned one.
2. **Evidence that the interferon module is pathogenic rather than protective in MVA.** That
   would reopen 4.3. It needs patient-derived cells, which this hackathon does not release.
3. **A tissue-selective way to set the tolerance dial** differently in premalignant and
   somatic cells. Section 4.1 asserts they cannot be separated; a marker that distinguishes
   them would falsify it.
4. **Functional data on p.Asn1002Lys.** Everything here treats it as a hypomorph on
   computational grounds. If the protein turns out to be null, the missense allele stops being
   the druggable one and section 1 changes.

---

## 7. What is actionable now, without a drug

Two outputs bear on care, and neither requires a medicine.

A secondary-finding survey. An unbiased scan of the whole VCF against all 346,384 ClinVar
Pathogenic/Likely-pathogenic records returned seven variants. **None is in an ACMG SF v3.2
gene, and none is reported here.**

A separate ClinVar Pathogenic/Likely-pathogenic hit in PRSS1 was **rejected**: it is present in
20.2% of gnomAD exomes, sits at VAF 0.147, carries MQRankSum -4.996, and lies inside the
somatically recombining TRB locus.

A tissue result. The released genome does not show the disease's defining phenotype.
Using 2,281,104 heterozygous PASS SNPs and a statistic robust to pericentromeric mapping
artifact, the 22 autosomes span 0.06495 to 0.06881, a spread of 0.00386. Simulation against the
observed depth distribution puts a 12% mosaic trisomy at a shift of 0.00440.
**Mosaic whole-chromosome gain in 12% or more of cells is excluded**; the depth channel agrees
at f <= 0.09. That floor lands on the 12% mosaic metaphases reported by karyotype in another
MVA1 patient (doi:10.1007/s10072-022-06247-w). The likely explanation is tissue: mosaicism is
assayed in cultured lymphocytes or fibroblasts, and aneuploid cells are selected against in the
proliferating lymphocyte pool that dominates bulk blood. **Future functional work on this child
should not use bulk blood.**

---

## 8. Limits

1. No candidate medication is proposed. That is the finding, and it does not answer the
   question Track 2 asked.
2. The closed-window argument (4.1) is reasoned from published mechanism, not demonstrated
   experimentally. I did not find it stated in the MVA or aneuploidy literature, but the
   searches were targeted rather than exhaustive.
3. p.Asn1002Lys has no functional data. The 1.5 kcal/mol figure is a point estimate from two
   predictors with no calibrated interval, computed on a predicted structure. **No experimental
   structure resolves the BUBR1 pseudokinase domain**: 6TLJ and 5KHU map to UniProt 1-1050
   because that is the construct, but model nothing past residue 345, at 3.8 and 4.8 A.
4. Phase is not demonstrated, but a hidden third allele is now excluded. The two coding
   variants are 10,911 bp apart with no parental data released, so *trans* configuration stays
   inferred from the recessive phenotype rather than shown. The sharper concern was that the
   real second hit might be a cryptic non-coding or structural allele, the recurring pattern in
   published MVA1 cases: one Kawakami allele was a 3.0 kb Alu/SVA insertion invisible to genome
   sequencing. That was screened directly and found absent. SpliceAI, calibrated on BUB1B to
   recover 3/3 published pathogenic MVA1 splice alleles at its recommended cutoff, scores every
   intronic variant negative (though it is only validated near the splice site, and the
   proband's intronic variants sit 905-3,653 bp deep). A k-mer copy-number profile over all
   1.07 billion reads, and a bwa alignment with Delly plus an Alu/L1/SVA homology screen, find
   no exon-level deletion, no structural variant and no mobile-element insertion, with read
   depth uniform across all 23 exons. Every variant class short of an event contained within a
   single read is excluded; indels up to 1,170 bp are already in the VCF. Detail and detection
   floors are in `splice_sv/findings/` (SPLICE_INTRONIC, BUB1B_KMER_CN, BUB1B_SV_ALIGNMENT).
5. Both ddG predictors are trained on small single-domain soluble proteins. An excised
   pseudokinase domain from a 1,050-residue multidomain protein is out of distribution.
6. LINCS is largely cancer cell lines at 6 and 24 hours. The neural progenitor context most
   relevant here is represented by three signatures.

---

## 9. Reproducibility and disclosure

The thirteen core analysis scripts are in `scripts/`, with run order, runtimes, fixed seeds
and which need gated data documented in `scripts/README.md`. Scripts 06 to 13 need no patient
data and can be run by anyone. The structural follow-up that closes Limit 4 is a separate
workstream in `splice_sv/` (scripts 14-22 plus the AlphaGenome positive control), documented in
`splice_sv/README.md`. Gated data never enters the repository; a missing input names its fetch
command rather than raising, and every genotype-bearing output stays gitignored.

Public resources used: Ensembl REST and VEP REST, gnomAD v4 GraphQL, ClinVar, UniProt,
AlphaFold DB, PDBe, SigCom-LINCS, Enrichr, ColabFold MMseqs2, ChEMBL, ClinicalTrials.gov,
PubMed, UCSC RepeatMasker, Dfam (Alu/L1/SVA consensus), and the hg38 no-alt analysis set.
Tools: bcftools, ThermoMPNN (doi:10.1073/pnas.2314853121), DDGun, Biopython, SpliceAI 1.3.1
(doi:10.1016/j.cell.2018.12.015; non-commercial), AlphaGenome (non-commercial, positive control
only, not run on proband data), bwa, samtools, Delly v2.6.0. MELT was not used: it requires a
registered download, and a bwa homology screen against the Alu/L1/SVA consensus covers the same
mobile-element signal.

AI assistant. Provider: Anthropic. Product: Claude (Opus-class models) via Claude Code. Tier:
a Claude Max subscription (Anthropic consumer tier). Settings: the model-training control ("Help
improve our AI models") was off for the account, and outputs were never rated, because
user-submitted Feedback is a stated exception under which content is used for training
regardless of the setting. No credit programme was used.

Applicable data-handling terms, documented per the organisers' response in Community discussion
#21, which noted the setting alone is not sufficient:

- **Ownership.** Under the Anthropic Consumer Terms (effective 2025-10-08), the user retains
  right, title and interest in Inputs, and Anthropic assigns its right, title and interest in
  Outputs to the user.
- **Training.** Consumer content is not used to train Anthropic's models unless the user opts
  in; with the control off, it is not. Two exceptions apply regardless: content submitted as
  Feedback (avoided here), and content flagged for safety review.
- **Retention.** With training off, deleted conversations are removed from back-end storage
  within 30 days (consumer retention policy effective 2026-07-01). Safety-flagged content is
  retained up to 2 years (inputs and outputs) and 7 years (classification scores).
- **Operational licence.** A licence for Anthropic to "provide, maintain, and improve the
  Services" persists independent of the training opt-out. It is a service-provision licence, not
  an assignment of rights; ownership remains as above.

Two conditions lie outside the user's control and are stated rather than claimed resolved:
safety-flagged content may be retained and used regardless of the setting, and the operational
"improve the Services" licence is broader than the model-training opt-out. Whether these
consumer-tier terms meet the requirement that the provider not train on or take rights in
inputs and outputs, with retention limited in time and purpose, is left for the organisers to
assess on the terms above; their guidance notes that an API or commercial tier carries explicit
data-handling terms should the consumer terms not be verifiable.

Every numeric claim here was produced by the code in the repository, not by the model. Two
model-generated claims were wrong on first pass and are corrected in the record: that BUB1B's
PDB entries cover the pseudokinase domain (they do not), and that episomal dilution closes gene
delivery (it closes AAV).
