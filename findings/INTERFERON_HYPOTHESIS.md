# The interferon hypothesis, tested

Run 2026-09-03, same treatment as TWO_ARMS.md: build the controls first, then look.

Hypothesis: BUB1B loss → chromosome missegregation → micronuclei → cGAS-STING → type-I
interferon; the interferon module is druggable with JAK inhibitors already approved for
children, without killing aneuploid cells.

**Verdict: three of four checks pass. The fourth, on safety direction, fails, and that
failure is disqualifying for this child.** The hypothesis survives as mechanism and as a candidate
biomarker. It cannot be proposed as a therapy here.

## Check 1: Is the interferon module specific to BUB1B, or generic knockout stress? PASSES

Metric: mean percentile-rank of the 224 Hallmark Interferon Alpha + Gamma Response genes
within each consensus knockout signature (0.5 = no shift). Background: 300 randomly sampled
consensus CRISPR-KO signatures from `l1000_mean_xpr` (7,502 available), seed fixed.

| | ISG percentile-rank |
|---|---|
| background mean (n=300) | 0.4939 |
| background sd | 0.1045 |
| **BUB1B** | **0.7061, 98.3rd percentile, z = +2.03** |

The background sits at 0.4939 against a null of 0.5, so the metric is unbiased. Interferon
activation is a specific consequence of BUB1B knockout, not what any knockout does.

## Check 2: Does it track with the mechanism? PARTIAL, and the split is informative

If the driver is checkpoint weakening → missegregation → micronuclei, other spindle-checkpoint
knockouts should behave the same way. They do not behave uniformly:

| Gene | ISG rank | Percentile | z |
|---|---|---|---|
| **BUB1** | 0.7157 | 98.7 | +2.12 |
| BUB1B | 0.7061 | 98.3 | +2.03 |
| ZWINT | 0.6760 | 94.7 | +1.74 |
| INCENP | 0.6503 | 90.0 | +1.50 |
| TTK | 0.6054 | 82.3 | +1.07 |
| CENPE | 0.5723 | 75.7 | +0.75 |
| PLK1 | 0.5319 | 66.3 | +0.36 |
| AURKB | 0.5165 | 63.0 | +0.22 |
| KIF11 | 0.4864 | 51.3 | −0.07 |
| SPC25 | 0.4461 | 37.0 | −0.46 |
| NDC80 | 0.4023 | 19.0 | −0.88 |
| NUF2 | 0.3481 | 7.0 | −1.39 |

**BUB1, which is BUB1B's paralog and its partner in the mitotic checkpoint complex, is the
single highest of all twelve.** That is strong corroboration.

The split has a coherent reading: the high genes (BUB1B, BUB1, ZWINT, TTK) are checkpoint
*signalling*; the low ones (NDC80, NUF2, SPC25) are the NDC80 attachment complex. Weakening the
checkpoint lets cells through mitosis missegregating, producing viable aneuploid daughters that
form micronuclei. Destroying attachment holds the checkpoint on and arrests cells, which never
reach the state that generates interferon.

**That reading is post-hoc.** It was constructed after seeing the split, on 4 genes versus 3.

**It was then put to a prospective test, and the test could not be run.** The design
(`scripts/10_sac_partition_test.py`) pre-registers two groups from protein annotation before
scoring anything, excludes all twelve genes whose scores were seen while forming the
hypothesis, and asks whether group A exceeds group B by label permutation. The script aborts
if a generating gene leaks into the test set.

It aborts for a different reason. The LINCS consensus CRISPR-KO library covers 7,502 genes and
holds **none** of the twenty-one held-out candidates: not MAD1L1, MAD2L1, BUB3, TRIP13, ZW10,
ZWILCH, KNTC1, KNL1/CASC5, CENPF or SPDL1 on the signalling side, and not SPC24, DSN1, NSL1,
MIS12, PMF1, SKA1, SKA2 or SKA3 on the attachment side. The shRNA library covers a few
(BUB3, ZW10, TRIP13, CENPF, and PMF1 and SKA2) but that leaves group B with two genes, and
would require building consensus signatures from individual shRNA profiles, a different
perturbation modality with its own seed-effect confounding.

The twelve genes that generated the hypothesis are very nearly the only spindle-checkpoint
genes this resource contains. **The split therefore cannot be upgraded from post-hoc to tested
using this data, and stays labelled post-hoc.** Testing it needs a perturbation resource with
broader kinetochore coverage, or targeted experiments.

Checks 1 and 3 are unaffected: specificity and the JAK direction control stand on their own
data.

## Check 3: Do JAK inhibitors actually reverse it? PASSES, with measurable contamination

Restricting the query to the 68 ISGs among the 600 most-up genes in BUB1B KO (plus 118 down
genes) and searching `l1000_cp`:

| Side | JAK-family hits |
|---|---|
| **reversers** | tofacitinib (r82), TG-101348 (r93, r262), ruxolitinib (r95, r187, r275), baricitinib (r260, r307), **8 signatures, 4 drugs, 6 cell lines** |
| mimickers | AZD-1480 (r135, r346, both TMD8), S-ruxolitinib (r451), **3 signatures** |

**8 reversers against 3 mimickers, 73% directional consistency.** Compare the control that
killed Arm A: BUB1B overexpression split 2 reversers / 2 mimickers, 50%. This is a materially
different result, and the reversers span six cell lines rather than clustering in one.

Two honest deductions from it:
- The JAK inhibitors sit at ranks 82–307, not the top. The very top of this list is still the
  same BRD-* and STK33 tool compounds that contaminated Arm A. THP1 is 7.6% of the top 500 here
  versus 13.0% before, reduced, not eliminated.
- S-ruxolitinib appearing on the mimicker side is consistent with it being the less active
  enantiomer, but that is a convenient reading of a single signature and is not evidence.

Note also that the Hallmark ISG sets are by construction JAK-STAT-dependent transcripts, so
"a JAK inhibitor lowers ISGs" is close to tautological. The check is worth something because it
could have failed, as it did for BUB1B overexpression, but it is a weak test of the
therapeutic idea, not a strong one.

## Check 4: Which direction is the interferon pointing? FAILS

This is the same question that MECHANISM.md §3 raised about aneuploidy-selective killing, and
it has the same answer here.

The child has **already had a rhabdomyosarcoma**, and carries a constitutional
cancer-predisposition syndrome. Type-I interferon and JAK-STAT signalling are central to tumour
immunosurveillance. The interferon module measured above is a *response to* the chromosomal
instability that causes the cancer predisposition. Whether it is pathology or defence is not
established by anything measured here, and the therapeutic proposal requires it to be
pathology.

The clinical evidence points the wrong way. In ORAL Surveillance, a randomised trial of
tofacitinib versus TNF inhibitor in 4,362 patients, the hazard ratio for malignancies excluding
non-melanoma skin cancer was **1.48 (95% CI 1.04–2.09)**, noninferiority was not met, and the
risk diverged after month 18 (HR 0.93 before, 1.93 after; interaction p=0.0469)
(Ytterberg et al. 2022, NEJM, doi:10.1056/NEJMoa2109927; Curtis et al. 2022, Ann Rheum Dis,
doi:10.1136/ard-2022-222543). This is the basis of the class boxed warning.

An elevated malignancy signal in adults without a cancer-predisposition syndrome is a different
proposition from the same drug class in a child whose defining risk is malignancy. **Proposing
chronic JAK inhibition here means suppressing tumour immunosurveillance in a patient selected
for tumour susceptibility.** Absent evidence that the interferon is pathogenic rather than
protective, that is the wrong direction.

## What survives, and what would change it

Surviving:
- **A specific, quantified mechanistic finding**: BUB1B knockout activates interferon at the
  98.3rd percentile of 300 random knockouts, corroborated by its paralog BUB1 at 98.7th, and
  tracking with checkpoint-signalling rather than attachment genes.
- **A candidate biomarker**: an ISG panel is a measurable readout of chromosomal instability in
  this patient, independent of any therapeutic use.

Not surviving: JAK inhibition as a proposed treatment for this child.

What would change the verdict: evidence separating pathogenic from protective interferon in
MVA specifically, whether the chronic ISG elevation contributes to the growth failure and
inflammatory burden independently of its surveillance role. That is not answerable from a WGS
VCF and a public perturbation atlas. It needs patient-derived cells, which this hackathon does
not release.

## Method note

The Hallmark gene sets came from Enrichr (`MSigDB_Hallmark_2020`). Consensus knockout
signatures came from `s3.amazonaws.com/lincs-dcic/LINCS-sigs-2021/means/xpr_mean_sigs/`;
312 downloaded (BUB1B + 11 SAC/mitotic + 300 random). Enrichment queries used SigCom-LINCS
`enrich/ranktwosided` against `l1000_cp`.
