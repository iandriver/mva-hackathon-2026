# Both repurposing arms return negatives

Run 2026-09-03. Signature reversal (METHOD.md) and structural docking (BUB1B_LOCUS.md) were
run in parallel. Neither produces a candidate list. Both failures were caught by controls
built before the results were seen.

---

# Arm A: LINCS signature reversal

## The machinery works: positive control passes

Consensus signature `LINCS_L1000_CRISPRKO_BUB1B_mean_sig`, 12,327 genes, top/bottom 250 taken
as the up/down sets (250 and 248 mapped to SigCom entity IDs).

Two independent validations that the signature is real:

1. **BUB1B itself is in its own bottom-15 down genes.** The knocked-out gene is down, as it
   must be.
2. Querying the CRISPR library (`l1000_xpr`, ~70,000 signatures) returns **BUB1B knockout at
   rank 2** (z-sum 24.85, A549) and again at **rank 14** (23.17, HT29), recovered in cell
   lines other than those driving the consensus.

The signature is also mechanistically coherent. Top up-genes: **IFI6, ISG15, NNMT, CA2, CLU,
IRF7, IFI44L, ZNF91, CDKN1A, CXCL14, EIF5A, ATF3, MX1**. That is a type-I interferon module
plus p21 (CDKN1A). Top down-genes include **PCNA, HMGB2, CDCA4**, which are replication and
proliferation genes. Loss of a spindle checkpoint gene producing interferon activation, p53/p21
arrest and reduced proliferation is exactly what the micronucleus/cGAS-STING literature
predicts, and it matches the p53 activation reported in BubR1-deficient cortical progenitors
(Sterling et al. 2023, doi:10.3389/fcell.2023.1282182).

## The direction control fails

METHOD.md §3 specified the test: drugs that reverse the knockout signature and drugs that mimic
the overexpression signature should agree. The prerequisite is that **BUB1B overexpression
should itself reverse the BUB1B knockout signature.** Querying the overexpression library
(`l1000_oe`, 17,326 signatures) with the knockout signature:

| BUB1B OE signature | Side | Rank | z-sum |
|---|---|---|---|
| PC3, 96 h | reverser | 218 | −12.74 |
| A375, 24 h | reverser | 1,204 | −9.34 |
| A375, 96 h | **mimicker** | 993 | **+9.61** |
| HA1E, 96 h | **mimicker** | 1,205 | **+9.17** |

**Two reverse, two mimic**, including the same cell line (A375) landing on opposite sides at
24 h and 96 h. Overexpressing BUB1B does not reliably produce the transcriptional opposite of
knocking it out.

This is not noise in the query, the same pipeline recovered BUB1B knockout at rank 2 out of
70,000. The asymmetry is biological, and in hindsight it is expected: knocking out a spindle
checkpoint gene causes missegregation, micronuclei, interferon and arrest, while overexpressing
it does not un-missegregate anything. **Loss and gain are not mirror images, so "reverse the
loss signature" has no validated meaning for this target.**

## The compound list is confounded independently

Top-200 compound reversers (`l1000_cp`):

- **THP1 cells are 13.0% of the list versus 2.06% of the library, 6.3× enrichment.** A single
  cell line at a single timepoint (6 h) drives the top of the ranking.
- The top 25 are dominated by uncharacterised BRD-* screening compounds and STK33 tool
  compounds (ML-281 appears 5 times, MW-STK33-1C, BG-STK33-55).
- The few approved drugs present are not usable: **troglitazone** (withdrawn for
  hepatotoxicity), phenytoin, cefixime, pinacidil, emetine, a list with no mechanistic
  relationship to BUBR1 and every appearance of rank noise.

Either failure alone would be disqualifying. Together they mean **Arm A yields no candidate.**

---

# Arm B: docking the pseudokinase ATP pocket

The hypothesis (BUB1B_LOCUS.md) was that the ATP pocket of the BUBR1 pseudokinase is a
*stability* element (Suijkerbuijk et al. 2012, doi:10.1016/j.devcel.2012.03.009), so a
small molecule occupying it might act as a pharmacological chaperone for N1002K.

Two geometric checks on the AlphaFold model, before docking anything:

**1. The pocket is closed.** Probing an 8 Å sphere at the centroid of the annotated ATP-binding
residues (772–780, 795) on a 1 Å grid: **63 of 2,109 grid points (3.0%) are more than 3.3 Å
from any protein atom**, roughly 63 Å³ of probe-accessible volume. A druggable kinase ATP site
is typically 300–600 Å³. This is an occluded pocket, which is what a pseudokinase that does not
productively bind nucleotide should look like.

**2. The variant is nowhere near it.** Residue 1002 sits **33.3 Å** from the ATP-site centroid.
Even a ligand that did bind the pocket has no direct structural route to the N1002K position;
any stabilisation would have to be global or allosteric.

**Arm B is not a well-posed docking experiment on the available structure**, and 3,311 approved
small molecules were not docked. Producing scores against a 63 Å³ occlusion 33 Å from the
variant would have generated a ranked list with no meaning.

Caveat in the honest direction: AlphaFold returns one conformation, and pseudokinases are often
modelled closed. A cryptic pocket could open under dynamics. That is a reason to withhold the
claim that the site is undruggable, not a reason to dock into a closed one and report the
numbers. There is no experimental structure of this domain to check against
(N1002K_STABILITY.md, "Correction").

Noted separately: PandaDock's CLI banner claims "R > 0.8 correlation with experimental binding
affinities". Its own repository reports Pearson 0.531 pooled on PDBbind v2020 and **within-target
median R = 0.31**. The banner figure should not be quoted.

---

# What survives

The negatives are bounded, not vague:

| Claim | Bound |
|---|---|
| Signature reversal gives no candidate for BUB1B | Positive control recovered BUB1B KO at rank 2/70,000, so the failure is the loss/gain asymmetry, not the method's sensitivity |
| Compound reversers are context-confounded | THP1 6.3× over background, quantified against the full library |
| ATP pocket not dockable on this model | 63 Å³ accessible vs 300–600 Å³ typical; variant 33.3 Å away |

**The constructive result is in Arm A's signature, not its reversal.** The consensus BUB1B
knockout signature's strongest up-genes are a type-I interferon module. IFI6, ISG15, IRF7,
IFI44L, MX1, alongside CDKN1A. Chromosomal instability generating micronuclei, cytosolic DNA
and cGAS-STING-driven interferon is well established, and it is the one consequence of this
child's genotype that is both downstream of an undruggable protein and reachable with drugs
already approved for children (JAK inhibitors in type-I interferonopathies).

That targets the inflammatory consequence rather than the cause. It does not require killing
aneuploid cells, so it does not re-enter the direction problem (MECHANISM.md §3). It is a
hypothesis generated from one consensus signature and carries none of the validation done
above; it is the next thing to test, not a finding.

## Reproducing

Queries: SigCom-LINCS `enrich/ranktwosided` (POST needs `curl -L --post301 --post302 --post303`;
urllib drops the body on the 307). Databases `l1000_xpr`, `l1000_cp`, `l1000_oe`.
Consensus signature from `s3.amazonaws.com/lincs-dcic/LINCS-sigs-2021/means/xpr_mean_sigs/`.
Structure and pocket probe: `structures/BUBR1_AF_kinase_766_1050.pdb`.
