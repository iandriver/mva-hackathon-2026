# ΔΔG of folding for BUBR1 p.Asn1002Lys

Run 2026-09-03. Closes the tool gap flagged in METHOD.md §6 and BUB1B_LOCUS.md caveat 3.

## Result

| Method | Input | ΔΔG | Reading |
|---|---|---|---|
| ThermoMPNN | AlphaFold, kinase domain 766–1050 | **+1.507** | destabilising ~1.5 kcal/mol |
| ThermoMPNN | AlphaFold, full-length 1–1050 | **+1.506** | destabilising ~1.5 kcal/mol |
| DDGun (sequence) | 661-sequence MSA | **−1.533** | destabilising ~1.5 kcal/mol |

The two methods use **opposite sign conventions**. ThermoMPNN positive = destabilising,
DDGun negative = destabilising. Both were checked internally: in each, proline substitution at
1002 is the most extreme of the 19, which is the expected behaviour and confirms the direction.

**p.Asn1002Lys is modestly destabilising, roughly 1.5 kcal/mol. It is not a catastrophic
substitution.**

## Calibration: the number alone means nothing

Site-saturation over all 285 kinase-domain positions, on the same run:

| | ThermoMPNN (domain) | DDGun |
|---|---|---|
| n substitutions | 5,415 | 5,377 |
| mean | 0.874 | −1.414 |
| median | 0.647 | −0.993 |
| **N1002K percentile** | **75.2%** | **62.9%** |

So N1002K is more destabilising than about two-thirds to three-quarters of substitutions in
this domain: above average, but well outside the extreme tail (ThermoMPNN 90th percentile is
2.29; N1002K is 1.51). Within its own position it ranks 8th of 19 possible substitutions.

Both predictors are skewed toward destabilising calls (ThermoMPNN 15.7% negative, DDGun 15.0%
positive), which is a known property of the method class. Quoting the raw value without the
percentile would overstate the finding.

## Two cross-checks

**Structural context.** Trimming to the kinase domain versus using the full-length model
changes the ThermoMPNN prediction by 0.001 kcal/mol (1.507 to 1.506). The low-confidence
disordered regions of the full-length model (whole-protein pLDDT 63.6) contribute nothing to
this site, so the trimming decision does not carry the result.

**Method independence.** ThermoMPNN is a supervised graph neural network transfer-learned from
the Tsuboyama megascale experimental dataset (Dieckhaus et al. 2024, PNAS,
doi:10.1073/pnas.2314853121). DDGun is untrained, an unsupervised combination of
substitution-matrix, hydrophobicity and conservation terms, so it cannot overfit a ΔΔG training
set. Different inductive biases converging on the same magnitude and direction is the
meaningful result; the agreement to three significant figures is coincidence.

## Orthogonal evidence from the MSA

From the same 661-sequence alignment used for DDGun, at residue 1002:

| Residue | Frequency |
|---|---|
| **N (wild type)** | **0.882** |
| gap | 0.115 |
| S | 0.003 |
| **K (the variant)** | **0.000** |

Asparagine is strongly conserved and lysine is never observed. This is evolutionary constraint,
independent of both ΔΔG predictors and of the clinical data.

## Why this matters for the repurposing strategy

The pseudokinase argument (BUB1B_LOCUS.md) says the damage from a kinase-domain missense in
BUBR1 should be conformational rather than catalytic, because human BUBR1's catalysis is
dispensable while its ATP-pocket residues are required for conformational stability
(Suijkerbuijk et al. 2012, doi:10.1016/j.devcel.2012.03.009). That argument predicted a
stability effect. A stability effect is what both methods find.

The *size* of the effect is the useful part. A ~1.5 kcal/mol destabilisation is:

- **large enough to matter**, sufficient to reduce steady-state cellular abundance through
  degradation of the partly-unfolded species;
- **small enough to be a plausible target**, recovering ~1.5 kcal/mol is within the range
  small-molecule stabilisers achieve, unlike a grossly misfolded protein.

This is the quantitative form of branch (iii). It supports the pharmacological-chaperone /
proteostasis hypothesis over both alternatives, restoring catalysis (there is none) and
aneuploidy-selective killing (wrong direction, MECHANISM.md §3).

It is also consistent with an attenuated presentation: a hypomorphic
allele retaining most of its fold, in trans with a null, predicts a milder presentation than
two severe alleles, toward the mild end of the published MVA1 spectrum.

## Correction to an earlier claim

MECHANISM.md §0 and METHOD.md §6 stated that BUB1B has nine experimental structures so the
structural arm needs no new structure determination. That is wrong for this residue.

Checking observed (not merely construct-mapped) residue ranges in every PDB entry for O60566:

| Entry | BubR1 chain | Observed residues |
|---|---|---|
| 6TLJ | S | 19–210, 221–233, 241–323, 341–345 |
| 5KHU | Q | 18–205, 224–230, 272–277, 291–308 |
| 5LCW | S | 1–560 (mapping) |
| 2WVI, 3SI5 | A | 57–220 |
| 5JJA | C, D | 647–720 |
| 5K6S | B | 663–681 |
| 5SWF | B | 668–676 |
| 4GGD | C, D | 20–42 |

**No experimental structure resolves the pseudokinase domain (766–1050).** 6TLJ and 5KHU map to
UniProt 1–1050 in PDBe's best-structures listing because that is the construct, but they model
nothing past residue 345, and they are 3.8 Å and 4.8 Å cryo-EM, unusable for side-chain
geometry regardless.

The structural arm therefore rests entirely on the AlphaFold model. That model is confident
where it needs to be, pLDDT 91.1 at residue 1002, 91.5 mean over 992–1012, 81.9 over the
kinase domain, but it is a prediction, and every downstream structural claim inherits that.

## Limits

1. Both predictors are trained or tuned on **small single-domain soluble proteins**. A
   pseudokinase domain excised from a 1050-residue multidomain protein is out of distribution.
2. **ΔΔG of folding is a proxy** for the quantity that actually matters, which is steady-state
   cellular abundance of BubR1 N1002K. The link between the two is empirical and imperfect.
3. Input is a **predicted structure**, not experimental (see the correction above).
4. No experimental validation of this variant exists. It remains a VUS with computational
   support, not a characterised hypomorph.
5. The ~1.5 kcal/mol figure is a point estimate with no uncertainty attached. Neither method
   returns a calibrated interval.

## Reproducing

```
tools/ddg-env          uv venv, ddgun 0.0.2 (+ biopython 1.81, setuptools<81)
tools/tmpnn-env        uv venv, torch 2.14, pytorch-lightning 2.6.5
tools/ThermoMPNN       cloned; local.yaml repointed to this repo, accel: cpu
structures/            AF-O60566 v6 full and 766-1050 trim, plus both SSM output CSVs
```
MSA: ColabFold MMseqs2 public API (`https://api.colabfold.com`), env mode, uniref a3m,
661 length-matched sequences after stripping a3m insertion columns.
