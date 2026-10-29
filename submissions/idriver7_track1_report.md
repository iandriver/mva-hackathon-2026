# Track 1: BUB1B compound heterozygote in WGS_EX2312012

Submitted by @idriver7. Reference GRCh38. Prediction: **one primary row**.

| | |
|---|---|
| Variant 1 | chr15:40209701 T>G, `ENST00000287598.11:c.2210T>G`, **p.Leu737Ter** |
| Variant 2 | chr15:40220612 T>G, `ENST00000287598.11:c.3006T>G`, **p.Asn1002Lys** |
| Gene | BUB1B (BUBR1), ENSG00000156970 |
| Disease | Mosaic variegated aneuploidy type 1, autosomal recessive |
| EPCR | 0.95 |

## 1. How the variants were reached

The released VCF is a single-sample GATK 4.2.4.0 callset, 5,012,204 records, reference
`GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_hs38d1`. **Contigs carry no `chr` prefix**;
the prefix is added when writing this submission, and its absence is a silent scoring failure
for anyone who does not check.

The proband phenotype, a childhood rhabdomyosarcoma with prenatal growth restriction, names a
childhood cancer-predisposition syndrome. In a recessive setting that combination points at the
spindle assembly checkpoint; biallelic BUB1B was established as the cause of MVA by Hanks et al.
2004 (doi:10.1038/ng1449). The challenge phenotype document is confidential, so this public copy
summarises it at the level the diagnosis already implies rather than listing it term by term.

Within the BUB1B gene span (15:40,160,984–40,221,137) the VCF holds 14 variants. Twelve are
intronic with gnomAD allele frequencies between 0.002 and 0.999. The remaining two are the
coding pair above.

This VCF retains `LowQual` calls rather than removing them: 271,414 of 5,012,204 records
carry a non-PASS FILTER, and one `LowQual;QD2` indel sits inside the BUB1B window at
15:40180642. A pipeline that filters to `PASS` by default discards those candidates without
reporting that it did. Both variants reported here are `PASS`.

## 2. Evidence

| | p.Leu737Ter | p.Asn1002Lys |
|---|---|---|
| Genotype | 0/1 | 0/1 |
| AD (ref,alt) / DP | 21,25 / 46 (VAF 0.54) | 15,13 / 28 (VAF 0.46) |
| QUAL / GQ | 708.77 / 99 | 344.77 / 99 |
| MQ, QD, FS, SOR | 60, 15.41, 1.115, 0.523 | 60, 12.31, 1.617, 0.330 |
| Exon | 17 / 23 | 23 / 23 |
| gnomAD v4 exome | 115 / 1,461,846 | **1 / 1,461,878** |
| gnomAD v4 genome | 5 / 152,174 | absent |
| ClinVar | **Pathogenic/Likely pathogenic**, multiple submitters, no conflicts | **no record for c.3006T>G** |
| CADD / REVEL / SIFT / PolyPhen | 38.0 / – / – / – | 26.4 / 0.472 / 0.00 / 0.997 |
| SpliceAI | 0.03 | 0.02 |

Both calls are clean at MQ 60, with no strand bias and balanced allele fractions consistent
with germline heterozygosity rather than mosaicism or artifact.

ClinVar has no entry for c.3006T**>G**.
It holds c.3006T**>A**, a different nucleotide producing the *same* substitution p.Asn1002Lys,
classified Uncertain significance. Asn is AAT/AAC and Lys is AAA/AAG, so both third-base
changes give Asn to Lys. Protein-level evidence transfers; the nucleotide-level record does not.
Describing this allele as "a ClinVar VUS" without that distinction is incorrect.

## 3. Interpretation: one null allele, one destabilised allele

UniProt O60566 is 1050 aa; residues 737 and 1002 were verified against the reference sequence.
The protein kinase domain spans 766–1050, with the active site at 882 and ATP-binding residues
at 772–780 and 795.

- **p.Leu737Ter truncates 29 residues before the kinase domain begins**, in exon 17 of 23. The
  transcript is NMD-competent; the expected outcome is a null allele.
- **p.Asn1002Lys sits inside the kinase domain**, 120 residues past the active site.

The domain is a **pseudokinase**: human BUBR1 retains the catalytic triad, but catalysis is
dispensable for error-free chromosome segregation, and the residues that contact ATP in
conventional kinases are instead required for conformational stability (Suijkerbuijk et al.
2012, doi:10.1016/j.devcel.2012.03.009). The domain is needed to promote KARD phosphorylation
and PP2A-B56 recruitment to the outer kinetochore (Gama Braga et al. 2020,
doi:10.1016/j.celrep.2020.108397). So a missense allele here should do **conformational**
damage, not catalytic.

We tested that prediction rather than asserting it. Two ΔΔG-of-folding predictors with
opposite sign conventions and unrelated inductive biases:

| Method | Input | ΔΔG | Percentile within the domain |
|---|---|---|---|
| ThermoMPNN (supervised, megascale) | AlphaFold, domain 766–1050 | +1.507 | 75.2% |
| ThermoMPNN | AlphaFold, full length | +1.506 | 72.1% |
| DDGun (untrained) | 661-sequence MSA | −1.533 | 62.9% |

Both put p.Asn1002Lys at **~1.5 kcal/mol destabilising**, above average for this domain but
outside the extreme tail (ThermoMPNN 90th percentile is 2.29). Trimming versus full-length
changes the prediction by 0.001 kcal/mol. Independently, in the same 661-sequence alignment,
asparagine is 88.2% conserved at residue 1002 and **lysine is never observed**.

A hypomorphic allele retaining most of its fold, in trans with a null, predicts a milder
presentation than two severe alleles, consistent with the attenuated end of the published
MVA1 phenotypic spectrum (García-Castillo et al. 2008, doi:10.1002/ajmg.a.32315).

No PDB entry resolves the pseudokinase domain. 6TLJ and 5KHU map to
UniProt 1–1050 in PDBe's best-structures listing because that is the construct, but they model
nothing past residue 345, at 3.8 Å and 4.8 Å. Every structural statement above rests on the
AlphaFold model, which is confident here (pLDDT 91.1 at residue 1002, 81.9 across the domain)
but is a prediction.

## 4. The defining phenotype is not visible in this genome

MVA is defined by mosaic variegated aneuploidy. We tested whether the released WGS shows it,
using 2,281,104 heterozygous PASS SNPs (depth 20–70) and a statistic robust to the
pericentromeric mapping artifacts that dominate a naive maximum: the median across windows of
2,000 het SNPs of mean |BAF − 0.5|.

| | |
|---|---|
| observed across 22 autosomes | 0.06495 (chr18) to 0.06881 (chr19) |
| between-chromosome spread | 0.00386 |
| whole-chromosome depth ratio | 1.000 – 1.045 |

A mosaic trisomy at fraction *f* splits BAF to (1+f)/(2+f) and 1/(2+f) and raises depth by
1+f/2. Calibrating by simulation against the observed depth distribution:

| mosaic fraction | shift in the statistic | exceeds the observed spread? |
|---|---|---|
| 0.06 | 0.00099 | no |
| 0.10 | 0.00349 | no |
| **0.12** | **0.00440** | **yes** |
| 0.20 | 0.01034 | yes |

**Detection floor: a mosaic whole-chromosome gain present in ≥12% of cells would move one
chromosome by more than the entire observed between-chromosome spread. No autosome does.**
The depth channel is consistent and slightly tighter (max ratio 1.045 to f ≤ 0.09), though that
figure includes GC coverage bias on chr16/19/22.

This is a bounded negative, not an absence of effort. It sits right at the level reported for
another MVA1 patient by karyotype, the 12% mosaic metaphases reported by Pavone et al. 2022
(doi:10.1007/s10072-022-06247-w). We can exclude mosaicism at that magnitude and not below it. The likely explanation is tissue: published MVA mosaicism is assayed in cultured
lymphocytes or fibroblasts, and aneuploid cells are selected against in the proliferating
lymphocyte pool that dominates bulk blood DNA.

**The diagnosis therefore rests on the genotype. The defining phenotype is not confirmable
from the released data**, and any pipeline that claimed to confirm MVA cytogenetically from
this VCF should be asked what its detection floor was.

## 5. Limits

1. **Phase is inferred, not demonstrated.** 10,911 bp apart, no parental data. *Trans* is
   parsimonious for a recessive disorder with this phenotype; this VCF cannot exclude *cis*.
2. **p.Asn1002Lys has no functional data.** It is essentially private and computationally
   predicted to be destabilising. It is not a characterised hypomorph.
3. **Secondary findings: surveyed, none reportable.** The whole VCF was intersected against
   all 346,384 ClinVar Pathogenic/Likely-pathogenic records (2026-08-29) and reportability
   judged against ACMG SF v3.2. Seven exact matches; **none in an SF gene**, and none is
   reported here. One ClinVar P/LP hit was **rejected**: PRSS1 7:142750561, which is in 20% of
   gnomAD exomes, sits at VAF 0.147 with MQRankSum -4.996, and lies inside the TRB locus.
   No structural or copy-number findings: this VCF contains no SV/CNV calls.
4. **The twelve intronic BUB1B variants were not run through splice prediction.** All are
   common in gnomAD, which is a frequency argument rather than a splice argument. Published
   MVA1 second alleles are repeatedly non-coding, including a 3.0 kb Alu/SVA insertion
   invisible to genome sequencing (Kawakami et al. 2023, doi:10.1038/s10038-023-01211-8), so
   this is where to look first if the trans hypothesis ever fails.
5. FASTQs were not processed; all calls are taken from the released VCF.

## 6. Tooling and AI disclosure

bcftools 1.x; Ensembl VEP REST; gnomAD v4 GraphQL; ClinVar E-utilities; UniProt REST;
AlphaFold DB (AF-O60566-F1-model_v6); ThermoMPNN (doi:10.1073/pnas.2314853121); DDGun 0.0.2;
ColabFold MMseqs2 public API for the MSA. Analysis code in the linked repository.

AI assistant: **Anthropic Claude, accessed through Claude Code on a Claude Max subscription
(consumer plan). "Help improve our AI models" is off, so inputs and outputs are not used to
train models.** Outputs were not rated at any point, per the organisers' caution that flagged
content can be used regardless of the account setting. No credit programme was used; if the
Hackathon's Claude credits are used for later work, the terms attached to them will be checked
separately rather than assumed to match this account.

Used for analysis orchestration, code, and drafting. Every numeric claim in this report was
produced by the code in the repository, not by the model. The model's structural claim about
available PDB coverage was wrong on first pass and was corrected against PDBe observed residue
ranges (see `findings/N1002K_STABILITY.md`).
