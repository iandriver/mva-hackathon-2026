# The BUB1B locus in WGS_EX2312012

Analysed 2026-09-03 from `WGS_EX2312012_HGWCNDSX7.vcf.gz` (GRCh38, 5,012,204 records).
Gene span from Ensembl: 15:40,160,984–40,221,137, + strand, ENSG00000156970,
canonical ENST00000287598.11. Fourteen variants fall in the gene body; twelve are intronic
with gnomAD AF between 0.002 and 0.999 and are not considered further here.

## The two coding variants

| | Variant 1 | Variant 2 |
|---|---|---|
| Position (GRCh38) | 15:40209701 T>G | 15:40220612 T>G |
| HGVS c. | ENST00000287598.11:c.2210T>G | ENST00000287598.11:c.3006T>G |
| HGVS p. | p.Leu737Ter | p.Asn1002Lys |
| Consequence | stop_gained (HIGH) | missense (MODERATE) |
| Exon | 17 / 23 | 23 / 23 |
| Genotype | 0/1 | 0/1 |
| AD (ref,alt) / DP | 21,25 / 46 (VAF 0.54) | 15,13 / 28 (VAF 0.46) |
| QUAL / GQ | 708.77 / 99 | 344.77 / 99 |
| MQ, QD, FS, SOR | 60, 15.41, 1.115, 0.523 | 60, 12.31, 1.617, 0.330 |
| dbSNP | rs759242053 | none |
| gnomAD v4 exome | 115 / 1,461,846 (7.87e-05) | **1 / 1,461,878 (6.84e-07)** |
| gnomAD v4 genome | 5 / 152,174 (3.29e-05) | absent |
| ClinVar | **Pathogenic/Likely pathogenic**, multiple submitters, no conflicts, 2024-10-09 | no record for c.3006T>G |

Both calls are clean: MQ 60, no strand bias, balanced allele fractions consistent with germline
heterozygosity rather than mosaicism or artifact. Neither is filtered.

### A nucleotide-level subtlety on variant 2

ClinVar has no record for **c.3006T>G**. It does hold **c.3006T>A**, a different nucleotide
change producing the *same* protein consequence p.Asn1002Lys, classified **Uncertain
significance** (single submitter, 2025-09-19). Asn = AAT/AAC; Lys = AAA/AAG; both T>G and T>A
at the third base give Asn to Lys.

So protein-level evidence transfers between the two records and nucleotide-level evidence does
not. The honest statement is that this child's missense allele is *absent from ClinVar*, and
the closest record, same amino acid substitution, different nucleotide, is a VUS. Reporting
it as "ClinVar VUS" without that distinction would be wrong.

## Position within the protein

UniProt O60566, 1050 aa. Verified against the reference sequence: residue 737 is L, residue
1002 is N.

| Domain | Residues |
|---|---|
| BUB1 N-terminal | 62–226 |
| KNL1-interaction region | 152–185 |
| **Protein kinase** | **766–1050** |
| Active site (proton acceptor) | 882 |
| ATP-binding site | 772–780, 795 |

- **p.Leu737Ter truncates 29 residues before the kinase domain begins.** Exon 17 of 23, so the
  transcript is NMD-competent; the expected outcome is a null allele.
- **p.Asn1002Lys sits inside the kinase domain**, 120 residues past the active site and 48 from
  the C-terminus. It is not at the active site or the annotated ATP-binding residues.

## This resolves the allele-architecture fork

MECHANISM.md §2 set out three branches. The genotype selects **branch (iii): missense with
residual protein**, one null allele plus one full-length allele carrying a kinase-domain
substitution. Not two nulls. Consequences:

- Dose restoration and fold stabilisation are coherent goals rather than speculation.
- The structural arm (Boltz, cross-checked by PandaDock) now has a defined input.
- The SigCom-LINCS overexpression-mimicry control (METHOD.md §3) is directly meaningful,
  because there is a protein whose abundance or stability could be raised.

## What the pseudokinase literature does to the interpretation

The BUBR1 C-terminal kinase domain is a **pseudokinase**, and this is not incidental:

- Human BUBR1 retains the catalytic triad, but **putative catalysis is dispensable for
  error-free chromosome segregation**. Instead, residues that contact ATP in conventional
  kinases are **essential for conformational stability** of BUBR1.
  Suijkerbuijk et al. 2012, Dev Cell. doi:10.1016/j.devcel.2012.03.009
- The pseudokinase domain is required to promote KARD phosphorylation; mutating or removing it
  reduces PP2A-B56 recruitment to the outer kinetochore, attenuates checkpoint silencing, and
  causes chromosome alignment errors through Aurora B imbalance.
  Gama Braga et al. 2020, Cell Rep. doi:10.1016/j.celrep.2020.108397
- KARD phosphorylation by PLK1 drives the BUBR1–PP2A-B56α interaction that opposes Aurora B at
  kinetochores. Suijkerbuijk et al. 2012, Dev Cell. doi:10.1016/j.devcel.2012.09.005

The consequence for repurposing: the damage done by a kinase-domain missense allele in BUBR1
is **most likely conformational, not catalytic**. There is no catalytic activity to restore.
That rules out one class of intervention and points at another, stabilising the fold, or
raising the abundance of the variant protein, which is a pharmacological-chaperone /
proteostasis hypothesis rather than a kinase-modulator one.

It also sharpens the structural question into something specific and testable: the ATP pocket
of this pseudokinase is a *stability* element (Suijkerbuijk 2012), so whether an approved drug
can occupy that degenerate pocket and stabilise the N1002K protein is a well-posed question for
the structural arm, with nine experimental structures plus AlphaFold available.

## Caveats, stated rather than assumed

1. **Phase is not demonstrated.** The two variants are 10,911 bp apart, far beyond short-read
   or read-pair phasing distance, and no parental data is in the release. *Trans*
   configuration is the parsimonious reading given a recessive disorder and the phenotype, but
   this VCF alone cannot exclude *cis*. Resolving it needs parental genotypes or long reads.
2. **p.Asn1002Lys has no functional data.** It is essentially private (1 allele in ~1.46M),
   which supports PM2, and the same substitution is a VUS in ClinVar. Its predicted
   destabilising effect is a prediction and has not been computed yet.
3. **Neither Boltz nor PandaDock predicts ΔΔG of folding**, the quantity that matters most for
   claim (§ pseudokinase). This gap was flagged in METHOD.md §6 and is now load-bearing rather
   than hypothetical.
4. **The intronic variants were evaluated for splice impact and the answer is inconclusive.**
   See [SPLICE_INTRONIC.md](../splice_sv/findings/SPLICE_INTRONIC.md). SpliceAI scores all fourteen gene-body
   variants below DS 0.03, but a saturation scan shows that at eleven of them no possible
   single-nucleotide change within ±50 bp reaches DS 0.2, so those zeros carry no information.
   The two *coding* variants are genuine negatives, since the model reads splice signal at
   those positions (neighbourhood maxima 0.76 and 0.92). The twelve intronic variants sit 905
   to 3653 bp from an exon and remain unadjudicated; the frequency argument is all that stands.
   Published MVA1 second alleles are repeatedly non-coding (MECHANISM.md §2), but the ones with
   experimental support are *near-splice*, and one was a 3.0 kb mobile element insertion no
   SNV-based predictor can see.
