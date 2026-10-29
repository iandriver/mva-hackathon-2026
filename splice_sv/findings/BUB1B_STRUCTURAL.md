# Structural variants in BUB1B: what this release can and cannot say

Run 2026-09-09. Reproduce with `scripts/16_bub1b_sv_assessment.py`.

[MECHANISM.md](../../track2/MECHANISM.md) §2 records that the Kawakami MVA1 second
allele was a 3.0 kb Alu/SVA insertion in intron 2, invisible to short-read genome
sequencing and found only by RNA-seq plus adaptive-sampling long reads.
[SPLICE_INTRONIC.md](SPLICE_INTRONIC.md) left structural variants as the untested
class. This is that test.

## Structural variants were never called

The VCF is GATK HaplotypeCaller output passed through VariantFiltration. It
declares no `SVTYPE`, `SVLEN`, `END`, `CIPOS`, `CIEND`, `IMPRECISE`, `MATEID` or
`CN` INFO field, defines no symbolic ALT alleles, and contains zero records with
a symbolic or breakend ALT across all 5,012,204 lines. There is no SV callset to
search.

The largest event the caller emitted anywhere in the genome is **1,170 bp**, with
1,095 indels over 100 bp out of 927,756 (0.12%). An insertion the size of the
Kawakami allele is above this caller's demonstrated range. That comparison is the
point: the one published MVA1 second allele of this class would not have appeared
in a file like this one.

## A large heterozygous deletion is mostly not excluded

A heterozygous deletion is visible in a small-variant callset without any SV
calling, because it halves read depth and removes heterozygous calls. Both tests
were run.

chr15 background depth: median 43, p10 25, p90 57 (n = 140,405). A heterozygous
deletion would put depth near 22.

Fourteen variants sit in the gene body, thirteen PASS with depth, median 40,
range 26–48.

- **Seven are heterozygous.** A heterozygous call cannot lie inside a
  heterozygous deletion, so each excludes one there regardless of depth.
- **Six are homozygous**, and could in principle be hemizygous. None are. The
  lowest sits at the 11th percentile of the chr15 depth distribution, which is
  ordinary coverage rather than half coverage.

So no deletion overlaps any of the thirteen informative positions. The problem is
how few positions that is.

| Window with no informative position | Span | Exons fully inside |
|---|---|---|
| 40168265–40181092 | 12,828 | 3 (exons 3–5) |
| 40198240–40209700 | 11,461 | 7 (exons 10–16) |
| 40160984–40168263 | 7,280 | 2 (exons 1–2) |
| 40209702–40216469 | 6,768 | 3 (exons 18–20) |
| 40187405–40192891 | 5,487 | none |
| 40193436–40198238 | 4,803 | 1 (exon 9) |
| 40216571–40220611 | 4,041 | 2 (exons 21–22) |
| 40182810–40186691 | 3,882 | 3 (exons 6–8) |
| 40181094–40182808 | 1,715 | none |

**Twenty-one of the twenty-three exons sit inside a window with no depth
information at all.** The VCF holds coverage only where it called a variant, and
in this gene that is thirteen points across 60 kb.

## Detection floor

**Excluded:** any heterozygous deletion overlapping one of the thirteen PASS
variant positions.

**Not excluded:** a deletion confined to any of the nine windows above, the
largest 12,828 bp and one of them containing seven consecutive exons. Also not
excluded: any insertion, any balanced event such as an inversion or
translocation, and any deletion the caller never attempted to find.

An exon-level deletion in BUB1B is therefore **untested**, not absent. Stating it
as absent would be inferring absence from a pipeline that never looked.

## What would resolve it

The release ships FASTQs, eight files across four lanes, about 85 GB, and no
BAM or CRAM. Any copy-number or structural analysis requires aligning from reads
first. Two routes, in increasing cost:

1. **Targeted depth.** Count BUB1B-specific k-mers across the FASTQs and compare
   coverage along the gene. Streams the reads once, no alignment, and fills the
   nine windows above with a continuous depth profile. Answers deletions, not
   insertions or balanced events. **Done — see [BUB1B_KMER_CN.md](BUB1B_KMER_CN.md):
   21 of 23 exons are diploid, no heterozygous exon-level deletion; exons 4 and 16
   are repeat-dense and remain uninterpretable.**
2. **Full alignment plus SV calling.** Align all four lanes, then run a
   split-read and discordant-pair caller (Delly) plus a mobile-element homology
   screen for the Alu/SVA/L1 class. **Done — see
   [BUB1B_SV_ALIGNMENT.md](BUB1B_SV_ALIGNMENT.md): no SV, no mobile-element
   insertion, uniform read depth (which also covers exons 4 and 16 that the k-mer
   method could not).**

Neither is blocked by the data-transfer agreement, since both run locally.

## Bearing on the case

This does not weaken the compound-heterozygous call. Two coding variants are in
hand, one a ClinVar pathogenic nonsense allele, and the Track 1 submission scored
F-max 1.000 against the confirmed answer. The relevant caveat is the one already
stated in [BUB1B_LOCUS.md](../../findings/BUB1B_LOCUS.md): phase is not demonstrated. A third
structural allele would matter for that question and for the several published
cases where exon screening found only one variant.
