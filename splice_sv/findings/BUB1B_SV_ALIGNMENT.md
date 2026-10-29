# Alignment-based structural-variant and mobile-element screen of BUB1B

Run 2026-09-10/11. Reproduce with `../scripts/20_targeted_align.sh`,
`../scripts/21_delly_call.sh`, `../scripts/22_me_screen.py`.

Targets the classes [BUB1B_KMER_CN.md](BUB1B_KMER_CN.md) left open: an insertion
(e.g. a Kawakami-type Alu/SVA mobile element) or a balanced event, neither
visible to a small-variant VCF or to unique-k-mer depth.

## Method

All four lanes were aligned to full hg38 with bwa mem; only read pairs anchored in
BUB1B +/- 5 kb were kept (26,615 records, ~6,600/lane). Full genome is required so
paralog and mobile-element mates map to their true copies rather than piling onto
BUB1B. Two callers: Delly (split-read + discordant-pair SV calling) and a
transparent screen for soft-clip breakpoints, discordant/one-end-anchored
clusters, and homology of breakpoint reads to Alu/L1/SVA consensus. MELT was not
used: it requires a registered download.

## Result: no structural variant or mobile-element insertion

**Delly: 0 calls** in or near BUB1B (deletions, insertions, duplications,
inversions, translocations).

**Coverage is uniform** across the gene: median 187 reads per 500 bp window, the
only sub-threshold windows being the two flank edges outside the gene body. This
cross-checks the k-mer copy-number profile and extends it: exons 4 and 16, which
were repeat-dense and uninterpretable by unique-k-mer depth, carry normal read
depth here (194 and 236 reads), so a deletion is now excluded there too.

**The raw screen flags are all native reference repeats, not a novel event.**
BUB1B contains 107 reference SINE/LINE elements. Every flagged signal coincides
with one:

| Signal | Position | Native element there |
|---|---|---|
| soft-clip cluster (22 reads) | chr15:40214833 | DNA/MER20B |
| soft-clip cluster (16 reads) | chr15:40180667 | SINE/AluY |
| soft-clip cluster (12 reads) | chr15:40193492 | SINE/AluJr |
| AluY-homologous reads (35) | six positions | all inside native AluSz/Sp/Y/Sq |
| L1/SVA-homologous reads (4) | four positions | inside native L1MEf / AluY |

All three soft-clip clusters are **one-sided** (right-clip only, zero left-clip),
carry **no polyA/T tail**, and have negligible discordant support (0–4 reads). A
polymorphic mobile-element insertion produces the opposite: clips on **both**
sides converging at one breakpoint, a polyA tail, and a stack of discordant reads
whose mates map to dispersed element copies. None of the three has that
signature. They are the expected mapping artifact where reads extend past a
reference-repeat boundary.

The 180 discordant reads across the 70 kb are background for an intron region this
repeat-rich; the two weak 5-read clusters sit at native repeats and are not
matched by a coincident breakpoint.

## Detection floor

A non-reference mobile-element insertion larger than ~300 bp (Alu and up) would
appear as a two-sided soft-clip pileup at one non-reference position, coincident
discordant reads, and ME-homologous clips, at >= 5 reads given the ~40x depth
here. None is present. A balanced translocation would appear as a discordant
cluster to another chromosome; none is. Not resolvable by this method: an event
fully contained within a single read, and small indels, which are already in the
VCF (largest indel called genome-wide was 1,170 bp).

## Bearing on the case

Four independent lines now say the BUB1B locus carries no third structural allele:
no SV called, uniform read depth, no k-mer copy-number change (21/23 exons, the
other two covered by read depth here), and no insertion breakpoint. This does not
change the compound-heterozygous call (Track 1, F-max 1.000); it closes the
standing caveat that phase is undemonstrated and a third allele was possible, for
every variant class short of a within-read event.

## Provenance

hg38 no-alt analysis set (chr naming). bwa mem, samtools, Delly v2.6.0 (bioconda
osx-64 under Rosetta). Alu/L1/SVA consensus from Dfam. Public reads, streamed and
deleted per lane. The screen summary is
`../alphagenome_control/results/bub1b_sv_screen.csv`; the BAM (proband reads) and
hg38 index stay in the gitignored results/ and on the SSD.
