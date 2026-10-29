# Alignment-free copy-number profile of BUB1B

Run 2026-09-10. Reproduce with `../scripts/17_kmer_panel.py`,
`../scripts/18_kmer_count.py`, `../scripts/19_kmer_cnv_profile.py`.

[BUB1B_STRUCTURAL.md](BUB1B_STRUCTURAL.md) left an exon-level deletion untested:
21 of 23 exons have no depth information in the VCF, and structural variants were
never called. This profiles copy number directly from the reads.

## Method

The 8 FASTQs (84.67 GB) were streamed from HuggingFace, gunzipped in flight, and
matched against 58,266 unique 31-mers tiling BUB1B in 500 bp windows. Reads were
never written to disk, so the release never landed locally (89 GB free on the SSD
would not have held it). 1,071,835,104 reads processed, ~40x genome-wide depth.

Copy number is normalized to the gene's own diploid level (median per-k-mer depth
of the exon windows). It is CN 2 because the VCF carries normal-depth
heterozygous SNVs across the gene, which excludes a uniform loss.

## Result

Among the 23 assessable exon windows, per-k-mer depth is diploid throughout:
CN min 1.63, median 1.99, max 4.10. **No window sits at the ~1.0 a heterozygous
deletion would produce.** The nearest is exon 23 at CN 1.63, still well above
hemizygous.

21 of 23 exons have at least one assessable window and are diploid. Exons 4 and
16 are uninterpretable: their windows are repeat-dense (per-k-mer depth 30-130x
the diploid level), so unique-k-mer depth cannot be read there.

Detection floor: the 3-sigma Poisson CN resolution is 0.041, and the empirical
spread of assessable exon windows is roughly CN 1.6-2.3. A heterozygous deletion
(CN 2 to 1) falls far below that band, so any assessable exon carrying one would
read near CN 1.0. None does.

## Direction of the method's bias

Per-window depth is a mean over the window's k-mers, so a few residual
repeat-contaminated k-mers can inflate a window (three exon windows read CN
3.5-4.1; these are single windows amid diploid neighbours, i.e. contamination,
not focal duplication). Contamination can only push depth **up**. A deletion
removes copies of every unique k-mer and pulls depth **down**. The negative for
deletions is therefore robust to this bias; the method cannot reliably call
duplications, which is not the question here.

## Two normalization facts, stated because they shaped the numbers

- The external control regions built into the panel (four 50 kb chr15 windows)
  are repeat-contaminated: their k-mers, filtered only for uniqueness within
  their own region, recur genome-wide and read 96x the diploid level. They are
  unusable as a baseline and serve here only as the diagnostic that exposed the
  problem. Normalization is self-referential instead.
- Roughly half the intronic windows are repeat-dense for the same reason. The
  gene-wide mean copy number is therefore uninformative; only the per-window,
  repeat-filtered exon values carry the signal.

## What this closes and what it does not

Closes: a heterozygous exon-level deletion across 21 of 23 BUB1B exons is now
excluded, where the VCF had no depth to say so.

Does not close: exons 4 and 16 (repeat-dense, uninterpretable by this method); any
insertion, including a Kawakami-type Alu/SVA mobile element, which adds sequence
rather than removing unique k-mers; and any balanced event. Those remain the
alignment plus Manta/Delly/MELT job described in BUB1B_STRUCTURAL.md.

This does not bear on the compound-heterozygous call (Track 1, F-max 1.000). It
addresses the standing caveat that phase is undemonstrated and a third allele
was possible.

## Provenance

Reference hg38 chr15 (UCSC). k=31, 500 bp windows. Public reads, streamed. The
CN profile table is `../alphagenome_control/results/bub1b_kmer_cn_profile.csv`
(no genotypes). Per-file k-mer counts carry no genotypes but are kept in the
gitignored results/ as intermediates.
