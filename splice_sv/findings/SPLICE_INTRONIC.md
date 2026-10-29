# Splice evaluation of the BUB1B gene-body variants

Run 2026-09-08. SpliceAI 1.3.1, five-model ensemble, hg38, `-D 50`, GRCh38
annotation, run locally. No proband data left this machine. Reproduce with
`scripts/14_splice_intronic.py` then `scripts/15_splice_detection_floor.py`.

Addresses caveat 4 of [BUB1B_LOCUS.md](../../findings/BUB1B_LOCUS.md): the intronic variants
were never evaluated for splice impact, and [MECHANISM.md](../../track2/MECHANISM.md) §2
records that published MVA1 second alleles are repeatedly non-coding and
splice-disrupting.

## Calibration first

SpliceAI was validated on this gene before being pointed at the proband, on a
582-variant panel of public gnomAD and literature variants
(`scripts/calibration/bub1b_spliceai_calibration.csv`). At its recommended
DS ≥ 0.5 it recovered 3/3 published pathogenic MVA1 splice alleles (Rio Frio
c.2386-11A>G, Kawakami c.1402-5A>G, Laberko c.1288+5G>A) while passing 6.9% of
background and 0% of common and deep intronic variants. Its directional calls
reproduced the published mechanisms unprompted, including the acceptor
gain plus acceptor loss signature of the Rio Frio allele.

**That validation covers near-splice variants only.** All three positives sit
within 11 bp of an exon boundary. It licenses no claim about deep intronic
positions.

## Result

Fourteen variants in the gene body: 9 SNVs and 5 indels, 2 exonic and 12
intronic. All fourteen were scored, and the script asserts the count so a
partial extraction fails rather than under-reporting.

**No variant reaches DS 0.5. The maximum observed is 0.03.**

| Zone | n | Max DS |
|---|---|---|
| exonic | 2 | 0.03 |
| intronic, 905–3653 bp from the nearest exon | 12 | 0.02 |

## The negative is mostly not a negative

A DS of 0.00 means something only if SpliceAI could have produced a higher
number at that position. `scripts/15_splice_detection_floor.py` scores every
possible SNV within ±50 bp of each variant, 3,936 in total, and asks what the
best attainable score in each neighbourhood is.

| Centre | Zone | Observed | Neighbourhood max | n ≥ 0.5 | Verdict |
|---|---|---|---|---|---|
| 40220612 | exonic (p.Asn1002Lys) | 0.02 | **0.92** | 8 | informative |
| 40209701 | exonic (p.Leu737Ter) | 0.03 | **0.76** | 5 | informative |
| 40187404 | donor +1762 | 0.02 | 0.37 | 0 | weak |
| 40168264 | acceptor −1798 | 0.01 | 0.09 | 0 | blind |
| 40198239 | acceptor −1376 | 0.01 | 0.05 | 0 | blind |
| 40216568 | acceptor −928 | 0.00 | 0.04 | 0 | blind |
| 40186692 | donor +1050 | 0.00 | 0.06 | 0 | blind |
| 40180642 | acceptor −3072 | 0.00 | 0.02 | 0 | blind |
| 40181093 | acceptor −2621 | 0.00 | 0.01 | 0 | blind |
| 40192892 | acceptor −3653 | 0.00 | 0.01 | 0 | blind |
| 40182809 | acceptor −905 | 0.00 | 0.00 | 0 | blind |
| 40193435 | acceptor −3110 | 0.00 | 0.00 | 0 | blind |
| 40216470 | acceptor −1026 | 0.00 | 0.00 | 0 | blind |
| 40216570 | acceptor −926 | 0.00 | 0.00 | 0 | blind |

Across the whole scan, 13 of 3,936 possible SNVs reach DS ≥ 0.5, and all 13 sit
in the two exonic neighbourhoods.

**Two findings and eleven non-results.**

1. **Neither coding variant disrupts splicing.** At both exonic positions the
   model demonstrably reads splice signal, since single base changes nearby
   reach 0.76 and 0.92. The observed 0.02 and 0.03 are real negatives. The
   p.Asn1002Lys allele in particular acts through the protein, not through the
   transcript, which is what the pseudokinase stability argument in
   BUB1B_LOCUS.md assumes.

2. **The twelve intronic variants are unevaluated, not cleared.** At eleven of
   them no possible single-nucleotide change within ±50 bp reaches DS 0.2. A
   score of 0.00 there is uninformative.

This experiment cannot separate two explanations for the blind windows: SpliceAI
may be insensitive in deep intron, or the sequence there may be too far from any
splice-like motif for one base change to create a site. Either way the
conclusion for this case is the same, and the opposite of what the raw scores
suggest on their own. Reporting "twelve intronic variants show no splice
disruption" would have been a plausible wrong number.

## What caveat 4 now says

It is narrowed, not closed. The frequency argument is unchanged, all twelve are
common enough in gnomAD to be unlikely drivers. What is added is that the
standard splice predictor cannot adjudicate them, because they all sit 905 to
3653 bp from an exon and the tool was only shown to work within 11 bp.

The published MVA1 second alleles that motivated the concern were themselves
near-splice, at −5 and −11, so a deep intronic second allele in this child is a
different hypothesis from the literature precedent rather than the same one. The
Kawakami case is the relevant warning: its second allele was a 3.0 kb Alu/SVA
insertion in intron 2, invisible to short-read sequencing and to any SNV-based
predictor including this one.

## Open

- AlphaGenome uses 100 kb of context and may behave differently in deep intron
  where SpliceAI is flat. It is API-only, so it remains blocked on the
  data-transfer and HIPAA questions with the organizers.
- Structural variants and mobile element insertions are not addressed by
  anything run so far and are the class that defeated exon screening in the
  published cases.
- Branch point and splicing-regulatory-element disruption are not what SpliceAI
  scores.

## Provenance

SpliceAI: Jaganathan et al. 2019, Cell, doi:10.1016/j.cell.2018.12.015. Source
under PolyForm Strict 1.0.0, models under CC BY-NC 4.0, non-commercial use only;
run unmodified and not redistributed. Not for clinical decision-making.
