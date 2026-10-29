# AlphaGenome recovers known pathogenic MVA1 splice alleles in BUB1B

Run 2026-09-08. AlphaGenome 0.9.0, 100 kb context, GRCh38, ENST00000287598.11.
582 variants scored on three splice heads. Public data only.

## Question

Does AlphaGenome rank experimentally confirmed BUB1B splice alleles above their
positionally matched neighbours? If not, the tool cannot close the non-coding
gap and there is no reason to ask anyone for permission to run it on patient
data.

## Panel

| Set | n | Definition | Expectation |
|---|---|---|---|
| `positive` | 3 | Published MVA1 alleles with reported mis-splicing | high |
| `canonical` | 40 | gnomAD SNVs at intron offset ±1,2 (the GT/AG dinucleotides) | high, validity anchor |
| `background` | 463 | gnomAD intronic SNVs at \|offset\| 3–20, positives removed | the matched null |
| `common` | 16 | gnomAD intronic SNVs at AF ≥ 1e-3 | low |
| `deep` | 60 | gnomAD intronic SNVs at \|offset\| > 100 | low, floor |

The positives:

| Allele | GRCh38 | Offset | gnomAD AF | Reported effect |
|---|---|---|---|---|
| c.2386-11A>G | 15:40212488 A>G | acceptor −11 | 2.11e-05 | de novo acceptor preferred over the authentic site (Rio Frio 2010, NEJM, doi:10.1056/NEJMoa1006565) |
| c.1402-5A>G | 15:40200239 A>G | acceptor −5 | 2.06e-06 | mis-splicing with 4-bp insertion, BUBR1 severely diminished on immunoblot (Kawakami 2023, J Hum Genet, doi:10.1038/s10038-023-01211-8) |
| c.1288+5G>A | 15:40196779 G>A | donor +5 | 1.86e-05 | compound het with c.498_505del, MDS with monosomy 7 (Laberko 2019, Orphanet J Rare Dis, doi:10.1186/s13023-019-1073-x) |

The background is matched on position, not on frequency. No common variant
exists within 10 bp of a BUB1B splice site in gnomAD v4 (0 at AF ≥ 1e-4 for
\|offset\| 1–2, 3 for \|offset\| 3–10), so a frequency-matched near-splice
negative set cannot be assembled. The background is therefore rare variants of
the same positional class, which is also the comparison a user actually makes:
this candidate against its positional peers.

## Result

**Validity anchor.** Canonical dinucleotide variants separate from background by
roughly sixtyfold on median max\|raw\| score: SPLICE_SITES 0.9914 vs 0.0156,
SPLICE_SITE_USAGE 0.9515 vs 0.0176, SPLICE_JUNCTIONS 6.87 vs 0.21. AUC 0.977 to
0.991. The model reads splice sites at this locus.

**The test.** All three positives land in the top 4% of the matched background on
every metric.

| Allele | SPLICE_SITES | SPLICE_SITE_USAGE | SPLICE_JUNCTIONS |
|---|---|---|---|
| c.1402-5A>G | rank 1/464 (top 0.22%) | 5/464 (1.08%) | 7/464 (1.51%) |
| c.1288+5G>A | 8/464 (1.72%) | 17/464 (3.66%) | 12/464 (2.59%) |
| c.2386-11A>G | 5/464 (1.08%) | 8/464 (1.72%) | 14/464 (3.02%) |

Joint one-sided rank probability that all three land at or above the worst
observed rank by chance: 5.1e-06 on SPLICE_SITES, 2.8e-05 on SPLICE_JUNCTIONS,
4.9e-05 on SPLICE_SITE_USAGE. Positives-vs-background AUC 0.992, 0.978, 0.981.

The three non-canonical alleles score as high as genuine GT/AG disruptions
(positives 0.876 to 0.996; canonical median 0.991), which is what their
published phenotypes predict.

## Detection floor

At SPLICE_SITES ≥ 0.4250, the 95th percentile of the matched background:

| Set | Passing |
|---|---|
| positive | 3/3 (100%) |
| canonical | 39/40 (97.5%) |
| background | 24/463 (5.2%) |
| common | 0/16 (0%) |
| deep | 0/60 (0%) |

On a candidate set of twelve rare intronic variants this threshold gives 0.6
expected false positives at 3/3 sensitivity on the known alleles. Some of the
24 passing background variants are likely real splice disruptors rather than
errors, since the background is rare variants that were never functionally
tested, so 5.2% is an upper bound on the false positive rate.

## Does it beat the simple baseline?

The positives sit at offsets −5, +5 and −11. Two are inside the extended splice
motif, so a proximity heuristic might rank them without any model. Ranked over
the same 463-variant background, with proximity ties counted in the baseline's
favour:

| Allele | Offset | AlphaGenome | Proximity (−\|offset\|) |
|---|---|---|---|
| c.1402-5A>G | −5 | **1** / 464 | 74 / 464 |
| c.1288+5G>A | +5 | **8** / 464 | 74 / 464 |
| c.2386-11A>G | −11 | **5** / 464 | 229 / 464 |

Removing position entirely, by ranking each positive only against background
variants at the *same* offset:

| Allele | Rank among same-offset variants | Their median | This score |
|---|---|---|---|
| c.1402-5A>G | 1 / 12 | 0.0830 | 0.9959 |
| c.2386-11A>G | 1 / 7 | 0.0200 | 0.9843 |
| c.1288+5G>A | 4 / 18 | 0.2738 | 0.8755 |

Two of three are the top-scoring variant at their exact distance from the splice
site, by roughly an order of magnitude over the stratum median. The signal is
sequence, not position. Offset alone carries some information but not much
(Spearman rho = −0.284 against score across the background).

The Rio Frio allele is the case that matters. At −11 it is outside the canonical
motif, a proximity filter ranks it 229/464, and AlphaGenome ranks it 5/464 and
first among the seven background variants at that same offset. That is the class
of allele a position or motif filter misses and the reason to use a model here.

The Laberko allele is the weak one: 4/18 within its stratum, with three untested
background variants scoring above it. One of three positives is not cleanly
separated from its positional peers.

## SpliceAI comparison

SpliceAI 1.3.1, all five ensemble models, hg38 chr15 from UCSC, `-D 50`, bundled
GRCh38 annotation. Same 582 variants, same background, same statistic. All 582
REF bases were checked against the reference FASTA before scoring, which is a
third independent confirmation of the positives' coordinates.

Licensing: source under PolyForm Strict 1.0.0, models under CC BY-NC 4.0. Any
noncommercial purpose is permitted, and personal research use is named
explicitly. Two conditions: the software may not be redistributed or modified,
and the models require attribution. It is run here unmodified and not
redistributed. Citation: Jaganathan et al. 2019, Cell, doi:10.1016/j.cell.2018.12.015.

**Ranks in the matched background (n = 463):**

| Allele | Offset | AlphaGenome | SpliceAI | Proximity |
|---|---|---|---|---|
| c.1402-5A>G | −5 | 1 / 464 | 3 / 464 | 74 / 464 |
| c.2386-11A>G | −11 | 5 / 464 | 12 / 464 | 229 / 464 |
| c.1288+5G>A | +5 | 8 / 464 | 29 / 464 | 74 / 464 |

AUC against the background: positives 0.992 (AlphaGenome) vs 0.971 (SpliceAI);
canonical 0.991 vs 0.990. Both put `common` and `deep` variants below background
(AUC 0.33 and 0.44 for AlphaGenome, 0.25 and 0.36 for SpliceAI).

**The difference is not significant.** AlphaGenome ranks all three positives
higher, but with three positives the best attainable exact sign-test p is
0.5³ = 0.125. Three variants cannot distinguish these tools, and no claim that
one is better is supported here.

At SpliceAI's own recommended threshold, DS ≥ 0.5, it recovers 3/3 positives at
6.9% background and 0% of common and deep variants, which is the same regime as
AlphaGenome at its 5% operating point. The single miss reported earlier (2/3)
came from imposing a background-percentile threshold of 0.64 on SpliceAI rather
than its published cutoff, and is an artifact of that choice.

**They are complementary at the tail.** Spearman 0.762 across the panel. At each
tool's own 5% operating point: 54 variants flagged by both, 12 by AlphaGenome
only, 11 by SpliceAI only, Jaccard 0.70.

**Independent mechanistic agreement.** SpliceAI's directional calls reproduce
the published mechanisms without being told them:

| Allele | Acceptor gain | Acceptor loss | Donor gain | Donor loss | Published effect |
|---|---|---|---|---|---|
| c.1402-5A>G | 1.00 | 0.20 | 0.00 | 0.00 | mis-splicing with 4-bp insertion |
| c.2386-11A>G | 0.84 | 0.69 | 0.01 | 0.00 | de novo acceptor **preferred over** the authentic one |
| c.1288+5G>A | 0.00 | 0.00 | 0.56 | 0.31 | donor +5 |

The Rio Frio allele is called as a simultaneous acceptor gain and acceptor loss,
which is the signature of a competing site outcompeting the native one. That is
what the NEJM paper reports, and neither model was given that information.

## The comparison changes the plan

SpliceAI runs entirely locally. The weights ship in the package, no network call
is made, and no data leaves the machine. That removes the blocking question:
scoring the proband's twelve intronic variants with SpliceAI is not a
third-party disclosure, so neither the data-transfer agreement's prohibition on
granting data access nor AlphaGenome's HIPAA clause is engaged.

AlphaGenome stays blocked, because it is API-only and the weights are not
released. It buys a second opinion on a partially disjoint tail, which is worth
having but is not worth waiting for.

## Negative result: the quantile score is unusable here

AlphaGenome's `quantile_score`, taken as the maximum across tracks, is saturated
in this locus and separates nothing. Background median 0.9974, background
minimum 0.9704, and the 95th percentile, 99th percentile and maximum are all
1.0000. Positives and canonical variants are also at 1.0000. Ranking positives
by max quantile puts them at 11th, 11th and 14th, worse than every raw-score
metric, and that ranking is produced almost entirely by tie-breaking noise.

Use raw scores. A max-over-tracks quantile reduction saturates because some
track, somewhere in a 600-track panel, is always extreme.

## Bearing on the hackathon

The detector works on the exact variant class that MECHANISM.md §2 identifies as
the recurring MVA1 second allele, at a locus-specific operating point with a
measured floor. That makes the open question worth putting to the organizers:
whether the data-transfer agreement permits submitting the proband's variant
coordinates to a third-party API, and whether the release is HIPAA-regulated,
since the AlphaGenome terms prohibit processing HIPAA-covered health information.

Until both are answered, nothing here is run on proband data.

## Caveats

1. Three positives is a small set. It is every published BUB1B splice allele with
   an experimentally reported effect that could be resolved to GRCh38, not a
   sample from a larger pool.
2. The background is presumed benign by rarity and by the absence of reports, not
   by function. It is not a validated negative set.
3. Scores are locus-specific. The 0.4250 threshold is the BUB1B background's 95th
   percentile and does not transfer to another gene.
4. Splice heads only. Expression, chromatin and TF-binding heads were not run, so
   this says nothing about regulatory-element variants of the kind reported by
   Qu et al. 2025.
5. Prediction, not function. A high score is a hypothesis that needs RNA.
6. Proximity is a weak baseline. It is the one available without a model, but a
   splice-motif PWM (MaxEntScan) or SpliceAI would be a harder comparison and
   was not run. The claim here is that AlphaGenome beats position, not that it
   beats the best available splice predictor.

Output and Derivatives here are subject to `LEGALLY_BINDING_TERMS_OF_USE.txt`.
Model citation: Avsec et al. 2026.
