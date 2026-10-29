# Splice and structural-variant evaluation of BUB1B: work index

Session 2026-09-08 to 09, consolidated under `splice_sv/` on 2026-09-10.
Entry point for a block of work that ran three tools against
the BUB1B locus to close two caveats left open in
[BUB1B_LOCUS.md](../findings/BUB1B_LOCUS.md): whether the twelve intronic variants disrupt
splicing (caveat 4), and whether a structural allele was missed. Each piece has
its own findings document; this ties them together and records what lives where.

## Motivation

[MECHANISM.md](../track2/MECHANISM.md) §2 records that published MVA1 second
alleles are repeatedly non-coding and splice-disrupting, and that the Kawakami
case's second allele was a 3.0 kb Alu/SVA insertion invisible to short-read
sequencing. The compound-heterozygous call rests on two coding variants with phase
undemonstrated, so a missed third allele of either class is the standing risk.

## What was done, in order

1. **AlphaGenome terms review.** A key was added to `~/Downloads/GEO_database/.env`
   (gitignored there, never in this repo). The terms permit individual
   non-commercial use but require attribution notices on published output and
   prohibit processing HIPAA-covered health information; it is API-only, so
   submitting proband variants is a third-party disclosure. Left blocked on the
   organizers, pending the data-transfer and HIPAA questions.

2. **AlphaGenome positive control**, public data only, in `alphagenome_control/`. Established that
   AlphaGenome recovers 3/3 published pathogenic MVA1 splice alleles in the top
   4% of a position-matched background, that its `quantile_score` is saturated at
   this locus and unusable, and that it beats a proximity baseline.

3. **SpliceAI comparison**, also in `alphagenome_control/`. SpliceAI matches
   AlphaGenome on the three positives (the difference is not significant at n=3)
   and, decisively, **runs entirely locally** — weights ship in the package, no
   network call. That removed the blocker for the proband scan.

4. **Proband intronic scan** → [SPLICE_INTRONIC.md](findings/SPLICE_INTRONIC.md),
   `scripts/14_splice_intronic.py` and `scripts/15_splice_detection_floor.py`.
   All fourteen gene-body variants score below DS 0.03. A saturation scan of
   3,936 neighbouring SNVs shows eleven of fourteen positions are BLIND (no
   possible change reaches DS 0.2), so those zeros are uninformative. The two
   *coding* variants are genuine negatives — the model reads splice signal there
   — which confirms the pseudokinase stability argument's assumption that
   p.Asn1002Lys acts through the protein, not the transcript.

5. **Structural-variant assessment** → [BUB1B_STRUCTURAL.md](findings/BUB1B_STRUCTURAL.md),
   `scripts/16_bub1b_sv_assessment.py`. The VCF is a GATK small-variant callset
   with no SV calling; largest event genome-wide is 1,170 bp, below the 3.0 kb
   Kawakami class. A large heterozygous deletion is excluded at the 13 informative
   positions but 21 of 23 exons sit inside windows with no depth information. An
   exon-level deletion is untested, not absent.

6. **Copy-number profile** → [findings/BUB1B_KMER_CN.md](findings/BUB1B_KMER_CN.md),
   `scripts/17_kmer_panel.py`, `scripts/18_kmer_count.py`, `scripts/19_kmer_cnv_profile.py`.
   All 8 FASTQs (84.67 GB) streamed from HuggingFace and matched against 58,266
   unique BUB1B 31-mers, reads never stored. 21 of 23 exons are diploid with no
   window at hemizygous depth; exons 4 and 16 are repeat-dense and uninterpretable.
   Excludes a heterozygous exon-level deletion where the VCF had no depth.

7. **Alignment + mobile-element screen** → [findings/BUB1B_SV_ALIGNMENT.md](findings/BUB1B_SV_ALIGNMENT.md),
   `scripts/20_targeted_align.sh`, `scripts/21_delly_call.sh`, `scripts/22_me_screen.py`.
   All four lanes aligned to hg38, BUB1B-anchored pairs kept. Delly called no SVs;
   coverage is uniform (covering exons 4 and 16 that the k-mer method could not);
   every soft-clip and Alu/L1/SVA-homology signal coincides with a native
   reference repeat, and the clip clusters are one-sided with no polyA. No
   insertion or balanced event.

## Net result

Both caveats are narrowed, neither is closed, and both negatives carry detection
floors rather than being asserted. Two things were actually settled: neither
coding variant disrupts splicing, and no deletion overlaps the covered positions.
Everything else — the twelve deep-intronic variants, any insertion, any balanced
event, the 21 uncovered exons — is bounded as untested with the power number
attached. Nothing changes the Track 1 result (F-max 1.000).

## Layout (everything under `splice_sv/`, git-tracked in place)

**Retainable — code and derived analysis, no genotypes:**

| Path | Content |
|---|---|
| `README.md` (this file) | index |
| `findings/SPLICE_INTRONIC.md` | proband splice result + detection floor |
| `findings/BUB1B_STRUCTURAL.md` | SV assessment + exon coverage map |
| `scripts/14_splice_intronic.py` | scores proband BUB1B variants with SpliceAI (local) |
| `scripts/15_splice_detection_floor.py` | saturation scan establishing the floor |
| `scripts/16_bub1b_sv_assessment.py` | SV capability + depth-based deletion test |
| `scripts/calibration/bub1b_spliceai_calibration.csv` | 582 public gnomAD/literature variants; the calibration panel |
| `alphagenome_control/` | the AlphaGenome positive control and SpliceAI comparison, public data only |

`alphagenome_control/` holds `FINDINGS.md` (the writeup), the pipeline
(`build_panel.py`, `run_scores.py`, `analyse.py`, `operating_point.py`,
`baseline.py`, `compare.py`, `make_vcf.py`), `data/` (the public gnomAD and
literature panel), `results/` (reduced per-variant scores and the SpliceAI VCF),
and `LEGALLY_BINDING_TERMS_OF_USE.txt`, the AlphaGenome attribution notice
required if that output is published. It touches no proband data.

**Gitignored, carry proband genotypes (repo-root `results/`, on the deletion
list):** `bub1b_proband_snvs.vcf`, `bub1b_proband_spliceai_D50.vcf`,
`bub1b_proband_splice_scores.csv`, `bub1b_floor_scan*.vcf`, `floor.log`.
Regenerate with scripts 14 and 15.

**Gitignored, regenerable or bulky:** `alphagenome_control/.venv*` (Python
environments), `alphagenome_control/results/shards/` (raw AlphaGenome score
frames; the reduced CSVs are tracked), `alphagenome_control/results/*.log`, and
the reference `/Volumes/IanSSD/mva-hackathon/ref/chr15.fa.bgz` (hg38 chr15 from
UCSC).

## Tool provenance

- SpliceAI 1.3.1 — Jaganathan et al. 2019, Cell, doi:10.1016/j.cell.2018.12.015.
  Source PolyForm Strict 1.0.0, models CC BY-NC 4.0. Run unmodified, not
  redistributed, non-commercial.
- AlphaGenome (API) — Avsec et al. 2026. Non-commercial terms; output attribution
  required. Not run on proband data.
- Neither tool is for clinical decision-making.

## Still open

- k-mer depth profile over the FASTQs to fill the 21 uncovered exons (alignment-
  free, local, not DTA-blocked). The proportionate next step for the SV question.
- Full alignment + Manta/Delly + MELT for the mobile-element class specifically.
- AlphaGenome on the proband variants, if the organizers clear the API disclosure
  and HIPAA question. Its 100 kb context may behave differently in deep intron
  where SpliceAI is flat.
