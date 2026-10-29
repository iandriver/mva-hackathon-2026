# AlphaGenome splice detector, positive control

The AlphaGenome positive control and SpliceAI comparison. Now lives under
`splice_sv/alphagenome_control/` in the hackathon repo (moved 2026-09-10; it was
built as a standalone directory). The AlphaGenome positive control itself uses
only public gnomAD and literature variants and does not read, write, or reference
the proband's data; it is kept apart from the proband analysis because it exercises
an API tool still under organizer review. Two files in `results/` are the
exception: `bub1b_kmer_cn_profile.csv` and `bub1b_sv_screen.csv` are proband-derived
aggregate summaries (per-window copy-number ratios and a structural-variant screen
signal count; no genotypes, no allele fractions) produced by the `splice_sv/` SV
workstream and written alongside this directory's outputs.

## Why this exists

`findings/BUB1B_LOCUS.md` in the hackathon repo left twelve intronic BUB1B
variants unevaluated for splice impact, and `track2/MECHANISM.md` §2 records
that published MVA1 second alleles are repeatedly non-coding and
splice-disrupting. AlphaGenome is a candidate tool for that gap.

Before pointing any tool at the child's variants, two things had to be settled.
One is whether the AlphaGenome terms and the data-transfer agreement permit it,
which is open with the organizers. The other is whether the tool works on this
specific problem, which is what this directory answers, using only public data.

## Inputs

- Three published MVA1 splice alleles, coordinates resolved through the Ensembl
  variant recoder and independently re-derived from the Ensembl exon table.
- gnomAD v4 SNVs in BUB1B, retrieved from the public GraphQL API.

No proband data is used as input to the AlphaGenome panel.

## Outputs

- `results/per_variant_scores.csv` per-variant reduced scores
- `results/positive_ranks.csv` ranks of the three positives in the background
- `results/shards/` raw AlphaGenome score frames
- `FINDINGS.md` the result

## Run

    python3 -m venv .venv && ./.venv/bin/pip install alphagenome requests
    export ALPHAGENOME_API_KEY=...        # not stored here
    ./.venv/bin/python build_panel.py
    ./.venv/bin/python run_scores.py
    ./.venv/bin/python analyse.py
    ./.venv/bin/python operating_point.py

The key lives in `~/Downloads/GEO_database/.env`, which is gitignored there.
AlphaGenome credentials are personal under the terms and must not be published
or shared, including inside an organization.

## Terms

Outputs here are AlphaGenome Output and Derivatives, subject to
`LEGALLY_BINDING_TERMS_OF_USE.txt`. Non-commercial use only. Not for clinical
decision-making or medical advice. Citations: Avsec et al. 2026 (model),
Cheng et al. 2026 (Atlas scores, not used here).
