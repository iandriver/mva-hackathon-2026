# Reproducing the analysis

Every number in `findings/` and in the Track 1 report is produced by a script here.
Nothing is transcribed from a session.

## Inputs

Gated data never enters this repository. Scripts read it from `$MVA_DATA`:

```bash
export MVA_DATA=/path/to/mva-hackathon/data     # default: /Volumes/IanSSD/mva-hackathon/data
export MVA_REF=/path/to/mva-hackathon/ref       # ClinVar; downloaded on first use
export MVA_RESULTS=./results                    # default
```

The dataset is gated. Obtaining it requires an accepted Data Transfer Agreement at
<https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data>, then:

```bash
hf download SageBio/mva-hackathon-2026-data --repo-type dataset --local-dir "$MVA_DATA" \
  WGS_EX2312012_HGWCNDSX7.vcf.gz WGS_EX2312012_HGWCNDSX7.vcf.gz.tbi
```

A missing input names its fetch command rather than raising a traceback (`scripts/env.sh`).

## Order

| Script | Produces | Needs gated data | Runtime |
|---|---|---|---|
| `01_vcf_facts.sh` | reference build, sample, contig naming, FILTER census, record count | yes | ~2 min |
| `02_bub1b_locus.py` | every variant in the BUB1B span, VEP + gnomAD + predictors | yes | ~2 min |
| `03_secondary_findings.sh` | proband ∩ all ClinVar P/LP (downloads ClinVar, ~190 MB) | yes | ~10 min first run |
| `04_acmg_triage.py` | ACMG SF v3.2 reportability, with frequency/VAF/mapping rejection | yes | ~1 min |
| `05_mosaicism.py` | BAF + depth by chromosome, simulated detection floor | yes | ~10 min |
| `06_pocket_geometry.py` | ATP-pocket volume and variant distance (Arm B) | no | ~1 min |
| `07_lincs_arm.py` | LINCS reversal with positive and direction controls (Arm A) | no | ~5 min |
| `08_interferon_specificity.py` | ISG enrichment vs 300 random KOs + SAC panel | no | ~15 min first run |
| `09_ddg_stability.sh` | ThermoMPNN + DDGun ddG, calibrated (calls `09b`, `09c`) | no | ~30 min first run |
| `make_submission.py` | the Track 1 CSV, from `results/`, with format validation | yes | seconds |

Scripts 06–09 need no patient data at all, they operate on public structures and public
perturbation atlases. They can be run by anyone.

## What the controls are for

Three scripts exist mainly to refuse conclusions, and they are the point rather than overhead:

- **`07_lincs_arm.py`** runs a positive control (does the BUB1B-KO signature recover BUB1B
  knockout? rank 2 of ~70,000) *and* a direction control (does BUB1B overexpression reverse
  it?). The second fails, which is why no repurposing candidate is claimed from this arm.
- **`06_pocket_geometry.py`** is the check run *before* docking. The pocket is occluded
  (~63 Å³) and the variant is 33 Å away, so no ligands were docked.
- **`04_acmg_triage.py`** rejects on population frequency, allele balance and mapping bias
  before ClinVar's label is allowed to carry a finding. It is what rejects the PRSS1 hit.

## Splice and structural-variant work

Moved to [../splice_sv/](../splice_sv/) on 2026-09-10 as a self-contained
workstream: scripts 14 (proband SpliceAI, local), 15 (saturation-scan detection
floor) and 16 (SV capability + depth-based deletion test), the SpliceAI
calibration panel, the two findings documents, and the AlphaGenome positive
control. See `splice_sv/README.md`.

## Environments

`09_ddg_stability.sh` builds its own virtualenvs under `tools/` (gitignored) and pins the two
compatibility fixes DDGun 0.0.2 needs (`biopython==1.81`, `setuptools<81`). Everything else
needs only Python 3 stdlib plus `bcftools`.

External services used, all public and unauthenticated: Ensembl REST, Ensembl VEP REST,
gnomAD GraphQL, NCBI E-utilities and ClinVar FTP, UniProt REST, AlphaFold DB, PDBe,
SigCom-LINCS, Enrichr, ColabFold MMseqs2.

## Determinism

`05_mosaicism.py` and `08_interferon_specificity.py` sample; both fix a seed (3 and 11) and
record it in their JSON output. Everything else is deterministic given the same input files.
Results depend on live database releases, ClinVar and gnomAD in particular, so each output
JSON records what it used.
