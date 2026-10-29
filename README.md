# MVA Hackathon 2026: BUB1B compound heterozygote, and why no repurposed drug fits

Entry by @idriver7 to *Rare Disease, Real Kid: The MVA Hackathon 2026*, organised by Sage
Bionetworks with the MVA Society, Hugging Face and BEACON.

## In brief

**Track 1, the cause.** The child is a compound heterozygote for **BUB1B**, the gene behind
mosaic variegated aneuploidy type 1. One allele carries p.Leu737Ter, a ClinVar
pathogenic/likely-pathogenic nonsense variant that truncates the protein before its kinase
domain and is expected to be null. The other carries p.Asn1002Lys, a missense change in the
pseudokinase domain found on 1 of 1,461,878 gnomAD v4 exome alleles. Two stability predictors
put it at about 1.5 kcal/mol destabilising. Phase is inferred, not demonstrated. Mosaic
whole-chromosome gain in 12% or more of cells is excluded in the released blood genome, so the
diagnosis rests on the genotype. A splice, copy-number, structural-variant and mobile-element
screen found no hidden third allele.

**Track 2, the drug question.** No repurposed approved drug is proposed, and the report argues
why. Five candidate classes were each refused against a control fixed in advance:
aneuploidy-selective compounds, LINCS signature reversal, pocket docking, nonsense readthrough,
and pharmacological chaperones. The refusals share one cause. The only level at which treating
the disease and lowering cancer risk point the same way is upstream, at the rate of chromosome
missegregation, and that level needs functional BUBR1, which has no druggable site. Section 6 of
the report states what would falsify this.

## Where to start

| Path | Contents |
|---|---|
| `submissions/idriver7_track1_report.md` | Track 1 report |
| `submissions/idriver7_track2_report.md` | Track 2 report |
| `submissions/idriver7_bub1b-compound-het.csv` | Track 1 prediction (primary row) |
| `findings/` | one document per analysis: locus, ΔΔG stability, mosaicism, both repurposing arms, interferon specificity, readthrough, chaperones, modality map, the closed tolerance window |
| `track2/MECHANISM.md`, `track2/METHOD.md` | mechanism groundwork and the signature-reversal design, both fixed before data access |
| `splice_sv/` | splice, k-mer copy-number, alignment, SV and mobile-element work closing the "hidden third allele" question; start at `splice_sv/README.md` |
| `structures/` | AlphaFold models and ThermoMPNN predictions for BUBR1 |
| `scripts/` | analysis scripts 01–13 plus `make_submission.py`; run order, runtimes and inputs in `scripts/README.md` |

## Reproducing

Scripts 01–05 and `make_submission.py` need the gated dataset, which requires an accepted Data
Transfer Agreement at <https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data>. Point
`MVA_DATA` at it (see `scripts/env.sh`). Scripts 06–13 use only public structures and
perturbation atlases, so anyone can run them. Analysis environments (ThermoMPNN, DDGun, docking)
were built locally under `tools/` and are not committed.

## Data handling

No primary Hackathon data is in this repository. `DATA_HANDLING.md` lists what patient-derived
content is here, what was removed before publication, and how the public copies differ from
the submitted versions. Please do not attempt to contact the child or the family.

## AI assistance

Analysis orchestration, code and drafting used Anthropic's Claude through Claude Code. The
provider, tier, data-handling settings and applicable terms are declared in Track 2 report §9.
Every number in the reports is produced by the scripts here, not by the model.

## Licence

Reports and findings are CC BY 4.0, code is MIT, and third-party model outputs (AlphaGenome,
SpliceAI, AlphaFold DB) keep their own terms. See `LICENSE`.

## Acknowledgements

With thanks to the child and family who opened this case to the research community, and to the
organisers: Sage Bionetworks, the MVA Society, Hugging Face and BEACON.
