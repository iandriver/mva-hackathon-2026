# Data handling

## What is not in this repository

No primary hackathon data. The VCF, its index, the FASTQs and the clinical phenotype document
live outside the repository on local encrypted storage and are excluded from git. They will be
deleted at the conclusion of the Hackathon and deletion confirmed by email to
MVAHackathon2026@synapse.org, per the Data Transfer Agreement (WCG IRB protocol #20252010).

## What is in this repository

Derived results about a real, living child, limited to the primary MVA question:

| File | Patient-derived content |
|---|---|
| `submissions/*.csv` | the two causal BUB1B variants |
| `findings/BUB1B_LOCUS.md` | genotypes and call quality at the BUB1B locus |
| `findings/MOSAICISM.md` | genome-wide allele-fraction summaries (aggregate, no genotypes) |
| `splice_sv/alphagenome_control/results/bub1b_kmer_cn_profile.csv`, `bub1b_sv_screen.csv` | aggregate copy-number and structural-variant screen summaries over the BUB1B locus (no genotypes, no allele fractions) |

Every genotype-bearing output beyond the two causal variants stays gitignored and is never
committed. The one exception is in the reports: they quote the allele fraction and mapping
metrics of a single **rejected** ClinVar call (PRSS1, inside the somatically recombining TRB
locus and present in about 20% of gnomAD exomes). That call is kept to document the filter that
rejected it, not as a finding.

## Removed before this repository was made public

The Hackathon requires the repository to be public once it ends, with submissions released
under CC BY 4.0. That requirement is about the submission and the method. Two categories of
additional private information about the child were removed from the working tree **and from
the git history** before the visibility change, so they are not recoverable from this
repository:

1. **An incidental secondary finding** — a tumour-predisposition allele unrelated to the
   primary MVA diagnosis, together with the carrier variants surfaced by the same unbiased
   ClinVar scan. Publishing the causal variants is the point of the exercise; publishing an
   unrelated predisposition allele is a separate disclosure about a child, so it was not
   published.
2. **The content of the challenge phenotype document**, which is marked
   *Confidential · Do not redistribute*. The working transcription was deleted. Where the
   submitted reports quoted that document directly (HPO identifiers, birth details, and family
   reproductive history), the public copies replace the quotation with what the MVA diagnosis
   itself implies: a childhood rhabdomyosarcoma with prenatal growth restriction.

The phenotype document also states that publications must not include information that could
re-identify the patient or family beyond what is already public through the family's own blog
posts. The public copies refer to the child's diagnosis only at that level.

## Public copies versus the submitted versions

The two reports and the Track 1 CSV here differ from the versions submitted to the organisers
only by the removals above: the incidental secondary-finding row and its discussion, and the
phenotype quotations. The primary BUB1B row, and every analysis and number bearing on it, are unchanged.

## Licence

See `LICENSE`. Reports and findings are CC BY 4.0, code is MIT, and the third-party outputs
listed there keep their own terms.
