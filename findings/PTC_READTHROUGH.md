# Nonsense readthrough for p.Leu737Ter: assessed, and weak

Run 2026-09-07. Reproduce with `scripts/11_ptc_readthrough_context.py`.

The null allele is a premature termination codon, so readthrough agents are an obvious
candidate class and one of the few Goal A routes with approved or near-approved drugs
(ataluren, aminoglycosides such as gentamicin and amikacin, ELX-02). Readthrough is strongly
sequence-context dependent, so this is answerable from the sequence rather than by assertion.

## The PTC and its context

`ENST00000287598` codon 737 is **TTA**, and c.2210T>G makes it **TGA** (opal).

```
   -9  -8  -7  -6  -5  -4  -3  -2  -1 |  +1  +2  +3  +4  +5  +6 ...
    C   T   A   C   C   A   G   A   G |   T   G   A   A   G   T
```

Scored against the positions Tarraso et al. 2020 identify as influential
(doi:10.1242/dmm.043281; the first base of the stop is +1):

| position | this PTC | favourable base | |
|---|---|---|---|
| stop codon | **TGA** | TGA | **favourable** — the most readthrough-permissive stop |
| −1 | **G** | G | **favourable** |
| −2 | **A** | not C | **favourable by absence** — C at −2 counteracts readthrough |
| −3 | G | C | unfavourable |
| −9 | C | G | unfavourable |
| +4 | **A** | C | **unfavourable** |

So: the best stop codon, one favourable flank, one counteracting base avoided, and three
positions that are not favourable, including +4.

## Two reasons this is not a strong candidate

**The closest tested comparator failed.** Tarraso et al. tested PYGM p.R50X, also a TGA, with
a context they judged partially favourable (G at −1 and −9, C at −3, counteracted by C at −2
and no C at +4). Across three cellular models and four agents — amlexanox, ataluren, RTC13 and
G418 — they found **no detectable readthrough**. Our context is comparable in quality: better
at −2, worse at −3 and −9, and identical in lacking C at +4.

**The transcript is NMD-competent.** The PTC sits in exon 17 of 23, 314 codons before the
native stop, well upstream of the last exon-exon junction. Nonsense-mediated decay degrades
the message, so there is less template for a readthrough agent to act on. This obstacle is
independent of sequence context and is the larger of the two.

Combining readthrough with NMD inhibition is the standard answer to the second problem, and it
is not available here: NMD inhibitors are not approved and carry substantial toxicity. That
takes the combination outside a repurposing brief.

## What this changes

It moves the target. **In this compound heterozygote the druggable allele is the hypomorph,
not the null.** p.Asn1002Lys produces full-length protein destabilised by roughly 1.5 kcal/mol
(N1002K_STABILITY.md); p.Leu737Ter produces a transcript that is largely degraded. Raising the
abundance or stability of the protein that exists is a better-posed problem than manufacturing
protein from a message that is not there.

Recorded as assessed-and-weak rather than untried. Aminoglycoside readthrough would in any
case be a hard sell in a child: ototoxicity and nephrotoxicity against a candidate with no
supporting context data and a failed comparator.
