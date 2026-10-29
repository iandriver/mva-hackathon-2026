# The released WGS does not show mosaic aneuploidy, with a floor

Run 2026-09-04 on `WGS_EX2312012_HGWCNDSX7.vcf.gz`. Input: 2,281,104 heterozygous PASS SNPs on
chr1-22/X, depth filtered to 20-70.

## Why the naive statistic fails

Per-chromosome maximum of window mean |BAF - 0.5| ranges 0.082-0.225 across autosomes, which
looks like strong mosaicism. It is not. The largest, chr7 at 58,067,345, is the chr7
centromere; the other extremes fall on chr1, chr9, chr16, chr20, chr21, chr22, pericentromeric
and acrocentric repeat regions where multi-mapping inflates allele-fraction deviation. The
maximum-window statistic measures mapping artifact, not biology.

The robust statistic is the **median across windows** of 2,000 het SNPs, which those artifacts
cannot move.

## Result

| | |
|---|---|
| median-window statistic, 22 autosomes | 0.06495 (chr18) – 0.06881 (chr19) |
| between-chromosome spread | **0.00386** |
| sd across chromosomes | 0.00080 |
| whole-chromosome depth ratio | 1.000 – 1.045 |

## Detection floor by simulation

A mosaic trisomy at cell fraction *f* gives BAF modes at (1+f)/(2+f) and 1/(2+f) and a depth
ratio of 1+f/2. Simulating that against the observed depth distribution:

| f | simulated statistic | shift vs null | exceeds observed spread |
|---|---|---|---|
| 0.02 | 0.06125 | 0.00037 | no |
| 0.06 | 0.06187 | 0.00099 | no |
| 0.10 | 0.06437 | 0.00349 | no |
| **0.12** | **0.06529** | **0.00440** | **yes** |
| 0.15 | 0.06689 | 0.00601 | yes |
| 0.20 | 0.07122 | 0.01034 | yes |

**A mosaic whole-chromosome gain in ≥12% of cells would shift one chromosome by more than the
entire observed between-chromosome spread. None does.** The depth channel agrees and is
slightly tighter (max ratio 1.045, f ≤ 0.09), though chr16/19/22 carry GC coverage bias.

MVA is *variegated*, different chromosomes in different cells, so the between-chromosome
spread is the right empirical null: a real mosaic gain would make one chromosome stand out
from the others, and nothing does.

## Reading

12% is not an arbitrary threshold. Pavone et al. reported **12% mosaic metaphases** on
karyotype in an MVA1 patient (doi:10.1007/s10072-022-06247-w). We can exclude mosaicism at
that magnitude in this sample and not below it.

The likely explanation is tissue rather than absence of disease. Published MVA mosaicism is
assayed in cultured lymphocytes or fibroblasts; aneuploid cells are selected against in the
proliferating lymphocyte pool that dominates bulk blood DNA, and the released WGS is bulk.

**Consequence: the diagnosis rests on the genotype, and the defining phenotype is not
confirmable from the released data.**

## Not done

- Segmental and arm-level events. The windowed scan would show them, but the artifact
  contamination at centromeres was not modelled well enough to call sub-chromosomal changes,
  and no such call is made.
- Mosaic *loss* was not separately simulated; the floor above is for gain.
- Uniparental disomy, which produces BAF deviation without depth change, was not tested.
