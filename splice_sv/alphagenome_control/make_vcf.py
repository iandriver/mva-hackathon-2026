"""Write the panel as a VCF for SpliceAI, checking every REF base against the
reference FASTA first. A silent REF mismatch would produce plausible wrong
scores rather than an error."""
import json, pathlib, pysam

ROOT = pathlib.Path(__file__).parent
panel = json.load(open(ROOT/'data'/'panel.json'))
fa = pysam.FastaFile(str(ROOT/'ref'/'chr15.fa.bgz'))

bad = []
for r in panel:
    obs = fa.fetch('chr15', r['pos']-1, r['pos']).upper()
    if obs != r['ref'].upper():
        bad.append((r['name'], r['pos'], r['ref'], obs))
assert not bad, f"REF mismatch vs hg38 chr15: {bad[:5]}"
print(f"all {len(panel)} REF bases match hg38 chr15")

out = ROOT/'data'/'panel.vcf'
with open(out, 'w') as f:
    f.write('##fileformat=VCFv4.2\n')
    f.write('##contig=<ID=chr15,length=101991189>\n')
    f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
    for r in sorted(panel, key=lambda r: r['pos']):
        f.write(f"chr15\t{r['pos']}\t{r['name']}\t{r['ref']}\t{r['alt']}\t.\t.\t.\n")
print("wrote", out)
