#!/usr/bin/env python3
"""Builds the DDGun profile for the BUBR1 kinase domain from a ColabFold MMseqs2 MSA."""
import json, os, re, subprocess, sys, tarfile, time, urllib.parse, urllib.request

REPO = os.path.join(os.path.dirname(__file__), "..")
WORK = os.environ.get("MVA_WORK", os.path.join(REPO, "tools", "ddgun-work"))
os.makedirs(WORK, exist_ok=True)
PROF = os.path.join(WORK, "bubr1_kinase.prof")
if os.path.exists(PROF):
    print(f"profile exists: {PROF}"); sys.exit(0)

seq = []
for l in open(os.path.join(REPO, "structures", "BUBR1_AF_kinase_766_1050.pdb")):
    if l.startswith("ATOM") and l[12:16].strip() == "CA":
        seq.append((int(l[22:26]), l[17:20].strip()))
m = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
     'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
     'TYR':'Y','VAL':'V'}
kin = "".join(m.get(a, "X") for _, a in sorted(seq))
print(f"kinase domain: {len(kin)} aa; residue 1002 = {kin[1002-766]}")

tid = json.load(urllib.request.urlopen(urllib.request.Request(
    "https://api.colabfold.com/ticket/msa",
    data=urllib.parse.urlencode({"q": f">q\n{kin}\n", "mode": "env"}).encode()), timeout=120))["id"]
for _ in range(90):
    s = json.load(urllib.request.urlopen(f"https://api.colabfold.com/ticket/{tid}", timeout=60))
    if s.get("status") in ("COMPLETE", "ERROR"): break
    time.sleep(10)
if s.get("status") != "COMPLETE":
    sys.exit(f"MSA server returned {s}")
tar = os.path.join(WORK, "msa.tar.gz")
urllib.request.urlretrieve(f"https://api.colabfold.com/result/download/{tid}", tar)
tarfile.open(tar).extractall(os.path.join(WORK, "msa_out"))

seqs, name, buf = [], None, []
for l in open(os.path.join(WORK, "msa_out", "uniref.a3m")):
    l = l.rstrip("\n")
    if l.startswith(">"):
        if name: seqs.append((name, "".join(buf)))
        name, buf = l[1:].split()[0], []
    elif l:
        buf.append(l)
if name: seqs.append((name, "".join(buf)))
seqs = [(n, re.sub("[a-z]", "", s)) for n, s in seqs]        # a3m lowercase = insertion
L = len(seqs[0][1]); seqs = [x for x in seqs if len(x[1]) == L]
ok = set("ACDEFGHIKLMNPQRSTVWY-")
aln = os.path.join(WORK, "bubr1_kinase.aln")
with open(aln, "w") as f:
    for n, s in seqs:
        f.write(f"{n} {''.join(c if c in ok else '-' for c in s.upper())}\n")
print(f"MSA: {len(seqs)} length-matched sequences")
subprocess.run([os.path.join(REPO, "tools", "ddg-env", "bin", "ddgun"),
                "mkprof", aln, "--output", PROF], check=True,
               stderr=subprocess.DEVNULL)
print(f"wrote {PROF}")
