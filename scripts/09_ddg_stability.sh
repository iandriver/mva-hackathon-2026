#!/usr/bin/env bash
# Reproduces findings/N1002K_STABILITY.md -- ddG of folding for BUBR1 p.Asn1002Lys.
#
# Two predictors with OPPOSITE sign conventions and unrelated inductive biases:
#   ThermoMPNN  supervised GNN, transfer-learned from the Tsuboyama megascale set
#               (positive = destabilising)
#   DDGun       untrained; substitution matrix + hydrophobicity + conservation
#               (negative = destabilising)
# A single number means nothing, so both are calibrated against site-saturation
# across the whole kinase domain and reported as a percentile.
source "$(dirname "$0")/env.sh"
TOOLS="$MVA_REPO/tools"
STRUCT="$MVA_REPO/structures"
mkdir -p "$TOOLS" "$STRUCT"
need_cmd uv "brew install uv"

# --- structure: AlphaFold, since no PDB entry resolves the pseudokinase domain ---
if [ ! -f "$STRUCT/BUBR1_AF_full.pdb" ]; then
  URL=$(curl -s https://alphafold.ebi.ac.uk/api/prediction/O60566 | python3 -c 'import json,sys;print(json.load(sys.stdin)[0]["pdbUrl"])')
  curl -sL -o "$STRUCT/BUBR1_AF_full.pdb" "$URL"
fi
python3 - "$STRUCT" <<'PY'
import sys, os
s = sys.argv[1]
out = os.path.join(s, "BUBR1_AF_kinase_766_1050.pdb")
if not os.path.exists(out):
    with open(out, "w") as o:
        for l in open(os.path.join(s, "BUBR1_AF_full.pdb")):
            if l.startswith(("ATOM", "HETATM")) and 766 <= int(l[22:26]) <= 1050:
                o.write(l)
        o.write("TER\nEND\n")
    print(f"wrote {out}")
PY

# --- ThermoMPNN ---
if [ ! -d "$TOOLS/ThermoMPNN" ]; then
  git clone -q --depth 1 https://github.com/Kuhlman-Lab/ThermoMPNN.git "$TOOLS/ThermoMPNN"
  # local.yaml ships with the authors' cluster paths and a gpu default
  python3 - "$TOOLS/ThermoMPNN" <<'PY'
import sys, pathlib
p = pathlib.Path(sys.argv[1]) / "local.yaml"
t = p.read_text().replace('accel: "gpu"', 'accel: "cpu"')
t = t.replace('thermompnn_dir: "/proj/kuhl_lab/ThermoMPNN"',
              f'thermompnn_dir: "{sys.argv[1]}"')
import re
t = re.sub(r'cache_dir: ".*"', f'cache_dir: "{sys.argv[1]}/.cache"', t)
p.write_text(t)
PY
fi
if [ ! -d "$TOOLS/tmpnn-env" ]; then
  uv venv --python 3.10 -q "$TOOLS/tmpnn-env"
  VIRTUAL_ENV="$TOOLS/tmpnn-env" uv pip install -q torch torchmetrics pytorch-lightning \
    biopython tqdm pandas numpy omegaconf joblib wandb
fi
mkdir -p "$TOOLS/ThermoMPNN/.cache"
for pdb in BUBR1_AF_kinase_766_1050 BUBR1_AF_full; do
  [ -f "$STRUCT/ThermoMPNN_inference_${pdb}.csv" ] && continue
  ( cd "$TOOLS/ThermoMPNN" && "$TOOLS/tmpnn-env/bin/python" analysis/custom_inference.py \
      --pdb "$STRUCT/${pdb}.pdb" --chain A \
      --model_path models/thermoMPNN_default.pt --out_dir "$STRUCT/" ) >/dev/null 2>&1
done

# --- DDGun (needs an MSA; ColabFold's public MMseqs2 server supplies one) ---
if [ ! -d "$TOOLS/ddg-env" ]; then
  uv venv --python 3.11 -q "$TOOLS/ddg-env"
  # ddgun 0.0.2 predates the Biopython three_to_one removal and needs pkg_resources
  VIRTUAL_ENV="$TOOLS/ddg-env" uv pip install -q ddgun "biopython==1.81" "setuptools<81"
fi
"$MVA_REPO/scripts/09b_ddgun_msa.py"

# --- report ---
"$MVA_REPO/scripts/09c_ddg_report.py"
