#!/usr/bin/env bash
# Run ColabFold / AlphaFold2-Multimer v3 on one target (or all) on the GPU. Resume-safe.
#
#   scripts/run_colabfold.sh 9ASS
#   scripts/run_colabfold.sh all
#
# colabfold_batch flags: https://github.com/sokrypton/ColabFold (colabfold/batch.py)
# Output layout: results/raw/colabfold/<id>/
#    <id>_unrelaxed_rank_001_alphafold2_multimer_v3_model_k_seed_000.pdb   (pLDDT in B-factor)
#    <id>_scores_rank_001_alphafold2_multimer_v3_model_k_seed_000.json     (plddt, pae, max_pae, ptm, iptm[, ipsae...])
#    <id>.done.txt                                                         (completion marker)
# Ranking: 'multimer' metric = 0.8*ipTM + 0.2*pTM.
# GPU memory: 8 GB is enough for ~550 residues with alphafold2_multimer_v3 (jax unified memory is
#   enabled by colabfold_batch by default). If OOM appears at >500 aa, add --max-msa 256:512 or
#   set TF_FORCE_UNIFIED_MEMORY=1 XLA_PYTHON_CLIENT_MEM_FRACTION=4.0 (already exported below).
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

TARGET="${1:?usage: run_colabfold.sh <PDB_ID|all>}"
CONDA_ENV="${CONDA_ENV:-colabfold}"
ENV_BIN="${CONDA_PREFIX_OVERRIDE:-$HOME/miniforge3/envs/$CONDA_ENV}/bin"
CFB="$ENV_BIN/colabfold_batch"
# config is read with the scoring env's python (yaml available there); fall back to the colabfold env
PY="$HOME/miniforge3/envs/cofold/bin/python"; [[ -x "$PY" ]] || PY="$ENV_BIN/python"
cfg() { "$PY" -c "import yaml,sys; c=yaml.safe_load(open('config.yaml')); print(c$1)"; }

IN_DIR="$(cfg "['paths']['colabfold_inputs']")"
OUT_ROOT="$(cfg "['paths']['colabfold_out']")"
MODEL_TYPE="$(cfg "['colabfold']['model_type']")"
NUM_RECYCLE="$(cfg "['colabfold']['num_recycle']")"
NUM_MODELS="$(cfg "['colabfold']['num_models']")"
NUM_SEEDS="$(cfg "['colabfold']['num_seeds']")"
MSA_MODE="$(cfg "['colabfold']['msa_mode']")"
PAIR_MODE="$(cfg "['colabfold']['pair_mode']")"
RANK="$(cfg "['colabfold']['rank']")"
SEED="$(cfg "['colabfold']['random_seed']")"
LOG_DIR="$OUT_ROOT/logs"; mkdir -p "$LOG_DIR"

if [[ "$TARGET" == "all" ]]; then
  IDS=$(ls "$IN_DIR"/*.fasta | xargs -n1 basename | sed 's/\.fasta$//')
else
  IDS="$TARGET"
fi

[[ -x "$CFB" ]] || { echo "ERROR: $CFB not found -- install colabfold in env '$CONDA_ENV' (pip install 'colabfold[alphafold]' + jax cuda)"; exit 2; }
nvidia-smi --query-gpu=name,memory.used,memory.total --format=csv,noheader 2>/dev/null || echo "WARN: nvidia-smi unavailable"
export TF_FORCE_UNIFIED_MEMORY=1 XLA_PYTHON_CLIENT_MEM_FRACTION=4.0

for ID in $IDS; do
  FASTA="$IN_DIR/$ID.fasta"
  OUT="$OUT_ROOT/$ID"
  if [[ -f "$OUT/$ID.done.txt" ]]; then echo "[skip] $ID already predicted ($OUT/$ID.done.txt)"; continue; fi
  mkdir -p "$OUT"
  echo "[run ] $ID  $(date -Is)  $MODEL_TYPE recycle=$NUM_RECYCLE models=$NUM_MODELS seeds=$NUM_SEEDS"
  "$CFB" "$FASTA" "$OUT" \
      --model-type "$MODEL_TYPE" \
      --num-recycle "$NUM_RECYCLE" --num-models "$NUM_MODELS" --num-seeds "$NUM_SEEDS" \
      --msa-mode "$MSA_MODE" --pair-mode "$PAIR_MODE" \
      --rank "$RANK" --random-seed "$SEED" \
      ${COLABFOLD_EXTRA:-} 2>&1 | tee "$LOG_DIR/$ID.log"
  if [[ -f "$OUT/$ID.done.txt" ]]; then echo "[done] $ID"; else echo "[FAIL] $ID -- see $LOG_DIR/$ID.log"; fi
done
