#!/usr/bin/env bash
# Run Boltz-2 on one target (or all) -- PHASE 2 (GPU). Resume-safe: skips finished targets.
#
#   scripts/run_boltz.sh 9ASS            # one target
#   scripts/run_boltz.sh all             # every data/inputs/boltz/*.yaml
#   BOLTZ_EXTRA="--use_potentials" scripts/run_boltz.sh 9ASS
#
# Input format / flags: https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md
# Output layout:  results/raw/boltz/<id>/boltz_results_<id>/predictions/<id>/
#                   <id>_model_{0..N-1}.cif, confidence_<id>_model_k.json (ptm, iptm, protein_iptm,
#                   complex_plddt, pair_chains_iptm ...).  Models are sorted by confidence_score.
# GPU memory: an RTX 4060 Laptop (8 GB) handles ~550 tokens with recycling_steps 3 and
#   diffusion_samples 5 (samples run sequentially: --max_parallel_samples 1 keeps peak memory low).
#   If CUDA OOM: lower --max_parallel_samples, --diffusion_samples, or add --no_kernels.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

TARGET="${1:?usage: run_boltz.sh <PDB_ID|all>}"
CONDA_ENV="${CONDA_ENV:-cofold}"
PY="${CONDA_PREFIX_OVERRIDE:-$HOME/miniforge3/envs/$CONDA_ENV}/bin/python"
BOLTZ="$(dirname "$PY")/boltz"
cfg() { "$PY" -c "import yaml,sys; c=yaml.safe_load(open('config.yaml')); print(c$1)"; }

IN_DIR="$(cfg "['paths']['boltz_inputs']")"
OUT_ROOT="$(cfg "['paths']['boltz_out']")"
RECYCLE="$(cfg "['boltz']['recycling_steps']")"
SAMPLES="$(cfg "['boltz']['diffusion_samples']")"
STEPS="$(cfg "['boltz']['sampling_steps']")"
FORMAT="$(cfg "['boltz']['output_format']")"
PAIRING="$(cfg "['boltz']['msa_pairing_strategy']")"
SEED="$(cfg "['boltz']['seed']")"
USE_POT="$(cfg "['boltz']['use_potentials']")"
LOG_DIR="$OUT_ROOT/logs"; mkdir -p "$LOG_DIR"

if [[ "$TARGET" == "all" ]]; then
  IDS=$(ls "$IN_DIR"/*.yaml | xargs -n1 basename | sed 's/\.yaml$//')
else
  IDS="$TARGET"
fi

[[ -x "$BOLTZ" ]] || { echo "ERROR: $BOLTZ not found -- install boltz in env '$CONDA_ENV' (pip install 'boltz[cuda]')"; exit 2; }
nvidia-smi --query-gpu=name,memory.used,memory.total --format=csv,noheader 2>/dev/null || echo "WARN: nvidia-smi unavailable"

EXTRA="${BOLTZ_EXTRA:-}"
[[ "$USE_POT" == "True" ]] && EXTRA="$EXTRA --use_potentials"

for ID in $IDS; do
  YAML="$IN_DIR/$ID.yaml"
  OUT="$OUT_ROOT/$ID"
  DONE="$OUT/boltz_results_$ID/predictions/$ID/confidence_${ID}_model_0.json"
  if [[ -f "$DONE" ]]; then echo "[skip] $ID already predicted ($DONE)"; continue; fi
  mkdir -p "$OUT"
  echo "[run ] $ID  $(date -Is)  recycling=$RECYCLE samples=$SAMPLES steps=$STEPS"
  # --no_kernels: the cuEquivariance CUDA kernels are optional and not installed in the cofold env;
  # the plain PyTorch path gives identical results, ~1.5-2x slower.
  # shellcheck disable=SC2086
  "$BOLTZ" predict "$YAML" \
      --out_dir "$OUT" \
      --use_msa_server --msa_pairing_strategy "$PAIRING" --no_kernels \
      --recycling_steps "$RECYCLE" --diffusion_samples "$SAMPLES" --sampling_steps "$STEPS" \
      --max_parallel_samples 1 \
      --output_format "$FORMAT" --seed "$SEED" --write_full_pae \
      $EXTRA 2>&1 | tee "$LOG_DIR/$ID.log"
  if [[ -f "$DONE" ]]; then echo "[done] $ID"; else echo "[FAIL] $ID -- see $LOG_DIR/$ID.log"; fi
done
