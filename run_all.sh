#!/usr/bin/env bash
# cofold-bench single entry point.
#   ./run_all.sh                 # all stages, GPU stages are skipped unless the predictor is installed
#   ./run_all.sh select download # only these stages
# Stages: select | download | inputs | test | predict-boltz | predict-af2 | score | report
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

ENV_DIR="${COFOLD_ENV_DIR:-$HOME/miniforge3/envs/cofold}"
PY="$ENV_DIR/bin/python"
[[ -x "$PY" ]] || { echo "ERROR: cofold env python not found at $PY (create with: mamba env create -f environment.yml)"; exit 2; }
"$PY" -c "import cofoldbench" 2>/dev/null || { echo "[setup] pip install -e . into $ENV_DIR"; "$ENV_DIR/bin/pip" install -q -e .; }

STAGES=("$@"); [[ ${#STAGES[@]} -eq 0 ]] && STAGES=(select download inputs test predict-boltz predict-af2 score report)
has() { for s in "${STAGES[@]}"; do [[ "$s" == "$1" ]] && return 0; done; return 1; }
banner() { echo; echo "=================== [$1] $(date -Is) ==================="; }

if has select; then
  banner select
  if [[ -f data/targets/targets.csv && -z "${FORCE_SELECT:-}" ]]; then
    echo "targets.csv exists -> skipping RCSB re-selection (set FORCE_SELECT=1 to redo)"
  else
    "$PY" -m cofoldbench.select_targets
  fi
fi
if has download; then banner download; "$PY" -m cofoldbench.download; fi
if has inputs;   then banner inputs;   "$PY" -m cofoldbench.make_inputs; fi
if has test;     then banner test;     "$PY" -m pytest -q; fi

if has predict-boltz; then
  banner predict-boltz
  if [[ -x "$ENV_DIR/bin/boltz" ]] && nvidia-smi >/dev/null 2>&1 && [[ -n "${RUN_GPU:-}" ]]; then
    scripts/run_boltz.sh all
  else
    echo "SKIPPED: Boltz-2 GPU stage (phase 2). Requires boltz in $ENV_DIR, a GPU, and RUN_GPU=1."
  fi
fi
if has predict-af2; then
  banner predict-af2
  if [[ -x "$HOME/miniforge3/envs/colabfold/bin/colabfold_batch" ]] && nvidia-smi >/dev/null 2>&1 && [[ -n "${RUN_GPU:-}" ]]; then
    scripts/run_colabfold.sh all
  else
    echo "SKIPPED: AF2-Multimer GPU stage (phase 2). Requires colabfold_batch in env 'colabfold', a GPU, and RUN_GPU=1."
  fi
fi
if has score;  then banner score;  "$PY" -m cofoldbench.score; fi
if has report; then banner report; "$PY" -m cofoldbench.plot; fi
echo; echo "run_all.sh finished: ${STAGES[*]}"
