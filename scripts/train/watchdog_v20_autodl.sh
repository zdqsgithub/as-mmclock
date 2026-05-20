#!/usr/bin/env bash
set -euo pipefail

ROOT="/root/mouse_methyl_work"
OUT="/root/autodl-tmp/mouse_methyl_work/results/autoresearch_v20_deep_learning"
LOG_DIR="/root/autodl-tmp/mouse_methyl_work/logs"
WATCHDOG_LOG="${LOG_DIR}/v20_watchdog.log"
WORKERS="${V20_WORKERS:-6}"

mkdir -p "${LOG_DIR}"

timestamp() {
  date -u +"%Y-%m-%dT%H:%M:%SZ"
}

is_done() {
  /root/mouse_methyl_work/.venv/bin/python - <<'PY'
from pathlib import Path
import sys
import pandas as pd

out = Path("/root/autodl-tmp/mouse_methyl_work/results/autoresearch_v20_deep_learning")
summary = out / "v20_dl_summary.csv"
lodo = out / "v20_lodo_summary.csv"
final_model = out / "final_all_data_model" / "final_deep_clock_model.pt"
if not (summary.exists() and lodo.exists() and final_model.exists()):
    sys.exit(1)
try:
    df = pd.read_csv(summary)
except Exception:
    sys.exit(1)
if "preprocess" not in df.columns or "status" not in df.columns:
    sys.exit(1)
robust_done = df[(df["preprocess"] == "robust") & (df["status"].fillna("completed") == "completed")]
if len(df) < 100 or robust_done.empty:
    sys.exit(1)
sys.exit(0)
PY
}

start_training() {
  local train_log
  train_log="${LOG_DIR}/v20_dl_watchdog_$(date +%Y%m%d_%H%M%S).log"
  cd "${ROOT}"
  screen -dmS v20_dl bash -lc \
    "/root/.local/bin/uv run python scripts/train/run_v20_dl_autoresearch_autodl.py --budget-hours 24 --top-lodo 12 --top-seed-repeats 3 --workers ${WORKERS} > ${train_log} 2>&1"
  echo "$(timestamp) restarted screen=v20_dl workers=${WORKERS} log=${train_log}" >> "${WATCHDOG_LOG}"
}

while true; do
  if is_done; then
    echo "$(timestamp) done summary_and_lodo_present" >> "${WATCHDOG_LOG}"
    exit 0
  fi

  if ! screen -ls | grep -q "[.]v20_dl"; then
    echo "$(timestamp) missing_screen restarting" >> "${WATCHDOG_LOG}"
    start_training
  else
    gpu_line="$(nvidia-smi --query-gpu=utilization.gpu,memory.used,power.draw,temperature.gpu --format=csv,noheader | head -1 || true)"
    rows="0"
    if test -s "${OUT}/v20_dl_summary_partial.csv"; then
      rows="$(wc -l < "${OUT}/v20_dl_summary_partial.csv")"
    fi
    procs="$(pgrep -fc train_deep_clock.py || true)"
    echo "$(timestamp) ok rows=${rows} procs=${procs} gpu=${gpu_line}" >> "${WATCHDOG_LOG}"
  fi

  sleep 300
done
