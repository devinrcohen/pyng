#!/usr/bin/env bash
set -euo pipefail

# ---------------- USER SETTINGS ----------------
PY_SCRIPT="${PY_SCRIPT:-smoketest_pyng.py}"      # default script to run
BUILD_DIR="${BUILD_DIR:-cmake-build-pyng-conda}" # where the .so lives
# -----------------------------------------------

echo "=== run_python.sh ==="
echo "PWD:        $PWD"
echo "BUILD_DIR:  $BUILD_DIR"
echo "PY_SCRIPT:  $PY_SCRIPT"
echo "PYTHON:     $(command -v python || true)"
echo "CONDA_PREFIX: ${CONDA_PREFIX:-<unset>}"

if [[ ! -d "$BUILD_DIR" ]]; then
  echo "[FAIL] Build directory not found: $BUILD_DIR"
  exit 1
fi

if [[ ! -f "$PY_SCRIPT" ]]; then
  echo "[FAIL] Python script not found: $PY_SCRIPT"
  exit 1
fi

# Ensure the build dir is on PYTHONPATH so "import pyng" resolves.
export PYTHONPATH="$PWD/$BUILD_DIR${PYTHONPATH:+:$PYTHONPATH}"

echo "[INFO] Running: python $PY_SCRIPT"
python "$PY_SCRIPT"
