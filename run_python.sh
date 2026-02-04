#!/usr/bin/env bash
set -euo pipefail

PY_SCRIPT="${1:-smoketest_pyng.py}"
BUILD_DIR="${BUILD_DIR:-cmake-build-pyng-conda}"

# Prefer active conda python if present, otherwise fall back to python3/python.
if [[ -n "${CONDA_PREFIX:-}" && -x "${CONDA_PREFIX}/bin/python" ]]; then
  PYTHON="${CONDA_PREFIX}/bin/python"
else
  PYTHON="$(command -v python3 || command -v python)"
fi

echo "=== run_python.sh ==="
echo "PWD:         $(pwd)"
echo "BUILD_DIR:   ${BUILD_DIR}"
echo "PY_SCRIPT:   ${PY_SCRIPT}"
echo "PYTHON:      ${PYTHON}"
echo "CONDA_PREFIX:${CONDA_PREFIX:-<unset>}"

# Ensure the extension module is discoverable.
# For a pybind11_add_module(pyng ...) built with Ninja, the .so ends up directly in BUILD_DIR.
export PYTHONPATH="$(pwd)/${BUILD_DIR}:$(pwd):${PYTHONPATH:-}"

echo "[INFO] Running: ${PYTHON} ${PY_SCRIPT}"
exec "${PYTHON}" "${PY_SCRIPT}"
