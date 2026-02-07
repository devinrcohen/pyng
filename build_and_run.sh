#!/usr/bin/env bash
set -euo pipefail

# ---------------- USER SETTINGS ----------------
PY_SCRIPT="${PY_SCRIPT:-smoketest_pyng.py}"      # default script to run
BUILD_DIR="${BUILD_DIR:-cmake-build-pyng-conda}" # CMake build directory

if [ ! -d "$BUILD_DIR" ]; then
  echo "[INFO] Build dir missing; configuring: $BUILD_DIR"
  cmake -S . -B "$BUILD_DIR" -G Ninja \
    -DPython3_EXECUTABLE="$CONDA_PREFIX/bin/python" \
    -DPython3_ROOT_DIR="$CONDA_PREFIX" \
    -DCMAKE_PREFIX_PATH="$CONDA_PREFIX"
fi

TARGET="${TARGET:-pyng}"                         # CMake target to build
# ----------------------------------------------

echo "=== build_and_run.sh ==="
echo "PWD:        $PWD"
echo "BUILD_DIR:  $BUILD_DIR"
echo "TARGET:     $TARGET"
echo "PY_SCRIPT:  $PY_SCRIPT"
echo "PYTHON:     $(command -v python || true)"
echo "CONDA_PREFIX: ${CONDA_PREFIX:-<unset>}"

if [[ ! -d "$BUILD_DIR" ]]; then
  echo "[FAIL] Build directory not found: $BUILD_DIR"
  echo "Tip: configure first, e.g. cmake -S . -B $BUILD_DIR -G Ninja ..."
  exit 1
fi

echo "[INFO] Building target '$TARGET'..."
cmake --build "$BUILD_DIR" -v --target "$TARGET"

# Sanity check: extension module should exist after build
if ! ls "$BUILD_DIR"/pyng.cpython-*.so >/dev/null 2>&1; then
  echo "[WARN] Could not find $BUILD_DIR/pyng.cpython-*.so"
  echo "Contents of build dir:"
  ls -la "$BUILD_DIR" || true
  echo "[WARN] Continuing anyway (your module name/build output might differ)."
fi

if [[ ! -f "$PY_SCRIPT" ]]; then
  echo "[FAIL] Python script not found: $PY_SCRIPT"
  exit 1
fi

# Ensure the build dir is on PYTHONPATH so "import pyng" resolves.
export PYTHONPATH="$PWD/$BUILD_DIR${PYTHONPATH:+:$PYTHONPATH}"

echo "[INFO] Running: python $PY_SCRIPT"
python "$PY_SCRIPT"
