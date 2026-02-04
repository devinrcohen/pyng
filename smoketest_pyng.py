#!/usr/bin/env python3
import os
import sys
import platform
import importlib
import traceback

def die(msg: str, code: int = 1) -> None:
    print(f"[FAIL] {msg}")
    sys.exit(code)

def ok(msg: str) -> None:
    print(f"[OK] {msg}")

def main() -> int:
    # Build dir containing pyng.cpython-311-darwin.so
    build_dir = os.environ.get("PYNG_BUILD_DIR", "cmake-build-pyng-conda")
    build_dir = os.path.abspath(build_dir)

    netlist = """
VDIVIDER.cir
V1 1 0 10
R1 1 2 R = {gauss(3k,0.1,10)}
R2 2 0 R = {gauss(7k,0.1,10)}
.end
"""

    print("=== pyng smoketest ===")
    print(f"Platform: {platform.platform()}")
    print(f"Python:   {sys.version.split()[0]}")
    print(f"Exe:      {sys.executable}")
    print(f"BuildDir: {build_dir}")

    if not os.path.isdir(build_dir):
        die(f"Build directory not found: {build_dir}")

    # Ensure we import the module from the build directory
    sys.path.insert(0, build_dir)

    # numpy check
    try:
        import numpy as np
        ok(f"numpy import ok (version {np.__version__})")
    except Exception as e:
        die(f"numpy import failed: {e}")

    # import extension
    try:
        pyng = importlib.import_module("pyng")
        ok("import pyng ok")
    except Exception:
        print(traceback.format_exc())
        die("import pyng failed (traceback above)")

    # basic API smoke
    try:
        info = pyng.info()
        if not isinstance(info, dict):
            die(f"pyng.info() expected dict, got {type(info)}")
        ok(f"pyng.info() returned dict: {info}")
    except Exception:
        print(traceback.format_exc())
        die("pyng.info() failed")

    try:
        v = pyng.add_ints(2, 3)
        if v != 5:
            die(f"pyng.add_ints(2,3) expected 5, got {v}")
        ok("pyng.add_ints ok")
    except Exception:
        print(traceback.format_exc())
        die("pyng.add_ints failed")

    # numpy array outputs + dtype checks
    try:
        a = pyng.linspace(5)
        if not hasattr(a, "dtype"):
            die("pyng.linspace did not return a numpy-like array")
        if a.shape != (5,):
            die(f"pyng.linspace shape expected (5,), got {a.shape}")
        if str(a.dtype) != "float64":
            die(f"pyng.linspace dtype expected float64, got {a.dtype}")
        if not np.allclose(a, np.array([0, 1, 2, 3, 4], dtype=np.float64)):
            die(f"pyng.linspace values unexpected: {a}")
        ok("pyng.linspace ok (float64, correct shape/values)")
    except Exception:
        print(traceback.format_exc())
        die("pyng.linspace failed")

    try:
        b = pyng.unit_phasors(3)
        if b.shape != (3,):
            die(f"pyng.unit_phasors shape expected (3,), got {b.shape}")
        if str(b.dtype) != "complex64":
            die(f"pyng.unit_phasors dtype expected complex64, got {b.dtype}")
        ok("pyng.unit_phasors ok (complex64, correct shape)")
    except Exception:
        print(traceback.format_exc())
        die("pyng.unit_phasors failed")

    try:
        pyng.ngspice_init()
    except Exception:
        print(traceback.format_exc())
        die("pyng.LoadNetlist failed")

    try:
        if pyng.load_netlist(netlist) == 1:
            pyng.run_command("reset")
            pyng.run_command("op")

    except Exception:
        print(traceback.format_exc())
        die("pyng.LoadNetlist failed")

    ok("All smoketests passed.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
