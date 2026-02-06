#!/usr/bin/env python3
import os
import sys
import platform
import importlib
import traceback
import pprint
from pickletools import uint8

import numpy as np
import matplotlib.pyplot as plt

def unpack_multirun_tuple(sim_results):
    n_runs, param_names, signal_names, run_ids, npts, param_vals, signals = sim_results

    n_signals = len(signal_names)
    assert len(signals) == n_runs * n_signals, "signals length doesn't match runs*signals"
    assert len(param_vals) == n_runs * len(param_names), "param_vals length doesn't match runs*params"

    # params: (runs, params)
    params = np.asarray(param_vals, dtype=np.complex128).reshape(n_runs, len(param_names))

    # signals: list of runs, each run is dict name->ndarray
    runs = []
    idx = 0
    for r in range(n_runs):
        run = {}
        for sname in signal_names:
            v = np.asarray(signals[idx], dtype=np.complex128)
            idx += 1
            run[sname] = v
        runs.append(run)

    return {
        "n_runs": n_runs,
        "param_names": param_names,
        "signal_names": signal_names,
        "run_ids": np.asarray(run_ids, dtype=np.int64),
        "npts": np.asarray(npts, dtype=np.int64),
        "params": params,
        "runs": runs,   # ragged per-run arrays
    }
def deep_getsizeof(o, seen=None):
    if seen is None:
        seen = set()
    oid = id(o)
    if oid in seen:
        return 0
    seen.add(oid)
    size = sys.getsizeof(o)

    if isinstance(o, dict):
        size += sum(deep_getsizeof(k, seen) + deep_getsizeof(v, seen) for k, v in o.items())
    elif isinstance(o, (list, tuple, set, frozenset)):
        size += sum(deep_getsizeof(i, seen) for i in o)
    return size

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
V1 1 0 10 AC=1
*R1 1 2 R = {gauss(3k,0.1,10)}
R1 1 2 R = {unif(3k,0.01)}
C1 2 0 C = {unif(1u,0.10)}
*R2 2 0 R = {gauss(7k,0.1,10)}
R2 2 0 R = {unif(7k,0.01)}
*R2 2 0 1e12
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
    # try:
    #     info = pyng.info()
    #     if not isinstance(info, dict):
    #         die(f"pyng.info() expected dict, got {type(info)}")
    #     ok(f"pyng.info() returned dict: {info}")
    # except Exception:
    #     print(traceback.format_exc())
    #     die("pyng.info() failed")

    # try:
    #     v = pyng.add_ints(2, 3)
    #     if v != 5:
    #         die(f"pyng.add_ints(2,3) expected 5, got {v}")
    #     ok("pyng.add_ints ok")
    # except Exception:
    #     print(traceback.format_exc())
    #     die("pyng.add_ints failed")

    # numpy array outputs + dtype checks
    # try:
    #     a = pyng.linspace(5)
    #     if not hasattr(a, "dtype"):
    #         die("pyng.linspace did not return a numpy-like array")
    #     if a.shape != (5,):
    #         die(f"pyng.linspace shape expected (5,), got {a.shape}")
    #     if str(a.dtype) != "float64":
    #         die(f"pyng.linspace dtype expected float64, got {a.dtype}")
    #     if not np.allclose(a, np.array([0, 1, 2, 3, 4], dtype=np.float64)):
    #         die(f"pyng.linspace values unexpected: {a}")
    #     ok("pyng.linspace ok (float64, correct shape/values)")
    # except Exception:
    #     print(traceback.format_exc())
    #     die("pyng.linspace failed")

    # try:
    #     b = pyng.unit_phasors(3)
    #     if b.shape != (3,):
    #         die(f"pyng.unit_phasors shape expected (3,), got {b.shape}")
    #     if str(b.dtype) != "complex64":
    #         die(f"pyng.unit_phasors dtype expected complex64, got {b.dtype}")
    #     ok("pyng.unit_phasors ok (complex64, correct shape)")
    # except Exception:
    #     print(traceback.format_exc())
    #     die("pyng.unit_phasors failed")

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

    # try:
    #     my_tuple = pyng.get_a_tuple()
    #     print(my_tuple)
    # except Exception:
    #     print(traceback.format_exc())
    #     die("pyng.get_a_tuple failed")

    try:
        # sim_results = pyng.multirun_proto(netlist, 100)
        # number_of_runs, param_names, signal_names, run_indices, num_datapoints, param_values_per_run, x_axes_per_run, signals = sim_results
        number_of_runs, param_names, signal_names, run_indices, num_datapoints, param_values_per_run, x_axes_per_run, x_label, signals = pyng.multirun_proto(netlist, 100)
        run_indices = np.array(run_indices, dtype=np.uint64)
        num_datapoints = np.array(num_datapoints, dtype=np.uint64)
        param_values_per_run = np.array(param_values_per_run, dtype=np.float64)
        x_axes_per_run = np.array(x_axes_per_run, dtype=np.float64)
        signals = np.array(signals, dtype=np.complex128)
        #some_array = np.block([run_indices, x_axes_per_run])
        print(f"[number_of_runs]: {deep_getsizeof(number_of_runs):,} Bytes")
        print(f"[param_names]: {deep_getsizeof(param_names):,} Bytes")
        print(f"[signal_names]: {deep_getsizeof(signal_names):,} Bytes")
        print(f"[run_indices]: {deep_getsizeof(run_indices):,} Bytes")
        print(f"[num_datapoints]: {deep_getsizeof(num_datapoints):,} Bytes")
        print(f"[param_values_per_run]: {deep_getsizeof(param_values_per_run):,} Bytes")
        print(f"[x_axes]: {deep_getsizeof(x_axes_per_run):,} Bytes")
        print(f"[x_label]: {deep_getsizeof(x_label)} Bytes")
        print(f"[signals]: {deep_getsizeof(signals):,} Bytes")
        # print(f"Total size: {deep_getsizeof(sim_results):,} Bytes")
        # freqs_30 = x_axes[30-1,:]
        # freq_30_51 = freqs_30[51-1]
        # print(f"{freq_30_51}")
    except Exception:
        print(traceback.format_exc())
        die("pyng.multirun_proto failed")

    ok("All smoketests passed.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
