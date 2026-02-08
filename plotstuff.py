#!/usr/bin/env python3
import os
import sys
import importlib
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# Load local-built extension first
# -----------------------------------------------------------------------------
build_dir = os.environ.get("PYNG_BUILD_DIR", "cmake-build-pyng-conda")
build_dir = os.path.abspath(build_dir)
sys.path.insert(0, build_dir)  # ensure our built .so wins

pyng = importlib.import_module("pyng")


# -----------------------------------------------------------------------------
# Plot helpers
# -----------------------------------------------------------------------------
def plot_overlay_mag(pkg, signal_name=None, run_slice=slice(None)):
    """
    Overlay magnitude(dB) of one signal across many runs.
    Works only for uniform-point packages (same n_pts each run).
    """
    assert pkg["uniform"]
    sig_names = pkg["signal_names"]
    sig = pkg["signals"]
    x = pkg["x"][run_slice]

    if signal_name is None:
        j = 0
        title = sig_names[0]
    else:
        # case-insensitive lookup
        sig_names_l = [s.lower() for s in sig_names]
        try:
            j = sig_names_l.index(signal_name.lower())
        except ValueError as e:
            raise RuntimeError(f"signal not found. available={sig_names}") from e
        title = sig_names[j]

    y = sig[run_slice, j, :]
    mag_db = 20 * np.log10(np.abs(y))

    plt.figure()
    for i in range(mag_db.shape[0]):
        plt.semilogx(x[i], mag_db[i])
    plt.xlabel(pkg["x_label"])
    plt.ylabel("Magnitude (dB)")
    plt.title(f"{title} (all runs)")
    plt.show()


# -----------------------------------------------------------------------------
# Unpack / reshape
# -----------------------------------------------------------------------------
def unpack_multirun(sim_results):
    """
    sim_results tuple:
      number_of_runs,
      param_names,
      signal_names,
      run_indices,
      num_datapoints,
      param_values_per_run,
      x_axes_per_run,
      x_label,
      signals

    Returns a dict with either dense arrays (uniform points) or ragged lists.
    """
    (n_runs,
     param_names,
     signal_names,
     run_indices,
     num_datapoints,
     param_values_per_run,
     x_axes_per_run,
     x_label,
     signals) = sim_results

    n_runs = int(n_runs)
    n_params = len(param_names)
    n_signals = len(signal_names)

    run_indices = np.asarray(run_indices, dtype=np.uint64)
    num_datapoints = np.asarray(num_datapoints, dtype=np.uint64)

    # Params: expect shape (n_runs, n_params)
    params = np.asarray(param_values_per_run, dtype=np.float64)
    if params.shape != (n_runs, n_params):
        raise ValueError(f"params shape {params.shape} != ({n_runs},{n_params})")

    # Decide uniform vs ragged
    uniform = bool(np.all(num_datapoints == num_datapoints[0]))
    if uniform:
        n_pts = int(num_datapoints[0])

        x = np.asarray(x_axes_per_run, dtype=np.float64)
        if x.shape != (n_runs, n_pts):
            raise ValueError(f"x shape {x.shape} != ({n_runs},{n_pts})")

        sig_flat = np.asarray(signals, dtype=np.complex128)

        # Normalize to (n_runs, n_signals, n_pts)
        if sig_flat.ndim == 2 and sig_flat.shape == (n_runs, n_pts * n_signals):
            # layout A: per-run row, concatenated signals
            sig = sig_flat.reshape(n_runs, n_signals, n_pts)

        elif sig_flat.ndim == 2 and sig_flat.shape == (n_runs * n_signals, n_pts):
            # layout B: stacked vectors:
            #   [run0_sig0, run0_sig1, ..., run1_sig0, run1_sig1, ...]
            sig = sig_flat.reshape(n_runs, n_signals, n_pts)

        elif sig_flat.ndim == 2 and sig_flat.shape == (n_runs, n_pts) and n_signals == 1:
            # layout C: single signal already per-run
            sig = sig_flat[:, None, :]  # add signal axis

        elif sig_flat.ndim == 3 and sig_flat.shape == (n_runs, n_signals, n_pts):
            # layout D: already correct
            sig = sig_flat

        else:
            raise ValueError(
                f"signals shape {sig_flat.shape} not understood for "
                f"n_runs={n_runs}, n_signals={n_signals}, n_pts={n_pts}"
            )

        return {
            "n_runs": n_runs,
            "param_names": list(param_names),
            "signal_names": list(signal_names),
            "run_indices": run_indices,
            "num_datapoints": num_datapoints,
            "params": params,     # (n_runs, n_params)
            "x": x,               # (n_runs, n_pts)
            "x_label": str(x_label),
            "signals": sig,       # (n_runs, n_signals, n_pts)
            "uniform": True,
        }

    # Ragged case (adaptive transient etc.)
    x_list = [np.asarray(x_axes_per_run[i], dtype=np.float64) for i in range(n_runs)]

    # signals as list-of-dicts keyed by signal name
    if isinstance(signals, (list, tuple)) and len(signals) == n_runs and isinstance(signals[0], (list, tuple)):
        # nested: signals[run][sig]
        runs_signal_dicts = [
            {signal_names[j]: np.asarray(signals[i][j], dtype=np.complex128) for j in range(n_signals)}
            for i in range(n_runs)
        ]
    else:
        # flat: [run0_sig0, run0_sig1, ..., run1_sig0, ...]
        flat = list(signals)
        if len(flat) != n_runs * n_signals:
            raise ValueError(f"ragged flat signals len={len(flat)} != n_runs*n_signals={n_runs*n_signals}")
        runs_signal_dicts = []
        idx = 0
        for i in range(n_runs):
            d = {}
            for j in range(n_signals):
                d[signal_names[j]] = np.asarray(flat[idx], dtype=np.complex128)
                idx += 1
            runs_signal_dicts.append(d)

    return {
        "n_runs": n_runs,
        "param_names": list(param_names),
        "signal_names": list(signal_names),
        "run_indices": run_indices,
        "num_datapoints": num_datapoints,
        "params": params,
        "x_list": x_list,                  # list[np.ndarray]
        "x_label": str(x_label),
        "signals_list": runs_signal_dicts,  # list[dict[str,np.ndarray]]
        "uniform": False,
    }


# -----------------------------------------------------------------------------
# Optional: deep sizeof helper (Python object overhead, not NumPy buffer bytes)
# -----------------------------------------------------------------------------
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


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main() -> int:
#     netlist = """
# VDIVIDER.cir
# V1 1 0 10 AC=1
# R1 1 2 R = {unif(3k,0.05)}
# C1 2 0 C = {unif(1u,0.47)}
# R2 2 0 R = {unif(7k,0.05)}
# .end
# """
    netlist = """
ADCREF.cir
X1 N003 x1n x1out ISL70444_FREQ
R1 Vout x1nn {unif(4.99k,0.05)}
R2 0 x1nn {unif(4.99k,0.05)}
R3 x1out N004 100
R4 N005 Vout {unif(10,0.024)}
R5 N001 N002 33.2
V1 N001 0 15
V2 x1n x1nn 0 AC 1
Q1 N002 N004 N005 0 NPN
V3 N003 0 2.5
R6 x1out N006 {pow(10, aunif(3,3))}
C1 N006 x1nn {pow(10, aunif(-7,2.5))}
*R6 x1out N006 {unif(50k,0.99)}
*C1 N006 x1nn {unif(1u,0.99999)}
C2 Vout 0 10u
*I1 Vout 0 10m
*I2 Vout 0 PULSE(0 10m 1m 10u 10u 50m)

* block symbol definitions
.subckt ISL70444_FREQ 1 2 3
*.param f1=1
*.param f2=19meg
*.param AOL=1000000
*.param rout=10 rin=1G
R1 N001 0 1
*C1 N001 0 {1/(2*3.14159*f1)}
C1 N001 0 0.159154943
R2 N002 0 1
*C2 N002 0 {1/(2*3.14159*f2)}
C2 N002 0 0.000000008
*G1 0 N001 1 2 {AOL}
G1 0 N001 1 2 1000000
G2 0 N002 N001 0 1
*R3 1 2 {rin}
R3 1 2 1G
R4 N003 0 1
*C3 N003 0 {1/(2*3.14159*f2)}
C3 N003 0 0.000000008
G3 0 N003 N002 0 1
E1 N004 0 N003 0 1
*R5 3 N004 {rout}
R5 3 N004 10
.ends ISL70444_FREQ

.model NPN NPN
.end
"""
    pyng.ngspice_init()

    # This is your C++ hard-coded prototype entry point:
    sim_results = pyng.multirun_proto(netlist, 100)

    pkg = unpack_multirun(sim_results)

    # ---- Example plot: H = v(2)/v(1), overlay all runs ----
    assert pkg["uniform"], "this plotting path assumes uniform AC sweep points"

    sig_names = pkg["signal_names"]
    sig_names_l = [s.lower() for s in sig_names]

    out_name = "v(x1nn)"
    in_name = "v(x1n)"

    try:
        j_out = sig_names_l.index(out_name.lower())
        j_in = sig_names_l.index(in_name.lower())
    except ValueError as e:
        raise RuntimeError(f"needed signals not found. available={sig_names}") from e

    x = pkg["x"]                    # (n_runs, n_pts)
    Vout = pkg["signals"][:, j_out, :]  # (n_runs, n_pts)
    Vin  = pkg["signals"][:, j_in, :]   # (n_runs, n_pts)

    H = Vout / Vin

    mag_db = 20 * np.log10(np.abs(H))
    ph_deg = np.angle(H, True)

    plt.figure()
    for i in range(pkg["n_runs"]):
        plt.semilogx(x[i], mag_db[i])
        plt.semilogx(x[i], ph_deg[i])
    plt.xlabel(pkg["x_label"])
    plt.ylabel("|Vout/Vin| (dB)")
    plt.title("H = v(2)/v(1) overlay (all runs)")
    plt.show()

    # ---- Optional: quick sanity overlay of one raw signal ----
    # plot_overlay_mag(pkg, signal_name="v(2)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
