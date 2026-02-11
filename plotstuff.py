#!/usr/bin/env python3
import os
import sys
import importlib
import numpy as np
import matplotlib.pyplot as plt
from typing import NamedTuple

# -----------------------------------------------------------------------------
# Load local-built extension first
# -----------------------------------------------------------------------------
build_dir = os.environ.get("PYNG_BUILD_DIR", "cmake-build-pyng-conda")
build_dir = os.path.abspath(build_dir)
sys.path.insert(0, build_dir)  # ensure our built .so wins

pyng = importlib.import_module("pyng")


def _normalize_angle_deg(a: np.ndarray) -> np.ndarray:
    """Map angles to (-180, 180]."""
    return (a + 180.0) % 360.0 - 180.0

def _unwrap_deg(ph_deg: np.ndarray) -> np.ndarray:
    """Unwrap phase in degrees."""
    return np.rad2deg(np.unwrap(np.deg2rad(ph_deg)))

def _interp_logf(f1, f2, t):
    """Interpolate frequency in log domain."""
    lf1 = np.log(f1)
    lf2 = np.log(f2)
    return float(np.exp(lf1 + t * (lf2 - lf1)))

def _first_gain_crossover(freq, mag_db, phase_deg_unwrapped):
    """
    Find FIRST (lowest-f) gain crossover where mag_db crosses 0.
    Returns (fc, phase_at_fc) or None.
    """
    n = len(freq)
    for i in range(n - 1):
        d1 = mag_db[i]
        d2 = mag_db[i + 1]
        f1 = freq[i]
        f2 = freq[i + 1]
        p1 = phase_deg_unwrapped[i]
        p2 = phase_deg_unwrapped[i + 1]

        # exact hits
        if d1 == 0.0 and np.isfinite(p1) and f1 > 0:
            return float(f1), float(p1)
        if d2 == 0.0 and np.isfinite(p2) and f2 > 0:
            return float(f2), float(p2)

        # crossing
        if (d1 > 0.0 and d2 < 0.0) or (d1 < 0.0 and d2 > 0.0):
            denom = (d2 - d1)
            if denom == 0.0:
                continue
            t = (0.0 - d1) / denom  # fraction from i -> i+1 in dB domain
            if not (0.0 <= t <= 1.0):
                continue
            if f1 <= 0.0 or f2 <= 0.0:
                continue
            fc = _interp_logf(f1, f2, t)
            phc = float(p1 + t * (p2 - p1))
            return fc, phc

    return None

def _first_phase_crossover(freq, phase_deg_unwrapped, mag_db, phase_cross_deg):
    """
    Find FIRST (lowest-f) phase crossover where phase crosses phase_cross_deg.
    Returns (fc, mag_db_at_fc) or None.
    """
    n = len(freq)
    for i in range(n - 1):
        p1 = phase_deg_unwrapped[i]
        p2 = phase_deg_unwrapped[i + 1]
        f1 = freq[i]
        f2 = freq[i + 1]
        d1 = mag_db[i]
        d2 = mag_db[i + 1]

        # exact hits
        if p1 == phase_cross_deg and np.isfinite(d1) and f1 > 0:
            return float(f1), float(d1)
        if p2 == phase_cross_deg and np.isfinite(d2) and f2 > 0:
            return float(f2), float(d2)

        # crossing
        if (p1 > phase_cross_deg and p2 < phase_cross_deg) or (p1 < phase_cross_deg and p2 > phase_cross_deg):
            denom = (p2 - p1)
            if denom == 0.0:
                continue
            t = (phase_cross_deg - p1) / denom  # fraction in phase domain
            if not (0.0 <= t <= 1.0):
                continue
            if f1 <= 0.0 or f2 <= 0.0:
                continue
            fc = _interp_logf(f1, f2, t)
            dbc = float(d1 + t * (d2 - d1))
            return fc, dbc

    return None

def _fit_line_y_vs_logx(x, y):
    """
    Fit y = a*log10(x) + b. Returns (x_sorted, y_fit, r2) or (None,None,None).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    m = np.isfinite(x) & np.isfinite(y) & (x > 0.0)
    x = x[m]
    y = y[m]
    if x.size < 2:
        return None, None, None

    lx = np.log10(x)
    # Guard against all identical x -> singular fit
    if np.allclose(lx, lx[0]):
        return None, None, None

    # Least squares on [lx, 1]
    A = np.vstack([lx, np.ones_like(lx)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    a, b = coef

    yhat = a * lx + b
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else float("nan")

    order = np.argsort(x)
    x_sorted = x[order]
    y_fit = (a * np.log10(x_sorted) + b)
    return x_sorted, y_fit, r2

def plot_stability_margins_vs_params(
        pkg: dict,
        out_signal: str,
        in_signal: str,
        phase_cross_deg: float = 0.0,
        use_first_crossover: bool = True,
):
    """
    Expects pkg from your unpack_multirun() (uniform=True), with:
      pkg["signals"] shape (n_runs, n_signals, n_pts)
      pkg["x"] shape (n_runs, n_pts)
      pkg["params"] shape (n_runs, n_params)
      pkg["param_names"], pkg["signal_names"], pkg["x_label"]

    Computes:
      PM (deg): phase_cross_deg - phase_at_gain_crossover  (wrapped to (-180,180])
      GM (dB): -mag_db_at_phase_crossover

    Skips runs with no crossover(s).
    Makes separate scatter plots (PM, GM) per parameter with x-axis log scale, fit line, R^2 annotation.
    """
    if not pkg.get("uniform", False):
        raise ValueError("This function currently expects uniform=True package (dense arrays).")

    sig_names = list(pkg["signal_names"])
    try:
        j_out = sig_names.index(out_signal)
        j_in  = sig_names.index(in_signal)
    except ValueError as e:
        raise ValueError(f"Signal not found. available={sig_names}") from e

    x_all = np.asarray(pkg["x"], dtype=float)                 # (n_runs, n_pts)
    sig   = np.asarray(pkg["signals"], dtype=np.complex128)   # (n_runs, n_sig, n_pts)
    params = np.asarray(pkg["params"], dtype=float)           # (n_runs, n_params)
    pnames = list(pkg["param_names"])
    n_runs = int(pkg["n_runs"])

    pm_vals = []
    gm_vals = []
    pm_runmask = np.zeros(n_runs, dtype=bool)
    gm_runmask = np.zeros(n_runs, dtype=bool)

    for r in range(n_runs):
        f = x_all[r, :]
        H = sig[r, j_out, :] / sig[r, j_in, :]

        # basic cleaning
        m = np.isfinite(f) & (f > 0.0) & np.isfinite(H.real) & np.isfinite(H.imag)
        f = f[m]
        H = H[m]
        if f.size < 2:
            continue

        mag_db = 20.0 * np.log10(np.abs(H))
        ph_deg = np.angle(H, deg=True)
        ph_deg = _unwrap_deg(ph_deg)

        # Gain crossover -> PM
        gc = _first_gain_crossover(f, mag_db, ph_deg)
        if gc is not None:
            _, ph_at_gc = gc
            pm = phase_cross_deg - ph_at_gc
            pm = float(_normalize_angle_deg(np.array(pm)))
            pm_vals.append((r, pm))
            pm_runmask[r] = True

        # Phase crossover -> GM
        pc = _first_phase_crossover(f, ph_deg, mag_db, phase_cross_deg)
        if pc is not None:
            _, mag_at_pc_db = pc
            gm = -mag_at_pc_db  # +dB margin if magnitude is below 0dB at phase crossover
            gm_vals.append((r, float(gm)))
            gm_runmask[r] = True

    # Convert to arrays aligned to runs (with NaNs for missing)
    pm_arr = np.full(n_runs, np.nan, dtype=float)
    gm_arr = np.full(n_runs, np.nan, dtype=float)
    for r, v in pm_vals:
        pm_arr[r] = v
    for r, v in gm_vals:
        gm_arr[r] = v

    # Plot per parameter
    for j, pname in enumerate(pnames):
        xparam = params[:, j]

        # --- Phase margin plot ---
        x_pm = xparam[pm_runmask]
        y_pm = pm_arr[pm_runmask]
        if x_pm.size >= 1:
            plt.figure()
            plt.scatter(x_pm, y_pm)
            plt.xscale("log")
            plt.xlabel(pname)
            plt.ylabel(f"Phase margin (deg) vs {phase_cross_deg:.0f}° crossing")
            plt.title(f"Phase margin vs {pname}")

            x_fit, y_fit, r2 = _fit_line_y_vs_logx(x_pm, y_pm)
            if x_fit is not None:
                plt.plot(x_fit, y_fit)
                plt.text(0.05, 0.95, f"$R^2$ = {r2:.3f}", transform=plt.gca().transAxes, va="top")

            plt.grid(True, which="both", linestyle=":")
            plt.show()

        # --- Gain margin plot ---
        x_gm = xparam[gm_runmask]
        y_gm = gm_arr[gm_runmask]
        if x_gm.size >= 1:
            plt.figure()
            plt.scatter(x_gm, y_gm)
            plt.xscale("log")
            plt.xlabel(pname)
            plt.ylabel(f"Gain margin (dB) at {phase_cross_deg:.0f}° phase crossover")
            plt.title(f"Gain margin vs {pname}")

            x_fit, y_fit, r2 = _fit_line_y_vs_logx(x_gm, y_gm)
            if x_fit is not None:
                plt.plot(x_fit, y_fit)
                plt.text(0.05, 0.95, f"$R^2$ = {r2:.3f}", transform=plt.gca().transAxes, va="top")

            plt.grid(True, which="both", linestyle=":")
            plt.show()

    return {
        "pm_deg": pm_arr, "pm_mask": pm_runmask,
        "gm_db": gm_arr, "gm_mask": gm_runmask,
    }

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
*R6 x1out N006 {pow(10, aunif(3,3))}
R6 x1out N006 50k
*R6 x1out N006 10
*C1 N006 x1nn {pow(10, aunif(-7,2.5))}
C1 N006 x1nn 1n
*C2 Vout 0 10u
*C2 Vout 0 {pow(10, aunif(-6, 1))}
C2 Vout 0 700n
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

.model NPN NPN(Is=14.34f Xti=3 Eg=1.11 Vaf=74.03 Bf=255.9 Ne=1.307
+ Ise=14.34f Ikf=.2847 Xtb=1.5 Br=6.092 Nc=2 Isc=0 Ikr=0 Rc=1 Cjc=7.306p
+ Mjc=.3416 Vjc=.75 Fc=.5 Cje=22.01p Mje=.377 Vje=.75 Tr=46.91n Tf=411.1p
+ Itf=.6Vtf=1.7 Xtf=3 Rb=10 Vceo=40)

.end
"""
    pyng.ngspice_init()

    # This is your C++ hard-coded prototype entry point:
    # sim_results = pyng.multirun_proto(netlist, 500)
    run_command = "ac dec 10 0.1 1meg"
    signal_labels = ["v(x1n)", "v(x1nn)"]
    params = [("r1", "r"), ("r2", "r"), ("r4", "r"), ("r6", "r"), ("c1", "c"), ("c2", "c")]
    # params = [("r6", "r")]
    sim_results = pyng.montecarlo(netlist, signal_labels, run_command, params, "frequency", 500)
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
    plt.title(f"H = {out_name}/{in_name} overlay (all runs)")
    plt.show()

    results = plot_stability_margins_vs_params(
        pkg,
        out_signal="v(x1nn)",
        in_signal="v(x1n)",
        phase_cross_deg=0.0,
    )

    # ---- Optional: quick sanity overlay of one raw signal ----
    # plot_overlay_mag(pkg, signal_name="v(2)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
