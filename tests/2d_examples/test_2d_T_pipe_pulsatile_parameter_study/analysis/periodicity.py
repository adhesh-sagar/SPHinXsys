"""Cycle-to-cycle periodicity: has the flow reached the periodic steady state?

This is the standard procedure in pulsatile cardiovascular CFD and it has a name:
the simulation is run for N cycles, the *start-up transient* is discarded, and
periodicity is DEMONSTRATED before the last cycles are analysed. It is not the
same thing as spatial convergence (see convergence.py).

The metric used here is the relative L2 difference between consecutive cycles of
a waveform f, folded onto a common phase grid:

    eps_N = || f_N(phi) - f_{N-1}(phi) ||_2 / || f_N(phi) ||_2

Periodicity is declared at eps_N < 1 %.

The theory that makes this quantitative: the slowest start-up mode of a plane
channel decays at lambda_1 = nu (pi/2R)^2 = (pi^2/4)/tau_visc with
tau_visc = R^2/nu, and tau_visc/T = alpha^2/(2 pi). So the number of cycles to
periodicity scales as ~0.4 alpha^2 - two cycles at alpha = 2, forty at alpha = 10.
Comparing measured eps_N curves against that scaling is the point of this module.
"""

from __future__ import annotations

import numpy as np

from io_utils import Case, interp_periodic


def cycle_waveforms(t, y, case: Case, n_phase: int = 200):
    """Fold EVERY complete pulsation cycle (not just the window) onto a phase grid.

    Returns ``(cycle_numbers, phase, Y)`` with ``Y`` of shape (n_cycles, n_phase).
    Unlike ``io_utils.resample_cycles`` this deliberately covers the whole run,
    because the transient is exactly what we want to look at.
    """
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)

    # drop the ramp: before t_ramp there is no pulsation to fold
    keep = t >= case["t_ramp"]
    t, y = t[keep], y[keep]

    cyc = case.cycle_index(t)
    phase_grid = np.linspace(0.0, 1.0, n_phase, endpoint=False)

    numbers, rows = [], []
    for c in np.unique(cyc):
        sel = cyc == c
        if sel.sum() < 8:
            continue
        ph = case.phase(t[sel])
        order = np.argsort(ph)
        numbers.append(int(c))
        rows.append(interp_periodic(phase_grid, ph[order], y[sel][order]))

    if not rows:
        return np.array([], dtype=int), phase_grid, np.zeros((0, n_phase))
    return np.asarray(numbers), phase_grid, np.vstack(rows)


def epsilon_series(t, y, case: Case, n_phase: int = 200):
    """Relative L2 difference between consecutive cycles.

    Returns ``(cycle_numbers[1:], eps)``: eps[i] compares cycle i+1 with cycle i.
    """
    numbers, _, Y = cycle_waveforms(t, y, case, n_phase)
    if Y.shape[0] < 2:
        return np.array([]), np.array([])

    eps = []
    for i in range(1, Y.shape[0]):
        denom = np.linalg.norm(Y[i])
        if denom <= 0.0:
            eps.append(np.nan)
        else:
            eps.append(np.linalg.norm(Y[i] - Y[i - 1]) / denom)
    return numbers[1:], np.asarray(eps)


def first_periodic_cycle(numbers, eps, tol: float = 0.01):
    """First cycle at which eps drops below ``tol`` AND stays below it."""
    if len(eps) == 0:
        return None
    below = eps < tol
    for i in range(len(below)):
        if below[i:].all():
            return int(numbers[i])
    return None


def assess(case: Case, tol: float = 0.01) -> dict:
    """Full periodicity report for one case, over three independent signals.

    Flow rate, pressure drop and total kinetic energy are checked because they
    probe different things: Q is a bulk kinematic quantity, dp is dominated by
    the near-wall stress, and the kinetic energy integrates the whole field.
    """
    out = {
        "case": case.name,
        "Re": case["Re"],
        "alpha": case["alpha"],
        "A": case["A"],
        "dp": case["dp"],
        "n_cycles": case["n_cycles"],
        "decay_cycles_predicted": case.get("decay_cycles", np.nan),
        "alpha_sq_scaling": 0.4 * case["alpha"] ** 2,
    }

    signals = {}
    t, qi, _, _ = case.all_flow_rates()
    signals["Q_inlet"] = (t, qi)
    t, dp_up = case.pressure_drop("upper")
    signals["dp_upper"] = (t, dp_up)
    t, ke = case.kinetic_energy()
    signals["kinetic_energy"] = (t, ke)

    out["signals"] = {}
    for name, (tt, yy) in signals.items():
        numbers, eps = epsilon_series(tt, yy, case)
        out["signals"][name] = {
            "cycles": numbers,
            "eps": eps,
            "eps_final": float(eps[-1]) if len(eps) else np.nan,
            "first_periodic": first_periodic_cycle(numbers, eps, tol),
        }

    finals = [s["eps_final"] for s in out["signals"].values()]
    out["eps_final_max"] = float(np.nanmax(finals)) if finals else np.nan
    out["periodic"] = bool(out["eps_final_max"] < tol)
    return out


def wcsph_checks(case: Case) -> dict:
    """The two weakly-compressible validity checks, over the analysis window.

    Mach number must stay <= 0.1 and the density deviation <= 1 %, otherwise the
    weak-compressibility assumption underlying WCSPH is violated and the pressure
    field should not be trusted quantitatively.
    """
    t, umax = case.max_speed()
    t2, drho = case.max_density_deviation()

    win = case.in_window(t)
    win2 = case.in_window(t2)
    mach = umax[win] / case["c_f"] if win.any() else umax / case["c_f"]
    dev = drho[win2] if win2.any() else drho

    return {
        "case": case.name,
        "mach_max": float(np.max(mach)),
        "mach_ok": bool(np.max(mach) <= 0.1),
        "density_deviation_max": float(np.max(dev)),
        "density_ok": bool(np.max(dev) <= 0.01),
        "delta_stokes_over_dp": case["delta_stokes"] / case["dp"],
        "stokes_resolved": bool(case.get("stokes_resolved", False)),
    }


def table(cases: dict[str, Case], tol: float = 0.01):
    """Periodicity + WCSPH check table for a whole study, as a list of dicts."""
    rows = []
    for name in sorted(cases):
        c = cases[name]
        try:
            a = assess(c, tol)
            w = wcsph_checks(c)
        except (FileNotFoundError, ValueError, AssertionError) as exc:
            print(f"  {name}: periodicity check skipped ({exc})")
            continue
        rows.append(
            {
                "case": name,
                "Re": c["Re"],
                "alpha": c["alpha"],
                "A": c["A"],
                "dp": c["dp"],
                "n_cycles": c["n_cycles"],
                "cycles_predicted": 0.4 * c["alpha"] ** 2,
                "eps_Q": a["signals"]["Q_inlet"]["eps_final"],
                "eps_dp": a["signals"]["dp_upper"]["eps_final"],
                "eps_KE": a["signals"]["kinetic_energy"]["eps_final"],
                "periodic": a["periodic"],
                "mach_max": w["mach_max"],
                "mach_ok": w["mach_ok"],
                "drho_max": w["density_deviation_max"],
                "density_ok": w["density_ok"],
                "delta_over_dp": w["delta_stokes_over_dp"],
            }
        )
    return rows
