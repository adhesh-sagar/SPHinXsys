"""Particle-resolution convergence: Richardson extrapolation and the GCI.

This is the standard CFD verification procedure (Roache 1994/1998; ASME V&V 20).
For a functional f computed on three resolutions with a constant refinement ratio
r = dp_coarse / dp_fine, the observed order of convergence is

    q = ln( (f_3 - f_2) / (f_2 - f_1) ) / ln r

(1 = finest), the Richardson-extrapolated exact value is

    f_exact ~ f_1 + (f_1 - f_2) / (r^q - 1),

and the Grid Convergence Index on the finest pair is

    GCI_fine = Fs * |(f_1 - f_2) / f_1| / (r^q - 1),   Fs = 1.25 for three grids.

GCI is reported as a percentage and is interpreted as a numerical uncertainty
band on the finest result.

Two functionals are tracked, chosen because they probe different parts of the
solution:

* ``Q_mean``   cycle-mean inlet flow rate - a bulk, well-behaved quantity
* ``tawss``    peak time-averaged wall shear stress - near-wall, the most
                resolution-sensitive quantity in the whole study

Plus a direct analytical error for the A = 0 case (steady Poiseuille), which does
not need extrapolation at all because the exact answer is known.
"""

from __future__ import annotations

import numpy as np

from io_utils import Case, resample_cycles
import wss as wss_mod


# ---------------------------------------------------------------------------
# functionals
# ---------------------------------------------------------------------------


def functional_Q_mean(case: Case) -> float:
    """Cycle-mean inlet volume flow rate over the analysis window."""
    t, Q = case.flow_rate("inlet")
    phase, Y = resample_cycles(t, Q, case)
    if Y.shape[0] == 0:
        return np.nan
    return float(Y.mean())


def functional_peak_tawss(case: Case) -> float:
    """Largest TAWSS anywhere on the walls."""
    try:
        res = wss_mod.wall_metrics(case)
    except (FileNotFoundError, ValueError) as exc:
        print(f"  {case.name}: peak TAWSS unavailable ({exc})")
        return np.nan
    vals = [np.nanmax(v["tawss"]) for v in res.values() if len(v["tawss"])]
    return float(np.nanmax(vals)) if vals else np.nan


def functional_split(case: Case) -> float:
    """Cycle-mean flow split Q_upper / (Q_upper + Q_lower)."""
    t, qi, qu, ql = case.all_flow_rates()
    _, YU = resample_cycles(t, qu, case)
    _, YL = resample_cycles(t, ql, case)
    if YU.shape[0] == 0 or YL.shape[0] == 0:
        return np.nan
    up, lo = YU.mean(), YL.mean()
    return float(up / (up + lo)) if (up + lo) != 0 else np.nan


def functional_tawss_p95(case: Case) -> float:
    """95th-percentile TAWSS over all wall stations.

    Preferred over the raw maximum for the convergence study. ``peak_tawss`` is a
    single-point maximum over a few hundred probes, so one noisy probe moves it:
    across the resolution series it came out 0.208 / 0.311 / 0.188, i.e.
    non-monotone, purely from which probe happened to be the outlier. A high
    percentile measures the same near-wall physics without being hostage to one
    sample, which is standard practice for reporting peak wall shear.
    """
    try:
        res = wss_mod.wall_metrics(case)
    except (FileNotFoundError, ValueError, AssertionError) as exc:
        print(f"  {case.name}: TAWSS percentile unavailable ({exc})")
        return np.nan
    vals = np.concatenate([v["tawss"] for v in res.values() if len(v["tawss"])])
    vals = vals[np.isfinite(vals)]
    return float(np.percentile(vals, 95)) if vals.size else np.nan


FUNCTIONALS = {
    "Q_mean": functional_Q_mean,
    "tawss_p95": functional_tawss_p95,
    "peak_tawss": functional_peak_tawss,
    "flow_split": functional_split,
}


# ---------------------------------------------------------------------------
# Richardson / GCI
# ---------------------------------------------------------------------------


def observed_order(f1: float, f2: float, f3: float,
                   r21: float, r32: float | None = None) -> float:
    """Observed order of convergence from three successively coarser values.

    ``f1`` is the finest. ``r21 = dp2/dp1`` and ``r32 = dp3/dp2``; if ``r32`` is
    omitted a constant ratio is assumed.

    For a constant ratio this reduces to the familiar
    ``q = ln(eps32/eps21) / ln(r)``. For UNEQUAL ratios it solves the standard
    ASME V&V 20 / Celik et al. (2008) relation by fixed-point iteration:

        q     = |ln|eps32/eps21| + Q(q)| / ln(r21)
        Q(q)  = ln( (r21^q - s) / (r32^q - s) ),    s = sign(eps32/eps21)

    This matters in practice: SPH resolutions are constrained to those where DH/dp
    is an integer, and the affordable integer ladder (60/40/30 layers) does not
    have a constant ratio. Forcing a constant-ratio formula onto it returns a
    negative, meaningless order.

    Returns NaN when the series cannot support Richardson extrapolation:
    oscillatory convergence (the differences change sign), a non-positive order
    (coarse-grid differences smaller than fine-grid ones), or an implausibly high
    order (> 4, not believable for a second-order scheme). Each of those means the
    series is outside the asymptotic range, which is a result worth reporting
    rather than a number to fabricate.
    """
    if r32 is None:
        r32 = r21
    e21, e32 = f2 - f1, f3 - f2
    if e21 == 0 or e32 == 0:
        return np.nan
    ratio = e32 / e21
    if ratio <= 0:
        return np.nan  # oscillatory: the sequence changes direction

    s = 1.0  # sign(ratio), and ratio > 0 here by the test above
    q = float(np.log(ratio) / np.log(r21))  # constant-ratio starting guess
    if abs(r21 - r32) > 1e-12:
        for _ in range(200):
            denom = r32**q - s
            if denom == 0 or (r21**q - s) / denom <= 0:
                return np.nan
            q_new = abs(np.log(abs(ratio)) + np.log((r21**q - s) / denom)) / np.log(r21)
            if not np.isfinite(q_new):
                return np.nan
            if abs(q_new - q) < 1e-10:
                q = q_new
                break
            q = q_new

    q = float(q)
    if not np.isfinite(q) or q <= 0.0 or q > 4.0:
        return np.nan
    return q


def richardson(f1: float, f2: float, r: float, q: float) -> float:
    """Extrapolated zero-spacing value from the two finest results."""
    if not np.isfinite(q):
        return np.nan
    return float(f1 + (f1 - f2) / (r**q - 1.0))


def gci(f1: float, f2: float, r: float, q: float, Fs: float = 1.25) -> float:
    """Grid Convergence Index on the finest pair, in percent."""
    if not np.isfinite(q) or f1 == 0:
        return np.nan
    return float(100.0 * Fs * abs((f1 - f2) / f1) / (r**q - 1.0))


def study(cases: dict[str, Case], series: list[str] | None = None):
    """Convergence report for a resolution series.

    ``series`` names the cases in the series; if omitted, every case sharing the
    modal (Re, alpha, A, branch_ratio, stenosis) combination is used, which picks
    out group C plus its group-A anchor automatically.
    """
    if series is None:
        # group by the physical case, ignoring dp
        buckets: dict[tuple, list[str]] = {}
        for name, c in cases.items():
            key = (c["Re"], c["alpha"], c["A"], c["branch_ratio"], c["stenosis"],
                   c["cfl_scale"])
            buckets.setdefault(key, []).append(name)
        # the resolution series is the largest bucket with more than one dp
        best, best_key = [], None
        for key, names in buckets.items():
            dps = {cases[n]["dp"] for n in names}
            if len(dps) > len(best):
                best, best_key = names, key
        series = best
        if len(series) < 3:
            print("  convergence: fewer than 3 resolutions found, skipping")
            return {
                "series": series,
                "dp": [cases[n]["dp"] for n in series],
                "functionals": {},
            }

    # sort fine -> coarse
    series = sorted(series, key=lambda n: cases[n]["dp"])
    dps = [cases[n]["dp"] for n in series]

    out = {"series": series, "dp": dps, "functionals": {}}
    for fname, fn in FUNCTIONALS.items():
        vals = [fn(cases[n]) for n in series]
        entry = {"values": vals, "dp": dps, "cases": series}

        if len(series) >= 3 and all(np.isfinite(v) for v in vals[:3]):
            f1, f2, f3 = vals[0], vals[1], vals[2]
            r_a = dps[1] / dps[0]
            r_b = dps[2] / dps[1]
            entry["r_fine_pair"] = r_a
            entry["r_coarse_pair"] = r_b
            entry["r_uniform"] = bool(abs(r_a - r_b) / r_a < 0.05)
            # The unequal-ratio (Celik) form is solved directly, so a non-constant
            # ladder no longer has to be approximated by a geometric-mean ratio.
            r = r_a
            entry["r"] = r
            q = observed_order(f1, f2, f3, r_a, r_b)
            entry["order"] = q
            entry["extrapolated"] = richardson(f1, f2, r, q)
            entry["gci_percent"] = gci(f1, f2, r, q)
            entry["monotone"] = bool(np.isfinite(q))
        out["functionals"][fname] = entry
    return out


# ---------------------------------------------------------------------------
# direct analytical error (steady case, no extrapolation needed)
# ---------------------------------------------------------------------------


def poiseuille_error(case: Case):
    """L2 / Linf error of the inlet profile against the exact parabola.

    Only meaningful for A = 0, where the exact steady solution is
    u(y) = 1.5 U (1 - (y/R)^2) with U the mean velocity. Because the exact answer
    is known, this is a validation error rather than an extrapolated estimate, and
    its decay with dp is the cleanest convergence evidence in the report.
    """
    if case["A"] != 0.0:
        raise ValueError(f"{case.name}: A = {case['A']}, not the steady case")

    R = case["R"]
    t, coord, u_along, _ = case.profile("inlet")
    win = case.in_window(t)
    u_sim = u_along[win].mean(axis=0)

    t2, Q = case.flow_rate("inlet")
    U = float(np.mean(Q[case.in_window(t2)])) / (2.0 * R)

    y = coord - R
    u_ana = 1.5 * U * np.maximum(0.0, 1.0 - (y / R) ** 2)

    denom = np.linalg.norm(u_ana)
    return {
        "case": case.name,
        "dp": case["dp"],
        "y_over_R": y / R,
        "u_sim": u_sim,
        "u_ana": u_ana,
        "U_mean": U,
        "L2_rel": float(np.linalg.norm(u_sim - u_ana) / denom) if denom else np.nan,
        "Linf_rel": float(
            np.max(np.abs(u_sim - u_ana)) / max(np.max(np.abs(u_ana)), 1e-30)
        ),
    }
