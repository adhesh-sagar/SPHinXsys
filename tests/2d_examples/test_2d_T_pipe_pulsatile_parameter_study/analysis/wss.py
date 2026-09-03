"""Wall shear stress and the standard hemodynamic wall metrics.

DEFINITIONS
-----------
The viscous traction on a wall with inward unit normal n is

    t = mu (grad u + grad u^T) . n

and the wall shear stress is its tangential part

    tau_w = t - (t . n) n.

From the WSS waveform over one cycle the three metrics used in essentially every
cardiovascular CFD paper follow:

    TAWSS = (1/T) int |tau_w| dt                      time-averaged WSS magnitude
    OSI   = 0.5 * (1 - |int tau_w dt| / int |tau_w| dt)   oscillatory shear index
    RRT   = 1 / ((1 - 2 OSI) * TAWSS)                  relative residence time

OSI needs the VECTOR WSS, not just its magnitude: it measures how much the shear
direction reverses over the cycle. It is 0 for unidirectional shear and 0.5 for
shear that is perfectly balanced forwards and backwards.

Low TAWSS combined with high OSI (equivalently, high RRT) marks the classical
atherosclerosis-prone regions - in a bifurcation these are the flow-divider apex
and the outer walls just downstream of the junction (Ku & Giddens 1983;
Zarins & Glagov 1985).

EXTRACTION FROM SPH
-------------------
The velocity gradient is interpolated at probe lines placed at two offsets
(1 dp and 2 dp) from each wall, then linearly extrapolated to the wall:

    tau_wall = 2 tau(1 dp) - tau(2 dp).

Evaluating exactly at the wall is avoided on purpose: SPH interpolation there is
one-sided (the observer only sees fluid particles), so it is biased. Two interior
offsets plus a linear extrapolation is the standard workaround and is stated as
such in the report.

The SPHinXsys ``VelocityGradient`` tensor is stored as G[a][b] = du_a/dx_b, but
only the symmetric combination G + G^T is used here, so the result does not
depend on that convention.
"""

from __future__ import annotations

import numpy as np

from io_utils import Case, trapezoid


def _traction(grad: np.ndarray, normal: np.ndarray, mu: float) -> np.ndarray:
    """mu (G + G^T) . n for stacked tensors.

    ``grad`` is (..., 2, 2) and ``normal`` is (..., 2); returns (..., 2).
    """
    sym = grad + np.swapaxes(grad, -1, -2)
    return mu * np.einsum("...ij,...j->...i", sym, normal)


def _tangential(vec: np.ndarray, normal: np.ndarray) -> np.ndarray:
    """Remove the normal component."""
    n_dot = np.einsum("...i,...i->...", vec, normal)
    return vec - n_dot[..., None] * normal


def wall_shear_stress(case: Case):
    """Extrapolated wall shear stress vector for every wall probe station.

    Returns a dict keyed by wall name, each with

    ``s``      (n_stations,)  arclength along the wall
    ``t``      (n_rows,)      time
    ``tau``    (n_rows, n_stations, 2)  WSS vector at the wall
    ``valid``  (n_stations,)  False where a probe had no fluid support
    """
    d = case.wall_probes()
    mu = case["mu_f"]

    offsets = d["offsets"]
    if len(offsets) < 2:
        raise ValueError("two probe offsets are required to extrapolate to the wall")
    d1, d2 = float(offsets[0]), float(offsets[1])

    # A probe with no fluid neighbours gets a Shepard-normalised value of exactly
    # zero in every component at every time. Real probes never do.
    support = ~np.all(d["grad"] == 0.0, axis=(0, 2, 3))

    out = {}
    for wi, wname in enumerate(d["wall_names"]):
        sel0 = (d["wall_id"] == wi) & (d["offset_id"] == 0)
        sel1 = (d["wall_id"] == wi) & (d["offset_id"] == 1)
        idx0 = np.flatnonzero(sel0)
        idx1 = np.flatnonzero(sel1)
        if len(idx0) == 0 or len(idx0) != len(idx1):
            continue

        s0, s1 = d["s"][idx0], d["s"][idx1]
        if not np.allclose(s0, s1):
            raise AssertionError(f"{wname}: probe stations differ between offsets")

        n = d["normal"][idx0]  # identical for both offsets by construction
        tau0 = _tangential(_traction(d["grad"][:, idx0], n, mu), n)
        tau1 = _tangential(_traction(d["grad"][:, idx1], n, mu), n)

        # linear extrapolation to zero offset
        w = d1 / (d2 - d1)
        tau_wall = tau0 + (tau0 - tau1) * w

        out[wname] = {
            "s": s0,
            "t": d["t"],
            "tau": tau_wall,
            "normal": n,
            "points": d["points"][idx0],
            "valid": support[idx0] & support[idx1],
        }
    return out


def wall_metrics(case: Case):
    """TAWSS, OSI and RRT along every wall, averaged over the analysis window.

    Integrals are taken with the trapezoidal rule on the actual (non-uniform)
    sample times rather than assuming uniform spacing, because the solver uses an
    adaptive time step.
    """
    wss = wall_shear_stress(case)
    out = {}
    for wname, w in wss.items():
        t = w["t"]
        keep = case.in_window(t)
        if keep.sum() < 3:
            keep = np.ones_like(t, dtype=bool)
        tt = t[keep]
        tau = w["tau"][keep]  # (n_rows, n_stations, 2)
        span = tt[-1] - tt[0]
        if span <= 0:
            continue

        mag = np.linalg.norm(tau, axis=-1)             # (n_rows, n_stations)
        int_mag = trapezoid(mag, tt, axis=0)         # (n_stations,)
        int_vec = trapezoid(tau, tt, axis=0)         # (n_stations, 2)

        tawss = int_mag / span
        with np.errstate(divide="ignore", invalid="ignore"):
            osi = 0.5 * (1.0 - np.linalg.norm(int_vec, axis=-1) / int_mag)
            rrt = 1.0 / ((1.0 - 2.0 * osi) * tawss)
        osi = np.where(int_mag > 0, osi, np.nan)
        rrt = np.where(np.isfinite(rrt), rrt, np.nan)

        invalid = ~w["valid"]
        tawss = np.where(invalid, np.nan, tawss)
        osi = np.where(invalid, np.nan, osi)
        rrt = np.where(invalid, np.nan, rrt)

        # fraction of the cycle during which the shear points backwards relative
        # to its own time-mean direction: a direct measure of flow reversal
        mean_dir = int_vec / np.maximum(
            np.linalg.norm(int_vec, axis=-1, keepdims=True), 1e-30
        )
        proj = np.einsum("tsi,si->ts", tau, mean_dir)
        reversal = np.mean(proj < 0.0, axis=0)
        reversal = np.where(invalid, np.nan, reversal)

        out[wname] = {
            "s": w["s"],
            "points": w["points"],
            "tawss": tawss,
            "osi": osi,
            "rrt": rrt,
            "reversal_fraction": reversal,
            "valid": w["valid"],
        }
    return out


def summary(cases: dict[str, Case]):
    """One row per case: the wall metrics condensed to scalars for the study table.

    ``peak_tawss`` locates the most strongly sheared point and ``max_osi`` the most
    oscillatory one; ``apex_*`` restricts attention to the divider wall opposite the
    inlet, which carries the stagnation point where the flow splits.
    """
    rows = []
    for name in sorted(cases):
        c = cases[name]
        try:
            m = wall_metrics(c)
        except (FileNotFoundError, ValueError, AssertionError) as exc:
            print(f"  {name}: WSS skipped ({exc})")
            continue
        if not m:
            continue

        all_tawss = np.concatenate([v["tawss"] for v in m.values()])
        all_osi = np.concatenate([v["osi"] for v in m.values()])
        all_rrt = np.concatenate([v["rrt"] for v in m.values()])

        row = {
            "case": name,
            "Re": c["Re"],
            "alpha": c["alpha"],
            "A": c["A"],
            "dp": c["dp"],
            "branch_ratio": c["branch_ratio"],
            "stenosis": c["stenosis"],
            "peak_tawss": float(np.nanmax(all_tawss)),
            "mean_tawss": float(np.nanmean(all_tawss)),
            "max_osi": float(np.nanmax(all_osi)),
            "mean_osi": float(np.nanmean(all_osi)),
            "max_rrt": float(np.nanmax(all_rrt)),
        }
        if "divider" in m:
            dv = m["divider"]
            row["apex_peak_tawss"] = float(np.nanmax(dv["tawss"]))
            row["apex_min_tawss"] = float(np.nanmin(dv["tawss"]))
            row["apex_max_osi"] = float(np.nanmax(dv["osi"]))
            # the stagnation point is where the shear is weakest on the divider
            i = int(np.nanargmin(dv["tawss"]))
            row["stagnation_s"] = float(dv["s"][i])
        rows.append(row)
    return rows


def nondimensional_tawss(case: Case, tawss: np.ndarray) -> np.ndarray:
    """Scale WSS by the reference viscous stress mu U / R so cases are comparable."""
    ref = case["mu_f"] * case["U_f"] / case["R"]
    return tawss / ref
