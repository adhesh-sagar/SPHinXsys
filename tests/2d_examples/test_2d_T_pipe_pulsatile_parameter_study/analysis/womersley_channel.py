"""Analysis of the dedicated Womersley channel benchmark (womersley_channel_validation.cpp).

WHY THIS IS A SEPARATE CASE FROM THE T-PIPE

The T-pipe cannot validate the Womersley solution and this module must not pretend
otherwise. In the T-pipe the inlet velocity profile is PRESCRIBED as a parabola by
the inflow buffer, and the inlet channel is only 3.5 long while the viscous
entrance length is ~0.04 Re D = 12 and the oscillatory development length is
U/omega = 3.0. The flow therefore arrives at the junction still carrying the
imposed profile: measured peakedness there is 1.41, against 1.50 for a pure
parabola and 1.21 for the Womersley prediction at alpha = 5.

The benchmark case removes both problems - a straight channel, periodic in x, with
no velocity boundary condition anywhere, driven only by an oscillatory body force
G cos(omega t). Everything in the interior is genuinely computed, so the comparison
below is a real validation of the scheme's unsteady viscous response.

The solver writes the analytical first harmonic at each probe into case_params.json,
so this module only has to extract the simulated harmonic and compare.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from io_utils import interp_periodic, read_dat


def load(folder: str | Path) -> dict:
    """Load one benchmark run and reduce it to first-harmonic profiles."""
    folder = Path(folder)
    with (folder / "case_params.json").open() as f:
        p = json.load(f)
    if p.get("kind") != "womersley_channel":
        raise ValueError(f"{folder} is not a womersley_channel run")

    y = np.asarray(p["y"], dtype=float)
    R, T = p["R"], p["T"]
    uh_ana = np.asarray([complex(a, b) for a, b in p["u_hat_analytic"]])

    t, d, _ = read_dat(folder / "ProfileObserver_Velocity.dat")
    u = d.reshape(len(t), len(y), 2)[:, :, 0]

    mask = t >= p["analysis_start"]
    tw, uw = t[mask], u[mask]

    grid = np.linspace(0.0, 1.0, 200, endpoint=False)
    cyc = np.floor(tw / T).astype(int)
    rows = []
    for c in np.unique(cyc):
        sel = cyc == c
        if sel.sum() < 20:
            continue
        ph = np.mod(tw[sel] / T, 1.0)
        order = np.argsort(ph)
        rows.append(
            np.array(
                [
                    interp_periodic(grid, ph[order], uw[sel][:, j][order])
                    for j in range(len(y))
                ]
            )
        )
    if not rows:
        raise ValueError(f"{folder}: no complete cycles in the analysis window")
    U = np.mean(rows, axis=0)  # (n_y, n_phase)

    # The forcing is G cos(omega t) starting at t = 0, so the simulation and the
    # analytical amplitude share a phase origin and can be compared directly.
    basis = np.exp(2j * np.pi * grid)
    uh_sim = 2.0 * np.mean((U - U.mean(axis=1, keepdims=True)) * np.conj(basis), axis=1)

    scale = float(np.max(np.abs(uh_ana)))
    return {
        "p": p,
        "alpha": p["alpha"],
        "y_over_R": (y - R) / R,
        "phase": grid,
        "u_phase": U,
        "u_hat_sim": uh_sim,
        "u_hat_ana": uh_ana,
        "L2_rel": float(
            np.linalg.norm(uh_sim - uh_ana) / np.linalg.norm(uh_ana)
        ),
        "Linf_rel": float(np.max(np.abs(uh_sim - uh_ana)) / scale),
        "amp_scale": scale,
        "delta_over_dp": p["delta_over_dp"],
    }


def load_all(root: str | Path = ".") -> dict[float, dict]:
    """Every womersley_channel run under ``root``, keyed by alpha."""
    root = Path(root)
    out = {}
    for folder in sorted(root.glob("output_*")):
        meta = folder / "case_params.json"
        if not meta.is_file():
            continue
        try:
            with meta.open() as f:
                if json.load(f).get("kind") != "womersley_channel":
                    continue
            r = load(folder)
        except (ValueError, KeyError, OSError) as exc:
            print(f"  skipping {folder.name}: {exc}")
            continue
        out[r["alpha"]] = r
    return out


def profiles_at_phases(res: dict, n_phase: int = 8):
    """Reconstruct simulated and analytical profiles at ``n_phase`` cycle phases."""
    phases = np.arange(n_phase) / n_phase
    sim, ana = [], []
    for ph in phases:
        e = np.exp(2j * np.pi * ph)
        sim.append(np.real(res["u_hat_sim"] * e))
        ana.append(np.real(res["u_hat_ana"] * e))
    return phases, np.asarray(sim), np.asarray(ana)


def table(root: str | Path = "."):
    """Error summary across the alpha sweep."""
    rows = []
    for alpha, r in sorted(load_all(root).items()):
        i_mid = int(np.argmin(np.abs(r["y_over_R"])))
        rows.append(
            {
                "case": r["p"]["case_name"],
                "alpha": alpha,
                "Re": r["p"]["Re"],
                "dp": r["p"]["dp"],
                "delta_over_dp": r["delta_over_dp"],
                "n_cycles": r["p"]["n_cycles"],
                "L2_rel_pct": 100.0 * r["L2_rel"],
                "Linf_rel_pct": 100.0 * r["Linf_rel"],
                "centre_amp_sim": float(abs(r["u_hat_sim"][i_mid])),
                "centre_amp_ana": float(abs(r["u_hat_ana"][i_mid])),
                "centre_phase_sim_deg": float(
                    np.degrees(np.angle(r["u_hat_sim"][i_mid]))
                ),
                "centre_phase_ana_deg": float(
                    np.degrees(np.angle(r["u_hat_ana"][i_mid]))
                ),
                # off-centre location of the amplitude peak: the Richardson
                # annular effect, which only appears once alpha is large
                "annular_peak_sim": float(
                    abs(r["y_over_R"][int(np.argmax(np.abs(r["u_hat_sim"])))])
                ),
                "annular_peak_ana": float(
                    abs(r["y_over_R"][int(np.argmax(np.abs(r["u_hat_ana"])))])
                ),
            }
        )
    return rows
