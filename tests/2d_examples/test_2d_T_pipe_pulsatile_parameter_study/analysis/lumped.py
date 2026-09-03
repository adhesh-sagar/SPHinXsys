"""Lumped-parameter (RL / Windkessel-style) reduction of the simulated bifurcation.

MODEL
-----
Over one branch the unsteady 1D momentum balance integrates to

    dp(t) = R_h Q(t) + L_h dQ/dt

where R_h is the hydraulic resistance and L_h the inertance (fluid inertia). This
is the resistive-inductive element that lumped circulation models are built from,
and fitting it to the simulation connects the distributed CFD result back to the
lumped description used in the lectures.

Analytical values for a plane channel of length L, width h, per unit depth:

    R_h = 12 mu L / h^3          (Poiseuille resistance)
    L_h = rho L / h              (inertance)

Their ratio defines a natural time scale L_h / R_h = rho h^2 / (12 mu), and the
dimensionless group that says which term dominates is

    omega L_h / R_h = omega rho h^2 / (12 mu) = alpha^2 (h/R)^2 / 12,

so the inertial term becomes important at exactly the same place the Womersley
number becomes large. The fit should therefore show L_h mattering more and more
as alpha grows - which is the lumped-model view of the same physics that the
phase lag in womersley.py measures.

FIT
---
R_h and L_h are obtained by ordinary least squares on the cycle-averaged
waveforms, using the analytic derivative of the first harmonic rather than a
finite difference of the noisy signal (WCSPH pressure carries acoustic noise, and
differentiating it directly would amplify that).
"""

from __future__ import annotations

import numpy as np

from io_utils import Case, resample_cycles
from womersley import first_harmonic


def analytic_R(case: Case, which: str) -> float:
    """12 mu L / h^3 for one outlet branch, per unit depth."""
    from conservation import branch_length

    h = case["DL"] - case["DL1"]
    return 12.0 * case["mu_f"] * branch_length(case, which) / h**3


def analytic_L(case: Case, which: str) -> float:
    """rho L / h for one outlet branch, per unit depth."""
    from conservation import branch_length

    h = case["DL"] - case["DL1"]
    return case["rho0_f"] * branch_length(case, which) / h


def fit(case: Case, branch: str = "upper") -> dict:
    """Least-squares fit of dp = R Q + L dQ/dt over the analysis window.

    Both dp and Q are first reduced to their cycle-averaged mean plus first
    harmonic. That makes dQ/dt exact (i omega times the harmonic) instead of a
    finite difference, and it discards the acoustic noise that WCSPH puts on the
    pressure signal at frequencies far above the pulsation.
    """
    omega = case["omega"]

    t, Q = case.flow_rate(branch)
    phQ, YQ = resample_cycles(t, Q, case)
    tp, dp = case.pressure_drop(branch)
    phP, YP = resample_cycles(tp, dp, case)
    if YQ.shape[0] == 0 or YP.shape[0] == 0:
        raise ValueError("no complete cycles in the analysis window")

    Q_mean, Q_hat = first_harmonic(phQ, YQ.mean(axis=0))
    P_mean, P_hat = first_harmonic(phP, YP.mean(axis=0))

    # Reconstruct on a dense phase grid and fit the two coefficients.
    phase = np.linspace(0.0, 1.0, 400, endpoint=False)
    e = np.exp(2j * np.pi * phase)
    q = Q_mean + np.real(Q_hat * e)
    dqdt = np.real(1j * omega * Q_hat * e)      # exact derivative of the harmonic
    p = P_mean + np.real(P_hat * e)

    A = np.column_stack([q, dqdt])
    coef, *_ = np.linalg.lstsq(A, p, rcond=None)
    R_fit, L_fit = float(coef[0]), float(coef[1])

    resid = p - A @ coef
    denom = np.sum((p - p.mean()) ** 2)
    r2 = float(1.0 - np.sum(resid**2) / denom) if denom > 0 else np.nan

    R_ana = analytic_R(case, branch)
    L_ana = analytic_L(case, branch)

    return {
        "case": case.name,
        "branch": branch,
        "Re": case["Re"],
        "alpha": case["alpha"],
        "A": case["A"],
        "dp": case["dp"],
        "R_fit": R_fit,
        "L_fit": L_fit,
        "R_analytic": R_ana,
        "L_analytic": L_ana,
        "R_ratio": R_fit / R_ana if R_ana else np.nan,
        "L_ratio": L_fit / L_ana if L_ana else np.nan,
        "r_squared": r2,
        # which term dominates: > 1 means inertia-dominated
        "omega_L_over_R": float(omega * L_fit / R_fit) if R_fit else np.nan,
        "omega_L_over_R_analytic": float(omega * L_ana / R_ana) if R_ana else np.nan,
        "phase": phase,
        "dp_measured": p,
        "dp_model": A @ coef,
        "Q": q,
    }


def table(cases: dict[str, Case], branch: str = "upper"):
    """RL fit for every case, as a list of scalar rows."""
    rows = []
    for name in sorted(cases):
        try:
            r = fit(cases[name], branch)
        except (FileNotFoundError, ValueError, AssertionError, np.linalg.LinAlgError) as exc:
            print(f"  {name}: RL fit skipped ({exc})")
            continue
        rows.append({k: v for k, v in r.items() if not isinstance(v, np.ndarray)})
    return rows
