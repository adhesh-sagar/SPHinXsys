"""Analytical Womersley solution for a plane channel, and comparison to the SPH result.

THEORY
------
For fully developed oscillatory flow between parallel plates at y = +-R driven by
a pressure gradient -dp/dx = Re{G e^{i w t}}, the momentum equation

    rho du/dt = -dp/dx + mu d2u/dy2

has the solution

    u(y,t) = Re{ (-i G / (rho w)) [ 1 - cosh(lam y) / cosh(lam R) ] e^{i w t} },
    lam = sqrt(i w / nu) = (1 + i) / delta,   delta = sqrt(2 nu / w),

where delta is the Stokes layer thickness. Writing lam R = (1+i) alpha / sqrt(2)
shows the whole shape depends on the Womersley number alpha = R sqrt(w/nu) alone.

Two limits, both of which the simulation should reproduce:

* alpha << 1 (quasi-steady): cosh expands to 1 + (lam y)^2/2, giving the
  instantaneous PARABOLA u = (G/2mu)(R^2 - y^2) in phase with the gradient.
* alpha >> 1 (inertia dominated): the core moves as a slug, u -> Re{-iG/(rho w)},
  which LAGS the pressure gradient by 90 degrees, and the shear is confined to a
  layer of thickness delta at the wall. The velocity maximum moves off the
  centreline (the Richardson annular effect).

The flow rate per unit depth follows by integrating across the channel:

    Q(t) = Re{ (-i G / (rho w)) * 2R * [1 - tanh(lam R) / (lam R)] e^{i w t} }

so the complex impedance Z = G / Q_hat gives the phase lag between the driving
gradient and the flow, which tends to 90 degrees as alpha -> infinity.

WHAT IS COMPARED, AND AN IMPORTANT CAVEAT
-----------------------------------------
The simulation prescribes the inlet VELOCITY, not the pressure gradient, so the
comparison is made the other way round: the measured flow-rate amplitude is used
to fix G, and the resulting analytical profile shape is compared to the measured
one.

CAVEAT: applied to the T-PIPE, this is NOT a validation of the Womersley solution,
and the report must not present it as one. The T-pipe inlet profile is imposed as a
parabola by the inflow buffer, and the inlet channel is only 3.5 long while

    viscous entrance length      ~ 0.04 Re D = 12
    oscillatory development len  ~ U / omega = 3.0

so the flow reaches the junction still carrying the imposed profile. Measured
peakedness at the inlet station is ~1.41, against 1.50 for a pure parabola and
1.21 for the Womersley prediction at alpha = 5 - i.e. much closer to the imposed
profile than to the theory, exactly as those length scales require.

What this module measures on the T-pipe is therefore how far the imposed profile
has RELAXED toward the Womersley solution by the measurement station, which is a
statement about the geometry, not about the numerical scheme.

The actual validation of the scheme's unsteady viscous response is done by the
separate periodic-channel benchmark - see womersley_channel.py and
womersley_channel_validation.cpp - where no velocity is prescribed anywhere and
the agreement is 1-3 % in L2 with the annular peak location matching exactly.
"""

from __future__ import annotations

import numpy as np

from io_utils import Case, resample_cycles


# ---------------------------------------------------------------------------
# analytical solution
# ---------------------------------------------------------------------------


def lam(omega: float, nu: float) -> complex:
    """The Womersley wave number lambda = sqrt(i omega / nu)."""
    return np.sqrt(1j * omega / nu)


def velocity_hat(y, R: float, omega: float, nu: float, rho: float, G: complex):
    """Complex velocity amplitude u_hat(y) for gradient amplitude G.

    ``y`` is measured from the channel centreline; the physical velocity is
    ``Re{u_hat(y) exp(i omega t)}``.
    """
    y = np.asarray(y, dtype=float)
    l = lam(omega, nu)
    # rho du/dt = G e^{iwt} + mu u''  =>  mu uh'' - i w rho uh = -G,
    # whose particular solution is uh_p = G / (i w rho) = -i G / (w rho).
    # The sign matters: getting it wrong flips the predicted phase by 180 deg.
    return (-1j * G / (rho * omega)) * (1.0 - np.cosh(l * y) / np.cosh(l * R))


def flow_rate_hat(R: float, omega: float, nu: float, rho: float, G: complex) -> complex:
    """Complex flow-rate amplitude per unit depth, integral of u_hat over the channel."""
    l = lam(omega, nu)
    return (-1j * G / (rho * omega)) * 2.0 * R * (1.0 - np.tanh(l * R) / (l * R))


def gradient_for_flow_rate(Q_hat: complex, R, omega, nu, rho) -> complex:
    """Invert ``flow_rate_hat``: the G that produces a given flow-rate amplitude."""
    l = lam(omega, nu)
    shape = (-1j / (rho * omega)) * 2.0 * R * (1.0 - np.tanh(l * R) / (l * R))
    return Q_hat / shape


def flow_phase_lag_deg(alpha: float) -> float:
    """Phase lag (degrees) of flow rate behind the driving pressure gradient.

    Depends on alpha alone, and is the single cleanest signature of the Womersley
    problem. Reference values (used as a self-test below):

        alpha = 0.5 ->  5.7 deg     quasi-steady, flow follows the gradient
        alpha = 2   -> 57.2 deg
        alpha = 5   -> 80.7 deg
        alpha = 10  -> 85.7 deg     inertia dominated, approaching the 90 deg limit
    """
    # Work in units R = nu = rho = 1, so omega = alpha^2 and G = 1.
    R = nu = rho = 1.0
    omega = alpha**2
    q = flow_rate_hat(R, omega, nu, rho, 1.0 + 0j)
    return float(np.degrees(-np.angle(q)))


def centreline_amplitude_ratio(alpha: float) -> float:
    """|u(0)| / |u_mean|: how peaked the profile is.

    3/2 in the quasi-steady (parabolic) limit, tending to 1 as alpha -> infinity
    because the core becomes a slug.
    """
    R = nu = rho = 1.0
    omega = alpha**2
    u0 = velocity_hat(0.0, R, omega, nu, rho, 1.0 + 0j)
    q = flow_rate_hat(R, omega, nu, rho, 1.0 + 0j)
    u_mean = q / (2.0 * R)
    return float(abs(u0) / abs(u_mean))


def annular_peak_offset(alpha: float, n: int = 2001) -> float:
    """|y|/R at which the velocity amplitude peaks (Richardson annular effect).

    Zero for small alpha; grows towards the wall as alpha increases.
    """
    R = nu = rho = 1.0
    omega = alpha**2
    y = np.linspace(0.0, R, n)
    amp = np.abs(velocity_hat(y, R, omega, nu, rho, 1.0 + 0j))
    return float(y[int(np.argmax(amp))] / R)


# ---------------------------------------------------------------------------
# harmonic extraction from the simulation
# ---------------------------------------------------------------------------


def first_harmonic(phase: np.ndarray, y: np.ndarray) -> tuple[float, complex]:
    """Mean and complex first-harmonic amplitude of a phase-folded waveform.

    Returns ``(mean, y_hat)`` such that ``y ~ mean + Re{y_hat exp(2 pi i phase)}``.
    Using a projection rather than an FFT keeps this correct for the non-uniform
    sample counts that come out of an adaptive time step.
    """
    phase = np.asarray(phase, dtype=float)
    y = np.asarray(y, dtype=float)
    mean = float(np.mean(y))
    basis = np.exp(2j * np.pi * phase)
    # y_hat = 2 * <y, e^{i phi}> for a real signal
    y_hat = 2.0 * np.mean((y - mean) * np.conj(basis))
    return mean, y_hat


def measured_flow_harmonic(case: Case, section: str = "inlet"):
    """Cycle-averaged mean and first harmonic of the flow rate at one section."""
    t, Q = case.flow_rate(section)
    phase, Y = resample_cycles(t, Q, case)
    if Y.shape[0] == 0:
        return np.nan, np.nan + 0j
    return first_harmonic(phase, Y.mean(axis=0))


def measured_profile_harmonic(case: Case, section: str = "inlet"):
    """Per-point mean and first harmonic of the velocity profile.

    Returns ``(coord, mean, u_hat)`` with one entry per profile point.
    """
    t, coord, u_along, _ = case.profile(section)
    means, hats = [], []
    for j in range(u_along.shape[1]):
        phase, Y = resample_cycles(t, u_along[:, j], case)
        if Y.shape[0] == 0:
            means.append(np.nan)
            hats.append(np.nan + 0j)
            continue
        m, h = first_harmonic(phase, Y.mean(axis=0))
        means.append(m)
        hats.append(h)
    return coord, np.asarray(means), np.asarray(hats)


# ---------------------------------------------------------------------------
# validation
# ---------------------------------------------------------------------------


def compare_inlet_profile(case: Case, n_phase: int = 8):
    """Simulated vs analytical inlet profile at ``n_phase`` phases of the cycle.

    Returns a dict with the centreline-relative coordinate, the simulated and
    analytical profiles at each phase, and the error norms.

    The analytical curve has ONE free parameter (the gradient amplitude G, fixed
    from the measured flow-rate harmonic); the shape and the phase distribution
    across the channel are predictions.
    """
    R = case["R"]
    nu = case["nu"]
    rho = case["rho0_f"]
    omega = case["omega"]

    coord, u_mean_sim, u_hat_sim = measured_profile_harmonic(case, "inlet")
    y = coord - R  # profile points are given in global y in [0, DH]

    Q_mean, Q_hat = measured_flow_harmonic(case, "inlet")
    G = gradient_for_flow_rate(Q_hat, R, omega, nu, rho)

    u_hat_ana = velocity_hat(y, R, omega, nu, rho, G)

    # The steady part is a plain parabola carrying the measured mean flow rate.
    U_mean = Q_mean / (2.0 * R)
    u_mean_ana = 1.5 * U_mean * np.maximum(0.0, 1.0 - (y / R) ** 2)

    phases = np.arange(n_phase) / n_phase
    sim, ana = [], []
    for ph in phases:
        e = np.exp(2j * np.pi * ph)
        sim.append(u_mean_sim + np.real(u_hat_sim * e))
        ana.append(u_mean_ana + np.real(u_hat_ana * e))
    sim = np.asarray(sim)
    ana = np.asarray(ana)

    scale = np.max(np.abs(ana)) if np.max(np.abs(ana)) > 0 else 1.0
    return {
        "case": case.name,
        "alpha": case["alpha"],
        "Re": case["Re"],
        "y_over_R": y / R,
        "phases": phases,
        "u_sim": sim,
        "u_ana": ana,
        "u_mean_sim": u_mean_sim,
        "u_mean_ana": u_mean_ana,
        "u_hat_sim": u_hat_sim,
        "u_hat_ana": u_hat_ana,
        "G": G,
        "L2_rel": float(np.linalg.norm(sim - ana) / np.linalg.norm(ana)),
        "Linf_rel": float(np.max(np.abs(sim - ana)) / scale),
        "L2_steady_rel": float(
            np.linalg.norm(u_mean_sim - u_mean_ana)
            / max(np.linalg.norm(u_mean_ana), 1e-30)
        ),
    }


def phase_lag_table(cases: dict[str, Case]):
    """Measured vs analytical phase lag and peakedness for every case.

    The measured lag is taken between the inlet pressure gradient (approximated by
    the inlet-to-upper-branch pressure drop, which is proportional to it for a
    fixed geometry) and the inlet flow rate. Comparing this against
    ``flow_phase_lag_deg(alpha)`` is the single most informative validation figure
    in the study.
    """
    rows = []
    for name in sorted(cases):
        c = cases[name]
        if c["A"] == 0.0:
            # no pulsation: the first harmonic is numerical noise and its phase
            # carries no information, so the case is excluded rather than plotted
            continue
        try:
            t, Q = c.flow_rate("inlet")
            phQ, YQ = resample_cycles(t, Q, c)
            tp, dp = c.pressure_drop("upper")
            phP, YP = resample_cycles(tp, dp, c)
            if YQ.shape[0] == 0 or YP.shape[0] == 0:
                continue
            _, Q_hat = first_harmonic(phQ, YQ.mean(axis=0))
            _, P_hat = first_harmonic(phP, YP.mean(axis=0))
        except (FileNotFoundError, ValueError, AssertionError) as exc:
            print(f"  {name}: phase lag skipped ({exc})")
            continue

        lag_measured = np.degrees(np.angle(P_hat) - np.angle(Q_hat))
        lag_measured = (lag_measured + 180.0) % 360.0 - 180.0

        coord, u_mean_sim, u_hat_sim = measured_profile_harmonic(c, "inlet")
        R = c["R"]
        y = coord - R
        i_mid = int(np.argmin(np.abs(y)))
        u_mean_bulk = Q.mean() / (2.0 * R) if len(Q) else np.nan

        rows.append(
            {
                "case": name,
                "Re": c["Re"],
                "alpha": c["alpha"],
                "A": c["A"],
                "dp": c["dp"],
                "lag_measured_deg": float(lag_measured),
                "lag_analytic_deg": flow_phase_lag_deg(c["alpha"]),
                "peakedness_measured": float(
                    abs(u_hat_sim[i_mid]) / max(abs(Q_hat) / (2.0 * R), 1e-30)
                ),
                "peakedness_analytic": centreline_amplitude_ratio(c["alpha"]),
                "annular_offset_analytic": annular_peak_offset(c["alpha"]),
                "annular_offset_measured": float(
                    abs(y[int(np.argmax(np.abs(u_hat_sim)))]) / R
                ),
                "u_mean_bulk": float(u_mean_bulk),
            }
        )
    return rows
