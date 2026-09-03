#!/usr/bin/env python3
"""Generate every report figure and the master results table.

Usage
-----
    python3 figures.py --root . --out figures

``--root`` is the directory holding the ``output_<case>/`` folders (i.e. where
run_study.sh was executed). Figures are written as PDF (vector, for LaTeX) and
every number that appears in a figure is also written to ``results/*.csv`` so the
tables in the report are traceable to the run that produced them.

Figures are grouped by report section so that they map onto the outline:

  verification/   convergence, periodicity, WCSPH validity
  validation/     Poiseuille, Womersley profiles, phase lag, mass conservation
  results/        Re x alpha maps, amplitude sweep, PI, WSS/OSI, RL fit
  extensions/     branch ratio, stenosis
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import conservation
import convergence
import lumped
import periodicity
import womersley
import womersley_channel
import wss as wss_mod
from io_utils import Case, load_study, resample_cycles

# ---------------------------------------------------------------------------
# plotting style: readable when dropped into a two-column report at 50 % width
# ---------------------------------------------------------------------------
plt.rcParams.update(
    {
        "figure.figsize": (5.2, 3.6),
        "figure.dpi": 130,
        "savefig.bbox": "tight",
        "font.size": 9,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "lines.linewidth": 1.4,
    }
)

CMAP = plt.get_cmap("viridis")


def _save(fig, out: Path, name: str):
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{name}.pdf"
    fig.savefig(path)
    fig.savefig(out / f"{name}.png", dpi=160)
    plt.close(fig)
    print(f"    {path.relative_to(path.parents[2]) if len(path.parents) > 2 else path}")


def _csv(rows, out: Path, name: str):
    if not rows:
        return
    out.mkdir(parents=True, exist_ok=True)
    keys = list(rows[0].keys())
    for r in rows:  # union of keys, stable order
        for k in r:
            if k not in keys:
                keys.append(k)
    path = out / f"{name}.csv"
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"    {path}")


def _group(cases: dict[str, Case], prefix: str) -> dict[str, Case]:
    return {k: v for k, v in cases.items() if k.startswith(prefix)}


# ===========================================================================
# verification
# ===========================================================================


def fig_periodicity(cases, out: Path, res: Path):
    """eps_N vs cycle number, plus the alpha^2 scaling that predicts it."""
    rows = periodicity.table(cases)
    _csv(rows, res, "periodicity_and_wcsph_checks")
    if not rows:
        return

    # (a) eps_N decay for a low-alpha and a high-alpha case
    picks = []
    for target in (2.0, 5.0, 10.0):
        for name in sorted(cases):
            c = cases[name]
            if c["alpha"] == target and c["Re"] == 100 and c["A"] == 0.5:
                picks.append(name)
                break
    if picks:
        fig, ax = plt.subplots()
        for i, name in enumerate(picks):
            c = cases[name]
            t, qi, _, _ = c.all_flow_rates()
            n, eps = periodicity.epsilon_series(t, qi, c)
            if len(n) == 0:
                continue
            ax.semilogy(
                n, eps, "o-", ms=3,
                color=CMAP(i / max(len(picks) - 1, 1)),
                label=rf"$\alpha$ = {c['alpha']:g}",
            )
            ax.axvline(
                0.4 * c["alpha"] ** 2, color=CMAP(i / max(len(picks) - 1, 1)),
                ls=":", lw=1,
            )
        ax.axhline(0.01, color="k", ls="--", lw=1)
        ax.text(0.02, 0.011, "1 % periodicity criterion", fontsize=7,
                transform=ax.get_yaxis_transform())
        ax.set_xlabel("pulsation cycle $N$")
        ax.set_ylabel(r"$\varepsilon_N$  (relative $L_2$ change in $Q(t)$)")
        ax.set_title("Approach to the periodic steady state\n"
                     r"dotted: predicted $0.4\,\alpha^2$ cycles", fontsize=8)
        ax.legend()
        _save(fig, out, "periodicity_epsilon_decay")

    # (b) cycles-to-periodicity vs alpha, measured against theory
    fig, ax = plt.subplots()
    al = np.array([r["alpha"] for r in rows], dtype=float)
    pred = np.array([r["cycles_predicted"] for r in rows], dtype=float)
    ax.scatter(al, pred, s=18, c="C0", label=r"$0.4\,\alpha^2$ (theory)")
    a_fine = np.linspace(min(al.min(), 1), al.max() * 1.1, 100)
    ax.plot(a_fine, 0.4 * a_fine**2, "C0-", lw=1)
    n_run = np.array([r["n_cycles"] for r in rows], dtype=float)
    ax.scatter(al, n_run, s=18, marker="s", c="C1", label="cycles actually run")
    ax.set_xlabel(r"Womersley number $\alpha$")
    ax.set_ylabel("cycles")
    ax.set_title("Why a fixed cycle count fails", fontsize=9)
    ax.legend()
    _save(fig, out, "periodicity_cycles_vs_alpha")

    # (c) WCSPH validity: Mach and density deviation
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.0))
    ma = np.array([r["mach_max"] for r in rows])
    dr = np.array([r["drho_max"] for r in rows]) * 100.0
    idx = np.arange(len(rows))
    a1.bar(idx, ma, color=np.where(ma <= 0.1, "C0", "C3"))
    a1.axhline(0.1, color="k", ls="--", lw=1)
    a1.set_ylabel("max Mach number")
    a2.bar(idx, dr, color=np.where(dr <= 1.0, "C0", "C3"))
    a2.axhline(1.0, color="k", ls="--", lw=1)
    a2.set_ylabel(r"max $|\rho-\rho_0|/\rho_0$  [%]")
    for a in (a1, a2):
        a.set_xticks(idx)
        a.set_xticklabels([r["case"] for r in rows], rotation=90, fontsize=5)
    fig.suptitle("Weak-compressibility validity checks", fontsize=9)
    _save(fig, out, "wcsph_validity_checks")


def fig_convergence(cases, out: Path, res: Path):
    """Resolution convergence: functionals vs dp, order, GCI, and the exact-error case."""
    rep = convergence.study(cases)
    series = rep.get("series", [])
    if len(series) < 2 or not rep.get("functionals"):
        print("  convergence: not enough resolutions for a series, skipping")
        # the steady-case comparison below needs no series at all, so fall through
        # to it rather than returning
        _fig_poiseuille(cases, out)
        return

    rows = []
    for fname, e in rep["functionals"].items():
        rows.append(
            {
                "functional": fname,
                "cases": " ".join(e["cases"]),
                "dp": " ".join(f"{d:g}" for d in e["dp"]),
                "values": " ".join(f"{v:.6g}" for v in e["values"]),
                "refinement_ratio": e.get("r", np.nan),
                "r_fine_pair": e.get("r_fine_pair", np.nan),
                "r_coarse_pair": e.get("r_coarse_pair", np.nan),
                "r_uniform": e.get("r_uniform", ""),
                "observed_order": e.get("order", np.nan),
                "extrapolated": e.get("extrapolated", np.nan),
                "gci_percent": e.get("gci_percent", np.nan),
                "monotone": e.get("monotone", False),
            }
        )
    _csv(rows, res, "convergence_gci")

    fig, ax = plt.subplots()
    for i, (fname, e) in enumerate(rep["functionals"].items()):
        v = np.asarray(e["values"], dtype=float)
        d = np.asarray(e["dp"], dtype=float)
        good = np.isfinite(v)
        if good.sum() < 2:
            continue
        ref = e.get("extrapolated", np.nan)
        if not np.isfinite(ref):
            ref = v[good][0]
        err = np.abs(v - ref) / max(abs(ref), 1e-30)
        ax.loglog(d[good], err[good], "o-", ms=4,
                  color=CMAP(i / max(len(rep["functionals"]) - 1, 1)),
                  label=f"{fname} (q = {e.get('order', float('nan')):.2f})")
    dd = np.array(rep["dp"], dtype=float)
    if len(dd) >= 2:
        ax.loglog(dd, (dd / dd[0]) ** 2 * 1e-2, "k--", lw=1, label="2nd order")
    ax.set_xlabel("particle spacing $dp$")
    ax.set_ylabel("relative error vs Richardson extrapolation")
    ax.set_title("Particle-resolution convergence", fontsize=9)
    ax.legend(fontsize=7)
    _save(fig, out, "convergence_functionals")

    _fig_poiseuille(cases, out)


def _fig_poiseuille(cases, out: Path):
    """Exact error against the analytic parabola for any steady (A = 0) case."""
    steady = [c for c in cases.values() if c["A"] == 0.0]
    if steady:
        fig, ax = plt.subplots()
        for i, c in enumerate(sorted(steady, key=lambda x: x["dp"])):
            try:
                r = convergence.poiseuille_error(c)
            except (FileNotFoundError, ValueError) as exc:
                print(f"  {c.name}: Poiseuille error skipped ({exc})")
                continue
            col = CMAP(i / max(len(steady) - 1, 1))
            ax.plot(r["u_sim"], r["y_over_R"], "o", ms=3, color=col,
                    label=f"SPH, dp = {r['dp']:g} ($L_2$ = {r['L2_rel']*100:.2f} %)")
            if i == 0:
                ax.plot(r["u_ana"], r["y_over_R"], "k-", lw=1.2,
                        label="analytic parabola")
        ax.set_xlabel("$u / U_f$")
        ax.set_ylabel("$y / R$")
        ax.set_title("Steady Poiseuille validation (A = 0)", fontsize=9)
        ax.legend(fontsize=7)
        _save(fig, out, "poiseuille_profile")


# ===========================================================================
# validation
# ===========================================================================


def fig_womersley(cases, out: Path, res: Path):
    """Profile overlays at 8 phases, plus phase lag and peakedness vs alpha."""
    rows = womersley.phase_lag_table(cases)
    _csv(rows, res, "tpipe_profile_relaxation")

    # (a) profile overlays for the three alpha values at Re = 100
    picks = []
    for target in (2.0, 5.0, 10.0):
        for name in sorted(cases):
            c = cases[name]
            if c["alpha"] == target and c["Re"] == 100 and c["A"] == 0.5:
                picks.append(name)
                break
    if picks:
        fig, axes = plt.subplots(1, len(picks), figsize=(2.6 * len(picks), 3.4),
                                 sharey=True)
        axes = np.atleast_1d(axes)
        for ax, name in zip(axes, picks):
            c = cases[name]
            try:
                r = womersley.compare_inlet_profile(c, n_phase=8)
            except (FileNotFoundError, ValueError) as exc:
                print(f"  {name}: Womersley overlay skipped ({exc})")
                continue
            for k in range(len(r["phases"])):
                col = CMAP(k / len(r["phases"]))
                ax.plot(r["u_ana"][k], r["y_over_R"], "-", color=col, lw=1)
                ax.plot(r["u_sim"][k], r["y_over_R"], "o", color=col, ms=2)
            ax.set_title(rf"$\alpha$ = {c['alpha']:g}" "\n"
                         rf"$L_2$ = {r['L2_rel']*100:.1f} %", fontsize=8)
            ax.set_xlabel("$u / U_f$")
        axes[0].set_ylabel("$y / R$")
        fig.suptitle(
            "T-pipe inlet station: imposed profile vs Womersley prediction\n"
            "NOT a validation - the inlet profile is prescribed and the channel is\n"
            "shorter than the entrance length (see womersley_channel_* figures)",
            fontsize=7)
        _save(fig, out, "tpipe_inlet_profile_relaxation")

    if not rows:
        return

    # (b) phase lag vs alpha against the analytic curve
    fig, ax = plt.subplots()
    a = np.array([r["alpha"] for r in rows], dtype=float)
    meas = np.array([r["lag_measured_deg"] for r in rows], dtype=float)
    a_fine = np.linspace(0.3, max(a.max() * 1.15, 12), 200)
    ana = np.array([womersley.flow_phase_lag_deg(x) for x in a_fine])
    ax.plot(a_fine, ana, "k-", lw=1.3, label="Womersley theory")
    ax.axhline(90, color="grey", ls=":", lw=1)
    ax.text(a_fine[-1], 88, r"$90^\circ$ limit", ha="right", va="top", fontsize=7)
    ax.scatter(a, meas, s=22, c="C3", zorder=3, label="SPH")
    ax.set_xlabel(r"Womersley number $\alpha$")
    ax.set_ylabel(r"phase lag of $Q$ behind $\Delta p$  [deg]")
    ax.set_title("Flow lags the driving pressure gradient as $\\alpha$ grows",
                 fontsize=9)
    ax.legend()
    _save(fig, out, "tpipe_phase_lag_vs_theory")

    # (c) profile peakedness and the Richardson annular effect
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.0))
    pk_m = np.array([r["peakedness_measured"] for r in rows])
    an_m = np.array([r["annular_offset_measured"] for r in rows])
    a_fine = np.linspace(0.3, max(a.max() * 1.15, 12), 120)
    a1.plot(a_fine, [womersley.centreline_amplitude_ratio(x) for x in a_fine],
            "k-", lw=1.3, label="theory")
    a1.scatter(a, pk_m, s=22, c="C3", zorder=3, label="SPH")
    a1.axhline(1.5, color="grey", ls=":", lw=1)
    a1.text(a_fine[-1], 1.51, "parabolic (3/2)", ha="right", fontsize=7)
    a1.set_xlabel(r"$\alpha$")
    a1.set_ylabel(r"$|u(0)| / |u_{\rm mean}|$")
    a1.set_title("Core flattening", fontsize=9)
    a1.legend(fontsize=7)

    a2.plot(a_fine, [womersley.annular_peak_offset(x) for x in a_fine],
            "k-", lw=1.3, label="theory")
    a2.scatter(a, an_m, s=22, c="C3", zorder=3, label="SPH")
    a2.set_xlabel(r"$\alpha$")
    a2.set_ylabel(r"$|y|_{\rm peak} / R$")
    a2.set_title("Richardson annular effect", fontsize=9)
    a2.legend(fontsize=7)
    _save(fig, out, "tpipe_profile_shape_vs_theory")


def fig_womersley_channel(cases, out: Path, res: Path, root: Path):
    """The dedicated Womersley benchmark: profiles at 8 phases against theory.

    This is the figure that actually validates the unsteady viscous scheme. The
    T-pipe figures cannot do it (prescribed inlet profile, inlet channel far
    shorter than the entrance length), so the two are kept visually distinct.
    """
    runs = womersley_channel.load_all(root)
    if not runs:
        print("    no womersley_channel runs found, skipping")
        return
    rows = womersley_channel.table(root)
    _csv(rows, res, "womersley_channel_validation")

    alphas = sorted(runs)
    fig, axes = plt.subplots(1, len(alphas), figsize=(2.7 * len(alphas), 3.5),
                             sharey=True)
    axes = np.atleast_1d(axes)
    for ax, a in zip(axes, alphas):
        r = runs[a]
        phases, sim, ana = womersley_channel.profiles_at_phases(r, 8)
        for k in range(len(phases)):
            col = CMAP(k / len(phases))
            ax.plot(ana[k], r["y_over_R"], "-", color=col, lw=1.1)
            ax.plot(sim[k], r["y_over_R"], "o", color=col, ms=2.2, mfc="none")
        ax.set_title(
            rf"$\alpha$ = {a:g}" "\n"
            rf"$L_2$ = {100*r['L2_rel']:.2f} %, $\delta/dp$ = {r['delta_over_dp']:.1f}",
            fontsize=8,
        )
        ax.set_xlabel("$u / U_{ref}$")
        ax.axhline(0, color="0.8", lw=0.6)
    axes[0].set_ylabel("$y / R$")
    fig.suptitle(
        "Womersley benchmark: oscillatory channel flow\n"
        "lines analytical, open circles SPH (colour = phase)", fontsize=8)
    _save(fig, out, "womersley_channel_profiles")

    # amplitude and phase across the channel, the quantitative version
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.1))
    for i, a in enumerate(alphas):
        r = runs[a]
        col = CMAP(i / max(len(alphas) - 1, 1))
        a1.plot(np.abs(r["u_hat_ana"]), r["y_over_R"], "-", color=col, lw=1.2,
                label=rf"$\alpha$ = {a:g}")
        a1.plot(np.abs(r["u_hat_sim"]), r["y_over_R"], "o", color=col, ms=2.5,
                mfc="none")
        a2.plot(np.degrees(np.angle(r["u_hat_ana"])), r["y_over_R"], "-",
                color=col, lw=1.2)
        a2.plot(np.degrees(np.angle(r["u_hat_sim"])), r["y_over_R"], "o",
                color=col, ms=2.5, mfc="none")
    a1.set_xlabel(r"$|\hat{u}|$"); a1.set_ylabel("$y / R$")
    a1.set_title("Amplitude", fontsize=9); a1.legend(fontsize=7)
    a2.set_xlabel(r"$\arg\,\hat{u}$  [deg]")
    a2.set_title("Phase (the Stokes layer leads the core)", fontsize=9)
    fig.suptitle("Womersley benchmark, first harmonic: lines theory, circles SPH",
                 fontsize=8)
    _save(fig, out, "womersley_channel_harmonic")


def fig_conservation(cases, out: Path, res: Path):
    """Q1 = Q2 + Q3 waveforms, closure error vs dp, and the flow split."""
    rows = conservation.table(cases)
    _csv(rows, res, "mass_conservation_and_split")
    if not rows:
        return

    # (a) waveforms for the baseline case
    base = None
    for name in sorted(cases):
        c = cases[name]
        if c["Re"] == 100 and c["alpha"] == 5 and c["A"] == 0.5 and c["stenosis"] == 0 \
                and c["branch_ratio"] == 1.0:
            base = c
            break
    if base is not None:
        try:
            a = conservation.assess(base)
            fig, ax = plt.subplots()
            ax.plot(a["phase"], a["Q_inlet"], "k-", label="$Q_1$ inlet")
            ax.plot(a["phase"], a["Q_upper"], "-", color="C0", label="$Q_2$ upper")
            ax.plot(a["phase"], a["Q_lower"], "-", color="C1", label="$Q_3$ lower")
            ax.plot(a["phase"], a["Q_upper"] + a["Q_lower"], "--", color="C3",
                    label="$Q_2+Q_3$")
            ax.set_xlabel("phase $t/T$")
            ax.set_ylabel("volume flow rate")
            ax.set_title(
                "Mass conservation over the cycle\n"
                f"rms closure error {a['closure_err_rms_pct']:.2f} %, "
                f"cycle-mean {a['closure_err_mean_pct']:.3f} %",
                fontsize=8,
            )
            ax.legend(fontsize=7, ncol=2)
            _save(fig, out, "mass_conservation_waveforms")
        except (FileNotFoundError, ValueError) as exc:
            print(f"  conservation waveform skipped ({exc})")

    # (b) closure error vs resolution
    by_dp = {}
    for r in rows:
        if r["branch_ratio"] == 1.0 and r["stenosis"] == 0 and r["alpha"] == 5 \
                and r["Re"] == 100 and r["A"] == 0.5:
            by_dp[r["dp"]] = r
    if len(by_dp) >= 2:
        d = np.array(sorted(by_dp))
        e_rms = np.array([by_dp[x]["closure_err_rms_pct"] for x in d])
        e_mean = np.array([by_dp[x]["closure_err_mean_pct"] for x in d])
        fig, ax = plt.subplots()
        ax.loglog(d, e_rms, "o-", ms=4, label="instantaneous (rms)")
        ax.loglog(d, e_mean, "s-", ms=4, label="cycle mean")
        ax.set_xlabel("particle spacing $dp$")
        ax.set_ylabel("closure error  [% of $Q_1$]")
        ax.set_title("Mass-conservation error shrinks with resolution", fontsize=9)
        ax.legend(fontsize=7)
        _save(fig, out, "mass_conservation_convergence")

    # (c) flow split vs branch ratio against the resistance network
    d_rows = [r for r in rows if r["stenosis"] == 0]
    if len({r["branch_ratio"] for r in d_rows}) >= 2:
        br = np.array([r["branch_ratio"] for r in d_rows])
        meas = np.array([r["split_measured"] for r in d_rows])
        order = np.argsort(br)
        b_fine = np.linspace(min(br.min(), 0.4), 1.0, 100)
        fig, ax = plt.subplots()
        ax.plot(b_fine, b_fine / (1.0 + b_fine), "k-", lw=1.3,
                label=r"resistance network  $b/(1+b)$")
        ax.scatter(br[order], meas[order], s=26, c="C3", zorder=3, label="SPH")
        ax.set_xlabel("lower-branch length ratio $b$")
        ax.set_ylabel(r"$Q_2 / (Q_2 + Q_3)$")
        ax.set_title("Flow division follows the lumped resistance model", fontsize=9)
        ax.legend(fontsize=7)
        _save(fig, out, "flow_split_vs_branch_ratio")


# ===========================================================================
# parameter study results
# ===========================================================================


def fig_parameter_maps(cases, out: Path, res: Path):
    """Re x alpha maps of the headline quantities, and the amplitude sweep."""
    wrows = wss_mod.summary(cases)
    _csv(wrows, res, "wall_metrics")
    crows = conservation.table(cases)
    lrows = lumped.table(cases)
    _csv(lrows, res, "lumped_RL_fit")

    fac = [r for r in wrows if r["A"] == 0.5 and r["branch_ratio"] == 1.0
           and r["stenosis"] == 0]
    if fac:
        Re_vals = sorted({r["Re"] for r in fac})
        al_vals = sorted({r["alpha"] for r in fac})
        if len(Re_vals) >= 2 and len(al_vals) >= 2:
            for key, label in (
                ("peak_tawss", "peak TAWSS"),
                ("max_osi", "max OSI"),
            ):
                M = np.full((len(al_vals), len(Re_vals)), np.nan)
                for r in fac:
                    M[al_vals.index(r["alpha"]), Re_vals.index(r["Re"])] = r[key]
                fig, ax = plt.subplots()
                im = ax.imshow(M, origin="lower", aspect="auto", cmap="viridis")
                ax.set_xticks(range(len(Re_vals)))
                ax.set_xticklabels([f"{v:g}" for v in Re_vals])
                ax.set_yticks(range(len(al_vals)))
                ax.set_yticklabels([f"{v:g}" for v in al_vals])
                ax.set_xlabel("Reynolds number $Re$")
                ax.set_ylabel(r"Womersley number $\alpha$")
                ax.set_title(f"{label} over the $Re \\times \\alpha$ factorial",
                             fontsize=9)
                ax.grid(False)
                for i in range(M.shape[0]):
                    for j in range(M.shape[1]):
                        if np.isfinite(M[i, j]):
                            ax.text(j, i, f"{M[i, j]:.3g}", ha="center",
                                    va="center", color="w", fontsize=7)
                fig.colorbar(im, ax=ax)
                _save(fig, out, f"factorial_{key}")

    # pulsatility index attenuation vs alpha
    pi_rows = [r for r in crows if "PI_inlet" in r and np.isfinite(r.get("PI_inlet", np.nan))]
    if pi_rows:
        fig, ax = plt.subplots()
        a = np.array([r["alpha"] for r in pi_rows])
        att = np.array([
            r["PI_upper"] / r["PI_inlet"] if r["PI_inlet"] else np.nan
            for r in pi_rows
        ])
        ax.scatter(a, att, s=24, c=[r["Re"] for r in pi_rows], cmap="viridis")
        ax.axhline(1.0, color="k", ls="--", lw=1)
        ax.set_xlabel(r"$\alpha$")
        ax.set_ylabel(r"PI$_{\rm branch}$ / PI$_{\rm inlet}$")
        ax.set_title("Pulsatility is damped through the bifurcation", fontsize=9)
        _save(fig, out, "pulsatility_attenuation")

    # amplitude sweep: flow reversal
    b = _group(cases, "B_")
    if b:
        fig, ax = plt.subplots()
        for i, name in enumerate(sorted(b, key=lambda n: b[n]["A"])):
            c = b[name]
            try:
                t, q = c.flow_rate("inlet")
                ph, Y = resample_cycles(t, q, c)
                if Y.shape[0] == 0:
                    continue
            except (FileNotFoundError, ValueError) as exc:
                print(f"  {name}: amplitude sweep skipped ({exc})")
                continue
            ax.plot(ph, Y.mean(axis=0), color=CMAP(i / max(len(b) - 1, 1)),
                    label=f"A = {c['A']:g}")
        ax.axhline(0.0, color="k", ls="--", lw=1)
        ax.set_xlabel("phase $t/T$")
        ax.set_ylabel("inlet flow rate")
        ax.set_title("Amplitude sweep: A > 1 drives flow reversal", fontsize=9)
        ax.legend(fontsize=7)
        _save(fig, out, "amplitude_sweep_flow_reversal")


def fig_wall_metrics(cases, out: Path, res: Path):
    """TAWSS / OSI distributions along the walls for the baseline case."""
    base = None
    for name in sorted(cases):
        c = cases[name]
        if c["Re"] == 100 and c["alpha"] == 5 and c["A"] == 0.5 \
                and c["stenosis"] == 0 and c["branch_ratio"] == 1.0:
            base = c
            break
    if base is None:
        return
    try:
        m = wss_mod.wall_metrics(base)
    except (FileNotFoundError, ValueError, AssertionError) as exc:
        print(f"  wall metric figure skipped ({exc})")
        return

    for key, label in (("tawss", "TAWSS"), ("osi", "OSI")):
        fig, ax = plt.subplots(figsize=(6.0, 3.4))
        for i, (wname, w) in enumerate(m.items()):
            y = w[key]
            if key == "tawss":
                y = wss_mod.nondimensional_tawss(base, y)
            ax.plot(w["s"], y, "-", color=CMAP(i / max(len(m) - 1, 1)), label=wname)
        ax.set_xlabel("arclength along wall")
        ax.set_ylabel(
            r"TAWSS / $(\mu U_f / R)$" if key == "tawss" else "OSI"
        )
        if key == "osi":
            ax.set_ylim(0, 0.5)
        ax.set_title(f"{label} distribution, {base.name}", fontsize=9)
        ax.legend(fontsize=7, ncol=3)
        _save(fig, out, f"wall_{key}_distribution")

    # spatial map: colour the probe points by TAWSS so the geometry is visible
    fig, ax = plt.subplots(figsize=(5.4, 4.6))
    xs, ys, vs = [], [], []
    for w in m.values():
        good = np.isfinite(w["tawss"])
        xs.append(w["points"][good, 0])
        ys.append(w["points"][good, 1])
        vs.append(wss_mod.nondimensional_tawss(base, w["tawss"][good]))
    sc = ax.scatter(np.concatenate(xs), np.concatenate(ys),
                    c=np.concatenate(vs), s=10, cmap="viridis")
    ax.set_aspect("equal")
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_title("Where the wall is sheared", fontsize=9)
    ax.grid(False)
    fig.colorbar(sc, ax=ax, label=r"TAWSS / $(\mu U_f/R)$")
    _save(fig, out, "wall_tawss_map")


def fig_extensions(cases, out: Path, res: Path):
    """Stenosis series: peak WSS and the throat jet."""
    e = _group(cases, "E_")
    if not e:
        return
    rows = []
    for name in sorted(e, key=lambda n: e[n]["stenosis"]):
        c = e[name]
        try:
            m = wss_mod.wall_metrics(c)
        except (FileNotFoundError, ValueError, AssertionError) as exc:
            print(f"  {name}: stenosis metrics skipped ({exc})")
            continue
        vals = [np.nanmax(v["tawss"]) for v in m.values() if np.isfinite(v["tawss"]).any()]
        rows.append({
            "case": name,
            "stenosis": c["stenosis"],
            "peak_tawss": float(np.nanmax(vals)) if vals else np.nan,
            "peak_tawss_nd": float(
                wss_mod.nondimensional_tawss(c, np.nanmax(vals))
            ) if vals else np.nan,
        })
    _csv(rows, res, "stenosis_series")
    if len(rows) < 2:
        return
    fig, ax = plt.subplots()
    ax.plot([r["stenosis"] * 100 for r in rows],
            [r["peak_tawss_nd"] for r in rows], "o-", ms=5)
    ax.set_xlabel("area reduction  [%]")
    ax.set_ylabel(r"peak TAWSS / $(\mu U_f/R)$")
    ax.set_title("Stenosis severity drives peak wall shear", fontsize=9)
    _save(fig, out, "stenosis_peak_wss")


# ===========================================================================


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=".", help="directory holding output_* folders")
    ap.add_argument("--out", default="figures", help="where to write figures")
    args = ap.parse_args()

    root = Path(args.root)
    out = Path(args.out)
    res = out.parent / "results" if out.name == "figures" else out / "results"

    print(f"loading study from {root.resolve()}")
    cases = load_study(root)
    if not cases:
        print("no completed cases found (looking for output_*/case_params.json)")
        return 1
    print(f"loaded {len(cases)} case(s): {', '.join(sorted(cases))}\n")

    print("  [womersley channel benchmark]")
    try:
        fig_womersley_channel(cases, out / "validation", res, root)
    except Exception as exc:
        print(f"    !! womersley channel failed: {type(exc).__name__}: {exc}")

    for label, fn, sub in (
        ("verification", fig_periodicity, "verification"),
        ("convergence", fig_convergence, "verification"),
        ("womersley validation", fig_womersley, "validation"),
        ("mass conservation", fig_conservation, "validation"),
        ("parameter maps", fig_parameter_maps, "results"),
        ("wall metrics", fig_wall_metrics, "results"),
        ("extensions", fig_extensions, "extensions"),
    ):
        print(f"  [{label}]")
        try:
            fn(cases, out / sub, res)
        except Exception as exc:  # a bad case must not kill the whole report
            print(f"    !! {label} failed: {type(exc).__name__}: {exc}")

    print(f"\nfigures in {out.resolve()}\nnumbers in {res.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
