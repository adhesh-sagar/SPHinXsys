"""Mass conservation, flow division, and the resistance-network check.

THREE THINGS ARE TESTED HERE
----------------------------
1. Instantaneous mass conservation  Q1 = Q2 + Q3.
   The closure error is reported as a percentage of the inlet flow rate and must
   shrink with the particle spacing. It is NOT expected to be exactly zero at any
   instant: WCSPH is weakly compressible, so the junction can store volume, and
   the imbalance is physically dV_stored/dt. The cycle-MEAN imbalance, by
   contrast, must vanish, because nothing can be stored over a full period once
   the flow is periodic. Both are reported, and the distinction between them is
   worth a paragraph in the report.

2. Flow division  Q2 / (Q2 + Q3).
   For the symmetric geometry this must be 0.5; any deviation is a direct measure
   of numerical asymmetry and bounds the trustworthiness of the split results.

3. The resistance-network prediction (extension D).
   Treating each outlet branch as a plane-channel resistor of hydraulic
   resistance (per unit depth)

       R_h = 12 mu L / h^3

   with L the branch length and h its width, two branches in parallel fed from a
   common junction divide the flow inversely to their resistances:

       Q_up / Q_lo = R_lo / R_up = (L_lo / L_up) * (h_up / h_lo)^3.

   Here both branches have the same width h = DL - DL1 and only the LENGTH is
   varied through --branch_ratio b, so the prediction collapses to

       Q_up / Q_lo = b,   i.e.  Q_up / (Q_up + Q_lo) = b / (1 + b).

   Comparing the measured split against this is a quantitative validation of the
   simulation against a lumped-parameter model - exactly the kind of link between
   distributed and lumped descriptions that the lecture course is about.

   The prediction is a QUASI-STEADY one: it neglects the inertance of each branch,
   so it is expected to hold at low alpha and to degrade as alpha grows. That
   deviation is a result, not a failure, and lumped.py quantifies the inertial
   part separately.
"""

from __future__ import annotations

import numpy as np

from io_utils import Case, resample_cycles


def branch_length(case: Case, which: str) -> float:
    """Length of an outlet branch measured from the junction to its outlet."""
    DH = case["DH"]
    if which == "upper":
        return DH            # y from DH to 2 DH
    if which == "lower":
        return case["branch_ratio"] * DH   # y from -b DH to 0
    raise ValueError(which)


def hydraulic_resistance(case: Case, which: str) -> float:
    """Plane-channel resistance per unit depth, R_h = 12 mu L / h^3."""
    h = case["DL"] - case["DL1"]
    return 12.0 * case["mu_f"] * branch_length(case, which) / h**3


def predicted_split(case: Case) -> float:
    """Q_upper / (Q_upper + Q_lower) from the parallel-resistance network."""
    r_up = hydraulic_resistance(case, "upper")
    r_lo = hydraulic_resistance(case, "lower")
    # parallel divider: the branch with the LOWER resistance takes more flow
    return r_lo / (r_up + r_lo)


def assess(case: Case) -> dict:
    """Conservation and split report for one case, over the analysis window."""
    t, qi, qu, ql = case.all_flow_rates()

    phase, QI = resample_cycles(t, qi, case)
    _, QU = resample_cycles(t, qu, case)
    _, QL = resample_cycles(t, ql, case)
    if QI.shape[0] == 0:
        raise ValueError("no complete cycles in the analysis window")

    qi_c = QI.mean(axis=0)
    qu_c = QU.mean(axis=0)
    ql_c = QL.mean(axis=0)

    imbalance = qi_c - (qu_c + ql_c)
    scale = np.max(np.abs(qi_c))
    scale = scale if scale > 0 else 1.0

    mean_i, mean_u, mean_l = qi_c.mean(), qu_c.mean(), ql_c.mean()
    split = mean_u / (mean_u + mean_l) if (mean_u + mean_l) != 0 else np.nan

    return {
        "case": case.name,
        "Re": case["Re"],
        "alpha": case["alpha"],
        "A": case["A"],
        "dp": case["dp"],
        "branch_ratio": case["branch_ratio"],
        "stenosis": case["stenosis"],
        "phase": phase,
        "Q_inlet": qi_c,
        "Q_upper": qu_c,
        "Q_lower": ql_c,
        "imbalance": imbalance,
        # instantaneous worst case: includes the compressible storage term
        "closure_err_max_pct": float(100.0 * np.max(np.abs(imbalance)) / scale),
        "closure_err_rms_pct": float(
            100.0 * np.sqrt(np.mean(imbalance**2)) / scale
        ),
        # cycle mean: storage must integrate to zero, so this is the strict test
        "closure_err_mean_pct": float(100.0 * np.abs(imbalance.mean()) / scale),
        "Q_mean_inlet": float(mean_i),
        "Q_mean_upper": float(mean_u),
        "Q_mean_lower": float(mean_l),
        "split_measured": float(split),
        "split_predicted": predicted_split(case),
        "split_error_pct": float(100.0 * (split - predicted_split(case))),
        # theoretical inlet flow rate: mean velocity U_f across the full height
        "Q_theory_inlet": float(case["U_f"] * case["DH"]),
        "Q_inlet_error_pct": float(
            100.0 * (mean_i - case["U_f"] * case["DH"]) / (case["U_f"] * case["DH"])
        ),
    }


def pulsatility_index(case: Case) -> dict:
    """PI = (Q_max - Q_min) / Q_mean at the inlet and both outlets.

    The attenuation of PI from the inlet into the branches is a standard clinical
    descriptor (Lecture 3) and is one of the clearest alpha-dependent effects in
    the study: at high alpha the branches damp the pulsation more strongly.
    """
    t, qi, qu, ql = case.all_flow_rates()
    out = {"case": case.name, "Re": case["Re"], "alpha": case["alpha"], "A": case["A"]}
    for label, q in (("inlet", qi), ("upper", qu), ("lower", ql)):
        _, Q = resample_cycles(t, q, case)
        if Q.shape[0] == 0:
            out[f"PI_{label}"] = np.nan
            continue
        qc = Q.mean(axis=0)
        m = qc.mean()
        out[f"PI_{label}"] = float((qc.max() - qc.min()) / m) if m != 0 else np.nan
    if np.isfinite(out.get("PI_inlet", np.nan)) and out["PI_inlet"] != 0:
        out["PI_attenuation_upper"] = out["PI_upper"] / out["PI_inlet"]
        out["PI_attenuation_lower"] = out["PI_lower"] / out["PI_inlet"]
    return out


def table(cases: dict[str, Case]):
    """Conservation + split + PI table for a whole study."""
    rows = []
    for name in sorted(cases):
        c = cases[name]
        try:
            a = assess(c)
            pi = pulsatility_index(c)
        except (FileNotFoundError, ValueError, AssertionError) as exc:
            print(f"  {name}: conservation skipped ({exc})")
            continue
        row = {k: v for k, v in a.items() if not isinstance(v, np.ndarray)}
        for k, v in pi.items():
            if k.startswith("PI"):
                row[k] = v
        rows.append(row)
    return rows
