"""Loaders for the pulsatile T-pipe study output.

Every run writes one ``output_<case>/`` folder containing

* ``case_params.json``  - all derived parameters plus the exact observer point tables
* ``*_<Quantity>.dat``  - whitespace-separated time series, one header line

The .dat header format differs slightly between the two SPHinXsys recorder types
(``ReducedQuantityRecording`` quotes ``run_time``, ``ObservedQuantityRecording``
does not), so the reader here is deliberately tolerant: it splits the header on
whitespace, strips quotes, and treats the first column as time.

Nothing in this package re-derives a physical parameter: everything comes from
case_params.json, so the Python side can never disagree with the solver.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# np.trapz was renamed np.trapezoid in NumPy 2.0; support both.
trapezoid = getattr(np, "trapezoid", None) or np.trapz

# ---------------------------------------------------------------------------
# low-level .dat reading
# ---------------------------------------------------------------------------

_QUOTE = re.compile(r'^"|"$')


def read_dat(path: str | Path) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Read a SPHinXsys recorder .dat file.

    Returns ``(t, data, names)`` where ``t`` has shape (n_rows,), ``data`` has
    shape (n_rows, n_cols - 1) and ``names`` labels the data columns.
    """
    path = Path(path)
    with path.open() as f:
        header = f.readline()
    names = [_QUOTE.sub("", tok) for tok in header.split()]

    raw = np.loadtxt(path, skiprows=1, ndmin=2)
    if raw.size == 0:
        raise ValueError(f"{path} contains no data rows")
    return raw[:, 0], raw[:, 1:], names[1:]


def read_observer(path: str | Path, n_points: int, comps: int):
    """Read an observer .dat and reshape to (n_rows, n_points, comps).

    ``comps`` is 1 for Real, 2 for Vecd, 4 for Matd (row-major 2x2 in 2D).
    """
    t, data, _ = read_dat(path)
    expected = n_points * comps
    if data.shape[1] != expected:
        raise ValueError(
            f"{path}: expected {expected} data columns "
            f"({n_points} points x {comps} components), found {data.shape[1]}"
        )
    return t, data.reshape(len(t), n_points, comps)


# ---------------------------------------------------------------------------
# one case
# ---------------------------------------------------------------------------


@dataclass
class Case:
    """One completed run: parameters plus lazily-loaded time series."""

    folder: Path
    p: dict

    # -- construction -------------------------------------------------------
    @classmethod
    def load(cls, folder: str | Path) -> "Case":
        folder = Path(folder)
        meta = folder / "case_params.json"
        if not meta.exists():
            raise FileNotFoundError(
                f"{meta} missing - the run probably did not finish "
                "(case_params.json is written last, on purpose)"
            )
        with meta.open() as f:
            p = json.load(f)
        return cls(folder=folder, p=p)

    # -- convenience accessors ---------------------------------------------
    @property
    def name(self) -> str:
        return self.p["case_name"]

    def __getitem__(self, key):
        return self.p[key]

    def get(self, key, default=None):
        """Tolerant accessor, for fields added in later schema versions."""
        return self.p.get(key, default)

    def phase(self, t: np.ndarray) -> np.ndarray:
        """Cycle phase in [0, 1) measured from the start of the pulsation.

        The inflow condition subtracts t_ramp before applying sin(omega t), so
        phase 0 is always "mean flow, rising", identically in every run.
        """
        return np.mod((t - self.p["t_ramp"]) / self.p["T"], 1.0)

    def cycle_index(self, t: np.ndarray) -> np.ndarray:
        """0-based pulsation cycle number; negative during the ramp."""
        return np.floor((t - self.p["t_ramp"]) / self.p["T"]).astype(int)

    def in_window(self, t: np.ndarray) -> np.ndarray:
        """Mask selecting the analysis window (the trailing cycles)."""
        return t >= self.p["analysis_start"]

    # -- global scalar histories -------------------------------------------
    def scalar(self, stem: str) -> tuple[np.ndarray, np.ndarray]:
        """A single-column body diagnostic, e.g. 'TotalKineticEnergy'."""
        t, d, _ = read_dat(self.folder / f"WaterBody_{stem}.dat")
        return t, d[:, 0]

    def kinetic_energy(self):
        return self.scalar("TotalKineticEnergy")

    def max_speed(self):
        return self.scalar("MaximumSpeed")

    def max_density_deviation(self):
        return self.scalar("MaxDensityDeviation")

    # -- volume flow rate ---------------------------------------------------
    def flow_rate(self, section: str) -> tuple[np.ndarray, np.ndarray]:
        """Volume flow rate through one cross-section, outflow positive.

        The solver reduces two sums over a slab of particles:
        ``[sum(sign * u_n * Vol), sum(Vol)]``. The flow rate is

            Q = duct_width * sum(sign * u_n * Vol) / sum(Vol)

        i.e. the duct width times the volume-weighted mean normal velocity.
        Normalising by the accumulated volume rather than the nominal slab
        thickness removes the lattice-layer-counting bias (a slab a few dp thick
        holds a whole number of particle layers, so the nominal thickness can be
        off by tens of percent).
        """
        sec = self.p["observers"]["sections"]
        i = sec["names"].index(section)
        width = sec["width"][i]

        t, d, _ = read_dat(self.folder / f"WaterBody_FluxAndVolume_{section}.dat")
        flux, vol = d[:, 0], d[:, 1]
        Q = np.full_like(flux, np.nan)
        good = vol > 0.0
        Q[good] = width * flux[good] / vol[good]
        return t, Q

    def all_flow_rates(self):
        """``(t, Q_inlet, Q_upper, Q_lower)`` on the inlet's time base.

        The three recorders are written in the same loop pass, so their time
        bases are already identical; this is asserted rather than interpolated.
        """
        t, qi = self.flow_rate("inlet")
        t2, qu = self.flow_rate("upper")
        t3, ql = self.flow_rate("lower")
        n = min(len(t), len(t2), len(t3))
        if not (np.allclose(t[:n], t2[:n]) and np.allclose(t[:n], t3[:n])):
            raise AssertionError("flow-rate recorders are not on a common time base")
        return t[:n], qi[:n], qu[:n], ql[:n]

    # -- point histories ----------------------------------------------------
    def centreline(self, which: str, quantity: str):
        """Velocity (n,2) or Pressure (n,1) at one centreline point."""
        body = {
            "inlet": "InletCentrelineObserver",
            "upper": "UpperOutletObserver",
            "lower": "LowerOutletObserver",
        }[which]
        comps = 2 if quantity == "Velocity" else 1
        t, a = read_observer(self.folder / f"{body}_{quantity}.dat", 1, comps)
        return t, a[:, 0, :]

    def pressure_drop(self, branch: str):
        """p_inlet - p_branch as a time series, on the inlet time base."""
        t, pi = self.centreline("inlet", "Pressure")
        tb, pb = self.centreline(branch, "Pressure")
        n = min(len(t), len(tb))
        return t[:n], (pi[:n, 0] - pb[:n, 0])

    # -- cross-section profiles --------------------------------------------
    def profile(self, section: str):
        """``(t, coord, u_along, u_across)`` for one cross-section.

        ``coord`` is the transverse coordinate of each profile point (y for the
        inlet section, x for the branch sections) and ``u_along`` is the velocity
        component along the duct axis.
        """
        prof = self.p["observers"]["profile"]
        i = prof["names"].index(section)
        b, n = prof["begin"][i], prof["count"][i]
        pts = np.asarray(prof["points"], dtype=float)

        n_total = len(pts)
        t, a = read_observer(
            self.folder / "ProfileObserver_Velocity.dat", n_total, 2
        )
        sl = slice(b, b + n)

        axis = self.p["observers"]["sections"]["axis"][
            self.p["observers"]["sections"]["names"].index(section)
        ]
        coord = pts[sl, 1 - axis]
        u_along = a[:, sl, axis]
        u_across = a[:, sl, 1 - axis]
        return t, coord, u_along, u_across

    # -- wall probes --------------------------------------------------------
    def wall_probes(self):
        """All near-wall probe data needed for WSS.

        Returns a dict with

        ``t``          (n_rows,)
        ``wall_id``    (n_probes,) index into ``wall_names``
        ``offset_id``  (n_probes,) index into ``offsets``
        ``s``          (n_probes,) arclength along the wall
        ``normal``     (n_probes, 2) inward unit normal at the wall foot point
        ``points``     (n_probes, 2) probe location
        ``vel``        (n_rows, n_probes, 2)
        ``grad``       (n_rows, n_probes, 2, 2) velocity gradient tensor
        ``pressure``   (n_rows, n_probes)
        """
        wp = self.p["observers"]["wall_probes"]
        n = len(wp["wall_id"])

        t, vel = read_observer(
            self.folder / "WallProbeObserver_Velocity.dat", n, 2
        )
        _, grad = read_observer(
            self.folder / "WallProbeObserver_VelocityGradient.dat", n, 4
        )
        _, pres = read_observer(
            self.folder / "WallProbeObserver_Pressure.dat", n, 1
        )

        return dict(
            t=t,
            wall_names=wp["wall_names"],
            offsets=np.asarray(wp["offsets"], dtype=float),
            wall_id=np.asarray(wp["wall_id"], dtype=int),
            offset_id=np.asarray(wp["offset_id"], dtype=int),
            s=np.asarray(wp["arclength"], dtype=float),
            normal=np.column_stack(
                [np.asarray(wp["normal_x"]), np.asarray(wp["normal_y"])]
            ).astype(float),
            points=np.asarray(wp["points"], dtype=float),
            vel=vel,
            grad=grad.reshape(len(t), n, 2, 2),
            pressure=pres[:, :, 0],
        )


# ---------------------------------------------------------------------------
# a whole study
# ---------------------------------------------------------------------------


def load_study(root: str | Path = ".") -> dict[str, Case]:
    """Load every ``output_*`` folder under ``root`` that has case_params.json.

    Folders from runs that crashed (no case_params.json) are skipped with a
    warning rather than aborting the whole analysis.
    """
    root = Path(root)
    cases: dict[str, Case] = {}
    for folder in sorted(root.glob("output_*")):
        if not folder.is_dir():
            continue
        try:
            case = Case.load(folder)
        except FileNotFoundError as exc:
            print(f"  skipping {folder.name}: {exc}")
            continue
        # The Womersley channel benchmark writes its own, much smaller metadata
        # into the same kind of folder. It is analysed by womersley_channel.py and
        # must not be handed to the T-pipe analyses, which would fail on the
        # missing keys.
        if case.get("kind", "t_pipe_study") != "t_pipe_study":
            continue
        cases[case.name] = case
    return cases


def interp_periodic(phase_grid, ph, yy):
    """Interpolate a periodic waveform onto ``phase_grid``.

    The samples of one cycle never bracket the full [0, 1) grid: the first sample
    sits a little after phase 0 and the last a little before 1. Rejecting such
    cycles (or letting np.interp clamp at the ends) is wrong, because the signal
    IS periodic - the value just before phase 1 is the value just before phase 0
    of the next cycle. So the data is wrapped by one sample at each end before
    interpolating.

    This matters most for the coarsely sampled signals (profiles and wall probes
    are written 50 times per cycle): with a naive bracketing test every single
    cycle gets discarded and the whole analysis silently returns NaN.
    """
    ph = np.concatenate(([ph[-1] - 1.0], ph, [ph[0] + 1.0]))
    yy = np.concatenate(([yy[-1]], yy, [yy[0]]))
    return np.interp(phase_grid, ph, yy)


def resample_cycles(t, y, case: Case, n_phase: int = 200, min_samples: int = 8):
    """Fold the analysis window onto a common phase grid.

    Returns ``(phase, Y)`` where ``phase`` has ``n_phase`` points in [0, 1) and
    ``Y`` has shape (n_cycles_in_window, n_phase). Cycles carrying fewer than
    ``min_samples`` samples are dropped as under-resolved; the rest are wrapped
    periodically (see ``interp_periodic``) rather than requiring them to bracket
    the grid.
    """
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = case.in_window(t)
    t, y = t[mask], y[mask]
    if t.size == 0:
        return np.linspace(0.0, 1.0, n_phase, endpoint=False), np.zeros((0, n_phase))

    cyc = case.cycle_index(t)
    phase_grid = np.linspace(0.0, 1.0, n_phase, endpoint=False)
    rows = []
    for c in np.unique(cyc):
        sel = cyc == c
        if sel.sum() < min_samples:
            continue
        ph = case.phase(t[sel])
        order = np.argsort(ph)
        rows.append(interp_periodic(phase_grid, ph[order], y[sel][order]))
    if not rows:
        return phase_grid, np.zeros((0, n_phase))
    return phase_grid, np.vstack(rows)


def cycle_average(t, y, case: Case, n_phase: int = 200):
    """Ensemble-average one waveform over the cycles in the analysis window."""
    phase, Y = resample_cycles(t, y, case, n_phase)
    if Y.shape[0] == 0:
        return phase, np.full(n_phase, np.nan)
    return phase, Y.mean(axis=0)
