"""Two report figures for Section 6 (Results), replacing the dimensional factorial
maps and the baseline-only OSI plot.

  figures/results/factorial_normalised.{pdf,png}
      Re x alpha maps of peak TAWSS*, mean TAWSS* (both / (mu U_f / R)) and max OSI.
      The dimensional maps fall with Re only because mu = rho U D / Re falls.
  figures/results/wall_maps.{pdf,png}
      Wall probes coloured by (a) TAWSS* for the baseline, (b) OSI at Re 200 alpha 5
      (separation-driven reversal), (c) OSI at Re 100 alpha 10 (Womersley reversal).

Usage (from bin/):  python3 analysis/results_clean.py
"""
import sys
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.filterwarnings("ignore")
from io_utils import Case  # noqa: E402
import wss  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "figures" / "results"

plt.rcParams.update({
    "font.size": 9, "savefig.bbox": "tight", "axes.spines.top": False,
    "axes.spines.right": False, "legend.frameon": False,
})

RES = (50, 100, 200)
ALPHAS = (2, 5, 10)


def load(name):
    return Case.load(ROOT / f"output_{name}")


def save(fig, name):
    fig.savefig(OUT / f"{name}.pdf")
    fig.savefig(OUT / f"{name}.png", dpi=200)
    plt.close(fig)
    print("wrote", OUT / f"{name}.{{pdf,png}}")


# ---------------------------------------------------------------- factorial
peak = np.zeros((3, 3))
mean = np.zeros((3, 3))
osi = np.zeros((3, 3))
for j, re in enumerate(RES):
    for i, al in enumerate(ALPHAS):
        c = load(f"A_Re{re}_al{al}")
        m = wss.wall_metrics(c)
        taw = np.concatenate([wss.nondimensional_tawss(c, w["tawss"]) for w in m.values()])
        o = np.concatenate([w["osi"] for w in m.values()])
        peak[i, j] = np.nanmax(taw)
        mean[i, j] = np.nanmean(taw)
        osi[i, j] = np.nanmax(o)

fig, axes = plt.subplots(1, 3, figsize=(10.0, 3.1))
panels = (
    (peak, r"(a) peak TAWSS / $(\mu U_f/R)$", "viridis", "{:.1f}", None),
    (mean, r"(b) wall-mean TAWSS / $(\mu U_f/R)$", "viridis", "{:.2f}", None),
    (osi, "(c) maximum OSI", "magma", "{:.2f}", (0, 0.5)),
)
for ax, (data, title, cmap, fmt, lim) in zip(axes, panels):
    vmin, vmax = lim if lim else (data.min(), data.max())
    im = ax.imshow(data, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    for i in range(3):
        for j in range(3):
            v = data[i, j]
            light = (v - vmin) / (vmax - vmin + 1e-30) < 0.6
            ax.text(j, i, fmt.format(v), ha="center", va="center", fontsize=8,
                    color="white" if light else "black")
    ax.set_xticks(range(3), [str(r) for r in RES])
    ax.set_yticks(range(3), [str(a) for a in ALPHAS])
    ax.set_xlabel("Reynolds number $Re$")
    ax.set_title(title, fontsize=9)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
axes[0].set_ylabel(r"Womersley number $\alpha$")
fig.tight_layout()
save(fig, "factorial_normalised")

# ---------------------------------------------------------------- wall maps
maps = (
    ("A_Re100_al5", "tawss", r"(a) TAWSS / $(\mu U_f/R)$" "\n" r"$Re$ 100, $\alpha$ 5", "viridis", None),
    ("A_Re200_al5", "osi", "(b) OSI\n" r"$Re$ 200, $\alpha$ 5", "magma", (0, 0.5)),
    ("A_Re100_al10", "osi", "(c) OSI\n" r"$Re$ 100, $\alpha$ 10", "magma", (0, 0.5)),
)
fig, axes = plt.subplots(1, 3, figsize=(10.0, 4.4), sharey=True)
for ax, (name, key, title, cmap, lim) in zip(axes, maps):
    c = load(name)
    m = wss.wall_metrics(c)
    xs, ys, vs = [], [], []
    for w in m.values():
        v = w[key]
        if key == "tawss":
            v = wss.nondimensional_tawss(c, v)
        good = np.isfinite(v)
        xs.append(w["points"][good, 0])
        ys.append(w["points"][good, 1])
        vs.append(v[good])
    kw = dict(vmin=lim[0], vmax=lim[1]) if lim else {}
    sc = ax.scatter(np.concatenate(xs), np.concatenate(ys), c=np.concatenate(vs),
                    s=12, cmap=cmap, edgecolors="0.6", linewidths=0.2, **kw)
    ax.set_aspect("equal")
    ax.set_xlabel("$x$")
    ax.set_title(title, fontsize=9)
    ax.set_facecolor("0.93")
    fig.colorbar(sc, ax=ax, fraction=0.05, pad=0.03)
axes[0].set_ylabel("$y$")
fig.tight_layout()
save(fig, "wall_maps")
