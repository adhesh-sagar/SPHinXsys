"""Report figure for Section 7 (Extensions), replacing flow_split_vs_branch_ratio,
whose title claimed agreement with the naive resistance model it actually refutes.

  figures/extensions/extensions_summary.{pdf,png}
      (a) upper-branch flow fraction vs lower-branch length ratio b: SPH, the naive
          Poiseuille network b/(1+b), and the network with a length-independent
          resistance R0 added to each branch (R0/R fitted, mean of the two cases).
      (b) stenosis series: upper-branch flow fraction and its pulsatility index.

Usage (from bin/):  python3 analysis/extensions_clean.py
"""
import csv
import sys
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.filterwarnings("ignore")
from io_utils import Case  # noqa: E402
import conservation  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "figures" / "extensions"

plt.rcParams.update({
    "font.size": 9, "savefig.bbox": "tight", "axes.grid": True, "grid.alpha": 0.25,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
})
cmap = plt.get_cmap("viridis")

rows = {r["case"]: r for r in csv.DictReader(open(ROOT / "results/mass_conservation_and_split.csv"))}
split = lambda name: float(rows[name]["split_measured"])  # noqa: E731

# (a) branch ratio
b_pts = np.array([0.5, 0.7, 1.0])
f_pts = np.array([split("D_br0.5"), split("D_br0.7"), split("A_Re100_al5")])
r0 = [(b - f * (1 + b)) / (2 * f - 1) for b, f in zip(b_pts[:2], f_pts[:2])]
r0_mean = float(np.mean(r0))
b = np.linspace(0.4, 1.0, 200)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.4, 3.4))
ax1.plot(b, b / (1 + b), "--", color="0.35", label=r"Poiseuille network $b/(1+b)$")
ax1.plot(b, (r0_mean + b) / (2 * r0_mean + 1 + b), "-", color=cmap(0.35),
         label=rf"with branch loss $R_0$ = {r0_mean:.1f} $R$")
ax1.plot(b_pts, f_pts, "o", color=cmap(0.8), ms=6, mec="k", mew=0.5, label="SPH")
ax1.set_xlabel("lower-branch length ratio $b$")
ax1.set_ylabel(r"upper-branch flow fraction $Q_2/Q_1$")
ax1.set_title("(a) asymmetric branches", fontsize=9)
ax1.legend(fontsize=8, loc="lower right")

# (b) stenosis
names = ["A_Re100_al5", "E_st0.3", "E_st0.5", "E_st0.7"]
st = np.array([0, 30, 50, 70])
frac = np.array([split(n) for n in names])
pi_up = np.array([conservation.pulsatility_index(Case.load(ROOT / f"output_{n}"))["PI_upper"]
                  for n in names])
ax2.plot(st, frac, "o-", color=cmap(0.15), label="flow fraction $Q_2/Q_1$")
ax2.set_xlabel("stenosis severity [% width reduction]")
ax2.set_ylabel(r"upper-branch flow fraction $Q_2/Q_1$", color=cmap(0.15))
ax2.set_ylim(0, 0.55)
ax2.set_xticks(st)
ax3 = ax2.twinx()
ax3.plot(st, pi_up, "s--", color=cmap(0.7), label="upper-branch PI")
ax3.set_ylabel("upper-branch pulsatility index", color=cmap(0.7))
ax3.set_ylim(0.9, 1.3)
ax3.grid(False)
ax3.spines["right"].set_visible(True)
ax2.set_title("(b) stenosis in the upper branch", fontsize=9)
h1, l1 = ax2.get_legend_handles_labels()
h2, l2 = ax3.get_legend_handles_labels()
ax2.legend(h1 + h2, l1 + l2, fontsize=8, loc="lower left")

fig.tight_layout()
OUT.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT / "extensions_summary.pdf")
fig.savefig(OUT / "extensions_summary.png", dpi=200)
print("wrote", OUT / "extensions_summary.{pdf,png}")
print("R0/R per case:", [round(x, 2) for x in r0], "mean", round(r0_mean, 2))
print("re-predicted:", [round((r0_mean + x) / (2 * r0_mean + 1 + x), 4) for x in b_pts[:2]],
      "measured:", [round(x, 4) for x in f_pts[:2]])
print("stenosis frac", frac.round(4), "PI_up", pi_up.round(3))
