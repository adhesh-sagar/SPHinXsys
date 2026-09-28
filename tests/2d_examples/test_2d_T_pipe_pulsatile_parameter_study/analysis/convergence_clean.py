"""Clean two-panel convergence figure for report Section 4.

(a) inlet flow rate: relative error vs the Richardson-extrapolated value, log-log,
    with a 2nd-order reference slope.
(b) wall shear: raw values vs dp, showing the non-monotone behaviour.
All numbers are read from results/convergence_gci.csv.
"""
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

root = Path(sys.argv[1])
rows = {r["functional"]: r for r in csv.DictReader(open(root / "results/convergence_gci.csv"))}

plt.rcParams.update({
    "font.size": 9, "axes.grid": True, "grid.alpha": 0.25,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "lines.linewidth": 1.4, "savefig.bbox": "tight",
})
cmap = plt.get_cmap("viridis")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.0, 3.3))

# (a) flow rate
q = rows["Q_mean"]
dp = np.array(q["dp"].split(), float)
val = np.array(q["values"].split(), float)
ext = float(q["extrapolated"])
order = float(q["observed_order"])
gci = float(q["gci_percent"])
err = np.abs(val - ext) / ext
ax1.loglog(dp, err, "o-", color=cmap(0.15), label=f"inlet flow rate $Q$ (observed order {order:.2f})")
ref = err[0] * (dp / dp[0]) ** 2
ax1.loglog(dp, ref, "k--", lw=1.0, label="2nd-order slope")
ax1.set_xlabel("particle spacing $dp$")
ax1.set_ylabel("relative error vs extrapolated $Q$")
ax1.set_title("(a) inlet flow rate", fontsize=9)
ax1.set_xticks(dp)
ax1.set_xticklabels([f"{d:g}" for d in dp])
ax1.minorticks_off()
ax1.legend(fontsize=8, loc="upper left")

# (b) wall shear
for key, lab, c in (("peak_tawss", "peak TAWSS", 0.35), ("tawss_p95", "TAWSS 95th percentile", 0.7)):
    r = rows[key]
    ax2.plot(np.array(r["dp"].split(), float), np.array(r["values"].split(), float),
             "o-", color=cmap(c), label=lab)
ax2.set_xlabel("particle spacing $dp$")
ax2.set_ylabel("wall shear stress")
ax2.set_title("(b) wall shear stress", fontsize=9)
ax2.set_xticks(dp)
ax2.set_ylim(bottom=0)
ax2.legend(fontsize=8, loc="lower center")

fig.tight_layout()
out = root / "figures/verification"
fig.savefig(out / "convergence_clean.pdf")
fig.savefig(out / "convergence_clean.png", dpi=200)
print("wrote", out / "convergence_clean.{pdf,png}")
