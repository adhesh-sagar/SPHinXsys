# Report guide — what to include, where it comes from, what to say

For a **15–20 page** Biofluid Mechanics report on the pulsatile T-pipe study.

Every figure below is produced automatically by
`python3 analysis/figures.py --root . --out figures` (PDF for LaTeX, PNG for Word),
except the ParaView screenshots, which are flagged **[ParaView]** and which you
make by hand — instructions in §4.

Numbers to quote are in `results/*.csv`. **Quote from the CSV, never from this
document** — the values here are from a partial run and will change once the full
sweep finishes.

---

## 1. Page budget

| § | Section | Pages | Figures | Tables |
|---|---|---|---|---|
| 1 | Introduction | 1.5 | F1 | — |
| 2 | Theory | 2.5 | F3 | T1 |
| 3 | Numerical method | 2 | F2 | — |
| 4 | Setup and run matrix | 1.5 | F2 | T2 |
| 5 | **Verification** | 3 | F4 F5 F6 | T3 |
| 6 | **Validation** | 2.5 | F7 F8 F9 | T4 |
| 7 | **Results** | 4 | F10 F11 F12 F13 F14 | T5 |
| 8 | Extensions | 1.5 | F15 F16 | T6 |
| 9 | Discussion + conclusions | 1.5 | — | — |
| | **Total** | **20** | 16 | 6 |

To land at **15 pages** instead: cut §2 to 1.5 (state results, cite the derivations),
§3 to 1.5, §7 to 3 (drop F13, F14), §8 to 1. Do **not** cut §5 or §6 — verification
and validation are where the marks are.

§§5–7 are the graded core. Roughly half the report should be there.

---

## 2. Figure inventory

### Section 1 — Introduction

**F1 — the geometry** **[ParaView]**
A single clean render of the T-bifurcation with the inlet, both outlets, and the
junction labelled, plus the coordinate origin. Annotate `DH`, `DL1`, `DL`, and the
two outlet branches. This is the reference every later figure points back to.
*Source:* `output_A_Re100_al5/WaterBody_0000000000.vtp` + `WallBoundary_*.vtp`.

### Section 2 — Theory

**F3 — analytical Womersley profiles at several α**
Pure theory, no simulation — shows what the reader should expect before any results.
*Source:* generate with three lines of Python:
```python
import sys; sys.path.insert(0, "analysis")
import womersley, numpy as np
# velocity_hat(y, R, omega, nu, rho, G) at alpha = 1, 2, 5, 10
```
Show the parabola collapsing to a plug and the phase lag growing.

**T1 — analytical reference values.** From `womersley.flow_phase_lag_deg` and
`centreline_amplitude_ratio`:

| α | phase lag | peakedness \|u(0)\|/\|u_mean\| | annular peak \|y\|/R |
|---|---|---|---|
| 0.5 | 5.7° | 1.50 | 0.00 |
| 2 | 57.2° | 1.48 | 0.00 |
| 5 | 80.7° | 1.21 | 0.33 |
| 10 | 85.7° | 1.07 | 0.68 |

### Section 3 — Numerical method

**F2 — particle discretisation close-up** **[ParaView]**
Zoom on the junction showing fluid particles, wall particles, and the inflow buffer,
coloured by `Indicator`. Demonstrates you understand what SPH particles *are* and
where the boundary treatment lives. Reuse the same render for §4.

### Section 4 — Setup

**T2 — the run matrix.** Reproduce the table from `README.md` §"The study design"
(Re × α × A, plus the extension groups). State that ω is *derived* from α, and give
the derived T and cycle count per case — read them from each
`output_<case>/case_params.json` (`omega`, `T`, `n_cycles`, `decay_cycles`).

### Section 5 — VERIFICATION (the biggest single grade lever)

**F4 — `verification/convergence_functionals.pdf`**
Error vs particle spacing on log-log with the 2nd-order reference slope.
*Say:* inlet flow rate converges at observed order ≈ 2.7 with GCI ≈ 0.04 %, against
the exact value `U·D = 3`. **Also say** that wall shear did *not* converge
monotonically and that the analysis reports NaN rather than a fabricated order —
see §5 discussion below.

**F5 — `verification/periodicity_epsilon_decay.pdf`**
ε_N vs cycle number, log scale, with the 1 % criterion and the predicted
`0.4α²` cycle count marked.
*Say:* this is *cycle-to-cycle periodicity* — reaching the periodic steady state and
discarding the start-up transient. Distinguish it explicitly from grid convergence.

**F6 — `verification/periodicity_cycles_vs_alpha.pdf`**
Cycles needed vs α, measured against the `0.4α²` law.
*Say:* the slowest channel start-up mode decays at `λ₁ = ν(π/2R)²`, and
`τ_visc/T = α²/2π`, so cycles-to-periodicity scales as `0.4α²` — 2 cycles at α = 2,
40 at α = 10. **This is why a fixed cycle count cannot work across the matrix.**
This figure is your own contribution and is worth a full paragraph.

**Also include** `verification/wcsph_validity_checks.pdf` as a small inset or
appendix figure: Mach < 0.1 and density deviation, per case.

**T3 — verification summary.** From `results/convergence_gci.csv` and
`results/periodicity_and_wcsph_checks.csv`: per functional, the values at each dp,
refinement ratios, observed order, Richardson extrapolation, GCI; and per case,
ε_Q, ε_KE, max Mach, max Δρ/ρ, δ/dp.

### Section 6 — VALIDATION

**F7 — `validation/womersley_channel_profiles.pdf`** ← *the single most important figure*
Simulated vs analytical profiles at 8 phases, for α = 2 / 5 / 10.
*Say:* this is the **dedicated benchmark**, not the T-pipe — a periodic channel with
no velocity boundary condition anywhere, driven only by `G·cos(ωt)`, which is exactly
the configuration Womersley's solution describes. Agreement is 0.96 / 1.63 / 3.25 %
in L2. Point out core flattening and the annular effect appearing as α grows.

**F8 — `validation/womersley_channel_harmonic.pdf`**
Amplitude and phase across the channel, first harmonic.
*Say:* centre amplitudes match within 0.2 %, phases within 1°, and the **Richardson
annular peak position matches exactly** (0.00 / 0.34 / 0.68).

**F9 — `validation/mass_conservation_waveforms.pdf`**
`Q₁` vs `Q₂+Q₃` over the cycle.
*Say:* rms closure error and, importantly, **why it is not identically zero
instantaneously** — WCSPH is weakly compressible, so the junction stores volume and
the imbalance is physically `dV_stored/dt`. The cycle-*mean* error must vanish, and
does. Support with `validation/mass_conservation_convergence.pdf`.

**T4 — validation summary.** From `results/womersley_channel_validation.csv` and
`results/mass_conservation_and_split.csv`.

**A paragraph you must include** (it is the difference between an honest report and
an overclaiming one): *the T-pipe cannot validate the Womersley solution.* Its inlet
profile is prescribed parabolic, and the inlet channel is 3.5 long while the viscous
entrance length is `0.04·Re·D` = 12 and the oscillatory development length is
`U/ω` = 3.0. Measured peakedness at the inlet station is 1.41 against 1.50 for a
parabola and 1.21 for Womersley — i.e. still essentially the imposed profile. That is
precisely why the separate benchmark exists. Use
`validation/tpipe_inlet_profile_relaxation.pdf` here, labelled as *relaxation*, not
validation.

### Section 7 — RESULTS (needs the full sweep)

**F10 — `results/factorial_peak_tawss.pdf`** — Re × α heat map of peak TAWSS
**F11 — `results/factorial_max_osi.pdf`** — Re × α heat map of max OSI
*These only appear once group A's 9 cases exist.*

**F12 — `results/wall_tawss_map.pdf`**
Wall probes coloured by TAWSS on the real geometry. The clearest single hemodynamics
figure in the study.
*Say:* identify the stagnation point on the divider wall (TAWSS minimum) and the
accelerating regions either side of it.

**F13 — `results/wall_osi_distribution.pdf`** + `wall_tawss_distribution.pdf`
TAWSS and OSI along each wall.
*Say:* low TAWSS + high OSI marks atherosclerosis-prone regions (Ku & Giddens 1983;
Zarins & Glagov 1985) — in a bifurcation, the flow-divider apex and the outer walls
just downstream of the junction. **This is the paragraph that makes it a biofluid
report rather than a generic CFD report.** Note OSI ≈ 0.001 for steady flow vs up to
0.22 pulsatile: the index correctly vanishes when nothing oscillates.

**F14 — `results/pulsatility_attenuation.pdf`** — PI damping through the bifurcation
**Also** `results/amplitude_sweep_flow_reversal.pdf` from group B.

**T5 — results summary.** From `results/wall_metrics.csv` (TAWSS, OSI, RRT, apex
values per case) and `results/lumped_RL_fit.csv` (fitted R and L vs the analytic
`12μL/h³` and `ρL/h`).

**Optional depth if you have the pages:** the lumped RL fit. R comes out within ~5 %
of the analytic Poiseuille resistance; L about 56 % above `ρL/h`, plausibly because
the junction and entrance region add inertance. Connects the distributed CFD result
back to the lumped models from the lectures.

### Section 8 — EXTENSIONS

**F15 — `validation/flow_split_vs_branch_ratio.pdf`**
Measured split vs the analytic resistance-network prediction `b/(1+b)`.
*Say:* two branches in parallel divide flow inversely to their resistances
`R_h = 12μL/h³`; with only the length varied this collapses to `Q_up/Q_lo = b`. A
quantitative validation against a lumped-parameter model.

**F16 — stenosis series** — `extensions/stenosis_peak_wss.pdf` plus a **[ParaView]**
velocity-magnitude render of the throat jet at 50 % or 70 % occlusion.
*Say:* peak WSS vs area reduction; locate the jet and any recirculation downstream.

**T6** — from `results/stenosis_series.csv`.

### Section 9 — Discussion

No new figures. Cover the limitations honestly — this reads as competence:
- **2D, not 3D** — no secondary/helical flow, which is significant in real bifurcations
- **Laminar** — Re ≤ 200 here; aortic Re is ~4000
- **Rigid walls** — no compliance, so no wave propagation
- **Newtonian** — no shear-thinning (see `--carreau`, reserved but unimplemented)
- **Unidirectional inlet** — A ≤ 1 only; genuine flow reversal needs a bidirectional buffer
- **WSS not converged** — report TAWSS/OSI as relative comparisons, not absolute values

---

## 3. The three negative results — include them

Markers reward knowing the limits of your own method. Each of these is a real finding:

1. **The T-pipe cannot validate Womersley** (§6) — diagnosed from entrance-length
   scaling, and the reason a second benchmark case was built.
2. **Wall shear did not converge monotonically** (§5) — dp = 0.075 sits ~35 % above
   its neighbours in *both* the peak and the 95th percentile, so it is systematic
   across the whole wall, not one noisy probe. The analysis reports NaN rather than
   fabricating an order.
3. **A > 1 is unsupported** (§9) — at A = 1.5 the inlet velocity reverses but the
   emitter/buffer is unidirectional; the run collapses to Dt ≈ 3 × 10⁻⁵. Bounded
   experimentally: A = 0.25 / 0.75 / 1.0 all run cleanly.

Two more worth a sentence each in §3 or §5, as evidence of careful numerics:
- The continuity step had to be switched to the **Riemann** variant; with the
  baseline's `NoRiemann` 1 run in 3 collapsed depending only on thread scheduling.
- **`DH/dp` must be an integer** or the leftover partial particle layer biases every
  volume-weighted average — the flow-rate error flipped sign between resolutions
  until this was fixed.

---

## 4. ParaView screenshots — how to make them

Four screenshots (F1, F2, F12-companion, F16), plus the phase montage below.

### First: the case must have been run with `--vtp_all=1`

By default VTPs are written only inside the analysis window, which is fine for
figures but makes a broken-looking animation — the series jumps from the t = 0
at-rest frame straight to fully developed flow. For anything you intend to animate
or screenshot, run a dedicated case with `--vtp_all=1` (see `HANDOFF.md` §4b).
`output_ANIM/` is already prepared this way: 142 continuous snapshots.

### Open the grouped `.vtp` series

SPHinXsys writes **one VTP per snapshot**, named by physical time padded to ten
digits — `WaterBody_0113107344.vtp` is t = 113.107. Same convention as the
`test_2d_dambreak` and `test_2d_T_shaped_pipe` benchmarks.

In ParaView's file dialog these collapse into one entry displayed as

```
WaterBody_0000000000.vtp*          <- the asterisk means "file series"
```

**Select that grouped entry** (not an individual file) and Apply — you get the whole
animation. Load `WallBoundary_0000000000.vtp` alongside it; the wall is static, so it
is genuinely a single file.

> Expanding the group and picking `WaterBody_0000000000.vtp` on its own gives the
> t = 0 frame, where the fluid is at rest and nothing appears to move. That is the
> most common way to conclude the simulation failed when it did not.

The one cost is that the time slider shows a bare frame **index**. To find the frames
you want:

```bash
python3 analysis/frames.py output_ANIM --quarters   # the four montage frames
python3 analysis/frames.py output_ANIM              # full index -> time -> phase table
```

`frames.py` also reports whether the snapshot series is continuous — a large max gap
means the case was run without `--vtp_all=1` and will animate badly.

(A `.pvd` collection would label the slider with real time, and `analysis/make_pvd.py`
can still write one, but ParaView did not read them reliably here. The grouped series
is the supported path.)

### Render settings

1. **Open** the grouped `WaterBody_0000000000.vtp*` **and**
   `WallBoundary_0000000000.vtp`, then **Apply**.
2. **Representation:** `Points`; raise **Point Size** to ~4–6 so the particles read as
   a filled body rather than a dot cloud.
3. **Colour by:** `Velocity` (Magnitude) for flow figures, `Pressure` for pressure
   figures, `Indicator` for the discretisation figure (F2). Also available:
   `Density`, `VelocityGradient`.
4. **Rescale the colour bar to a FIXED range** across any snapshots you intend to
   compare — use *Rescale to Custom Data Range*, not *Rescale to Data Range*.
   Otherwise every frame rescales itself, the flow looks identical at every phase,
   and the comparison is meaningless.
5. **Camera:** View ▸ Camera ▸ `-Z` for a straight-on 2D view. Turn off the
   orientation axes and the ParaView logo for a clean figure.
6. **Save:** File ▸ Save Screenshot, PNG, at least 1600 × 1200.

### The phase montage — the most convincing figure in the report

Use a case run with `--vtp_all=1`; `output_ANIM/` is already prepared. Ask
`frames.py` for the four frames of the last full cycle:

```
$ python3 analysis/frames.py output_ANIM --quarters
   panel   frame   time        phase
       1     121     113.107     4.001
       2     126     117.812     4.250
       3     131     122.533     4.501
       4     136     127.243     4.750
```

Type each **frame** number into ParaView's Time box, colour by velocity magnitude on
the same fixed range, save each, and lay them out 2 × 2 labelled by phase
(0.00 T, 0.25 T, 0.50 T, 0.75 T). That single figure shows the pulsatile cycle better
than any line plot.

---

## 5. Practical notes

- **PDF for LaTeX, PNG for Word** — both are written for every figure.
- Every figure has its underlying numbers in `results/*.csv`, so tables and figures
  are guaranteed consistent, and you can state that in the methods section.
- **Report everything dimensionlessly** — `u/U_f`, `t/T`, `p/(ρU_f²)`,
  `τ_w/(μU_f/R)`. That is what makes it a parameter study rather than a set of runs.
- If a figure is missing after a run, the case it needs probably failed — check
  `manifest.csv` and the matching `log_<case>.txt`.
- Cite SPHinXsys, and note the case derives from the `T_shaped_pipe` benchmark.
