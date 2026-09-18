# Handoff / work log — pulsatile T-pipe parameter study

**Date:** 2026-09-19 · **Status:** COMPLETE - all 21 study runs executed, 0 failures

This document records what was planned, what was actually executed, what was
found, and what remains. It is written so that (a) a new chat session can pick the
work up cold, and (b) the findings can go straight into the report.

Read `README.md` first for how to build and run. This file is the *log*.

**Companion documents:**
- `RUN_EVERYTHING.sh` — one resumable command that produces every number and figure
- `REPORT_GUIDE.md` — what to put in the report, figure by figure, with sources

---

## 0. One-paragraph summary

Your original single-case file was extended into a scriptable, dimensionless
parameter study plus a separate analytical benchmark. **Your original file was not
modified.** The C++ solver, the sweep driver, and a 10-module Python analysis
package are complete and working end to end. The physics has been validated
against the analytical Womersley solution to 1–3 %. Four of the twenty-one study
runs have been executed (enough to verify every code path and to produce a
complete convergence study); the remaining seventeen are compute time, not work.

---

## 1. What was planned vs. what was executed

### 1.1 Code — all planned items delivered

| Planned | Status | File |
|---|---|---|
| CLI-driven case definition (Re, α, A, dp, …) | **Done** | `case_params.h` (489 ln) |
| Study solver with flow-rate reduction, WSS probes, decoupled sampling | **Done** | `T_pipe_pulsatile_parameter_study.cpp` (969 ln) |
| Asymmetric branch (extension D) | **Done, smoke-tested** | same |
| Cosine stenosis (extension E) | **Done, smoke-tested** | same |
| 22-run sweep driver | **Done** (now 21, see §3.4) | `run_study.sh` (257 ln) |
| Python analysis + figures | **Done** | `analysis/*.py` (2 054 ln) |
| Womersley benchmark | **Added — not in original plan** | `womersley_channel_validation.cpp` (375 ln) |

Total ≈ 4 850 lines. Everything compiles clean; all four CMake targets build.

### 1.2 Deviations from the original plan, and why

1. **A separate Womersley benchmark executable was added.** The plan assumed the
   T-pipe inlet profile could be compared against the Womersley solution. It
   cannot — see §3.1. This was the single largest addition and it is what gives
   the report a real validation section.
2. **Non-Newtonian (Carreau) rheology was dropped.** You selected asymmetric
   branch + stenosis as the two extensions, so this was never built. A `--carreau`
   flag is parsed and reserved but does nothing. The API exists
   (`CarreauViscosity`, `ShearRateDependentViscosity`,
   `NonNewtonianViscousForceWithWall`) — see
   `tests/2d_examples/test_2d_lid_driven_cavity_non_newtonian/lid_driven_cavity.cpp:137,161-163`.
3. **The run matrix shrank from 22 to 21 cases**: A = 1.5 was removed (§3.4).
4. **Production resolution changed from dp = 0.15 to dp = 0.05**, forced by the
   Stokes-layer criterion (§2.2).

### 1.3 Runs actually executed - ALL 21, zero failures

| Group | Cases | What it gives |
|---|---|---|
| A | 9 | Re (50/100/200) x alpha (2/5/10) factorial at A = 0.5, dp = 0.05 |
| B | 4 | amplitude sweep A = 0, 0.25, 0.75, 1.0 |
| C | 2 (+A anchor) | resolution series dp = 0.05 / 0.075 / 0.10 |
| D | 2 | asymmetric branch, b = 0.5 and 0.7 |
| E | 3 | stenosis, 30 / 50 / 70 % |
| F | 1 | CFL insensitivity (half time step) |
| benchmark | 3 | Womersley channel, alpha = 2 / 5 / 10 |

Slowest case: `A_Re200_al2` at 22 302 s (its period is T = 236, so one run covers
~1650 time units). Total sweep ~29 h CPU, run 4-way parallel.

**Auxiliary runs live in `bin/aux/`** - the animation case (`ANIM`), the two
geometry checks (`GEOM_*`) and the early smoke tests (`S_*`). They were moved out
of the top level because `load_study` globs `output_*` and they were contaminating
the analysis: `ANIM` duplicated `C_dp0.10`'s resolution in the convergence series,
and `ANIM`/`S_cfl`/`C_*` all collided on the (Re 100, alpha 5) cell of the
factorial. The analysis is now also hardened against this (one case per resolution;
factorial restricted to the production dp), but keep auxiliary runs in `aux/`.

### 1.4 Not executed

Nothing from the study matrix. Optional extras only:

- **T4** `SectionMeanPressure` reducer - would fix the noisy `eps_dp` diagnostic
- **T5** wall-shear convergence - still unresolved, see §3.3
- **T6** Carreau non-Newtonian rheology - `--carreau` reserved but unimplemented
- **T7** bidirectional inlet for A > 1 - a project, not a task

## 2. Design decisions worth defending in the report

### 2.1 ω must be derived from α, never prescribed

`Re` and `α` are coupled through ν, so they cannot be varied independently by
editing ω:

```
ν = U·DH/Re        ω = ν(α/R)²        T = 2π/ω
```

Your original file set ω directly, which is why the α values in its header comment
do not match its own formula.

### 2.2 Run length scales as α², and resolution is set by the Stokes layer

The slowest start-up mode of a plane channel decays at `λ₁ = ν(π/2R)² =
(π²/4)/τ_visc`, and `τ_visc/T = α²/2π`, so

> **cycles to periodicity ≈ 0.4 α²** — 2 cycles at α = 2, **40 at α = 10**

Your original `end_time = t_ramp + 5·T_pulse` is far too short at high α. The
solver now derives the cycle count per case.

The Stokes layer `δ = √(2ν/ω) = R√2/α` is **independent of Re** and is the
smallest physical length in the problem. Resolving it with ≥ 4 particles requires
`dp ≤ 0.795/α`, so the largest α sets dp for the whole factorial: **dp = 0.05**
gives δ/dp = 4.2 at α = 10. Group A must use a single dp or resolution becomes
confounded with α — the very thing the factorial exists to separate.

### 2.3 Two constraints on convergence-study resolutions

1. **`DH/dp` must be an integer.** The lattice fills the channel with whole
   particle layers; a partial layer biases every volume-weighted average. At
   dp = 0.1125 (26.67 layers) and 0.16875 (17.78) the inlet flow-rate error jumped
   to +2.07 % and −2.98 %, *flipping sign*. At 60/40/30 layers it is −0.32 %,
   −0.38 %, −0.49 %. The solver now warns when `DH/dp` is not an integer.
2. **Every member must be inside the asymptotic range.** dp = 0.20 gives δ/dp = 2.1
   and an implied order of 4.6 — not believable for a second-order scheme.

Those two constraints leave an *unequal* ratio ladder (1.5 then 1.333), so
`convergence.py` solves the **ASME V&V 20 / Celik unequal-ratio relation** by
fixed-point iteration rather than assuming a constant ratio. Self-tested: it
recovers order 2.0000 exactly on synthetic 2nd-order data for all three ladders.

---

## 3. Findings — the substance for the report

### 3.1 The T-pipe cannot validate the Womersley solution

This is the most important finding and it changes what the report can claim.

The T-pipe inlet profile is **prescribed** as a parabola by the inflow buffer, and
the inlet channel is only 3.5 long while

- viscous entrance length ≈ `0.04·Re·D` = **12**
- oscillatory development length ≈ `U/ω` = **3.0**

so the flow reaches the junction still carrying the imposed profile. Measured
peakedness at the inlet station is **1.41**, against 1.50 for a pure parabola and
1.21 for the Womersley prediction at α = 5 — i.e. much closer to the imposed
profile than to theory, exactly as those length scales require.

**Consequence:** the `tpipe_*_vs_theory` figures are labelled as profile
*relaxation*, not validation, and the analysis module carries the caveat in its
docstring. Do not present them as validation.

### 3.2 The dedicated benchmark validates the scheme properly

`womersley_channel_validation.cpp` — straight channel, periodic in x, **no velocity
BC anywhere**, driven only by `G·cos(ωt)`. This is exactly the configuration
Womersley's solution describes, so everything in the interior is genuinely computed.

| α | δ/dp | L2 | centre \|û\| sim / theory | centre phase sim / theory | annular peak sim / theory |
|---|---|---|---|---|---|
| 2 | 14.1 | 0.96 % | 1.4846 / 1.4830 | −61.1° / −60.9° | 0.00 / 0.00 |
| 5 | 5.7 | 1.63 % | 1.2093 / 1.2114 | −91.4° / −91.2° | 0.34 / 0.34 |
| 10 | 2.8 | 3.25 % | 1.0714 / 1.0717 | −90.8° / −89.9° | 0.68 / 0.68 |

Core flattening, the phase approaching the 90° limit, and the **Richardson annular
effect** are all reproduced — the annular peak position matches exactly. The error
growing as δ/dp falls independently confirms the ≥ 4 criterion of §2.2.

Analytical reference values (`womersley.flow_phase_lag_deg`), useful as a
report table: lag = 5.7° (α=0.5), 57.2° (α=2), 80.7° (α=5), 85.7° (α=10).

### 3.3 Verification and results — FINAL NUMBERS (full 21-case study)

**Convergence** (dp 0.05 / 0.075 / 0.10, Celik unequal-ratio):

| Functional | values | order | extrapolated | GCI |
|---|---|---|---|---|
| `Q_mean` | 2.99051, 2.98858, 2.98519 | **2.69** | **2.9915** (exact = 3) | **0.04 %** |
| `tawss_p95` | 0.1813, 0.2422, 0.1703 | NaN | — | — |
| `flow_split` | 0.49969, 0.49952, 0.49866 | NaN | — | — |

**Periodicity and WCSPH validity** across all 21 cases:

- `eps_Q` 0.0069–0.0145, `eps_KE` 0.0025–0.0062 → 14/21 formally periodic at the
  1 % criterion, the rest marginally above. Both are an order of magnitude better
  than `eps_dp` (0.3–0.9), which is why the verdict now uses Q and KE only.
- Mach ≤ 0.096 everywhere **except the stenosis cases**: `E_st0.5` reaches 0.113
  and `E_st0.7` reaches 0.131, because the throat accelerates the flow. Those two
  marginally violate weak compressibility — state it as a limitation.
- Density deviation < 2 % everywhere except `A_Re50_al10` at 3.11 %, the
  worst-resolved case in the matrix.

**Conservation, across all 21 cases:**

- inlet flow-rate error **−0.14 % to −0.66 %** against the exact `U·D = 3`
- rms closure error 0.50–0.79 % (one outlier: `A_Re50_al10` at 6.86 %)
- flow split 0.4995–0.5006 for every symmetric case

**Pulsatility index tracks theory exactly** — this is a clean validation figure:

| A | 0 | 0.25 | 0.5 | 0.75 | 1.0 |
|---|---|---|---|---|---|
| PI measured | 0.020 | 0.507 | 1.012 | 1.516 | 2.021 |
| PI = 2A | 0 | 0.5 | 1.0 | 1.5 | 2.0 |

**Re × α factorial** (peak TAWSS, dimensional): TAWSS falls with Re
(0.33 → 0.21 → 0.135 for Re 50 → 100 → 200) and is nearly independent of α.
Max OSI is ~0 at α = 2 and 5 for Re 50–100, but jumps to 0.34–0.50 at α = 10 and
at Re = 200 — i.e. **flow reversal near the walls appears at high α and high Re**,
which is the main physical result of the factorial.

**Still unresolved: wall shear does not converge.** Both the raw peak
(0.208 / 0.311 / 0.188) and the 95th percentile (0.181 / 0.242 / 0.170) put
dp = 0.075 ~35 % above its neighbours, so it is systematic across the whole wall,
not one noisy probe. The analysis reports NaN rather than fabricating an order.
Report TAWSS/OSI as **relative comparisons between cases at the production
resolution**, never as converged absolute values.

### 3.3b Extension D: the lumped model fails, and the failure is quantitative

The naive resistance-network prediction `split = b/(1+b)` is **wrong by 15 %**:

| b | measured | naive prediction | R₀/R_branch | re-predicted with R₀ |
|---|---|---|---|---|
| 0.50 | 0.4809 | 0.3333 | 5.80 | 0.4797 |
| 0.70 | 0.4872 | 0.4118 | 4.99 | 0.4880 |

The naive form assumes the only resistance in each path is the fully developed
Poiseuille resistance `12μL/h³` of that branch. But the branches have
**L/h = DH/(DL−DL1) = 2** — far too short for Poiseuille flow to dominate — so much
of the pressure drop happens in the junction and entrance and is *common to both
paths*. Adding a shared series resistance R₀,

```
split = (R0 + b R) / (2 R0 + (1 + b) R)
```

and inverting at each measured split gives **R₀ ≈ 5 R from two independent cases**
(agreeing to 14 %), and re-predicting with the mean R₀ reproduces both splits to
~0.2 %. So the junction carries roughly five times the resistance of a branch.

This is a **better report result than the original prediction would have been**:
it shows when a lumped model applies (L/h ≫ 1) and what to do when it does not.
Computed by `conservation.series_resistance_fit`.

### 3.3c Extension E: stenosis diverts flow, strongly and monotonically

| stenosis | Q_upper fraction | PI_upper |
|---|---|---|
| 0 % | 0.4997 | 1.024 |
| 30 % | 0.4746 | 1.049 |
| 50 % | 0.4116 | 1.102 |
| 70 % | 0.2706 | 1.214 |

A 70 % occlusion pushes the upper branch from half the flow down to 27 %, and the
flow that still gets through is *more* pulsatile (PI rises 1.02 → 1.21).

### 3.4 A > 1 is not supported — flow reversal needs a bidirectional inlet

At A > 1 the prescribed inlet velocity goes **negative** over part of the cycle,
i.e. fluid must flow back out through the inlet. The emitter + inflow-buffer
combination inherited from `T_shaped_pipe` is unidirectional.

Tested at dp = 0.10: A = 0.25, 0.75, 1.0 all run cleanly (Dt ≈ 0.015–0.027);
**A = 1.5 collapses to Dt ≈ 2.7 × 10⁻⁵ and stalls at t/T = 2.4.**

The solver now refuses A > 1 with an explanatory message, and group B uses
A ∈ {0, 0.25, 0.75, 1.0}. Studying genuine reversal needs the bidirectional buffer
in `tests/extra_source_and_tests/extra_src/shared/pressure_boundary/` — separate
work, not a parameter change.

### 3.5 Numerical robustness problems found and fixed

| Problem | Symptom | Fix |
|---|---|---|
| **Flaky runs** | 1 in 3 identical runs collapsed to dt ≈ 5e-5 and exhausted the particle buffer, depending only on thread scheduling | Switched the continuity step from `Integration2ndHalfWithWallNoRiemann` to `…WithWallRiemann`. 3/3 stable; rogue particles gone (max speed 2.37 vs 70–418) |
| **`MaximumSpeed` is misleading** | Reduces over buffer/disposer particles, reporting 70–418 in a flow of order 1. Sizing `c_f` from it inflated it 2× and *destabilised* the run (2/3 failures at c_f = 45, 0/3 at 22.5) | Added `MaxBulkSpeed` (interior particles only); `CORNER_FACTOR` calibrated to 1.25 |
| **Density check measured the wrong thing** | 7.9–67 % deviation, dominated by free-surface particles with truncated kernel support | `MaxDensityDeviation` restricted to `Indicator == 0`; now 0.8–1.9 % |
| **Parallel runs raced** | `Abort trap: 6` at start-up; every process creates `./output`, `./restart`, `./sphinxsys.log` in the CWD | `run_study.sh` gives each case its own `.work_<case>/` dir |
| **`GROUPS` is a bash builtin array** | `-g` selection silently produced an empty matrix | Renamed `SEL_GROUPS` |
| **Asymmetric measurement stations** | Upper/lower flow-rate stations were at different distances from the junction, biasing the `Q₁ = Q₂ + Q₃` imbalance | Both now at the same fraction (0.50) of their branch |
| **Womersley sign error** | Analytical lag came out −99°, must lie in [0°, 90°] | Particular solution is `G/(iωρ) = −iG/(ωρ)`, not `+i` |
| **Cycle-folding rejected every cycle** | Coarsely-sampled signals (50/cycle) never bracket the phase grid → silent NaN everywhere | `interp_periodic` wraps the waveform periodically |
| **Time-dependent forcing applied once** | Benchmark had zero oscillation; first harmonic was noise | `apply_forcing.exec()` moved inside the time loop |
| **Wall-coincident probes** | Observer exactly on the wall reads \|u\| = 0.17 where theory says 0 (one-sided stencil), inflating L∞ from 0.3 % to 14 % | Probes inset by 0.5 dp |

---

## 4. Verified-correct things NOT to "fix"

- `position[1]` in `InflowVelocity` is in the **inlet-box local frame**
  (`InflowVelocityCondition::update` calls `transform_.shiftBaseStationToFrame`),
  so `1 − y²/R²` peaks on the centreline. Correct in your original too.
- `createProfileObserverPoints` uses **global** coordinates — also correct.
- `NaN` in `convergence_gci.csv` is deliberate: it means the series does not
  support Richardson extrapolation.
- `flow_split` reporting a NaN order is correct — it is already converged to
  ~0.1 % at every resolution, so its differences are noise.

---

## 4b. OUTPUT CONVENTION — follow this for every case, always

Benchmarked against the two stock SPHinXsys cases that animate correctly in
ParaView: `build/tests/2d_examples/test_2d_dambreak/bin/output/` and
`test_2d_T_shaped_pipe/bin/output/`. Any new case in this project must produce the
same file set, for the same reasons.

### The canonical structure

```
output[_<case>]/
    ShapeSPHSystemDomain.vtp        domain bounding shape, written once
    WallBoundary_0000000000.vtp     static body -> exactly ONE file
    WaterBody_0000000000.vtp        \
    WaterBody_0000536250.vtp         |  the time series: one VTP per snapshot,
    WaterBody_0001072500.vtp         |  named by physical time x 1e6, zero-padded
    ...                             /   to 10 digits
    WaterBody_<Quantity>.dat        reduced scalar histories
    <Observer>_<Quantity>.dat       observer time series
```

**Rules, each with the failure it prevents:**

1. **One VTP per snapshot, named `<Body>_<t x 1e6 padded to 10>.vtp`.** This is what
   `BodyStatesRecording::writeToFile()` (no argument) does. Do not pass an iteration
   number — the `writeToFile(size_t)` overload names files `ite_*` instead, which
   breaks the time series. `test_2d_T_shaped_pipe` writes 197 files spanning
   t = 0 → 100.06 this way.

2. **A static body gets exactly one file.** `WallBoundary_0000000000.vtp` is written
   once. Do not re-emit it per snapshot.

3. **Snapshots must be CONTINUOUS from t = 0.** This is the rule most easily broken
   here. Writing VTPs only inside the analysis window produces a series that jumps
   from the t = 0 at-rest frame straight to a fully developed flow — ParaView shows
   two disconnected states and the animation looks broken. The study defaults to
   window-only to keep 21 cases' disk use sane, so **any case you intend to animate
   or screenshot must be run with `--vtp_all=1`**, which restores continuous
   coverage (142 snapshots, max gap 0.97, for a 5-cycle case).

4. **Emit `ShapeSPHSystemDomain.vtp`.** The `SPHSystem` constructor writes it into
   the default `./output`, so any case that calls `resetOutputFolder()` destroys it
   and must call `sph_system.writeSystemDomainShapeToVtp()` afterwards.

5. **Open the grouped `.vtp` series in ParaView — not a `.pvd`.** In the file
   dialog the snapshots collapse into a single entry shown as
   `WaterBody_0000000000.vtp*`; selecting that loads the whole series. This is how
   the stock benchmarks are used and it is the supported path here.

   `analysis/make_pvd.py` exists and can write `.pvd` collections that would label
   the time slider with physical time or cycle phase, but **ParaView did not read
   them reliably in practice**, so it is no longer run by `RUN_EVERYTHING.sh` and
   nothing depends on it. The grouped series works; leave it at that.

   The cost of the grouped series is that the slider shows a bare frame INDEX.
   `analysis/frames.py` prints the index -> time -> phase mapping:

   ```bash
   python3 analysis/frames.py output_ANIM --quarters   # the four montage frames
   python3 analysis/frames.py output_ANIM              # full table
   ```
   It also reports whether the series is continuous, which is the check in rule 3.

6. **Per-case output folders** (`output_<case>/`) are a deliberate deviation from
   the stock single `output/`, because a 21-case sweep needs isolation. The file set
   *inside* each folder is identical to the stock convention.

### Making an animation case

```bash
./test_2d_T_pipe_pulsatile_parameter_study \
    --Re=100 --alpha=5 --A=0.5 --dp=0.10 --cycles=5 --vtp_all=1 --case=ANIM
python3 analysis/make_pvd.py output_ANIM
```
Then open `output_ANIM/WaterBody_phase.pvd` in ParaView.

### Checking a case conforms

```bash
ls output_<case>/ | sed 's/_[0-9]\{10\}\.vtp/_<time>.vtp/' | sort -u
```
Expect `ShapeSPHSystemDomain.vtp`, one `WallBoundary_<time>.vtp`, one
`WaterBody_<time>.vtp`, plus `.pvd` and `.dat` files. Then confirm continuity:

```bash
python3 -c "
import glob; f=sorted(glob.glob('output_<case>/WaterBody_*.vtp'))
t=[int(x[-14:-4])/1e6 for x in f]
print(len(f),'snapshots, t =',t[0],'->',t[-1],'max gap',max(b-a for a,b in zip(t,t[1:])))"
```
A max gap comparable to the mean spacing means the series is continuous. A gap of
hundreds of time units means VTPs were confined to the analysis window — re-run
with `--vtp_all=1`.

---

## 5. How to continue (for a fresh session)

> **The actionable queue is §6.** This section is the context behind it; §7 covers
> how to phrase the prompts.

1. **Read** `README.md` (build/run) then this file (state and rationale).
2. **Run the study:**
   ```bash
   cd <build>/tests/2d_examples/test_2d_T_pipe_pulsatile_parameter_study/bin
   ./run_study.sh -n          # confirm 21 cases
   ./run_study.sh -j 4        # ~2-3 h
   python3 analysis/figures.py --root . --out figures
   ```
   Check `manifest.csv` for failures; each case has a `log_<case>.txt`.
3. **Expect** group A's α = 2 cases to be the slowest (T = 235 at Re 200, α 2).

### Known risks / open items

- **Groups D and E have never had a production run**, only 3-cycle smoke tests.
  The stenosis geometry is built by polygon subtraction; **dump the wall body to
  VTP and look at it** before trusting a full run. Wall probes inside the stenosis
  footprint are dropped automatically (they appear as `NaN`).
- **WSS convergence is unresolved** (§3.3). If the report needs converged absolute
  WSS, try a kernel-corrected velocity gradient
  (`VelocityGradientWithWall<LinearGradientCorrection>` plus
  `LinearGradientCorrectionMatrixComplex`) instead of `NoKernelCorrection`.
- **`eps_dp` (pressure-drop periodicity) is 0.3–0.9**, far worse than `eps_Q` and
  `eps_KE` (both < 2 %). Point pressure in WCSPH is noise-dominated. Use `eps_Q`
  and `eps_KE` as the periodicity criteria. A `SectionMeanPressure` reducer
  (mirroring `SectionFlowRate`) would fix this and also improve the RL fit — a
  clean ~20-line addition that was scoped but not built.
- **The RL lumped fit** gives R within 5 % of `12μL/h³` but L about 56 % above
  `ρL/h`; plausibly real (the junction and entrance add inertance) but unverified.
- **Carreau rheology** is reserved but unimplemented.

---

## 6. TASK BACKLOG — the ordered queue

Numbered, self-contained, dependency-ordered. Quote a number and a fresh session
knows exactly what to do: *"Read HANDOFF.md, then do tasks T1-T4 from §6."*

Tasks marked **[parallel-safe]** can be done while T2's sweep is running.
Effort is Claude's working time, not your waiting time.

---

> **T1-T3 are now automated by `RUN_EVERYTHING.sh`** (resumable, one command).
> The descriptions below remain the reference for what each step is for.

### T1 — Visually verify the stenosis and asymmetric-branch geometry
**Blocks:** T2 groups D and E · **Effort:** small · **Risk if skipped:** high

Groups D and E have only ever had 3-cycle smoke tests. The stenosis is built by
polygon subtraction from the fluid body and addition to the wall body; if the
throat is malformed or the wall is not closed, a full sweep wastes an hour per case.

```bash
cd <build>/tests/2d_examples/test_2d_T_pipe_pulsatile_parameter_study/bin
./test_2d_T_pipe_pulsatile_parameter_study --stenosis=0.7 --dp=0.05 --cycles=3 --case=geomchk
./test_2d_T_pipe_pulsatile_parameter_study --branch_ratio=0.5 --dp=0.05 --cycles=3 --case=geomchk2
```
Open `output_geomchk/WallBoundary_*.vtp` and `WaterBody_0000000000.vtp` in ParaView.
**Done when:** the throat is open, the wall is continuous around it, the throat is
≥ 6 particles across, and the shortened lower branch still contains its disposer.

---

### T2 — Run the full 21-case study
**Depends on:** T1 · **Effort:** small to launch, ~2–3 h wall time with `-j 4`

```bash
./run_study.sh -n            # confirm 21 cases
./run_study.sh -j 4
```
Group A's α = 2 cases are slowest (T = 235 at Re 200, α 2).
**Done when:** `manifest.csv` shows 21 ok, 0 failed. On failure read
`log_<case>.txt`; a collapsed `Dt` means the instability of §3.5, not a bad parameter.

---

### T3 — Generate all figures and tables
**Depends on:** T2 · **Effort:** small

```bash
python3 analysis/figures.py --root . --out figures
```
**Done when:** no `!!` lines; `figures/` and `results/` populated. Then sanity-check
against §3.3: inlet Q error < 1 % at dp = 0.05, split ≈ 0.5, PI ≈ 2A.

---

### T4 — Add a `SectionMeanPressure` reducer  **[parallel-safe]**
**Depends on:** nothing · **Effort:** medium · **Value:** high

Fixes the one weak diagnostic in the study. `eps_dp` is 0.3–0.9 because pressure is
sampled at a *single* centreline point and WCSPH point pressure is noise-dominated.
A slab-averaged pressure over the same cross-sections used for flow rate would cut
that noise by ~√N and simultaneously improve the RL fit of §5's open items.

Mirror `SectionFlowRate` in `T_pipe_pulsatile_parameter_study.cpp` (search for
`class SectionFlowRate`): same slab geometry, reduce `ReduceSum<Vecd>` of
`[p_i·Vol_i, Vol_i]`, divide in Python exactly as `io_utils.flow_rate` does. Record
it over the **whole run** (not just the analysis window) so it feeds periodicity.
Then use it in `periodicity.assess` and `lumped.fit` instead of
`case.pressure_drop`.

**Done when:** `eps_dp` drops below ~0.05 and `lumped.py`'s r² stays ≥ 0.99.

---

### T5 — Resolve the wall-shear convergence problem  **[parallel-safe]**
**Depends on:** nothing (re-runs 3 cases) · **Effort:** medium · **Value:** high

§3.3's one unresolved result. WSS is twice-derived and dp = 0.075 sits ~35 % above
its neighbours. Two things to try, in order:

1. **Kernel-corrected gradient.** Swap
   `VelocityGradientWithWall<NoKernelCorrection>` for
   `<LinearGradientCorrection>`, which additionally requires constructing
   `InteractionWithUpdate<LinearGradientCorrectionMatrixComplex>` and calling it
   before the gradient each step — see
   `tests/2d_examples/test_2d_lid_driven_cavity_non_newtonian/lid_driven_cavity.cpp:155,161`.
2. **Probe offsets.** They are currently 1 dp and 2 dp
   (`probes::buildWallProbes`). Try 1.5 dp and 2.5 dp so both sit fully inside the
   kernel support, and compare the extrapolation.

Re-run only the group C series (dp = 0.05 / 0.075 / 0.10) and re-check
`convergence_gci.csv`.
**Done when:** `tawss_p95` is monotone and returns a finite order, **or** you have
evidence it cannot be made monotone — which is itself a reportable result. Do not
force a number.

---

### T6 — Non-Newtonian (Carreau) blood rheology  **[parallel-safe]** — OPTIONAL
**Depends on:** nothing · **Effort:** medium · **Value:** extra report section

Only if the report needs more material. `--carreau` is parsed and reserved but does
nothing. The API exists: `addMaterialProperty<CarreauViscosity>(...)`,
`SimpleDynamics<fluid_dynamics::ShearRateDependentViscosity>`,
`NonNewtonianViscousForceWithWall<AngularConservative>` — see
`lid_driven_cavity.cpp:137,161-163`. Note `AdvectionViscousTimeStep` must then use
the maximum viscosity.
**Done when:** a Carreau run at matched Re shows the expected shear-thinning
reduction in near-wall viscosity vs the Newtonian baseline.

---

### T7 — Bidirectional inlet for genuine flow reversal — OPTIONAL, LARGE
**Depends on:** nothing · **Effort:** large · **Value:** only if reversal is required

§3.4: A > 1 needs an inlet that accepts backflow. The machinery is in
`tests/extra_source_and_tests/extra_src/shared/pressure_boundary/` (bidirectional
buffers, Windkessel BCs) but it is a different framework and lives outside the main
test tree. **This is a project, not a task.** Do not start it unless the report
specifically needs flow reversal — the A ≤ 1 sweep already covers the physiological
range for this geometry.

---

## 7. How to prompt from here

### First: most of what is left is NOT a prompt

The sweep is a terminal command, not a conversation. Do this yourself:

```bash
cd <build>/tests/2d_examples/test_2d_T_pipe_pulsatile_parameter_study/bin
nohup ./RUN_EVERYTHING.sh &
tail -f run_everything.log
```

Resumable — if the laptop sleeps or you Ctrl-C, run it again and it skips whatever
finished. Asking a chat session to babysit a 3-hour job wastes the session; use
Claude for the parts that need judgement.

### The one rule for a fresh chat

A new session has **no memory of this work** and will not find these files on its
own. Every opening prompt needs: **the path**, **which document**, **which section**.

```
Read HANDOFF.md and REPORT_GUIDE.md in
/Users/adheshsagar/CodeFiles/bfm11/sphinxsys/tests/2d_examples/test_2d_T_pipe_pulsatile_parameter_study/
then <what you want>.
```

### Prompts by where you are

**A. Sweep finished, want to know what you got**

> Read `HANDOFF.md` §3 and every CSV in `bin/results/`. Compare the actual numbers
> against the ones quoted in `HANDOFF.md` §3.3 and `REPORT_GUIDE.md` — those came
> from a partial run. Tell me what changed, what is now better or worse, and
> whether any case failed. Update §1.3 and §3.3 of `HANDOFF.md` with the real
> numbers.

**B. Some cases failed**

> `manifest.csv` shows these cases FAILED: <paste>. Read the matching
> `log_<case>.txt`, diagnose, and fix. Check `HANDOFF.md` §3.5 first — this may be
> a known failure mode rather than a new one.

**C. Writing a report section**

> Read `REPORT_GUIDE.md` §2. I am writing section <N>. Using the real numbers in
> `bin/results/*.csv`, draft that section: what the figures show, what the numbers
> are, and what to conclude. Quote only numbers that are actually in the CSVs -
> do not carry over any value from the guide without checking it.

That last sentence matters. The numbers in `REPORT_GUIDE.md` are from a partial
run and *will* be stale.

**D. Stuck on ParaView**

> Read `REPORT_GUIDE.md` §4. I want <the 2x2 phase montage / the stenosis throat
> render>. Walk me through it for the files actually in `output_<case>/`, and tell
> me which snapshot numbers to use.

**E. Remaining code work (T4/T5)**

> Read `HANDOFF.md` §6 and do T4. Do not touch anything listed in §4.

For T5, bound it explicitly — it is the one genuine rabbit hole:

> Do T5, but try only the two fixes listed. If neither makes `tawss_p95` monotone,
> stop and report that, rather than going deeper. A negative result is acceptable
> here and is already written into the report plan.

### Phrasing that works

- **Name the document and section**, not "the next few things". *"T1-T4 in §6 of
  HANDOFF.md"* beats *"the next 4 things"*.
- **Batch by dependency, not by count.** T4 and T5 are independent and marked
  parallel-safe; T1→T2→T3 is one chain.
- **Ask for the log to be updated** at the end of any session that changes state,
  or this document goes stale and the next session inherits a wrong picture.
- **Bound the rabbit holes** — say what "done" and "give up" both look like.
- **Say "do not fabricate numbers"** when asking for report text. It is the single
  most useful sentence once real results exist.
- **Carry the negative results forward.** §3.1, §3.3 and §3.4 are things a fresh
  session will otherwise rediscover the hard way — or silently "fix" back into a
  broken state. §4 lists what must not be touched.

## 8. Suggested report structure (≈20 pages)

| § | Content | Pages |
|---|---|---|
| 1 | Introduction: bifurcation hemodynamics, why pulsatility matters | 1.5 |
| 2 | Theory: Poiseuille, Womersley, α/Re/A groups, Stokes layer, TAWSS/OSI, lumped RL | 3 |
| 3 | Numerical method: WCSPH, Riemann fluxes, transport-velocity correction, inflow/outflow buffers | 2.5 |
| 4 | Setup: geometry, BCs, observers, run matrix | 1.5 |
| 5 | **Verification**: convergence + GCI (§3.3), periodicity and the 0.4α² law (§2.2), WCSPH validity, lattice/asymptotic constraints (§2.3) | 3 |
| 6 | **Validation**: Womersley benchmark (§3.2), steady Poiseuille, mass conservation; and *why the T-pipe cannot do this* (§3.1) | 2.5 |
| 7 | **Results**: Re × α factorial, amplitude sweep, PI attenuation, TAWSS/OSI, RL fit | 4 |
| 8 | **Extensions**: branch ratio vs resistance network, stenosis series | 1.5 |
| 9 | Discussion (2D, laminar, rigid-wall, Newtonian, unidirectional inlet §3.4) and conclusions | 1.5 |

§§5–7 are the marked core. The negative results (§3.3 WSS, §3.4 A > 1) are worth
real space — showing you found the limits of your own method reads as competence,
not failure.
