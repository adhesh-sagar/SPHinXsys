# Handoff / work log — pulsatile T-pipe parameter study

**Date:** 2026-08-18 · **Status:** framework complete and verified; 4 of 21 study runs executed

This document records what was planned, what was actually executed, what was
found, and what remains. It is written so that (a) a new chat session can pick the
work up cold, and (b) the findings can go straight into the report.

Read `README.md` first for how to build and run. This file is the *log*.

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

### 1.3 Runs actually executed

| Case | Group | Parameters | Wall time |
|---|---|---|---|
| `A_Re100_al5` | A (1/9) | Re 100, α 5, A 0.5, dp 0.05 | 286 s |
| `C_dp0.075` | C (complete) | Re 100, α 5, A 0.5, dp 0.075 | 107 s |
| `C_dp0.10` | C (complete) | Re 100, α 5, A 0.5, dp 0.10 | 74 s |
| `wom_al2` / `wom_al5` / `wom_al10` | benchmark (complete) | α 2 / 5 / 10 | — |

**Group C (convergence) is complete** — its three members are dp = 0.05 (the
group-A anchor), 0.075 and 0.10.

Additionally smoke-tested (3 cycles, dp = 0.10, results discarded): stenosis 0.5
and 0.7, branch_ratio 0.5, cfl_scale 0.5, and A = 0.25 / 0.75 / 1.0 / 1.5.

### 1.4 Not executed

- **Group A: 8 of 9 cases.** The (Re × α) factorial — the core of the study.
- **Group B: all 4.** Amplitude sweep.
- **Group D: both.** Asymmetric branch (code smoke-tested, no production run).
- **Group E: all 3.** Stenosis (code smoke-tested, no production run).
- **Group F: 1.** CFL check (smoke-tested).

Run them with `./run_study.sh -j 4`. Estimated ~2–3 h wall time with `-j 4`.

---

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

### 3.3 Verification results (Re 100, α 5, A 0.5; dp 0.05/0.075/0.10)

- **Inlet flow rate** → exact `U·D = 3`: 2.99051, 2.98858, 2.98519. Observed order
  **2.69**, Richardson extrapolation **2.9915**, **GCI 0.04 %**.
- **Flow split** 0.4987–0.4997 vs symmetric prediction 0.5.
- **Mass conservation** `Q₁ = Q₂ + Q₃`: rms closure error 0.78 → 0.69 → 0.53 %.
- **Pulsatility index** 1.012 vs theoretical `2A` = 1.0.
- **OSI** ≈ 0.001 steady vs up to 0.22 pulsatile — correctly vanishes with nothing
  oscillating.
- **Mach** 0.084–0.087; bulk density deviation 1.28 / 1.05 / 0.80 %.

**Honest negative result:** *wall shear did not converge monotonically.* Both the
raw peak (0.208 / 0.311 / 0.188) and the 95th percentile (0.181 / 0.242 / 0.170)
put dp = 0.075 ~35 % above its neighbours, so the anomaly is systematic across the
whole wall, not one bad probe. WSS here is twice-derived (SPH velocity gradient →
interpolated to probes at 1 dp and 2 dp → linearly extrapolated to the wall), and
near-wall WCSPH quantities are known to converge slowly. The analysis reports
`NaN` rather than a fabricated order. **Report TAWSS/OSI as relative comparisons
between cases at the production resolution, not as converged absolute values.**

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

## 5. How to continue (for a fresh session)

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

## 6. Suggested report structure (≈20 pages)

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
