# Pulsatile T-pipe parameter study

A dimensionless parameter study of pulsatile flow through a 2D T-shaped
bifurcation, built on the SPHinXsys `T_shaped_pipe` benchmark, plus a dedicated
Womersley benchmark that validates the unsteady viscous scheme.

**Your original single-case file is untouched.** It still lives in
`../test_2d_T_pipe_pulsatile_poiseuille_study/` and still builds and runs. This
folder is entirely additional.

---

## What is here

| File | What it is |
|---|---|
| `T_pipe_pulsatile_parameter_study.cpp` | the study solver (one run per invocation) |
| `case_params.h` | case definition, CLI parsing, metadata output |
| `womersley_channel_validation.cpp` | the separate Womersley benchmark (see below) |
| `run_study.sh` | drives the full 22-run matrix |
| `analysis/` | Python post-processing; produces every report figure and table |

Two executables are built: `test_2d_T_pipe_pulsatile_parameter_study` and
`womersley_channel`.

## Build and run

```bash
cd <build dir> && cmake . && ninja test_2d_T_pipe_pulsatile_parameter_study \
                                  womersley_channel \
                                  test_2d_T_pipe_pulsatile_parameter_study_tools
cd tests/2d_examples/test_2d_T_pipe_pulsatile_parameter_study/bin

./test_2d_T_pipe_pulsatile_parameter_study --help_case   # all options
./run_study.sh -n                                        # dry run: print the matrix
./run_study.sh -j 4                                      # the whole study
./womersley_channel --alpha=5 --cycles=14 --case=wom_al5 # one benchmark run

python3 -m pip install -r analysis/requirements.txt
python3 analysis/figures.py --root . --out figures
```

Each run writes `output_<case>/` containing `.dat` time series and
`case_params.json`. Figures land in `figures/` (PDF + PNG), and every number that
appears in a figure is also written to `results/*.csv` so the report's tables are
traceable to the run that produced them.

---

## The study design in one page

The case is defined by **three dimensionless groups and nothing else**:

| Group | Definition | Values |
|---|---|---|
| `Re` | ρ U D / μ | 50, 100, 200 |
| `α` (Womersley) | R √(ω/ν) | 2, 5, 10 |
| `A` (pulsatility) | ΔU / Ū | 0, 0.25, 0.5, 1.0, 1.5 |

`Re` and `α` are **coupled through ν**, so `ω` is never set directly. It is
*derived*: `ν = U D / Re`, then `ω = ν (α/R)²`. Setting `ω` by hand — as the
original file did — means each run lands on whatever `α` happens to result, which
is why the α values in the original header comment did not match its own formula.

### Two derived quantities that drive everything

**Run length.** The slowest start-up mode of a plane channel decays at
`λ₁ = ν(π/2R)² = (π²/4)/τ_visc`, and `τ_visc/T = α²/2π`, so

```
cycles to periodicity ≈ 0.4 α²      →   α=2: 2 cycles,  α=10: 40 cycles
```

A fixed cycle count cannot work across the matrix. `--cycles=0` (the default)
derives it per case.

**Resolution.** The smallest physical length is the Stokes layer
`δ = √(2ν/ω) = R√2/α`, which is **independent of Re**. Resolving it with ≥ 4
particles needs `dp ≤ 0.795/α`, so the largest α in the matrix sets `dp` for the
whole factorial: `dp = 0.05` gives `δ/dp = 4.2` at α = 10. Group A must use a
single `dp`, otherwise resolution is confounded with α — which is exactly what the
factorial exists to separate.

The benchmark confirms this criterion independently: its L2 error is 0.96 % at
δ/dp = 14, 1.63 % at 5.7, and 3.25 % at 2.8.

### Choosing resolutions for the convergence study

Two constraints must hold **simultaneously**, and violating either silently
destroys the convergence result:

1. **Constant refinement ratio.** Roache's Richardson/GCI procedure assumes one.
   A mixed series (0.075 / 0.10 / 0.15, ratios 1.33 and 1.5) produced a *negative*
   observed order for the mean flow rate — meaningless — even though the values
   were converging monotonically.

2. **`DH/dp` must be an integer.** The lattice generator fills the channel with
   whole particle layers; a leftover partial layer biases every volume-weighted
   average. At dp = 0.1125 (26.67 layers) and 0.16875 (17.78 layers) the inlet
   flow-rate error jumped to +2.07 % and −2.98 %, *flipping sign* between them. At
   60, 40 and 30 layers the same error is −0.32 %, −0.38 %, −0.53 %. The solver
   now prints a warning when `DH/dp` is not an integer.

Group C therefore uses **dp = 0.025 / 0.05 / 0.10** (120 / 60 / 30 layers,
constant ratio 2), with the production resolution deliberately in the **middle**
so the GCI computed on the finest pair bounds the runs the study actually reports.
Going one step coarser instead (0.05 / 0.10 / 0.20) puts the coarsest member
outside the asymptotic range — δ/dp = 2.1 there, below this study's own criterion —
and the extracted order comes out at 4.6, which is not believable for a
second-order scheme.

---

## Why there are two executables

**The T-pipe cannot validate the Womersley solution.** Its inlet profile is
*prescribed* as a parabola by the inflow buffer, and the inlet channel is 3.5 long
while the viscous entrance length is `0.04 Re D = 12` and the oscillatory
development length is `U/ω = 3.0`. The flow reaches the junction still carrying the
imposed profile — measured peakedness is 1.41, against 1.50 for a pure parabola and
1.21 for the Womersley prediction at α = 5.

So `womersley_channel_validation.cpp` does the validation properly: a straight
channel, periodic in x, **no velocity boundary condition anywhere**, driven only by
an oscillatory body force `G cos(ωt)` — precisely the configuration Womersley's
solution describes. Measured against theory:

| α | δ/dp | L2 | centre \|û\| sim / theory | centre phase sim / theory | annular peak sim / theory |
|---|---|---|---|---|---|
| 2  | 14.1 | 0.96 % | 1.4846 / 1.4830 | −61.1° / −60.9° | 0.00 / 0.00 |
| 5  | 5.7  | 1.63 % | 1.2093 / 1.2114 | −91.4° / −91.2° | 0.34 / 0.34 |
| 10 | 2.8  | 3.25 % | 1.0714 / 1.0717 | −90.8° / −89.9° | 0.68 / 0.68 |

Core flattening, the phase approaching the 90° limit, and the Richardson annular
effect are all reproduced — the annular peak position matches exactly.

Use the `womersley_channel_*` figures for the report's validation section. The
`tpipe_*_vs_theory` figures are labelled as profile *relaxation*, not validation.

---

## Results that check out

Measured at Re = 100, α = 5, A = 0.5 on the dp = 0.05 / 0.075 / 0.10 series
(60 / 40 / 30 layers). See `results/*.csv` for the full tables.

- **Inlet flow rate** converges monotonically to the exact value `U·D = 3`:
  2.99051, 2.98858, 2.98519 (errors −0.32 %, −0.38 %, −0.49 %). Observed order of
  convergence **2.69**, Richardson extrapolation **2.9915** against an exact 3,
  **GCI 0.04 %**.
- **Flow split** 0.4987–0.4997 against the symmetric prediction of 0.5.
- **Mass conservation** `Q₁ = Q₂ + Q₃`: rms closure error 0.78 → 0.69 → 0.53 % as
  dp refines; the cycle-mean error is smaller still.
- **Pulsatility index** 1.012 at the inlet, against the theoretical `2A = 1.0`.
- **OSI** ≈ 0.001 for steady flow and up to 0.22 for pulsatile — the oscillatory
  shear index correctly vanishes when there is nothing oscillating.
- **Mach number** 0.084–0.087 (< 0.1); bulk density deviation 1.28 % at dp = 0.05,
  1.05 % at 0.075, 0.80 % at 0.10.
- **Womersley benchmark** reproduces theory to 0.96 / 1.63 / 3.25 % in L2 at
  α = 2 / 5 / 10, with centre amplitudes within 0.2 %, phases within 1°, and the
  Richardson annular peak position matching exactly.

## Things worth knowing before you trust a number

- **The analysis reports `NaN` rather than a fabricated convergence order** when a
  series is not in the asymptotic range. `observed_order` rejects a non-positive
  order (coarse-grid differences smaller than fine-grid ones) and anything above 4
  (not believable for a second-order scheme). A `NaN` in `convergence_gci.csv`
  means "this series does not support Richardson extrapolation", which is a result
  worth reporting, not a bug to work around.
- **Wall shear did NOT converge monotonically on the tested series**, and this is
  a genuine verification finding rather than a tooling problem. Both the raw peak
  (0.208 / 0.311 / 0.188 at dp = 0.05 / 0.075 / 0.10) and the more robust 95th
  percentile (0.181 / 0.242 / 0.170) put dp = 0.075 about 35 % above its
  neighbours, so the anomaly is systematic across the whole wall, not one noisy
  probe. WSS here is a *twice-derived* quantity — an SPH velocity gradient,
  interpolated to probes at 1 dp and 2 dp, then linearly extrapolated to the wall —
  and near-wall quantities in WCSPH are well known to converge slowly and
  irregularly. Report TAWSS/OSI as *relative* comparisons between cases at the
  production resolution, not as converged absolute values, and say so explicitly.
  If you need a converged absolute WSS, that needs a finer series and probably a
  kernel-corrected gradient.
- **`flow_split` reports a NaN order** because it is already converged to ~0.1 % at
  every resolution, so its differences are noise; the implied order (~4) is
  rejected as not believable for a second-order scheme.
- **`eps_dp` (pressure-drop periodicity) is 0.3–0.9**, far worse than the flow-rate
  and kinetic-energy criteria (both < 2 %). Point pressure in WCSPH is noise
  dominated. Use `eps_Q` and `eps_KE` as the periodicity criteria and report
  `eps_dp` with that caveat.
- **`flow_split` GCI is `NaN`** because the quantity is already converged to
  ~5·10⁻⁴ of 0.5 at every resolution, so the differences are noise. That is the
  correct behaviour, not a failure.
- Wall probes falling inside the stenosis are dropped automatically (no fluid
  support); they show up as `NaN`, not as zeros.

## Problems that were found and fixed

Recorded here because both are easy to reintroduce.

**The run was flaky.** With the baseline's `Integration2ndHalfWithWallNoRiemann`,
identical runs either finished cleanly or collapsed to `dt ~ 5e-5` and exhausted
the particle buffer, depending only on thread scheduling — 1 in 3 failed. The
baseline is validated only over a short steady run (t ≈ 33); this study runs to
t ≈ 150 with a pulsating inflow, i.e. an order of magnitude more exposure to the
sharp re-entrant corners at the junction. Switching the continuity step to
`Integration2ndHalfWithWallRiemann` fixed it: 3/3 stable, and the rogue particles
disappeared (max speed 2.37 instead of 70–418).

**Parallel runs need separate working directories.** Every SPHinXsys process
creates `./output`, `./restart`, `./reload` and `./sphinxsys.log` in its working
directory, and the solver then renames `./output` out of the way. Run several in
one directory and they race on those paths — the loser aborts during start-up
(`Abort trap: 6`) before printing a single time step. `run_study.sh` gives each
case its own `.work_<case>/` directory and moves the results back afterwards.

**`GROUPS` is a bash builtin array.** The sweep driver originally used it for the
`-g` selection; assigning to it as a scalar silently does nothing, so every group
test failed and the matrix came out empty. It is now `SEL_GROUPS`.

**Do not size `c_f` from `MaximumSpeed`.** The library's `MaximumSpeed` reduces
over *every* particle, including those in the inflow buffer and the outflow
disposers, where it reports 70–418 in a flow of order 1. Sizing the sound speed
against that inflates it 2× and measurably *destabilises* the run. `MaxBulkSpeed`
in the solver restricts to interior particles and is the number to use.
