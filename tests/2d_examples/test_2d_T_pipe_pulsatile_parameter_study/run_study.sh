#!/usr/bin/env bash
#
# run_study.sh - drive the full pulsatile T-pipe parameter study.
#
# The study is a designed matrix in the three dimensionless groups (Re, alpha, A)
# plus a resolution series and two geometric extensions. omega is never set here:
# it is derived from alpha inside the solver, because Re and alpha are coupled
# through the kinematic viscosity.
#
#   Group A  full factorial Re x alpha at A = 0.5, dp = 0.05        9 runs
#   Group B  amplitude sweep at Re = 100, alpha = 5                 4 runs
#   Group C  particle-resolution convergence at Re = 100, alpha = 5 3 runs
#   Group D  asymmetric lower branch (resistance network check)     2 runs
#   Group E  stenosis in the upper branch                           3 runs
#   Group F  CFL insensitivity check                                1 run
#                                                                  --------
#                                                                  22 runs
#
# WHY dp = 0.05 IS THE PRODUCTION RESOLUTION
#   The smallest physical length in the problem is the oscillatory (Stokes) layer
#       delta = sqrt(2 nu / omega) = R * sqrt(2) / alpha,
#   which is INDEPENDENT of Re. Resolving it with at least 4 particles requires
#       dp <= 0.53 * R / alpha = 0.795 / alpha,
#   so the largest alpha in the matrix sets dp for the whole factorial:
#       alpha = 10  ->  dp <= 0.0795  ->  dp = 0.05 gives delta/dp = 4.2.
#   Group A must use a SINGLE dp, otherwise resolution is confounded with alpha,
#   which is exactly what the factorial is meant to separate. Group C then shows
#   convergence towards that resolution.
#
# Expected cost: ~4600 units of physical time in total for group A, roughly
# 1 s of wall time per unit at dp = 0.05 on a laptop, i.e. ~1.5 h sequential or
# ~30 min with -j 4. Groups B-F add a comparable amount.
#
# Usage (from the build tree, next to the binary):
#   ./run_study.sh                 run everything sequentially
#   ./run_study.sh -j 4            run 4 cases at a time
#   ./run_study.sh -g A -g C       run only groups A and C
#   ./run_study.sh -n              dry run: print the commands only
#   ./run_study.sh -- --log_level=4   pass extra flags through to the solver
#
# Each case writes ./output_<case>/ containing the .dat time series and
# case_params.json. A manifest.csv summarising all runs is written alongside.

set -u -o pipefail

BIN="./test_2d_T_pipe_pulsatile_parameter_study"
MANIFEST="manifest.csv"
JOBS=1
DRY=0
SEL_GROUPS=""
EXTRA_ARGS=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -j) JOBS="$2"; shift 2 ;;
        -g) SEL_GROUPS="${SEL_GROUPS}$2 "; shift 2 ;;
        -n) DRY=1; shift ;;
        --) shift; EXTRA_ARGS="$*"; break ;;
        -h|--help) sed -n '2,30p' "$0"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 1 ;;
    esac
done

if [[ ! -x "$BIN" ]]; then
    echo "error: $BIN not found or not executable." >&2
    echo "Run this script from the directory containing the built binary." >&2
    exit 1
fi

# NB: do not name this variable GROUPS - that is a bash builtin array.
[[ -z "$SEL_GROUPS" ]] && SEL_GROUPS="A B C D E F"

wanted() { [[ " $SEL_GROUPS " == *" $1 "* ]]; }

#----------------------------------------------------------------------
# Build the case list. One line per case:
#   <group>|<case_name>|<solver arguments>
#----------------------------------------------------------------------
# Production resolution: set by the Stokes layer at the largest alpha (see header).
DP_PROD=0.05
CASES=()

# ---- Group A: full factorial in (Re, alpha) -------------------------
# The core of the study. A and dp are held fixed so that the two factors are
# cleanly separated; alpha = 2 is quasi-steady, alpha = 10 is inertia dominated.
if wanted A; then
    for Re in 50 100 200; do
        for AL in 2 5 10; do
            CASES+=("A|A_Re${Re}_al${AL}|--Re=${Re} --alpha=${AL} --A=0.5 --dp=${DP_PROD}")
        done
    done
fi

# ---- Group B: pulsatility amplitude --------------------------------
# A = 0     steady Poiseuille validation case (no pulsation at all)
# A = 1.0   sectional mean velocity just touches zero at trough - the largest
#           amplitude this inlet treatment supports
#
# A > 1 IS DELIBERATELY EXCLUDED. At A = 1.5 the prescribed inlet velocity becomes
# NEGATIVE over part of the cycle, i.e. fluid would have to flow back out through
# the inlet. The emitter + inflow-buffer combination inherited from T_shaped_pipe is
# unidirectional: it can inject but not accept backflow. Tested at dp = 0.10, A =
# 1.5 collapses to Dt ~ 2.7e-5 and stalls at t/T = 2.4, while A = 0.25, 0.75 and
# 1.0 all run cleanly with Dt ~ 0.015-0.027.
#
# Studying genuine flow reversal needs a bidirectional buffer - see
# tests/extra_source_and_tests/extra_src/shared/pressure_boundary/ - which is a
# separate piece of work, not a parameter change.
if wanted B; then
    for AMP in 0 0.25 0.75 1.0; do
        CASES+=("B|B_A${AMP}|--Re=100 --alpha=5 --A=${AMP} --dp=${DP_PROD}")
    done
fi

# ---- Group C: particle resolution convergence ----------------------
#     dp    = 0.05  ->  0.075  ->  0.10        (0.05 comes from A_Re100_al5)
#     DH/dp =  60   ->   40    ->   30         layers across the channel
#
# The finest member IS the production resolution, so the GCI bounds the runs the
# study actually reports.
#
# TWO constraints drive this choice, and both were learned the hard way:
#
#  1. DH/dp must be an INTEGER, so the particle lattice fits the channel height
#     exactly. With dp = 0.1125 (26.67 layers) and dp = 0.16875 (17.78 layers) the
#     leftover partial layer biases the volume-weighted mean velocity and the inlet
#     flow-rate error jumps to +2.07 % and -2.98 %, with the SIGN flipping between
#     them - which destroys monotonicity. At 60, 40 and 30 layers the same error is
#     -0.32 %, -0.38 %, -0.53 % and well behaved. The solver warns when DH/dp is
#     not an integer.
#
#  2. Every member must be inside the asymptotic range. Going coarser to get a
#     constant ratio (0.05/0.10/0.20) puts dp = 0.20 outside it: the Stokes layer
#     spans only 2.1 particles there, below this study's own >= 4 criterion, the
#     flow-rate error jumps to -5.7 % and the extracted order comes out at 4.6.
#
# Those two constraints leave an UNEQUAL ratio ladder (1.5 then 1.333). That is
# fine: convergence.py solves the ASME V&V 20 / Celik unequal-ratio relation
# rather than assuming a constant ratio. Going finer to force a constant ratio
# (dp = 0.025, 120 layers, 40 800 particles) costs ~9 h for this single case and
# buys nothing.
if wanted C; then
    for DP in 0.075 0.10; do
        CASES+=("C|C_dp${DP}|--Re=100 --alpha=5 --A=0.5 --dp=${DP}")
    done
fi

# ---- Group D: asymmetric branch resistance -------------------------
# Shortening the lower branch lowers its hydraulic resistance, so the flow split
# should follow the analytic resistance-network prediction Q_up/Q_lo = R_lo/R_up.
# branch_ratio = 1.0 is already covered by A_Re100_al5.
if wanted D; then
    for BR in 0.7 0.5; do
        CASES+=("D|D_br${BR}|--Re=100 --alpha=5 --A=0.5 --dp=${DP_PROD} --branch_ratio=${BR}")
    done
fi

# ---- Group E: stenosis in the upper branch -------------------------
# At 70 % the throat is 1.5 * 0.3 = 0.45 wide, i.e. 9 particles at dp = 0.05, so
# the production resolution already resolves the worst case comfortably.
if wanted E; then
    for ST in 0.3 0.5 0.7; do
        CASES+=("E|E_st${ST}|--Re=100 --alpha=5 --A=0.5 --dp=${DP_PROD} --stenosis=${ST}")
    done
fi

# ---- Group F: time-step (CFL) insensitivity ------------------------
# Halving the acoustic time step must move the reported functionals by far less
# than the spatial discretisation error established in group C.
if wanted F; then
    CASES+=("F|F_cfl0.5|--Re=100 --alpha=5 --A=0.5 --dp=${DP_PROD} --cfl_scale=0.5")
fi

if [[ ${#CASES[@]} -eq 0 ]]; then
    echo "no cases selected (groups: $SEL_GROUPS). Valid groups are A B C D E F." >&2
    exit 1
fi
echo "cases selected: ${#CASES[@]}  (groups: $SEL_GROUPS)"

#----------------------------------------------------------------------
# Runner for a single case.
#----------------------------------------------------------------------
run_one() {
    local line="$1"
    local group="${line%%|*}"
    local rest="${line#*|}"
    local name="${rest%%|*}"
    local args="${rest#*|}"
    local log="log_${name}.txt"

    echo "[$group] $name : $args"
    if [[ "$DRY" == "1" ]]; then
        return 0
    fi

    # Each case runs in its own working directory.
    #
    # This is required for -j to work at all. Every SPHinXsys process creates
    # ./output, ./restart, ./reload and ./sphinxsys.log in its working directory,
    # and the solver then renames ./output out of the way. Run several in one
    # directory and they race on those paths: the loser aborts during start-up
    # (Abort trap: 6) before it prints a single time step.
    #
    # The per-case output folder is moved back up afterwards so the analysis still
    # sees a flat set of output_<case>/ directories.
    local work=".work_${name}"
    rm -rf "$work"; mkdir -p "$work"

    local t0 t1 status
    t0=$(date +%s)
    # shellcheck disable=SC2086
    if ( cd "$work" && "../$BIN" $args --case="$name" $EXTRA_ARGS ) > "$log" 2>&1; then
        status="ok"
    else
        status="FAILED"
        echo "  !! $name failed, see $log" >&2
    fi
    t1=$(date +%s)

    if [[ -d "$work/output_$name" ]]; then
        rm -rf "output_$name"
        mv "$work/output_$name" "output_$name"
    fi
    rm -rf "$work"

    printf '%s,%s,%s,%s,%s\n' "$group" "$name" "$status" "$((t1 - t0))" "$args" \
        >> "$MANIFEST.part"
}
export -f run_one
export BIN DRY MANIFEST EXTRA_ARGS

#----------------------------------------------------------------------
# Execute.
#----------------------------------------------------------------------
rm -f "$MANIFEST.part"
: > "$MANIFEST.part"

if [[ "$JOBS" -gt 1 && "$DRY" == "0" ]]; then
    # Each case is single-threaded enough that running several at once is the
    # cheapest way to get through the matrix on a multi-core laptop.
    printf '%s\n' ${CASES[@]+"${CASES[@]}"} | xargs -P "$JOBS" -I{} bash -c 'run_one "$@"' _ {}
else
    for c in ${CASES[@]+"${CASES[@]}"}; do run_one "$c"; done
fi

if [[ "$DRY" == "1" ]]; then
    rm -f "$MANIFEST.part"
    echo "dry run complete, nothing executed."
    exit 0
fi

{
    echo "group,case,status,wall_time_s,args"
    sort "$MANIFEST.part"
} > "$MANIFEST"
rm -f "$MANIFEST.part"

n_ok=$(grep -c ',ok,' "$MANIFEST" || true)
n_bad=$(grep -c ',FAILED,' "$MANIFEST" || true)
echo
echo "wrote $MANIFEST : $n_ok ok, $n_bad failed, of ${#CASES[@]} cases"

if [[ "$n_bad" != "0" ]]; then
    exit 1
fi

echo
echo "next step:"
echo "  python3 analysis/figures.py --root . --out figures"
