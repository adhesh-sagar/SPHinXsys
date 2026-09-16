#!/usr/bin/env bash
#
# RUN_EVERYTHING.sh - produce every number and figure the report needs, in one go.
#
#   1. geometry check runs for the two extensions (writes VTPs for you to eyeball)
#   2. the Womersley validation benchmark, alpha = 2 / 5 / 10
#   3. the full 21-case parameter sweep
#   4. all figures, result tables and ParaView .pvd collections
#
# It is RESUMABLE and safe to re-run: anything already finished is skipped, so if
# the laptop sleeps or you Ctrl-C, just start it again and it picks up where it
# stopped. (A case is "finished" only once its case_params.json exists, which the
# solver writes last - a half-finished case is correctly re-run.)
#
# USAGE, from the directory containing the built binaries:
#
#     ./RUN_EVERYTHING.sh              # foreground, 4 parallel jobs
#     ./RUN_EVERYTHING.sh -j 2         # fewer jobs if the laptop gets hot
#     nohup ./RUN_EVERYTHING.sh &      # survives closing the terminal
#     tail -f run_everything.log       # watch progress
#
# Expect roughly 2-4 h total with -j 4. The slowest cases are group A at alpha = 2
# (the pulse period there is T = 235, so a single case covers ~1650 time units).

set -u -o pipefail

JOBS=4
while [[ $# -gt 0 ]]; do
    case "$1" in
        -j) JOBS="$2"; shift 2 ;;
        -h|--help) sed -n '2,26p' "$0"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 1 ;;
    esac
done

LOG="run_everything.log"
exec > >(tee -a "$LOG") 2>&1

say() { echo; echo "=============================================================="; \
        echo "  $*"; echo "=============================================================="; }

for b in ./test_2d_T_pipe_pulsatile_parameter_study ./womersley_channel ./run_study.sh; do
    [[ -x "$b" ]] || { echo "error: $b not found. Run this from the build bin/ directory." >&2; exit 1; }
done

START=$(date +%s)
say "START $(date).  Jobs: $JOBS.  Log: $LOG"

#----------------------------------------------------------------------
say "STEP 1/4  Geometry checks for the stenosis and asymmetric branch"
#----------------------------------------------------------------------
# Short runs whose only purpose is to write the wall/fluid VTPs so the geometry can
# be inspected before hours are spent on the full extension sweeps.
geom() {
    local name="$1"; shift
    if [[ -f "output_$name/case_params.json" ]]; then echo "  $name already done"; return 0; fi
    rm -rf ".w_$name"; mkdir -p ".w_$name"
    ( cd ".w_$name" && ../test_2d_T_pipe_pulsatile_parameter_study "$@" --case="$name" \
        > "../log_$name.txt" 2>&1 ) || true
    [[ -d ".w_$name/output_$name" ]] && { rm -rf "output_$name"; mv ".w_$name/output_$name" .; }
    rm -rf ".w_$name"
    echo "  $name done"
}
geom GEOM_stenosis0.7 --Re=100 --alpha=5 --A=0.5 --dp=0.05 --stenosis=0.7 --cycles=3 --analysis_cycles=1
geom GEOM_branch0.5   --Re=100 --alpha=5 --A=0.5 --dp=0.05 --branch_ratio=0.5 --cycles=3 --analysis_cycles=1

cat <<'EOT'

  >>> ACTION FOR YOU (can be done later, does not block the rest):
      Open these in ParaView and confirm the geometry is sane -
          output_GEOM_stenosis0.7/WallBoundary_0000000000.vtp
          output_GEOM_stenosis0.7/WaterBody_0000000000.vtp
          output_GEOM_branch0.5/WallBoundary_0000000000.vtp
      Check: the stenosis throat is open and >= 6 particles across, the wall is
      continuous around it, and the shortened lower branch still reaches its outlet.
      These same screenshots go in the report (see REPORT_GUIDE.md, Figure 2).

EOT

#----------------------------------------------------------------------
say "STEP 2/4  Womersley validation benchmark (alpha = 2, 5, 10)"
#----------------------------------------------------------------------
wom() {
    local a="$1" c="$2"
    if [[ -f "output_wom_al${a}/case_params.json" ]]; then echo "  alpha=$a already done"; return 0; fi
    rm -rf ".w_wom$a"; mkdir -p ".w_wom$a"
    ( cd ".w_wom$a" && ../womersley_channel --alpha="$a" --Re=100 --n_across=40 \
        --cycles="$c" --case="wom_al${a}" > "../log_wom_al${a}.txt" 2>&1 ) || true
    [[ -d ".w_wom$a/output_wom_al${a}" ]] && { rm -rf "output_wom_al${a}"; mv ".w_wom$a/output_wom_al${a}" .; }
    rm -rf ".w_wom$a"
    echo "  alpha=$a done"
}
# cycle counts follow the 0.4*alpha^2 periodicity rule plus the analysis window
wom 2 8
wom 5 14
wom 10 45

#----------------------------------------------------------------------
say "STEP 3/4  Full 21-case parameter sweep (this is the long one)"
#----------------------------------------------------------------------
./run_study.sh -r -j "$JOBS" || echo "  (some cases failed - see manifest.csv and log_<case>.txt)"

#----------------------------------------------------------------------
say "STEP 4/4  Figures, result tables and ParaView collections"
#----------------------------------------------------------------------
python3 analysis/figures.py --root . --out figures

# ParaView .pvd collections: without these the animation is indexed by frame
# number rather than by physical time or cycle phase.
python3 analysis/make_pvd.py

END=$(date +%s)
say "DONE in $(( (END-START)/60 )) min"

echo "  cases completed : $(ls -d output_*/ 2>/dev/null | wc -l | tr -d ' ')"
echo "  figures         : $(find figures -name '*.pdf' 2>/dev/null | wc -l | tr -d ' ')"
echo "  result tables   : $(ls results/*.csv 2>/dev/null | wc -l | tr -d ' ')"
echo
echo "  Sweep status (manifest.csv):"
sed 's/^/    /' manifest.csv 2>/dev/null | head -25
echo
echo "  Next: follow REPORT_GUIDE.md - it maps every figure and table to a report section."
