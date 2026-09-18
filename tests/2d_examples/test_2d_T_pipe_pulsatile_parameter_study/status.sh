#!/usr/bin/env bash
#
# status.sh - one-shot answer to "is anything running, and how far along is it?"
#
#     ./status.sh          print status once
#     ./status.sh -w       refresh every 20 s until everything finishes
#
# Progress per case is (current physical time) / (end_time), both scraped from
# that case's log. It is a fair progress bar because the solver advances physical
# time roughly uniformly - but note the time STEP shrinks in hard cases, so the
# last few percent can take disproportionately long.

WATCH=0
[[ "${1:-}" == "-w" ]] && WATCH=1

show() {
    local running done_n
    running=$(pgrep -f "test_2d_T_pipe_pulsatile_parameter_study --" 2>/dev/null | wc -l | tr -d ' ')
    running=$(( running + $(pgrep -f "womersley_channel --" 2>/dev/null | wc -l | tr -d ' ') ))
    done_n=$(ls -d output_*/ 2>/dev/null | wc -l | tr -d ' ')

    echo "======================================================================"
    date "+  %H:%M:%S"
    echo "  completed cases: $done_n        processes computing: $running"
    echo "======================================================================"

    if [[ "$running" == "0" ]]; then
        if pgrep -f RUN_EVERYTHING.sh >/dev/null 2>&1; then
            echo "  RUN_EVERYTHING.sh is alive but no solver is running"
            echo "  (probably between steps, or generating figures)"
        else
            echo "  NOTHING IS RUNNING."
            if [[ -f manifest.csv ]]; then
                echo
                echo "  Last sweep result:"
                sed 's/^/    /' manifest.csv
            else
                echo "  No manifest.csv - the sweep has not completed a full pass."
                echo "  Start or resume it with:   nohup ./RUN_EVERYTHING.sh &"
            fi
        fi
        return
    fi

    printf "  %-18s %9s %9s %7s   %s\n" CASE TIME END PROGRESS STEP
    for f in log_*.txt; do
        [[ -f "$f" ]] || continue
        # only cases whose last line is a progress line are still active
        local last end t pct name
        last=$(tail -1 "$f")
        [[ "$last" == N=* ]] || continue
        name=${f#log_}; name=${name%.txt}
        end=$(grep -m1 "end_time" "$f" | awk '{print $3}')
        t=$(sed -n 's/.*  t=\([0-9.]*\) .*/\1/p' <<< "$last")
        [[ -z "$t" || -z "$end" ]] && continue
        pct=$(awk -v a="$t" -v b="$end" 'BEGIN{printf "%.1f", 100*a/b}')
        printf "  %-18s %9.2f %9.2f %6s%%   %s\n" "$name" "$t" "$end" "$pct" \
            "$(grep -o 'dt=[0-9.e-]*' <<< "$last" | tail -1)"
    done
}

if [[ "$WATCH" == "1" ]]; then
    while true; do
        clear; show
        pgrep -f "RUN_EVERYTHING.sh" >/dev/null 2>&1 || { echo; echo "  pipeline finished."; break; }
        sleep 20
    done
else
    show
fi
