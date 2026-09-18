#!/usr/bin/env python3
"""Print the ParaView frame-index -> time -> cycle-phase map for a case.

WHY
---
In ParaView, selecting the grouped entry `WaterBody_0000000000.vtp*` loads the
whole snapshot series, but the time slider then reads a plain frame INDEX
(0, 1, 2, ...) with no indication of what physical time or cycle phase each frame
corresponds to. This prints that mapping, so you can jump straight to the frames
you want for a figure.

SPHinXsys names each snapshot by physical time multiplied by 1e6 and zero-padded
to ten digits, so the time is recoverable from the filename alone - nothing else
is needed. Cycle phase additionally uses `T` and `t_ramp` from case_params.json.

USAGE
-----
    python3 analysis/frames.py output_ANIM              # full table
    python3 analysis/frames.py output_ANIM --quarters   # just the montage frames
    python3 analysis/frames.py output_ANIM --body WallProbeObserver
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

_STAMP = re.compile(r"_(\d{10})\.vtp$")


def load_series(folder: Path, body: str):
    """Return [(frame_index, physical_time, filename)] sorted by time."""
    out = []
    for f in folder.glob(f"{body}_*.vtp"):
        m = _STAMP.search(f.name)
        if m:
            out.append((int(m.group(1)) / 1.0e6, f.name))
    out.sort()
    return [(i, t, n) for i, (t, n) in enumerate(out)]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder")
    ap.add_argument("--body", default="WaterBody")
    ap.add_argument("--quarters", action="store_true",
                    help="print only the four quarter-cycle frames of the last full cycle")
    args = ap.parse_args()

    folder = Path(args.folder)
    series = load_series(folder, args.body)
    if not series:
        print(f"no {args.body}_*.vtp files in {folder}")
        return 1

    T = t_ramp = None
    meta = folder / "case_params.json"
    if meta.is_file():
        try:
            p = json.loads(meta.read_text())
            T, t_ramp = p.get("T"), p.get("t_ramp", 0.0)
        except (ValueError, OSError):
            pass

    def phase(t):
        return None if not T else (t - t_ramp) / T

    gaps = [b[1] - a[1] for a, b in zip(series, series[1:])]
    print(f"{folder}/{args.body}: {len(series)} frames, "
          f"t = {series[0][1]:.3f} -> {series[-1][1]:.3f}")
    if gaps:
        print(f"  spacing: mean {sum(gaps)/len(gaps):.3f}, max {max(gaps):.3f}", end="")
        # A max gap far above the mean means snapshots are not continuous - usually
        # a case run without --vtp_all, which animates badly.
        print("   [NOT CONTINUOUS - re-run with --vtp_all=1]"
              if max(gaps) > 5 * (sum(gaps) / len(gaps)) else "   [continuous]")
    print()

    if args.quarters:
        if not T:
            print("cycle phase unavailable (no case_params.json)")
            return 1
        last = int(phase(series[-1][1]))
        if phase(series[-1][1]) - last < 0.75:
            last -= 1
        print(f"  montage frames, cycle {last}:\n")
        print("   panel   frame   time        phase")
        for k in range(4):
            target = last + k * 0.25
            i, t, _ = min(series, key=lambda r: abs(phase(r[1]) - target))
            print(f"   {k+1:5d}   {i:5d}   {t:9.3f}   {phase(t):7.3f}")
        return 0

    print("   frame       time     phase")
    for i, t, _ in series:
        ph = phase(t)
        mark = ""
        if ph is not None and i and abs(ph - round(ph)) < 1e-2:
            mark = f"   <- cycle {int(round(ph))}"
        if i == 0:
            mark = "   <- t=0, fluid at rest"
        print(f"   {i:5d}   {t:9.3f}   " + (f"{ph:7.3f}" if ph is not None else "    -") + mark)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
