#!/usr/bin/env python3
"""Write ParaView .pvd collection files for the VTP series in each output folder.

WHY THIS IS NEEDED
------------------
SPHinXsys writes one VTP per snapshot, named by physical time padded to ten
digits: ``WaterBody_0226205018.vtp`` means t = 226.205018. (This is the same
convention as the test_2d_dambreak benchmark - nothing unusual.)

ParaView will group those into a time series on its own, but it has no way to know
what the numbers MEAN, so the animation is indexed 0, 1, 2, ... instead of by time.
A .pvd file is a small XML index that maps each file to its true timestep, which
gives:

  * a time slider showing real physical time instead of frame numbers
  * correct playback spacing when snapshots are unevenly spaced in time
  * one thing to open instead of picking files out of a list

For the pulsatile study the phase within the cycle matters more than absolute
time, so a SECOND collection is written per body with the time axis expressed in
cycles since the pulsation started, ``<body>_phase.pvd``. Opening that one makes
the time slider read 0.00, 0.25, 0.50 ... in units of T, which is exactly what the
report figures are labelled by.

USAGE
-----
    python3 analysis/make_pvd.py                 # every output_* folder
    python3 analysis/make_pvd.py output_A_Re100_al5
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from xml.sax.saxutils import escape

# WaterBody_0226205018.vtp -> body "WaterBody", stamp 226205018
_PATTERN = re.compile(r"^(?P<body>.+)_(?P<stamp>\d{10})\.vtp$")

# The solver pads physical time to 10 digits after multiplying by this factor.
TIME_SCALE = 1.0e6


def collect(folder: Path) -> dict[str, list[tuple[float, str]]]:
    """Group the VTP files in ``folder`` by body name, sorted by physical time."""
    series: dict[str, list[tuple[float, str]]] = {}
    for f in sorted(folder.glob("*.vtp")):
        m = _PATTERN.match(f.name)
        if not m:
            continue  # e.g. SPHSystemDomainProxy.vtp, which is not a time series
        body = m.group("body")
        t = int(m.group("stamp")) / TIME_SCALE
        series.setdefault(body, []).append((t, f.name))
    for v in series.values():
        v.sort()
    return series


def write_pvd(path: Path, entries: list[tuple[float, str]]) -> None:
    lines = [
        '<?xml version="1.0"?>',
        '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">',
        "  <Collection>",
    ]
    for t, name in entries:
        lines.append(
            f'    <DataSet timestep="{t:.9g}" group="" part="0" '
            f'file="{escape(name)}"/>'
        )
    lines += ["  </Collection>", "</VTKFile>", ""]
    path.write_text("\n".join(lines))


def process(folder: Path) -> int:
    series = collect(folder)
    if not series:
        return 0

    # Cycle phase needs the period and the ramp duration; both are in the metadata.
    T = t_ramp = None
    meta = folder / "case_params.json"
    if meta.is_file():
        try:
            p = json.loads(meta.read_text())
            T = p.get("T")
            t_ramp = p.get("t_ramp", 0.0)
        except (ValueError, OSError):
            pass

    written = 0
    for body, entries in series.items():
        if len(entries) < 2:
            continue  # a single snapshot is not a series
        write_pvd(folder / f"{body}.pvd", entries)
        written += 1
        if T:
            cycles = [((t - t_ramp) / T, n) for t, n in entries]
            write_pvd(folder / f"{body}_phase.pvd", cycles)
            written += 1

    label = f"  {folder.name}: "
    parts = [f"{b} ({len(e)} steps)" for b, e in series.items() if len(e) >= 2]
    print(label + ", ".join(parts) + (" [+ phase]" if T else ""))
    return written


def main(argv: list[str]) -> int:
    targets = [Path(a) for a in argv[1:]] or sorted(Path(".").glob("output_*"))
    targets = [t for t in targets if t.is_dir()]
    if not targets:
        print("no output_* folders found (run this from the build bin/ directory)")
        return 1

    total = sum(process(t) for t in targets)
    print(f"\nwrote {total} .pvd files")
    print("Open the .pvd in ParaView, not the individual .vtp files:")
    print("  <body>.pvd        time slider in physical time")
    print("  <body>_phase.pvd  time slider in cycles since the pulsation started")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
