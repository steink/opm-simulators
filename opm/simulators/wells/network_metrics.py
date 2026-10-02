#!/usr/bin/env python3
"""Compare flow runs on the metrics that matter for the well/group/network work.

Usage:
    network_metrics.py RUN [RUN ...] [--ref RUN] [--csv]

Each RUN is a directory holding one case's output (CASE.PRT, CASE.DBG and the
summary files), as written by flow with --output-dir=RUN, or a path to the
.PRT/.DBG file itself. Runs are printed as columns. With --ref, cumulative
production is also given as a relative difference to that run. Summary
vectors are read with the `summary` tool from opm-common when it is on PATH.

Most counters come from debug lines in the .DBG file, so they are only as
good as those messages; a zero can mean "did not happen" or "not logged by
this build".
"""

import argparse
import collections
import glob
import os
import re
import shutil
import subprocess
import sys

# (label, pattern) counted over the .DBG file.
DBG_COUNTERS = [
    ("timestep chops", r"Timestep chopped"),
    ("network: solves", r"Network: solved the production network"),
    ("network: give-ups", r"is not possible at report step"),
    ("group-tree: rounds", r"GroupTreeWorkflow: report step \d+ round done"),
    ("group-tree: cliff stops", r"converged past its own tubing-curve cliff"),
    ("group-tree: reopens", r"Network: reopening well"),
    ("group-tree: reopen failed", r"reopened by the network but stopped"),
    ("local: stopped", r"gets STOPPED during iteration"),
    ("local: re-opened", r"is re-opened after being stopped"),
    ("local: switch limit hit", r"is oscillating between open and stop"),
    ("local: not converged", r"did not converge in \d+ inner iterations"),
    ("network: max outer its", r"Maximum of \d+ network iterations"),
]

# Status events per well, for the "busiest wells" line.
WELL_EVENTS = [
    re.compile(r"well (\S+) gets STOPPED during iteration"),
    re.compile(r"(\S+) is re-opened after being stopped"),
    re.compile(r"well (\S+) converged past its own tubing-curve cliff"),
    re.compile(r"Network: reopening well (\S+)"),
]

PRT_FIELDS = [
    ("timesteps", r"Number of timesteps:\s+(\d+)", int),
    ("simulation time [s]", r"Simulation time:\s+([\d.]+) s", float),
    ("well assembly [s]", r"Well assembly:\s+([\d.]+) s", float),
    ("linearizations", r"Overall Linearizations:\s+(\d+)", int),
    ("newton iterations", r"Overall Newton Iterations:\s+(\d+)", int),
    ("wasted newton", r"Overall Newton Iterations:\s+\d+\s+\(Wasted:\s+(\d+)", int),
    ("wasted newton [%]", r"Overall Newton Iterations:\s+\d+\s+\(Wasted:\s+\d+;\s+([\d.]+)%", float),
]

SUMMARY_VECTORS = ["FOPT", "FGPT", "FWPT"]


def case_files(run):
    """Return (name, prt, dbg, summary_base) for a run directory or file."""
    if os.path.isdir(run):
        prts = glob.glob(os.path.join(run, "*.PRT"))
        if len(prts) != 1:
            sys.exit(f"{run}: expected exactly one .PRT file, found {len(prts)}")
        base = prts[0][:-4]
    else:
        base = os.path.splitext(run)[0]
    name = run.rstrip("/")
    return name, base + ".PRT", base + ".DBG", base


def read(path):
    try:
        with open(path, errors="replace") as f:
            return f.read()
    except OSError:
        return ""


def summary_end_values(base):
    """Last values of SUMMARY_VECTORS, or {} if the summary tool is missing."""
    if shutil.which("summary") is None:
        return {}
    try:
        out = subprocess.run(["summary", base] + SUMMARY_VECTORS, capture_output=True,
                             text=True, check=False).stdout
    except OSError:
        return {}
    rows = [line.split() for line in out.splitlines() if line.strip()]
    rows = [r for r in rows if len(r) == len(SUMMARY_VECTORS)]
    if not rows:
        return {}
    try:
        return {k: float(v) for k, v in zip(SUMMARY_VECTORS, rows[-1])}
    except ValueError:
        return {}


def metrics(run):
    name, prt_path, dbg_path, base = case_files(run)
    prt, dbg = read(prt_path), read(dbg_path)
    out = collections.OrderedDict()
    for label, pattern, conv in PRT_FIELDS:
        m = re.findall(pattern, prt)
        out[label] = conv(m[-1]) if m else None
    for label, pattern in DBG_COUNTERS:
        out[label] = len(re.findall(pattern, dbg)) if dbg else None
    per_well = collections.Counter()
    for rx in WELL_EVENTS:
        per_well.update(rx.findall(dbg))
    out["busiest wells (status events)"] = (
        ", ".join(f"{w}:{n}" for w, n in per_well.most_common(3)) if per_well else "-")
    out.update(summary_end_values(base))
    return name, out


def fmt(value):
    if value is None:
        return "-"
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("runs", nargs="+")
    parser.add_argument("--ref", help="run to compare cumulative production against")
    parser.add_argument("--csv", action="store_true", help="comma-separated output")
    args = parser.parse_args()

    runs = [metrics(r) for r in args.runs]
    ref = metrics(args.ref)[1] if args.ref else None
    if ref:
        for _, m in runs:
            for v in SUMMARY_VECTORS:
                if v in m and v in ref and ref[v]:
                    m[f"{v} vs ref [%]"] = 100.0 * (m[v] - ref[v]) / ref[v]

    labels = []
    for _, m in runs:
        labels += [k for k in m if k not in labels]
    names = [os.path.basename(n) or n for n, _ in runs]

    if args.csv:
        print(",".join(["metric"] + names))
        for label in labels:
            print(",".join([label] + [fmt(m.get(label)) for _, m in runs]))
        return

    width = max(len(l) for l in labels)
    colw = [max(len(n), *(len(fmt(m.get(l))) for l in labels)) for n, (_, m) in zip(names, runs)]
    print(" " * width + "  " + "  ".join(n.rjust(w) for n, w in zip(names, colw)))
    for label in labels:
        cells = [fmt(m.get(label)).rjust(w) for (_, m), w in zip(runs, colw)]
        print(label.ljust(width) + "  " + "  ".join(cells))


if __name__ == "__main__":
    main()
