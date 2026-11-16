#!/usr/bin/env python3
# Plot the force traces from CALIBRATE_PRESSURE_ADVANCE
#
# Copyright (C) 2026
#
# This file may be distributed under the terms of the GNU GPLv3 license.
import argparse
import collections
import csv

import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("Agg")

SEGMENT_COLORS = {
    "idle": "#dddddd",
    "lead": "#cfe3f7",
    "slow": "#bcd9f2",
    "fast": "#f7d9bc",
    "slow2": "#bcd9f2",
    "out": "#cfe3f7",
}


def read_csv(path):
    passes = collections.OrderedDict()
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            advance = float(row["advance"])
            passes.setdefault(advance, []).append(
                (
                    float(row["time"]),
                    float(row["counts"]),
                    float(row["force_g"]),
                    row["segment"],
                )
            )
    return passes


def segment_spans(samples):
    spans = []
    start = samples[0][0]
    label = samples[0][3]
    for when, _counts, _force, seg in samples:
        if seg != label:
            spans.append((label, start, when))
            start = when
            label = seg
    spans.append((label, start, samples[-1][0]))
    return spans


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", help="CSV written by CALIBRATE_PRESSURE_ADVANCE")
    parser.add_argument("-o", "--output", default="pa_calibrate.png")
    parser.add_argument(
        "--counts",
        action="store_true",
        help="plot raw counts instead of grams",
    )
    args = parser.parse_args()

    passes = read_csv(args.csv)
    if not passes:
        parser.error("no rows in %s" % (args.csv,))
    column = 1 if args.counts else 2
    label = "counts" if args.counts else "force (g)"

    fig, axes = plt.subplots(
        len(passes), 1, sharex=True, figsize=(10, 2.0 * len(passes))
    )
    if len(passes) == 1:
        axes = [axes]
    for axis, (advance, samples) in zip(axes, passes.items()):
        for seg, start, end in segment_spans(samples):
            color = SEGMENT_COLORS.get(seg)
            if color:
                axis.axvspan(start, end, color=color, linewidth=0)
        axis.plot(
            [s[0] for s in samples],
            [s[column] for s in samples],
            color="#204060",
            linewidth=0.8,
        )
        axis.set_ylabel("%.4f" % (advance,), rotation=0, ha="right")
        axis.grid(True, alpha=0.3)
    axes[-1].set_xlabel("time (s)")
    fig.suptitle("Pressure advance sweep - %s" % (label,))
    fig.tight_layout()
    fig.savefig(args.output, dpi=110)
    print("wrote %s" % (args.output,))


if __name__ == "__main__":
    main()
