#!/usr/bin/env python3

from pathlib import Path
import csv
import sys

root = Path(sys.argv[1])

depths = ("5x", "10x", "25x", "50x", "100x")

print(
    "seed\tdepth\traw_rows\tlcsv_rows\t"
    "raw_XY\tlcsv_XY\tlion_events"
)

for seed in range(101, 110):
    base = root / f"seed_{seed}/results"

    for depth in depths:
        tag = f"official_simulated_{depth}_noise0.1"

        raw_file = (
            base /
            "LIONS_official_by_depth" /
            depth /
            f"{tag}.LIONS.log"
        )

        out = (
            base /
            "LIONS_simuOptimal_by_depth" /
            depth
        )

        lcsv_file = out / f"{tag}.lcsv"
        lion_file = out / f"{tag}.lion"

        raw = []

        with raw_file.open(errors="replace") as fh:
            for line in fh:
                row = line.rstrip("\n").split("\t")

                if len(row) == 25 and ":" in row[0]:
                    raw.append(row)

        with lcsv_file.open(newline="") as fh:
            lcsv = list(
                csv.DictReader(fh, delimiter="\t")
            )

        with lion_file.open(errors="replace") as fh:
            lion_lines = [
                line for line in fh
                if line.strip()
            ]

        lion_events = max(0, len(lion_lines) - 1)

        raw_xy = sum(
            row[12].removeprefix("chr") in ("X", "Y")
            for row in raw
        )

        lcsv_xy = sum(
            row["Chromosome"].removeprefix("chr") in ("X", "Y")
            for row in lcsv
        )

        if len(raw) != len(lcsv):
            raise SystemExit(
                f"Row mismatch: seed={seed}, depth={depth}, "
                f"raw={len(raw)}, lcsv={len(lcsv)}"
            )

        if raw_xy != lcsv_xy:
            raise SystemExit(
                f"X/Y mismatch: seed={seed}, depth={depth}, "
                f"raw_XY={raw_xy}, lcsv_XY={lcsv_xy}"
            )

        print(
            seed,
            depth,
            len(raw),
            len(lcsv),
            raw_xy,
            lcsv_xy,
            lion_events,
            sep="\t"
        )

print(
    "\nPASS: all 45 LCSV files retain all "
    "raw interactions, including X/Y."
)
