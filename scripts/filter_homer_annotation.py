#!/usr/bin/env python3
"""Filter HOMER annotatePeaks output by absolute distance to TSS."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def find_column(header: list[str], predicate, description: str) -> int:
    for index, value in enumerate(header):
        if predicate(value.strip().lower()):
            return index
    raise ValueError(f"HOMER annotation 缺少列: {description}; header={header}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotation", required=True, type=Path)
    parser.add_argument("--summits", required=True, type=Path)
    parser.add_argument("--max-distance", required=True, type=int)
    parser.add_argument("--bed", required=True, type=Path)
    parser.add_argument("--summit-bed", required=True, type=Path)
    args = parser.parse_args()

    selected: dict[str, tuple[str, str, str]] = {}
    with args.annotation.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        id_col = 0
        chr_col = find_column(header, lambda x: x in {"chr", "chrom", "chromosome"}, "Chr")
        start_col = find_column(header, lambda x: x == "start", "Start")
        end_col = find_column(header, lambda x: x == "end", "End")
        distance_col = find_column(header, lambda x: "distance to tss" in x, "Distance to TSS")
        for row in reader:
            if len(row) <= max(id_col, chr_col, start_col, end_col, distance_col):
                continue
            try:
                distance = int(float(row[distance_col]))
            except ValueError:
                continue
            if abs(distance) <= args.max_distance:
                selected[row[id_col]] = (row[chr_col], row[start_col], row[end_col])

    args.bed.parent.mkdir(parents=True, exist_ok=True)
    with args.bed.open("w", encoding="utf-8") as bed_handle:
        for peak_id, (chrom, start, end) in selected.items():
            bed_handle.write(f"{chrom}\t{start}\t{end}\t{peak_id}\n")

    with args.summits.open("r", encoding="utf-8") as source, args.summit_bed.open("w", encoding="utf-8") as target:
        for line in source:
            if line.startswith("track") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 4 and fields[3] in selected:
                target.write(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
