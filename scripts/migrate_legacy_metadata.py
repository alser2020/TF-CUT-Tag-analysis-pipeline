#!/usr/bin/env python3
"""Convert legacy samples.csv/comparisons.csv into the release TSV format."""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path, PurePosixPath


def join_path(root: str, *parts: str) -> str:
    return str(PurePosixPath(root, *parts))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--comparisons", required=True, type=Path)
    parser.add_argument("--rawdata", required=True, help="Raw-data root written into r1/r2 paths")
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    samples_out = args.outdir / "samples.tsv"
    comparisons_out = args.outdir / "comparisons.tsv"
    for output in (samples_out, comparisons_out):
        if output.exists():
            raise SystemExit(f"Refusing to overwrite: {output}")

    counters: defaultdict[str, int] = defaultdict(int)
    with args.samples.open("r", encoding="utf-8-sig", newline="") as source, samples_out.open(
        "w", encoding="utf-8", newline=""
    ) as target:
        reader = csv.DictReader(source)
        required = {"sample_id", "merge_group"}
        if not required <= set(reader.fieldnames or []):
            raise SystemExit(f"Legacy samples file must contain: {sorted(required)}")
        writer = csv.DictWriter(
            target,
            delimiter="\t",
            lineterminator="\n",
            fieldnames=[
                "sample_id", "group_id", "replicate", "label", "r1", "r2",
                "batch", "pair_id", "description",
            ],
        )
        writer.writeheader()
        for row in reader:
            sample_id = (row.get("sample_id") or "").strip()
            group_id = (row.get("merge_group") or "").strip()
            counters[group_id] += 1
            writer.writerow(
                {
                    "sample_id": sample_id,
                    "group_id": group_id,
                    "replicate": counters[group_id],
                    "label": (row.get("sample_name") or sample_id).strip(),
                    "r1": join_path(args.rawdata, sample_id, f"{sample_id}_1.fq.gz"),
                    "r2": join_path(args.rawdata, sample_id, f"{sample_id}_2.fq.gz"),
                    "batch": "",
                    "pair_id": "",
                    "description": (row.get("description") or "").strip(),
                }
            )

    with args.comparisons.open("r", encoding="utf-8-sig", newline="") as source, comparisons_out.open(
        "w", encoding="utf-8", newline=""
    ) as target:
        reader = csv.DictReader(source)
        required = {"name", "treatment_group", "control_group"}
        if not required <= set(reader.fieldnames or []):
            raise SystemExit(f"Legacy comparisons file must contain: {sorted(required)}")
        writer = csv.DictWriter(
            target,
            delimiter="\t",
            lineterminator="\n",
            fieldnames=["comparison_id", "treatment_group", "control_group", "description"],
        )
        writer.writeheader()
        for row in reader:
            writer.writerow(
                {
                    "comparison_id": (row.get("name") or "").strip(),
                    "treatment_group": (row.get("treatment_group") or "").strip(),
                    "control_group": (row.get("control_group") or "").strip(),
                    "description": "",
                }
            )

    print(samples_out)
    print(comparisons_out)
    print("Replicate numbers were inferred from row order within each legacy merge_group; verify them manually.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
