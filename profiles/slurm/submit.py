#!/usr/bin/env python3
"""Submit one Snakemake 7 legacy-cluster job to SLURM."""
from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

from snakemake.utils import read_job_properties


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)[:120]


def main() -> int:
    if len(sys.argv) != 2:
        print(f"usage: {sys.argv[0]} JOBSCRIPT", file=sys.stderr)
        return 2
    jobscript = Path(sys.argv[1]).resolve()
    properties = read_job_properties(str(jobscript))
    cluster = properties.get("cluster", {}) or {}
    rule = properties.get("rule", "job")
    jobid = properties.get("jobid", "na")
    wildcards = properties.get("wildcards", {}) or {}
    wildcard_text = ".".join(f"{key}-{value}" for key, value in sorted(wildcards.items()))
    name = safe_name(".".join(part for part in ("ct", rule, wildcard_text, str(jobid)) if part))

    log_dir = Path(cluster.get("log_dir") or "slurm_logs")
    if not log_dir.is_absolute():
        log_dir = Path.cwd() / log_dir
    log_dir.mkdir(parents=True, exist_ok=True)

    command = [
        "sbatch",
        "--parsable",
        f"--job-name={name}",
        f"--cpus-per-task={max(1, int(properties.get('threads', 1)))}",
        f"--output={log_dir}/%x.%j.out",
        f"--error={log_dir}/%x.%j.err",
    ]
    optional = {
        "account": "--account",
        "partition": "--partition",
        "time": "--time",
    }
    for key, option in optional.items():
        value = cluster.get(key)
        if value not in (None, "", "null"):
            command.append(f"{option}={value}")
    mem_mb = cluster.get("mem_mb")
    if mem_mb not in (None, "", "null"):
        command.append(f"--mem={int(mem_mb)}M")
    command.append(str(jobscript))

    result = subprocess.run(command, check=True, text=True, capture_output=True)
    scheduler_id = result.stdout.strip().split(";", 1)[0]
    if not scheduler_id:
        print(result.stderr, file=sys.stderr)
        return 1
    print(scheduler_id)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
