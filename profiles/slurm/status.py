#!/usr/bin/env python3
"""Report a SLURM job state to Snakemake legacy cluster executor."""
from __future__ import annotations

import subprocess
import sys

RUNNING = {
    "PENDING", "RUNNING", "CONFIGURING", "COMPLETING", "SUSPENDED", "RESIZING",
    "REQUEUED", "REQUEUE_FED", "REQUEUE_HOLD", "SIGNALING", "STAGE_OUT",
}
SUCCESS = {"COMPLETED"}
FAILED = {
    "BOOT_FAIL", "CANCELLED", "DEADLINE", "FAILED", "NODE_FAIL", "OUT_OF_MEMORY",
    "PREEMPTED", "REVOKED", "SPECIAL_EXIT", "STOPPED", "TIMEOUT",
}


def command_output(command: list[str]) -> str:
    try:
        result = subprocess.run(command, check=False, text=True, capture_output=True)
    except OSError:
        return ""
    return result.stdout.strip()


def normalize(value: str) -> str:
    # sacct can return values such as CANCELLED by 12345 or COMPLETED+.
    return value.strip().split()[0].rstrip("+").split("|", 1)[0]


def main() -> int:
    if len(sys.argv) != 2:
        print("running")
        return 0
    jobid = sys.argv[1]
    output = command_output(
        ["sacct", "-X", "-j", jobid, "--format=State", "--noheader", "--parsable2"]
    )
    states = [normalize(line) for line in output.splitlines() if line.strip()]
    state = states[0] if states else ""
    if not state:
        state = normalize(command_output(["squeue", "-h", "-j", jobid, "-o", "%T"]))
    if state in SUCCESS:
        print("success")
    elif state in FAILED:
        print("failed")
    else:
        # Includes accounting delay and all active states.
        print("running")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
