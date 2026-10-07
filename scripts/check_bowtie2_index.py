#!/usr/bin/env python3
"""Check whether a complete small or large Bowtie2 index exists."""
from __future__ import annotations

import argparse
from pathlib import Path

SMALL = (".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
LARGE = tuple(x + "l" for x in SMALL)


def index_kind(prefix: str | Path) -> str | None:
    p = str(prefix)
    if all(Path(p + suffix).is_file() for suffix in SMALL):
        return "bt2"
    if all(Path(p + suffix).is_file() for suffix in LARGE):
        return "bt2l"
    return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("prefix")
    args = parser.parse_args()
    kind = index_kind(args.prefix)
    if kind:
        print(kind)
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
