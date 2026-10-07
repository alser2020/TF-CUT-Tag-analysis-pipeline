#!/usr/bin/env bash
set -euo pipefail
PIPELINE_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
DEST=${1:?"Usage: $0 PROJECT_DIRECTORY"}
mkdir -p "$DEST"
for target in config.yaml samples.tsv comparisons.tsv; do
  if [[ -e "$DEST/$target" ]]; then
    echo "Refusing to overwrite existing file: $DEST/$target" >&2
    exit 1
  fi
done
cp "$PIPELINE_ROOT/config/config.example.yaml" "$DEST/config.yaml"
cp "$PIPELINE_ROOT/config/samples.example.tsv" "$DEST/samples.tsv"
cp "$PIPELINE_ROOT/config/comparisons.example.tsv" "$DEST/comparisons.tsv"
echo "Initialized project: $(cd "$DEST" && pwd)"
echo "Next: edit config.yaml, samples.tsv and comparisons.tsv"
