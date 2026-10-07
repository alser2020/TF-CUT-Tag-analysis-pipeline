#!/usr/bin/env bash
set -euo pipefail

PIPELINE_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
CONFIG=""
PROFILE="$PIPELINE_ROOT/profiles/slurm"
CHECK_ONLY=false
DRY_RUN=false
LOCAL=false
CORES=4
EXTRA_ARGS=()

usage() {
  cat <<USAGE
Usage: $0 --config PATH [options] [-- extra Snakemake arguments]

Options:
  --config PATH       Project config YAML (required)
  --check-only        Run validation only
  --dry-run           Run Snakemake dry-run after validation
  --profile PATH      SLURM profile directory (default: bundled profile)
  --local             Run locally instead of submitting SLURM jobs
  --cores N           Cores for --local mode (default: 4)
  -h, --help          Show this message
USAGE
}

while (($#)); do
  case "$1" in
    --config)
      CONFIG=${2:?"--config requires a path"}
      shift 2
      ;;
    --check-only)
      CHECK_ONLY=true
      shift
      ;;
    --dry-run)
      DRY_RUN=true
      shift
      ;;
    --profile)
      PROFILE=${2:?"--profile requires a path"}
      shift 2
      ;;
    --local)
      LOCAL=true
      shift
      ;;
    --cores)
      CORES=${2:?"--cores requires an integer"}
      shift 2
      ;;
    --)
      shift
      EXTRA_ARGS+=("$@")
      break
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

[[ -n "$CONFIG" ]] || { echo "--config is required" >&2; usage >&2; exit 2; }
CONFIG=$(python3 -c 'import pathlib,sys; print(pathlib.Path(sys.argv[1]).expanduser().resolve())' "$CONFIG")
[[ -f "$CONFIG" ]] || { echo "Config not found: $CONFIG" >&2; exit 2; }
PROJECT_DIR=$(dirname "$CONFIG")

python3 "$PIPELINE_ROOT/scripts/preflight.py" --config "$CONFIG"
$CHECK_ONLY && exit 0

command -v snakemake >/dev/null 2>&1 || {
  echo "snakemake is not available in PATH; activate the Snakemake 7 environment first." >&2
  exit 127
}

SNAKEMAKE_ARGS=(
  --snakefile "$PIPELINE_ROOT/workflow/Snakefile"
  --configfile "$CONFIG"
  --directory "$PROJECT_DIR"
)

if $LOCAL; then
  SNAKEMAKE_ARGS+=(--cores "$CORES" --use-conda --rerun-incomplete --latency-wait 60 --printshellcmds)
else
  PROFILE=$(python3 -c 'import pathlib,sys; print(pathlib.Path(sys.argv[1]).expanduser().resolve())' "$PROFILE")
  [[ -f "$PROFILE/config.yaml" && -f "$PROFILE/cluster.yaml" && -x "$PROFILE/submit.py" && -x "$PROFILE/status.py" ]] || {
    echo "Invalid SLURM profile: $PROFILE" >&2
    exit 2
  }
  SNAKEMAKE_ARGS+=(
    --profile "$PROFILE"
    --cluster-config "$PROFILE/cluster.yaml"
    --cluster "python $PROFILE/submit.py"
    --cluster-status "python $PROFILE/status.py"
  )
fi

$DRY_RUN && SNAKEMAKE_ARGS+=(--dry-run)
if ((${#EXTRA_ARGS[@]})); then
  SNAKEMAKE_ARGS+=("${EXTRA_ARGS[@]}")
fi
exec snakemake "${SNAKEMAKE_ARGS[@]}"
