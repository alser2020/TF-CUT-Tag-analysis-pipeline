#!/usr/bin/env bash
#SBATCH --job-name=cuttag_driver
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=2-00:00:00
#SBATCH --output=cuttag_driver.%j.out
#SBATCH --error=cuttag_driver.%j.err

set -euo pipefail

resolve_pipeline_root() {
  if [[ -n ${CUTTAG_PIPELINE_ROOT:-} ]]; then
    printf '%s\n' "$CUTTAG_PIPELINE_ROOT"
    return
  fi
  if [[ -n ${SLURM_JOB_ID:-} ]] && command -v scontrol >/dev/null 2>&1; then
    local job_info command_path
    job_info=$(scontrol show job "$SLURM_JOB_ID" -o)
    command_path=${job_info#* Command=}
    command_path=${command_path%% *}
    if [[ -n "$command_path" && -f "$command_path" ]]; then
      cd -- "$(dirname -- "$command_path")" && pwd
      return
    fi
  fi
  cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd
}

PIPELINE_ROOT=$(resolve_pipeline_root)
CONFIG=${1:?"Usage: sbatch [site options] submit.slurm /absolute/path/config.yaml"}

# 在此之前确保 Snakemake 7 已位于 PATH；也可按集群实际情况加载 module/conda。
command -v snakemake >/dev/null 2>&1 || {
  echo "snakemake is not available in PATH" >&2
  exit 127
}
[[ -x "$PIPELINE_ROOT/run.sh" ]] || {
  echo "run.sh not found under resolved pipeline root: $PIPELINE_ROOT" >&2
  exit 1
}

exec "$PIPELINE_ROOT/run.sh" --config "$CONFIG"
