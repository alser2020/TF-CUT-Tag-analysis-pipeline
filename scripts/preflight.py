#!/usr/bin/env python3
"""Validate CUT&Tag configuration and cross-file metadata before submission."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

import yaml

ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
SAMPLE_REQUIRED = ("sample_id", "group_id", "replicate", "r1", "r2")
COMPARISON_REQUIRED = ("comparison_id", "treatment_group", "control_group")
SMALL_INDEX = (".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
LARGE_INDEX = tuple(x + "l" for x in SMALL_INDEX)


class Report:
    def __init__(self) -> None:
        self.errors: list[str] = []
        self.warnings: list[str] = []
        self.info: list[str] = []

    def error(self, message: str) -> None:
        self.errors.append(message)

    def warn(self, message: str) -> None:
        self.warnings.append(message)

    def note(self, message: str) -> None:
        self.info.append(message)

    def emit(self) -> None:
        for message in self.info:
            print(f"[INFO] {message}")
        for message in self.warnings:
            print(f"[WARN] {message}", file=sys.stderr)
        for message in self.errors:
            print(f"[ERROR] {message}", file=sys.stderr)
        print(
            f"Preflight summary: {len(self.errors)} error(s), "
            f"{len(self.warnings)} warning(s)."
        )


def nested(data: dict[str, Any], keys: Iterable[str], report: Report) -> Any:
    current: Any = data
    path: list[str] = []
    for key in keys:
        path.append(key)
        if not isinstance(current, dict) or key not in current:
            report.error(f"config 缺少必需键: {'.'.join(path)}")
            return None
        current = current[key]
    return current


def resolve_path(value: str, root: Path) -> Path:
    expanded = os.path.expandvars(os.path.expanduser(str(value)))
    path = Path(expanded)
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def read_tsv(path: Path, required: tuple[str, ...], label: str, report: Report) -> list[dict[str, str]]:
    if not path.is_file():
        report.error(f"{label} 文件不存在: {path}")
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        missing = [column for column in required if column not in fields]
        if missing:
            report.error(f"{label} 缺少列: {', '.join(missing)}")
            return []
        rows = []
        for line_no, raw in enumerate(reader, start=2):
            row = {key: (value or "").strip() for key, value in raw.items()}
            row["__line__"] = str(line_no)
            rows.append(row)
    if not rows:
        report.error(f"{label} 没有数据行: {path}")
    return rows


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def fasta_contigs(path: Path) -> set[str]:
    contigs: set[str] = set()
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].strip().split()[0]
                if name:
                    contigs.add(name)
    return contigs


def gtf_contigs(path: Path) -> set[str]:
    contigs: set[str] = set()
    with open_text(path) as handle:
        for line in handle:
            if line and not line.startswith("#"):
                fields = line.rstrip("\n").split("\t", 1)
                if fields and fields[0]:
                    contigs.add(fields[0])
    return contigs


def index_exists(prefix: Path) -> bool:
    value = str(prefix)
    return all(Path(value + suffix).is_file() for suffix in SMALL_INDEX) or all(
        Path(value + suffix).is_file() for suffix in LARGE_INDEX
    )


def inspect_index_contigs(prefix: Path, report: Report) -> set[str] | None:
    executable = shutil.which("bowtie2-inspect")
    if not executable:
        report.warn("未找到 bowtie2-inspect，跳过 Bowtie2 index contig 一致性检查")
        return None
    try:
        result = subprocess.run(
            [executable, "-n", str(prefix)],
            check=True,
            text=True,
            capture_output=True,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        report.warn(f"bowtie2-inspect 检查失败，未验证 index contig: {exc}")
        return None
    return {line.strip() for line in result.stdout.splitlines() if line.strip()}


def validate(config_path: Path, skip_file_checks: bool = False) -> tuple[Report, dict[str, Any]]:
    report = Report()
    root = config_path.parent.resolve()
    try:
        config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    except (OSError, yaml.YAMLError) as exc:
        report.error(f"无法读取 YAML: {exc}")
        return report, {}
    if not isinstance(config, dict):
        report.error("config 顶层必须是 mapping")
        return report, {}

    schema_path = Path(__file__).resolve().parents[1] / "schemas" / "config.schema.yaml"
    try:
        import jsonschema
    except ImportError:
        report.warn("Python 环境缺少 jsonschema，跳过 config schema 校验")
    else:
        schema = yaml.safe_load(schema_path.read_text(encoding="utf-8"))
        validator = jsonschema.Draft202012Validator(schema)
        for error in sorted(
            validator.iter_errors(config),
            key=lambda item: tuple(str(part) for part in item.absolute_path),
        ):
            location = ".".join(str(part) for part in error.absolute_path) or "<root>"
            report.error(f"config schema [{location}]: {error.message}")

    output_value = nested(config, ("project", "output_dir"), report)
    sample_value = nested(config, ("project", "samples"), report)
    comparison_value = nested(config, ("project", "comparisons"), report)
    assembly = nested(config, ("reference", "assembly"), report)
    fasta_value = nested(config, ("reference", "fasta"), report)
    gtf_value = nested(config, ("reference", "gtf"), report)
    genome_size = nested(config, ("reference", "effective_genome_size"), report)
    chromosomes = nested(config, ("reference", "chromosomes"), report)

    output_dir = resolve_path(output_value, root) if output_value else None
    samples_path = resolve_path(sample_value, root) if sample_value else None
    comparisons_path = resolve_path(comparison_value, root) if comparison_value else None
    fasta_path = resolve_path(fasta_value, root) if fasta_value else None
    gtf_path = resolve_path(gtf_value, root) if gtf_value else None

    named_paths = {
        "project root": root,
        "output_dir": output_dir,
        "samples": samples_path,
        "comparisons": comparisons_path,
        "FASTA": fasta_path,
        "GTF": gtf_path,
    }
    for label, path in named_paths.items():
        if path is not None and any(character.isspace() for character in str(path)):
            report.error(f"当前版本不支持路径中含空白字符 ({label}): {path}")

    if not isinstance(assembly, str) or not assembly.strip():
        report.error("reference.assembly 必须是非空字符串")
    if not isinstance(genome_size, int) or isinstance(genome_size, bool) or genome_size <= 0:
        report.error("reference.effective_genome_size 必须是正整数")
    if not isinstance(chromosomes, list) or not chromosomes or not all(isinstance(x, str) and x for x in chromosomes):
        report.error("reference.chromosomes 必须是非空字符串列表")

    samples = read_tsv(samples_path, SAMPLE_REQUIRED, "samples", report) if samples_path else []
    comparisons = read_tsv(comparisons_path, COMPARISON_REQUIRED, "comparisons", report) if comparisons_path else []

    sample_ids: set[str] = set()
    group_replicates: set[tuple[str, int]] = set()
    groups: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    fastq_paths: list[Path] = []
    for row in samples:
        line = row["__line__"]
        sample_id = row["sample_id"]
        group_id = row["group_id"]
        if not ID_RE.fullmatch(sample_id):
            report.error(f"samples 第 {line} 行 sample_id 非法: {sample_id!r}")
        if not ID_RE.fullmatch(group_id):
            report.error(f"samples 第 {line} 行 group_id 非法: {group_id!r}")
        if sample_id in sample_ids:
            report.error(f"samples sample_id 重复: {sample_id}")
        sample_ids.add(sample_id)
        try:
            replicate = int(row["replicate"])
            if replicate <= 0:
                raise ValueError
        except ValueError:
            report.error(f"samples 第 {line} 行 replicate 必须是正整数: {row['replicate']!r}")
            continue
        key = (group_id, replicate)
        if key in group_replicates:
            report.error(f"samples 中 (group_id, replicate) 重复: {key}")
        group_replicates.add(key)
        groups[group_id].append(row)
        r1 = resolve_path(row["r1"], root)
        r2 = resolve_path(row["r2"], root)
        if r1 == r2:
            report.error(f"samples 第 {line} 行 R1 与 R2 指向同一文件")
        fastq_paths.extend((r1, r2))
        for read_name, path in (("R1", r1), ("R2", r2)):
            if any(character.isspace() for character in str(path)):
                report.error(f"samples 第 {line} 行 {read_name} 路径含空白字符: {path}")
        if not skip_file_checks:
            for read_name, path in (("R1", r1), ("R2", r2)):
                if not path.is_file():
                    report.error(f"samples 第 {line} 行 {read_name} 不存在: {path}")
        if output_dir and (output_dir == r1 or output_dir == r2 or output_dir in r1.parents or output_dir in r2.parents):
            report.error(f"FASTQ 不应位于 output_dir 内: {r1 if output_dir in r1.parents else r2}")

    comparison_ids: set[str] = set()
    idr_mode = nested(config, ("analysis", "idr", "mode"), report)
    if idr_mode not in {"auto", "off"}:
        report.error("analysis.idr.mode 只能是 auto 或 off")
    for row in comparisons:
        line = row["__line__"]
        comp = row["comparison_id"]
        treatment = row["treatment_group"]
        control = row["control_group"]
        if not ID_RE.fullmatch(comp):
            report.error(f"comparisons 第 {line} 行 comparison_id 非法: {comp!r}")
        if comp in comparison_ids:
            report.error(f"comparisons comparison_id 重复: {comp}")
        comparison_ids.add(comp)
        if treatment not in groups:
            report.error(f"comparison {comp} 的 treatment_group 不存在: {treatment}")
        if control not in groups:
            report.error(f"comparison {comp} 的 control_group 不存在: {control}")
        if treatment == control:
            report.error(f"comparison {comp} 的 treatment_group 与 control_group 相同")
        n_reps = len(groups.get(treatment, []))
        if idr_mode == "auto":
            if n_reps == 2:
                report.note(f"comparison {comp}: treatment 有 2 个重复，将运行 IDR")
            elif n_reps < 2:
                report.warn(f"comparison {comp}: treatment 只有 {n_reps} 个重复，跳过 IDR")
            else:
                report.warn(f"comparison {comp}: treatment 有 {n_reps} 个重复，当前版本不自动截取，跳过 IDR")

    macs = nested(config, ("analysis", "macs3"), report)
    if isinstance(macs, dict):
        qvalue = macs.get("qvalue")
        if not isinstance(qvalue, (int, float)) or isinstance(qvalue, bool) or not 0 < qvalue <= 1:
            report.error("analysis.macs3.qvalue 必须在 (0, 1] 范围内")
        if macs.get("broad") is not False:
            report.error("当前发行版输出 narrowPeak/summits，因此 analysis.macs3.broad 必须为 false")
        if macs.get("call_summits") is not True:
            report.error("motif/annotation 依赖 summits，因此 analysis.macs3.call_summits 必须为 true")
        keep_dup = macs.get("keep_dup")
        if keep_dup != "auto" and (
            not isinstance(keep_dup, int) or isinstance(keep_dup, bool) or keep_dup <= 0
        ):
            report.error("analysis.macs3.keep_dup 必须是 auto 或正整数")

    alignment = nested(config, ("analysis", "alignment"), report)
    if isinstance(alignment, dict):
        try:
            insert_min = int(alignment.get("insert_min"))
            insert_max = int(alignment.get("insert_max"))
            if insert_min < 0 or insert_max <= insert_min:
                raise ValueError
        except (TypeError, ValueError):
            report.error("alignment 要求 0 <= insert_min < insert_max")
        mapq = alignment.get("mapq_min")
        if not isinstance(mapq, int) or isinstance(mapq, bool) or mapq < 0:
            report.error("analysis.alignment.mapq_min 必须是非负整数")

    annotation = nested(config, ("analysis", "annotation"), report)
    motif = nested(config, ("analysis", "motif"), report)
    if isinstance(annotation, dict) and annotation.get("method") != "homer":
        report.error("当前 annotation.method 仅支持 homer")
    if isinstance(motif, dict):
        methods = motif.get("methods")
        if not isinstance(methods, list) or not methods or not set(methods) <= {"homer", "memechip"}:
            report.error("analysis.motif.methods 必须是 homer/memechip 的非空列表")
        if "memechip" in (methods if isinstance(methods, list) else []):
            report.note("MEME-ChIP 使用 treatment peak summit FASTA；当前不使用 external negative/control FASTA")
        scope = motif.get("input_scope")
        if scope not in {"all", "promoter"}:
            report.error("analysis.motif.input_scope 只能是 all 或 promoter")
        if scope == "promoter" and not (isinstance(annotation, dict) and annotation.get("enabled") is True):
            report.error("motif.input_scope=promoter 时必须启用 annotation")
        memechip = motif.get("memechip")
        if isinstance(memechip, dict):
            min_width = memechip.get("min_width")
            max_width = memechip.get("max_width")
            if isinstance(min_width, int) and isinstance(max_width, int) and min_width > max_width:
                report.error("analysis.motif.memechip 要求 min_width <= max_width")
            for key in ("flank", "center_cut", "meme_p", "meme_nmotifs", "streme_nmotifs"):
                value = memechip.get(key)
                if not isinstance(value, int) or isinstance(value, bool):
                    report.error(f"analysis.motif.memechip.{key} 必须是整数")
            if isinstance(memechip.get("flank"), int) and memechip["flank"] < 1:
                report.error("analysis.motif.memechip.flank 必须 >= 1")
            if isinstance(memechip.get("center_cut"), int) and memechip["center_cut"] < 0:
                report.error("analysis.motif.memechip.center_cut 必须 >= 0")
            if isinstance(memechip.get("meme_p"), int) and memechip["meme_p"] < 1:
                report.error("analysis.motif.memechip.meme_p 必须 >= 1")
            if isinstance(memechip.get("meme_nmotifs"), int) and memechip["meme_nmotifs"] < 1:
                report.error("analysis.motif.memechip.meme_nmotifs 必须 >= 1")
            if isinstance(memechip.get("streme_nmotifs"), int):
                if memechip["streme_nmotifs"] < 0:
                    report.error("analysis.motif.memechip.streme_nmotifs 必须 >= 0")
                elif memechip["streme_nmotifs"] > 0:
                    report.warn("memechip.streme_nmotifs > 0 将启用 MEME-ChIP 内部 STREME；当前推荐设置为 0")
            db = memechip.get("known_motif_db")
            if db is not None and not isinstance(db, str):
                report.error("analysis.motif.memechip.known_motif_db 必须是路径字符串或 null")
        else:
            report.error("analysis.motif.memechip 配置缺失")


    if not skip_file_checks:
        for label, path in (("FASTA", fasta_path), ("GTF", gtf_path)):
            if path and not path.is_file():
                report.error(f"reference {label} 不存在: {path}")
        if fasta_path and fasta_path.is_file() and gtf_path and gtf_path.is_file():
            fasta_names = fasta_contigs(fasta_path)
            gtf_names = gtf_contigs(gtf_path)
            if not fasta_names:
                report.error(f"FASTA 未解析到 contig: {fasta_path}")
            missing_gtf = sorted(gtf_names - fasta_names)
            if missing_gtf:
                preview = ", ".join(missing_gtf[:10])
                report.error(f"GTF 中有 {len(missing_gtf)} 个 seqname 不在 FASTA 中: {preview}")
            missing_chr = sorted(set(chromosomes or []) - fasta_names)
            if missing_chr:
                report.error(f"reference.chromosomes 不在 FASTA 中: {', '.join(missing_chr)}")

            index_value = config.get("reference", {}).get("bowtie2_index")
            if index_value:
                prefix = resolve_path(index_value, root)
                if any(character.isspace() for character in str(prefix)):
                    report.error(f"当前版本不支持 Bowtie2 index 路径含空白字符: {prefix}")
                if not index_exists(prefix):
                    report.error(f"Bowtie2 index 不完整: {prefix}")
                else:
                    index_names = inspect_index_contigs(prefix, report)
                    if index_names is not None and index_names != fasta_names:
                        missing = sorted(fasta_names - index_names)
                        extra = sorted(index_names - fasta_names)
                        report.error(
                            "Bowtie2 index 与 FASTA contig 不一致; "
                            f"index 缺少 {missing[:5]}, index 额外 {extra[:5]}"
                        )
            else:
                report.note("未提供预建 Bowtie2 index，workflow 将根据本次 FASTA 构建")

    normalized = {
        "config": str(config_path.resolve()),
        "project_root": str(root),
        "output_dir": str(output_dir) if output_dir else None,
        "samples": str(samples_path) if samples_path else None,
        "comparisons": str(comparisons_path) if comparisons_path else None,
        "sample_count": len(samples),
        "comparison_count": len(comparisons),
        "groups": {key: len(value) for key, value in sorted(groups.items())},
        "errors": report.errors,
        "warnings": report.warnings,
    }
    return report, normalized


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--skip-file-checks", action="store_true")
    parser.add_argument("--json-report", type=Path)
    args = parser.parse_args()
    config_path = args.config.expanduser().resolve()
    if not config_path.is_file():
        print(f"[ERROR] config 文件不存在: {config_path}", file=sys.stderr)
        return 2
    report, normalized = validate(config_path, args.skip_file_checks)
    report.emit()
    if args.json_report:
        args.json_report.parent.mkdir(parents=True, exist_ok=True)
        args.json_report.write_text(json.dumps(normalized, indent=2, ensure_ascii=False) + "\n")
    return 1 if report.errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
