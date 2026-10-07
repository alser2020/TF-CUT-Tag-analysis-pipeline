rule multiqc:
    input:
        fastp=expand(os.path.join(DIRS["trimmed"], "{sample}", "{sample}.fastp.json"), sample=samples),
        flagstat=expand(os.path.join(DIRS["merged"], "{group}", "{group}.flagstat.txt"), group=groups),
        idxstats=expand(os.path.join(DIRS["merged"], "{group}", "{group}.idxstats.txt"), group=groups),
    output:
        html=os.path.join(DIRS["qc"], "multiqc_report.html"),
        zip=os.path.join(DIRS["qc"], "multiqc_report_data.zip"),
    conda:
        CORE_ENV
    threads: 1
    log:
        os.path.join(DIRS["logs"], "multiqc.log"),
    shell:
        """
        mkdir -p {DIRS[qc]} $(dirname {log})
        multiqc \
          {input.fastp} {input.flagstat} {input.idxstats} \
          -o {DIRS[qc]} \
          --filename multiqc_report.html \
          --zip-data-dir \
          > {log} 2>&1
        """

rule package_results:
    input:
        upstream=analysis_targets,
        config=str(PROJECT_DIR / "config.yaml"),
        samples=SAMPLE_SHEET,
        comparisons=COMPARISON_SHEET,
        reference_sizes=REFERENCE_SIZES,
        reference_fai=REFERENCE_FAI,
    output:
        archive=PACKAGE_ARCHIVE,
        manifest=PACKAGE_MANIFEST,
        checksum=PACKAGE_CHECKSUM,
    params:
        project_dir=str(PROJECT_DIR),
    threads: 1
    log:
        os.path.join(DIRS["logs"], "package_results.log"),
    benchmark:
        os.path.join(DIRS["benchmarks"], "package_results.txt"),
    run:
        if not PACKAGE_ENABLED:
            raise ValueError("package_results requested but outputs.package.enabled is false")
        project = Path(params.project_dir)
        archive = Path(output.archive)
        manifest = Path(output.manifest)
        checksum = Path(output.checksum)
        archive.parent.mkdir(parents=True, exist_ok=True)
        manifest.parent.mkdir(parents=True, exist_ok=True)
        log_path = Path(log[0])
        log_path.parent.mkdir(parents=True, exist_ok=True)
        package_manifest = project / ".cuttag_package_manifest.txt"
        included = [
            Path("config.yaml"),
            Path(input.samples).relative_to(project),
            Path(input.comparisons).relative_to(project),
            Path(input.reference_sizes).relative_to(project),
            Path(input.reference_fai).relative_to(project),
        ]
        for root in ("results/04.deeptools", "results/05.peaks", "results/06.motif", "results/07.qc", "results/logs", "results/benchmarks"):
            root_path = project / root
            if root_path.exists():
                included.extend(path.relative_to(project) for path in root_path.rglob("*") if path.is_file())
        included = sorted(set(included), key=str)
        package_manifest.write_text(
            "# CUT&Tag result package\n"
            f"# project: {project.name}\n"
            f"# archive: {archive.name}\n"
            + "\n".join(str(path) for path in included)
            + "\n",
            encoding="utf-8",
        )
        try:
            import subprocess
            tar_files = [str(path) for path in included] + [package_manifest.relative_to(project).as_posix()]
            subprocess.run(
                ["tar", "-czf", str(archive), "-C", str(project), *tar_files],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=log_path.open("w", encoding="utf-8"),
            )
            manifest.write_text(package_manifest.read_text(encoding="utf-8"), encoding="utf-8")
            with checksum.open("w", encoding="utf-8") as handle:
                subprocess.run(["sha256sum", str(archive)], check=True, stdout=handle)
            log_path.write_text(f"Created {archive}\nIncluded files: {len(included)}\n", encoding="utf-8")
        finally:
            package_manifest.unlink(missing_ok=True)
