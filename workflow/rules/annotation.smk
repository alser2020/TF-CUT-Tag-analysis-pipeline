rule annotate_peaks:
    input:
        summits=os.path.join(DIRS["peaks"], "{comp}", "{comp}_summits.bed"),
        fasta=REFERENCE_LINK,
        gtf=str(REFERENCE_GTF),
    output:
        table=os.path.join(DIRS["peaks"], "{comp}", "{comp}.annotation.tsv"),
    params:
        clean=os.path.join(DIRS["peaks"], "{comp}", "{comp}.summits.clean.bed"),
    conda:
        HOMER_ENV
    threads: 1
    log:
        os.path.join(DIRS["logs"], "annotation", "{comp}.log"),
    benchmark:
        os.path.join(DIRS["benchmarks"], "annotation.{comp}.txt"),
    shell:
        """
        mkdir -p $(dirname {output.table}) $(dirname {log}) {DIRS[benchmarks]}
        grep -v '^track' {input.summits} > {params.clean}
        HOMER_PREFIX=$(dirname "$(dirname "$(command -v homerTools)")")
        export LD_LIBRARY_PATH="$HOMER_PREFIX/lib:${{LD_LIBRARY_PATH:-}}"
        annotatePeaks.pl {params.clean} {input.fasta} -gtf {input.gtf} \
          > {output.table} 2> {log}
        if grep -Eq 'CXXABI_|Illegal division by zero' {log}; then
            cat {log} >&2
            exit 1
        fi
        rm -f {params.clean}
        """


rule filter_promoter_peaks:
    input:
        annotation=os.path.join(DIRS["peaks"], "{comp}", "{comp}.annotation.tsv"),
        summits=os.path.join(DIRS["peaks"], "{comp}", "{comp}_summits.bed"),
    output:
        bed=os.path.join(DIRS["peaks"], "{comp}", "{comp}.promoter.bed"),
        summits=os.path.join(DIRS["peaks"], "{comp}", "{comp}.promoter_summits.bed"),
    params:
        script=str(SCRIPT_DIR / "filter_homer_annotation.py"),
    conda:
        CORE_ENV
    threads: 1
    log:
        os.path.join(DIRS["logs"], "promoter", "{comp}.log"),
    shell:
        """
        mkdir -p $(dirname {output.bed}) $(dirname {log})
        python {params.script} \
          --annotation {input.annotation} \
          --summits {input.summits} \
          --max-distance {ANNOTATION[promoter_max_abs_tss_distance]} \
          --bed {output.bed} \
          --summit-bed {output.summits} \
          > {log} 2>&1
        """
