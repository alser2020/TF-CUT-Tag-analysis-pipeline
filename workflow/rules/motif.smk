rule homer_prepare_genome:
    input:
        fasta=REFERENCE_LINK,
    output:
        marker=HOMER_PREPARSE_MARKER,
    params:
        size=MOTIF["homer"]["size"],
        preparsed_dir=os.path.join(DIRS["reference"], "preparsed"),
        prefix=lambda w: os.path.join(
            DIRS["reference"], "preparsed", f"genome.fa.{MOTIF['homer']['size']}"
        ),
        benchmark_dir=DIRS["benchmarks"],
    conda:
        HOMER_ENV
    threads: 1
    log:
        os.path.join(DIRS["logs"], "homer_preparse.log"),
    benchmark:
        os.path.join(DIRS["benchmarks"], "homer_preparse.txt"),
    shell:
        """
        mkdir -p {params.preparsed_dir} $(dirname {log}) {params.benchmark_dir}
        HOMER_PREFIX=$(dirname "$(dirname "$(command -v homerTools)")")
        export LD_LIBRARY_PATH="$HOMER_PREFIX/lib:${{LD_LIBRARY_PATH:-}}"
        rm -f {params.prefix}.*
        preparseGenome.pl {input.fasta} -size {params.size} > {log} 2>&1
        if grep -Eq 'CXXABI_|Illegal division by zero|Might have something wrong' {log}; then
            cat {log} >&2
            exit 1
        fi
        test -s {params.prefix}.gcbins
        touch {output.marker}
        """


rule homer_motif:
    input:
        summits=motif_input_bed,
        fasta=REFERENCE_LINK,
        preparsed=HOMER_PREPARSE_MARKER,
    output:
        directory(os.path.join(DIRS["motif"], "{comp}", "HOMER")),
    conda:
        HOMER_ENV
    threads: 16
    log:
        os.path.join(DIRS["logs"], "motif_homer", "{comp}.log"),
    benchmark:
        os.path.join(DIRS["benchmarks"], "motif_homer.{comp}.txt"),
    shell:
        """
        rm -rf {output}
        mkdir -p {output} $(dirname {log}) {DIRS[benchmarks]}
        grep -v '^track' {input.summits} > {output}/summits.clean.bed
        HOMER_PREFIX=$(dirname "$(dirname "$(command -v homerTools)")")
        export LD_LIBRARY_PATH="$HOMER_PREFIX/lib:${{LD_LIBRARY_PATH:-}}"
        findMotifsGenome.pl \
          {output}/summits.clean.bed {input.fasta} {output} \
          -size {MOTIF[homer][size]} \
          -len {HOMER_LENGTHS} \
          -p {threads} \
          > {log} 2>&1
        if grep -Eq 'CXXABI_|Illegal division by zero|Might have something wrong' {log}; then
            cat {log} >&2
            exit 1
        fi
        test -s {output}/knownResults.txt
        """


rule motif_fasta:
    input:
        summits=motif_input_bed,
        fasta=REFERENCE_LINK,
        fai=REFERENCE_FAI,
        sizes=REFERENCE_SIZES,
    output:
        fasta=os.path.join(DIRS["motif"], "{comp}", "summits.fasta"),
        bed=os.path.join(DIRS["motif"], "{comp}", "summits.slop.bed"),
    params:
        clean=os.path.join(DIRS["motif"], "{comp}", "summits.clean.bed"),
    conda:
        CORE_ENV
    threads: 1
    log:
        os.path.join(DIRS["logs"], "motif_fasta", "{comp}.log"),
    shell:
        """
        mkdir -p $(dirname {output.fasta}) $(dirname {log})
        grep -v '^track' {input.summits} > {params.clean}
        bedtools slop \
          -i {params.clean} -g {input.sizes} \
          -b {MOTIF[memechip][flank]} \
          > {output.bed} 2> {log}
        bedtools getfasta \
          -fi {input.fasta} -bed {output.bed} -fo {output.fasta} \
          >> {log} 2>&1
        rm -f {params.clean}
        """




rule memechip:
    input:
        fasta=os.path.join(DIRS["motif"], "{comp}", "summits.fasta"),
    output:
        complete=os.path.join(DIRS["motif"], "{comp}", "MEME-ChIP", ".complete"),
    params:
        output_dir=lambda wc: os.path.join(DIRS["motif"], wc.comp, "MEME-ChIP"),
        center_cut=lambda wc: MOTIF["memechip"]["center_cut"],
        min_width=lambda wc: MOTIF["memechip"]["min_width"],
        max_width=lambda wc: MOTIF["memechip"]["max_width"],
        meme_p=lambda wc: MOTIF["memechip"]["meme_p"],
        meme_nmotifs=lambda wc: MOTIF["memechip"]["meme_nmotifs"],
        streme_nmotifs=lambda wc: MOTIF["memechip"]["streme_nmotifs"],
        db_flag=lambda wc: (
            f'-db {shlex.quote(str(project_path(MOTIF["memechip"]["known_motif_db"]))) }'
            if MOTIF["memechip"].get("known_motif_db")
            else ""
        ),
    conda:
        MEME_ENV
    threads: 16
    log:
        os.path.join(DIRS["logs"], "motif_memechip", "{comp}.log"),
    benchmark:
        os.path.join(DIRS["benchmarks"], "motif_memechip.{comp}.txt"),
    shell:
        """
        OUT={params.output_dir}
        rm -rf "$OUT"
        mkdir -p "$OUT" $(dirname {log}) {DIRS[benchmarks]}
        meme-chip {input.fasta} \\
          -ccut {params.center_cut} \\
          -meme-p {params.meme_p} \\
          -meme-nmotifs {params.meme_nmotifs} \\
          -streme-nmotifs {params.streme_nmotifs} \\
          -minw {params.min_width} \\
          -maxw {params.max_width} \\
          {params.db_flag} \\
          -oc "$OUT" \\
          > {log} 2>&1
        test -s "$OUT/meme-chip.html" || find "$OUT" -maxdepth 2 -type f -name '*.html' -size +0c -print -quit | grep -q .
        test -s "$OUT/combined.meme" || test -s "$OUT/summary.tsv" || find "$OUT" -maxdepth 2 -type f \\
          \( -name 'meme.txt' -o -name 'centrimo.txt' -o -name 'tomtom.txt' \) -size +0c -print -quit | grep -q .
        touch {output.complete}
        """

