rule macs3_peak:
    input:
        unpack(macs3_inputs),
    output:
        peaks=os.path.join(DIRS["peaks"], "{comp}", "{comp}_peaks.narrowPeak"),
        summits=os.path.join(DIRS["peaks"], "{comp}", "{comp}_summits.bed"),
    params:
        outdir=os.path.join(DIRS["peaks"], "{comp}"),
        name="{comp}",
        call_summits=CALL_SUMMITS_FLAG,
        trackline=TRACKLINE_FLAG,
    conda:
        CORE_ENV
    threads: 1
    log:
        os.path.join(DIRS["logs"], "macs3", "{comp}.log"),
    benchmark:
        os.path.join(DIRS["benchmarks"], "macs3.{comp}.txt"),
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log}) {DIRS[benchmarks]}
        macs3 callpeak \
          -t {input.treatment} -c {input.control} \
          -f {MACS3[format]} \
          -g {REFERENCE[effective_genome_size]} \
          -n {params.name} \
          --outdir {params.outdir} \
          -q {MACS3[qvalue]} \
          --keep-dup {MACS3[keep_dup]} \
          {params.call_summits} {params.trackline} \
          > {log} 2>&1
        """


rule macs3_idr_replicate:
    input:
        treatment=os.path.join(DIRS["mapped"], "{sample}", "{sample}.filtered.bam"),
        control=idr_control_bam,
    output:
        peaks=os.path.join(DIRS["peaks"], "{comp}", "idr", "{sample}_peaks.narrowPeak"),
    params:
        outdir=os.path.join(DIRS["peaks"], "{comp}", "idr"),
        call_summits=CALL_SUMMITS_FLAG,
        trackline=TRACKLINE_FLAG,
    conda:
        CORE_ENV
    threads: 1
    log:
        os.path.join(DIRS["logs"], "idr_macs3", "{comp}.{sample}.log"),
    benchmark:
        os.path.join(DIRS["benchmarks"], "idr_macs3.{comp}.{sample}.txt"),
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log}) {DIRS[benchmarks]}
        macs3 callpeak \
          -t {input.treatment} -c {input.control} \
          -f {MACS3[format]} \
          -g {REFERENCE[effective_genome_size]} \
          -n {wildcards.sample} \
          --outdir {params.outdir} \
          -q {MACS3[qvalue]} \
          --keep-dup {MACS3[keep_dup]} \
          {params.call_summits} {params.trackline} \
          > {log} 2>&1
        """


rule idr_analysis:
    input:
        peaks=idr_replicate_peaks,
    output:
        table=os.path.join(DIRS["peaks"], "{comp}", "idr", "idr.tsv"),
        plot=os.path.join(DIRS["peaks"], "{comp}", "idr", "idr.tsv.png"),
    params:
        outdir=os.path.join(DIRS["peaks"], "{comp}", "idr"),
        sort_column=IDR_SORT_COLUMN,
    conda:
        IDR_ENV
    threads: 1
    log:
        os.path.join(DIRS["logs"], "idr", "{comp}.log"),
    benchmark:
        os.path.join(DIRS["benchmarks"], "idr.{comp}.txt"),
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log}) {DIRS[benchmarks]}
        test $(echo {input.peaks} | wc -w) -eq 2
        grep -v '^track' {input.peaks[0]} \
          | sort -k{params.sort_column},{params.sort_column}nr -k1,1V -k2,2n \
          > {params.outdir}/rep1.sorted.narrowPeak
        grep -v '^track' {input.peaks[1]} \
          | sort -k{params.sort_column},{params.sort_column}nr -k1,1V -k2,2n \
          > {params.outdir}/rep2.sorted.narrowPeak
        idr \
          --samples {params.outdir}/rep1.sorted.narrowPeak {params.outdir}/rep2.sorted.narrowPeak \
          --input-file-type narrowPeak \
          --rank {IDR[rank]} \
          --idr-threshold {IDR[threshold]} \
          --output-file {output.table} \
          --plot \
          > {log} 2>&1
        test -s {output.plot}
        """
