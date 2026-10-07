rule prepare_reference:
    input:
        fasta=str(REFERENCE_FASTA),
    output:
        fasta=REFERENCE_LINK,
        fai=REFERENCE_FAI,
        sizes=REFERENCE_SIZES,
    conda:
        CORE_ENV
    threads: 1
    log:
        os.path.join(DIRS["logs"], "prepare_reference.log"),
    shell:
        """
        mkdir -p $(dirname {output.fasta}) $(dirname {log})
        ln -sfn {input.fasta} {output.fasta}
        samtools faidx {output.fasta} 2> {log}
        cut -f1,2 {output.fai} > {output.sizes}
        """


rule bowtie2_index:
    input:
        fasta=REFERENCE_LINK,
        fai=REFERENCE_FAI,
    output:
        ready=BOWTIE2_INDEX_READY,
    params:
        prefix=BOWTIE2_INDEX_PREFIX,
        allow_build="true" if BOWTIE2_INDEX_ALLOW_BUILD else "false",
        checker=str(SCRIPT_DIR / "check_bowtie2_index.py"),
    conda:
        CORE_ENV
    threads: 16
    log:
        os.path.join(DIRS["logs"], "bowtie2_index.log"),
    benchmark:
        os.path.join(DIRS["benchmarks"], "bowtie2_index.txt"),
    shell:
        """
        mkdir -p $(dirname {output.ready}) $(dirname {params.prefix}) $(dirname {log}) {DIRS[benchmarks]}
        if ! python {params.checker} {params.prefix} >/dev/null; then
            if [[ {params.allow_build} != true ]]; then
                echo "Configured Bowtie2 index is missing or incomplete: {params.prefix}" > {log}
                exit 1
            fi
            bowtie2-build --threads {threads} {input.fasta} {params.prefix} > {log} 2>&1
        fi
        python {params.checker} {params.prefix} >> {log} 2>&1
        touch {output.ready}
        """


rule fastp:
    input:
        r1=lambda w: sample_fastq(w, "r1"),
        r2=lambda w: sample_fastq(w, "r2"),
    output:
        r1=(os.path.join(DIRS["trimmed"], "{sample}", "{sample}.R1.fastq.gz") if KEEP_TRIMMED else temp(os.path.join(DIRS["trimmed"], "{sample}", "{sample}.R1.fastq.gz"))),
        r2=(os.path.join(DIRS["trimmed"], "{sample}", "{sample}.R2.fastq.gz") if KEEP_TRIMMED else temp(os.path.join(DIRS["trimmed"], "{sample}", "{sample}.R2.fastq.gz"))),
        unpaired1=temp(os.path.join(DIRS["trimmed"], "{sample}", "{sample}.unpaired.R1.fastq.gz")),
        unpaired2=temp(os.path.join(DIRS["trimmed"], "{sample}", "{sample}.unpaired.R2.fastq.gz")),
        failed=temp(os.path.join(DIRS["trimmed"], "{sample}", "{sample}.failed.fastq.gz")),
        json=os.path.join(DIRS["trimmed"], "{sample}", "{sample}.fastp.json"),
        html=os.path.join(DIRS["trimmed"], "{sample}", "{sample}.fastp.html"),
    params:
        cut_front="--cut_front" if FASTP["cut_front"] else "",
        cut_tail="--cut_tail" if FASTP["cut_tail"] else "",
    conda:
        CORE_ENV
    threads: 8
    log:
        os.path.join(DIRS["logs"], "fastp", "{sample}.log"),
    benchmark:
        os.path.join(DIRS["benchmarks"], "fastp.{sample}.txt"),
    shell:
        """
        mkdir -p $(dirname {output.r1}) $(dirname {log}) {DIRS[benchmarks]}
        fastp \
          -i {input.r1} -I {input.r2} \
          -o {output.r1} -O {output.r2} \
          --unpaired1 {output.unpaired1} --unpaired2 {output.unpaired2} \
          --failed_out {output.failed} \
          --json {output.json} --html {output.html} \
          --thread {threads} \
          --qualified_quality_phred {FASTP[qualified_quality_phred]} \
          --length_required {FASTP[length_required]} \
          {params.cut_front} {params.cut_tail} \
          --cut_window_size {FASTP[cut_window_size]} \
          --cut_mean_quality {FASTP[cut_mean_quality]} \
          --average_qual {FASTP[average_qual]} \
          > {log} 2>&1
        """


rule bowtie2_align:
    input:
        index=BOWTIE2_INDEX_READY,
        r1=os.path.join(DIRS["trimmed"], "{sample}", "{sample}.R1.fastq.gz"),
        r2=os.path.join(DIRS["trimmed"], "{sample}", "{sample}.R2.fastq.gz"),
    output:
        bam=temp(os.path.join(DIRS["mapped"], "{sample}", "{sample}.unsorted.bam")),
    params:
        prefix=BOWTIE2_INDEX_PREFIX,
        flags=BOWTIE2_FLAGS,
        rg="@RG\\tID:{sample}\\tSM:{sample}\\tLB:CUTTag\\tPL:ILLUMINA",
        bowtie_log=os.path.join(DIRS["logs"], "bowtie2", "{sample}.summary.log"),
    conda:
        CORE_ENV
    threads: 16
    log:
        os.path.join(DIRS["logs"], "bowtie2", "{sample}.pipeline.log"),
    benchmark:
        os.path.join(DIRS["benchmarks"], "bowtie2.{sample}.txt"),
    shell:
        """
        mkdir -p $(dirname {output.bam}) $(dirname {log}) {DIRS[benchmarks]}
        bowtie2 {params.flags} --phred33 \
          -I {ALIGNMENT[insert_min]} -X {ALIGNMENT[insert_max]} \
          -p {threads} -x {params.prefix} \
          -1 {input.r1} -2 {input.r2} \
          2> {params.bowtie_log} \
        | samtools view -@ {threads} -b -F 0x04 - \
        | samtools addreplacerg -r '{params.rg}' -O bam -o {output.bam} - \
          > {log} 2>&1
        """


rule clean_bam:
    input:
        bam=os.path.join(DIRS["mapped"], "{sample}", "{sample}.unsorted.bam"),
    output:
        bam=os.path.join(DIRS["mapped"], "{sample}", "{sample}.filtered.bam"),
        bai=os.path.join(DIRS["mapped"], "{sample}", "{sample}.filtered.bam.bai"),
    params:
        chromosomes=CHROMOSOMES,
    conda:
        CORE_ENV
    threads: 8
    log:
        os.path.join(DIRS["logs"], "clean_bam", "{sample}.log"),
    benchmark:
        os.path.join(DIRS["benchmarks"], "clean_bam.{sample}.txt"),
    shell:
        """
        mkdir -p $(dirname {output.bam}) $(dirname {log}) {DIRS[benchmarks]}
        tmp_sorted={output.bam}.sorted.tmp.bam
        tmp_filtered={output.bam}.filtered.tmp.bam
        trap 'rm -f "$tmp_sorted" "$tmp_sorted.bai" "$tmp_filtered"' EXIT
        samtools sort -@ {threads} -o "$tmp_sorted" {input.bam} 2> {log}
        samtools index -@ {threads} "$tmp_sorted" >> {log} 2>&1
        samtools view -@ {threads} -b -q {ALIGNMENT[mapq_min]} -F 0x04 \
          -o "$tmp_filtered" "$tmp_sorted" {params.chromosomes} >> {log} 2>&1
        mv "$tmp_filtered" {output.bam}
        samtools index -@ {threads} {output.bam} >> {log} 2>&1
        rm -f "$tmp_sorted" "$tmp_sorted.bai"
        trap - EXIT
        """


rule samtools_merge:
    input:
        bams=group_bams,
    output:
        bam=os.path.join(DIRS["merged"], "{group}", "{group}.bam"),
        bai=os.path.join(DIRS["merged"], "{group}", "{group}.bam.bai"),
    conda:
        CORE_ENV
    threads: 8
    log:
        os.path.join(DIRS["logs"], "merge", "{group}.log"),
    benchmark:
        os.path.join(DIRS["benchmarks"], "merge.{group}.txt"),
    shell:
        """
        mkdir -p $(dirname {output.bam}) $(dirname {log}) {DIRS[benchmarks]}
        samtools merge -f -@ {threads} {output.bam} {input.bams} 2> {log}
        samtools index -@ {threads} {output.bam} >> {log} 2>&1
        """


rule samtools_stats:
    input:
        bam=os.path.join(DIRS["merged"], "{group}", "{group}.bam"),
        bai=os.path.join(DIRS["merged"], "{group}", "{group}.bam.bai"),
    output:
        flagstat=os.path.join(DIRS["merged"], "{group}", "{group}.flagstat.txt"),
        idxstats=os.path.join(DIRS["merged"], "{group}", "{group}.idxstats.txt"),
    conda:
        CORE_ENV
    threads: 1
    log:
        os.path.join(DIRS["logs"], "stats", "{group}.log"),
    shell:
        """
        mkdir -p $(dirname {output.flagstat}) $(dirname {log})
        samtools flagstat {input.bam} > {output.flagstat} 2> {log}
        samtools idxstats {input.bam} > {output.idxstats} 2>> {log}
        """
