rule bam_to_bigwig_group:
    input:
        bam=os.path.join(DIRS["merged"], "{group}", "{group}.bam"),
        bai=os.path.join(DIRS["merged"], "{group}", "{group}.bam.bai"),
    output:
        bw=os.path.join(DIRS["deeptools"], "group_bigwig", "{group}.bw"),
    params:
        ignore_duplicates=IGNORE_DUPLICATES_FLAG,
    conda:
        CORE_ENV
    threads: 16
    log:
        os.path.join(DIRS["logs"], "bigwig_group", "{group}.log"),
    benchmark:
        os.path.join(DIRS["benchmarks"], "bigwig_group.{group}.txt"),
    shell:
        """
        mkdir -p $(dirname {output.bw}) $(dirname {log}) {DIRS[benchmarks]}
        bamCoverage -b {input.bam} -o {output.bw} \
          --binSize {COVERAGE[bin_size]} \
          --normalizeUsing {COVERAGE[normalize_using]} \
          --effectiveGenomeSize {REFERENCE[effective_genome_size]} \
          --extendReads {COVERAGE[extend_reads]} \
          --smoothLength {COVERAGE[smooth_length]} \
          {params.ignore_duplicates} \
          --numberOfProcessors {threads} \
          > {log} 2>&1
        """


rule bam_to_bigwig_sample:
    input:
        bam=os.path.join(DIRS["mapped"], "{sample}", "{sample}.filtered.bam"),
        bai=os.path.join(DIRS["mapped"], "{sample}", "{sample}.filtered.bam.bai"),
    output:
        bw=os.path.join(DIRS["deeptools"], "sample_bigwig", "{sample}.bw"),
    params:
        ignore_duplicates=IGNORE_DUPLICATES_FLAG,
    conda:
        CORE_ENV
    threads: 16
    log:
        os.path.join(DIRS["logs"], "bigwig_sample", "{sample}.log"),
    benchmark:
        os.path.join(DIRS["benchmarks"], "bigwig_sample.{sample}.txt"),
    shell:
        """
        mkdir -p $(dirname {output.bw}) $(dirname {log}) {DIRS[benchmarks]}
        bamCoverage -b {input.bam} -o {output.bw} \
          --binSize {COVERAGE[bin_size]} \
          --normalizeUsing {COVERAGE[normalize_using]} \
          --effectiveGenomeSize {REFERENCE[effective_genome_size]} \
          --extendReads {COVERAGE[extend_reads]} \
          --smoothLength {COVERAGE[smooth_length]} \
          {params.ignore_duplicates} \
          --numberOfProcessors {threads} \
          > {log} 2>&1
        """


rule sample_correlation:
    input:
        bigwigs=comparison_sample_bigwigs,
    output:
        matrix=os.path.join(DIRS["deeptools"], "{comp}", "correlation.npz"),
        heatmap=os.path.join(DIRS["deeptools"], "{comp}", "correlation_heatmap.png"),
        table=os.path.join(DIRS["deeptools"], "{comp}", "correlation_matrix.tsv"),
        pca=os.path.join(DIRS["deeptools"], "{comp}", "PCA.png"),
    params:
        labels=comparison_labels,
        outdir=os.path.join(DIRS["deeptools"], "{comp}"),
    conda:
        CORE_ENV
    threads: 16
    log:
        os.path.join(DIRS["logs"], "correlation", "{comp}.log"),
    benchmark:
        os.path.join(DIRS["benchmarks"], "correlation.{comp}.txt"),
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log}) {DIRS[benchmarks]}
        multiBigwigSummary bins \
          -b {input.bigwigs} \
          --binSize {QC[summary_bin_size]} \
          --numberOfProcessors {threads} \
          -o {output.matrix} \
          > {log} 2>&1
        plotCorrelation \
          -in {output.matrix} \
          --corMethod {QC[correlation_method]} \
          --whatToPlot heatmap \
          --removeOutliers \
          --labels {params.labels} \
          --outFileCorMatrix {output.table} \
          -o {output.heatmap} \
          >> {log} 2>&1
        plotPCA \
          -in {output.matrix} \
          --labels {params.labels} \
          -o {output.pca} \
          >> {log} 2>&1
        """


rule plot_fingerprint:
    input:
        bams=comparison_sample_bams,
        bais=comparison_sample_bais,
    output:
        png=os.path.join(DIRS["deeptools"], "{comp}", "fingerprint.png"),
        counts=os.path.join(DIRS["deeptools"], "{comp}", "fingerprint_counts.tsv"),
    conda:
        CORE_ENV
    threads: 16
    log:
        os.path.join(DIRS["logs"], "fingerprint", "{comp}.log"),
    benchmark:
        os.path.join(DIRS["benchmarks"], "fingerprint.{comp}.txt"),
    shell:
        """
        mkdir -p $(dirname {output.png}) $(dirname {log}) {DIRS[benchmarks]}
        plotFingerprint \
          -b {input.bams} \
          --smartLabels \
          --outRawCounts {output.counts} \
          --numberOfSamples {QC[fingerprint_samples]} \
          --skipZeros \
          --binSize {QC[fingerprint_bin_size]} \
          --numberOfProcessors {threads} \
          -o {output.png} \
          > {log} 2>&1
        """


rule tss_enrichment:
    input:
        bigwigs=comparison_sample_bigwigs,
        gtf=str(REFERENCE_GTF),
    output:
        matrix=os.path.join(DIRS["deeptools"], "{comp}", "TSS_matrix.gz"),
        heatmap=os.path.join(DIRS["deeptools"], "{comp}", "TSS_heatmap.png"),
        profile=os.path.join(DIRS["deeptools"], "{comp}", "TSS_profile.png"),
    params:
        labels=comparison_labels,
    conda:
        CORE_ENV
    threads: 16
    log:
        os.path.join(DIRS["logs"], "tss", "{comp}.log"),
    benchmark:
        os.path.join(DIRS["benchmarks"], "tss.{comp}.txt"),
    shell:
        """
        mkdir -p $(dirname {output.matrix}) $(dirname {log}) {DIRS[benchmarks]}
        computeMatrix reference-point \
          --referencePoint TSS \
          -S {input.bigwigs} \
          -R {input.gtf} \
          --beforeRegionStartLength {QC[tss_before]} \
          --afterRegionStartLength {QC[tss_after]} \
          --binSize {QC[tss_bin_size]} \
          --skipZeros \
          --samplesLabel {params.labels} \
          --numberOfProcessors {threads} \
          -o {output.matrix} \
          > {log} 2>&1
        plotHeatmap \
          -m {output.matrix} \
          --colorMap viridis \
          --zMin 0 \
          -o {output.heatmap} \
          >> {log} 2>&1
        plotProfile \
          -m {output.matrix} \
          --perGroup \
          -o {output.profile} \
          >> {log} 2>&1
        """
