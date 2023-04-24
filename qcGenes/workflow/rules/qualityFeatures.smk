# =============================================================================
# Snakemake workflow to compute quality features
# =============================================================================

# ------------------------------------------------------------------------------
rule qualityFeatures:
    input:
        expand("output/main/scores/{dataid}.scores.txt", dataid=DATAIDS_MAIN),
        expand("output/batched/scores/{dataid}.scores.txt", dataid=DATAIDS_BATCH)
#        get_all_aggregate_scores_files,
#	get_all_batched_aggregate_scores_files
#	,        get_multiqc_files,
#	get_batched_multiqc_files

rule allFastqc:
    input:
        expand("output/qc/{dataid}/fastqc/{sample}.fastqc.txt", zip, dataid=SAMPLES_ALL_GEO_COL, sample=SAMPLES_ALL_RUN_COL)

#rule allMultiqcsalmon:
#    input:
#        get_multiqc_salmon_files

rule allBowtie:
    input:
        expand("output/qc/{dataid}/bowtie2/{sample}.txt", zip, dataid=SAMPLES_ALL_GEO_COL, sample=SAMPLES_ALL_RUN_COL)

rule allScores:
    input:
        expand("output/qc/{dataid}/scores/{sample}.score.txt", zip, dataid=SAMPLES_ALL_GEO_COL, sample=SAMPLES_ALL_RUN_COL)


# ------------------------------------------------------------------------------
rule fastqc:
    input:
#	path=rules.download_sample.output
        path="data/datasets/{dataid}/{sample}_1.fastq.gz"
    output:
        arch="output/qc/{dataid}/fastqc/{sample}_fastqc.zip",
        html="output/qc/{dataid}/fastqc/{sample}_fastqc.html",
        summ="output/qc/{dataid}/fastqc/{sample}.fastqc.txt"
    group:
        "qualityFeatures"
    log:
        "output/logs/qc/{dataid}.{sample}.fastqc.log"
    conda:
        "../envs/fastqc.yaml"
    params:
        in_arch_path="{sample}_fastqc/summary.txt",
        out_dir="output/qc/{dataid}/fastqc/",
        extracted_file_path="output/qc/{dataid}/fastqc/summary.txt",
        tmp_input="data/datasets/{dataid}/{sample}.fastq.gz"
#	,
#        tmp_arch="output/qc/{dataid}/fastqc/{sample}_fastqc.zip",
#        tmp_html="output/qc/{dataid}/fastqc/{sample}_fastqc.html"
    threads:
        1
    shell:
        "ln -frs {input} {params.tmp_input} > {log} 2>&1 && "
        "fastqc -q --noextract -t {threads} --outdir {params.out_dir} {params.tmp_input} >> {log} 2>&1 && "
        "rm -f {params.tmp_input} && "
        "unzip -j -q -o '{output.arch}' '{params.in_arch_path}' -d '{params.out_dir}' >> {log} 2>&1 && "
        "mv '{params.extracted_file_path}' '{output.summ}' >> {log} 2>&1 "

# ------------------------------------------------------------------------------
rule bowtie:
    input:
        path="data/datasets/{dataid}/{sample}_1.fastq.gz"
#	path=rules.download_sample.output
    output:
        bam="data/datasets/{dataid}/bowtie2/{sample}.bam",
        txt="output/qc/{dataid}/bowtie2/{sample}.txt"
    log:
        "output/logs/qc/{dataid}.{sample}.bowtie.log"
    benchmark:
        "output/logs/qc/{dataid}.{sample}.bowtie.bench"
    conda:
        "../envs/bowtie.yaml"
    params:
        idx=config["BOWTIE_IDX"],
        sam="data/datasets/{dataid}/bowtie2/{sample}.sam",
        unsorted_sam="data/datasets/{dataid}/bowtie2/{sample}.unsorted.sam",
        tmp_fq="output/qc/{dataid}/bowtie2/{sample}.tmp.fastq.gz",
        tmp_txt="output/qc/{dataid}/bowtie2/{sample}.txt.tmp",
        max_reads=config['BOWTIE2_MAX_READS']
    threads:
        config['BOWTIE2_THREADS']
    priority: 50
    shell:
        "if [ {params.max_reads} -gt 0 ] ; then seqtk sample -s 100 {input} {params.max_reads} | gzip -1 > {params.tmp_fq} ; fi && "
        "if [ {params.max_reads} -lt 0 ] ; then cp -f {input} {params.tmp_fq} ; fi && "
        "bowtie2 -x {params.idx} -q {params.tmp_fq} -S {params.sam} -p {threads} --mm 2> {params.tmp_txt} > {log} && "
        "grep -v Warning {params.tmp_txt} > {output.txt} 2>> {log} && "
        "samtools view -bS {params.sam} -o {params.unsorted_sam} >> {log} 2>&1 && "
        "samtools sort {params.unsorted_sam} -o {output.bam} >> {log} 2>&1 && "
        "rm -f {params.sam} {params.unsorted_sam} {params.tmp_txt} {params.tmp_fq}"

# ------------------------------------------------------------------------------
rule reads_all_anno:
    input:
        rules.bowtie.output.bam
#        "data/datasets/{dataid}/bowtie2/{sample}.bam"
    output:
        ann="output/qc/{dataid}/reads_anno/{sample}.readsannot.txt",
        tss="output/qc/{dataid}/tss_anno/{sample}.tss.txt"
    log:
        "output/logs/qc/{dataid}.{sample}.reads_all_anno.log"
    benchmark:
        "output/logs/qc/{dataid}.{sample}.reads_all_anno.bench"
    conda:
        "../envs/readsAnno.yaml"
    params:
        script=config["WRAPPER_SCRIPT_READS_ALL_ANNO_PY"],
        script_ann_r=config["WRAPPER_SCRIPT_READS_ALL_ANNO_ANN_R"],
        script_tss_r=config["WRAPPER_SCRIPT_READS_ALL_ANNO_TSS_R"],
        script_ann_r_name=config["WRAPPER_SCRIPT_READS_ALL_ANNO_ANN_R_NAME"],
        script_tss_r_name=config["WRAPPER_SCRIPT_READS_ALL_ANNO_TSS_R_NAME"],
        species=config["SPECIES_COMMON_NAME"]
    shell:
        "python {params.script} {input} {params.species} {output.ann} {output.tss} > {log} 2>&1 &&"
        "rm -f Rplots.pdf"

# ------------------------------------------------------------------------------
rule scorer:
    input:
        raw=rules.fastqc.output.summ,
        map=rules.bowtie.output.txt,
        loc=rules.reads_all_anno.output.ann,
        tss=rules.reads_all_anno.output.tss
#        raw="output/qc/{dataid}/fastqc/{sample}.fastqc.txt",
#        map="output/qc/{dataid}/bowtie2/{sample}.txt",
#        loc="output/qc/{dataid}/reads_anno/{sample}.readsannot.txt",
#        tss="output/qc/{dataid}/tss_anno/{sample}.tss.txt"
    output:
        "output/qc/{dataid}/scores/{sample}.score.txt"
    group:
        "qualityFeatures"
    log:
        "output/logs/qc/{dataid}.{sample}.scorer.log"
    conda:
        "../envs/seqQscorer.yaml"
    params:
#        raw="../../output/qc/{dataid}/fastqc/{sample}.fastqc.txt",
#        map="../../output/qc/{dataid}/bowtie2/{sample}.txt",
#        loc="../../output/qc/{dataid}/reads_anno/{sample}.readsannot.txt",
#        tss="../../output/qc/{dataid}/tss_anno/{sample}.tss.txt",
#        out="../../output/qc/{dataid}/scores/{sample}.score.txt",
#        out="../../output/qc/{dataid}/scores/{sample}.score.txt",
        species=config["SPECIES_COMMON_NAME"]
    shell:
        "cd lib/seqQscorer.git/ && "
#        "python seqQscorer.py --spec {params.species} --assay RNA-seq --rt single-ended --raw {params.raw} --map {params.map} --loc {params.loc} --tss {params.tss} --out {params.out} > ../../{log} 2>&1 && "
        "python seqQscorer.py --spec None --assay None --rt None --raw ../../{input.raw} --map ../../{input.map} --loc ../../{input.loc} --tss ../../{input.tss} --out ../../{output} > ../../{log} 2>&1 && "
        "cd ../.."

# ------------------------------------------------------------------------------


rule aggregate_scores:
    input:
#        expand("output/qc/{dataid}/scores/{sample}.score.txt", dataid="{dataid}", sample=get_samples_by_analysis_and_dataid_BIS("{{analysis}}", "{{dataid}}"))
#        expand("output/qc/{dataid}/scores/{sample}.score.txt", dataid="{dataid}", sample=get_samples_by_analysis_and_dataid_BIS("main", "GSE99816"))
#        expand("output/qc/{dataid}/scores/{sample}.score.txt", dataid="{dataid}", sample=get_samples_by_analysis_and_dataid_BIS("main", "GSE99816"))
#        expand("output/qc/{dataid}/scores/{sample}.score.txt", dataid="{dataid}", sample=get_samples_by_analysis_and_dataid)
        get_samples_scores_files_by_analysis_and_GEOSeries
    output:
        "output/{analysis}/scores/{dataid}.scores.txt"
    group:
        "qualityFeatures"
    log:
        "output/logs/{analysis}/qc/{dataid}.aggregate_scores.log"
    conda:
        "../envs/linux.yaml"
    shell:
        "cat {input} > {output} 2> {log}"


# ------------------------------------------------------------------------------
#rule aggregate_batched_scores:
#    input:
#        get_batched_samples_scores_files_by_GEOSeries
#    output:
#        "output/batched/scores/{dataid}.scores.txt"
#    group:
#        "qualityFeatures"
#    log:
#        "output/logs/batched/qualityFeatures/{dataid}.aggregate_scores.log"
#    conda:
#        "../envs/linux.yaml"
#    shell:
#        "cat {input} > {output} 2> {log}"



# ------------------------------------------------------------------------------
rule picard_rrna:
    input:
        gtf=config["ANNOT_PATH"],
        bam=rules.bowtie.output.bam
#        bam="data/datasets/{dataid}/bowtie2/{sample}.bam"
    output:
        temp("output/qc/{dataid}/picard/{dataid}.{sample}.rRNA.bed")
    log:
        "output/logs/qc/{dataid}.{sample}.picard_rrna.log"
    benchmark:
        "output/logs/qc/{dataid}.{sample}.picard_rrna.bench"
    conda:
        "../envs/samtools.yaml"
    shell:
        """samtools view -H {input.bam} > {output} && """
        """zcat {input.gtf} | grep rRNA | awk 'BEGIN{{OFS="\\t"}}{{print $1, $4-1, $5, $7, NR}}' >> {output} 2> {log}"""

rule picard:
    input:
        bam=rules.bowtie.output.bam,
        bed=rules.picard_rrna.output,
        ref=config["PICARD_REFFLAT_FILE"],
#        bam="data/datasets/{dataid}/bowtie2/{sample}.bam",
#        bed="output/qc/{dataid}/picard/{dataid}.{sample}.rRNA.bed"
    output:
        RNA_Metrics="output/qc/{dataid}/picard/{sample}.RNA_Metrics",
        ALN_Metrics="output/qc/{dataid}/picard/{sample}.ALN_Metrics"
    log:
        "output/logs/qc/{dataid}.{sample}.picard.log"
    benchmark:
        "output/logs/qc/{dataid}.{sample}.picard.bench"
    conda:
        "../envs/picard.yaml"
    shell:
        "picard CollectRnaSeqMetrics I={input.bam} O={output.RNA_Metrics} REF_FLAT={input.ref} STRAND=NONE RIBOSOMAL_INTERVALS={input.bed} > {log} 2>&1 && "
        "picard CollectAlignmentSummaryMetrics INPUT={input.bam} OUTPUT={output.ALN_Metrics} >> {log} 2>&1"

# ------------------------------------------------------------------------------
rule multiqc:
    input:
        expand("output/qc/{{dataid}}/picard/{sample}.RNA_Metrics", sample=get_all_samples_by_dataid),
        expand("data/datasets/{{dataid}}/bowtie2/{sample}.bam", sample=get_all_samples_by_dataid),
        expand("output/qc/{{dataid}}/fastqc/{sample}.fastqc.txt", sample=get_all_samples_by_dataid)
    output:
        "output/qc/{dataid}/multiqc/{dataid}_multiqc_report.html"
    log:
        "output/logs/qc/{dataid}.multiqc.log"
    benchmark:
        "output/logs/qc/{dataid}.multiqc.bench"
    conda:
        "../envs/multiqc.yaml"
    params:
        outdir="output/qc/{dataid}/multiqc",
        picarddir="output/qc/{dataid}/picard",
        bowtiedir="output/qc/{dataid}/bowtie2",
        fastqcdir="output/qc/{dataid}/fastqc"
    shell:
        "multiqc -f -i {wildcards.dataid} -o {params.outdir} {params.picarddir} {params.bowtiedir} {params.fastqcdir} > {log} 2>&1"

# ------------------------------------------------------------------------------
#rule multiqc_salmon:
#    input:
#        "output/main/scores/{dataid}.scores.txt",
#        "output/main/RASflowResults/{dataid}/trans/tpmFile/all.samples.tpm.tsv",
#        get_picard_files_by_GEOSeries,
#	get_batched_picard_files_by_GEOSeries
#    output:
#        "output/main/qc/{dataid}/multiqc_salmon/{dataid}_multiqc_report.html"
#    log:
#        "output/logs/main/qc/{dataid}.multiqc_salmon.log"
#    benchmark:
#        "output/logs/main/qc/{dataid}.multiqc_salmon.bench"
#    conda:
#        "../envs/multiqc.yaml"
#    params:
#        outdir="output/main/qc/{dataid}/multiqc_salmon",
#        picarddir="output/qc/{dataid}/picard",
#        bowtiedir="output/qc/{dataid}/bowtie2",
#        fastqcdir="output/qc/{dataid}/fastqc",
#        salmondir="output/main/RASflowResults/{dataid}/trans/quant"
#    shell:
#        "multiqc -f -i {wildcards.dataid} -o {params.outdir} {params.picarddir} {params.bowtiedir} {params.fastqcdir} {params.salmondir} > {log} 2>&1"

# ------------------------------------------------------------------------------
#rule batched_multiqc_salmon:
#    input:
#        "output/batched/scores/{dataid}.scores.txt",
#        "output/batched/RASflowResults/{dataid}/trans/tpmFile/all.samples.tpm.tsv",
#        get_picard_files_by_GEOSeries,
#	get_batched_picard_files_by_GEOSeries
#    output:
#        "output/batched/qc/{dataid}/multiqc_salmon/{dataid}_multiqc_report.html"
#    log:
#        "output/logs/batched/qc/{dataid}.multiqc_salmon.log"
#    benchmark:
#        "output/logs/batched/qc/{dataid}.multiqc_salmon.bench"
#    conda:
#        "../envs/multiqc.yaml"
#    params:
#        outdir="output/batched/qc/{dataid}/multiqc_salmon",
#        picarddir="output/qc/{dataid}/picard",
#        bowtiedir="output/qc/{dataid}/bowtie2",
#        fastqcdir="output/qc/{dataid}/fastqc",
#        salmondir="output/batched/RASflowResults/{dataid}/trans/quant"
#    shell:
#        "multiqc -f -i {wildcards.dataid} -o {params.outdir} {params.picarddir} {params.bowtiedir} {params.fastqcdir} {params.salmondir} > {log} 2>&1"

# ==============================================================================
# BATCH ANALYSIS
# ==============================================================================

# ------------------------------------------------------------------------------
#rule batchedQualityFeatures:
#    input:
#        get_all_batched_aggregate_scores_files
#        expand("output/batched/scores/{dataid}.scores.txt", dataid=DATAIDS_BATCH)

