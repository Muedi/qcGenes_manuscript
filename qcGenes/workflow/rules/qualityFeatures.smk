# =============================================================================
# Snakemake workflow to compute quality features
# =============================================================================

# ------------------------------------------------------------------------------
rule qualityFeatures:
    input:
        expand("output/main/scores/{dataid}.scores.txt", dataid=DATAIDS_MAIN)

rule allFastqc:
    input:
        expand("output/qc/{dataid}/fastqc/{sample}.fastqc.txt", zip, dataid=SAMPLES_MAIN_GEO_COL, sample=SAMPLES_MAIN_RUN_COL)

rule allBowtie:
    input:
        expand("output/qc/{dataid}/bowtie2/{sample}.txt", zip, dataid=SAMPLES_MAIN_GEO_COL, sample=SAMPLES_MAIN_RUN_COL)

rule allScores:
    input:
        expand("output/qc/{dataid}/scores/{sample}.score.txt", zip, dataid=SAMPLES_MAIN_GEO_COL, sample=SAMPLES_MAIN_RUN_COL)

rule allMultiQC:
    input:
        expand("output/qc/{dataid}/multiqc/{dataid}_multiqc_report.html", zip, dataid=DATAIDS_MAIN)


# ------------------------------------------------------------------------------
rule fastqc:
    input:
        path="data/datasets/{dataid}/{sample}_1.fastq.gz"
    output:
        arch="output/qc/{dataid}/fastqc/{sample}_fastqc.zip",
        html="output/qc/{dataid}/fastqc/{sample}_fastqc.html",
        summ="output/qc/{dataid}/fastqc/{sample}.fastqc.txt"
    group:
        "qualityFeatures"
    log:
        "output/logs/qc/{dataid}.{sample}.fastqc.log"
    params:
        in_arch_path="{sample}_fastqc/summary.txt",
        out_dir="output/qc/{dataid}/fastqc/",
        extracted_file_path="output/qc/{dataid}/fastqc/summary.txt",
        tmp_input="data/datasets/{dataid}/{sample}.fastq.gz"
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
        fastq_path="data/datasets/{dataid}/{sample}_1.fastq.gz",
        bowtie_index=rules.extract_bowtie_idx.output
    output:
        bam="data/datasets/{dataid}/bowtie2/{sample}.bam",
        txt="output/qc/{dataid}/bowtie2/{sample}.txt"
    log:
        "output/logs/qc/{dataid}.{sample}.bowtie.log"
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
        "if [ {params.max_reads} -gt 0 ] ; then seqtk sample -s 100 {input.fastq_path} {params.max_reads} | gzip -1 > {params.tmp_fq} ; fi && "
        "if [ {params.max_reads} -lt 0 ] ; then cp -f {input.fastq_path} {params.tmp_fq} ; fi && "
        "bowtie2 -x {params.idx} -q {params.tmp_fq} -S {params.sam} -p {threads} --mm 2> {params.tmp_txt} > {log} && "
        "grep -v Warning {params.tmp_txt} > {output.txt} 2>> {log} && "
        "samtools view -bS {params.sam} -o {params.unsorted_sam} >> {log} 2>&1 && "
        "samtools sort {params.unsorted_sam} -o {output.bam} >> {log} 2>&1 && "
        "rm -f {params.sam} {params.unsorted_sam} {params.tmp_txt} {params.tmp_fq}"

# ------------------------------------------------------------------------------
rule reads_all_anno:
    input:
        rules.bowtie.output.bam
    output:
        ann="output/qc/{dataid}/reads_anno/{sample}.readsannot.txt",
        tss="output/qc/{dataid}/tss_anno/{sample}.tss.txt"
    log:
        "output/logs/qc/{dataid}.{sample}.reads_all_anno.log"
    params:
        script=config["WRAPPER_SCRIPT_READS_ALL_ANNO_PY"],
        script_ann_r=config["WRAPPER_SCRIPT_READS_ALL_ANNO_ANN_R"],
        script_tss_r=config["WRAPPER_SCRIPT_READS_ALL_ANNO_TSS_R"],
        species=config["SPECIES_COMMON_NAME"],
        max_reads=config["READS_ANNOT_MAX_READS"]
    threads:
        config['BOWTIE2_THREADS']
    shell:
        "python3 {params.script} {input} {params.species} {output.ann} {output.tss} {params.script_ann_r} {params.script_tss_r} {params.max_reads} > {log} 2>&1 &&"
        "rm -f Rplots.pdf"

# ------------------------------------------------------------------------------
rule scorer:
    input:
        raw=rules.fastqc.output.summ,
        map=rules.bowtie.output.txt,
        loc=rules.reads_all_anno.output.ann,
        tss=rules.reads_all_anno.output.tss
    output:
        "output/qc/{dataid}/scores/{sample}.score.txt"
    group:
        "qualityFeatures"
    log:
        "output/logs/qc/{dataid}.{sample}.scorer.log"
    params:
        species=config["SPECIES_COMMON_NAME"]
    shell:
        "cd lib/seqQscorer.git/ && "
        "python3 seqQscorer.py --spec None --assay None --rt None --raw ../../{input.raw} --map ../../{input.map} --loc ../../{input.loc} --tss ../../{input.tss} --out ../../{output} > ../../{log} 2>&1 && "
        "cd ../.."

# ------------------------------------------------------------------------------
rule aggregate_scores:
    input:
        get_samples_scores_files_by_analysis_and_GEOSeries
    output:
        "output/{analysis}/scores/{dataid}.scores.txt"
    group:
        "qualityFeatures"
    log:
        "output/logs/{analysis}/qc/{dataid}.aggregate_scores.log"
    shell:
        "cat {input} > {output} 2> {log}"

# ------------------------------------------------------------------------------
rule picard_rrna:
    input:
        gtf=config["ANNOT_PATH"],
        bam=rules.bowtie.output.bam
    output:
        temp("output/qc/{dataid}/picard/{dataid}.{sample}.rRNA.bed")
    log:
        "output/logs/qc/{dataid}.{sample}.picard_rrna.log"
    shell:
        """samtools view -H {input.bam} > {output} && """
        """zcat {input.gtf} | grep rRNA | awk 'BEGIN{{OFS="\\t"}}{{print $1, $4-1, $5, $7, NR}}' >> {output} 2> {log}"""

rule picard:
    input:
        bam=rules.bowtie.output.bam,
        bed=rules.picard_rrna.output,
        ref=config["PICARD_REFFLAT_FILE"],
    output:
        RNA_Metrics="output/qc/{dataid}/picard/{sample}.RNA_Metrics",
        ALN_Metrics="output/qc/{dataid}/picard/{sample}.ALN_Metrics"
    log:
        "output/logs/qc/{dataid}.{sample}.picard.log"
    threads:
        config["SALMON_INDEX_THREADS"]
    shell:
        "PicardCommandLine CollectRnaSeqMetrics I={input.bam} O={output.RNA_Metrics} REF_FLAT={input.ref} STRAND=NONE RIBOSOMAL_INTERVALS={input.bed} > {log} 2>&1 && "
        "PicardCommandLine CollectAlignmentSummaryMetrics INPUT={input.bam} OUTPUT={output.ALN_Metrics} >> {log} 2>&1"

# ------------------------------------------------------------------------------
rule multiqc:
    input:
        # get_picard_files_by_GEOSeries,
        get_bowtie2_files_by_GEOSeries,
        get_fastqc_files_by_GEOSeries
    output:
        "output/qc/{dataid}/multiqc/{dataid}_multiqc_report.html"
    log:
        "output/logs/qc/{dataid}.multiqc.log"
    params:
        outdir="output/qc/{dataid}/multiqc",
        # picarddir="output/qc/{dataid}/picard",
        bowtiedir="data/datasets/{dataid}/bowtie2",
        fastqcdir="output/qc/{dataid}/fastqc"
    shell:
        "multiqc -f -i {wildcards.dataid} -o {params.outdir} {params.bowtiedir} {params.fastqcdir} > {log} 2>&1"