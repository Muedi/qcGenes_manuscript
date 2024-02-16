# =============================================================================
# Snakemake workflow to download all datasets
# =============================================================================

# ------------------------------------------------------------------------------
rule downloadDatasets:
    input:
        expand("data/datasets/{dataid}/{sample}_1.fastq.gz", zip, dataid=SAMPLES_ALL_GEO_COL, sample=SAMPLES_ALL_RUN_COL)

# ------------------------------------------------------------------------------
rule download_sample:
    output:
        "data/datasets/{dataid}/{sample}_1.fastq.gz"
    group:
        "downloadDatasets"
    log:
        "output/logs/main/datasets/{dataid}/{sample}.download_sample.log"
    benchmark:
        "output/logs/main/datasets/{dataid}/{sample}.download_sample.bench"
    conda:
        "../envs/sratools.yaml"
    params:
        max_reads=config['MAX_READS']
    shell:
        "fastq-dump -X {params.max_reads} --outdir data/datasets/{wildcards.dataid} --gzip --skip-technical --split-files {wildcards.sample} > {log} 2>&1"
