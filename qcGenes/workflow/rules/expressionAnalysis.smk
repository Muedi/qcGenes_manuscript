# =============================================================================
# Snakemake workflow to perform gene expression analysis
# =============================================================================

# ------------------------------------------------------------------------------
rule expressionAnalysis:
    input:
        expand("data/output/{dataid}/salmon/{sampleid}_cts.tsv", zip, dataid=SAMPLES_MAIN_GEO_COL, sampleid=SAMPLES_MAIN_RUN_COL),
        expand("data/output/{dataid}/salmon/{sampleid}_tpm.tsv", zip, dataid=SAMPLES_MAIN_GEO_COL, sampleid=SAMPLES_MAIN_RUN_COL),
        expand("output/main/gene_expression/{dataid}/{dataid}.samples.tpm.tsv.gz", dataid=DATAIDS_MAIN),
        expand("output/main/gene_expression/{dataid}/{dataid}.samples.cts.tsv.gz", dataid=DATAIDS_MAIN),
        expand("output/main/gene_expression/{dataid}/{dataid}.diff.genes.tsv.gz", dataid=DATAIDS_MAIN)

        
rule quantification:
    input:
        # expand("data/output/{dataid}/salmon/{sample}/quant.sf", zip, dataid=SAMPLES_ALL_GEO_COL, sample=SAMPLES_ALL_RUN_COL),
        expand("data/output/{dataid}/salmon/{sample}_tpm.tsv", zip, dataid=SAMPLES_ALL_GEO_COL, sample=SAMPLES_ALL_RUN_COL)

# rule pipeline:
#     input:
#         expand("output/main/pipelines/{dataid}/main.py", dataid=DATAIDS_MAIN)

# rule metadata:
#     input:
#         expand("output/main/pipelines/{dataid}/configs/metadata.tsv", dataid=DATAIDS_MAIN)

# rule run_all_pipelines:
#     input:
#         expand("output/main/RASflowResults/{dataid}/trans/dea/countGroup/tx2gene.RData", dataid=DATAIDS_MAIN)

# # ------------------------------------------------------------------------------
# rule unzip_gtf:
#     input:
#         gtfgz = config["ANNOT_PATH"]
#     output:
#         gtf = temp(config["ANNOT_PATH_UNZIPPED"])
#     group:
#         "expressionAnalysis"
#     log:
#         "output/logs/expression/expressionAnalysis/unzip_gtf.log"
#     benchmark:
#         "output/logs/expression/expressionAnalysis/unzip_gtf.bench"
#     shell:
#         """
#         gzip -d -c {input.gtfgz} > {output.gtf} 2> {log}
#         """

# ------------------------------------------------------------------------------
rule salmon_index:
    input:
        trans = config["TRANS_PATH"],
        genom = config["GENOM_PATH"]
    output:
        index = directory(config["SALMON_INDEX_DIR"]),
	info  = os.path.join(config["SALMON_INDEX_DIR"], "info.json")
    group:
        "expressionAnalysis"
    log:
        "output/logs/expression/expressionAnalysis/salmon_index.log"
    benchmark:
        "output/logs/expression/expressionAnalysis/salmon_index.bench"
    conda:
        "../envs/salmon.yaml"
    threads:
        config["SALMON_INDEX_THREADS"]
    params:
        decoy=os.path.join(config["SALMON_INDEX_DIR"], "decoys.txt"),
	gentrome=config["GENTROME_PATH"]
    shell:
        "mkdir -p {output.index} > {log} 2>&1 && "
        "grep '^>' <(gunzip -c {input.genom}) | cut -d ' ' -f 1 > {params.decoy} 2>> {log} && "
	"sed -i.bak -e 's/>//g' {params.decoy} >> {log} 2>&1 && "
	"cat {input.trans} {input.genom} > {params.gentrome} 2>> {log} 2>&1 && "
        "salmon index -t {params.gentrome} -d {params.decoy} -i {output.index} -p {threads} >> {log} 2>&1"
#        "salmon index -t {input.trans} -i {output.index} --type quasi -k 31 -p {threads} > {log} 2>&1"

# # ------------------------------------------------------------------------------
# rule quantify_sample:
#     input:
#         read1="data/datasets/{dataid}/{sampleid}_1.fastq.gz",
#         index=rules.salmon_index.output.index,
#         info=rules.salmon_index.output.info
#     output:
#          tpm="data/output/{dataid}/salmon/{sampleid}_tpm.tsv"
#     group:
#         "expressionAnalysis"
#     log:
#         "output/logs/expression/{dataid}.{sampleid}.quantify_sample.log"
#     benchmark:
#         "output/logs/expression/{dataid}.{sampleid}.quantify_sample.bench"
#     conda:
#         "../envs/salmon.yaml"
#     threads:
#         config['SALMON_THREADS']
#     params:
#         quant_dir="data/output/{dataid}/salmon/{sampleid}",
#         read2="data/datasets/{dataid}/{sampleid}_2.fastq.gz",
# 	gtf=config["ANNOT_PATH"],
#         max_reads=config['SALMON_MAX_READS']
#     shell:
#         """
#         if [ -f "{params.read2}" ]; then
#             salmon quant -l A -1 {input.read1} -2 {params.read2} -o {params.quant_dir} -i {input.index} -g {params.gtf} -p {threads} --seqBias --gcBias --posBias --reduceGCMemory 2> /dev/null
#             awk 'NR==1{{next}}{{print $1"\\t"$4}}' {params.quant_dir}/quant.sf > {output.tpm}
#         else 
#             salmon quant -l A -r {input.read1} -o {params.quant_dir} -i {input.index} -g {params.gtf} -p {threads} --seqBias --gcBias --posBias --reduceGCMemory 2> /dev/null
#             awk 'NR==1{{next}}{{print $1"\\t"$4}}' {params.quant_dir}/quant.sf > {output.tpm}
#         fi
#         """


# # ------------------------------------------------------------------------------
rule quantify_sample:
    input:
        read1="data/datasets/{dataid}/{sampleid}_1.fastq.gz",
        index=rules.salmon_index.output.index,
        info=rules.salmon_index.output.info,
	gtf=config["ENSIDVER_PATH"]
    output:
         tpm="data/output/{dataid}/salmon/{sampleid}_tpm.tsv"
#        quant_dir=directory("output/{analysis}/RASflowResults/{dataid}/trans/quant/{sampleid}"),
#        tpm="output/{analysis}/RASflowResults/{dataid}/trans/tpmFile/{sampleid}_tpm.tsv"
    group:
        "expressionAnalysis"
    log:
        "output/logs/expression/{dataid}.{sampleid}.quantify_sample.log"
    benchmark:
        "output/logs/expression/{dataid}.{sampleid}.quantify_sample.bench"
    conda:
       "../envs/salmon.yaml"
    threads:
       config['SALMON_THREADS']
    params:
        quant_dir="data/output/{dataid}/salmon/{sampleid}",
        read2="data/datasets/{dataid}/{sampleid}_2.fastq.gz",
	    tmp_dir="data/output/{dataid}/reads",
        tmp_fq_1="data/output/{dataid}/reads/{sampleid}.1.fastq.gz",
        tmp_fq_2="data/output/{dataid}/reads/{sampleid}.2.fastq.gz",
#        gtf=config["ANNOT_PATH"],
        max_reads=config['SALMON_MAX_READS']
    shell:
        """
        mkdir -p {params.tmp_dir}
        if [ {params.max_reads} -gt 0 ] ; then seqtk sample -s 100 {input.read1} {params.max_reads} 2> {log} | gzip -1 > {params.tmp_fq_1} 2>> {log}; fi 
        if [ {params.max_reads} -lt 0 ] ; then cp -f {input.read1} {params.tmp_fq_1} >> {log} 2>&1; fi
        if [ -f "{params.read2}" ]; then
            if [ {params.max_reads} -gt 0 ] ; then seqtk sample -s 100 {params.read2} {params.max_reads} 2>> {log} | gzip -1 > {params.tmp_fq_2} 2>> {log}; fi
            if [ {params.max_reads} -lt 0 ] ; then cp -f {params.read2} {params.tmp_fq_2} >> {log} 2>&1; fi
#            salmon quant -l A -1 {params.tmp_fq_1} -2 {params.tmp_fq_2} -o {params.quant_dir} -i {input.index} -g {input.gtf} -p {threads} --seqBias --gcBias --posBias >> {log} 2> /dev/null
            salmon quant -l A -1 {params.tmp_fq_1} -2 {params.tmp_fq_2} -o {params.quant_dir} -i {input.index} -g {input.gtf} -p {threads} --seqBias --gcBias --posBias --reduceGCMemory >> {log} 2> /dev/null
#            salmon quant -l A -1 {input.read1} -2 {params.read2} -o {params.quant_dir} -i {input.index} -g {input.gtf} -p {threads} --seqBias --gcBias --posBias --reduceGCMemory >> {log} 2> /dev/null
            awk 'NR==1{{next}}{{print $1"\\t"$4}}' {params.quant_dir}/quant.sf > {output.tpm}
        else 
#            salmon quant -l A -r {params.tmp_fq_1} -o {params.quant_dir} -i {input.index} -g {input.gtf} -p {threads} --seqBias --gcBias --posBias >> {log} 2> /dev/null
            salmon quant -l A -r {params.tmp_fq_1} -o {params.quant_dir} -i {input.index} -g {input.gtf} -p {threads} --seqBias --gcBias --posBias --reduceGCMemory >> {log} 2> /dev/null
#            salmon quant -l A -r {input.read1} -o {params.quant_dir} -i {input.index} -g {input.gtf} -p {threads} --seqBias --gcBias --posBias --reduceGCMemory 2> /dev/null
            awk 'NR==1{{next}}{{print $1"\\t"$4}}' {params.quant_dir}/quant.sf > {output.tpm}
        fi
        rm -f {params.tmp_fq_1} {params.tmp_fq_2} >> {log} 2>&1
        """


rule extract_counts_from_quant:
    input:
        rules.quantify_sample.output.quant
    output:
        "data/output/{dataid}/salmon/{sampleid}_cts.tsv"
    group:
        "expressionAnalysis"
    log:
         "output/logs/expressionAnalysis/{dataid}.{sampleid}.extract_counts_from_quant.log"
    threads:
        1
    shell:
        """awk 'NR==1{{next}}{{print $1"\\t"$5}}' {input} > {output} 2> {log}"""

rule extract_tpm_from_quant:
    input:
        rules.quantify_sample.output.quant
    output:
        "data/output/{dataid}/salmon/{sampleid}_tpm.tsv"
    group:
        "expressionAnalysis"
    log:
         "output/logs/expressionAnalysis/{dataid}.{sampleid}.extract_tpm_from_quant.log"
    threads:
        1
    shell:
        """awk 'NR==1{{next}}{{print $1"\\t"$4}}' {input} > {output} 2> {log}"""


# ------------------------------------------------------------------------------
rule merge_tpm:
    input:
        ids=config["ENSIDFULLVER_PATH"],
        files=get_samples_tpm_files_by_analysis_and_GEOSeries
    output:
        "output/{analysis}/gene_expression/{dataid}/{dataid}.samples.tpm.tsv.gz"
    group:
        "expressionAnalysis"
    log:
         "output/logs/{analysis}/expressionAnalysis/{dataid}.merge_tpm.log"
    threads:
        1
    shell:
        "Rscript workflow/scripts/merge_tpm.R {wildcards.dataid} {input.ids} {output} {input.files} > {log} 2>&1"


# ------------------------------------------------------------------------------
rule merge_cts:
    input:
        ids=config["ENSIDFULLVER_PATH"],
        files=get_samples_cts_files_by_analysis_and_GEOSeries
    output:
        "output/{analysis}/gene_expression/{dataid}/{dataid}.samples.cts.tsv.gz"
    group:
        "expressionAnalysis"
    log:
         "output/logs/{analysis}/expressionAnalysis/{dataid}.merge_cts.log"
    threads:
        1
    shell:
        "Rscript workflow/scripts/merge_cts.R {wildcards.dataid} {input.ids} {output} {input.files} > {log} 2>&1"


# ------------------------------------------------------------------------------
rule deseq2:
    input:
        ids=config["ENSIDVER_PATH"],
        scores=rules.aggregate_scores.output,
        files=get_samples_quant_files_by_GEOSeries
    output:
        rlog="output/{analysis}/gene_expression/{dataid}/{dataid}.samples.rlg.tsv.gz",
        diff="output/{analysis}/gene_expression/{dataid}/{dataid}.diff.genes.tsv.gz"
    group:
        "expressionAnalysis"
    log:
         "output/logs/{analysis}/expressionAnalysis/{dataid}.diff_analyze.log"
    threads:
        config["DESEQ2_THREADS"]
    params:
        datasets=config["DATASETS_FILE"]
    shell:
        #"Rscript workflow/scripts/diff_analyze.R {wildcards.dataid} {input.datasets} {input.samples} {input.ids} {output.rlog} {output.diff} {input.files} > {log} 2>&1"
        "Rscript workflow/scripts/diff_analyze.R {wildcards.dataid} {params.datasets} {input.ids} {input.scores} {output.rlog} {output.diff} {input.files} > {log} 2>&1"