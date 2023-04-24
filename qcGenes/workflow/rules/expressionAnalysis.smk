# =============================================================================
# Snakemake workflow to perform gene expression analysis
# =============================================================================

# ------------------------------------------------------------------------------
rule expressionAnalysis:
    input:
#        get_multiqc_salmon_files,
#        get_batched_multiqc_salmon_files,
#        get_all_pipeline_ack_files,
#        get_all_batched_pipeline_ack_files
#        "config/metadata/Mega_SraRunTable.csv",
#        "config/metadata/Datasets.csv"
        expand("output/main/RASflowResults/{dataid}/trans/tpmFile/all.samples.tpm.tsv", dataid=DATAIDS_MAIN),
        expand("output/batched/RASflowResults/{dataid}/trans/tpmFile/all.samples.tpm.tsv", dataid=DATAIDS_BATCH)

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

# ------------------------------------------------------------------------------
rule create_pipeline_dir:
    input:
        rules.salmon_index.output,
	get_samples_tpm_files_by_analysis_and_GEOSeries
    output:
        "output/{analysis}/pipelines/{dataid}/main.py"
    group:
        "expressionAnalysis"
    log:
        "output/logs/{analysis}/expressionAnalysis/{dataid}.create_pipeline_dir.log"
    params:
        control=getControlGroup,
        treat=getTreatGroup,
        trim="no",
        end=getEndType,
        samplesPairing=getSamplesPairing,
        threads=config["MAX_THREADS"],
        salmon_threads=config["RASFLOW_THREADS"],
	salmon_index_dir=config["SALMON_INDEX_DIR"],
        genomPath=config["GENOM_PATH"],
        annotPath=config["ANNOT_PATH"],
        transPath=config["TRANS_PATH"],
        ensidPath=config["ENSID_PATH"],
        ensemblDataset="hsapiens_gene_ensembl",
        subdir="{analysis}"
    conda:
        "../envs/linux.yaml"
    shell:
        "rm -rf output/{params.subdir}/pipelines/{wildcards.dataid}/ && "
        "git clone -q lib/RASflow.git output/{params.subdir}/pipelines/{wildcards.dataid} > {log} 2>&1 && "
        "rm -rf output/{params.subdir}/pipelines/{wildcards.dataid}/data output/{params.subdir}/pipelines/{wildcards.dataid}/output output/{params.subdir}/pipelines/{wildcards.dataid}/configs/metadata.tsv && "
        "ln -s ../../../../data output/{params.subdir}/pipelines/{wildcards.dataid}/ >> {log} 2>&1 && "
	"rm -rf output/{params.subdir}/pipelines/{wildcards.dataid}/data/output/{wildcards.dataid}/transcripts_index >> {log} 2>&1 && "
	"mkdir -p output/{params.subdir}/pipelines/{wildcards.dataid}/data/output/{wildcards.dataid} >> {log} 2>&1 && "
	"ln -s ../../../{params.salmon_index_dir} output/{params.subdir}/pipelines/{wildcards.dataid}/data/output/{wildcards.dataid}/transcripts_index >> {log} 2>&1 && "
	"rm -rf output/{params.subdir}/RASflowResults/{wildcards.dataid}/trans/quant >> {log} 2>&1 && "
	"mkdir -p output/{params.subdir}/RASflowResults/{wildcards.dataid}/trans >> {log} 2>&1 && "
	"ln -s ../../../../../data/output/{wildcards.dataid}/salmon output/{params.subdir}/RASflowResults/{wildcards.dataid}/trans/quant >> {log} 2>&1 && "
        "cp output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.template.yaml output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
        "sed -i \"s/VAL_PROJECT_NAME/{wildcards.dataid}/\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
        "sed -i \"s/VAL_TRIMMED/{params.trim}/\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
        "sed -i \"s#VAL_READSPATH#data/datasets/{wildcards.dataid}#\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
        "sed -i \"s#VAL_METAFILE#configs/metadata.tsv#\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
        "sed -i \"s/VAL_END/{params.end}/\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
        "sed -i \"s/VAL_NCORE_SALMON/{params.salmon_threads}/\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
        "sed -i \"s/VAL_NCORE/{params.threads}/\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
        "sed -i \"s#VAL_FINALOUTPUT#../../RASflowResults#\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
        "sed -i \"s#VAL_GENOME#{params.genomPath}#\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
        "sed -i \"s#VAL_ANNOTATION#{params.annotPath}#\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
        "sed -i \"s#VAL_BIOMART_ENS_IDS#{params.ensidPath}#\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
        "sed -i \"s#VAL_TRANS#{params.transPath}#\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
        "sed -i \"s/VAL_PAIR/{params.samplesPairing}/\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
        "sed -i \"s/VAL_CONTROL/['{params.control}']/\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
        "sed -i \"s/VAL_TREAT/['{params.treat}']/\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
        "sed -i \"s/VAL_EnsemblDataSet/{params.ensemblDataset}/\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 "

# ------------------------------------------------------------------------------
rule create_metadata_file:
    input:
#        "output/{analysis}/pipelines/{dataid}/main.py"
        rules.create_pipeline_dir.output
    output:
        "output/{analysis}/pipelines/{dataid}/configs/metadata.tsv"
    group:
        "expressionAnalysis"
    log:
        "output/logs/{analysis}/expressionAnalysis/{dataid}.create_metadata_file.log"
    conda:
        "../envs/tidyverse.yaml"
    shell:
        "Rscript workflow/scripts/create_metadata_file.R {wildcards.dataid} {wildcards.analysis} {output} > {log} 2>&1"

#    run:
#        """
#        data = pd.read_csv("config/metadata/Mega_SraRunTable.csv")
#        data = data.rename(columns={"Run": "sample"})
#	if params[1] == "main":
#	    data["my_selection"] = data["Selected"]
#	elif params[1] == "batched":
#	    data["my_selection"] = data["BatchAnalysisSelection"]
#	else:
#	    data["my_selection"] = data["Selected"]
#        data = data[data.my_selection.eq(1) & data.GEO_Series.eq(params[0])][["sample", "group", "subject", "batch"]]
#        data.to_csv(output[0], sep='\t', index=False)
#        """

# ------------------------------------------------------------------------------
rule run_pipeline:
    input:
#        "output/main/pipelines/{dataid}/main.py",
#        "output/main/pipelines/{dataid}/configs/metadata.tsv",
#	get_samples_fastq_files_by_GEOSeries,
        rules.create_pipeline_dir.output,
	rules.create_metadata_file.output,
#	get_samples_tpm_files_by_GEOSeries,
#        expand("output/{analysis}/RASflowResults/{dataid}/trans/tpmFile/{sample}_tpm.tsv", analysis="{analysis}", dataid="{dataid}", sample=get_samples_by_analysis_and_dataid),
        get_samples_tpm_files_by_analysis_and_GEOSeries,
        config["GENOM_PATH"],
        config["ANNOT_PATH"],
        config["TRANS_PATH"]
    output:
        protected("output/{analysis}/RASflowResults/{dataid}/trans/dea/countGroup/tx2gene.RData")
    group:
        "expressionAnalysis"
    log:
        "output/logs/{analysis}/expressionAnalysis/{dataid}.run_pipeline.log"
    benchmark:
        "output/logs/{analysis}/expressionAnalysis/{dataid}.run_pipeline.bench"
    conda:
        "../envs/rasflow.yaml"
#    threads:
#        config['MAX_THREADS']
    shell:
        "cd output/{wildcards.analysis}/pipelines/{wildcards.dataid} && (python main.py) > ../../../../{log} 2>&1 && cd ../../../.."

# ------------------------------------------------------------------------------
rule post_pipeline:
    input:
#        "output/main/pipelines/{dataid}/main.py",
#        "output/main/pipelines/{dataid}/configs/metadata.tsv",
#        "output/main/RASflowResults/{dataid}/trans/dea/countGroup/tx2gene.RData"
        rules.create_pipeline_dir.output,
	rules.create_metadata_file.output,
	rules.run_pipeline.output
    output:
        "output/{analysis}/RASflowResults/{dataid}/trans/tpmFile/all.samples.tpm.tsv"
    group:
        "expressionAnalysis"
    log:
        "output/logs/{analysis}/expressionAnalysis/{dataid}.post_pipeline.log"
    params:
        "output/{analysis}/pipelines/{dataid}"
    benchmark:
        "output/logs/{analysis}/expressionAnalysis/{dataid}.post_pipeline.bench"
    conda:
        "../envs/postExpression.yaml"
    threads:
        1
    shell:
        "cd {params} && "
	    "Rscript scripts/merge_tpm.R     > ../../../../{log} 2>&1 && "
        "Rscript scripts/merge_deseq2.R >> ../../../../{log} 2>&1 && cd ../../../.."


# ==============================================================================
# BATCH ANALYSIS
# ==============================================================================

# ------------------------------------------------------------------------------
#rule batchedExpressionAnalysis:
#    input:
#        get_all_batched_pipeline_ack_files
#        expand("output/batched/RASflowResults/{dataid}/trans/tpmFile/all.samples.tpm.tsv", dataid=DATAIDS_BATCH)

# ------------------------------------------------------------------------------
#rule create_batched_pipeline_dir:
#    input:
#        rules.salmon_index.output
#    output:
#        "output/batched/pipelines/{dataid}/main.py"
#    group:
#        "expressionAnalysis"
#    log:
#        "output/logs/batched/expressionAnalysis/{dataid}.create_pipeline_dir.log"
#    params:
#        control=getControlGroup,
#        treat=getTreatGroup,
#        trim="no",
#        end=getEndType,
#        samplesPairing=getSamplesPairing,
#        threads=config["RASFLOW_THREADS"],
#        salmon_threads=config["SALMON_THREADS"],
#	salmon_index_dir=config["SALMON_INDEX_DIR"],
#        genomPath=config["GENOM_PATH"],
#        annotPath=config["ANNOT_PATH"],
#        transPath=config["TRANS_PATH"],
#        ensidPath=config["ENSID_PATH"],
#        ensemblDataset="hsapiens_gene_ensembl",
#        subdir="batched"
#    conda:
#        "../envs/linux.yaml"
#    shell:
#        "rm -rf output/{params.subdir}/pipelines/{wildcards.dataid}/ && "
#        "git clone -q lib/RASflow.git output/{params.subdir}/pipelines/{wildcards.dataid} > {log} 2>&1 && "
#        "rm -rf output/{params.subdir}/pipelines/{wildcards.dataid}/data output/{params.subdir}/pipelines/{wildcards.dataid}/output output/{params.subdir}/pipelines/{wildcards.dataid}/configs/metadata.tsv && "
#        "ln -s ../../../../data output/{params.subdir}/pipelines/{wildcards.dataid}/ >> {log} 2>&1 && "
#	"rm -rf output/{params.subdir}/pipelines/{wildcards.dataid}/data/output/{wildcards.dataid}/transcripts_index >> {log} 2>&1 && "
#	"mkdir -p output/{params.subdir}/pipelines/{wildcards.dataid}/data/output/{wildcards.dataid} >> {log} 2>&1 && "
#	"ln -s ../../../{params.salmon_index_dir} output/{params.subdir}/pipelines/{wildcards.dataid}/data/output/{wildcards.dataid}/transcripts_index >> {log} 2>&1 && "
#        "cp output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.template.yaml output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
#        "sed -i \"s/VAL_PROJECT_NAME/{wildcards.dataid}/\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
#        "sed -i \"s/VAL_TRIMMED/{params.trim}/\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
#        "sed -i \"s#VAL_READSPATH#data/datasets/{wildcards.dataid}#\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
#        "sed -i \"s#VAL_METAFILE#configs/metadata.tsv#\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
#        "sed -i \"s/VAL_END/{params.end}/\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
#        "sed -i \"s/VAL_NCORE_SALMON/{params.salmon_threads}/\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
#        "sed -i \"s/VAL_NCORE/{params.threads}/\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
#        "sed -i \"s#VAL_FINALOUTPUT#../../RASflowResults#\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
#        "sed -i \"s#VAL_GENOME#{params.genomPath}#\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
#        "sed -i \"s#VAL_ANNOTATION#{params.annotPath}#\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
#        "sed -i \"s#VAL_BIOMART_ENS_IDS#{params.ensidPath}#\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
#        "sed -i \"s#VAL_TRANS#{params.transPath}#\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
#        "sed -i \"s/VAL_PAIR/{params.samplesPairing}/\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
#        "sed -i \"s/VAL_CONTROL/['{params.control}']/\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
#        "sed -i \"s/VAL_TREAT/['{params.treat}']/\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 && "
#        "sed -i \"s/VAL_EnsemblDataSet/{params.ensemblDataset}/\" output/{params.subdir}/pipelines/{wildcards.dataid}/configs/config_main.yaml >> {log} 2>&1 "

# ------------------------------------------------------------------------------
#rule create_batched_metadata_file:
#    input:
#        "output/batched/pipelines/{dataid}/main.py",
#        rules.create_batched_pipeline_dir.output
#    output:
#        "output/batched/pipelines/{dataid}/configs/metadata.tsv"
#    group:
#        "expressionAnalysis"
#    log:
#        "output/logs/batched/expressionAnalysis/{dataid}.create_metadata_file.log"
#    params:
#        dataid="{dataid}"
#    run:
#        data = pd.read_csv("config/metadata/Mega_SraRunTable.csv")
#        data = data.rename(columns={"Run": "sample"})
#        data = data[data.BatchAnalysisSelection.eq(1) & data.GEO_Series.eq(params[0])][["sample", "group", "subject", "batch"]]
#        data.to_csv(output[0], sep='\t', index=False)

# ------------------------------------------------------------------------------
#rule run_batched_pipeline:
#    input:
#        "output/batched/pipelines/{dataid}/main.py",
#        "output/batched/pipelines/{dataid}/configs/metadata.tsv",
#	get_batched_samples_fastq_files_by_GEOSeries,
#        rules.create_batched_pipeline_dir.output,
#        rules.create_batched_metadata_file.output,
#        get_batched_samples_tpm_files_by_GEOSeries,
#        config["GENOM_PATH"],
#        config["ANNOT_PATH"],
#        config["TRANS_PATH"]
#    output:
#        "output/batched/RASflowResults/{dataid}/trans/dea/countGroup/tx2gene.RData"
#    group:
#        "expressionAnalysis"
#    log:
#        "output/logs/batched/expressionAnalysis/{dataid}.run_pipeline.log"
#    benchmark:
#        "output/logs/batched/expressionAnalysis/{dataid}.run_pipeline.bench"
#    conda:
#        "../envs/rasflow.yaml"
#    threads:
#        config['MAX_THREADS']
#    shell:
#        "cd output/batched/pipelines/{wildcards.dataid} && (python main.py <<< y) > ../../../../{log} 2>&1 && cd ../../../.."

# ------------------------------------------------------------------------------
#rule post_batched_pipeline:
#    input:
#        "output/batched/pipelines/{dataid}/main.py",
#        "output/batched/pipelines/{dataid}/configs/metadata.tsv",
#        "output/batched/RASflowResults/{dataid}/trans/dea/countGroup/tx2gene.RData"
#        rules.create_batched_pipeline_dir.output,
#        rules.create_batched_metadata_file.output,
#	rules.run_batched_pipeline.output
#    output:
#        "output/batched/RASflowResults/{dataid}/trans/tpmFile/all.samples.tpm.tsv"
#    group:
#        "expressionAnalysis"
#    log:
#        "output/logs/batched/expressionAnalysis/{dataid}.post_pipeline.log"
#    params:
#        "output/batched/pipelines/{dataid}"
#    benchmark:
#        "output/logs/batched/expressionAnalysis/{dataid}.post_pipeline.bench"
#    conda:
#        "../envs/postExpression.yaml"
#    threads:
#        1
#    shell:
#        "cd {params} && "
#        "Rscript scripts/merge_tpm.R     > ../../../../{log} 2>&1 && "
#        "Rscript scripts/merge_deseq2.R >> ../../../../{log} 2>&1 && cd ../../../.."
