# =============================================================================
# Common Functions for the Snakemake workflow
# =============================================================================
import pandas as pd
import numpy as np
import glob
from os import path

# Get sample file paths
# ------------------------------------------------------------------------------

def get_main_samples_by_dataid(wildcards):
    return(SAMPLES_MAIN_TABLE[SAMPLES_MAIN_TABLE.GEO_Series.eq(wildcards.dataid)])

# def get_batch_samples_by_dataid(wildcards):
#     return(SAMPLES_BATCH_TABLE[SAMPLES_BATCH_TABLE.GEO_Series.eq(wildcards.dataid)])

# def get_all_samples_by_dataid(wildcards):
#     return(SAMPLES_ALL_TABLE[SAMPLES_ALL_TABLE.GEO_Series.eq(wildcards.dataid)])


# qualityFeatures
# ------------------------------------------------------------------------------
def get_samples_scores_files_by_analysis_and_GEOSeries(wildcards):
    sampleids=SAMPLES_MAIN_TABLE
    # if wildcards.analysis == "batched":
    #     sampleids=SAMPLES_BATCH_TABLE
    sampleids = sampleids[sampleids.GEO_Series.eq(wildcards.dataid)][["GEO_Series", "Run"]]
    sampleids = sampleids.apply(lambda x : "output/qc/{}/scores/{}.score.txt".format(x[0], x[1]), axis=1)
    return(np.unique(sampleids))


# expressionAnalysis
# ------------------------------------------------------------------------------
def get_samples_quant_files_by_GEOSeries(wildcards):
    table=SAMPLES_MAIN_TABLE
    # if wildcards.analysis == "batched":
    #     table=SAMPLES_BATCH_TABLE
    table = table[table.GEO_Series.eq(wildcards.dataid)][["GEO_Series", "Run"]]
    sampleids = table.apply(lambda x : "data/output/{}/salmon/{}/quant.sf".format(x[0], x[1]), axis=1)
    return(np.unique(sampleids))

def get_samples_cts_files_by_analysis_and_GEOSeries(wildcards):
    table=SAMPLES_MAIN_TABLE
    # if wildcards.analysis == "batched":
    #     table=SAMPLES_BATCH_TABLE
    table = table[table.GEO_Series.eq(wildcards.dataid)][["GEO_Series", "Run"]]
    sampleids = table.apply(lambda x : "data/output/{}/salmon/{}_cts.tsv".format(x[0], x[1]), axis=1)
    return(np.unique(sampleids))

def get_samples_tpm_files_by_analysis_and_GEOSeries(wildcards):
    table=SAMPLES_MAIN_TABLE
    # if wildcards.analysis == "batched":
    #     table=SAMPLES_BATCH_TABLE
    table = table[table.GEO_Series.eq(wildcards.dataid)][["GEO_Series", "Run"]]
#    sampleids = sampleids.apply(lambda x : "output/{}/RASflowResults/{}/trans/tpmFile/{}_tpm.tsv".format(wildcards.analysis, x[0], x[1]), axis=1)
    sampleids = table.apply(lambda x : "data/output/{}/salmon/{}_tpm.tsv".format(x[0], x[1]), axis=1)
    return(np.unique(sampleids))

def getControlGroup(wildcards):
    data = pd.read_csv("config/metadata/Datasets.csv")
    return(data[data.GEO_Series.eq(wildcards.dataid)]["Control"].to_string(index=False).strip())

def getTreatGroup(wildcards):
    data = pd.read_csv("config/metadata/Datasets.csv")
    return(data[data.GEO_Series.eq(wildcards.dataid)]["Treat"].to_string(index=False).strip())

def getEndType(wildcards):
    data = pd.read_csv("config/metadata/Datasets.csv")
    if(data[data.GEO_Series.eq(wildcards.dataid)]["LibraryLayout"].to_string(index=False).strip() == "SINGLE"):
      return("single")
    else:
      return("pair")

def getSamplesPairing(wildcards):
    data = pd.read_csv("config/metadata/Datasets.csv")
    if(data[data.GEO_Series.eq(wildcards.dataid)]["SamplesPairing"].to_string(index=False).strip() == "1"):
      return("TRUE")
    else:
      return("FALSE")

# MultiQC
# ------------------------------------------------------------------------------
# def get_multiqc_files(wildcards):
#     data = pd.read_csv("config/metadata/Mega_SraRunTable.csv")
#     datasets = pd.read_csv("config/metadata/Datasets.csv")[["GEO_Series", "select"]]
#     data = data.merge(datasets, how="inner", on="GEO_Series")[["GEO_Series", "Selected", "select"]]
#     data = data[data.select.eq(1)]
#     GEO_Series_ids = data[data.Selected.eq(1)][["GEO_Series"]]
#     GEO_Series_ids = GEO_Series_ids.apply(lambda x : 'output/qc/{}/multiqc/{}_multiqc_report.html'.format(x[0],x[0]), axis=1)
#     return(np.unique(GEO_Series_ids))
# 
# def get_batched_multiqc_files(wildcards):
#     data = pd.read_csv("config/metadata/Mega_SraRunTable.csv")
#     datasets = pd.read_csv("config/metadata/Datasets.csv")[["GEO_Series", "batches"]]
#     data = data.merge(datasets, how="inner", on="GEO_Series")[["GEO_Series", "BatchAnalysisSelection", "batches"]]
#     data = data[data.batches.eq(1)]
#     GEO_Series_ids = data[data.BatchAnalysisSelection.eq(1)][["GEO_Series"]]
#     GEO_Series_ids = GEO_Series_ids.apply(lambda x : 'output/qc/{}/multiqc/{}_multiqc_report.html'.format(x[0],x[0]), axis=1)
#     return(np.unique(GEO_Series_ids))

# def get_multiqc_salmon_files(wildcards):
#     data = pd.read_csv("config/metadata/Mega_SraRunTable.csv")
#     datasets = pd.read_csv("config/metadata/Datasets.csv")[["GEO_Series", "select"]]
#     data = data.merge(datasets, how="inner", on="GEO_Series")[["GEO_Series", "Selected", "select"]]
#     data = data[data.select.eq(1)]
#     GEO_Series_ids = data[data.Selected.eq(1)][["GEO_Series"]]
#     GEO_Series_ids = GEO_Series_ids.apply(lambda x : 'output/main/qc/{}/multiqc/{}_multiqc_report.html'.format(x[0],x[0]), axis=1)
#     return(GEO_Series_ids)

# def get_batched_multiqc_salmon_files(wildcards):
#     data = pd.read_csv("config/metadata/Mega_SraRunTable.csv")
#     datasets = pd.read_csv("config/metadata/Datasets.csv")[["GEO_Series", "batches"]]
#     data = data.merge(datasets, how="inner", on="GEO_Series")[["GEO_Series", "BatchAnalysisSelection", "batches"]]
#     data = data[data.batches.eq(1)]
#     GEO_Series_ids = data[data.BatchAnalysisSelection.eq(1)][["GEO_Series"]]
#     GEO_Series_ids = GEO_Series_ids.apply(lambda x : 'output/batched/qc/{}/multiqc/{}_multiqc_report.html'.format(x[0],x[0]), axis=1)
#     return(GEO_Series_ids)

# Picard
# ------------------------------------------------------------------------------
def get_picard_files_by_GEOSeries(wildcards):
    table=SAMPLES_MAIN_TABLE
    # if wildcards.analysis == "batched":
    #     table=SAMPLES_BATCH_TABLE
    table = table[table.GEO_Series.eq(wildcards.dataid)][["GEO_Series", "Run"]]
    sampleids = table.apply(lambda x : "output/qc/{}/picard/{}.RNA_Metrics".format(x[0], x[1]), axis=1)
    return(np.unique(sampleids))

# def get_picard_files_by_GEOSeries(wildcards):
#     data = pd.read_csv("config/metadata/Mega_SraRunTable.csv")
#     datasets = pd.read_csv("config/metadata/Datasets.csv")[["GEO_Series", "select"]]
#     data = data.merge(datasets, how="inner", on="GEO_Series")[["GEO_Series", "Run", "Selected", "select"]]
#     data = data[data.select.eq(1)]
#     data = data[data.Selected.eq(1)]
#     #sampleids = data[data.Selected.eq(1)][["GEO_Series", "Run"]]
#     sampleids = data[data.GEO_Series.eq(wildcards.dataid)][["GEO_Series","Run"]]
#     sampleids = sampleids.apply(lambda x : 'output/qc/{}/picard/{}.RNA_Metrics'.format(x[0],x[1]), axis=1)
#     return(sampleids)
# 
# def get_batched_picard_files_by_GEOSeries(wildcards):
#     data = pd.read_csv("config/metadata/Mega_SraRunTable.csv")
#     datasets = pd.read_csv("config/metadata/Datasets.csv")[["GEO_Series", "batches"]]
#     data = data.merge(datasets, how="inner", on="GEO_Series")[["GEO_Series", "Run", "BatchAnalysisSelection", "batches"]]
#     data = data[data.batches.eq(1)]
#     data = data[data.Selected.eq(1)]
#     #sampleids = data[data.Selected.eq(1)][["GEO_Series", "Run"]]
#     sampleids = data[data.GEO_Series.eq(wildcards.dataid)][["GEO_Series","Run"]]
#     sampleids = sampleids.apply(lambda x : 'output/qc/{}/picard/{}.RNA_Metrics'.format(x[0],x[1]), axis=1)
#     return(sampleids)


# Bowtie2
# ------------------------------------------------------------------------------
# def get_bowtie2_files_by_GEOSeries(wildcards):
#     data = pd.read_csv("config/metadata/Mega_SraRunTable.csv")
#     datasets = pd.read_csv("config/metadata/Datasets.csv")[["GEO_Series", "select"]]
#     data = data.merge(datasets, how="inner", on="GEO_Series")[["GEO_Series", "Run", "Selected", "select"]]
#     data = data[data.select.eq(1)]
#     data = data[data.Selected.eq(1)]
#     #sampleids = data[data.Selected.eq(1)][["GEO_Series", "Run"]]
#     sampleids = data[data.GEO_Series.eq(wildcards.dataid)][["GEO_Series","Run"]]
#     sampleids = sampleids.apply(lambda x : 'data/datasets/{}/bowtie2/{}.bam'.format(x[0],x[1]), axis=1)
#     return(sampleids)

def get_bowtie2_files_by_GEOSeries(wildcards):
    table=SAMPLES_MAIN_TABLE
    # if wildcards.analysis == "batched":
    #     table=SAMPLES_BATCH_TABLE
    table = table[table.GEO_Series.eq(wildcards.dataid)][["GEO_Series", "Run"]]
    sampleids = table.apply(lambda x : "data/datasets/{}/bowtie2/{}.bam".format(x[0], x[1]), axis=1)
    return(np.unique(sampleids))

  
# Fastqc
# ------------------------------------------------------------------------------
# def get_fastqc_files_by_GEOSeries(wildcards):
#     data = pd.read_csv("config/metadata/Mega_SraRunTable.csv")
#     datasets = pd.read_csv("config/metadata/Datasets.csv")[["GEO_Series", "select"]]
#     data = data.merge(datasets, how="inner", on="GEO_Series")[["GEO_Series", "Run", "Selected", "select"]]
#     data = data[data.select.eq(1)]
#     data = data[data.Selected.eq(1)]
#     #sampleids = data[data.Selected.eq(1)][["GEO_Series", "Run"]]
#     sampleids = data[data.GEO_Series.eq(wildcards.dataid)][["GEO_Series","Run"]]
#     sampleids = sampleids.apply(lambda x : 'output/qc/{}/fastqc/{}.fastqc.txt'.format(x[0],x[1]), axis=1)
#     return(sampleids)

def get_fastqc_files_by_GEOSeries(wildcards):
    table=SAMPLES_MAIN_TABLE
    # if wildcards.analysis == "batched":
    #     table=SAMPLES_BATCH_TABLE
    table = table[table.GEO_Series.eq(wildcards.dataid)][["GEO_Series", "Run"]]
    sampleids = table.apply(lambda x : "output/qc/{}/fastqc/{}.fastqc.txt".format(x[0], x[1]), axis=1)
    return(np.unique(sampleids))


# quality vs expression
# ------------------------------------------------------------------------------
def get_q_vs_expr_output_files(wildcards):
    table=DATASETS[DATASETS.select.eq(1)][["GEO_Series"]]
    # table = table[["GEO_Series"]]
    file_paths = table.apply(lambda x : "output/{}/qualityVsExp/{}/{}.cor001.png".format(wildcards.analysis, x[0], x[0]), axis=1)
    return(np.unique(file_paths))

# DEseq2
# ------------------------------------------------------------------------------
def get_deseq2_output_files(wildcards):
    table=DATASETS[DATASETS.select.eq(1)][["GEO_Series"]]
    # table = table[["GEO_Series"]]
    file_paths = table.apply(lambda x : "output/{}/gene_expression/{}/{}.diff.genes.tsv.gz".format(wildcards.analysis, x[0], x[0]), axis=1)
    return(np.unique(file_paths))