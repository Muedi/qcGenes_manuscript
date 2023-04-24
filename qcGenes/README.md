# qcGenes: identifying quality-biased genes in RNA-seq datasets (beta!)

Several RNA-seq datasets and selected samples are described in a metadata table
that is processed by several Snakemake workflows to download, process and 
analyse the data. This is a beta version (not stable). 

---

## Changes

### Current to do list

* Major
  * Data
    * add more datasets
  * Analysis
    * Top quality-correlated genes: show additional stats (e.g. # of datasets where differential and biased)
    * compare with batch corrected data (https://www.bioconductor.org/packages/release/bioc/html/sva.html)
    * PCA after correction (incl. or not outliers removal), predicted batch vs quality, ...
  * Snakemake
    * post_pipeline rule for batched datasets
    * !! create and publish docker image
    * compare bowtie results (% mapping) on different fastq file sizes (1M, 2M, 10M ...)
  * RASflow:
    * create index with decoy sequence, see: https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
    * should raise error if failed (exit code)
* Minor
* Low-priority
  * Summary in PDF (may be difficult to have a single md file for both html and pdf output)
  * limit parallelization of dataset downloads (gnu-parallel?)

### Changes log
  * 2021-12-09
    * DONE Dataset Quality vs Size plot (cor.size.png): include dataset with minimal number of samples
    * DONE Plot Dataset Quality vs AvgSpotLength
    * DONE Use a Salmon index that already exists (avoid creation for each dataset)
    * DISCARDED call directly sub workflows instead of main.py (https://snakemake.readthedocs.io/en/v5.16.0/snakefiles/modularization.html#snakefiles-sub-workflows)
    * DONE by choosing generic model => tune seqQscorer parameters using Dataset description (assay type and run type)
    * DISCARDED: Linear regressions: add plot showing substration of DIF from ALL, and DIS from ALL
    * DISCARDED: define config files to have data input and output dir under same user defined root directory (merge data and output dirs)

  * 2020-10-23 (and earlier)
    * support paired designs (PAIR: FALSE by default in config/config.yaml)
    * set package versions in conda envs yaml files (tidyverse and readsAnno)
    * include human and mouse Biomart examples + do not use them anymore
  * 2020-07-13
    * Declined (overkill + not better than clustering evaluation on many datasets): compare outlier removal to random removal of same number of samples (n times) (aggregate results only for datasets with many samples)
  * 2020-06-09
    * Batched-datasets analysis: rules and new column in Mega metadata for file selection
    * Reorganising output dir (main and batched subdirs); moving pipelines dir to output dir
    * Declined (not practical if other features have to be computed): rename features files with raw.txt, map.txt, loc.txt and tss.txt extensions
  * 2020-06-05
    * Finalise report (figure captions)
    * P vs Batch plot: geom_jitter() instead?
    * add column in Datasets.tsv to easily include/exclude a dataset from the analysis (only for final aggregation analysis)
    * More datasets with batch annotation
    * Batch vs P_low: plots + correlation
    * Bowtie memory benchmarks
  * 2020-05-30
    * Optimize Bowtie2 run time and memoty imprint (manual benchmark)
    * replace ack files by real output files
  * 2020-05-29
    * add support of paired-end datasets (see Wrappers)
    * write input function to replace dependency to Mega Metadata file
    * SeqQscorer could compute Bowtie stats on 1M reads only
    * Bigger text size in figures
  * 2020-05-15
    * Rmarkdown Summary
  * 2020-05-13
    * Collect and plot correlation coefficients (e.g. quality vs deg FDR) (rule qualityCorGenes)
  * 2020-05-12
    * Compute enrichment of disease (and also differentially-expressed) genes in 2 gene sets: quality-correlated and non-quality-correlated (will summarise the related plots .cor001.png) => use scatter plots and correlation coefficients
  * 2020-05-11
    * Fix bug on plot Quality vs Diff genes top 10% that shows genes in extreme categories where no differential genes should be (see plots above) (qualityVsExp e.g. GSE108643) 
    * 2020-05-06
    * Draft report
    * conda env with gnu-linux tools (wget, tar, ...)
  * 2020-05-05
    * set groups for rules
    * Rule createPipelineDir: use quiet version (progression in log file)
  * 2020-05-04
    * rename rules files with .smk extension
    * commands to run on the Mogon cluster
    * revise download options (see https://edwards.sdsu.edu/research/fastq-dump/) => added --skip-technical
  * 2020-04-30
    * correlation calculated with or without outliers? => Keep all samples
    * make conda env with linux tools for rules without specific env (linux.yaml wget). => Done but not all tools available.
  * 2020-04-29
    * implement log files in RASflow's rules to reduce scatter on terminal (e.g. ChIPseeker summary on stdout, and RASflow)
    * use correlation of sample's group and sample's quality to enable dataset selection/filtering
    * (Failed) Change getReads rule to use original data path or symbolic links instead of scp
  * 2020-04-27
    * Analysis / Analyse aggregated quality vs expression statistics (e.g. correlation, dunn, entropy)
    * RASflow / generate expression table with all samples and all genes (or transcripts)
    * delete or save plot produced by ChIPpeakAnno.binOverFeature function to appropriate directory (e.g. use functions pdf / dev.off)
  * 2020-04-18: revise usage of ack files for acknowledging steps of the pipeline (use real output files in rules' input/output)
    * see directory() and touch() directives at https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#log-files

---

## Requirements

### Software

* Linux OS (e.g. Ubuntu 16.04)
* Anaconda or Miniconda installed (e.g. conda 4.8.3)
* git (e.g. 2.7.4)

### Storage

* qcGenes source code: 1-2 MB
* github dependencies: 1.2 GB
* conda environments: 10.1 GB
* reference files: 7.8 GB
* Example analysis of 17 datasets (226 fastq files of 2M reads): 
  * downloaded data: 37 GB
  * analysis output dir: 6 GB
  * RASflow pipelines: 225 MB

---

## Installation

### Installation with Docker as Root user
The host computer is your "normal" computer (here linux). You must install Docker to create and use containers that are like light-weight virtual machines on the host. You will use a terminal on the host to install and start containers. You will use a terminal within a container for computations using pre-installed software in the container.

#### Preparations
```bash
# Create 2 full paths to dirs on the host (outside docker)
DIR=/absolute/path/to/input_output/folder # to be changed to your needs
HOST_DATA_DIR=$DIR/data                   # For the downloaded data
HOST_OUTPUT_DIR=$DIR/output               # For computing results
HOST_SNAKELOG_DIR=$DIR/snakemake_logs     # For Snakemake logs

mkdir -p $HOST_DATA_DIR
mkdir -p $HOST_OUTPUT_DIR
mkdir -p $HOST_SNAKELOG_DIR

# Create parameter strings for docker cli
MAP_DIR1="-v $HOST_DATA_DIR:/home/qcgenes/data"
MAP_DIR2="-v $HOST_OUTPUT_DIR:/home/qcgenes/output"
MAP_DIR3="-v $HOST_SNAKELOG_DIR:/home/qcgenes/.snakemake/log"

# Create parameter strings for snakemake cli (set more cores if possible)
OPT1="--latency-wait 15 --restart-times 2 --use-conda --conda-frontend mamba --cores 4"
OPT2="--latency-wait 15 --restart-times 2 --use-conda --conda-frontend mamba --cores 8"
```

#### Create Docker image if not available for download
```bash
# Clone Main GIT Repo
git clone --branch docker.v1 git@github.com:JFK24/qcGenes.git qcGenes.git
cd qcGenes.git

# Clone and set GIT dependencies
git clone --branch docker.v1 git@github.com:JFK24/seqQscorer.git lib/seqQscorer.git
git clone --branch docker.v1 git@github.com:JFK24/RASflow.git lib/RASflow.git
cd lib/RASflow.git
git checkout mainz
cd ../..

# Build docker image (including required Snakemake's conda envs)
docker build -t jfk24/qcgenesenvs .
```

#### Use interactive terminal within Docker container
```bash
# Start interactive bash session in container
docker run -ti $MAP_DIR1 $MAP_DIR2 $MAP_DIR3 jfk24/qcgenesenvs bash

# We are now inside the container!!

# Test manually if conda envs have been properly created (>2 files in subfolders)
ll .snakemake/conda/*

# Snakemake parameters for use inside container
OPT1="--latency-wait 15 --restart-times 2 --use-conda --conda-frontend mamba --cores 4"
OPT2="--latency-wait 15 --restart-times 2 --use-conda --conda-frontend mamba --cores 8"

# OPTION 1 = Run the full pipeline with a defined number of cores (MAY NOT WORK / BROKEN LINKS BETWEEN RULES)
snakemake $OPT2

# OPTION 2 = Run the pipeline with less cores for download and more for computations
snakemake $OPT1 downloadRefFiles downloadDatasets && snakemake $OPT2
```

#### Use Non-Interactive Docker container
```bash
# OPTION 1 = Run the full pipeline with defined number of cores (set more cores if possible) (MAY NOT WORK / BROKEN LINKS BETWEEN RULES)
docker run $MAP_DIR1 $MAP_DIR2 $MAP_DIR3 jfk24/qcgenesenvs snakemake $OPT1

# OPTION 2 = Run the pipeline with less cores for download and more for computations (set more cores if possible)
docker run $MAP_DIR1 $MAP_DIR2 $MAP_DIR3 jfk24/qcgenesenvs /bin/bash -c "snakemake $OPT1 downloadRefFiles downloadDatasets && snakemake $OPT1"
```

#### Use custom config files and use Docker container
```bash
# ABSOLUTE PATHS ON THE HOST TO CUSTOMIZED COPIES OF 3 CONFIG FILES (CHANGE TO YOUR NEEDS)
HOST_CONFIG=$DIR/config/host.config.yaml
HOST_DATASETS=$DIR/config/host.Datasets.csv
HOST_SRA_TABLE=$DIR/config/host.Mega_SraRunTable.csv

# Create parameter strings for docker cli
MAP_FILE1="-v $HOST_CONFIG:/home/qcgenes/config/config.yaml"
MAP_FILE2="-v $HOST_DATASETS:/home/qcgenes/config/metadata/Datasets.csv"
MAP_FILE3="-v $HOST_SRA_TABLE:/home/qcgenes/config/metadata/Mega_SraRunTable.csv"

# START INTERACTIVE DOCKER SESSION
docker run -ti $MAP_DIR1 $MAP_DIR2 $MAP_DIR3 $MAP_FILE1 $MAP_FILE2 $MAP_FILE3 jfk24/qcgenesenvs bash

# START NON INTERACTIVE DOCKER SESSION
docker run $MAP_DIR1 $MAP_DIR2 $MAP_DIR3 $MAP_FILE1 $MAP_FILE2 $MAP_FILE3 jfk24/qcgenesenvs /bin/bash -c "snakemake $OPT1 downloadRefFiles downloadDatasets && snakemake $OPT1"
```

#### Some useful docker commands

##### Start a stopped container 
```bash
docker start -ai  <CONTAINER ID>
```

##### Start a stopped interactive container as a non interactive one with command 
```bash
CONTAINER_ID=e1a58550d3d4
NEW_IMAGE_NAME=test2

docker commit $CONTAINER_ID $NEW_IMAGE_NAME
docker run -v $HOST_DATA_DIR:/home/qcgenes/data -v $HOST_OUTPUT_DIR:/home/qcgenes/output $NEW_IMAGE_NAME snakemake $OPT2
# docker run -v $HOST_DATA_DIR:/home/qcgenes/data -v $HOST_OUTPUT_DIR:/home/qcgenes/output $NEW_IMAGE_NAME "snakemake $OPT1 downloadRefFiles downloadDatasets && snakemake $OPT2" # NOT WORKING

# Remove image and container when finished
docker image rm $NEW_IMAGE_NAME
#docker rm <NEW CONTAINER> 

```


### Installation as normal user (no docker)
```bash
# Clone Main Repo
git clone git@github.com:JFK24/qcGenes.git qcGenes.git
cd qcGenes.git

# Clone dependencies
git clone git@github.com:JFK24/seqQscorer.git lib/seqQscorer.git
# git clone -b mainz git@ github.com:JFK24/RASflow.git lib/RASflow.git # to be tested
git clone git@github.com:JFK24/RASflow.git lib/RASflow.git
cd lib/RASflow.git
git checkout mainz
cd ../..

# Create the conda environments
conda install -y -c conda-forge mamba
mamba create -y -c conda-forge -c bioconda -n qcgenes snakemake=5.26.1
#conda activate qcgenes
#conda install -y -c conda-forge graphviz=2.42.3 # for visualizing the DAG
#conda install -y -c conda-forge cookiecutter=1.7.2 # for using slurm profile

# OPTIONAL: Only for testing some tidyverse scripts (e.g. rasflow)
mamba create -y -c r -n tidyverse r-tidyverse=1.2.1 r-yaml r-factominer r-factoextra r-ggpubr r-ggrepel r-fpc
```

---

## Configuration

Note that changing config files may make snakemake rerun a large part of the pipeline. 

| File                                 | Description                                                                                                                    |
| ------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------ |
| config/config.yaml                   | Local paths and variables (e.g. max reads in Fastq files)                                                                      |
| config/metadata/Mega_SraRunTable.csv | Selected samples from RNA-seq experiments to process (SRA repository, important columns are Selected, group, subject, and Run) | 
| config/metadata/Datasets.csv         | Description of the Datasets detailed in Mega_SraRunTable.csv (e.g. groups of samples)                                          |

### config/config.yaml

* Bowtie index: See example config files for human and mouse
* Bowtie index must be compatible with R packages used in rule reads_all_anno (see also R scripts)

```bash
# HUMAN BOWTIE INDEX CONFIG
BOWTIE_IDX_DIR: data/ref/bowtie2
BOWTIE_IDX_PATH: data/ref/bowtie2/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz
BOWTIE_IDX_FILE1: data/ref/bowtie2/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.1.bt2
BOWTIE_IDX_URL: ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz
BOWTIE_IDX: data/ref/bowtie2/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index

# MOUSE BOWTIE INDEX CONFIG
BOWTIE_IDX_DIR: data/ref/bowtie2
BOWTIE_IDX_PATH: data/ref/bowtie2/mm10.zip
BOWTIE_IDX_FILE1: data/ref/bowtie2/mm10.1.bt2
BOWTIE_IDX_URL: https://genome-idx.s3.amazonaws.com/bt/mm10.zip
BOWTIE_IDX: data/ref/bowtie2/mm10
```

---

## Run the analysis

### Download data separately
This is the only way to restrict the simultaneous number of downloads. 
If not required, you can run the all pipeline directly (see below)

```bash
# Define sankemake options and activate conda environment.
# You may change the number of cores as available (after --cores)
OPT="--latency-wait 15 --restart-times 2 --use-conda --conda-frontend mamba --cores 4"
conda activate qcgenes
# Download reference files (7.8 GB)
snakemake $OPT downloadRefFiles
# Get datasets (e.g. 100 files 2M reads = 7.4GB)
snakemake $OPT downloadDatasets
```
 
### Run the pipeline 
Use following command for first-time analysis (it will also download the data) 
or after the separate downloads. 
It can be used after selecting new files in the Mega metadata file. 
But modifying already selected files in the metadata will request to force a rerun. 
```bash
snakemake $OPT
```

### Create the report
```bash
snakemake --report report.html
sed -i "s/max-width:800px/max-width:90%/" report.html
```

---

## Run on Mogon

### Adapt some settings
```bash
sed -i "s/MAX_READS: 2000000/MAX_READS: 10000000/" config/config.yaml
sed -i "s/MAX_THREADS: 4/MAX_THREADS: 8/" config/config.yaml
sed -i "s/BOWTIE2_THREADS: 4/BOWTIE2_THREADS: 8/" config/config.yaml
```

### Download datasets
```bash
# 1. Only download datasets (8 threads)
SBATCH_OPT='-J SMK -o SMK_DW.log -A jgu-cbdm -p andrade -c 8 --mem=8G --time=24:00:00 --mail-type=FAIL,END --mail-user=jfontain@uni-mainz.de'
OPT="--latency-wait 15 --restart-times 3 --use-conda --cores 8"
sbatch $SBATCH_OPT --wrap="snakemake $OPT --snakefile workflow/rules/downloadDatasets.smk"
```

### Run as full-node job
```bash
# 2. Full workflow (1 full 64-core node)
SBATCH_OPT='-J SMK -o SMK.log -A jgu-cbdm -p andrade -N 1 -c 64 --mem=110G --time=48:00:00 --mail-type=FAIL,END --mail-user=jfontain@uni-mainz.de'
OPT="--latency-wait 15 --restart-times 3 --use-conda --cores 64"
sbatch $SBATCH_OPT --dependency=singleton --wrap="snakemake $OPT" 
```

Report can be then created as described above.

---

## Results

Structure of the output directory:
```
- output
-- main              # main analysis
--- qualityCorGenes/ # quality correlated genes across datasets
--- qualityFeatures/ # quality features for each dataset (SeqQscorer tool)
--- qualityVsExp/    # plots and tables of correlated genes for each dataset
--- RASflowResults/  # gene expression analysis for each dataset (RASflow tool)
--- scores/          # quality scores of the samples for each dataset
--- statistics/      # statistical summaries
-- batched           # partly similar to main but limited to datasets with batch
```

---

## Notes

### Popular diseases

The script **popular_diseases_by_genes.r** generates a table of diseases sorted 
by the number of associated genes (source: GS2D). A copy of this table, 
including additional columns for manual notes, is in the **config/metadata** 
subdirectory to keep track of investigated datasets on GEO (**config/metadata/popular_diseases_by_genes annotated.tsv**).

### Adding a new dataset

* Search on GEO a human RNA-seq dataset with SRA Run selector link (note the **GEOSeriesID**, e.g. GSE12345)
* On SRA Run Selector page, check if the dataset use single-end RNA-seq (e.g. scroll down the first table or click on a RUN id)
* Download metadata file from SRA Run Selector save as **config/metadata/GEOSeriesID_SraRunTable.csv**
* Integrate manually this data to **config/metadata/Mega_SraRunTable.csv**, and fill first 3 columns to select samples and groups (2 different groups only for each dataset)
* Update manually the **config/metadata/Datasets.csv** file (1 row per dataset)
* Then you may run the following commands on a 4-core computer:

```bash
OPT="--latency-wait 15 --restart-times 2 --use-conda --conda-frontend mamba --cores 4"
conda activate qcgenes &&
snakemake $OPT downloadDatasets &&
snakemake $OPT &&
snakemake --report report.html &&
sed -i "s/max-width:800px/max-width:90%/" report.html

```


### Commands

```bash
# Secure downloaded data from accidental deletion
chmod a-w data/datasets/**/*.fastq.gz data/ref/

# Open all PCA plots
eog output/main/qualityVsExp/**/*.pca001.png

# Open all Scores plots
eog output/main/qualityVsExp/**/*.scores001.png

# Open all genes vs quality plots
eog output/main/qualityVsExp/**/*.cor001.png

# Open all volcano plots
gvfs-open output/main/RASflowResults/**/trans/dea/visualization/volcano*.pdf
xdg-open output/main/RASflowResults/GSE142582/trans/dea/visualization/volcano_plot_healthy_psoriasis.pdf 
gvfs-open output/main/RASflowResults/GSE142582/trans/dea/visualization/volcano_plot_healthy_psoriasis.pdf 
gio open output/main/RASflowResults/**/trans/dea/visualization/volcano*.pdf
evince output/main/RASflowResults/**/trans/dea/visualization/volcano*.pdf
gnome-open output/main/RASflowResults/**/trans/dea/visualization/volcano*.pdf

# Get data size
find data/datasets/ -type f -name '*.gz' -exec du -ch {} + | grep total$

# Clear output dir but not the downloads
rm -rf output data/output
rm -rf data/datasets/**/download_ack data/datasets/**/bowtie2

# Delete snakemake's conda envs (remove -n to really run this command)
snakemake -n --conda-cleanup-envs --cores 1
```

### Commands for batch analysis

```bash

# for id in GSE82177 GSE120099 GSE61491
for id in GSE61491
do
rm -rf output/main/RASflowResults/$id
rm -rf output/main/qualityVsExp/$id
rm -rf output/main/qualityFeatures/$id
rm -rf output/main/scores/$id.scores.txt
rm -rf output/main/pipelines/$id
snakemake $OPT 
done

eog output/main/qualityVsExp/**/*.batch.png
eog output/main/qualityVsExp/{GSE82177,GSE120099,GSE61491}/*.cor001.png
eog output/main/qualityVsExp/{GSE82177,GSE120099,GSE61491}/*.pca001.png
eog output/main/qualityVsExp/{GSE82177,GSE120099,GSE61491}/*.scores001.png


```

### Alternative runs

Remove -n or --dry-run option to really run the pipeline.


#### Global scope
```bash
# Rerun PCA plots from Quality Analysis
rm -rf output/main/qualityVsExp/*
snakemake $OPT --force qualityAnalysis

# Rerun PCA plots from Gene Expression pipelines
rm -f output/Rmain/ASflowResults/**/pipeline_run_ack
rm -f output/main/RASflowResults/**/trans/dea/visualization/*.pdf
snakemake $OPT expressionAnalysis

# Rerun all Gene Expression pipelines
# rm -rf output/RASflowResults/** pipelines/** output/qualityCorGenes
rm -rf output/main/RASflowResults/** pipelines/**
snakemake $OPT expressionAnalysis

# Runs the quality analysis for the aggregation of all datasets results
snakemake -n $OPT --force qualityAnalysis

```

#### On a selected dataset or sample
```bash
DATAID=GSE92724
SAMPLE=SRR5125028

# Clear all output related to a Dataset and rerun snakemake
rm -rf output/main/**/$DATAID output/main/scores/$DATAID.scores.txt output/main/qualityCorGenes/ data/output/$DATAID data/datasets/$DATAID/bowtie2 output/main/pipelines/$DATAID summary.html report.html
snakemake $OPT

# Runs gene expression pipeline for the selected dataset and quality analysis
rm -rf output/main/RASflowResults/$DATAID output/main/pipelines/$DATAID output/main/qualityVsExp/$DATAID output/main/qualityCorGenes/
snakemake -n $OPT

# Runs gene expression pipeline for the selected dataset
rm -rf output/main/RASflowResults/$DATAID output/main/pipelines/$DATAID
snakemake -n $OPT output/main/RASflowResults/$DATAID/pipeline_run_ack

# Runs quality vs expression  for the selected dataset
snakemake -n $OPT --force output/main/qualityVsExp/$DATAID/$DATAID.cor001.png

# Deleting several Gene Expression Analyses
for DATAID in GSE108643 GSE122619 GSE76220 GSE99816 GSE117875 GSE85567 GSE82177
do
rm -rf output/main/RASflowResults/$DATAID output/main/pipelines/$DATAID
done

# Generates quality features for the selected dataset
# rm output/qualityFeatures/GSE92724/scores/SRR5125028.score.txt
# rm output/scores/GSE92724.scores.txt
snakemake -n $OPT --force output/main/scores/$DATAID.scores.txt

# Runs the quality analysis for the selected dataset
snakemake -n $OPT --force output/main/qualityVsExp/$DATAID/$DATAID.cor001.png

# Runs Fastqc on a selected sample
snakemake -n $OPT --force output/main/qualityFeatures/$DATAID/fastqc/$SAMPLE.fastqc.txt

# Run Bowtie on a selected sample
snakemake -n $OPT --force data/datasets/$DATAID/bowtie2/$SAMPLE.bam

# Run Reads annot on a selected sample
snakemake -n $OPT --force output/main/qualityFeatures/$DATAID/tss_anno/$SAMPLE.tss.txt

# Run Scorer on a selected sample
snakemake -n $OPT --force output/main/qualityFeatures/$DATAID/scores/$SAMPLE.score.txt
```

#### Force all rules

Note that downloading again all the datasets, or re-running all the computations
may take a while.

```bash
snakemake -n $OPT --forceall --snakefile workflow/rules/downloadRefFiles.smk downloadRefFiles
snakemake -n $OPT --forceall --snakefile workflow/rules/downloadDatasets.smk downloadDatasets
snakemake -n $OPT --forceall qualityAnalysis
```

<!-- ----------------------------------------------------------------------- -->
### Slurm profile **[DO NOT WORK SO FAR]**
Create a config file (change as required MY_EMAIL, MY_ACCOUNT and MY_PARTITION):
```bash
WDIR=$(pwd)
EMAIL=MY_EMAIL
SLURM_ACCOUNT=MY_ACCOUNT
SLURM_PARTITION=MY_PARTITION
cat > config/cluster_config.yaml <<EOF
__default__:
  account: $MY_ACCOUNT
  partition: $MY_PARTITION
  time: 01:00:00
  mem: 4G

run_pipeline:
  cpus-per-task: 8
  mem: 16G

quality_vs_exp:
  mem: 16G
  
bowtie:
  cpus-per-task: 8
  mem: 8G

qualityAnalysis:
  mail-type=FAIL,END
  mail-user=$EMAIL

downloadDatasets:
  mail-type=FAIL,END
  mail-user=$EMAIL

downloadRefFiles:
  mail-type=FAIL,END
  mail-user=$EMAIL

expressionAnalysis:
  mail-type=FAIL,END
  mail-user=$EMAIL

qualityAnalysis:
  mail-type=FAIL,END
  mail-user=$EMAIL

qualityFeatures:
  mail-type=FAIL,END
  mail-user=$EMAIL
EOF
```
Download and install the profile:
```bash
# Install Slurm profile
mkdir -p ~/.config/snakemake/slurm
cp config/cluster_config.yaml ~/.config/snakemake/config.yaml
cd ~/.config/snakemake
rm -rf slurm

git clone git@github.com:Snakemake-Profiles/slurm.git slurm.git
slurm.git/
git checkout d33a404be9585eee8d8d9b7e0bdbc2670bd49722
cd ..

cp slurm.git/\{\{cookiecutter.profile_name\}\}/* slurm/

# sed -i "s///" slurm/slurm-submit.py
# sed -i "s///" slurm/slurm-submit.py
# sed -i "s///" slurm/slurm-submit.py
# 
# SBATCH_DEFAULTS = """{{cookiecutter.sbatch_defaults}}"""
# CLUSTER_CONFIG = "{{cookiecutter.cluster_config}}"
# ADVANCED_ARGUMENT_CONVERSION = {"yes": True, "no": False}["{{cookiecutter.advanced_argument_conversion}}"]
# 
# SBATCH_DEFAULTS = """"""
# CLUSTER_CONFIG = "../config.yaml"
# ADVANCED_ARGUMENT_CONVERSION = {"yes": True, "no": False}["no"]

# DO NOT WORK
cookiecutter https://github.com/Snakemake-Profiles/slurm.git <<EOF

slurm

../config.yaml

EOF
cd $WDIR

```

#### Run

#### Adapt some settings
```bash
sed -i "s/MAX_READS: 2000000/MAX_READS: 10000000/" config/config.yaml
sed -i "s/MAX_THREADS: 4/MAX_THREADS: 8/" config/config.yaml
sed -i "s/BOWTIE2_THREADS: 4/BOWTIE2_THREADS: 8/" config/config.yaml
```
#### Run across nodes as multiple batch jobs **[not working]**
```bash
cd my/scratch/working/dir
OPT="--latency-wait 15 --restart-times 3 --use-conda --cores 200  --immediate-submit --notemp --max-jobs-per-second 1 --max-status-checks-per-second 1"
snakemake $OPT --profile slurm

# OPT="--latency-wait 15 --restart-times 3 --use-conda --cores 160"
# nohup snakemake $OPT --profile slurm &

# snakemake $OPT --profile slurm downloadDatasets
# snakemake $OPT --profile slurm downloadRefFiles
```

#### Run without profile but on a full node as interactive job **[not recommended]**
```bash
srun --pty -p andrade -A jgu-cbdm --time=72:00:00 -N 1 --mem=110G bash -i
OPT="--latency-wait 15 --restart-times 2 --use-conda --cores 64"
cd my/scratch/working/dir
snakemake $OPT
```

---

<!-- ----------------------------------------------------------------------- -->
## Methods
<!-- ----------------------------------------------------------------------- -->

### PCA clustering validation (R package fpc)
* The Dunn Index is the ratio of the smallest distance between observations not in the same cluster to the largest intra-cluster distance (minimum separation / maximum diameter). The Dunn Index has a value between zero and ∞, and should be maximized. Dunn index, see Halkidi et al. (2002)
* pearsongamma: correlation between distances and a 0-1-vector where 0 means same cluster, 1
means different clusters. "Normalized gamma" in Halkidi et al. (2001)
  * Halkidi, M., Batistakis, Y., Vazirgiannis, M. (2001) On Clustering Validation Techniques, Journal
of Intelligent Information Systems, 17, 107-145.

---

## Benchmarks

### Bowtie2

* On a 4-core 3.30Ghz 4th-gen. Core i5, 16GB RAM 1333Mhz
* time snakemake $OPT --force data/datasets/GSE92724/bowtie2/SRR5125028.bam data/datasets/GSE92724/bowtie2/SRR5125029.bam data/datasets/GSE92724/bowtie2/SRR5125030.bam data/datasets/GSE92724/bowtie2/SRR5125031.bam

| Times            | RAM (15.6 max)    | Setting (seqtk 2M->1M reads) |
|------------------|-------------------|------------------------------|
| 1:34, 1:37, 1:37 | 2.14 -> 2.36GB    | bowtie2 (1 cores; -mm)       |
| 1:42, 1:44, 1:49 | 1.98 -> 2.07GB    | bowtie2 (2 cores; -mm)       |
| 1:55, 1:57, 1:58 | 1.98 -> 2.04GB    | bowtie2 (4 cores; -mm)       |
| 1:41, 1:47, 2:01 | 2.14 -> 15.1GB    | bowtie2 (1 cores)            |
| 1:44, 1:46, 2:16 | 2.14 -> 8.62GB    | bowtie2 (2 cores)            |
| 2:00, 2:02, 2:26 | 1.98 -> 5.24GB    | bowtie2 (4 cores)            |

## Troubleshouting

### FGSEA
```
Warning message:
Unknown or uninitialised column: `SYMBOL`. 
Warning message:
Unknown or uninitialised column: `SYMBOL`. 
Warning message:
In fgsea(pathways = pathways, stats = ranks, nperm = 1000) :
  There are ties in the preranked stats (98.91% of the list).
The order of those tied genes will be arbitrary, which may produce unexpected results.
```

### Datasets

#### GSE106899, GSE73189, GSE123496, GSE99816, GSE92724, GSE120852, GSE117875, GSE63744, ...
```
rule indexTrans:
    input: data/ref/GRCh38.99/Homo_sapiens.GRCh38.cdna.all.fa.gz
    output: data/output/GSE106899/trans/transcripts_index
    jobid: 22
[2020-04-25 22:52:36.229] [jointLog] [warning] Entry with header [ENST00000434970.2], had length less than the k-mer length of 31 (perhaps after poly-A clipping)
[2020-04-25 22:52:43.328] [jointLog] [warning] Removed 12034 transcripts that were sequence duplicates of indexed transcripts.
```


#### GSE142582

Review samples selection to choose only AvgSpotLen>200


### Fix Snakemake's conda.py file 

This problem should not occur if using a recent version of snakemake such as 
the version used above in the installation section. 

Otherwise, if you are using the following software:

* Snakemake 5.3.0
* conda 4.8.3

Snakemake 5.3.0 uses deprecated "source activate" command instead of "conda activate" (https://bitbucket.org/snakemake/snakemake/issues/1115/cannot-activate-conda-enironment-using) It fails to activate conda's environments from rules:

Methods described below may fix the problem.

#### Method 1 (tested successfully)

Run the following code to fix the problem.

```bash
conda activate rasflow
CONDASH=$(conda info --base)/etc/profile.d/conda.sh
SRCFILE=$(python <<EOF                              
from snakemake import conda
print(conda.__file__)
EOF)
sed "s|source activate|source $CONDASH && conda activate|" $SRCFILE
conda deactivate
```

```bash
# Example paths for $SRCFILE
echo $HOME/opt/miniconda3/envs/rasflow/lib/python3.6/site-packages/snakemake/conda.py
# Example paths for $CONDASH
echo $HOME/opt/miniconda3/etc/profile.d/conda.sh
```

#### Method 2 (not tested)

Put `shell.executable("/bin/bash")` at the beginning of the Snakefile 

#### Method 3 (not tested)

Install Snakemake 5.4.3 or newer (not tested)

---

## Bibliography / References (to annotate)

### Data

* GEO
* SRA
* http://sra.dbcls.jp/index.html (SRA search by keywords)
* ArrayExpress https://www.ebi.ac.uk/arrayexpress/help/programmatic_access.html

### Batch-annotated datasets

Search in GEO Datasets (2020-06-02): 

* batch effect 
* Filters: Organism=Homo Sapiens, Study type=Expression profiling by high throughput sequencing
* (batch[All Fields] AND effect[All Fields]) AND ("Homo sapiens"[Organism] AND "Expression profiling by high throughput sequencing"[Filter])

Datasets:

* GSE82177 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE82177
  * HCV Infection and HCC
  * 27 samples, 2 batches
* GSE115643 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115643
  * EVI1-mediated modulation of gene expression, HEK293 cells (H2O2 treatment)
  * 24 samples, 3 batches
* GSE120099 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120099
  * RNA-seq proiles of 3 stages of human SMC differentiation starting from iPSC
  * 92 samples, 3 batches
* GSE61491 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61491
  * RNA-seq in NPCs treated with shRNAs targeting CHD8. For controls, NPCs were treated with shRNAs targeting GFP and LacZ. Infection and sequencing was carried out in two separate batches, with one GFP and one LacZ sample in each batch. All samples were sequenced in two technical replicates.
  * 54 samples, 2 batches
  * replicate samples not distributed across batches

Unused datasets:

* GSE91019 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91019
  * Treatment of Detroit 562 cells with 5 stimuli: LPS, TNFα, Pam2CSK4, Poly I:C or M tri-DAP followed by gene expression analysis by RNA-seq. Two biological duplicates were analyzed for each condition.
  * 14 samples, 4 batches
  * 4 controls but 2 samples max for a treatment
* GSE111203 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111203
  * MCF7 cells treated by chemicals (PBDE, PPT, MPP)
  * 18 samples, batch 2 with only 1 sample
* GSE97471 Bad design
* GSE136340 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136340
  * Effects of darunavir upon gene expression in kidney tubular cells after transduction with HIV or EGFP-control lentivirus
  * 18 samples, 3 batches
* GSE136864 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136864
  * A549 cells were not-transduced or transduced with empty vector, ELF1-WT, or ELF1-R8A mutant. Indicated cultures were stimulated with interferon beta for 6 hours or 48 hours. All samples were harvested simultaneously and analyzed by RNA-Seq.
  * 17 samples, 3 batches
  * multiple factors
* GSE119345 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119345  
  * RNA was isolated from HuT 78 cells transfected with siRNA for integrins (si-ITGAV/B3); TRA (si-TR) or non-target sequences (si-CT) followed by treatment with bexarotene in the presence of physiological concentrations of TH, for 6 hours using TRI-REAGENT
  * 60 samples, 2 batches
  * sample groups not distributed across batches

### Genes related to quality and ML classification for quality

* Ma A., Zhu Z., Ye M., Wang F. (2019) EnsembleKQC: An Unsupervised Ensemble Learning Method for Quality Control of Single Cell RNA-seq Sequencing Data. In: Huang DS., Jo KH., Huang ZK. (eds) Intelligent Computing Theories and Application. ICIC 2019. Lecture Notes in Computer Science, vol 11644. Springer, Cham. https://doi.org/10.1007/978-3-030-26969-2_47
  * TO READ
* Tomislav Ilicic, Jong Kyoung Kim, Aleksandra A. Kolodziejczyk, Frederik Otzen Bagger, Davis James McCarthy, John C. Marioni & Sarah A. Teichmann. Classification of low quality cells from single-cell RNA-seq data. Genome Biology volume 17, Article number: 29 (2016). https://doi.org/10.1186/s13059-016-0888-1
  * "With microfluidic capture methods visual inspection under the microscope allows identification of wells containing broken, empty, and multiple cells to be found."
  * Unclear definition of specificity in the main text (better in methods). (Fig 4: "Sensitivity is defined as the proportion of correctly identified low quality cells. Specificity is defined as the proportion of correctly identified high quality cells. F-score is defined as the harmonic mean between sensitivity and specificity")
  * No precision shown or discussed (reverse engineered precision = 0.79 and 0.41 for top 2 best cases)
  * section Application to diverse cell types and protocols: 
  * ensemble classifier is better when using all features (not only common features) although those features are shown to be dataset-specific
  * ensemble classifier fails miserably on 1 dataset out of 3 (sensitivity on low-quality cells at 100% and sensitivity for high-quality cells at 0% = probably all cells were predicted low-quality)
  * Why is the following expected? "As expected, classification failed in other cell types and protocols since all cells are considered as high quality (zero sensitivity), due to training the model with cell-type specific features (Fig. 2e)"
  * How did they get annotation for published datasets missing prior annotation? "We further applied the PCA-based method on two datasets containing published human cancer cell lines without prior quality annotation"
  * "The mouse SVM model on human cancer cells performed best (65 % accuracy based on prior feature-based PCA annotation) when excluding genes relating to Cytoplasm as a feature. [...] This means that an SVM model trained on mouse cells cannot be directly applied to human cancer cell lines."
  * What is a "PCA-based version and our SVM model"?
  * R package: https://github.com/ti243/cellity 
  * Python pipeline: https://github.com/ti243/celloline

Reverse engineering precision for dataset Shalek 2014 (mentioned as working very well)
===========================================================================
paper: sen_low=84%, sen_high=90%
data (suppl. file): n_low=345, n_high=745

Stats to detect low-quality cells
TP = 345 * 0.84 = 289.8
FP = 745 * (1-0.90) = 74.5
TN = 745 * 0.90 = 670.5
FN = 345 * (1-0.84) = 55.2

sensitivity_low=0.84
precision_low=TP/(TP+FP) = 289.8 / (289.8+74.5) = 289.8 / 364.3 = **0.795**

Reverse engineering precision for dataset CD4+ T cells (mentioned as working very well)
===========================================================================

paper: sen_low=84%, sen_high=62%
data (suppl. file): n_low=449, n_high=1375 (total = 1824)

Stats to detect low-quality cells
TP = 449 * 0.84 = 377.16
FP = 1375 * (1-0.62) = 522.5
TN = 1375 * 0.62 = 852.5
FN = 449 * (1-0.84) = 71.84

sensitivity_low=0.84
precision_low=TP/(TP+FP) = 377.16 / (377.16+522.5) = 377.16 / 899.66 = **0.419**


### Gene co-expression across datasets

* Rhodes et al. (2004). ONCOMINE: A Cancer Microarray Database and Integrated Data-Mining Platform. https://doi.org/10.1016/s1476-5586(04)80047-2

* Adler et al. (2009). (MEM - Multi Experiment Matrix; data?). Mining for coexpression across hundreds of datasets using novel rank aggregation and visualization methods. https://doi.org/10.1186/gb-2009-10-12-r139

* Wolfe CJ, Kohane IS, Butte AJ: Systematic survey reveals general applicability of "guilt-by-association" within gene coexpression networks. BMC Bioinformatics. 2005, 6: 227 https://doi.org/10.1186/1471-2105-6-227

* Hughes TR, Marton MJ, Jones AR, Roberts CJ, Stoughton R, Armour CD, Bennett HA, Coffey E, Dai H, He YD, Kidd MJ, King AM, Meyer MR, Slade D, Lum PY, Stepaniants SB, Shoemaker DD, Gachotte D, Chakraburtty K, Simon J, Bard M, Friend SH: Functional discovery via a compendium of expression profiles. Cell. 2000, 102: 109-126. https://doi.org/10.1016/S0092-8674(00)00015-5

* Stuart JM, Segal E, Koller D, Kim SK: A gene-coexpression network for global discovery of conserved genetic modules. Science. 2003, 302: 249-255. https://doi.org/10.1126/science.1087447

* Wilson BJ, Giguère V: Identification of novel pathway partners of p68 and p72 RNA helicases through Oncomine meta-analysis. BMC Genomics. 2007, 8: 419 https://doi.org/10.1186/1471-2164-8-419

* Basso K, Margolin AA, Stolovitzky G, Klein U, Dalla-Favera R, Califano A: Reverse engineering of regulatory networks in human B cells. Nat Genet. 2005, 37: 382-390. https://doi.org/10.1038/ng1532

* Rhodes DR, Tomlins SA, Varambally S, Mahavisno V, Barrette T, Kalyana-Sundaram S, Ghosh D, Pandey A, Chinnaiyan AM: Probabilistic model of the human protein-protein interaction network. Nat Biotechnol. 2005, 23: 951-959. https://doi.org/10.1038/nbt1103

* Kemmeren P, van Berkum NL, Vilo J, Bijma T, Donders R, Brazma A, Holstege FCP: Protein interaction verification and functional annotation by integrated analysis of genome-scale data. Mol Cell. 2002, 9: 1133-1143. https://doi.org/10.1016/S1097-2765(02)00531-2

* Pennacchio LA, Loots GG, Nobrega MA, Ovcharenko I: Predicting tissue-specific enhancers in the human genome. Genome Res. 2007, 17: 201-211. https://doi.org/10.1101/gr.5972507

* Brazma A, Jonassen I, Vilo J, Ukkonen E: Predicting gene regulatory elements in silico on a genomic scale. Genome Res. 1998, 8: 1202-1215.

* Adler P, Peterson H, Agius P, Reimand J, Vilo J: Ranking genes by their co-expression to subsets of pathway members. Ann NY Acad Sci. 2009, 1158: 1-13. https://doi.org/10.1111/j.1749-6632.2008.03747.x

### Outlier detection

* Drobin et al.: Molecular Profiling for Predictors of Radiosensitivity in Patients with Breast or Head-and-Neck Cancer. Cancers (Basel). 2020 Mar; 12(3): 753. https://doi.org/10.3390/cancers12030753
  * Outlier detection based on PCA
  * Combat to integrate different datasets
  
### Batch effect correction

* Zhou et al.: Influence of batch effect correction methods on drug induced differential gene expression profiles. BMC Bioinformatics. 2019; 20: 437. https://doi.org/10.1186/s12859-019-3028-6
  * **correction has strong positive impact** (LIMMA was fine; LEAPP not effective) 
  * Limma: **sample size must be sufficient (>40)** and include 2-3 principal components
  * smaller total sample size: results should be interpreted with caution
* Conesa et al. A survey of best practices for RNA-seq data analysis. Genome Biol. 2016; 17: 13. https://doi.org/10.1186/s13059-016-0881-8
  * The **NOISeq R package [20]** contains a wide variety of diagnostic plots to identify sources of biases in RNA-seq data and to apply appropriate normalization procedures in each case. 
  * Finally, despite these sample-specific normalization methods, batch effects may still be present in the data. These **effects can be minimized by appropriate experimental design** [51] or, alternatively, removed by batch-correction methods such as **COMBAT [52] or ARSyN**   
* Chao et al.: Systematic evaluation of RNA-Seq preparation protocol performance. BMC Genomics. 2019; 20: 571. https://doi.org/10.1186/s12864-019-5953-1
  * To exclude the possibility that these differences stemmed from batch effects, such as different set of libraries being prepared at different times, **we included additional technical replicates, prepared at different times**, for the TruSeq Stranded Total RNA and mRNA protocols (1 μg).
  * Taken together, these results demonstrate that the **variability among these library preparation protocols was not primarily due to batch effects**.
  
* Soneson et al.: Batch Effect Confounding Leads to Strong Bias in Performance Estimates Obtained by Cross-Validation. PLoS One. 2014; 9(6): e100335. https://doi.org/10.1371/journal.pone.0100335
  * Simulated data from real gene expression data as the basis for our simulation
  * The current study focuses on the impact of batch effects on the ability to build and evaluate the performance of a classifier based on gene expression data. Construction of classifiers, with the aim to assign samples to groups or predict some other trait of interest, is one of the most common goal in gene expression studies.
  * Studies comparing several of these tools have led to the conclusion that the **performance of most approaches is similar [8]**, [25].
  * Batch effects can be a big obstacle when **combining data sets**, and their characterization and potential elimination have recently received much attention in the literature (e.g., [7], [8], [23], [25]).
  * We have shown that in the presence of a batch effect with **at least moderate level of confounding with the main grouping variable, the performance estimates obtained by cross-validation are highly biased**.
  * The presence of a **batch effect that was completely non-confounded with the signal of interest did not introduce bias** in the performance estimates obtained by cross-validation.
  * However, **the bias in the cross-validation performance estimates is not eliminated by the batch effect removal**, and consequently the cross-validation performance estimates obtained after batch effect elimination are not more reliable measures of the true performance than those obtained without batch effect elimination.
  * In other words, batch effect removal methods should not be trusted blindly as a ‘post-experimental’ way of rescuing a badly designed experiment.
  * We have shown that **in the presence of batch effects that are confounded with the signal of interest, many of the highly ranked variables are associated only with the batch effect and not truly differentially expressed between the interesting groups**.
  * Similarly, we expect the results to generalize to other batch effect removal methods and classifiers.
* Fatai and Gamieldien: A 35-gene signature discriminates between rapidly- and slowly-progressing glioblastoma multiforme and predicts survival in known subtypes of the cancer. BMC Cancer. 2018; 18: 377. https://doi.org/10.1186/s12885-018-4103-5
  * As gene expression of the TCGA samples was profiled in batches which could introduce bias in classification analysis [12, https://www.ncbi.nlm.nih.gov/pubmed/24967636/], the statistical significance of batch effect was assessed as a function of the selected genes using **guided Principal Component Analysis (gPCA) from the R package gPCA** [13, https://www.ncbi.nlm.nih.gov/pubmed/23958724/]. **The approach used by TCGA** (2008) [14, https://www.ncbi.nlm.nih.gov/pubmed/18772890/; **Distance Weighted Discrimination (DWD) method** is applied to data for batchcorrection (https://dx.doi.org/10.1093%2Fbioinformatics%2Fbts096).] and Verhaak et al. (2011) [15, https://www.ncbi.nlm.nih.gov/pubmed/20129251/] was employed to generate gene-centric expression data.

* Müller et al.: Removing Batch Effects from Longitudinal Gene Expression - Quantile Normalization Plus ComBat as Best Approach for Microarray Transcriptome Data. PLoS One. 2016; 11(6): e0156594. https://doi.org/10.1371/journal.pone.0156594
  * Gutenberg Health Study
  * Various correction methods on microarrays
  * biological variability as validation
* Larsen et al.: Microarray-Based RNA Profiling of Breast Cancer: Batch Effect Removal Improves Cross-Platform Consistency. Biomed Res Int. 2014; 2014: 651751. https://doi.org/10.1155/2014/651751
  * Batch adjustment was found to be particularly **valuable in the detection of more delicate differences** in gene expression. 
  * Furthermore, our results show that prober adjustment is **essential for integration of gene expression data obtained from multiple sources**. 
  * We show that **high-variance genes are highly reproducibly expressed across platforms** making them particularly well suited as biomarkers and for building gene signatures, exemplified by prediction of estrogen-receptor status and molecular subtypes
* Zhang et al. ComBat-Seq: batch effect adjustment for RNA-Seq count data. https://doi.org/10.1101/2020.01.13.904730
  * negative binomial regression
  * However, **batch effects in composition**, i.e. the level of expression of genes scaled by the total expression (coverage) in each sample, **cannot be fully corrected with normalization**. Leek et al. (2010) provided an example of composition batch effects in microarray data
  * For heterogeneity from **unknown sources, SVASeq (Leek, 2014) and RUVSeq** (Risso et al., 2014b) are commonly used.
  * For differential expression, many common methods or procedures (e.g. edgeR (Robinson et al., 2010) and DESeq2 (Love et al., 2014)) suggest to **include batch variables as covariates in the linear models** behind these methods to account for the impact of batch.
  * We used the **polyester** R package (Frazee et al., 2015) to simulate realistic RNA-Seq studies
  * More specifically, batch 1 contains 5 replicates of cells overexpressing HER2, and 12 replicates for GFP controls (GEO accession **GSE83083**); batch 2 contains 6 replicates of each for EGFR and its corresponding controls (GEO accession **GSE59765**); batch 3 consists of 9 replicates of each for wild type KRAS and GFP controls (GEO accession GSE83083).
  * In this case **[when data contain no dispersion batch effect], ComBat-Seq which assumes separate dispersion across batch may be redundant and lead to higher false positives**
  * ComBat-Seq controls false positives and shows benefits in increased true positive rates only when a true dispersion batch effect is present in the data
  * This is consistent with the intuition of batch effect adjustment, that **modifying the data in any way comes with a risk of jeopardizing biological signals in the data**
  * Therefore, **batch effects should only be adjusted when they are present and result in unfavorable impact on downstream analysis**
  *  We focused primarily on addressing the unwanted impact of batch effect on downstream differential expression
  * Our ComBat-Seq method is based on a gene-wise negative binomial regression model, which, similar to other (generalized) linear models, **may not work well on data with severely or even completely confounded study designs**. 
  * However, **batch correction in confounded designs is challenging for most if not all the state-of-the-art batch adjustment methods**, and **careful experimental design has been widely advised** to mitigate the unfavorable impact of batch effects.
  * SVA - Leek et al.: The sva package for removing batch effects and other unwanted variation in high-throughput experiments. Bioinformatics. 2012 Mar 15; 28(6): 882–883. https://doi.org/10.1093/bioinformatics/bts034
  * SVA tool. Leek JT, Storey JD (2007) Capturing Heterogeneity in Gene Expression Studies by Surrogate Variable Analysis. PLoS Genet 3(9): e161. https://doi.org/10.1371/journal.pgen.0030161
  * Removing batch effects and using surrogate variables in differential expression analysis **have been shown to reduce dependence, stabilize error rate estimates, and improve reproducibility** (see [2, 3, 4] for more detailed information).
  * Show (hidden) batch effect in 9 public data. Leek, J., Scharpf, R., Bravo, H. et al. Tackling the widespread and critical impact of batch effects in high-throughput data. Nat Rev Genet 11, 733–739 (2010). https://doi.org/10.1038/nrg2825
  * Papiez et al.: BatchI: Batch effect Identification in high-throughput screening data using a dynamic programming algorithm. Bioinformatics. 2019; 35(11): p1885-1892. https://doi.org/10.1093/bioinformatics/bty900

### QC tools

* Tarazona S, Furió-Tarí P, Turrà D, Pietro AD, Nueda MJ, Ferrer A, et al. Data quality aware analysis of differential expression in RNA-seq with NOISeq R/Bioc package. Nucleic Acids Res. 2015;43:e140.
* Risso D, Schwartz K, Sherlock G, Dudoit S. GC-content normalization for RNA-seq data. BMC Bioinformatics. 2011;12:480.
* multiqc
* qualimap
* dupradar

### Normalization

* Evans C, Hardin J, Stoebel DM. Selecting between-sample RNA-Seq normalization methods from the perspective of their assumptions. Brief Bioinform. 2018 Sep 28;19(5):776-792. doi: 10.1093/bib/bbx008.
  * normalization methods fail to produce results comparable between samples for different mRNA/cell or global expression level shift
* Zhao S, Ye Z, Stanton R. Misuse of RPKM or TPM normalization when comparing across samples and sequencing protocols. RNA. 2020 Aug;26(8):903-909. doi: 10.1261/rna.074922.120. Epub 2020 Apr 13.
  * different sample composition (% of ribo RNA, mito RNA, ...) prevent from comparing samples meaningfully, especially samples from different studies (e.g. polyA enrichment vs ribo depletion) 
  * rna compartmentalization affect tpm values (cytosolic vs nuclear)
  * strandness has substantial impact
  * RPKM and TPM may not be meaningful for samples with varying mRNA levels

### Reproducibility

* Holcomb et al. 2021. Benchmarking Single-Cell mRNA–Sequencing Technologies Uncovers Differences in Sensitivity and Reproducibility in Cell Types With Low RNA Content. J Biomol Tech. 2021; 32(4): 3fc1f5fe.dbeabb2a. https://doi.org/10.7171%2F3fc1f5fe.dbeabb2a
  * Takara Bio USA
  * New kit improves in convenience, sensitivity, gene identification, and reproducibility
* Li et al. Genes expressed at low levels raise false discovery rates in RNA samples contaminated with genomic DNA. BMC Genomics
. 2022 Aug 3;23(1):554. https://doi.org/10.1186/s12864-022-08785-1
  * TO READ
  * contamination of RNA-seq samples with genomic DNA impact gene expression
* Leech et al. Incomplete reporting of manual therapy interventions and a lack of clinician and setting diversity in clinical trials for neck pain limits replication and real-world translation. A scoping review. J Man Manip Ther
. 2022 Sep 1;1-10. doi: 10.1080/10669817.2022.2113295.
  * clinical data and poor metadata limits reproducibility
* Eran Elhaik.  Principal Component Analyses (PCA)‑based fndings in population genetic studies are highly biased and must be reevaluated. Scientifc Reports. 2022; 12:14683. https://doi.org/10.1038/s41598-022-14395-4
  * PCA calculations instability question published results
