# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("tximport")

# ______________________________________________________________________________
# INIT
# ______________________________________________________________________________
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(DESeq2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(tximport, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(rstatix, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(yaml, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

args = commandArgs(trailingOnly=TRUE)

config.path                 <- file.path(".", "config", "config.yaml")
config                      <- yaml.load_file(config.path)
FILTER_OUT_QUALITY_OUTLIERS <- config$FILTER_OUT_QUALITY_OUTLIERS
QUALITY_OUTLIERS_COEF       <- config$QUALITY_OUTLIERS_COEF
P_LOW_AS_CONFOUNDER         <- config$P_LOW_AS_CONFOUNDER
samples.path                <-  config$SAMPLES_FILE


# TEST VALUES
# ______________________________________________________________________________
# test.dataid <- "GSE85567"
# args <- c(
#   test.dataid,
#   "config/metadata/Datasets.csv",
#   "output/data/ref/GRCh38.101/Homo_sapiens.GRCh38.101.gene2tx.ver.tsv", # ENSIDVER_PATH
#   file.path("output", "main", "scores", paste(test.dataid, "scores", "txt", sep=".")),
#   file.path("output", "main", "gene_expression", test.dataid, paste(test.dataid, "samples", "rlg", "tsv", "gz", sep=".")),
#   file.path("output", "main", "gene_expression", test.dataid, paste(test.dataid, "diff", "genes", "tsv", "gz", sep=".")),
#   file.path(
#     file.path("output/data/output", test.dataid, "salmon"),
#     grep(pattern="quant.sf",
#          x=dir(file.path("output/data/output", test.dataid, "salmon"), recursive=T),
#          value=TRUE)
#   )
# )


# VARIABLES
# ______________________________________________________________________________

dataid           <- args[1]
datasets.path    <- args[2]
tx2gene.path     <- args[3] # ENSIDVER_PATH
scores.path      <- args[4]
rlog.output.path <- args[5]
diff.output.path <- args[6]
files            <- args[7:length(args)]
names(files)     <- basename(dirname(files))
files <- sort(files)
# length(files)


# REF AND DATA FILES
# ______________________________________________________________________________

tx2gene <- read_tsv(tx2gene.path, show_col_types = FALSE)

datasets.metadata <- read_csv(datasets.path, show_col_types = FALSE) %>% 
  mutate(Control = str_replace_all(Control, "-", "_")) %>% 
  mutate(Treat = str_replace_all(Treat, "-", "_")) %>% 
  filter(GEO_Series==dataid) %>% 
  select(GEO_Series, batches, LibraryLayout, SamplesPairing, Control, Treat, `mesh_terms (pipe-separated list)`, `5-YEAR-IF2020`)
PAIRING <- datasets.metadata$SamplesPairing

samples <- read_csv(samples.path, show_col_types = FALSE) %>% 
  mutate(group = str_replace_all(group, "-", "_")) %>% 
  filter(GEO_Series==dataid) %>% 
  filter(Selected==1) %>% 
  filter(group %in% c(datasets.metadata$Control, datasets.metadata$Treat)) %>% 
  mutate(group = factor(group, levels=c(datasets.metadata$Control, datasets.metadata$Treat))) %>% 
  # mutate(group = relevel(group, ref=datasets.metadata$Control)) %>%
  mutate(subject = as.factor(subject)) %>% 
  mutate(gender = as.factor(gender)) %>% 
  select(subject, Run, group, LibraryLayout, Age, gender, batch, sequenced_molecule, Important_notes) %>% 
  arrange(Run)

# READ SCORES FILE
scores.metadata <- read_tsv(scores.path, col_names=c("Run", "P_low", "note"), show_col_types = F) %>% 
  select(-note) %>% 
  left_join(samples %>% select(Run, group), by="Run") %>% 
  group_by(group) %>% 
  mutate(outlier=is_outlier(P_low, coef = QUALITY_OUTLIERS_COEF)) %>% 
  ungroup()

samples <- samples %>% 
  left_join(scores.metadata %>% select(Run, P_low, outlier), by="Run")

# FILTER OUT QUALITY OUTLIERS
no.outlier.samples <- samples$Run
if(FILTER_OUT_QUALITY_OUTLIERS){
  if( datasets.metadata$SamplesPairing ){
    no.outlier.samples <- samples %>% 
      group_by(subject) %>%
      mutate(outlier = max(as.numeric(outlier))) %>% 
      filter(!outlier) %>% 
      pull(Run)
  } else{
    no.outlier.samples <- scores.metadata %>% 
      filter(!outlier) %>% 
      pull(Run)
  }
}
remaining.files <- base::intersect(names(files), no.outlier.samples) 
files <- files[remaining.files]
samples <- samples %>% 
  filter(Run %in% remaining.files) %>% 
  mutate(subject=as.factor(as.character(subject)))


# Create DDS Object
# ______________________________________________________________________________
txi <- tximport(files, type="salmon", tx2gene=tx2gene)

# length(files)

# ddsTxi <- DESeqDataSetFromTximport(
#   txi,
#   colData = samples ,
#   design = ~ group)

if(PAIRING==0 & P_LOW_AS_CONFOUNDER){
  # design(ddsTxi) <- formula(~ P_low + group)
  print("formula: ~ P_low + group")
  ddsTxi <- DESeqDataSetFromTximport(
    txi,
    colData = samples ,
    design = ~ P_low + group)
} else if(PAIRING==1 & !P_LOW_AS_CONFOUNDER){
  # # design(ddsTxi) <- formula(~ group + subject) # ~batch + condition
  # design(ddsTxi) <- formula(~ subject + group)
  print("formula: ~ subject + group")
  ddsTxi <- DESeqDataSetFromTximport(
    txi,
    colData = samples ,
    design = ~ subject + group)
} else if(PAIRING==1 & P_LOW_AS_CONFOUNDER){
  # design(ddsTxi) <- formula(~ subject + P_low + group)
  print("formula: ~ subject + P_low + group")
  ddsTxi <- DESeqDataSetFromTximport(
    txi,
    colData = samples ,
    design = ~ subject + P_low + group)
} else {
  # design(ddsTxi) <- formula(~ group)
  print("formula: ~ group")
  ddsTxi <- DESeqDataSetFromTximport(
    txi,
    colData = samples ,
    design = ~ group)
}


# Filtering
# ______________________________________________________________________________
# smallestGroupSize <- 3
# keep <- rowSums(counts(ddsTxi) >= 10) >= smallestGroupSize
# ddsTxi <- ddsTxi[keep,]


# Differential Analysis
# ______________________________________________________________________________
# ddsTxi <- DESeq(ddsTxi)
# ddsTxi <- DESeq(ddsTxi, sfType = "iterate")
ddsTxi <- DESeq(ddsTxi, test = "Wald", fitType = "parametric", sfType = "poscounts")
rlog.table <- as.data.frame(assay(rlog(ddsTxi, blind=TRUE))) %>%
  rownames_to_column(var = "geneid")

ddsTxi.res <- results(ddsTxi)
ddsTxi.res <- ddsTxi.res[order(ddsTxi.res$pvalue),]
ddsTxi.res.export.table <- as.data.frame(ddsTxi.res) %>% rownames_to_column("geneid")


# Write output
# ______________________________________________________________________________
write_tsv(rlog.table, rlog.output.path)
write_tsv(ddsTxi.res.export.table, diff.output.path)
