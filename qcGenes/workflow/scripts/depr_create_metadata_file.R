#!/usr/bin/env Rscript
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))


# ARGUMENTS
# ------------------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
if (length(args)<1 | length(args)>3) {
  stop("Please provide GEO_Series ID, analysis type (subdir) and output path", call.=FALSE)
}
dataid <- args[1]
subdir <- args[2]
outpath <- args[3]

#dataid <- "GSE111523"
#subdir <- "main"
#outpath <- "output/main/pipelines/GSE111523/configs/metadata.tsv"

sra.path <- file.path(".", "config", "metadata", "Mega_SraRunTable.csv")
sra.table <- read_csv(sra.path) %>%
  rename(sample=Run) %>%
  filter(GEO_Series==dataid)

output.table <- sra.table %>% filter(Selected==1) %>% select("sample", "group", "subject", "batch")

if(subdir=="batched"){
    output.table <- sra.table %>% filter(BatchAnalysisSelection==1) %>% select("sample", "group", "subject", "batch")
}

write_tsv(output.table, file=outpath)




