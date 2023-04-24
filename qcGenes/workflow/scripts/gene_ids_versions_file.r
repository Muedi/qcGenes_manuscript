#!/usr/bin/env Rscript
# ==============================================================================
# Create Gene-Transcript IDs file 
# ==============================================================================
# Rscript workflow/scripts/gene_ids_versions_file.r

suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(yaml, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

config.path <- file.path(".", "config", "config.yaml")
config.data <- yaml.load_file(config.path)
input.file  <- config.data$ANNOT_PATH
output.file <- config.data$ENSIDVER_PATH

ensembl_ids <- as_tibble(read.delim(input.file, header=F,  comment.char = "#", stringsAsFactors=F)) %>% 
  filter(V3=="transcript") %>% 
  mutate(gene_id=str_extract(V9, regex("(?<=gene_id )(.+?)(?=;)", dotall = TRUE))) %>% 
  mutate(gene_ver=str_extract(V9, regex("(?<=gene_version )(.+?)(?=;)", dotall = TRUE))) %>%
  mutate(gene_id=paste0(gene_id, ".", gene_ver)) %>%
  mutate(transcript_id=str_extract(V9, regex("(?<=transcript_id )(.+?)(?=;)", dotall = TRUE))) %>% 
  mutate(transcript_ver=str_extract(V9, regex("(?<=transcript_version )(.+?)(?=;)", dotall = TRUE))) %>% 
  mutate(transcript_id=paste0(transcript_id, ".", transcript_ver)) %>%
  select(transcript_id, gene_id)

write_tsv(ensembl_ids, output.file)
