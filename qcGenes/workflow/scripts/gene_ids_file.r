#!/usr/bin/env Rscript
# ==============================================================================
# Create Gene-Transcript IDs file 
# ==============================================================================
# Rscript workflow/scripts/gene_ids_file.r

suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(yaml, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

config.path <- file.path(".", "config", "config.yaml")
config.data <- yaml.load_file(config.path)
input.file  <- config.data$ANNOT_PATH
output.file <- config.data$ENSID_PATH

ensembl_ids <- as_tibble(read.delim(input.file, header=F,  comment.char = "#", stringsAsFactors=F)) %>% 
  filter(V3=="transcript") %>% 
  mutate(gene_id=str_extract(V9, regex("(?<=gene_id )(.+?)(?=;)", dotall = TRUE))) %>% 
  mutate(gene_name=str_extract(V9, regex("(?<=gene_name )(.+?)(?=;)", dotall = TRUE))) %>% 
  mutate(trans_id=str_extract(V9, regex("(?<=transcript_id )(.+?)(?=;)", dotall = TRUE))) %>% 
  mutate(trans_name=str_extract(V9, regex("(?<=transcript_name )(.+?)(?=;)", dotall = TRUE))) %>% 
  select(gene_id, trans_id, gene_name, trans_name)

write_tsv(ensembl_ids, output.file)
