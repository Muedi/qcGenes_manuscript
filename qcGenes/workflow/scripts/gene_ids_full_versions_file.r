suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
# suppressMessages(library(yaml, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

# config.path <- file.path(".", "config", "config.yaml")
# config.data <- yaml.load_file(config.path)
# input.file  <- config.data$ANNOT_PATH
# output.file <- config.data$ENSIDFULLVER_PATH

args = commandArgs(trailingOnly=TRUE)
input.file  <- args[1]
output.file <- args[2]

biomart <- as_tibble(read.delim(input.file, header=F,  comment.char = "#", stringsAsFactors=F)) %>% 
  filter(V3=="transcript") %>% 
  mutate(gene_id=str_extract(V9, regex("(?<=gene_id )(.+?)(?=;)", dotall = TRUE))) %>% 
  mutate(gene_ver=str_extract(V9, regex("(?<=gene_version )(.+?)(?=;)", dotall = TRUE))) %>% 
  mutate(gene_id_ver=paste0(gene_id, ".", gene_ver)) %>% 
  mutate(gene_name=str_extract(V9, regex("(?<=gene_name )(.+?)(?=;)", dotall = TRUE))) %>% 
  mutate(trans_id=str_extract(V9, regex("(?<=transcript_id )(.+?)(?=;)", dotall = TRUE))) %>% 
  mutate(trans_ver=str_extract(V9, regex("(?<=transcript_version )(.+?)(?=;)", dotall = TRUE))) %>% 
  mutate(trans_id_ver=paste0(trans_id, ".", trans_ver)) %>% 
  mutate(trans_name=str_extract(V9, regex("(?<=transcript_name )(.+?)(?=;)", dotall = TRUE))) %>% 
  select(trans_id_ver, trans_id, trans_name, gene_id_ver, gene_id, gene_name)

write_tsv(biomart, output.file)
