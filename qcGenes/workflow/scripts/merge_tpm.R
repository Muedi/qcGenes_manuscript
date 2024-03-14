suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

args = commandArgs(trailingOnly=TRUE)

# args <-
# c(
#   "GSE123496",
#   # "output/data/ref/GRCh38.101/Homo_sapiens.GRCh38.101.gtf.gz",
#   "output/data/ref/GRCh38.101/Homo_sapiens.GRCh38.101.gene2tx.full.ver.tsv",
#   "output/main/gene_expression/GSE123496/tpm/GSE123496.samples.tpm.tsv.gz",
#   "output/data/output/GSE123496/salmon/SRR8307929_tpm.tsv",
#   "output/data/output/GSE123496/salmon/SRR8307934_tpm.tsv",
#   "output/data/output/GSE123496/salmon/SRR8307939_tpm.tsv",
#   "output/data/output/GSE123496/salmon/SRR8307944_tpm.tsv",
#   "output/data/output/GSE123496/salmon/SRR8307949_tpm.tsv",
#   "output/data/output/GSE123496/salmon/SRR8307954_tpm.tsv",
#   "output/data/output/GSE123496/salmon/SRR8307959_tpm.tsv",
#   "output/data/output/GSE123496/salmon/SRR8307964_tpm.tsv",
#   "output/data/output/GSE123496/salmon/SRR8307969_tpm.tsv",
#   "output/data/output/GSE123496/salmon/SRR8307974_tpm.tsv"
# )

data.id <- args[1]
annotation.file <- args[2]
output.path <- args[3]
tpm.file.paths <- args[4:length(args)]

f <- tpm.file.paths[1]
file.name <- stringr::str_remove(basename(f), "_tpm.tsv") 
a <- read_tsv(f, col_names = c("trans_id_ver", file.name), show_col_types = FALSE)
 
for (f in tpm.file.paths[2:length(tpm.file.paths)]) {
  file.name <- stringr::str_remove(basename(f), "_tpm.tsv") 
  b <- read_tsv(f, col_names = c("trans_id_ver", file.name), show_col_types = FALSE)
  a <- a %>% left_join(b, by="trans_id_ver")
}

biomart <- read_tsv(annotation.file, show_col_types = FALSE) %>% 
  select(trans_id_ver, trans_name, gene_id_ver, gene_name)

a <- a %>% 
  left_join(biomart, by="trans_id_ver") %>% 
  select(trans_id_ver, trans_name, gene_id_ver, gene_name, everything()) 

a %>% write_tsv(output.path)
