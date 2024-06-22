suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(org.Hs.eg.db, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Please provide input and output file paths", call.=FALSE)
}
input.path <- args[1]
output.path <- args[2]

disease.genes <- readr::read_tsv(input.path, show_col_types = F) %>%
  dplyr::rename(entrezid=gene_id) %>%
  dplyr::filter(fdr<=0.05)

disease.genes$geneid <- mapIds(
  org.Hs.eg.db, 
  keys = as.character(disease.genes$entrezid), 
  keytype="ENTREZID", 
  column = "ENSEMBL",
  multiVals = "first")

readr::write_tsv(disease.genes, output.path)
