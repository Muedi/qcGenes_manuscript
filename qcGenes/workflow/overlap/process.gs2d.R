# setwd("~/projects/qcGenes")
# install.packages("BiocManager")
# BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
library(magrittr)

ref.data.dir <- file.path("data", "ref", "gs2d")
dir.create(ref.data.dir, showWarnings = F, recursive = T)

disease.genes.url <- file.path("http://cbdm-01.zdv.uni-mainz.de/~jfontain/downloads/20180416_gene2mesh.tsv")
output.path <- file.path(ref.data.dir, "20180416_gene2mesh.ens.tsv.gz")

disease.genes <- readr::read_tsv(disease.genes.url, show_col_types = F) %>%
  dplyr::rename(entrezid=gene_id) %>%
  dplyr::filter(fdr<=0.05)

disease.genes$geneid <- mapIds(org.Hs.eg.db, 
                  keys = as.character(disease.genes$entrezid), 
                  keytype="ENTREZID", 
       column = "ENSEMBL",
       multiVals = "first")

readr::write_tsv(disease.genes, output.path)
