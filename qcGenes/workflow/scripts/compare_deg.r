# Returns a union of 2 tables of differential genes (Rasflow deg files)
# Additional boolean columns identify differential genes in each input table
#
library(tidyverse)
# setwd("/media/jf/hdd4/backups/projects/qcGenes.truong.bopp.git")

args <- commandArgs(trailingOnly = TRUE)
file.x   = args[1] # First table
file.y   = args[2] # Second table
out.dir  = args[3] # Output directory
out.file = file.path(args[3], args[4]) # Output file path
gtf.path = args[5] # GTF file path

# file.x   = "output/main/RASflowResults/project1A/trans/dea/DEA/gene-level/deg_WTwoHCl_KOwoHCl.tsv"
# file.y   = "output/main/RASflowResults/project1B/trans/dea/DEA/gene-level/deg_WTwiHCl_KOwiHCl.tsv"
# out.dir  = "output/main/degComparison"
# out.file = "output/main/degComparison/IcerKO_vs_Gpr65KO.woHcl.tsv"
# gtf.path = "data/ref/GRCm38.101/Mus_musculus.GRCm38.101.gtf.gz"

ensembl.ids.table  <- as_tibble(read.delim(gtf.path, header=F,  comment.char = "#", stringsAsFactors=F)) %>%
  filter(V3=="transcript") %>%
  mutate(ensemblid=str_extract(V9, regex("(?<=gene_id )(.+?)(?=;)", dotall = TRUE))) %>%
  mutate(symbol=str_extract(V9, regex("(?<=gene_name )(.+?)(?=;)", dotall = TRUE))) %>%
  select(ensemblid, symbol) %>%
  unique()

table.x = read.delim(file.x, sep="\t", row.names=1)
table.y = read.delim(file.y, sep="\t", row.names=1)
table.x = table.x %>% rownames_to_column() %>% rename(ensemblid=rowname) %>% select(-baseMean, -stat, -pvalue, -lfcSE)
table.y = table.y %>% rownames_to_column() %>% rename(ensemblid=rowname) %>% select(-baseMean, -stat, -pvalue, -lfcSE)
table.x = tibble(table.x)
table.y = tibble(table.y)

dir.create(out.dir)

p <- table.x %>% 
  full_join(table.y, by="ensemblid") %>% 
  mutate(deg.x=ifelse(!is.na(padj.x) & padj.x<0.05, 1, 0)) %>% 
  mutate(deg.y=ifelse(!is.na(padj.y) & padj.y<0.05, 1, 0)) %>%
  mutate(deg.xy=ifelse(deg.x==1 & deg.y==1, 1, 0)) %>%
  left_join(ensembl.ids.table, by="ensemblid") %>%
  select(ensemblid, symbol, log2FoldChange.x, padj.x, log2FoldChange.y, padj.y, deg.x, deg.y, deg.xy)
p %>% 
  write_delim(out.file, delim="\t")

openxlsx::write.xlsx(p, file = paste0(out.file, ".xlsx"))

n.table.x <- sum(p$deg.x)
n.table.y <- sum(p$deg.y)
n.inter <- sum(p$deg.xy)
print(paste("X genes: ", n.table.x, "; Y genes: ", n.table.y, "; Intersect: ", n.inter))
