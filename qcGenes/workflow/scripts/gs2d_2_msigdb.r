suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(yaml, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Please provide input and output file paths", call.=FALSE)
}
input.path <- args[1]
output.path <- args[2]

d<- read_tsv(input.path)

d2 <- d %>% 
  filter(fdr<0.05 & fold_change>1.5 & count_name_in_gene_set>5) %>%
  group_by(name) %>%
  summarize(symbols = paste(symbol, collapse = '\t'), n=n()) %>%
  mutate(url="http://cbdm-01.zdv.uni-mainz.de/~jfontain/cms/?page_id=592") %>%
  filter(n>10) %>%
  select(name, url, symbols)

write_tsv(d2, output.path, col_names=FALSE)
