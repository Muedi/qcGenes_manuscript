suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(yaml, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

config.path                   <- file.path(".", "config", "config.yaml")
config.data                   <- yaml.load_file(config.path)
gs2d.path                     <- file.path(".", config.data$GS2D_PATH)
gs2d.msigdb.path              <- file.path(".", config.data$GS2D_MSIGDB_PATH)

d<- read_tsv(gs2d.path)

d2 <- d %>% 
  filter(fdr<0.05 & fold_change>1.5 & count_name_in_gene_set>5) %>%
  group_by(name) %>%
  summarize(symbols = paste(symbol, collapse = '\t'), n=n()) %>%
  mutate(url="http://cbdm-01.zdv.uni-mainz.de/~jfontain/cms/?page_id=592") %>%
  filter(n>10) %>%
  select(name, url, symbols)

write_tsv(d2, gs2d.msigdb.path, col_names=FALSE)
