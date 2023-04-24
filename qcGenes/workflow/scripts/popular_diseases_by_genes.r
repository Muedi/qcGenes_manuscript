#!/usr/bin/env Rscript
# ==============================================================================
# Generate a table of popular disease. 
# The more the disease genes (in gs2d), the more popular a disease.
# ==============================================================================

# LIBRARIES AND SYSTEM OPTIONS
# ------------------------------------------------------------------------------
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(yaml, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
options(readr.num_columns = 0) # to mute readr functions

# GLOBAL VARIABLES
# ------------------------------------------------------------------------------
config.path                   <- file.path(".", "config", "config.yaml")
config.data                   <- yaml.load_file(config.path)
gs2d.path                     <- file.path(".", config.data$GS2D_PATH)
gs2d.table                    <- read_tsv(gs2d.path) # "GeneSet 2 Disease" data
out.dir.path                  <- file.path(".", "output", "statistics")
dir.create(out.dir.path, showWarnings=FALSE, recursive=TRUE)

# CALCULATE AND SAVE STATISTICS
# ------------------------------------------------------------------------------
gs2d_top <- gs2d.table %>% 
  group_by(name) %>% 
  summarise(genes=n()) %>%
  arrange(desc(genes)) %>% 
  rename(disease=name) %>% 
  select(disease, genes)
out.file.path <- file.path(out.dir.path, "popular_diseases_by_genes.tsv")
write.table(gs2d_top, out.file.path, quote=FALSE, sep="\t", row.names = FALSE)
