library(tidyverse)
# ==============================================================================
datasets.path <- file.path("config", "metadata", "Datasets.csv")
samples.path <- file.path("config", "metadata", "Mega_SraRunTable.csv")
stats.path <- file.path("output", "main", "qualityCorGenes", "datasets.stats.tsv")
# ==============================================================================

# DATASETS METADATA
# d <- read_csv(datasets.path) %>% # , show_col_types = FALSE) %>%
d <- read_csv(datasets.path) %>% 
# select(GEO_Series, select, batches, SamplesPairing, Control, Treat, Critical, library_purification) %>% 
  filter(select==1) %>% 
  filter(GEO_Series!="GSE100925")
d %>% head()

# SAMPLES METADATA
samples <- read_csv(samples.path) %>% #, show_col_types = FALSE) %>% 
  filter(Selected==1) %>% 
  filter(!is.na(Age)) %>% 
  filter(!is.na(gender)) %>% 
  mutate(Age=round(as.numeric(Age), 0)) %>% 
  mutate(gender=str_to_lower(gender)) %>% 
  filter(Age>40) %>% 
#  filter(Age>40 & Age<60) %>% 
#  filter(Age>20 & Age<60) %>% 
  group_by(GEO_Series) %>% 
  mutate(n.male=sum(gender=="male")) %>% 
  mutate(n.female=sum(gender=="female")) %>% 
  mutate(majority.gender=if_else(n.male>n.female, "male", "female")) %>%
  filter(gender==majority.gender) %>% 
  select(-n.male, -n.female, -majority.gender) %>% 
  mutate(n.batches=length(unique(batch))) %>% 
  filter(n.batches==1) %>% 
  select(-n.batches) %>% 
  mutate(n.groups=length(unique(group))) %>% 
  filter(n.groups==2) %>% 
  mutate(n=n()) %>% 
  filter(n>30) %>% 
  arrange(n) 

# SEE NUMBER OF SELECTED SAMPLES PER GROUP
samples %>% 
  # select(GEO_Series, Run, group, subject, Age, gender) %>% 
  select(GEO_Series, group, n) %>% 
  group_by(GEO_Series, group) %>% 
  mutate(n=n()) %>% 
  distinct()

q(save="no")

# DATASETS STATISTICS
stats.table <- read_tsv(stats.path) %>% #, show_col_types = FALSE) %>% 
  select(dataset, p_group_cor, n_samples, n_batches, n_degs, n_degs_fdr_only, selected, SamplesPairing, mesh_terms) %>% 
  rename(GEO_Series=dataset) %>% 
  filter(GEO_Series %in% unique(samples$GEO_Series))
# stats.table %>% head()


# CREATE DATASET SUBSETS OF 5 SAMPLES PER GROUP
set.seed(12345678)
fake.dataset.table <- tibble()
# for (dataset in c("GSE116250")) {
for (dataset in unique(samples$GEO_Series)) {
  print(dataset)
  dataset.samples <- samples %>% filter(GEO_Series==dataset)
  groups <- unique(dataset.samples$group)
  scores.path <- file.path("output", "main", "scores", paste(dataset, "scores", "txt", sep="."))
  scores.table <- read_tsv(scores.path, col_names = c("sample", "p", "dummy")) %>% 
    select(-dummy)
  # scores.table %>% head()
  combin.list <- list()
  min.n.combin <- 10000000
  for (i in 1:length(groups)) {
    my.group <- groups[i]
    my.group.samples.ids <- (dataset.samples %>% filter(group==my.group))$Run
    combin.list[[i]] <- combn(my.group.samples.ids, 5)
    min.n.combin <- min(dim(combin.list[[i]])[2], min.n.combin)
  }
  # COMBINATION OF COLUMNS
  all.col.pairs.to.combine <- combn(min.n.combin, 2)
  print(dim(all.col.pairs.to.combine))
  # all.col.pairs.to.combine.index <- round(
  #   quantile(1:ncol(all.col.pairs.to.combine), probs = seq(0, 1, 0.05)), 
  #   0)
  all.col.pairs.to.combine.index <- 1:ncol(all.col.pairs.to.combine)
  if(ncol(all.col.pairs.to.combine)>200){
    all.col.pairs.to.combine.index <- sample(1:ncol(all.col.pairs.to.combine), 200)
  }
  all.col.pairs.to.combine <- all.col.pairs.to.combine[, as.numeric(all.col.pairs.to.combine.index)]
  print(dim(all.col.pairs.to.combine))
  
  for (col.pair.to.combine.index in 1:ncol(all.col.pairs.to.combine)) {
    fake.dataset.id <- paste(dataset, sprintf("%02d", col.pair.to.combine.index), sep="_")         # FAKE DATASET ID
    col.pair.to.combine <- all.col.pairs.to.combine[,col.pair.to.combine.index]
    group1.sample.ids <- combin.list[[1]][, col.pair.to.combine[1]]               # SAMPLE IDS FOR GROUP 1
    group2.sample.ids <- combin.list[[2]][, col.pair.to.combine[2]]               # SAMPLE IDS FOR GROUP 2
    fake.dataset.quality.table <- bind_rows(
      data.frame(dataset=dataset, fake.dataset=fake.dataset.id, sample=group1.sample.ids, group=0),
      data.frame(dataset=dataset, fake.dataset=fake.dataset.id, sample=group2.sample.ids, group=1)) %>% 
      left_join(scores.table, by="sample")
    p_group_cor <- abs(cor(fake.dataset.quality.table$p, fake.dataset.quality.table$group))
    fake.dataset.quality.table$p_group_cor <- p_group_cor
    # fake.dataset.quality.table %>% head()
    fake.dataset.table <- bind_rows(
      fake.dataset.table,
      fake.dataset.quality.table
    )
  }
}

# TABLE OF DATASET SUBSETS
fake.dataset.table %>% head()

# TABLE OF RANKS BY QUALITY
ranks.table <- fake.dataset.table %>%
  select(dataset, fake.dataset, p_group_cor) %>% 
  distinct() %>%
  group_by(dataset) %>% 
  mutate(rank=round(100*rank(p_group_cor)/n(), 0))
# sort(unique(ranks.table$rank))

ranks.table %>% 
  ggplot(aes(x=p_group_cor))+
  geom_histogram()+
  facet_wrap(vars(dataset))

# SELECTION OF A FEW DATASETS BY RANKS
fake.dataset.selection <- ranks.table %>% 
  # filter(rank %in% c(1,2,3, 20,21,22, 41,42,43, 61,62,63, 81,82,83, 98,99,100)) %>% 
  filter(rank %in% c(1,4,8, 21,24,28, 41,44,48, 61,62,63, 81,82,83, 98,99,100)) %>% 
  group_by(dataset, rank) %>% 
  slice_sample(n=1) %>% 
  ungroup()
# fake.dataset.selection %>% print(n=40)
fake.dataset.selection %>% head()

# CREATE DATASETS TABLE
final.datasets.table <- fake.dataset.selection %>%
  rename(GEO_Series=dataset) %>% 
  left_join(d, by="GEO_Series") %>% 
  select(-GEO_Series) %>% 
  rename(GEO_Series=fake.dataset)
final.datasets.table %>% head()
write_csv(final.datasets.table, file.path("config", "metadata", "Datasets.subsets.csv"))


# CREATE SAMPLES TABLE
final.samples.table <- fake.dataset.table %>% 
  rename(Run=sample) %>% 
  rename(source_dataset=dataset) %>% 
  select(-group, -p, -p_group_cor) %>% 
  left_join(samples, by="Run") %>% 
  # select(-GEO_Series) %>% 
  mutate(GEO_Series=fake.dataset) %>% 
  select(-fake.dataset)
final.samples.table %>% head()
write_csv(final.samples.table, file.path("config", "metadata", "Mega_SraRunTable.subsets.csv"))

# CREATE BASH LINKS FOR THE DATASETS FILES

commands <- fake.dataset.selection %>% 
  # select(Run, source_dataset, GEO_Series) %>% 
  mutate(cmd0=paste0("mkdir ", fake.dataset)) %>%
  mutate(cmd1=paste0("ln -s ../../../qcGenes.git/data/datasets/", dataset, " ", fake.dataset)) %>%
  mutate(cmd2=paste0("ln -s ../../../qcGenes.git/data/output/", dataset, " ", fake.dataset)) %>%
  mutate(cmd3=paste0("ln -s ../../../../qcGenes.git/output/main/scores/", dataset, ".scores.txt ", fake.dataset, ".scores.txt"))

write_csv(commands %>% select(cmd0), file.path("config", "metadata", "commands.for.data.dirs.csv"), col_names = F)
write_csv(commands %>% select(cmd1), file.path("config", "metadata", "commands.for.data.datasets.links.csv"), col_names = F)
write_csv(commands %>% select(cmd2), file.path("config", "metadata", "commands.for.data.output.links.csv"), col_names = F)
write_csv(commands %>% select(cmd3), file.path("config", "metadata", "commands.for.output.main.scores.links.csv"), col_names = F)
