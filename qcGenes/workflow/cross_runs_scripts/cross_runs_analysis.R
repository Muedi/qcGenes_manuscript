# ______________________________________________________________________________
# CROSS-RUN ANALYSIS
# ______________________________________________________________________________
library(tidyverse)

# run.0.path <- file.path("output", "main_C2") # Reference Run
# run.1.path <- file.path("output", "main_C3")
# run.2.path <- file.path("output", "main_C4")
# out.dir <- file.path("output", "cross_run_C")

run.0.path <- file.path("output", "main_P1") # Reference Run
run.1.path <- file.path("output", "main_P1_plow_conf")
run.2.path <- file.path("output", "main_P3")
run.3.path <- file.path("output", "main_P4") #_P4_default_DESEQ")
out.dir <- file.path("output", "cross_run_P")

# run.0.path <- file.path("output", "main_S2") # Reference Run
# run.1.path <- file.path("output", "main_S3")
# run.2.path <- file.path("output", "main_S4")
# out.dir <- file.path("output", "cross_run_S")


# ______________________________________________________________________________
# GET N OUTLIERS ====
# ______________________________________________________________________________
n_outliers <- read_tsv(file.path(run.0.path, "metrics", "n.outliers.tsv"), show_col_types = F)


# ______________________________________________________________________________
# GET MARKER GENES RESULTS FROM REFERENCE RUN (MARKER LIST, PATHWAYS)
# ______________________________________________________________________________

# ______________________________________________________________________________
# COMPARE PCAs (Dunn index) ====
# ______________________________________________________________________________
run.0.stats.path <- file.path(run.0.path, "qualityCorGenes", "datasets.stats.tsv")
run.1.stats.path <- file.path(run.1.path, "qualityCorGenes", "datasets.stats.tsv")
run.2.stats.path <- file.path(run.2.path, "qualityCorGenes", "datasets.stats.tsv")
run.3.stats.path <- file.path(run.2.path, "qualityCorGenes", "datasets.stats.tsv")

run.0.stats <- read_tsv(run.0.stats.path, show_col_types = F) %>% 
  mutate(confo_run="0") %>% 
  select(confo_run, dataset, p_group_cor, p_group_test, dunn_all, n_samples, n_degs, 
         n_degs_fdr_only, SamplesPairing, no.bias.group, low.bias.group, high.bias.group)
run.1.stats <- read_tsv(run.1.stats.path, show_col_types = F) %>% 
  mutate(confo_run="1") %>% 
  select(confo_run, dataset, p_group_cor, p_group_test, dunn_all, n_samples, n_degs, 
         n_degs_fdr_only, SamplesPairing, no.bias.group, low.bias.group, high.bias.group)
run.2.stats <- read_tsv(run.2.stats.path, show_col_types = F) %>% 
  mutate(confo_run="2") %>% 
  select(confo_run, dataset, p_group_cor, p_group_test, dunn_all, n_samples, n_degs, 
         n_degs_fdr_only, SamplesPairing, no.bias.group, low.bias.group, high.bias.group)
run.3.stats <- read_tsv(run.2.stats.path, show_col_types = F) %>% 
  mutate(confo_run="2") %>% 
  select(confo_run, dataset, p_group_cor, p_group_test, dunn_all, n_samples, n_degs, 
         n_degs_fdr_only, SamplesPairing, no.bias.group, low.bias.group, high.bias.group)

run.stats <- bind_rows(run.0.stats, run.1.stats, run.2.stats,run.3.stats)

dunn_table <- run.0.stats %>%
  select(dataset, n_samples, p_group_cor, p_group_test, dunn_all_0=dunn_all) %>% 
  left_join(run.1.stats %>% select(dataset, dunn_all_1=dunn_all), by="dataset") %>% 
  left_join(run.2.stats %>% select(dataset, dunn_all_2=dunn_all), by="dataset") %>% 
  left_join(run.3.stats %>% select(dataset, dunn_all_3=dunn_all), by="dataset") %>% 
  left_join(n_outliers, by="dataset")


dunn_plot <- dunn_table %>%
  filter(outliers>0) %>% 
  mutate(big_data=if_else(n_samples<=10, "small data", "medium and big data")) %>%
  mutate(bias=if_else(p_group_cor<0.3 & p_group_test>0.05, "low bias", "higher bias")) %>%
  ggplot(aes(x=dunn_all_0, y=dunn_all_2, size=n_samples, col=outliers)) + 
  geom_point() + 
  geom_abline(slope = 1) + 
  ylab("Dunn Index (no outliers)") +
  xlab("Dunn Index (all samples)") +
  ggtitle("Dunn Index") + 
  # facet_wrap(vars(big_data, bias))
  facet_wrap(vars(bias))
# dunn_plot
ggsave(file.path(out.dir, "dunn.png"), dunn_plot, height = 3, width = 7)


# ______________________________________________________________________________
# COMPARE PATHWAYS (%) ====
# ______________________________________________________________________________
run.0.pathways.path <- file.path(run.0.path, "qualityCorGenes", "pathways.comparison.tsv")
run.1.pathways.path <- file.path(run.1.path, "qualityCorGenes", "pathways.comparison.tsv")
run.2.pathways.path <- file.path(run.2.path, "qualityCorGenes", "pathways.comparison.tsv")
run.3.pathways.path <- file.path(run.2.path, "qualityCorGenes", "pathways.comparison.tsv")

run.0.pathways <- read_tsv(run.0.pathways.path, show_col_types = F) %>% 
  mutate(confo_run="0") %>% 
  select(confo_run, everything()) %>% 
  filter(n.pathways.dataset>10) %>% 
  filter(n.pathways.markers>0) %>% 
  select(type, lfc, dataset, perc.overlap.A=perc.overlap)
run.1.pathways <- read_tsv(run.1.pathways.path, show_col_types = F) %>% 
  mutate(confo_run="1") %>% 
  select(confo_run, everything()) %>% 
  filter(n.pathways.dataset>10) %>% 
  filter(n.pathways.markers>0) %>% 
  select(type, lfc, dataset, perc.overlap.B=perc.overlap)
run.2.pathways <- read_tsv(run.2.pathways.path, show_col_types = F) %>% 
  mutate(confo_run="2") %>% 
  select(confo_run, everything()) %>% 
  filter(n.pathways.dataset>10) %>% 
  filter(n.pathways.markers>0) %>% 
  select(type, lfc, dataset, perc.overlap.C=perc.overlap)
run.3.pathways <- read_tsv(run.2.pathways.path, show_col_types = F) %>% 
  mutate(confo_run="2") %>% 
  select(confo_run, everything()) %>% 
  filter(n.pathways.dataset>10) %>% 
  filter(n.pathways.markers>0) %>% 
  select(type, lfc, dataset, perc.overlap.D=perc.overlap)

comparison.table <- 
  bind_rows(
    run.0.pathways %>% left_join(run.1.pathways, by = join_by(type, lfc, dataset)) %>% mutate(comparison="DESeq 2"),
    run.0.pathways %>% left_join(run.2.pathways, by = join_by(type, lfc, dataset)) %>% mutate(comparison="Outliers out"),
    run.0.pathways %>% left_join(run.3.pathways, by = join_by(type, lfc, dataset)) %>% mutate(comparison="Outliers out & DESeq 2")
  ) %>% 
  filter(type %in% c("cell", "cure", "regu")) %>% 
  mutate(type=if_else(type=="cell", "Cell Type Signature Gene Sets", type)) %>% 
  mutate(type=if_else(type=="cure", "Curated Gene Sets", type)) %>% 
  mutate(type=if_else(type=="regu", "Regulatory Gene Sets", type))

plt <- comparison.table %>% 
  ggplot(aes(x=perc.overlap.A, y=perc.overlap.B, col=comparison)) +
  geom_point() + 
  geom_abline(slope=1) +
  facet_grid(rows=vars(type), cols=vars(lfc)) +
  xlab("Percentage of overlaping pathways before correction") + 
  ylab("Percentage of overlaping pathways after correction") + 
  ggtitle("DEG and Low-Quality Marker Enriched pathways")

ggsave(file.path(out.dir, "scatter-perc.png"), plt, height = 3, width = 7)


comparison.summary <- comparison.table %>%
  filter(!is.na(perc.overlap.A) & !is.na(perc.overlap.B)) %>%
  filter(perc.overlap.A>0) %>% 
  mutate(ratio = (perc.overlap.B-perc.overlap.A) / (perc.overlap.A + perc.overlap.B)  ) %>% # >0 if new higher than old
  mutate(worse = ratio>0.15) %>% 
  mutate(better = ratio<(-0.15)) %>% 
  mutate(irrelevant = abs(ratio)<0.15) %>% 
  group_by(comparison, type, lfc) %>% 
  summarise(
    n_datasets = n(),
    worse = sum(worse),
    no_change = sum(irrelevant),
    better = sum(better),
    # perc = round(better / n, 2),
    .groups = "drop"
  ) %>% 
  rename(
    method=comparison,
    gene_sets=type,
    deg_sign=lfc
  )

  
comparison.summary %>%
  write_tsv(file = file.path(out.dir, "pathway.comparison.summary.deseq.tsv"))
  
comparison.summary <- comparison.table %>%
  filter(!is.na(perc.overlap.A) & !is.na(perc.overlap.C)) %>%
  filter(perc.overlap.A>0) %>% 
  mutate(ratio = (perc.overlap.C-perc.overlap.A) / (perc.overlap.A + perc.overlap.C)  ) %>% # >0 if new higher than old
  mutate(worse = ratio>0.15) %>% 
  mutate(better = ratio<(-0.15)) %>% 
  mutate(irrelevant = abs(ratio)<0.15) %>% 
  group_by(comparison, type, lfc) %>% 
  summarise(
    n_datasets = n(),
    worse = sum(worse),
    no_change = sum(irrelevant),
    better = sum(better),
    # perc = round(better / n, 2),
    .groups = "drop"
  ) %>% 
  rename(
    method=comparison,
    gene_sets=type,
    deg_sign=lfc
  )

  
comparison.summary %>%
  write_tsv(file = file.path(out.dir, "pathway.comparison.summary.outlier.out.tsv"))
  
comparison.summary <- comparison.table %>%
  filter(!is.na(perc.overlap.A) & !is.na(perc.overlap.D)) %>%
  filter(perc.overlap.A>0) %>% 
  mutate(ratio = (perc.overlap.D-perc.overlap.A) / (perc.overlap.A + perc.overlap.D)  ) %>% # >0 if new higher than old
  mutate(worse = ratio>0.15) %>% 
  mutate(better = ratio<(-0.15)) %>% 
  mutate(irrelevant = abs(ratio)<0.15) %>% 
  group_by(comparison, type, lfc) %>% 
  summarise(
    n_datasets = n(),
    worse = sum(worse),
    no_change = sum(irrelevant),
    better = sum(better),
    # perc = round(better / n, 2),
    .groups = "drop"
  ) %>% 
  rename(
    method=comparison,
    gene_sets=type,
    deg_sign=lfc
  )

  
comparison.summary %>%
  write_tsv(file = file.path(out.dir, "pathway.comparison.summary.deseq.outlier.out.tsv"))
  

  

# # ______________________________________________________________________________
# # COMPARE PATHWAYS (N) ====
# # ______________________________________________________________________________
# run.0.pathways <- read_tsv(run.0.pathways.path, show_col_types = F) %>% 
#   mutate(confo_run="0") %>% 
#   select(confo_run, everything()) %>% 
#   filter(n.pathways.dataset>10) %>% 
#   filter(n.pathways.markers>0) %>% 
#   select(type, lfc, dataset, n.pathways.dataset.A=n.pathways.dataset)
# run.1.pathways <- read_tsv(run.1.pathways.path, show_col_types = F) %>% 
#   mutate(confo_run="1") %>% 
#   select(confo_run, everything()) %>% 
#   filter(n.pathways.dataset>10) %>% 
#   filter(n.pathways.markers>0) %>% 
#   select(type, lfc, dataset, n.pathways.dataset.B=n.pathways.dataset)
# run.2.pathways <- read_tsv(run.2.pathways.path, show_col_types = F) %>% 
#   mutate(confo_run="2") %>% 
#   select(confo_run, everything()) %>% 
#   filter(n.pathways.dataset>10) %>% 
#   filter(n.pathways.markers>0) %>% 
#   select(type, lfc, dataset, n.pathways.dataset.C=n.pathways.dataset)

# comparison.table <- 
#   bind_rows(
#     run.0.pathways %>% left_join(run.1.pathways, by = join_by(type, lfc, dataset)) %>% mutate(comparison="Outliers out"),
#     run.0.pathways %>% left_join(run.2.pathways, by = join_by(type, lfc, dataset)) %>% mutate(comparison="Outliers out & DESeq 2")
    
#   ) %>% 
#   filter(type %in% c("cell", "cure", "regu")) %>% 
#   mutate(type=if_else(type=="cell", "Cell Type Pathways", type)) %>% 
#   mutate(type=if_else(type=="cure", "Curated Pathways", type)) %>% 
#   mutate(type=if_else(type=="regu", "Regulation Pathways", type))

# comparison.table %>% 
#   ggplot(aes(x=n.pathways.dataset.A, y=n.pathways.dataset.B, col=comparison)) +
#   geom_point() + 
#   geom_abline(slope=1) +
#   facet_grid(rows=vars(type), cols=vars(lfc)) +
#   xlab("Number of pathways before correction") + 
#   ylab("Number of pathways after correction")
  

