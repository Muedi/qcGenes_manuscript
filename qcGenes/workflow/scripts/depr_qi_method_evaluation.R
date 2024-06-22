# install.packages("np")
# install.packages("coin")
suppressMessages(library(yaml, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(rstatix, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(broom, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(np, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))


# ______________________________________________________________________________
# Init ====
# ______________________________________________________________________________

args = commandArgs(trailingOnly=TRUE)
# args <- c(
#   "output/main/scores",
#   "output/main/metrics"
# )
# args <- c(
#   "output/Folders_3M/main_P1/scores",
#   "output/Folders_3M/main_P1/metrics"
# )
# args <- c(
#   "../qcGenes/big_runs_2024/main_P2_default_DESEQ/scores",
#   "../qcGenes/big_runs_2024/main_P2_default_DESEQ/metrics_new"
# )
# args <- c(
#   "../qcGenes/big_runs_2024/main_S1_default_DESEQ/scores",
#   "../qcGenes/big_runs_2024/main_S1_default_DESEQ/metrics_new"
# )

config.path <- file.path(".", "config", "config.yaml")
config <- yaml.load_file(config.path)

metadata.datasets.path <- config$DATASETS_FILE
metadata.files.path <- config$SAMPLES_FILE

qc.dir <- args[1]
out.dir <- args[2]
dir.create(out.dir, showWarnings = F, recursive = T)

N_PERMUTATIONS <- config$N_PERMUTATIONS
OUTLIER_COEF <- config$QUALITY_OUTLIERS_COEF

set.seed(12345)


# ______________________________________________________________________________
# FUNCTIONS
# ______________________________________________________________________________

gmd <- function(x){
  gmd <- 0
  for(i in 1:length(x)){
    for(j in 1:length(x)){
      gmd <- gmd + abs(x[i]-x[j])
    }
  }
  return(gmd/length(x)^2)
}
# gmd(1:8)

ctdiff <- function(x, y){  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7514071/
  (median(x)-median(y)) / 
    sqrt( (gmd(x)^2 + gmd(y)^2) /  2 )
}

ctdiff_value_class <- function(x, binary_class){  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7514071/
  # x <- sample(20)
  # binary_class <- c(rep(1, 10), rep(2,10))
  classes <- levels(factor(binary_class))
  difference <- abs(ctdiff(x[binary_class==classes[1]], x[binary_class==classes[2]]))
  return(difference)
}

# ctdiff(10:20, 12:22)

ctdiff.permut.test <- function(x, binary_class, permutations=1000){
  # x <- sample(20)
  # binary_class <- c(rep(1, 10), rep(2,10))
  classes <- levels(factor(binary_class))
  difference <- abs(ctdiff(x[binary_class==classes[1]], x[binary_class==classes[2]]))
  random.diffs <- c()
  for (i in 1:permutations) {
    binary_class <- sample(binary_class)
    random.diffs <- c(
      random.diffs, 
      abs(ctdiff(x[binary_class==classes[1]], x[binary_class==classes[2]]))
    )
  }  
  n.better.values <-length(random.diffs[random.diffs>difference])
  pvalue <- n.better.values / permutations
  return(pvalue)
}
# ctdiff.permut.test(scores.table$p, scores.table$group_numeric)
# ctdiff.permut.test(scores.table$p, scores.table$group_numeric)

permut_value_class <- function(x, binary_class, myfunc=mean){  
  # x <- sample(20)
  # binary_class <- c(rep(1, 10), rep(2,10))
  classes <- levels(factor(binary_class))
  difference <- abs(myfunc(x[binary_class==classes[1]])-myfunc(x[binary_class==classes[2]]))
  return(difference)
}

permut.test.myfunc <- function(x, binary_class, permutations=1000, myfunc=mean){
  # x <- sample(20)
  # binary_class <- c(rep(1, 10), rep(2,10))
  classes <- levels(factor(binary_class))
  difference <- abs(myfunc(x[binary_class==classes[1]])-myfunc(x[binary_class==classes[2]]))
  random.diffs <- c()
  for (i in 1:permutations) {
    binary_class <- sample(binary_class)
    random.diffs <- c(
      random.diffs,
      abs(myfunc(x[binary_class==classes[1]])-myfunc(x[binary_class==classes[2]]))
    )
  }
  n.better.values <-length(random.diffs[random.diffs>difference])
  pvalue <- n.better.values / permutations
  return(pvalue)
}
# # permut.test.myfunc(scores.table$p, scores.table$group_numeric, myfunc=mean)
# # permut.test.myfunc(scores.table$p, scores.table$group_numeric, myfunc=median)

permut.test.corfunc <- function(x, binary_class, permutations=1000, mymethod="pearson"){
  # x <- sample(20)
  # binary_class <- c(rep(1, 10), rep(2,10))
  classes <- levels(factor(binary_class))
  cor.coef <- abs(cor(x, binary_class, method=mymethod))
  random.cors <- c()
  for (i in 1:permutations) {
    binary_class <- sample(binary_class)
    random.cors <- c(
      random.cors,
      abs(cor(x, binary_class, method=mymethod))
    )
  }
  n.better.values <-length(random.cors[random.cors>cor.coef])
  pvalue <- n.better.values / permutations
  return(pvalue)
}
# permut.test.corfunc(scores.table$p, scores.table$group_numeric, myfunc=mean)
# permut.test.corfunc(scores.table$p, scores.table$group_numeric, myfunc=median)

# ______________________________________________________________________________
# ______________________________________________________________________________



# ______________________________________________________________________________
# Data ====
# ______________________________________________________________________________

# Datasets table
metadata.datasets <- read_csv(metadata.datasets.path, show_col_types = F) %>%
  # filter(select_for_article==1) %>% 
  filter(select==1) %>% 
  select(GEO_Series, SamplesPairing, Control, Treat)

# Samples table merged with Datasets table
metadata.files <- read_csv(metadata.files.path, show_col_types = F) %>% 
  filter(Selected==1) %>% 
  filter(GEO_Series %in% unique(metadata.datasets$GEO_Series)) %>% 
  select(Selected, group, Run, GEO_Series, Age, gender, batch) %>% 
  left_join(metadata.datasets, by = join_by(GEO_Series))

# Load all scores files
qc.table <- bind_rows(
  read_tsv(
    file.path(qc.dir, dir(file.path(qc.dir))[dir(file.path(qc.dir))!="simulations"]), 
    col_names=c("Run", "prob", "note"),
    show_col_types = F
  ) %>% select(-note)
)

# Merge QC and metadata tables
data <- qc.table %>% 
  left_join(metadata.files, by = join_by(Run)) %>% 
  rename(sample=Run, selected=Selected, dataset=GEO_Series, age=Age, 
         is_paired=SamplesPairing, control=Control, treat=Treat) %>% 
  filter(!is.na(dataset))

# Preparing analysis data (incl. outlier annotation)
analysis.data <- data %>% 
  mutate(is_control = if_else(group==control, 1, 0)) %>% 
  select(-selected, -group, -control, -treat) %>% 
  group_by(dataset, is_control) %>%
  mutate(outlier=is_outlier(prob, coef = OUTLIER_COEF)) %>% 
  mutate(extreme=is_extreme(prob)) %>% 
  ungroup() %>% 
  group_by(dataset) %>% 
  mutate(
    n_samples=n(),
    n_samples_nooutliers=n_samples-sum(outlier)
    ) %>% 
  ungroup() %>% 
  select(dataset, sample, is_control, prob, everything())

analysis.data %>% 
  group_by(dataset) %>% 
  summarise(outliers=sum(outlier)) %>% 
  write_tsv(file.path(out.dir, "n.outliers.tsv"))


# ______________________________________________________________________________
# Metrics ====
# ______________________________________________________________________________

# Empty table to collect below metrics from Pearson, Spearman and npsigtest
metrics <- tibble()


# ______________________________________________________________________________
## Correlations with and without outliers ====

# Adding Pearson metrics
# unique(analysis.data$dataset)
cor.pearson.res <- analysis.data %>% 
  group_by(dataset) %>% 
  cor_test(vars = prob, vars2 = is_control) %>% 
  rename(value=cor) %>% 
  mutate(no_outlier=FALSE) %>% 
  select(dataset, value, p, no_outlier, method)
metrics <- bind_rows(metrics, cor.pearson.res)

cor.pearson.res.noutliers <- analysis.data %>%
  filter(!outlier) %>% 
  group_by(dataset) %>% 
  cor_test(vars = prob, vars2 = is_control) %>% 
  rename(value=cor) %>% 
  mutate(no_outlier=TRUE) %>% 
  select(dataset, value, p, no_outlier, method)
metrics <- bind_rows(metrics, cor.pearson.res.noutliers)

# Adding Spearman metrics
cor.spearman.res <- analysis.data %>% 
  group_by(dataset) %>% 
  cor_test(vars = prob, vars2 = is_control, method = "spearman") %>% 
  rename(value=cor) %>% 
  mutate(no_outlier=FALSE) %>% 
  select(dataset, value, p, no_outlier, method)
metrics <- bind_rows(metrics, cor.spearman.res)

cor.spearman.res.noutliers <- analysis.data %>% 
  filter(!outlier) %>% 
  group_by(dataset) %>% 
  cor_test(vars = prob, vars2 = is_control, method = "spearman") %>% 
  rename(value=cor) %>% 
  mutate(no_outlier=TRUE) %>% 
  select(dataset, value, p, no_outlier, method)
metrics <- bind_rows(metrics, cor.spearman.res.noutliers)


# ______________________________________________________________________________
## npsigtest with and without outliers ====

# npsigtest.tidy <- function(np_object){
#   return(data.frame(p=np_object$P, IN=np_object$In))
# }
# 
# npsigtest.res <- analysis.data %>% 
#   nest(data=-dataset) %>% 
#   mutate(fit = map(data, ~ npsigtest(npreg(prob ~ is_control, regtype = "ll", bwmethod = "cv.aic", gradients = TRUE, data = .))),
#          results = map(fit, npsigtest.tidy)) %>% 
#   unnest(results) %>% 
#   rename(value=IN) %>% # bw R2 MSE
#   mutate(method="npsigtest") %>%
#   mutate(no_outlier=FALSE) %>% 
#   select(dataset, value, p, no_outlier, method)
# metrics <- bind_rows(metrics, npsigtest.res)
# 
# npsigtest.res.nooutliers <- analysis.data %>% 
#   filter(!outlier) %>% 
#   nest(data=-dataset) %>% 
#   mutate(fit = map(data, ~ npsigtest(npreg(prob ~ is_control, regtype = "ll", bwmethod = "cv.aic", gradients = TRUE, data = .))),
#          results = map(fit, npsigtest.tidy)) %>% 
#   unnest(results) %>% 
#   rename(value=IN) %>% # bw R2 MSE
#   mutate(method="npsigtest") %>%
#   mutate(no_outlier=TRUE) %>% 
#   select(dataset, value, p, no_outlier, method)
# metrics <- bind_rows(metrics, npsigtest.res.nooutliers)



# ______________________________________________________________________________
## CTDiff with and without outliers ====

ctdiff.res <- analysis.data %>% 
  group_by(dataset) %>% 
  summarise(
    value=ctdiff_value_class(prob, is_control),
    p=ctdiff.permut.test(prob, is_control, permutations=N_PERMUTATIONS),
    no_outlier=FALSE,
    method="ctdiff"
  )
metrics <- bind_rows(metrics, ctdiff.res)

ctdiff.res.nooutliers <- analysis.data %>% 
  filter(!outlier) %>% 
  group_by(dataset) %>% 
  summarise(
    value=ctdiff_value_class(prob, is_control),
    p=ctdiff.permut.test(prob, is_control, permutations=N_PERMUTATIONS),
    no_outlier=TRUE,
    method="ctdiff"
  )
metrics <- bind_rows(metrics, ctdiff.res.nooutliers)


# ______________________________________________________________________________
## Permut Pearson Test with and without outliers ====
  
permut.res <- analysis.data %>% 
  group_by(dataset) %>% 
  summarise(
    value=cor(prob, is_control, method="pearson"),
    # p=permut.test.myfunc(prob, is_control, permutations=N_PERMUTATIONS, myfunc=median),
    p=permut.test.corfunc(prob, is_control, permutations=N_PERMUTATIONS, mymethod="pearson"),
    no_outlier=FALSE,
    method="permut_pearson"
  )
metrics <- bind_rows(metrics, permut.res)

permut.res.nooutliers <- analysis.data %>% 
  filter(!outlier) %>% 
  group_by(dataset) %>% 
  summarise(
    value=cor(prob, is_control, method="pearson"),
    # p=permut.test.myfunc(prob, is_control, permutations=N_PERMUTATIONS, myfunc=median),
    p=permut.test.corfunc(prob, is_control, permutations=N_PERMUTATIONS, mymethod="pearson"),
    no_outlier=TRUE,
    method="permut_pearson"
  )
metrics <- bind_rows(metrics, permut.res.nooutliers)


# ______________________________________________________________________________
## Permut Spearman Test with and without outliers ====

permut.res <- analysis.data %>% 
  group_by(dataset) %>% 
  summarise(
    value=cor(prob, is_control, method="spearman"),
    # p=permut.test.myfunc(prob, is_control, permutations=N_PERMUTATIONS, myfunc=median),
    p=permut.test.corfunc(prob, is_control, permutations=N_PERMUTATIONS, mymethod="spearman"),
    no_outlier=FALSE,
    method="permut_spearman"
  )
metrics <- bind_rows(metrics, permut.res)

permut.res.nooutliers <- analysis.data %>% 
  filter(!outlier) %>% 
  group_by(dataset) %>% 
  summarise(
    value=cor(prob, is_control, method="spearman"),
    # p=permut.test.myfunc(prob, is_control, permutations=N_PERMUTATIONS, myfunc=median),
    p=permut.test.corfunc(prob, is_control, permutations=N_PERMUTATIONS, mymethod="spearman"),
    no_outlier=TRUE,
    method="permut_spearman"
  )
metrics <- bind_rows(metrics, permut.res.nooutliers)


# ______________________________________________________________________________
## Permut Means Test with and without outliers ====
  
permut.res <- analysis.data %>% 
  group_by(dataset) %>% 
  summarise(
    value=permut_value_class(prob, is_control, myfunc=mean),
    p=permut.test.myfunc(prob, is_control, permutations=N_PERMUTATIONS, myfunc=mean),
    no_outlier=FALSE,
    method="permut_means"
  )
metrics <- bind_rows(metrics, permut.res)

permut.res.nooutliers <- analysis.data %>% 
  filter(!outlier) %>% 
  group_by(dataset) %>% 
  summarise(
    value=permut_value_class(prob, is_control, myfunc=mean),
    p=permut.test.myfunc(prob, is_control, permutations=N_PERMUTATIONS, myfunc=mean),
    no_outlier=TRUE,
    method="permut_means"
  )
metrics <- bind_rows(metrics, permut.res.nooutliers)


# ______________________________________________________________________________
## Permut Medians Test with and without outliers ====
  
permut.res <- analysis.data %>% 
  group_by(dataset) %>% 
  summarise(
    value=permut_value_class(prob, is_control, myfunc=median),
    p=permut.test.myfunc(prob, is_control, permutations=N_PERMUTATIONS, myfunc=median),
    no_outlier=FALSE,
    method="permut_medians"
  )
metrics <- bind_rows(metrics, permut.res)

permut.res.nooutliers <- analysis.data %>% 
  filter(!outlier) %>% 
  group_by(dataset) %>% 
  summarise(
    value=permut_value_class(prob, is_control, myfunc=median),
    p=permut.test.myfunc(prob, is_control, permutations=N_PERMUTATIONS, myfunc=median),
    no_outlier=TRUE,
    method="permut_medians"
  )
metrics <- bind_rows(metrics, permut.res.nooutliers)


# ______________________________________________________________________________
# Saving all metrics to file ====
# ______________________________________________________________________________

metrics %>% 
  write_tsv(file.path(out.dir, "metrics.tsv")) #%>% 
  # head()


# ______________________________________________________________________________
## Number of outliers ====

# unique(metrics$method)
# metrics %>% 
#   filter(method=="Pearson") %>% 



# ______________________________________________________________________________
## Biased datasets tables ====

metrics %>% 
  filter(method %in% c("Pearson", "Spearman")) %>%
  mutate(no_outlier_label = if_else(no_outlier, "No.Out.", "Normal")) %>% 
  group_by(method, no_outlier_label) %>% 
  summarise(
    prop.low.val = sum(abs(value)<0.18) / n(),
    prop.high.val = sum(abs(value)>0.3) / n(),
    .groups = "drop"
  ) %>% 
  pivot_wider(names_from = no_outlier_label, values_from = c(prop.low.val, prop.high.val)) %>% 
  write_tsv(file.path(out.dir, "metrics.values.tsv")) #%>%
  # head()
  

# Selection with Spearman
metrics %>% 
  filter(method %in% c("Pearson", "Spearman", "npsigtest")) %>%
  select(-value) %>% 
   left_join(
     metrics %>% 
       filter(method %in% c("Pearson", "Spearman")) %>%
       mutate(value = abs(value)) %>% 
       select(-p) %>%
       pivot_wider(names_from = method, values_from = value)
     ,
     by=c("dataset", "no_outlier")
   ) %>% 
  mutate(Pearson = if_else(Pearson>0.3, 1, 0)) %>% 
  mutate(Spearman = if_else(Spearman>0.3, 1, 0)) %>% 
  mutate(p = if_else(p<0.05, 1, 0)) %>% 
  mutate(intersect_pearson = as.numeric(Pearson+p>1)) %>% 
  mutate(intersect_spearman = as.numeric(Spearman+p>1)) %>% 
  mutate(intersect_ref_spearman = as.numeric(intersect_spearman+Pearson>1)) %>% 
  group_by(no_outlier, method) %>% 
  summarise(
    n = n(),
    n_pos_pval = sum(p),
    perc_pos_pval = n_pos_pval/n(),

    pears_n_pos_val = sum(Pearson),
    pears_perc_pos_val = pears_n_pos_val/n(),
    pears_n_pos_inter = sum(intersect_pearson),
    pears_perc_pos_inter = pears_n_pos_inter/n(),
    
    spear_n_pos_value = sum(Spearman),
    spear_perc_pos_value = spear_n_pos_value/n(),
    spear_n_pos_inter = sum(intersect_spearman),
    spear_perc_pos_inter = spear_n_pos_inter/n(),

    spear_n_ref_inter = sum(intersect_ref_spearman),
    spear_perc_ref_inter = spear_n_ref_inter/n(),

    .groups = "drop"
  ) %>% 
  write_tsv(file.path(out.dir, "metrics.percent.highQI.tsv")) #%>%
  # head()


# Selection with Spearman
metrics %>% 
  filter(method %in% c("Pearson", "Spearman", "npsigtest")) %>%
  select(-value) %>% 
   left_join(
     metrics %>% 
       filter(method %in% c("Pearson", "Spearman")) %>%
       mutate(value = abs(value)) %>% 
       select(-p) %>%
       pivot_wider(names_from = method, values_from = value)
     ,
     by=c("dataset", "no_outlier")
   ) %>% 
  mutate(Pearson = if_else(Pearson<0.18, 1, 0)) %>% 
  mutate(Spearman = if_else(Spearman<0.18, 1, 0)) %>% 
  mutate(p = if_else(p>0.05, 1, 0)) %>% 
  mutate(intersect_pearson = as.numeric(Pearson+p>1)) %>% 
  mutate(intersect_spearman = as.numeric(Spearman+p>1)) %>% 
  mutate(intersect_ref_spearman = as.numeric(intersect_spearman+Pearson>1)) %>% 
  group_by(no_outlier, method) %>% 
  summarise(
    n = n(),
    n_pos_pval = sum(p),
    perc_pos_pval = n_pos_pval/n(),

    pears_n_pos_val = sum(Pearson),
    pears_perc_pos_val = pears_n_pos_val/n(),
    pears_n_pos_inter = sum(intersect_pearson),
    pears_perc_pos_inter = pears_n_pos_inter/n(),
    
    spear_n_pos_value = sum(Spearman),
    spear_perc_pos_value = spear_n_pos_value/n(),
    spear_n_pos_inter = sum(intersect_spearman),
    spear_perc_pos_inter = spear_n_pos_inter/n(),

    spear_n_ref_inter = sum(intersect_ref_spearman),
    spear_perc_ref_inter = spear_n_ref_inter/n(),
    
    
        .groups = "drop"
  ) %>% 
  write_tsv(file.path(out.dir, "metrics.percent.lowQI.tsv")) #%>%
  # head()


# ______________________________________________________________________________
## Plots ====

# Outliers
outliers.plot <- analysis.data %>% 
  select(dataset, n_samples, n_samples_nooutliers) %>% 
  distinct() %>% 
  mutate(diff=n_samples-n_samples_nooutliers) %>% 
  ggplot(aes(x=diff)) + 
  geom_bar() + 
  ylab("datasets") + 
  ggtitle("Number of Outliers")
ggsave(file.path(out.dir, "outliers.png"), outliers.plot, height = 3, width=3)


# P-values
p.plot <- metrics %>% 
  filter(method %in% c("Pearson", "Spearman", "ctdiff", "permut_medians", "permut_means")) %>%
  mutate(no_outlier_label = if_else(no_outlier, "No.Out.", "Normal")) %>% 
  select(-no_outlier, -value) %>%
  group_by(method) %>% 
  pivot_wider(names_from = no_outlier_label, values_from = c(p)) %>% 
  ggplot(aes(x=Normal, y=`No.Out.`)) + 
  geom_jitter() + 
  geom_abline(slope = 1) +
  facet_wrap(vars(method)) +
  ggtitle("Outliers impact on test p-values") + 
  xlab("P-value before outlier removal") +
  ylab("P-value after outlier removal")
# p.plot
ggsave(file.path(out.dir, "pvalues.png"), p.plot, height = 4, width=6)


# Values
val.plot <- metrics %>% 
  filter(method %in% c("Pearson", "Spearman", "ctdiff", "permut_medians", "permut_means")) %>%
  mutate(value = abs(value)) %>% 
  mutate(no_outlier_label = if_else(no_outlier, "No.Out.", "Normal")) %>% 
  select(-no_outlier, -p) %>%
  group_by(method) %>% 
  pivot_wider(names_from = no_outlier_label, values_from = c(value)) %>% 
  ggplot(aes(x=Normal, y=`No.Out.`)) + 
  geom_jitter() + 
  geom_abline(slope = 1) +
  facet_wrap(vars(method), scales="free") + 
  ggtitle("Outliers impact on test statistics or effect sizes") + 
  xlab("Statistics before outlier removal") +
  ylab("Statistics after outlier removal")
# val.plot
ggsave(file.path(out.dir, "values.png"), val.plot, height = 4, width=6)


# pcc.p.plot <- metrics %>% 
#   filter(p>=0) %>%
#   left_join(
#     metrics %>%
#       filter(method=="Pearson") %>% 
#       select(-p, -method) %>% 
#       rename(val.pearson=value) %>% 
#       mutate(val.pearson=abs(val.pearson))
#     , 
#     by=c("dataset", "no_outlier")
#   ) %>% 
#   filter(!method=="Pearson") %>% 
#   mutate(no_outlier_label = if_else(no_outlier, "No.Out.", "Normal")) %>% 
#   ggplot(aes(x=val.pearson, y=p)) + 
#   geom_point() +
#   geom_hline(yintercept = 0.5) +
#   geom_vline(xintercept = 0.18) +
#   geom_vline(xintercept = 0.3) +
#   facet_wrap(vars(method)) +
#   facet_grid(rows=vars(method), cols=vars(no_outlier_label)) +
#   ggtitle("PCC vs P-values")
# 
# ggsave(file.path(out.dir, "pcc.p.png"), pcc.p.plot, height = 8, width=6)


# valpcc.val.plot <- metrics %>% 
#   # filter(p>=0) %>%
# mutate(value = if_else(method=="Spearman", abs(value), value)) %>% 
#     left_join(
#     metrics %>%
#       filter(method=="Pearson") %>% 
#       select(-p, -method) %>% 
#       rename(val.pearson=value) %>% 
#       mutate(val.pearson=abs(val.pearson))
#     # select(-value, -method) %>% 
#     # rename(p.pearson=p)
#     , 
#     by=c("dataset", "no_outlier")
#   ) %>% 
#   filter(!method=="Pearson") %>% 
#   mutate(no_outlier_label = if_else(no_outlier, "No.Out.", "Normal")) %>% 
#   select(-no_outlier, -p) %>%
#   group_by(method) %>% 
#   ggplot(aes(x=val.pearson, y=value)) + 
#   geom_point() +
#   geom_vline(xintercept = 0.18) +
#   geom_vline(xintercept = 0.3) +
#   facet_wrap(vars(method)) +
#   facet_grid(rows=vars(method), cols=vars(no_outlier_label), scales="free_y") +
#   ggtitle("PCC vs Values")
# # valpcc.val.plot
# 
# ggsave(file.path(out.dir, "metrics.val.vs.valpcc.png"), valpcc.val.plot, height = 8, width=6)


# # PCC vs SCC
# # unique(metrics$method)
# pcc.scc.plot <- metrics %>%
#   # filter(method %in% c("Pearson", "Spearman", "npsigtest")) %>% 
#   pivot_wider(names_from = method, values_from = c(p, value)) %>% 
#   # mutate(p_test=p_permut_test) %>%
#   # mutate(p_test=p_npsigtest) %>%
#   # mutate(p_test=p_wilcoxon) %>%
#   mutate(p_test=p_Spearman) %>%
#   mutate(value_Pearson = abs(value_Pearson)) %>% 
#   mutate(value_Spearman = abs(value_Spearman)) %>% 
#   # mutate(QI=if_else(value_Pearson>=0.3, "high", "medium")) %>%
#   # mutate(QI=if_else(value_Pearson<0.18, "low", QI)) %>% 
#   # mutate(QI=if_else(p_test<=0.05, "sign.", "medium")) %>%
#   # mutate(QI=if_else(p_test>=0.5, "low", QI)) %>%
#   mutate(QI=if_else(p_test<=0.05, "significant", "non sign.")) %>%
#   mutate(no_outlier_label = if_else(no_outlier, "No Outliers", "Original")) %>% 
#   select(-no_outlier) %>%
#   left_join(
#     analysis.data %>% 
#       select(dataset, n_samples, n_samples_nooutliers) %>% 
#       distinct(),
#     by="dataset"
#   ) %>% 
#   mutate(n_samples=if_else(no_outlier_label=="Original", n_samples, n_samples_nooutliers)) %>% 
#   select(-n_samples_nooutliers) %>% 
#   mutate(small_data=if_else(n_samples<=15, TRUE, FALSE)) %>% 
#   ggplot(aes(x=value_Pearson, y=value_Spearman, col=QI)) + 
#   geom_point() +
#   geom_hline(yintercept = 0.2) +
#   geom_hline(yintercept = 0.3) +
#   geom_vline(xintercept = 0.2) +
#   geom_vline(xintercept = 0.3) +
#   facet_wrap(vars(no_outlier_label)) +
#   # facet_grid(rows=vars(small_data), cols=vars(no_outlier_label)) +
#   ggtitle("Correlation Coefficients vs Spearman's P-value as QI")
# # pcc.scc.plot
# 
# ggsave(file.path(out.dir, "pcc.scc.png"), pcc.scc.plot, height = 4, width=7)


# ______________________________________________________________________________
# SUMMARY
# ______________________________________________________________________________

# 1. Which measure to choose  (lower RMSE is better)?

## a. robustness to outliers => Spearman (p-value and correlation coefficient)
# see pvalues.png 
# see values.png
# see rmse.tsv (below)
bind_rows(
  metrics %>% 
    mutate(no_outlier_label = if_else(no_outlier, "NoOut", "Original")) %>% 
    select(-no_outlier, -p) %>% 
    # pivot_wider(names_from = no_outlier_label, values_from = c(value, p))
    pivot_wider(names_from = no_outlier_label, values_from = value)  %>% 
    mutate(squared_error = (Original - NoOut)^2) %>% 
    group_by(method) %>% 
    summarise(
      mse = mean(squared_error),
      rmse = sqrt(mse)
    ) %>% 
    mutate(value_type="value"),
  metrics %>% 
    mutate(no_outlier_label = if_else(no_outlier, "NoOut", "Original")) %>% 
    select(-no_outlier, -value) %>% 
    # pivot_wider(names_from = no_outlier_label, values_from = c(value, p))
    pivot_wider(names_from = no_outlier_label, values_from = p)  %>% 
    mutate(squared_error = (Original - NoOut)^2) %>% 
    group_by(method) %>% 
    summarise(
      mse = mean(squared_error),
      rmse = sqrt(mse)
    ) %>% 
    mutate(value_type="p-value")
) %>% 
  filter(!(method=="npsigtest" & value_type=="value")) %>% 
  select(method, value_type, mse, rmse) %>% 
  write_tsv(file.path(out.dir, "rmse.tsv")) %>%
  arrange(rmse) #%>% 
  # head()


# 2. How many datasets in the high-QI group ?
# see metrics.percent.highQI.tsv

metrics.percent.highQI <- read_tsv(file.path(out.dir, "metrics.percent.highQI.tsv"), show_col_types = F)

# using only p-value => 0.35
# metrics.percent.highQI %>% 
#   filter(method=="Spearman" & no_outlier==TRUE) %>% 
#   select(!starts_with("pears")) %>% 
#   select(n_pos_pval, perc_pos_pval)

# using p-value and correl coef => 0.3
# metrics.percent.highQI %>% 
#   filter(method=="Spearman" & no_outlier==TRUE) %>% 
#   select(!starts_with("pears")) %>% 
#   select(spear_n_pos_inter, spear_perc_pos_inter)

# using p-value, correl coef and reference set => 0.3
# metrics.percent.highQI %>% 
#   filter(method=="Spearman" & no_outlier==TRUE) %>% 
#   select(pears_n_pos_val, spear_n_ref_inter)


# 3. How many datasets in the low-QI group ?
# see metrics.percent.lowQI.tsv

metrics.percent.lowQI <- read_tsv(file.path(out.dir, "metrics.percent.lowQI.tsv"), show_col_types = F)

# using only p-value => 0.65
# metrics.percent.lowQI %>% 
#   filter(method=="Spearman" & no_outlier==TRUE) %>% 
#   select(!starts_with("pears")) %>% 
#   select(n_pos_pval, perc_pos_pval)

# using p-value and correl coef => 0.35
# metrics.percent.lowQI %>% 
#   filter(method=="Spearman" & no_outlier==TRUE) %>% 
#   select(!starts_with("pears")) %>% 
#   select(spear_n_pos_inter, spear_perc_pos_inter)

# using p-value, correl coef and reference set => 14, including 11 reference datasets
# metrics.percent.lowQI %>% 
#   filter(method=="Spearman" & no_outlier==TRUE) %>% 
#   select(pears_n_pos_val, spear_n_pos_inter, spear_n_ref_inter)