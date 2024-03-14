# ______________________________________________________________________________
# INIT ====
# ______________________________________________________________________________
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(ggpmisc, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(yaml, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(stringi, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

options(ggplot2.discrete.colour= c("#295D8A", "#A41720", "#4E4459"))
options(ggplot2.discrete.fill= c("#295D8A", "#A41720", "#4E4459"))

config.path <- file.path(".", "config", "config.yaml")
config <- yaml.load_file(config.path)

MIN_DATASET_N_DEG <- config$MIN_DATASET_N_DEG

MIN_N_DEGS <- config$OVERLAP_MIN_N_DEGS
# MIN_N_DEGS_QUALITY_MARKERS_OVERLAP <- config$OVERLAP_MIN_N_DEGS
# MIN_N_DEGS_DISEASE_MARKERS_OVERLAP <- 300

DIS_GENES <- config$OVERLAP_DIS_GENES

N_SAMPLES_BIG <- config$OVERLAP_N_SAMPLES_BIG
MAX_SAMPLES_FIG6 <- config$OVERLAP_MAX_SAMPLES_FIG6

QI_SIGNIFICANCE_TEST <- config$QI_SIGNIFICANCE_TEST
QI_SIGNIFICANCE_TEST_CUTOFF <- config$QI_SIGNIFICANCE_TEST_CUTOFF
bias.low.cutoff  <- config$BIAS_LOW_CUTOFF
bias.high.cutoff <- config$BIAS_HIGH_CUTOFF
samples.cutoff   <- config$MIN_SAMPLES
# disease.genes.path <- config$GS2D_ENSEMBL_PATH
# DATASETS_FILE <- config$DATASETS_FILE


# ______________________________________________________________________________
# ARGUMENTS ====
# ______________________________________________________________________________
args = commandArgs(trailingOnly=TRUE)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# test.dataids <- read_csv(config$DATASETS_FILE, show_col_types = F) %>% filter(select==1) %>% pull(GEO_Series)
# args <- c(
#   file.path("output", "main", "qualityCorGenes", "pos.cor.genes.noBias.tsv"),
#   file.path("output", "main", "qualityCorGenes", "neg.cor.genes.noBias.tsv"),
#   file.path("output", "main", "qualityCorGenes", "datasets.stats.tsv"),
#   file.path("output", "data", "ref", "gs2d", "20180416_gene2mesh.ens.tsv.gz"),
#   file.path("output", "main", "overlap"),
#   file.path("output", "main", "gene_expression", test.dataids, paste0(test.dataids, ".diff.genes.tsv.gz"))
# )
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

pos.cor.genes.path <- args[1]
neg.cor.genes.path <- args[2]
datasets.stats.path <- args[3]
disease.genes.path <- args[4]
output.dir <- args[5]
gene.list.dir.paths <- args[6:length(args)]

dir.create(output.dir, showWarnings = F, recursive = T)
output.table.path <- file.path(output.dir, "disease_deg_cor.tsv")
n_dis_genes.table.path <- file.path(output.dir, "n_dis_genes.tsv")
dis.genes.fdr.vs.quality.plot.path <- file.path(output.dir, "figure-7-dis.genes.fdr.vs.quality.plot.png")

dis_pos_neg_intersections.path <- file.path(output.dir, "dis_pos_neg_intersections.tsv")
datasets.vs.overlap.path <- file.path(output.dir, "datasets.vs.overlap.png")
datasets.vs.overlap.2.path <- file.path(output.dir, "datasets.vs.overlap.small.data.png")
datasets.vs.overlap.3.path <- file.path(output.dir, "datasets.vs.overlap.dataSize.png")

datasets.vs.overlap.dis.path <- file.path(output.dir, "datasets.vs.overlap.dis.png")
datasets.vs.overlap.dis.2.path <- file.path(output.dir, "datasets.vs.overlap.dis.small.data.png")
datasets.vs.overlap.dis.3.path <- file.path(output.dir, "datasets.vs.overlap.dis.dataSize.png")


# ______________________________________________________________________________
# FUNCTIONS
# ______________________________________________________________________________
conf.mat <- function(inter, n.dis.genes=300, n.quality.genes){
  Confusion <-
    matrix(c(20000-(n.quality.genes-inter + n.dis.genes-inter + inter), n.quality.genes-inter, n.dis.genes-inter, inter),
           nrow = 2,
           dimnames = list(Quality_p = c("no", "yes"),
                           Disease = c("no", "yes")))
  return(Confusion)
}

my.fisher <- function(inter, n.dis.genes=300, n.quality.genes){
  values <- c()
  for (i in inter) {
    my.p <- fisher.test(conf.mat(i, n.dis.genes, n.quality.genes))$p
    values <- c(values, my.p)
  }
  return(values)
}


# ______________________________________________________________________________
# DATASETS STATS 
# ______________________________________________________________________________
datasets.stats <- read_tsv(datasets.stats.path, show_col_types = F) %>% 
  mutate(nickname=paste(dataset, round(p_group_cor, 2), sep="_")) %>% 
  select(dataset, p_group_cor, p_group_test, nickname, n_batches, n_samples, n_degs, n_degs_fdr_only, 
         mesh_terms, SamplesPairing, samples_homogeneity) %>% 
  filter(n_degs>=MIN_DATASET_N_DEG) %>%
  # filter(n_degs>=MIN_N_DEGS) %>%
  mutate(mesh_terms=if_else(str_starts(mesh_terms, "Glioma") , "Glioma", mesh_terms)) %>% 
  mutate(mesh_terms=if_else(str_starts(mesh_terms, "Heart") , "Heart Failure", mesh_terms)) %>% 
  mutate(mesh_terms=if_else(str_starts(mesh_terms, "Liver") , "Liver Neoplasms", mesh_terms))


if(!QI_SIGNIFICANCE_TEST){
  datasets.stats <- datasets.stats %>% 
    mutate(no.bias.group=if_else(p_group_cor<=bias.low.cutoff & n_samples>=samples.cutoff, 1, 0)) %>% 
    mutate(low.bias.group=if_else(p_group_cor<=bias.high.cutoff & n_samples>=samples.cutoff, 1, 0)) %>% 
    mutate(high.bias.group=if_else(p_group_cor>bias.high.cutoff & n_samples>=samples.cutoff, 1, 0)) #%>% 
  # select(matches("bias"), everything())
} else{
  # mutate(no.bias.group=if_else(p_group_test>QI_SIGNIFICANCE_TEST_CUTOFF, 1, 0)) %>% 
  datasets.stats <- datasets.stats %>% 
    mutate(no.bias.group=if_else(p_group_cor<=bias.low.cutoff & n_samples>=samples.cutoff & p_group_test>QI_SIGNIFICANCE_TEST_CUTOFF, 1, 0)) %>% 
    mutate(low.bias.group=if_else(p_group_cor<=bias.high.cutoff & n_samples>=samples.cutoff & p_group_test>QI_SIGNIFICANCE_TEST_CUTOFF, 1, 0)) %>% 
    mutate(high.bias.group=if_else(p_group_cor>bias.high.cutoff & n_samples>=samples.cutoff & p_group_test<=QI_SIGNIFICANCE_TEST_CUTOFF, 1, 0)) #%>% 
  # select(matches("bias"), everything())
}

# datasets.stats <- datasets.stats %>% 
#   filter(no.bias.group==1 | high.bias.group==1)

# ______________________________________________________________________________
# OPEN REFERENCE GENE LISTS
# ______________________________________________________________________________
print("OPEN REFERENCE GENE LISTS")
pos.cor.genes <- read_tsv(pos.cor.genes.path, show_col_types = F) %>% 
  separate(geneid, into = c("geneid", "gene_version")) %>% 
  select(-gene_version) %>% 
  filter(n>3)
neg.cor.genes <- read_tsv(neg.cor.genes.path, show_col_types = F) %>% 
  separate(geneid, into = c("geneid", "gene_version")) %>% 
  select(-gene_version) %>% 
  filter(n>3)


# ______________________________________________________________________________
# OPEN DISEASE GENES
# ______________________________________________________________________________
# All disease genes
disease.genes <- read_tsv(disease.genes.path, show_col_types = F) %>% 
  left_join(pos.cor.genes %>% rename(n_pos_cor=n) %>% select(-genes), by="geneid") %>%
  left_join(neg.cor.genes %>% rename(n_neg_cor=n) %>% select(-genes), by="geneid") %>%
  mutate(is_pos_cor=if_else(!is.na(n_pos_cor), 1, 0)) %>% 
  mutate(is_neg_cor=if_else(!is.na(n_neg_cor), 1, 0)) %>% 
  mutate(cor_genes="none") %>% 
  mutate(cor_genes=if_else(is_pos_cor==1, "positive", cor_genes)) %>% 
  mutate(cor_genes=if_else(is_neg_cor==1, "negative", cor_genes)) %>% 
  mutate(cor_genes=if_else(is_pos_cor==1 & is_neg_cor==1, "both", cor_genes)) %>% 
  filter(name!="Genetic Predisposition to Disease") %>% 
  filter(name!="Disease") %>% 
  filter(name!="Neoplasms") %>%
  group_by(name) %>%
  mutate(name_freq_in_table=n()) %>% 
  ungroup()

# Genes of Diseases with most genes
disease.genes.1 <- disease.genes %>%
  filter(name_freq_in_table>=DIS_GENES)
# length(unique(disease.genes.1$name))
  
# Genes of Diseases in our study
disease.genes.2 <- disease.genes %>% 
  filter(name %in% unique(datasets.stats$mesh_terms)) # %>% 
# length(unique(disease.genes.2$name))
write_tsv(disease.genes.2, n_dis_genes.table.path)

# Diseases in our Study with known disease genes
diseases.table <- disease.genes.2 %>%
  group_by(name) %>%
  summarise(n_genes=n()) %>% 
  ungroup() %>% 
  arrange(desc(n_genes))

# diseases.table
write_tsv(diseases.table, n_dis_genes.table.path)

# Top genes of diseases in our study
disease.genes.3 <- disease.genes.2 %>% 
  filter(name_freq_in_table>=DIS_GENES) %>%
  group_by(name) %>% 
  slice_min(n=DIS_GENES, order_by=p_value, with_ties = F)
# length(unique(disease.genes.3$name))
# disease.genes.3 %>%
#   group_by(name) %>%
#   summarise(n=n())


# ______________________________________________________________________________
# DISEASE VS RELATION TO QUALITY - PLOTS WITH MANY DISEASES
# ______________________________________________________________________________
my.dataset.selection <- sample(unique(disease.genes.1$name), size=20)
dis.genes.fdr.vs.quality.plot <- disease.genes.1 %>% 
  filter(count_name_in_gene_set>3) %>% 
  filter(name %in% my.dataset.selection) %>% 
  # group_by(name) %>% 
  # slice_head(n=200) %>% 
  ggplot(aes(x=cor_genes, y=fdr)) +
  geom_boxplot() +
  geom_jitter(alpha=0.05, color="#295D8A") +
  facet_wrap(vars(name), scales = "free") +
  xlab("Type of correlation with sample low quality") +
  ylab("False Discovery Rate (FDR)") +
  ggtitle("Disease Genes Significance and Relation to Quality")
# dis.genes.fdr.vs.quality.plot
ggsave(dis.genes.fdr.vs.quality.plot.path, dis.genes.fdr.vs.quality.plot, width=14, height=9)


# ______________________________________________________________________________
# DISEASE GENES VS QUALITY MARKERS - ENRICHMENT OF QUALITY GENES IN KNOWN DISEASE GENES
# ______________________________________________________________________________
dis_pos_neg_intersections <- data.frame()
for (disease in unique(disease.genes.3$name)) {
  my.dis <- disease.genes.3 %>% filter(name==disease)
  dis.n <- nrow(my.dis)
  # dis.pos.n <- length(intersect(my.dis$geneid, pos.cor.genes$geneid))
  dis.pos.n <- sum(my.dis$is_pos_cor)
  dis.pos.perc <- dis.pos.n / dis.n
  # dis.neg.n <- length(intersect(my.dis$geneid, neg.cor.genes$geneid))
  dis.neg.n <- sum(my.dis$is_neg_cor)
  dis.neg.perc <- dis.neg.n / dis.n
  dis_pos_neg_intersections <- rbind(dis_pos_neg_intersections, data.frame(
    mesh=disease, 
    dis_genes=dis.n, 
    dis_pos_inter_n=dis.pos.n, 
    dis_pos_inter_perc=dis.pos.perc,
    dis_neg_inter_n=dis.neg.n, 
    dis_neg_inter_perc=dis.neg.perc
    ))
  }
dis_pos_neg_intersections <- dis_pos_neg_intersections %>% 
  mutate(dis_pos_fisher=my.fisher(dis_pos_inter_n, DIS_GENES, nrow(pos.cor.genes))) %>% 
  mutate(dis_neg_fisher=my.fisher(dis_neg_inter_n, DIS_GENES, nrow(neg.cor.genes)))

# dis_pos_neg_intersections
write_tsv(dis_pos_neg_intersections, dis_pos_neg_intersections.path)


# ______________________________________________________________________________
# DEGS VS QUALITY MARKERS
# ______________________________________________________________________________

N_DEGS_LIST <- c(50, 100, 250, 500)
N_DIS_GENES_LIST <- c(50, 100, 250)

print("DEGS VS QUALITY MARKERS")
degs_pos_neg_intersections <- data.frame()
for (i in 1:nrow(datasets.stats) ) {
  
  my.info <- datasets.stats[i,]
  dataset.id <- my.info$dataset
  gene.list.path <- grep(paste0(dataset.id, "/", dataset.id,".diff.genes.tsv.gz"), gene.list.dir.paths, value=T)
  
  dataset.degs <- read.delim(gene.list.path) %>%
    separate(geneid, into = c("geneid", "version")) %>% 
    tibble() %>%
    filter(abs(log2FoldChange)>1 & padj<0.05)
  dataset.deg.n <- nrow(dataset.degs)
  
  dataset.dis.genes <- disease.genes.2 %>% 
    filter(name==my.info$mesh_terms) 
  dataset.dis.genes.n <- nrow(dataset.dis.genes)
  
  if(!identical(gene.list.path, character(0))){
    
    my.list <- list()
    for (max_genes in N_DEGS_LIST) {
      
      if(dataset.deg.n>=max_genes){
        my.degs <- dataset.degs %>% 
          slice_min(order_by = pvalue, n=max_genes)
        deg.n <- nrow(my.degs)
        deg.pos.n <- length(intersect(my.degs$geneid, pos.cor.genes$geneid))
        deg.pos.perc <- deg.pos.n / deg.n
        deg.neg.n <- length(intersect(my.degs$geneid, neg.cor.genes$geneid))
        deg.neg.perc <- deg.neg.n / deg.n
        deg.both.n <- length(intersect(my.degs$geneid, union(pos.cor.genes$geneid, neg.cor.genes$geneid)))
        deg.both.perc <- deg.both.n / deg.n
        my.list[[paste0("deg_pos_inter_perc_", max_genes, "_", 0)]] <- deg.pos.perc
        my.list[[paste0("deg_neg_inter_perc_", max_genes, "_", 0)]] <- deg.neg.perc
        my.list[[paste0("deg_both_inter_perc_", max_genes, "_", 0)]] <- deg.both.perc
        for(n_disease_genes in N_DIS_GENES_LIST){
          if(dataset.dis.genes.n>=n_disease_genes){
            my.dis <- dataset.dis.genes %>% 
              filter(name_freq_in_table>=n_disease_genes) %>%
              group_by(name) %>% 
              slice_min(n=n_disease_genes, order_by=p_value, with_ties = F) 
            dis.genes.n <- nrow(my.dis)
            deg.dis.n <- length(intersect(my.degs$geneid, my.dis$geneid))
            deg.dis.perc <- deg.dis.n / deg.n
            my.list[[paste0("deg_dis_inter_perc_", max_genes, "_", n_disease_genes)]] <- deg.dis.perc
          } else {
            my.list[[paste0("deg_dis_inter_perc_", max_genes, "_", n_disease_genes)]] <- NA
          }
        }
      } else {
        my.list[[paste0("deg_pos_inter_perc_", max_genes, "_", 0)]] <- NA
        my.list[[paste0("deg_neg_inter_perc_", max_genes, "_", 0)]] <- NA
        my.list[[paste0("deg_both_inter_perc_", max_genes, "_", 0)]] <- NA
        for(n_disease_genes in N_DIS_GENES_LIST){
          my.list[[paste0("deg_dis_inter_perc_", max_genes, "_", n_disease_genes)]] <- NA
        }
      }
    }
    
    degs_pos_neg_intersections <- rbind(
      degs_pos_neg_intersections, 
      cbind(
        data.frame(
          dataset=my.info$dataset,
          p_group_cor=my.info$p_group_cor,
          n_samples=my.info$n_samples,
          SamplesPairing=my.info$SamplesPairing,
          dataset_genes=dataset.deg.n,
          mesh_genes=dataset.dis.genes.n,
          mesh=my.info$mesh_terms),
        as.data.frame(my.list)
      )
    )
  }
}
# degs_pos_neg_intersections

merged.tables <- degs_pos_neg_intersections %>% 
  left_join(dis_pos_neg_intersections %>% select(mesh, dis_genes, dis_pos_inter_perc, dis_neg_inter_perc), by="mesh") 
# merged.tables %>% arrange(p_group_cor)

write_tsv(merged.tables, output.table.path)


# ______________________________________________________________________________
# ______________________________________________________________________________
# ______________________________________________________________________________
# REFORMAT merged.tables FOR BETTER PLOTTING
# ______________________________________________________________________________
# ______________________________________________________________________________
# ______________________________________________________________________________

names(degs_pos_neg_intersections)
plot.data <- degs_pos_neg_intersections %>%
  pivot_longer(cols=ends_with("0"), names_sep = "_", names_to = c("genes", "markers", "inter", "perc", "n_genes", "dis_genes")) %>% 
  mutate(genes_label=if_else(dis_genes==0, n_genes, paste(n_genes, dis_genes, sep="_"))) %>% 
  mutate(genes_label=paste(markers, genes_label, sep="_")) %>% 
  mutate(n_genes=as.numeric(n_genes)) %>% 
  mutate(size=if_else(n_samples>=N_SAMPLES_BIG, "big", "small")) %>%
  mutate(design=if_else(SamplesPairing==1, "paired", "unpaired")) %>% 
  select(dataset, p_group_cor, n_samples, dataset_genes, mesh_genes, SamplesPairing, mesh, genes_label, everything()) %>% 
  select(-inter, -perc) %>% 
  # mutate(SamplesPairing=as.factor(SamplesPairing)) %>% 
  tibble()
# plot.data
# sort(unique(plot.data$n_samples))
# plot.data$value

masterplot.degsVmarkers <- plot.data %>% 
  mutate(n_genes=factor(paste(n_genes, "DEGs"), levels=unique(paste(sort(n_genes), "DEGs")))) %>% 
  mutate(markers=if_else(markers=="both", "Quality Markers", markers)) %>% 
  mutate(markers=if_else(markers=="pos", "Low-Quality Markers", markers)) %>% 
  mutate(markers=if_else(markers=="neg", "High-Quality Markers", markers)) %>% 
  filter(!is.na(value)) %>% 
  filter(markers!="dis") %>% 
  filter(size=="small") %>%
  filter(p_group_cor<3) %>% 
  ggplot(aes(x=p_group_cor, y=value, col=design)) + 
  geom_point() + 
  facet_grid(cols=vars(n_genes), rows=vars(markers), scales = "free") +
  geom_smooth(formula=y~x, method = "glm") +
  xlab("Quality Imbalance Index") + 
  ylab("Percent Quality Markes in Differential Genes List") +
  ggtitle("Differential vs Quality Markers") +
  # labs(subtitle = paste0("by Design and Size: small<", N_SAMPLES_BIG, " and big\u2265 ", N_SAMPLES_BIG, " samples")) +
  stat_poly_eq(formula = y~x,
               aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")),
               # aes(label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep = "~~~")),
               # aes(label = after_stat(rr.label)),
               size = rel(2),
               parse = TRUE) +
  theme_minimal()
# masterplot.degsVmarkers
my.path <- file.path("output/main/overlap/masterplot.degsVmarkers.small.datasets.png")
ggsave(my.path, masterplot.degsVmarkers, width=8, height=6)

masterplot.degsVgs2d <- plot.data %>% 
  mutate(n_genes=factor(paste(n_genes, "DEGs"), levels=unique(paste(sort(n_genes), "DEGs")))) %>% 
  mutate(dis_genes=factor(paste(dis_genes, "Disease Genes"), levels=unique(paste(sort(dis_genes), "Disease Genes")))) %>% 
  filter(!is.na(value)) %>% 
  filter(markers=="dis") %>% 
  filter(size=="small") %>%
  # filter(p_group_cor<3) %>% 
  ggplot(aes(x=p_group_cor, y=value, col=design)) + 
  geom_point() + 
  # facet_grid(cols=vars(size, n_genes), rows=vars(dis_genes), scales = "free") +
  facet_grid(cols=vars(n_genes), rows=vars(dis_genes), scales = "free") +
  geom_smooth(formula=y~x, method = "glm") +
  xlab("Quality Imbalance Index") + 
  ylab("Percent Disease Genes in Differential Genes List") +
  ggtitle("Differential vs Disease Genes") +
  # labs(subtitle = paste0("by Design and Size: small<", N_SAMPLES_BIG, " and big\u2265 ", N_SAMPLES_BIG, " samples")) +
  stat_poly_eq(formula = y~x,
               aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")),
               # aes(label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep = "~~~")),
               # aes(label = after_stat(rr.label)),
               size = rel(2),
               parse = TRUE) +
  theme_minimal()
# masterplot.degsVgs2d
my.path <- file.path("output/main/overlap/masterplot.degsVgs2d.small.datasets.png")
ggsave(my.path, masterplot.degsVgs2d, width=8, height=6)
dir.create("output/main/overlap/", recursive = T)
#


# ______________________________________________________________________________
# ______________________________________________________________________________
# ______________________________________________________________________________
# BOXPLOTS REFORMAT
# ______________________________________________________________________________
# ______________________________________________________________________________
# ______________________________________________________________________________

masterplot.degsVmarkers.box <- plot.data %>% 
  left_join(datasets.stats %>% select(dataset, no.bias.group, high.bias.group), by="dataset") %>% 
  filter(no.bias.group + high.bias.group>=1) %>% 
  mutate(group=ifelse(no.bias.group==1, "low", "high")) %>% 
  mutate(group=factor(group, ordered = FALSE)) %>% 
  mutate(group=relevel(group, ref="low")) %>% 
  mutate(n_genes=factor(paste(n_genes, "DEGs"), levels=unique(paste(sort(n_genes), "DEGs")))) %>% 
  mutate(markers=if_else(markers=="both", "Quality Markers", markers)) %>% 
  mutate(markers=if_else(markers=="pos", "Low-Quality Markers", markers)) %>% 
  mutate(markers=if_else(markers=="neg", "High-Quality Markers", markers)) %>% 
  filter(!is.na(value)) %>% 
  filter(markers!="dis") %>% 
  filter(size=="small") %>%
  filter(p_group_cor<3) %>% 
  ggplot(aes(x=group, y=value, col=group)) + 
  geom_boxplot() + 
  facet_grid(cols=vars(design, n_genes), rows=vars(markers), scales = "free") +
  # geom_smooth(formula=y~x, method = "glm") +
  xlab("Quality Imbalance") + 
  ylab("Percent Quality Markes in Differential Genes List") +
  ggtitle("Differential vs Quality Markers") +
  theme_minimal()
masterplot.degsVmarkers.box
my.path <- file.path("output/main/overlap/masterplot.degsVmarkers.box.small.datasets.png")
ggsave(my.path, masterplot.degsVmarkers.box, width=8, height=6)

masterplot.degsVgs2d.box <- plot.data %>% 
  left_join(datasets.stats %>% select(dataset, no.bias.group, high.bias.group), by="dataset") %>% 
  filter(no.bias.group + high.bias.group>=1) %>% 
  mutate(group=ifelse(no.bias.group==1, "low", "high")) %>% 
  mutate(group=factor(group, ordered = FALSE)) %>% 
  mutate(group=relevel(group, ref="low")) %>% 
  mutate(n_genes=factor(paste(n_genes, "DEGs"), levels=unique(paste(sort(n_genes), "DEGs")))) %>% 
  mutate(dis_genes=factor(paste(dis_genes, "Disease Genes"), levels=unique(paste(sort(dis_genes), "Disease Genes")))) %>% 
  filter(!is.na(value)) %>% 
  filter(markers=="dis") %>% 
  filter(size=="small") %>%
  # filter(p_group_cor<3) %>% 
  ggplot(aes(x=group, y=value, col=group)) + 
  geom_boxplot() + 
  # facet_grid(cols=vars(size, n_genes), rows=vars(dis_genes), scales = "free") +
  facet_grid(cols=vars(design, n_genes), rows=vars(dis_genes), scales = "free") +
  # geom_smooth(formula=y~x, method = "glm") +
  xlab("Quality Imbalance Index") + 
  ylab("Percent Disease Genes in Differential Genes List") +
  ggtitle("Differential vs Disease Genes") +
  # labs(subtitle = paste0("by Design and Size: small<", N_SAMPLES_BIG, " and big\u2265 ", N_SAMPLES_BIG, " samples")) +
  # stat_poly_eq(formula = y~x,
  #              aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")),
  #              # aes(label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep = "~~~")),
  #              # aes(label = after_stat(rr.label)),
  #              size = rel(2),
  #              parse = TRUE) +
  theme_minimal()
# masterplot.degsVgs2d.box
my.path <- file.path("output/main/overlap/masterplot.degsVgs2d.box.small.datasets.png")
ggsave(my.path, masterplot.degsVgs2d.box, width=8, height=6)

#



# PLOT BY DESIGN AND SIZE
# MIN_N_DEGS
# MIN_GENES
# MIN_N_DEGS <- 500
# names(merged.tables)
my.variable.1 <- "deg_both_inter_perc_1000"
my.variable.2 <- "deg_dis_inter_perc_1000"
if(MIN_N_DEGS<1000){my.variable.1 <- "deg_both_inter_perc_500_0"; my.variable.2 <-"deg_dis_inter_perc_500"}
if(MIN_N_DEGS<500) {my.variable.1 <- "deg_both_inter_perc_250_0"; my.variable.2 <-"deg_dis_inter_perc_250"}
if(MIN_N_DEGS<250) {my.variable.1 <- "deg_both_inter_perc_100_0"; my.variable.2 <-"deg_dis_inter_perc_100"}
if(MIN_N_DEGS<100) {my.variable.1 <- "deg_both_inter_perc_50_0" ; my.variable.2 <-"deg_dis_inter_perc_50"}

dis_genes_suffix <- ifelse(DIS_GENES %in% N_DIS_GENES_LIST, paste0("_", DIS_GENES), "_50")
my.variable.2 <- paste0(my.variable.2, dis_genes_suffix)

# 
# my.variable.1 <- "deg_neg_inter_perc_1000"
# my.variable.2 <- "deg_dis_inter_perc_1000"
# if(MIN_N_DEGS<1000){my.variable.1 <- "deg_neg_inter_perc_500"; my.variable.2 <-"deg_dis_inter_perc_500"}
# if(MIN_N_DEGS<500) {my.variable.1 <- "deg_neg_inter_perc_250"; my.variable.2 <-"deg_dis_inter_perc_250"}
# if(MIN_N_DEGS<250) {my.variable.1 <- "deg_neg_inter_perc_100"; my.variable.2 <-"deg_dis_inter_perc_100"}
# if(MIN_N_DEGS<100) {my.variable.1 <- "deg_neg_inter_perc_50" ; my.variable.2 <-"deg_dis_inter_perc_50"}
names(merged.tables)
datasets.vs.overlap.plot <- merged.tables %>% 
  mutate(label="datasets") %>% 
  # filter(deg_both_inter_perc_1000>0) %>% 
  filter(get(my.variable.1)>0) %>% 
  mutate(size=if_else(n_samples>=N_SAMPLES_BIG, "big", "small")) %>%
  mutate(design=if_else(SamplesPairing==1, "paired", "unpaired")) %>% 
  # ggplot(aes(y=deg_both_inter_perc_1000, x=p_group_cor, col=design)) +
  ggplot(aes(y=get(my.variable.1), x=p_group_cor, col=design)) +
  theme_gray(base_size = 15) +
  geom_point() +
  geom_smooth(formula=y~x, method = "glm") +
  facet_grid(rows=vars(size), cols=vars(label)) +
  xlab("Quality Imbalance Index") + 
  ylab("Percent Quality Markes in DEGs") +
  ggtitle("Differential Genes vs Quality in Datasets") +
  labs(subtitle = paste0("by Design and Size: small<", N_SAMPLES_BIG, " and big\u2265 ", N_SAMPLES_BIG, " samples")) +
  stat_poly_eq(formula = y~x, 
               aes(label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep = "~~~")),
               # aes(label = after_stat(rr.label)), 
               size = rel(4),
               parse = TRUE) + 
  theme_minimal()
# datasets.vs.overlap.plot
ggsave(datasets.vs.overlap.path, datasets.vs.overlap.plot, width=8, height=7)


# PLOT BY DESIGN
datasets.vs.overlap.2.plot <- merged.tables %>%
  tibble() %>% 
  # filter(p_group_cor<0.6) %>% 
  mutate(label="datasets") %>% 
  # filter(deg_both_inter_perc_1000>0) %>% 
  filter(get(my.variable.1)>0) %>% 
  filter(n_samples<MAX_SAMPLES_FIG6) %>%
  # filter(n_samples<30) %>%
  # mutate(size=if_else(n_samples>=N_SAMPLES_BIG, "big", "small")) %>%
  # filter(size=="small") %>%
  mutate(design=if_else(SamplesPairing==1, "paired", "unpaired")) %>% 
  # ggplot(aes(y=deg_both_inter_perc_1000, x=p_group_cor)) +
  ggplot(aes(y=get(my.variable.1), x=p_group_cor)) +
  theme_gray(base_size = 15) +
  geom_point(color="#295D8A") +
  geom_smooth(formula=y~x, method = "glm", color="#295D8A") +
  facet_grid(rows=vars(design), cols=vars(label)) +
  # ylim(0, 1) +
  xlab("Quality Imbalance Index") + 
  ylab("Percent Quality Markes in DEGs") +
  ggtitle("Differential Genes vs Quality in Datasets") +
  labs(subtitle = paste0("by Design; Samples<", MAX_SAMPLES_FIG6)) +
  stat_poly_eq(formula = y~x, 
               aes(label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep = "~~~")),
               # aes(label = after_stat(rr.label)), 
               size = rel(4),
               parse = TRUE,
               color="#295D8A") +
  theme_minimal()
# datasets.vs.overlap.2.plot
ggsave(datasets.vs.overlap.2.path, datasets.vs.overlap.2.plot, width=6, height=6)


# PLOT BY SIZE
datasets.vs.overlap.3.plot <- merged.tables %>% 
  mutate(label="datasets") %>% 
  # filter(deg_both_inter_perc_1000>0) %>% 
  filter(get(my.variable.1)>0) %>% 
  mutate(size=if_else(n_samples>=N_SAMPLES_BIG, "big", "small")) %>%
  mutate(design=if_else(SamplesPairing==1, "paired", "unpaired")) %>% 
  # ggplot(aes(y=deg_both_inter_perc_1000, x=p_group_cor, col=design)) +
  ggplot(aes(y=get(my.variable.1), x=p_group_cor)) +
  theme_gray(base_size = 15) +
  geom_point() +
  geom_smooth(formula=y~x, method = "glm") +
  facet_grid(rows=vars(size), cols=vars(label)) +
  xlab("Quality Imbalance Index") + 
  ylab("Percent Quality Markes in DEGs") +
  ggtitle("Differential Genes vs Quality in Datasets") +
  labs(subtitle = paste0("by Size: small<", N_SAMPLES_BIG, " and big\u2265 ", N_SAMPLES_BIG, " samples")) +
  stat_poly_eq(formula = y~x, 
               aes(label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep = "~~~")),
               # aes(label = after_stat(rr.label)), 
               size = rel(4),
               parse = TRUE) + 
  theme_minimal()
# datasets.vs.overlap.3.plot
ggsave(datasets.vs.overlap.3.path, datasets.vs.overlap.3.plot, width=8, height=7)


# ______________________________________________________________________________
# DEGS VS DISEASE MARKERS
# ______________________________________________________________________________

# DISEASE PLOT BY DESIGN AND SIZE
datasets.vs.overlap.dis.plot <- merged.tables %>%
  mutate(label="datasets") %>% 
  # filter(deg_dis_inter_perc_1000>0) %>% 
  filter(get(my.variable.2)>0) %>% 
  mutate(data.size=if_else(n_samples>=N_SAMPLES_BIG, "big", "small")) %>%
  mutate(design=if_else(SamplesPairing==1, "paired", "unpaired")) %>% 
  # ggplot(aes(y=deg_dis_inter_perc_1000, x=p_group_cor, col=design)) +
  ggplot(aes(y=get(my.variable.2), x=p_group_cor, col=design)) +
  theme_gray(base_size = 15) +
  geom_point() +
  geom_smooth(formula=y~x, method = "glm") +
  # facet_grid(rows=vars(data.size), cols=vars(label)) +
  facet_grid(rows=vars(data.size), cols=vars(label)) +
  xlab("Quality Imbalance Index") + 
  ylab("Percent Known Disease Genes in DEGs") +
  ggtitle("Differential vs Known Genes in Datasets") +
  # labs(subtitle = "by Design and Size. Data size: small<30 and big\u2265 30 samples") +
  labs(subtitle = paste0("by Design and Size. Data size: small<", N_SAMPLES_BIG, " and big\u2265 ", N_SAMPLES_BIG, " samples")) +
  stat_poly_eq(formula = y~x, 
               aes(label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep = "~~~")),
               # aes(label = after_stat(rr.label)), 
               size = rel(4),
               parse = TRUE) + 
  theme_minimal()
# datasets.vs.overlap.dis.plot
ggsave(datasets.vs.overlap.dis.path, datasets.vs.overlap.dis.plot, width=8, height=7)

# DISEASE PLOT BY DESIGN
datasets.vs.overlap.dis.2.plot <- merged.tables %>% 
  mutate(label="datasets") %>% 
  # filter(deg_dis_inter_perc_1000>0) %>% 
  filter(get(my.variable.2)>0) %>% 
  filter(n_samples<MAX_SAMPLES_FIG6) %>%
  # mutate(size=if_else(n_samples>=N_SAMPLES_BIG, "big", "small")) %>%
  # filter(size=="small") %>% 
  mutate(design=if_else(SamplesPairing==1, "paired", "unpaired")) %>% 
  # ggplot(aes(y=deg_dis_inter_perc_1000, x=p_group_cor)) +
  ggplot(aes(y=get(my.variable.2), x=p_group_cor)) +
  theme_gray(base_size = 15) +
  geom_point(color="#295D8A") +
  geom_smooth(formula=y~x, method = "glm", color="#295D8A") +
  facet_grid(rows=vars(design), cols=vars(label)) +
  xlab("Quality Imbalance Index") + 
  ylab("Percent Known Disease Genes in DEGs") +
  ggtitle("Differential vs Known Genes in Datasets") +
  # labs(subtitle = "by Design; Samples<50") +
  labs(subtitle = paste0("by Design; Samples<", MAX_SAMPLES_FIG6)) +
  stat_poly_eq(formula = y~x, 
               aes(label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep = "~~~")),
               # aes(label = after_stat(rr.label)), 
               size = rel(4),
               parse = TRUE,
               color="#295D8A") + 
  theme_minimal()
# datasets.vs.overlap.dis.2.plot
ggsave(datasets.vs.overlap.dis.2.path, datasets.vs.overlap.dis.2.plot, width=6, height=6)

# DISEASE PLOT BY SIZE
datasets.vs.overlap.dis.3.plot <- merged.tables %>%
  mutate(label="datasets") %>% 
  # filter(deg_dis_inter_perc_1000>0) %>% 
  filter(get(my.variable.2)>0) %>% 
  mutate(data.size=if_else(n_samples>=N_SAMPLES_BIG, "big", "small")) %>%
  mutate(design=if_else(SamplesPairing==1, "paired", "unpaired")) %>% 
  # ggplot(aes(y=deg_dis_inter_perc_1000, x=p_group_cor, col=design)) +
  ggplot(aes(y=get(my.variable.2), x=p_group_cor)) +
  theme_gray(base_size = 15) +
  geom_point() +
  geom_smooth(formula=y~x, method = "glm") +
  facet_grid(rows=vars(data.size), cols=vars(label)) +
  # facet_grid(data.size ~ design, scales="free") +
  xlab("Quality Imbalance Index") + 
  ylab("Percent Known Disease Genes in DEGs") +
  ggtitle("Differential vs Known Genes in Datasets") +
  # labs(subtitle = "by Design and Size. Data size: small<30 and big\u2265 30 samples") +
  labs(subtitle = paste0("by Size. Data size: small<", N_SAMPLES_BIG, " and big\u2265 ", N_SAMPLES_BIG, " samples")) +
  stat_poly_eq(formula = y~x, 
               aes(label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep = "~~~")),
               # aes(label = after_stat(rr.label)), 
               size = rel(4),
               parse = TRUE) + 
  theme_minimal()
# datasets.vs.overlap.dis.3.plot
ggsave(datasets.vs.overlap.dis.3.path, datasets.vs.overlap.dis.3.plot, width=8, height=7)


# figure 6 manuscript
figure <- ggarrange(
  datasets.vs.overlap.2.plot,
  datasets.vs.overlap.dis.2.plot,
  labels = c("A", "B"),
  font.label = list(size = 20, color = "black", face = "bold", family = NULL),
  common.legend = TRUE, legend = "bottom",
  ncol = 2, nrow = 1)
path <- file.path(output.dir, "figure.6.png")
# figure
ggsave(filename=path, plot=figure, width=10, height=7)


# figure 6_ALT manuscript
figure <- ggarrange(
  datasets.vs.overlap.plot,
  datasets.vs.overlap.dis.plot,
  labels = c("A", "B"),
  font.label = list(size = 20, color = "black", face = "bold", family = NULL),
  common.legend = TRUE, legend = "bottom",
  ncol = 2, nrow = 1)
path <- file.path(output.dir, "figure.6.size.and.design.png")
# figure
ggsave(filename=path, plot=figure, width=10, height=7)


# figure 6_ALT manuscript
figure <- ggarrange(
  datasets.vs.overlap.3.plot,
  datasets.vs.overlap.dis.3.plot,
  labels = c("A", "B"),
  font.label = list(size = 20, color = "black", face = "bold", family = NULL),
  common.legend = TRUE, legend = "bottom",
  ncol = 2, nrow = 1)
path <- file.path(output.dir, "figure.6.by.size.png")
# figure
ggsave(filename=path, plot=figure, width=10, height=7)

