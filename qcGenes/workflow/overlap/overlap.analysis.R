# ==============================================================================
# INIT
# ==============================================================================
# setwd("~/projects/qcGenes")

# install.packages("ggVennDiagram")
# install.packages("devtools")
# devtools::install_github("yanlinlin82/ggvenn")
# remove.packages("UpSetR")
# devtools::install_github("jonocarroll/UpSetR@678b2bc")
# install.packages("ComplexUpset")
# library(UpSetR)
library(tidyverse)
library(ggpmisc)
library(ggpubr)
library(ComplexUpset)
# library(ggVennDiagram)
# library(ggvenn)
options(ggplot2.discrete.colour= c("#295D8A", "#A41720", "#4E4459"))
options(ggplot2.discrete.fill= c("#295D8A", "#A41720", "#4E4459"))

ref.data.dir <- file.path("data", "ref", "gs2d")
base.dir <- file.path("output", "main")
deg.base.dir <- file.path(base.dir, "RASflowResults")

pos.cor.genes.path <- file.path(base.dir, "qualityCorGenes", "pos.cor.genes.noBias.tsv")
neg.cor.genes.path <- file.path(base.dir, "qualityCorGenes", "neg.cor.genes.noBias.tsv")
datasets.stats.path <- file.path(base.dir, "qualityCorGenes", "datasets.stats.tsv")
disease.genes.path <- file.path(ref.data.dir, "20180416_gene2mesh.ens.tsv.gz")

output.dir <- file.path(base.dir, "overlap")
dir.create(output.dir, showWarnings = F, recursive = T)
output.table.path <- file.path(output.dir, "disease_deg_cor.tsv")
n_dis_genes.table.path <- file.path(output.dir, "n_dis_genes.tsv")
dis.genes.fdr.vs.quality.plot.path <- file.path(output.dir, "figure-7-dis.genes.fdr.vs.quality.plot.pdf")

dis_pos_neg_intersections.path <- file.path(output.dir, "dis_pos_neg_intersections.tsv")
datasets.vs.overlap.path <- file.path(output.dir, "datasets.vs.overlap.pdf")
datasets.vs.overlap.2.path <- file.path(output.dir, "datasets.vs.overlap.small.data.pdf")

datasets.vs.overlap.dis.path <- file.path(output.dir, "datasets.vs.overlap.dis.pdf")
datasets.vs.overlap.dis.2.path <- file.path(output.dir, "datasets.vs.overlap.dis.small.data.pdf")

upset.1.path <- file.path(output.dir, "upset.1.pdf")
upset.2.path <- file.path(output.dir, "upset.2.pdf")
upset.3.path <- file.path(output.dir, "upset.3.pdf")

# n_samples.big    <- config.data$N_SAMPLES_BIG
n_samples.big <- 30


# ==============================================================================
# FUNCTIONS
# ==============================================================================
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


# ==============================================================================
# DATASETS STATS 
# ==============================================================================
datasets.stats <- read_tsv(datasets.stats.path, show_col_types = F) %>% 
  mutate(nickname=paste(dataset, round(p_group_cor, 2), sep="_")) %>% 
  select(dataset, p_group_cor, nickname, n_batches, n_samples, n_degs, n_degs_fdr_only, 
         mesh_terms, SamplesPairing, samples_homogeneity) %>% 
  filter(n_degs>=1000) %>% 
  mutate(mesh_terms=if_else(str_starts(mesh_terms, "Glioma") , "Glioma", mesh_terms)) %>% 
  mutate(mesh_terms=if_else(str_starts(mesh_terms, "Heart") , "Heart Failure", mesh_terms)) %>% 
  mutate(mesh_terms=if_else(str_starts(mesh_terms, "Liver") , "Liver Neoplasms", mesh_terms))
head(datasets.stats) 


# ==============================================================================
# OPEN REFERENCE GENE LISTS
# ==============================================================================
pos.cor.genes <- read_tsv(pos.cor.genes.path, show_col_types = F) %>% 
  filter(n>3)
neg.cor.genes <- read_tsv(neg.cor.genes.path, show_col_types = F) %>% 
  filter(n>3)


# ==============================================================================
# OPEN DISEASE GENES
# ==============================================================================
# All disease genes
disease.genes <- read_tsv(disease.genes.path, show_col_types = F)%>% 
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
length(unique(disease.genes$name))

# Diseases with most genes
disease.genes.1 <- disease.genes %>%
  filter(name_freq_in_table>=300)
length(unique(disease.genes.1$name))

# Diseases in our study
disease.genes.2 <- disease.genes %>% 
  filter(name %in% datasets.stats$mesh_terms) # %>% 
length(unique(disease.genes.2$name))
write_tsv(disease.genes.2, n_dis_genes.table.path)

disease.genes.table <- disease.genes.2 %>%
  group_by(name) %>%
  summarise(n=n())
disease.genes.table
write_tsv(disease.genes.table, n_dis_genes.table.path)

# Diseases in our study with a min number of genes
disease.genes.3 <- disease.genes.2 %>% 
  group_by(name) %>% 
  filter(name_freq_in_table>=300) %>% 
  slice_min(n=300, order_by=p_value, with_ties = F)
length(unique(disease.genes.3$name))
disease.genes.3 %>%
  group_by(name) %>%
  summarise(n=n())


# ==============================================================================
# DISEASE VS RELATION TO QUALITY - PLOTS WITH MANY DISEASES
# ==============================================================================
dis.genes.fdr.vs.quality.plot <- disease.genes.1 %>% 
  filter(count_name_in_gene_set>3) %>% 
  mutate(cor_genes=ifelse(
    cor_genes=="positive", "low", ifelse(
      cor_genes=="negative", "high", "none"
    ))) %>%
  mutate(cor_genes=factor(cor_genes, levels=c("high", "none", "low"))) %>%
  ggplot(aes(x=cor_genes, y=fdr)) +
  geom_boxplot() +
  geom_jitter(alpha=0.05, color="#295D8A") +
  facet_wrap(vars(name), scales = "free") +
  xlab("Quality marker genes") +
  ylab("False Discovery Rate (FDR)") +
  ggtitle("Disease Genes Significance and Relation to Quality")
# dis.genes.fdr.vs.quality.plot
ggsave(dis.genes.fdr.vs.quality.plot.path, dis.genes.fdr.vs.quality.plot, width=14, height=9)


# ==============================================================================
# DISEASE VS QUALITY - GENE LISTS INTERSECTION
# ==============================================================================
dis_pos_neg_intersections <- data.frame()
for (disease in unique(disease.genes.3$name)) {
  my.dis <- disease.genes.3 %>% filter(name==disease)
  dis.n <- nrow(my.dis)
  dis.pos.n <- length(intersect(my.dis$geneid, pos.cor.genes$geneid))
  dis.pos.perc <- dis.pos.n / dis.n
  dis.neg.n <- length(intersect(my.dis$geneid, neg.cor.genes$geneid))
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
  mutate(dis_pos_fisher=my.fisher(dis_pos_inter_n, 300, nrow(pos.cor.genes))) %>% 
  mutate(dis_neg_fisher=my.fisher(dis_neg_inter_n, 300, nrow(neg.cor.genes)))

dis_pos_neg_intersections
write_tsv(dis_pos_neg_intersections, dis_pos_neg_intersections.path)


# ==============================================================================
# DEGS VS QUALITY MARKERS
# ==============================================================================
degs_pos_neg_intersections <- data.frame()
for (i in 1:nrow(datasets.stats) ) {
  
  my.info <- datasets.stats[i,]
  dataset.id <- my.info$dataset
  my.dis <- disease.genes.3 %>% filter(name==my.info$mesh_terms)
  gene.list.dir.path <- file.path(deg.base.dir, dataset.id, "trans", "dea", "DEA", "gene-level")
  gene.list.path <- grep("deg.*", dir(gene.list.dir.path), perl=T, value=T)
  
  my.list <- list()
  for (max_genes in c(50, 100, 250, 500, 1000)) {
    my.degs <- read.delim(file.path(gene.list.dir.path, gene.list.path)) %>%
      rownames_to_column("geneid") %>%
      tibble() %>%
      filter(abs(log2FoldChange)>1) %>% 
      slice_min(order_by = pvalue, n=max_genes)
    deg.n <- nrow(my.degs)
    deg.pos.n <- length(intersect(my.degs$geneid, pos.cor.genes$geneid))
    deg.pos.perc <- deg.pos.n / deg.n
    deg.neg.n <- length(intersect(my.degs$geneid, neg.cor.genes$geneid))
    deg.neg.perc <- deg.neg.n / deg.n
    deg.both.n <- length(intersect(my.degs$geneid, union(pos.cor.genes$geneid, neg.cor.genes$geneid)))
    deg.both.perc <- deg.both.n / deg.n
    deg.dis.n <- length(intersect(my.degs$geneid, my.dis$geneid))
    deg.dis.perc <- deg.dis.n / deg.n
    my.list[[paste0("deg_pos_inter_perc_", max_genes)]] <- deg.pos.perc
    my.list[[paste0("deg_neg_inter_perc_", max_genes)]] <- deg.neg.perc
    my.list[[paste0("deg_both_inter_perc_", max_genes)]] <- deg.both.perc
    my.list[[paste0("deg_dis_inter_perc_", max_genes)]] <- deg.dis.perc
  }
  
  degs_pos_neg_intersections <- rbind(
    degs_pos_neg_intersections, 
    cbind(
      data.frame(
        dataset=my.info$dataset,
        p_group_cor=my.info$p_group_cor,
        n_samples=my.info$n_samples,
        SamplesPairing=my.info$SamplesPairing,
        deg_genes=deg.n, 
        mesh=my.info$mesh_terms),
      as.data.frame(my.list)
    )
  )
}
# degs_pos_neg_intersections
merged.tables <- degs_pos_neg_intersections %>% 
  left_join(dis_pos_neg_intersections %>% select(mesh, dis_genes, dis_pos_inter_perc, dis_neg_inter_perc), by="mesh") 
merged.tables %>% arrange(p_group_cor)

write_tsv(merged.tables, output.table.path)


# PLOT BY DESIGN AND SIZE
datasets.vs.overlap.plot <- merged.tables %>% 
  filter(deg_both_inter_perc_1000>0) %>% 
  mutate(size=if_else(n_samples>=n_samples.big, "big", "small")) %>%
  mutate(design=if_else(SamplesPairing==1, "paired", "unpaired")) %>% 
  ggplot(aes(y=deg_both_inter_perc_1000, x=p_group_cor, col=design)) +
  theme_gray(base_size = 15) +
  geom_point() +
  geom_smooth(formula=y~x, method = "glm") +
  facet_grid(rows=vars(size)) +
  xlab("Quality Imbalance Index") + 
  ylab("Percent Quality Markes in DEGs") +
  ggtitle("Differential Genes vs Quality in Datasets by Design and Size") +
  labs(subtitle = "Data size: small<30 and big\u2265 30 samples") +
  stat_poly_eq(formula = y~x, 
              #  aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
               aes(label = after_stat(rr.label)), 
               size = rel(4),
               parse = TRUE)
# datasets.vs.overlap.plot
ggsave(datasets.vs.overlap.path, datasets.vs.overlap.plot, width=8, height=7)


# PLOT BY DESIGN
datasets.vs.overlap.2.plot <- merged.tables %>% 
  filter(deg_both_inter_perc_1000>0) %>% 
  mutate(size=if_else(n_samples>=50, "big", "small")) %>%
  filter(size=="small") %>%
  mutate(design=if_else(SamplesPairing==1, "paired", "unpaired")) %>% 
  ggplot(aes(y=deg_both_inter_perc_1000, x=p_group_cor)) +
  theme_gray(base_size = 15) +
  geom_point(color="#295D8A") +
  geom_smooth(formula=y~x, method = "glm", color="#295D8A") +
  facet_grid(rows=vars(design)) +
  ylim(0, 0.2) +
  xlab("Quality Imbalance Index") + 
  ylab("Percent Quality Markes in DEGs") +
  ggtitle("Differential Genes vs Quality in Datasets") +
  labs(subtitle = "by Design; Samples<50") +
  stat_poly_eq(formula = y~x, 
               aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
               #aes(label = after_stat(rr.label)), 
               size = rel(4),
               parse = TRUE,
               color="#295D8A") +
  theme_minimal()
# datasets.vs.overlap.2.plot
ggsave(datasets.vs.overlap.2.path, datasets.vs.overlap.2.plot, width=6, height=6)



# ==============================================================================
# DEGS VS QUALITY MARKERS
# ==============================================================================

# DISEASE PLOT BY DESIGN AND SIZE
datasets.vs.overlap.dis.plot <- merged.tables %>% 
  filter(deg_dis_inter_perc_1000>0) %>% 
  mutate(size=if_else(n_samples>=n_samples.big, "big", "small")) %>%
  mutate(design=if_else(SamplesPairing==1, "paired", "unpaired")) %>% 
  ggplot(aes(y=deg_dis_inter_perc_1000, x=p_group_cor, col=design)) +
  theme_gray(base_size = 15) +
  geom_point() +
  geom_smooth(formula=y~x, method = "glm") +
  facet_grid(rows=vars(size)) +
  xlab("Quality Imbalance Index") + 
  ylab("Percent Known Disease Genes in DEGs") +
  ggtitle("Differential vs Known Genes in Datasets") +
  labs(subtitle = "by Design and Size. Data size: small<30 and big\u2265 30 samples") +
  stat_poly_eq(formula = y~x, 
               aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
               # aes(label = after_stat(rr.label)), 
               size = rel(4),
               parse = TRUE)
# datasets.vs.overlap.dis.plot
ggsave(datasets.vs.overlap.dis.path, datasets.vs.overlap.dis.plot, width=8, height=7)


# DISEASE PLOT BY DESIGN
datasets.vs.overlap.dis.2.plot <- merged.tables %>% 
  filter(deg_dis_inter_perc_1000>0) %>% 
  mutate(size=if_else(n_samples>=50, "big", "small")) %>%
  filter(size=="small") %>% 
  mutate(design=if_else(SamplesPairing==1, "paired", "unpaired")) %>% 
  ggplot(aes(y=deg_dis_inter_perc_1000, x=p_group_cor)) +
  theme_gray(base_size = 15) +
  geom_point(color="#295D8A") +
  geom_smooth(formula=y~x, method = "glm", color="#295D8A") +
  facet_grid(rows=vars(design)) +
  xlab("Quality Imbalance Index") + 
  ylab("Percent Known Disease Genes in DEGs") +
  ggtitle("Differential vs Known Genes in Datasets") +
  labs(subtitle = "by Design; Samples<50") +
  stat_poly_eq(formula = y~x, 
               aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
               # aes(label = after_stat(rr.label)), 
               size = rel(4),
               parse = TRUE,
               color="#295D8A") + 
  theme_minimal() 
# datasets.vs.overlap.dis.2.plot
ggsave(datasets.vs.overlap.dis.2.path, datasets.vs.overlap.dis.2.plot, width=6, height=6)

# figure 6 manuscript
figure <- ggarrange(
  datasets.vs.overlap.2.plot,
  datasets.vs.overlap.dis.2.plot,
  labels = c("A", "B"),
  font.label = list(size = 20, color = "black", face = "bold", family = NULL),
  common.legend = TRUE, legend = "bottom",
  ncol = 2, nrow = 1)
path <- file.path(output.dir, "figure-6.pdf")
ggsave(filename=path, plot=figure, width=14, height=7)

# ==============================================================================
# MESH TERMS WITH MULTIPLE DATASETS
# ==============================================================================
selected.terms <- (
  datasets.stats %>%
    filter(n_degs>=1000) %>%
    group_by(mesh_terms) %>%
    summarise(n=n()) %>%
    arrange(desc(n)) %>%
    filter(n>1)
)$mesh_terms

selected.terms <- intersect(selected.terms, disease.genes.3$name)
selected.terms


# ==============================================================================
# DEGs UPSET PLOTS
# ==============================================================================
my.list <- list()
my.plots.list <- list()
for (term in selected.terms) {
  # print(term)
  my.term.datasets <- datasets.stats %>% filter(mesh_terms==term)
  # my.list[[term]][["pos.cor"]] <- pos.cor.genes$geneid
  # my.list[[term]][["neg.cor"]] <- neg.cor.genes$geneid
  my.list[[term]][["quality"]] <- c(pos.cor.genes$geneid, neg.cor.genes$geneid)
  # my.list[[term]][["disease"]] <- disease.genes.3[disease.genes.3$name==term,]$geneid
  for (i in 1:nrow(my.term.datasets)) {
    dataset.id <- my.term.datasets$dataset[i]
    dataset.name <- my.term.datasets$nickname[i]
    # print(paste0("  _ ", dataset.id))
    gene.list.dir.path <- file.path(deg.base.dir, dataset.id, "trans", "dea", "DEA", "gene-level")
    gene.list.path <- grep("deg.*", dir(gene.list.dir.path), perl=T, value=T)
    gene.list <- read.delim(file.path(gene.list.dir.path, gene.list.path)) %>%
      rownames_to_column("geneid") %>%
      tibble() %>%
      slice_min(order_by = pvalue, n=1000)
    genes <- gene.list$geneid
    my.list[[term]][[dataset.name]] <- genes
  }
  my.table <- tibble(geneid=unique(unlist(my.list[[term]])))
  for (item in names(my.list[[term]])) {
    # print(item)
    my.table <- my.table %>% 
      mutate({{item}}:=if_else(geneid %in% my.list[[term]][[item]], 1, 0))
  }
  my.upset <- upset(my.table, intersect = names(my.list[[term]]), wrap=T, min_size=10) +
    ggtitle(term)
  my.plots.list[[term]] <- my.upset
}
# my.plots.list[[1]] / my.plots.list[[2]] / my.plots.list[[3]]

ggsave(upset.1.path, my.plots.list[[1]], width=20, height=10)
ggsave(upset.2.path, my.plots.list[[2]], width=12, height=8)
ggsave(upset.3.path, my.plots.list[[3]], width=12, height=8)


# ==============================================================================
# DISEASE VENN AND UPSET DIAGRAMS
# ==============================================================================

# my.plots.list <- list()
# for(term in selected.terms){
#   print(term)
#   x <- list(PosCor=pos.cor.genes$geneid,
#             NegCor=neg.cor.genes$geneid)
#   x[[term]] <- disease.genes.3[disease.genes.3$name==term,]$geneid
#   my.upset <- upset(fromList(x), order.by = "freq", nsets = length(x), title = term)
#   my.plots.list[[term]] <- my.upset
# }
# my.plots.list # Display UpsetPlots per disease
