#!/usr/bin/env Rscript
# ==============================================================================
# Collect and analyse genes correlated with quality in different datasets 
# Output directory: output/qualityCorGenes/
# ==============================================================================

# TO DO: 
#  * FGSEA: 1 dataset per mesh term
# ==============================================================================
# LIBRARIES AND SYSTEM OPTIONS
# ==============================================================================
suppressMessages(library(fgsea, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(yaml, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(ggpmisc, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

options(readr.num_columns = 0)
options(ggplot2.discrete.colour= c("#295D8A", "#A41720", "#4E4459"))
options(ggplot2.discrete.fill= c("#295D8A", "#A41720", "#4E4459"))


# ==============================================================================
# GLOBAL VARIABLES  
# ==============================================================================
config.path <- file.path(".", "config", "config.yaml")
config.data <- yaml.load_file(config.path)
#msigdb.kegg.path <- config.data$MSIGDB_PATH
msigdb.hallmark.path <- config.data$MSIGDB_HALLMARK_URL
msigdb.posi.path <- config.data$MSIGDB_POSITIONAL_PATH
msigdb.curated.path <- config.data$MSIGDB_CURATED_PATH
msigdb.curated.cp.path <- config.data$MSIGDB_CURATED_CP_PATH
msigdb.regulatory.path <- config.data$MSIGDB_REGULATORY_PATH
msigdb.cells.path <- config.data$MSIGDB_CELLS_PATH
msigdb.gs2d.path <- config.data$GS2D_MSIGDB_PATH

sra.path    <- file.path(".", "config", "metadata", "Mega_SraRunTable.csv")
datasets.path <- file.path(".", "config", "metadata", "Datasets.csv")
out.dir.path <- file.path("output", "main", "qualityCorGenes")
fgsea.dir.path <- file.path(out.dir.path, "fgsea")
pca.dir.path <- file.path(out.dir.path, "pca")
gsea.input_genes <- config.data$GSEA_INPUT_GENES
gsea.pathway.perc.cutoff <- config.data$GSEA_PATHWAY_REPRESENTATION_PERCENTAGE_CUTOFF
gsea.input.perc.cutoff <- config.data$GSEA_INPUTGENES_REPRESENTATION_PERCENTAGE_CUTOFF
gsea.table.cutoff <- config.data$GSEA_TABLE_PADJ_CUTOFF
gsea.plot.cutoff <- config.data$GSEA_PLOT_PADJ_CUTOFF2
dataset_p_group_cor_cutoff <- config.data$MAX_DATASET_P_VS_GROUP_CORRELATION
dataset_n_deg_cutoff <- config.data$MIN_DATASET_N_DEG

dir.create(out.dir.path, showWarnings=FALSE, recursive=TRUE)
dir.create(fgsea.dir.path, showWarnings=FALSE, recursive=TRUE)
dir.create(pca.dir.path, showWarnings=FALSE, recursive=TRUE)

bias.low.cutoff  <- config.data$BIAS_LOW_CUTOFF
bias.high.cutoff <- config.data$BIAS_HIGH_CUTOFF
#bias.low.cutoff <- 0.20
#bias.high.cutoff <- 0.40
#bias.low.cutoff <- 0.05
#bias.high.cutoff <- 0.1
#samples.cutoff <- 10
#n_samples.medium <- 15
#n_samples.big <- 30
samples.cutoff   <- config.data$MIN_SAMPLES
n_samples.medium <- config.data$N_SAMPLES_MEDIUM
n_samples.big    <- config.data$N_SAMPLES_BIG
near.zero.cutoff <- 0.05

# ==============================================================================
# GLOBAL DATA
# ==============================================================================
sra.table <- read_csv(sra.path, show_col_types=F) %>% 
  rename(sample=Run)
# dataids.list <- unique(sra.table$GEO_Series)
datasets.table <- read_csv(datasets.path, show_col_types=F) %>% 
  rename(mesh_terms=`mesh_terms (pipe-separated list)`)
dataids.list <- (datasets.table %>% filter(select==1 & batches==0))$GEO_Series
#dataids.list <- (datasets.table %>% filter(select==1 & batches==0 & curators=="JF"))$GEO_Series
# remove subsets
dataids.list <- dataids.list[!grepl("sub", dataids.list)]
# datasets deg
# ------------------------------------------------------------------------------
# deg.stats <- c()
# for(DATAID in dataids.list){
#   # print(paste0("quality_correlated_gemes.r: dataset ", DATAID))
#   cotrl <- (datasets.table %>% filter(GEO_Series==DATAID))$Control
#   treat <- (datasets.table %>% filter(GEO_Series==DATAID))$Treat
#   file <- file.path("output", "main", "RASflowResults", DATAID, "trans", "dea", "DEA", 
#                     "gene-level", paste0("deg_", cotrl, "_", treat, ".tsv") )
#   if (file.exists(file)){
#     deg.stats <- append(deg.stats, length(readLines(file)))
#   }
#   else{
#     deg.stats <- append(deg.stats, -1)
#   }
# }
# deg.stats.table <- as_tibble(cbind(dataset=dataids.list, deg=deg.stats)) %>% 
#   mutate(deg=as.numeric(deg)-1)

# datasets stats
# ------------------------------------------------------------------------------
datasets.stats0 <- tibble()
for(DATAID in dataids.list){
  file <- file.path("output", "main", "qualityVsExp", DATAID, paste0(DATAID, ".cor.tsv"))
  if (file.exists(file)){
    datasets.stats0 <- bind_rows(list(datasets.stats0, read_tsv(file, show_col_types=F)))
  }
}
# datasets.stats0 <- datasets.stats
datasets.stats <- datasets.stats0 %>% 
  # left_join(deg.stats.table, by="dataset") %>% 
  mutate(n_degs=if_else(is.na(n_degs), 0, n_degs)) %>%
  mutate(n_degs_fdr_only=if_else(is.na(n_degs_fdr_only), 0, n_degs_fdr_only)) %>%
  filter(n_degs>0) %>% 
  # mutate(selected=!(n_degs<dataset_n_deg_cutoff | p_group_cor>=dataset_p_group_cor_cutoff))
  mutate(selected=n_degs>=dataset_n_deg_cutoff & p_group_cor<dataset_p_group_cor_cutoff) %>% 
  left_join(datasets.table %>% rename(dataset=GEO_Series) %>% select(dataset, mesh_terms, balanced_groups, SamplesPairing, samples_homogeneity), by = "dataset")

path <- file.path(out.dir.path, paste0("datasets.stats.tsv"))
write_tsv(datasets.stats, path)  

#bias.low.cutoff  <- round(quantile(datasets.stats$p_group_cor, 0.25), 2)
#bias.high.cutoff <- round(quantile(datasets.stats$p_group_cor, 0.50), 2)

# correlation QI vs IF
# IF <- datasets.table %>% filter(GEO_Series %in% dataids.list) %>% select(GEO_Series, `2-YEAR-IF2020`)
# IF <- IF %>% left_join(datasets.stats %>% select(dataset, p_group_cor), by=c("GEO_Series"="dataset"))
# cor(IF$p_group_cor, IF$`2-YEAR-IF2020`, use="pairwise.complete.obs")

# ==============================================================================
# CORRELATED GENES IN LOW BIAS DATASETS
# ==============================================================================
#dataids.list <- (datasets.stats %>% filter(selected==TRUE))$dataset
dataids.list <- (datasets.stats %>% 
                   filter(p_group_cor<=bias.high.cutoff) %>% 
                   filter(n_samples>=samples.cutoff) %>% 
                   arrange(p_group_cor) %>% 
                   group_by(mesh_terms) %>% 
                   slice_head(n=1)
)$dataset


# POSITIVELY CORRELATED GENES - LOW BIAS
# ==============================================================================
d <- c()
for(DATAID in dataids.list){
  # print(paste0("quality_correlated_gemes.r: dataset ", DATAID))
  file <- file.path("output", "main", "qualityVsExp", DATAID, paste0(DATAID, ".cor_pos_genes.tsv"))
  if (file.exists(file)){
    d <- append(d, read_tsv(file, show_col_types=F)$SYMBOL)
  }
}
pos_genes <- tibble(genes=d) %>% 
  group_by(genes) %>% 
  summarise(n=n()) %>% 
  filter(n>=2 & genes!="NA") %>% 
  arrange(desc(n))
out.file1 <- file.path(out.dir.path, "pos.cor.genes.tsv")
write_tsv(pos_genes, out.file1)

path <- file.path(out.dir.path, paste0("pos.cor.genes.pdf"))
pos.cor.genes <- pos_genes %>% 
  rename(gene=genes, datasets=n) %>% 
  # top_n(20,datasets) %>%
  arrange(desc(datasets)) %>% 
  slice_head(n=25) %>% 
  ggplot(aes(x=reorder(gene, datasets), y=datasets)) +
  geom_bar(stat='identity', fill="#295D8A") +
  labs(title="Top pos. correlated genes in low-bias datasets",
       subtitle=paste0("Datasets selection (n=", length(dataids.list), "): bias\u22640.2 and samples\u2265", samples.cutoff), 
       x="genes") +
  coord_flip() +
  # theme_minimal(base_size = 20)
  theme_minimal()
ggsave(filename=path, plot=pos.cor.genes, device='pdf', width=5, height=5)

path <- file.path(out.dir.path, paste0("pos.cor.all.genes.pdf"))
pos.cor.all.genes <- pos_genes %>% 
  rename(gene=genes, datasets=n) %>% 
  ggplot(aes(x=reorder(gene, datasets), y=datasets)) +
  geom_bar(stat="identity", fill="#295D8A") +
  labs(title=paste("Distribution of", nrow(pos_genes) ,"pos. correlated genes in low-bias datasets"), 
       subtitle=paste0("Datasets selection (n=", length(dataids.list), "): bias\u22640.2 and samples\u2265", samples.cutoff), 
       x="genes") +
  coord_flip() +
  theme_minimal(base_size = 20) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
ggsave(filename=path, plot=pos.cor.all.genes, device='pdf', width=10, height=10)


# NEGATIVELY CORRELATED GENES - LOW BIAS
# ==============================================================================
d <- c()
for(DATAID in dataids.list){
  file <- file.path("output", "main", "qualityVsExp", DATAID, paste0(DATAID, ".cor_neg_genes.tsv"))
  if (file.exists(file)){
    d <- append(d, read_tsv(file, show_col_types=F)$SYMBOL)
  }
}
neg_genes <- tibble(genes=d) %>% 
  group_by(genes) %>% 
  summarise(n=n()) %>% 
  filter(n>=2 & genes!="NA") %>% 
  arrange(desc(n))
out.file2 <- file.path(out.dir.path, "neg.cor.genes.tsv")
write_tsv(neg_genes, out.file2)

path <- file.path(out.dir.path, paste0("neg.cor.genes.pdf"))
neg.cor.genes <- neg_genes %>% 
  rename(gene=genes, datasets=n) %>% 
  # top_n(20,datasets) %>%
  arrange(desc(datasets)) %>% 
  slice_head(n=25) %>% 
  # arrange(datasets) %>% 
  ggplot(aes(x=reorder(gene, datasets), y=datasets)) +
  geom_bar(stat='identity', fill="#295D8A") +
  labs(title="Top neg. correlated genes in low-bias datasets", 
       subtitle=paste0("Datasets selection (n=", length(dataids.list), "): bias\u22640.2 and samples\u2265", samples.cutoff), 
       x="genes") +
  coord_flip() +
  theme_minimal(base_size = 20)
ggsave(filename=path, plot=neg.cor.genes, device='pdf', width=10, height=10)

path <- file.path(out.dir.path, paste0("neg.cor.all.genes.pdf"))
neg.cor.all.genes <- neg_genes %>% 
  rename(gene=genes, datasets=n) %>% 
  ggplot(aes(x=reorder(gene, datasets), y=datasets)) +
  geom_bar(stat='identity', fill="#295D8A") +
  labs(title=paste("Distribution of", nrow(neg_genes) ,"neg. correlated genes in low-bias datasets"), 
       subtitle=paste0("Datasets selection (n=", length(dataids.list), "): bias\u22640.2 and samples\u2265", samples.cutoff), 
       x="genes") +
  coord_flip() +
  theme_minimal(base_size = 20) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
ggsave(filename=path, plot=neg.cor.all.genes, device='pdf', width=12, height=10)


# ==============================================================================
# CORRELATED GENES IN NO BIAS DATASETS
# ==============================================================================
#dataids.list <- (datasets.stats %>% filter(selected==TRUE))$dataset
# dataids.list <- (datasets.stats %>% filter(p_group_cor<=bias.low.cutoff) %>% filter(n_samples>=samples.cutoff))$dataset
dataids.list <- (datasets.stats %>% 
                   filter(p_group_cor<=bias.low.cutoff) %>% 
                   filter(n_samples>=samples.cutoff) %>% 
                   arrange(p_group_cor) %>% 
                   group_by(mesh_terms) %>% 
                   slice_head(n=1)
)$dataset


# POSITIVELY CORRELATED GENES - NO BIAS
# ==============================================================================
d <- c()
d.ids <- c()
for(DATAID in dataids.list){
  # print(paste0("quality_correlated_gemes.r: dataset ", DATAID))
  file <- file.path("output", "main", "qualityVsExp", DATAID, paste0(DATAID, ".cor_pos_genes.tsv"))
  if (file.exists(file)){
    d <- append(d, read_tsv(file, show_col_types=F)$SYMBOL)
    d.ids <- append(d.ids, read_tsv(file, show_col_types=F)$geneid)
  }
}
pos_genes_noBias <- tibble(genes=d, geneid=d.ids) %>% 
  group_by(genes, geneid) %>% 
  summarise(n=n()) %>% 
  filter(n>=2 & genes!="NA") %>% 
  arrange(desc(n)) %>% 
  ungroup()
out.file1 <- file.path(out.dir.path, "pos.cor.genes.noBias.tsv")
write_tsv(pos_genes_noBias, out.file1)

path <- file.path(out.dir.path, paste0("pos.cor.genes.noBias.pdf"))
pos.cor.genes.noBias <- pos_genes_noBias %>% 
  select(-geneid) %>%
  rename(gene=genes, datasets=n) %>% 
  # top_n(20,datasets) %>%
  arrange(desc(datasets)) %>% 
  slice_head(n=25) %>% 
  ggplot(aes(x=reorder(gene, datasets), y=datasets)) +
  geom_bar(stat='identity', fill="#295D8A") +
  labs(title="Top low-quality marker genes",
      #  subtitle=paste0("Datasets selection (n=", length(dataids.list), "): bias\u2264", bias.low.cutoff," and samples\u2265", samples.cutoff), 
       x="Quality Imbalance Index") +
  coord_flip() +
  theme_minimal()+
  scale_y_continuous(breaks=c(0,2,4,6,8,10))
ggsave(filename=path, plot=pos.cor.genes.noBias, device='pdf', width=5, height=5)

path <- file.path(out.dir.path, paste0("pos.cor.all.genes.noBias.pdf"))
pos.cor.all.genes.noBias <- pos_genes_noBias %>% 
  rename(gene=genes, datasets=n) %>% 
  ggplot(aes(x=reorder(gene, datasets), y=datasets)) +
  geom_bar(stat="identity", fill="#295D8A") +
  labs(title=paste("Distribution of", nrow(pos_genes_noBias) ,"pos. correlated genes in no-bias datasets"), 
      #  subtitle=paste0("Datasets selection (n=", length(dataids.list), "): bias\u2264", bias.low.cutoff," and samples\u2265", samples.cutoff), 
       x="genes") +
  coord_flip() +
  theme_minimal(base_size = 20) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  scale_y_continuous(breaks=c(0,2,4,6,8,10))
ggsave(filename=path, plot=pos.cor.all.genes.noBias, device='pdf', width=12, height=10)


# NEGATIVELY CORRELATED GENES - NO BIAS
# ==============================================================================
d <- c()
d.ids <- c()
for(DATAID in dataids.list){
  file <- file.path("output", "main", "qualityVsExp", DATAID, paste0(DATAID, ".cor_neg_genes.tsv"))
  if (file.exists(file)){
    d <- append(d, read_tsv(file, show_col_types=F)$SYMBOL)
    d.ids <- append(d.ids, read_tsv(file, show_col_types=F)$geneid)
  }
}
neg_genes_noBias <- tibble(genes=d, geneid=d.ids) %>% 
  group_by(genes, geneid) %>% 
  summarise(n=n()) %>% 
  filter(n>=2 & genes!="NA") %>% 
  arrange(desc(n))%>% 
  ungroup()
out.file2 <- file.path(out.dir.path, "neg.cor.genes.noBias.tsv")
write_tsv(neg_genes_noBias, out.file2)

path <- file.path(out.dir.path, paste0("neg.cor.genes.noBias.pdf"))
neg.cor.genes.noBias <- neg_genes_noBias %>% 
  rename(gene=genes, datasets=n) %>% 
  # top_n(20,datasets) %>%
  arrange(desc(datasets)) %>% 
  slice_head(n=25) %>% 
  arrange(datasets) %>% 
  ggplot(aes(x=reorder(gene, datasets), y=datasets)) +
  geom_bar(stat='identity', fill="#295D8A") +
  labs(title="Top high-quality marker genes", 
      #  subtitle=paste0("Datasets selection (n=", length(dataids.list), "): bias\u2264", bias.low.cutoff," and samples\u2265", samples.cutoff), 
       x="Quality Imbalance Index") +
  coord_flip() +
  theme_minimal()
ggsave(filename=path, plot=neg.cor.genes.noBias, device='pdf', width=5, height=5)

path <- file.path(out.dir.path, paste0("neg.cor.all.genes.noBias.pdf"))
neg.cor.all.genes.noBias <- neg_genes_noBias %>% 
  rename(gene=genes, datasets=n) %>% 
  ggplot(aes(x=reorder(gene, datasets), y=datasets)) +
  geom_bar(stat='identity', fill="#295D8A") +
  labs(title=paste("Distribution of", nrow(neg_genes_noBias) ,"neg. correlated genes in no-bias datasets"), 
      #  subtitle=paste0("Datasets selection (n=", length(dataids.list), "): bias\u2264", bias.low.cutoff," and samples\u2265", samples.cutoff), 
       x="genes") +
  coord_flip() +
  theme_minimal(base_size = 20) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
ggsave(filename=path, plot=neg.cor.all.genes.noBias, device='pdf', width=12, height=10)



# ==============================================================================
# CORRELATED GENES IN HIGH BIAS DATASETS
# ==============================================================================
#dataids.list <- (datasets.stats %>% filter(selected==TRUE))$dataset
# dataids.list <- (datasets.stats %>% filter(p_group_cor>bias.high.cutoff) %>% filter(n_samples>=samples.cutoff))$dataset
dataids.list <- (datasets.stats %>% 
                   filter(p_group_cor>bias.high.cutoff) %>% 
                   filter(n_samples>=samples.cutoff) %>% 
                   arrange(p_group_cor) %>% 
                   group_by(mesh_terms) %>% 
                   slice_head(n=1)
)$dataset


# POSITIVELY CORRELATED GENES - HIGH BIAS
# ==============================================================================
d <- c()
for(DATAID in dataids.list){
  # print(paste0("quality_correlated_gemes.r: dataset ", DATAID))
  file <- file.path("output", "main", "qualityVsExp", DATAID, paste0(DATAID, ".cor_pos_genes.tsv"))
  if (file.exists(file)){
    d <- append(d, read_tsv(file, show_col_types=F)$SYMBOL)
  }
}
pos_genes_highBias <- tibble(genes=d) %>% 
  group_by(genes) %>% 
  summarise(n=n()) %>% 
  filter(n>=2 & genes!="NA") %>% 
  arrange(desc(n))
out.file1 <- file.path(out.dir.path, "pos.cor.genes.highBias.tsv")
write_tsv(pos_genes_highBias, out.file1)

path <- file.path(out.dir.path, paste0("pos.cor.genes.highBias.pdf"))
pos.cor.genes.highBias <- pos_genes_highBias %>% 
  rename(gene=genes, datasets=n) %>% 
  # top_n(20,datasets) %>%
  arrange(desc(datasets)) %>% 
  slice_head(n=25) %>% 
  ggplot(aes(x=reorder(gene, datasets), y=datasets)) +
  geom_bar(stat='identity', fill="#295D8A") +
  labs(title="Top pos. correlated genes in high-bias datasets",
      #  subtitle=paste0("Datasets selection (n=", length(dataids.list), "): bias>", bias.high.cutoff," and samples\u2265", samples.cutoff), 
       x="genes") +
  coord_flip() +
  theme_minimal()
ggsave(filename=path, plot=pos.cor.genes.highBias, device='pdf', width=5, height=5)

path <- file.path(out.dir.path, paste0("pos.cor.all.genes.highBias.pdf"))
pos.cor.all.genes.highBias <- pos_genes_highBias %>% 
  rename(gene=genes, datasets=n) %>% 
  ggplot(aes(x=reorder(gene, datasets), y=datasets)) +
  geom_bar(stat="identity", fill="#295D8A") +
  labs(title=paste("Distribution of", nrow(pos_genes_highBias) ,"pos. correlated genes in high-bias datasets"), 
      #  subtitle=paste0("Datasets selection (n=", length(dataids.list), "): bias>", bias.high.cutoff," and samples\u2265", samples.cutoff), 
       x="genes") +
  coord_flip() +
  theme_minimal(base_size = 20) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
ggsave(filename=path, plot=pos.cor.all.genes.highBias, device='pdf', width=12, height=10)


# NEGATIVELY CORRELATED GENES - HIGH BIAS
# ==============================================================================
d <- c()
for(DATAID in dataids.list){
  file <- file.path("output", "main", "qualityVsExp", DATAID, paste0(DATAID, ".cor_neg_genes.tsv"))
  if (file.exists(file)){
    d <- append(d, read_tsv(file, show_col_types=F)$SYMBOL)
  }
}
neg_genes_highBias <- tibble(genes=d) %>% 
  group_by(genes) %>% 
  summarise(n=n()) %>% 
  filter(n>=2 & genes!="NA") %>% 
  arrange(desc(n))
out.file2 <- file.path(out.dir.path, "neg.cor.genes.highBias.tsv")
write_tsv(neg_genes_highBias, out.file2)

path <- file.path(out.dir.path, paste0("neg.cor.genes.highBias.pdf"))
neg.cor.genes.highBias <- neg_genes_highBias %>% 
  rename(gene=genes, datasets=n) %>% 
  # top_n(20,datasets) %>%
  arrange(desc(datasets)) %>% 
  slice_head(n=25) %>% 
  arrange(datasets) %>% 
  ggplot(aes(x=reorder(gene, datasets), y=datasets)) +
  geom_bar(stat='identity', fill="#295D8A") +
  labs(title="Top neg. correlated genes in high-bias datasets", 
      #  subtitle=paste0("Datasets selection (n=", length(dataids.list), "): bias>", bias.high.cutoff," and samples\u2265", samples.cutoff), 
       x="genes") +
  coord_flip() +
  theme_minimal()
ggsave(filename=path, plot=neg.cor.genes.highBias, device='pdf', width=5, height=5)

path <- file.path(out.dir.path, paste0("neg.cor.all.genes.highBias.pdf"))
neg.cor.all.genes.highBias <- neg_genes_highBias %>% 
  rename(gene=genes, datasets=n) %>% 
  ggplot(aes(x=reorder(gene, datasets), y=datasets)) +
  geom_bar(stat='identity', fill="#295D8A") +
  labs(title=paste("Distribution of", nrow(neg_genes_highBias) ,"neg. correlated genes in high-bias datasets"), 
      #  subtitle=paste0("Datasets selection (n=", length(dataids.list), "): bias>", bias.high.cutoff," and samples\u2265", samples.cutoff), 
       x="genes") +
  coord_flip() +
  theme_minimal(base_size = 20) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
ggsave(filename=path, plot=neg.cor.all.genes.highBias, device='pdf', width=12, height=10)

# figure 4 manuscript
figure <- ggarrange(
  pos.cor.genes.noBias,
  neg.cor.genes.noBias,
  labels = c("A", "B"),
  font.label = list(size = 20, color = "black", face = "bold", family = NULL),
  common.legend = TRUE, legend = "bottom",
  ncol = 2, nrow = 1)
path <- file.path(out.dir.path, "figure-4.pdf")
ggsave(filename=path, plot=figure, width=10, height=5)


# ==============================================================================
# GENE SET ENRICHMENT ANALYSIS
# ==============================================================================

# FUNCTION
# Computes gene set enrichment analysis (GSEA). Produces 2 tables and 1 plot
#
# gene_rank_table: table (col1:gene symbol; col2:rank)
# out.file.prefix: output file prefix
# out.dir: output directory
# top_n_to_plot: max number of enrichment (or gene sets / pathways) terms to plot
# table.cutoff: Adjusted p-value cutoff to limit output table
# plot.cutoff: Adjusted p-value cutoff to color enriched terms on plot
# genesets: local path to MSigDB genesets file (download at https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/)
gsea_analysis <- function(gene_rank_table,
                          out.file.prefix="mygsea",
                          out.dir=".",
                          top_n_to_plot=40,
                          table.cutoff=0.15,
                          plot.cutoff=0.15,
                          genesets) {

  # gene_rank_table = pos_genes
  # out.file.prefix="pos"
  # filename="pos"
  # out.dir=out.dir.path
  # top_n_to_plot=40
  # table.cutoff=0.5
  # plot.cutoff=0.15
  # genesets = msigdb.kegg.path

  path0 <- file.path(out.dir, paste0(out.file.prefix, ".gsea.sets.tsv"))
  path1 <- file.path(out.dir, paste0(out.file.prefix, ".gsea.plot.pdf"))
  path2 <- file.path(out.dir, paste0(out.file.prefix, ".gsea.genes.annot.tsv"))

  ranks <- deframe(gene_rank_table)
  pathways <- gmtPathways(genesets)
  fgseaRes <- fgsea(pathways=pathways, stats=ranks, nperm=1000) # Examples of fgsea at https://stephenturner.github.io/deseq-to-fgsea/
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
    filter(padj<table.cutoff) %>%
    mutate(significant=padj<plot.cutoff)

  if(dim(fgseaResTidy)[1]>0){

    write_tsv(fgseaResTidy, path0)

    gsea.plot <- ggplot(fgseaResTidy %>% top_n(top_n_to_plot, desc(abs(NES))), aes(reorder(pathway, NES), NES)) +
      geom_col(aes(fill=significant)) +
      coord_flip() +
      labs(x="Top pathways", y="Normalized Enrichment Score",
           title=paste("Gene Set Enrichment Analysis")) +
      theme_minimal(base_size = 20)
    ggsave(filename=path1, plot=gsea.plot, width=15.2, height=12.8)

    fgseaGenesAnnot <-  pathways %>%
      enframe("pathway", "SYMBOL") %>%
      unnest(cols = c(SYMBOL)) %>%
      inner_join(gene_rank_table, by=c("SYMBOL"="genes"))  %>%
      filter(pathway %in% fgseaResTidy$pathway)
    write_tsv(fgseaGenesAnnot, path2)
    } else{
      print("quality_correlated_genes.r: no correlated genes passed the cutoff")
      file.create(path0)
      file.create(path1)
      file.create(path2)
    }
}

#gsea_analysis(pos_genes, "pos.kegg", out.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.kegg.path)
#gsea_analysis(neg_genes, "neg.kegg", out.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.kegg.path)
#gsea_analysis(pos_genes, "pos.posi", out.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.posi.path)
#gsea_analysis(neg_genes, "neg.posi", out.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.posi.path)

# FUNCTION
# Computes gene set enrichment analysis (GSEA). If significant results are found:
# writes results in 2 tables and 1 plot on disk, and returns the plot as ggplot object. 
#
# genes: vector of strings as gene symbols (foreground set)
# universe: vector of strings as gene symbols (background set)
# out.file.prefix: string - file name prefix (not including file extension)
# out.dir: output directory
# top_n_to_plot: max number of enrichment (or gene sets / pathways) terms to plot
# enrich.cutoff: minumin percentage of pathway genes in the input gene list
# table.cutoff: Adjusted p-value cutoff to limit output table
# plot.cutoff: Adjusted p-value cutoff to color enriched terms on plot
# genesets: local path to MSigDB genesets file (download at https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/)
# create.default.files: if TRUE it creates empty files when the analysis has no significant results, if FALSE it creates non-empty files only for significant results

# RETURNS: ggplot object
ora_analysis <- function(plot.title="Gene Set Enrichment Analysis",
                         plot.subtitle="",
                         genes,
                         universe,
                          out.file.prefix="mygsea",
                          out.dir=".",
                          top_n_to_plot=40,
                          table.cutoff=0.15,
                          plot.cutoff=0.15,
                          genesets,
                          path.enrich.cutoff=0.01,
                          genes.enrich.cutoff=0.1,
                          create.default.files=FALSE) {

#   genes = pos_genes$genes[1:100]
#   universe = pos_genes$genes
#   out.file.prefix="pos.posi"
#   out.dir=out.dir.path
#   top_n_to_plot=40
#   table.cutoff=gsea.table.cutoff
#   plot.cutoff=gsea.plot.cutoff
#   genesets = msigdb.posi.path
#   path.enrich.cutoff=gsea.pathway.perc.cutoff
#   genes.enrich.cutoff=gsea.input.perc.cutoff
#   create.default.files=FALSE

  path0 <- file.path(out.dir, paste0(out.file.prefix, ".gsea.sets.tsv"))
  path1 <- file.path(out.dir, paste0(out.file.prefix, ".gsea.plot.pdf"))
  path2 <- file.path(out.dir, paste0(out.file.prefix, ".gsea.genes.annot.tsv"))
  pathways <- gmtPathways(genesets)
  pathways.universe <- unique(unlist(pathways))
  pathways.universe <- pathways.universe[ pathways.universe %in% universe ]

  pathways_tibble <- enframe(pathways) %>% 
    unnest_longer(value)
  pathways_tibble <- pathways_tibble[ pathways_tibble$value %in% universe, ]
  pathways_tibble <- pathways_tibble %>%
    group_by(name) %>%
    filter(n()>=5)
  pathways2 <- deframe(pathways_tibble %>% nest(data=c(value)))
  pathways2 <- lapply(pathways2, function(x) x$value )
  # print(paste0(genes))
  foraRes <- fora(pathways2, genes , universe)
  foraResTidy <- foraRes %>%
    as_tibble() %>%
    mutate(input_genes=length(genes)) %>%
    mutate(universe_genes=length(universe)) %>%
    filter(padj<table.cutoff) %>%
    filter((overlap/input_genes)>genes.enrich.cutoff) %>%
    filter((overlap/size)>path.enrich.cutoff) %>%
    mutate(significant=padj<plot.cutoff) %>%
    mutate(fold_change= round((overlap/input_genes) / (size/universe_genes), 2)) %>%
    arrange(desc(fold_change))
  # print(foraResTidy)
  gsea.plot <- ggplot() +
    labs(x="", y="Fold change",
         title=paste(plot.title),
         subtitle=paste(plot.subtitle)
    ) + theme_void()
  if(dim(foraResTidy)[1]>0){

    write_tsv(foraResTidy, path0)

    gsea.plot <- foraResTidy %>% 
      filter(padj<plot.cutoff) %>%
      slice_max(fold_change, n=top_n_to_plot) %>%
      # ggplot(aes(reorder(pathway, fold_change), fold_change)) +
      mutate(pathway=str_to_lower(str_wrap(str_replace_all(str_trunc(pathway, 60), "_", " "), 60))) %>% 
      ggplot(aes(reorder(pathway, fold_change), fold_change)) +
      geom_col(aes(fill=-log10(padj))) +
      coord_flip() +
      labs(x="", y="Fold change",
           title=paste(plot.title),
           subtitle=paste(plot.subtitle)
           ) +
      # theme_minimal(base_size = 10)
    theme_minimal()
    ggsave(filename=path1, plot=gsea.plot, width=8, height=8)

    fgseaGenesAnnot <-  pathways %>%
      enframe("pathway", "SYMBOL") %>%
      unnest(cols = c(SYMBOL)) %>%
      inner_join(data.frame(genes=genes), by=c("SYMBOL"="genes"))  %>%
      filter(pathway %in% foraResTidy$pathway)
    write_tsv(fgseaGenesAnnot, path2)
    } else{
      if(create.default.files==TRUE){
        print("quality_correlated_genes.r: no correlated genes passed the cutoff")
        file.create(path0)
        file.create(path1)
        file.create(path2)  
      }
    }
  return(gsea.plot)
}

# ora_analysis("Chromosomes enrichment", "Top pos. cor. genes | low-bias datasets", top_pos_genes, pos_genes$genes, "pos.posi", out.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.posi.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff, create.default.files=T)
# ora_analysis("Canonical pathways enrichment", "Top pos. cor. genes | low-bias datasets", top_pos_genes, pos_genes$genes, "pos.path", out.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.curated.cp.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff, create.default.files=T)
# ora_analysis("Chromosomes enrichment", "Top neg. cor. genes | low-bias datasets", top_neg_genes, neg_genes$genes, "neg.posi", out.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.posi.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff, create.default.files=T)
# ora_analysis("Canonical pathways enrichment", "Top neg. cor. genes | low-bias datasets", top_neg_genes, neg_genes$genes, "neg.path", out.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.curated.cp.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff, create.default.files=T)

# NO BIAS
top_pos_genes_noBias <- pos_genes_noBias$genes[1:gsea.input_genes]
top_neg_genes_noBias <- neg_genes_noBias$genes[1:gsea.input_genes]
print(top_pos_genes_noBias)
print(pos_genes_noBias)
gsea.pos.noBias.hall.plot <- ora_analysis("Hallmark pathways enrichment", "Low Quality Markers", top_pos_genes_noBias, pos_genes_noBias$genes, "fgsea.pos.hall.noBias", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.hallmark.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)
gsea.pos.noBias.posi.plot <- ora_analysis("Chromosomes enrichment", "Low Quality Markers", top_pos_genes_noBias, pos_genes_noBias$genes, "fgsea.pos.posi.noBias", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.posi.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff, create.default.files=F)
gsea.pos.noBias.cure.plot <- ora_analysis("Curated pathways enrichment", "Pathways of low-quality markers", top_pos_genes_noBias, pos_genes_noBias$genes, "fgsea.pos.cure.noBias", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.curated.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff, create.default.files=T)
gsea.pos.noBias.path.plot <- ora_analysis("Canonical pathways enrichment", "Low Quality Markers", top_pos_genes_noBias, pos_genes_noBias$genes, "fgsea.pos.path.noBias", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.curated.cp.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff, create.default.files=F)
gsea.pos.noBias.regu.plot <- ora_analysis("Regulatory target genes enrichment", "Regulators of low-quality markers", top_pos_genes_noBias, pos_genes_noBias$genes, "fgsea.pos.regu.noBias", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.regulatory.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)
gsea.pos.noBias.cell.plot <- ora_analysis("Cell type signature genes enrichment", "Low Quality Markers", top_pos_genes_noBias, pos_genes_noBias$genes, "fgsea.pos.cell.noBias", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.cells.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)
gsea.pos.noBias.gs2d.plot <- ora_analysis("Disease genes enrichment", "Low Quality Markers", top_pos_genes_noBias, pos_genes_noBias$genes, "fgsea.pos.gs2d.noBias", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.gs2d.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)

gsea.neg.noBias.hall.plot <- ora_analysis("Hallmark pathways enrichment", "High Quality Markers", top_neg_genes_noBias, neg_genes_noBias$genes, "fgsea.neg.hall.noBias", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.hallmark.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)
gsea.neg.noBias.posi.plot <- ora_analysis("Chromosomes enrichment", "High Quality Markers", top_neg_genes_noBias, neg_genes_noBias$genes, "fgsea.neg.posi.noBias", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.posi.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff, create.default.files=F)
gsea.neg.noBias.cure.plot <- ora_analysis("Curated pathways enrichment", "Pathways of high-quality markers", top_neg_genes_noBias, neg_genes_noBias$genes, "fgsea.neg.cure.noBias", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.curated.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff, create.default.files=T)
gsea.neg.noBias.path.plot <- ora_analysis("Canonical pathways enrichment", "High Quality Markers", top_neg_genes_noBias, neg_genes_noBias$genes, "fgsea.neg.path.noBias", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.curated.cp.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff, create.default.files=F)
gsea.neg.noBias.regu.plot <- ora_analysis("Regulatory target genes enrichment", "Regulators of high-quality markers", top_neg_genes_noBias, neg_genes_noBias$genes, "fgsea.neg.regu.noBias", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.regulatory.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)
gsea.neg.noBias.cell.plot <- ora_analysis("Cell type signature genes enrichment", "High Quality Markers", top_neg_genes_noBias, neg_genes_noBias$genes, "fgsea.neg.cell.noBias", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.cells.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)
gsea.neg.noBias.gs2d.plot <- ora_analysis("Disease genes enrichment", "High Quality Markers", top_neg_genes_noBias, neg_genes_noBias$genes, "fgsea.neg.gs2d.noBias", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.gs2d.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)


# LOW BIAS
top_pos_genes <- pos_genes$genes[1:gsea.input_genes]
top_neg_genes <- neg_genes$genes[1:gsea.input_genes]
gsea.pos.lowBias.hall.plot <- ora_analysis("Hallmark pathways enrichment", "Low Quality Markers | low-bias datasets", top_pos_genes, pos_genes$genes, "fgsea.pos.hall", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.hallmark.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)
gsea.pos.lowBias.posi.plot <- ora_analysis("Chromosomes enrichment", "Low Quality Markers | low-bias datasets", top_pos_genes, pos_genes$genes, "fgsea.pos.posi", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.posi.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff, create.default.files=F)
gsea.pos.lowBias.cure.plot <- ora_analysis("Curated pathways enrichment", "Low Quality Markers | low-bias datasets", top_pos_genes, pos_genes$genes, "fgsea.pos.cure", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.curated.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff, create.default.files=T)
gsea.pos.lowBias.path.plot <- ora_analysis("Canonical pathways enrichment", "Low Quality Markers | low-bias datasets", top_pos_genes, pos_genes$genes, "fgsea.pos.path", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.curated.cp.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff, create.default.files=F)
gsea.pos.lowBias.regu.plot <- ora_analysis("Regulatory target genes enrichment", "Low Quality Markers | low-bias datasets", top_pos_genes, pos_genes$genes, "fgsea.pos.regu", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.regulatory.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)
gsea.pos.lowBias.cell.plot <- ora_analysis("Cell type signature genes enrichment", "Low Quality Markers | low-bias datasets", top_pos_genes, pos_genes$genes, "fgsea.pos.cell", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.cells.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)
gsea.pos.lowBias.gs2d.plot <- ora_analysis("Disease genes enrichment", "Low Quality Markers | low-bias datasets", top_pos_genes, pos_genes$genes, "fgsea.pos.gs2d", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.gs2d.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)

gsea.neg.lowBias.hall.plot <- ora_analysis("Hallmark pathways enrichment", "High Quality Markers | low-bias datasets", top_neg_genes, neg_genes$genes, "fgsea.neg.hall", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.hallmark.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)
gsea.neg.lowBias.posi.plot <- ora_analysis("Chromosomes enrichment", "High Quality Markers | low-bias datasets", top_neg_genes, neg_genes$genes, "fgsea.neg.posi", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.posi.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff, create.default.files=F)
gsea.neg.lowBias.cure.plot <- ora_analysis("Curated pathways enrichment", "High Quality Markers | low-bias datasets", top_neg_genes, neg_genes$genes, "fgsea.neg.cure", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.curated.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff, create.default.files=T)
gsea.neg.lowBias.path.plot <- ora_analysis("Canonical pathways enrichment", "High Quality Markers | low-bias datasets", top_neg_genes, neg_genes$genes, "fgsea.neg.path", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.curated.cp.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff, create.default.files=F)
gsea.neg.lowBias.regu.plot <- ora_analysis("Regulatory target genes enrichment", "High Quality Markers | low-bias datasets", top_neg_genes, neg_genes$genes, "fgsea.neg.regu", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.regulatory.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)
gsea.neg.lowBias.cell.plot <- ora_analysis("Cell type signature genes enrichment", "High Quality Markers | low-bias datasets", top_neg_genes, neg_genes$genes, "fgsea.neg.cell", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.cells.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)
gsea.neg.lowBias.gs2d.plot <- ora_analysis("Disease genes enrichment", "High Quality Markers | low-bias datasets", top_neg_genes, neg_genes$genes, "fgsea.neg.gs2d", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.gs2d.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)

# HIGH BIAS
top_pos_genes_highBias <- pos_genes_highBias$genes[1:gsea.input_genes]
top_neg_genes_highBias <- neg_genes_highBias$genes[1:gsea.input_genes]
gsea.pos.highBias.hall.plot <- ora_analysis("Hallmark pathways enrichment", "Top pos. cor. genes | high-bias datasets", top_pos_genes_highBias, pos_genes_highBias$genes, "fgsea.pos.hall.highBias", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.hallmark.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)
gsea.pos.highBias.posi.plot <- ora_analysis("Chromosomes enrichment", "Top pos. cor. genes | high-bias datasets", top_pos_genes_highBias, pos_genes_highBias$genes, "fgsea.pos.posi.highBias", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.posi.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)
gsea.pos.highBias.cure.plot <- ora_analysis("Curated pathways enrichment", "Top pos. cor. genes | high-bias datasets", top_pos_genes_highBias, pos_genes_highBias$genes, "fgsea.pos.cure.highBias", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.curated.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff, create.default.files=T)
gsea.pos.highBias.path.plot <- ora_analysis("Canonical pathways enrichment", "Top pos. cor. genes | high-bias datasets", top_pos_genes_highBias, pos_genes_highBias$genes, "fgsea.pos.path.highBias", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.curated.cp.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)
gsea.pos.highBias.regu.plot <- ora_analysis("Regulatory target genes enrichment", "Top pos. cor. genes | high-bias datasets", top_pos_genes_highBias, pos_genes_highBias$genes, "fgsea.pos.regu.highBias", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.regulatory.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)
gsea.pos.highBias.cell.plot <- ora_analysis("Cell type signature genes enrichment", "Top pos. cor. genes | high-bias datasets", top_pos_genes_highBias, pos_genes_highBias$genes, "fgsea.pos.cell.highBias", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.cells.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)
gsea.pos.highBias.gs2d.plot <- ora_analysis("Disease genes enrichment", "Top pos. cor. genes | high-bias datasets", top_pos_genes_highBias, pos_genes_highBias$genes, "fgsea.pos.gs2d.highBias", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.gs2d.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)

gsea.neg.highBias.hall.plot <- ora_analysis("Hallmark pathways enrichment", "Top neg. cor. genes | high-bias datasets", top_neg_genes_highBias, neg_genes_highBias$genes, "fgsea.neg.hall.highBias", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.hallmark.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)
gsea.neg.highBias.posi.plot <- ora_analysis("Chromosomes enrichment", "Top neg. cor. genes | high-bias datasets", top_neg_genes_highBias, neg_genes_highBias$genes, "fgsea.neg.posi.highBias", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.posi.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)
gsea.neg.highBias.cure.plot <- ora_analysis("Curated pathways enrichment", "Top neg. cor. genes | high-bias datasets", top_neg_genes_highBias, neg_genes_highBias$genes, "fgsea.neg.cure.highBias", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.curated.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff, create.default.files=T)
gsea.neg.highBias.path.plot <- ora_analysis("Canonical pathways enrichment", "Top neg. cor. genes | high-bias datasets", top_neg_genes_highBias, neg_genes_highBias$genes, "fgsea.neg.path.highBias", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.curated.cp.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)
gsea.neg.highBias.regu.plot <- ora_analysis("Regulatory target genes enrichment", "Top neg. cor. genes | high-bias datasets", top_neg_genes_highBias, neg_genes_highBias$genes, "fgsea.neg.regu.highBias", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.regulatory.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)
gsea.neg.highBias.cell.plot <- ora_analysis("Cell type signature genes enrichment", "Top neg. cor. genes | high-bias datasets", top_neg_genes_highBias, neg_genes_highBias$genes, "fgsea.neg.cell.highBias", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.cells.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)
gsea.neg.highBias.gs2d.plot <- ora_analysis("Disease genes enrichment", "Top neg. cor. genes | high-bias datasets", top_neg_genes_highBias, neg_genes_highBias$genes, "fgsea.neg.gs2d.highBias", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.gs2d.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)


figure.width <- 14
figure.height <- 16

figure <- ggarrange(
  gsea.pos.noBias.hall.plot,
  gsea.neg.noBias.hall.plot,
  gsea.pos.lowBias.hall.plot,
  gsea.neg.lowBias.hall.plot,
  gsea.pos.highBias.hall.plot,
  gsea.neg.highBias.hall.plot,
  # labels = c("A", "B", "C", "D", "E", "F"),
  common.legend = TRUE, legend = "bottom",
  ncol = 2, nrow = 3)
path <- file.path(out.dir.path, "fgsea.hall.pdf")
ggsave(filename=path, plot=figure, width=figure.width, height=figure.height)

figure <- ggarrange(
  gsea.pos.noBias.posi.plot,
  gsea.neg.noBias.posi.plot,
  gsea.pos.lowBias.posi.plot,
  gsea.neg.lowBias.posi.plot,
  gsea.pos.highBias.posi.plot,
  gsea.neg.highBias.posi.plot,
  # labels = c("A", "B", "C", "D", "E", "F"),
  common.legend = TRUE, legend = "bottom",
  ncol = 2, nrow = 3)
path <- file.path(out.dir.path, "fgsea.posi.pdf")
ggsave(filename=path, plot=figure, width=figure.width, height=figure.height)

figure <- ggarrange(
  gsea.pos.noBias.cure.plot,
  gsea.neg.noBias.cure.plot,
  gsea.pos.lowBias.cure.plot,
  gsea.neg.lowBias.cure.plot,
  gsea.pos.highBias.cure.plot,
  gsea.neg.highBias.cure.plot,
  # labels = c("A", "B", "C", "D", "E", "F"),
  common.legend = TRUE, legend = "bottom",
  ncol = 2, nrow = 3)
path <- file.path(out.dir.path, "fgsea.cure.pdf")
ggsave(filename=path, plot=figure, width=figure.width, height=figure.height)

figure <- ggarrange(
  gsea.pos.noBias.path.plot,
  gsea.neg.noBias.path.plot,
  gsea.pos.lowBias.path.plot,
  gsea.neg.lowBias.path.plot,
  gsea.pos.highBias.path.plot,
  gsea.neg.highBias.path.plot,
  # labels = c("A", "B", "C", "D", "E", "F"),
  common.legend = TRUE, legend = "bottom",
  ncol = 2, nrow = 3)
path <- file.path(out.dir.path, "fgsea.path.pdf")
ggsave(filename=path, plot=figure, width=figure.width, height=figure.height)

figure <- ggarrange(
  gsea.pos.noBias.regu.plot,
  gsea.neg.noBias.regu.plot,
  gsea.pos.lowBias.regu.plot,
  gsea.neg.lowBias.regu.plot,
  gsea.pos.highBias.regu.plot,
  gsea.neg.highBias.regu.plot,
  # labels = c("A", "B", "C", "D", "E", "F"),
  common.legend = TRUE, legend = "bottom",
  ncol = 2, nrow = 3)
path <- file.path(out.dir.path, "fgsea.regu.pdf")
ggsave(filename=path, plot=figure, width=figure.width, height=figure.height)

figure <- ggarrange(
  gsea.pos.noBias.cell.plot,
  gsea.neg.noBias.cell.plot,
  gsea.pos.lowBias.cell.plot,
  gsea.neg.lowBias.cell.plot,
  gsea.pos.highBias.cell.plot,
  gsea.neg.highBias.cell.plot,
  # labels = c("A", "B", "C", "D", "E", "F"),
  common.legend = TRUE, legend = "bottom",
  ncol = 2, nrow = 3)
path <- file.path(out.dir.path, "fgsea.cell.pdf")
ggsave(filename=path, plot=figure, width=figure.width, height=figure.height)

figure <- ggarrange(
  gsea.pos.noBias.gs2d.plot,
  gsea.neg.noBias.gs2d.plot,
  gsea.pos.lowBias.gs2d.plot,
  gsea.neg.lowBias.gs2d.plot,
  gsea.pos.highBias.gs2d.plot,
  gsea.neg.highBias.gs2d.plot,
  # labels = c("A", "B", "C", "D", "E", "F"),
  common.legend = TRUE, legend = "bottom",
  ncol = 2, nrow = 3)
path <- file.path(out.dir.path, "fgsea.gs2d.pdf")
ggsave(filename=path, plot=figure, width=figure.width, height=figure.height)

# figure 5 of the manuscript
figure <- ggarrange(
  gsea.pos.noBias.regu.plot + ggtitle(NULL),
  gsea.neg.noBias.regu.plot + ggtitle(NULL),
  gsea.pos.noBias.cure.plot + ggtitle(NULL),
  gsea.neg.noBias.cure.plot + ggtitle(NULL),
  labels = c("A", "B", "C", "D"),
  font.label = list(size = 20, color = "black", face = "bold", family = NULL),
  common.legend = TRUE, legend = "bottom",
  ncol = 2, nrow = 2)
path <- file.path(out.dir.path, "figure-5.noBias.pdf")
ggsave(filename=path, plot=figure, width=14, height=14)

figure <- ggarrange(
  gsea.pos.lowBias.regu.plot + ggtitle(NULL),
  gsea.neg.lowBias.regu.plot + ggtitle(NULL),
  gsea.pos.lowBias.cure.plot + ggtitle(NULL),
  gsea.neg.lowBias.cure.plot + ggtitle(NULL),
  labels = c("A", "B", "C", "D"),
  font.label = list(size = 20, color = "black", face = "bold", family = NULL),
  common.legend = TRUE, legend = "bottom",
  ncol = 2, nrow = 2)
path <- file.path(out.dir.path, "figure-5.lowBias.pdf")
ggsave(filename=path, plot=figure, width=14, height=14)


# ==============================================================================
# LINEAR REGRESSIONS
# ==============================================================================
# d <- datasets.stats %>% filter(selected==TRUE)
d <- ""
datasets.stats <- datasets.stats %>% mutate(`Exceeds Threshold`=ifelse(p_group_cor >= 0.3, T, F)) 
# Quality vs genes' FDR
# ------------------------------------------------------------------------------
path <- file.path(out.dir.path, "cor.q.fdr.pdf")
# plot <- datasets.stats %>% 
#   select(dataset, cor_q_fdr_allgenes, cor_q_fdr_disgenes, cor_q_fdr_difgenes, selected, deg) %>%
#   gather(cor_q_fdr_allgenes, cor_q_fdr_disgenes, cor_q_fdr_difgenes, key = "feature", value = "value") %>%
#   ggplot(aes(x=reorder(dataset, 100*selected + value), y=value, fill=selected)) +
#   geom_bar(stat = "identity") +
#   geom_text(aes(label=round(value,2))) +
#   coord_flip() +
#   facet_wrap( ~ feature) +
#   theme_minimal() +
#   ggsave(path, width=15.2, height=12.8)
# 
plot <- datasets.stats %>% 
  select(dataset, selected, n_degs, `Exceeds Threshold`,
         cor_q_fdr_allgenes, cor_q_fdr_disgenes, cor_q_fdr_difgenes,
         cor_q_fdr_pos_allgenes, cor_q_fdr_pos_disgenes, cor_q_fdr_pos_difgenes,
         cor_q_fdr_neg_allgenes, cor_q_fdr_neg_disgenes, cor_q_fdr_neg_difgenes
         ) %>%
  rename(ALL_ALL_Rqe_FDRde=cor_q_fdr_allgenes, 
         DIS_ALL_Rqe_FDRtm=cor_q_fdr_disgenes,
         DIF_ALL_Rqe_FDRde=cor_q_fdr_difgenes,
         ALL_POS_Rqe_FDRde=cor_q_fdr_pos_allgenes, 
         DIS_POS_Rqe_FDRtm=cor_q_fdr_pos_disgenes, 
         DIF_POS_Rqe_FDRde=cor_q_fdr_pos_difgenes,
         ALL_NEG_Rqe_FDRde=cor_q_fdr_neg_allgenes, 
         DIS_NEG_Rqe_FDRtm=cor_q_fdr_neg_disgenes, 
         DIF_NEG_Rqe_FDRde=cor_q_fdr_neg_difgenes) %>% 
  # gather(cor_q_fdr_allgenes, cor_q_fdr_disgenes, cor_q_fdr_difgenes,
  #      cor_q_fdr_pos_allgenes, cor_q_fdr_pos_disgenes, cor_q_fdr_pos_difgenes,
  #      cor_q_fdr_neg_allgenes, cor_q_fdr_neg_disgenes, cor_q_fdr_neg_difgenes,
  #      key = "feature", value = "value") %>%
  gather(ALL_ALL_Rqe_FDRde, 
         DIS_ALL_Rqe_FDRtm, 
         DIF_ALL_Rqe_FDRde, 
         ALL_POS_Rqe_FDRde,
         DIS_POS_Rqe_FDRtm, 
         DIF_POS_Rqe_FDRde, 
         ALL_NEG_Rqe_FDRde,
         DIS_NEG_Rqe_FDRtm, 
         DIF_NEG_Rqe_FDRde, key = "feature", value = "value") %>%
  ggplot(aes(x=reorder(dataset, -10*selected + value), y=value, fill=`Exceeds Threshold`)) +
  geom_bar(stat = "identity") +
  # geom_text(aes(label=round(value,2)), size=4) +
  coord_flip() +
  facet_wrap( ~ feature, labeller = label_parsed) +
  labs(title=paste("Dataset Quality vs FDR"), 
       x="Datasets", 
       y="Linear regression coefficients") + 
  theme_minimal(base_size = 18)
ggsave(filename=path, plot=plot, width=15.2, height=20)


# ==============================================================================
# GLOBAL STATISTICS
# ==============================================================================
# d <- datasets.stats

# Quality vs Data Size plot
# ------------------------------------------------------------------------------
path <- file.path(out.dir.path, "cor.size.pdf")
plot <- datasets.stats %>% 
  select(dataset, p_bases_cor, p_bases_cor_ctrl, p_bases_cor_treat, p_bytes_cor, 
         p_bytes_cor_ctrl, p_bytes_cor_treat, n_degs, selected, `Exceeds Threshold`) %>% 
  gather(p_bases_cor, p_bases_cor_ctrl, p_bases_cor_treat, p_bytes_cor, 
         p_bytes_cor_ctrl, p_bytes_cor_treat, key="Feature", value="value") %>% 
  ggplot(aes(x=reorder(dataset, 100*selected + value), y=value, fill=`Exceeds Threshold`)) +
  geom_bar(stat='identity') +
  ylim(c(-1.1, 1.1)) +
  geom_text(aes(label=n_degs), size=5) +
  facet_wrap(~Feature) +
  coord_flip() +
  labs(title=paste("Dataset Quality vs Size"), 
       x="Datasets", 
       y="Correlation coefficient (Low-Quality Probability vs Size in Bases or Bytes)") + 
  theme_minimal(base_size = 20)
ggsave(filename=path, plot=plot, width=15.2, height=30)


# Dataset Selection Plot
# ------------------------------------------------------------------------------
path <- file.path(out.dir.path, "figure-1-cor.group.pdf")
bias_deg.all.cor <- round(cor(datasets.stats$p_group_cor, datasets.stats$n_degs, use="pairwise.complete.obs"), 2)
bias_deg.selection.cor <- round(cor(datasets.stats[datasets.stats$selected==TRUE,]$p_group_cor, datasets.stats[datasets.stats$selected==TRUE,]$n_degs), 2)

plot.cor <- datasets.stats %>%
  select(dataset, p_group_cor, `Exceeds Threshold`, n_degs) %>% 
  mutate(`Quality Imbalance` = ifelse(
    p_group_cor >= bias.high.cutoff, "high", ifelse(
      p_group_cor <= bias.low.cutoff, "low", "medium"
  ))) %>%
  ggplot(aes(x=reorder(dataset, -p_group_cor), y=p_group_cor, fill=`Quality Imbalance`)) +
  geom_bar(stat='identity') +
  ylim(c(0, 1.1)) +
  geom_text(aes(label=n_degs), hjust=-0.25, size=5) +
  # geom_hline(yintercept = dataset_p_group_cor_cutoff, col="gray") +
  coord_flip() +
  labs(#title=paste("Datasets Selection"), 
      #  subtitle = paste0("dataset bias<", dataset_p_group_cor_cutoff, 
      #                   " &  diff. expressed genes (degs)\u2265", dataset_n_deg_cutoff,
      #                   "\nbias-degs cor=", bias_deg.all.cor,
      #                   "| bias-degs cor (selected datasets)=", bias_deg.selection.cor),
       x="Datasets", 
       y="Quality Imbalance") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "bottomleft", legend.direction = "horizontal") +
  scale_fill_manual(values=list("high"= "#A41720",
                          "low" = "#295D8A",
                          "medium" = "#4E4459"))

ggsave(filename=path, plot=plot.cor, width=13, height=17)

# options(ggplot2.discrete.fill= c("#295D8A", "#A41720", "#4E4459"))

# plot statistics
meta_inp <- datasets.table %>% 
  filter(!grepl("sub", GEO_Series)) %>%
  filter(select == 1) %>%
  mutate(mesh_terms_simp = str_split(mesh_terms, "\\|", simplify = TRUE)[, 1])  


sapply(strsplit(meta_inp$mesh_terms, "\\|"), function(x) ifelse(length(x) > 1, x[1], x))

path <- file.path(out.dir.path, "test_bar_samples_pairing.pdf")
bar_samples_pairing <- meta_inp %>%
  mutate(SamplesPairing=ifelse(SamplesPairing==1, "TRUE", "FALSE")) %>%
  ggplot(aes(x=SamplesPairing)) +
  geom_bar(fill="#4E4459") +
  labs(x="Samples Pairing", y=NULL) +
    theme_minimal(base_size = 16) 

ggsave(filename=path, plot=bar_samples_pairing, width=2, height=4)

path <- file.path(out.dir.path, "test_bar_single-paired-end.pdf")
bar_reads <- meta_inp %>%
  ggplot(aes(x=LibraryLayout)) +
  geom_bar(fill="#4E4459") +
  labs(x="Library Layout", y=NULL) +
    theme_minimal(base_size = 16)

ggsave(filename=path, plot=bar_reads, width=2, height=4)


path <- file.path(out.dir.path, "test_dist-IF.pdf")
if_dens <- meta_inp %>%
  ggplot(aes(x=`2-YEAR-IF2020`)) +
  # geom_histogram() +
  geom_density(color="#4E4459") +
  labs(x="2 Year Impact Factor 2020", y=NULL) +
    theme_minimal(base_size = 16) 

ggsave(filename=path, plot=if_dens, width=2, height=4)

path <- file.path(out.dir.path, "test_bar_diseases.pdf")
# bar_disease <- meta_inp %>%
#   ggplot(aes(x="", y=mesh_terms_simp, fill=mesh_terms_simp)) +
#   geom_bar(width=1, stat="identity") +
#   #coord_polar("y", start=0) +
#   guides(fill="none")
bar_disease <- meta_inp %>%
  group_by(mesh_terms_simp) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  arrange(count) %>%
  ggplot(aes(x = "", y = count,  label = mesh_terms_simp)) +
  geom_bar(stat = "identity",fill = "#4E4459",  color = "#295D8A") +
  geom_text(aes(label = mesh_terms_simp), color='lightgrey', position = position_stack(vjust = 0.5)) +
  labs(x = "Distribution of diseases", y = NULL, fill = "Mesh terms") + 
    theme_minimal(base_size = 16) +
  guides(fill="none") 

ggsave(filename=path, plot=bar_disease, width=4, height=8)

# dist DEGs
path <- file.path(out.dir.path, "test_hist_degs.pdf")
degs <- datasets.stats  %>%
  ggplot(aes(x=n_degs)) +
  geom_histogram(fill="#4E4459") +
  labs(x="Number of DEGs", y=NULL) +
  theme_minimal(base_size = 16) 

ggsave(filename=path, plot=degs, width=4, height=8)


# figure-1-composit
figure <- ggarrange(plot.cor, bar_disease,
          ggarrange(bar_reads, bar_samples_pairing, degs, if_dens, ncol = 1, nrow = 4, labels =c("C", "D", "E", "F")), 
          labels = c("A", "B"),
          widths = c(2, 1.5, 1.5),
          heights = c(100, 1, 1),
          ncol = 3, nrow = 3,
          common.legend = TRUE, legend = "none")

path <- file.path(out.dir.path, "figure-1-composit.pdf")
ggsave(filename=path, plot=figure, width=14, height=14)

# path <- file.path(out.dir.path, "cor.group.pval.pdf")
# bias_deg.all.cor <- round(cor(datasets.stats$p_group_cor, datasets.stats$n_degs, use="pairwise.complete.obs"), 2)
# bias_deg.selection.cor <- round(cor(datasets.stats[datasets.stats$selected==TRUE,]$p_group_cor, datasets.stats[datasets.stats$selected==TRUE,]$n_degs), 2)

# plot <- datasets.stats %>% mutate(`Exceeds Threshold`=ifelse(p_group_cor >= 0.3, T, F)) %>%
#   select(dataset, p_group_cor, `Exceeds Threshold`, p_group_cor_pval, n_degs) %>% 
#   ggplot(aes(x=reorder(dataset, -p_group_cor), y=p_group_cor, fill=`Exceeds Threshold`)) +
#   geom_bar(stat='identity') +
#   ylim(c(0, 1.1)) +
#   geom_text(aes(label=round(p_group_cor_pval, 3)), hjust=-0.25, size=7) +
#   # geom_hline(yintercept = dataset_p_group_cor_cutoff, col="gray") +
#   coord_flip() +
#   labs(#title=paste("Datasets Selection"), 
#       #  subtitle = paste0("dataset bias<", dataset_p_group_cor_cutoff, 
#       #                   " &  diff. expressed genes (degs)\u2265", dataset_n_deg_cutoff,
#       #                   "\nbias-degs cor=", bias_deg.all.cor,
#       #                   "| bias-degs cor (selected datasets)=", bias_deg.selection.cor),
#        x="Datasets", 
#        y="Quality Imbalance") +
#   theme_minimal(base_size = 20)
# ggsave(filename=path, plot=plot, width=13, height=20)

# DEGs vs n.samples
# ------------------------------------------------------------------------------
# swithc colors to mach "imbalance"
options(ggplot2.discrete.colour= c("#A41720","#295D8A", "#4E4459"))

path <- file.path(out.dir.path, "deg.samples.pdf")
my.formula <- y ~ x

deg.samples <-  datasets.stats %>% 
  filter(n_degs>=dataset_n_deg_cutoff) %>%
  filter(n_samples>=samples.cutoff) %>%
  select(dataset, p_group_cor, selected, n_degs, n_samples) %>% 
  rename(bias=p_group_cor) %>% 
  mutate(Imbalance=as.factor(if_else(bias>bias.high.cutoff, "true", "false"))) %>% 
  ggplot(aes(y=n_degs, x=n_samples, size=bias, color=Imbalance)) +
  geom_point() +
  geom_smooth(method="glm", formula=my.formula <- y ~ x) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               size = rel(6),
               parse = TRUE) +  
  labs(title=paste("FDR and fold change"), 
      #  subtitle = paste0(
      #    "Diff. genes defined by FDR and fold change",
      #   #  "\nDatasets selection: diff genes\u2265", dataset_n_deg_cutoff,
      #   #  " and samples\u2265", samples.cutoff,
      #    "\nQuality Imbalance: low \u2264 ", bias.high.cutoff, "< high imbalance"
      #  ),
       y="Differential genes", 
       x="Samples",
       size="QI Index",
       color="High QI") +
  theme_minimal(base_size = 20) +
  scale_color_manual(values=list("false"="#295D8A", "true"="#A41720"))

ggsave(filename=path, plot=deg.samples, width=12, height=10)


# DEGs vs n.samples  (SamplesPairing only)
# ------------------------------------------------------------------------------
my.formula <- y ~ x
path <- file.path(out.dir.path, "deg.paired.samples.pdf")

deg.paired.samples <-  datasets.stats %>% 
  filter(n_degs>=dataset_n_deg_cutoff) %>%
  filter(n_samples>=samples.cutoff) %>%
  filter(SamplesPairing==1) %>%
  select(dataset, p_group_cor, selected, n_degs, n_samples) %>% 
  rename(bias=p_group_cor) %>% 
  mutate(Imbalance=as.factor(if_else(bias>bias.high.cutoff, "true", "false"))) %>% 
  ggplot(aes(y=n_degs, x=n_samples, size=bias, color=Imbalance)) +
  geom_point() +
  geom_smooth(method="glm", formula=my.formula <- y ~ x) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               size = rel(6),
               parse = TRUE) +  
  labs(title=paste("FDR and fold change \nPaired Samples"), 
      #  subtitle = paste0(
      #    "Diff. genes defined by FDR and fold change",
      #   #  "\nDatasets selection: diff genes\u2265", dataset_n_deg_cutoff,
      #   #  " and samples\u2265", samples.cutoff,
      #    "\nQuality Imbalance: low \u2264 ", bias.high.cutoff, "< high imbalance"
      #  ),
       y="Differential genes", 
       x="Samples",
       size="QI Index",
       color="High QI") +
  theme_minimal(base_size = 20)+
  scale_color_manual(values=list("false"="#295D8A", "true"="#A41720"))

ggsave(filename=path, plot=deg.paired.samples, width=12, height=10)


# DEGs vs n.samples  (NO SamplesPairing only)
# ------------------------------------------------------------------------------
my.formula <- y ~ x
path <- file.path(out.dir.path, "deg.no.paired.samples.pdf")

deg.no.paired.samples <-  datasets.stats %>% 
  filter(n_degs>=dataset_n_deg_cutoff) %>%
  filter(n_samples>=samples.cutoff) %>%
  filter(SamplesPairing==0) %>%
  select(dataset, p_group_cor, selected, n_degs, n_samples) %>% 
  rename(bias=p_group_cor) %>% 
  mutate(Imbalance=as.factor(if_else(bias>bias.high.cutoff, "true", "false"))) %>% 
  ggplot(aes(y=n_degs, x=n_samples, size=bias, color=Imbalance)) +
  geom_point() +
  geom_smooth(method="glm", formula=my.formula <- y ~ x) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               size = rel(6),
               parse = TRUE) +  
  labs(title=paste("FDR and fold change \nNo Paired Samples"), 
      #  subtitle = paste0(
      #    "Diff. genes defined by FDR and fold change",
      #   #  "\nDatasets selection: diff genes\u2265", dataset_n_deg_cutoff,
      #   #  " and samples\u2265", samples.cutoff,
      #    "\nQuality Imbalance: low \u2264 ", bias.high.cutoff, "< high imbalance"
      #  ),
       y="Differential genes", 
       x="Samples",
       size="QI Index",
       color="High QI") +
  theme_minimal(base_size = 20)+
  scale_color_manual(values=list("false"="#295D8A", "true"="#A41720"))

ggsave(filename=path, plot=deg.no.paired.samples, width=12, height=10)


# DEGs vs n.samples (FDR only)
# ------------------------------------------------------------------------------
my.formula <- y ~ x
path <- file.path(out.dir.path, "deg.fdrOnly.samples.pdf")
# bias_deg.all.cor <- round(cor(datasets.stats$p_group_cor, datasets.stats$deg), 2)
# bias_deg.selection.cor <- round(cor(d$p_group_cor, d$deg), 2)
# sample_deg.all.cor <- round(cor(datasets.stats$n_samples, datasets.stats$deg), 2)
# sample_deg.selection.cor <- round(cor(d$n_samples, d$deg), 2)
# bias.medium <- 0.15
# bias.high.cutoff <- 0.20
# samples.cutoff <- 10

deg.fdrOnly.samples <-  datasets.stats %>% 
  filter(n_degs>=dataset_n_deg_cutoff) %>%
  filter(n_samples>=samples.cutoff) %>%
  select(dataset, p_group_cor, selected, n_degs_fdr_only, n_samples, SamplesPairing) %>% 
  rename(bias=p_group_cor) %>% 
  mutate(Imbalance=as.factor(if_else(bias>bias.high.cutoff, "true", "false"))) %>% 
#  mutate(Imbalance=if_else(bias>bias.high.cutoff, "true", "false")) %>% 
#  mutate(SamplesPairing=if_else(SamplesPairing==1, "paired samples", "non-paired samples")) %>% 
#  mutate(Imbalance=as.factor(paste(SamplesPairing, data, sep=" | "))) %>% 
  ggplot(aes(y=n_degs_fdr_only, x=n_samples, size=bias, color=Imbalance)) +
  geom_point() +
  geom_smooth(method="glm", formula=my.formula <- y ~ x) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               size = rel(6),
               parse = TRUE) +  
  labs(title=paste("FDR only"), 
      #  subtitle = paste0(
      #    #"Max correlation coeff.:", dataset_p_group_cor_cutoff, 
      #    "Diff. genes defined by FDR only",
      #   #  "\nDatasets selection: diff genes\u2265", dataset_n_deg_cutoff,
      #   #  " and samples\u2265", samples.cutoff,
      #     "\nQuality Imbalance: low \u2264 ", bias.high.cutoff, "< high imbalance"
      #    # ,
      #    #                        "\nsamples-genes cor: ", sample_deg.all.cor,
      #    # "samples-genes cor: ", sample_deg.all.cor,
      #    # "| samples-genes cor (selected datasets): ", sample_deg.selection.cor
      #  ),
       y="Differential genes", 
       x="Samples",
       size="QI Index",
       color="High QI") +
  theme_minimal(base_size = 20)+
  scale_color_manual(values=list("false"="#295D8A", "true"="#A41720"))

ggsave(filename=path, plot=deg.fdrOnly.samples, width=12, height=10)


# DEGs vs n.samples (FDR only + SamplesPairing only)
# ------------------------------------------------------------------------------
my.formula <- y ~ x
path <- file.path(out.dir.path, "deg.fdrOnly.paired.samples.pdf")
# bias_deg.all.cor <- round(cor(datasets.stats$p_group_cor, datasets.stats$deg), 2)
# bias_deg.selection.cor <- round(cor(d$p_group_cor, d$deg), 2)
# sample_deg.all.cor <- round(cor(datasets.stats$n_samples, datasets.stats$deg), 2)
# sample_deg.selection.cor <- round(cor(d$n_samples, d$deg), 2)
# bias.medium <- 0.15
# bias.high.cutoff <- 0.20
# samples.cutoff <- 10

deg.fdrOnly.paired.samples <-  datasets.stats %>% 
  filter(n_degs>=dataset_n_deg_cutoff) %>%
  filter(n_samples>=samples.cutoff) %>%
  filter(SamplesPairing==1) %>%
  select(dataset, p_group_cor, selected, n_degs_fdr_only, n_samples, SamplesPairing) %>% 
  rename(bias=p_group_cor) %>% 
  mutate(Imbalance=as.factor(if_else(bias>bias.high.cutoff, "true", "false"))) %>% 
  mutate(SamplesPairing=as.factor(SamplesPairing)) %>% 
  ggplot(aes(y=n_degs_fdr_only, x=n_samples, size=bias, color=Imbalance)) +
  geom_point() +
  geom_smooth(method="glm", formula=my.formula <- y ~ x) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               size = rel(6),
               parse = TRUE) +  
  labs(title=paste("FDR only \nPaired Samples"), 
      #  subtitle = paste0(
      #    #"Max correlation coeff.:", dataset_p_group_cor_cutoff, 
      #    "Diff. genes defined by FDR only",
      #   #  "\nDatasets selection: diff genes\u2265", dataset_n_deg_cutoff,
      #   #  " and samples\u2265", samples.cutoff,
      #     "\nQuality Imbalance: low \u2264 ", bias.high.cutoff, "< high imbalance"
      #    # ,
      #    #                        "\nsamples-genes cor: ", sample_deg.all.cor,
      #    # "samples-genes cor: ", sample_deg.all.cor,
      #    # "| samples-genes cor (selected datasets): ", sample_deg.selection.cor
      #  ),
       y="Differential genes", 
       x="Samples",
       size="QI Index",
       color="High QI") +
  theme_minimal(base_size = 20)+
  scale_color_manual(values=list("false"="#295D8A", "true"="#A41720"))

ggsave(filename=path, plot=deg.fdrOnly.paired.samples, width=12, height=10)


# DEGs vs n.samples (FDR only + NO SamplesPairing only)
# ------------------------------------------------------------------------------
my.formula <- y ~ x
path <- file.path(out.dir.path, "deg.fdrOnly.no.paired.samples.pdf")
# bias_deg.all.cor <- round(cor(datasets.stats$p_group_cor, datasets.stats$deg), 2)
# bias_deg.selection.cor <- round(cor(d$p_group_cor, d$deg), 2)
# sample_deg.all.cor <- round(cor(datasets.stats$n_samples, datasets.stats$deg), 2)
# sample_deg.selection.cor <- round(cor(d$n_samples, d$deg), 2)
# bias.medium <- 0.15
# bias.high.cutoff <- 0.20
# samples.cutoff <- 10

deg.fdrOnly.no.paired.samples <-  datasets.stats %>% 
  filter(n_degs>=dataset_n_deg_cutoff) %>%
  filter(n_samples>=samples.cutoff) %>%
  filter(SamplesPairing==0) %>%
  select(dataset, p_group_cor, selected, n_degs_fdr_only, n_samples, SamplesPairing) %>% 
  rename(bias=p_group_cor) %>% 
  mutate(Imbalance=as.factor(if_else(bias>bias.high.cutoff, "true", "false"))) %>% 
  mutate(SamplesPairing=as.factor(SamplesPairing)) %>% 
  ggplot(aes(y=n_degs_fdr_only, x=n_samples, size=bias, color=Imbalance)) +
  geom_point() +
  geom_smooth(method="glm", formula=my.formula <- y ~ x) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               size = rel(6),
               parse = TRUE) +  
  labs(title=paste("FDR only \nNo Paired Samples"), 
      #  subtitle = paste0(
      #    #"Max correlation coeff.:", dataset_p_group_cor_cutoff, 
      #    "Diff. genes defined by FDR only",
      #   #  "\nDatasets selection: diff genes\u2265", dataset_n_deg_cutoff,
      #   #  " and samples\u2265", samples.cutoff,
      #     "\nQuality Imbalance: low \u2264 ", bias.high.cutoff, "< high imbalance"
      #    # ,
      #    #                        "\nsamples-genes cor: ", sample_deg.all.cor,
      #    # "samples-genes cor: ", sample_deg.all.cor,
      #    # "| samples-genes cor (selected datasets): ", sample_deg.selection.cor
      #  ),
       y="Differential genes",
       x="Samples",
       size="QI Index",
       color="High QI") +
  theme_minimal(base_size = 20)+
  scale_color_manual(values=list("false"="#295D8A", "true"="#A41720"))

ggsave(filename=path, plot=deg.fdrOnly.no.paired.samples, width=12, height=10)


### composite plots
# possible plots:
# deg.samples
# deg.paired.samples
# deg.no.paired.samples
# deg.fdrOnly.samples
# deg.fdrOnly.paired.samples
# deg.fdrOnly.no.paired.samples
# Figure 3
figure <- ggarrange(
  deg.fdrOnly.samples,
  deg.samples,
  labels = c("A", "B"),
  font.label = list(size = 20, color = "black", face = "bold", family = NULL),
  common.legend = TRUE, legend = "bottom",
  ncol = 2, nrow = 1)
path <- file.path(out.dir.path, "figure-3.pdf")
ggsave(filename=path, plot=figure, width=13, height=6.5)
# suppl figure S2
figure <- ggarrange(
  deg.fdrOnly.no.paired.samples,
  deg.no.paired.samples,
  deg.fdrOnly.paired.samples,
  deg.paired.samples,
  labels = c("A", "B", "C", "D"),
  font.label = list(size = 20, color = "black", face = "bold", family = NULL),
  common.legend = TRUE, legend = "bottom",
  ncol = 2, nrow = 2)
path <- file.path(out.dir.path, "suppl-figure-S2.pdf")
ggsave(filename=path, plot=figure, width=13, height=13)


# DEGs vs BIAS
# ------------------------------------------------------------------------------
path <- file.path(out.dir.path, "deg.plow.balanced.pdf")
# n_samples.median <- median(datasets.stats$n_samples)
# n_samples.median <- 30

plot <- datasets.stats %>%
  filter(n_degs>=dataset_n_deg_cutoff) %>%
  filter(n_samples>=samples.cutoff) %>%
  filter(balanced_groups==1) %>%
  select(dataset, p_group_cor, selected, n_degs, n_samples) %>% 
  mutate(data=as.factor(if_else(n_samples>n_samples.big, "big data", if_else(n_samples>=n_samples.medium, "medium data", "small data")))) %>% 
  rename(samples=n_samples) %>% 
  ggplot(aes(y=n_degs, x=p_group_cor, size=samples, color=data)) +
  geom_point() +
  geom_smooth(method="glm", formula=my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               size = rel(6),
               parse = TRUE) +  
  labs(title=paste("Differential Genes vs Dataset Bias - balanced_groups"), 
       subtitle = paste0("Datasets selection: diff genes\u2265", dataset_n_deg_cutoff,
                        " and samples\u2265", samples.cutoff,
                        "\nData size: ", "small<", n_samples.medium, " | big>", n_samples.big, " samples"
       ),
       y="Differential genes", 
       x="Dataset Bias") +
  theme_minimal(base_size = 20)
ggsave(filename=path, plot=plot, width=12, height=10)


path <- file.path(out.dir.path, "deg.plow.homogen.pdf")
plot <- datasets.stats %>%
  filter(n_degs>=dataset_n_deg_cutoff) %>%
  filter(n_samples>=samples.cutoff) %>%
  filter(samples_homogeneity==TRUE) %>%
  select(dataset, p_group_cor, selected, n_degs, n_samples) %>% 
  mutate(data=as.factor(if_else(n_samples>n_samples.big, "big data", if_else(n_samples>=n_samples.medium, "medium data", "small data")))) %>% 
  rename(samples=n_samples) %>% 
  ggplot(aes(y=n_degs, x=p_group_cor, size=samples, color=data)) +
  geom_point() +
  geom_smooth(method="glm", formula=my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               size = rel(6),
               parse = TRUE) +  
  labs(title=paste("Differential Genes vs Dataset Bias - samples_homogeneity"), 
       subtitle = paste0("Datasets selection: diff genes\u2265", dataset_n_deg_cutoff,
                        " and samples\u2265", samples.cutoff,
                        "\nData size: ", "small<", n_samples.medium, " | big>", n_samples.big, " samples"
       ),
       y="Differential genes", 
       x="Dataset Bias") +
  theme_minimal(base_size = 20)
ggsave(filename=path, plot=plot, width=12, height=10)


path <- file.path(out.dir.path, "deg.plow.regress.all.pdf")
plot <- datasets.stats %>%
  filter(n_degs>=dataset_n_deg_cutoff) %>%
  filter(n_samples>=samples.cutoff) %>%
  select(dataset, p_group_cor, selected, n_degs, n_samples) %>% 
  mutate(data=as.factor(if_else(n_samples>n_samples.big, "big data", if_else(n_samples>=n_samples.medium, "medium data", "small data")))) %>% 
  rename(samples=n_samples) %>% 
  ggplot(aes(y=n_degs, x=p_group_cor)) +
  geom_point() +
  geom_smooth(method="glm", formula=my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               size = rel(6),
               parse = TRUE) +  
  labs(title=paste("Differential Genes vs Dataset Bias - Regress All"), 
       subtitle = paste0("Datasets selection: diff genes\u2265", dataset_n_deg_cutoff,
                        " and samples\u2265", samples.cutoff,
                        "\nData size: ", "small<", n_samples.medium, " | big>", n_samples.big, " samples"
       ),
       y="Differential genes", 
       x="Dataset Bias") +
  theme_minimal(base_size = 20)
ggsave(filename=path, plot=plot, width=12, height=10)

path <- file.path(out.dir.path, "deg.plow.pdf")
plot <- datasets.stats %>%
  filter(n_degs>=dataset_n_deg_cutoff) %>%
  filter(n_samples>=samples.cutoff) %>%
  select(dataset, p_group_cor, selected, n_degs, n_samples) %>% 
  mutate(data=as.factor(if_else(n_samples>n_samples.big, "big data", if_else(n_samples>=n_samples.medium, "medium data", "small data")))) %>% 
  rename(samples=n_samples) %>% 
  ggplot(aes(y=n_degs, x=p_group_cor, size=samples, color=data)) +
  geom_point() +
  geom_smooth(method="glm", formula=my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               size = rel(6),
               parse = TRUE) +  
  labs(title=paste("Differential Genes vs Dataset Bias"), 
       subtitle = paste0("Datasets selection: diff genes\u2265", dataset_n_deg_cutoff,
                        " and samples\u2265", samples.cutoff,
                        "\nData size: ", "small<", n_samples.medium, " | big>", n_samples.big, " samples"
       ),
       y="Differential genes", 
       x="Dataset Bias") +
  theme_minimal(base_size = 20)
ggsave(filename=path, plot=plot, width=12, height=10)


# switch color back
options(ggplot2.discrete.colour= c("#295D8A","#A41720", "#4E4459"))


# SAMPLES vs BIAS
# # ------------------------------------------------------------------------------
# path <- file.path(out.dir.path, "samples.plow.pdf")
# n_deg.median <- median(datasets.stats$n_degs)
# # n_deg.median <- 20
# my.formula <- y ~ x
# plot <- datasets.stats %>%
#   filter(n_degs>dataset_n_deg_cutoff) %>%
#   select(dataset, p_group_cor, selected, n_degs, n_samples) %>% 
#   rename(samples=n_samples) %>% 
#   rename(`differential genes`=n_degs) %>% 
#   mutate(data=as.factor(if_else(`differential genes`>n_deg.median, "many degs", "few degs"))) %>% 
#   ggplot(aes(x=samples, y=p_group_cor, size=`differential genes`)) +
#   geom_point() +
#   geom_smooth(method="glm", formula=my.formula) +
#   stat_poly_eq(formula = my.formula, 
#                aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
#                size = rel(6),
#                parse = TRUE) +  
#   labs(title=paste("Samples vs Dataset Bias"), 
#        subtitle = paste("Datasets selection: diff genes>", dataset_n_deg_cutoff
#                         # ,
#                         # "| Big data definition: samples>", n_samples.median
#                         # ,
#                         # "\np-genes cor (big data): ", sample_deg.bigdata.cor,
#                         # "| p-genes cor (small data): ", sample_deg.smalldata.cor
#        ),
#        x="Samples", 
#        y="Dataset Bias") +
#   theme_minimal(base_size = 20)
# ggsave(filename=path, plot=plot, width=12, height=10)
# 

# ==============================================================================
# CLUSTERING QUALITY
# ==============================================================================
median.p <- median(datasets.stats$p_group_cor)
dtmp <- datasets.stats%>% 
  rename(quality_bias=p_group_cor) %>%
  mutate(diff_dunn=dunn_noout-dunn_all) %>%
  mutate(diff_wbrat=wbrat_noout-wbrat_all) %>%  
  mutate(diff_pgamma=pgamma_noout-pgamma_all) %>% 
  mutate(diff_dunn_nogenes=dunn_noout_nocor-dunn_noout_norand) %>%  
  mutate(diff_wbrat_nogenes=wbrat_noout_nocor-wbrat_noout_norand) %>% 
  mutate(diff_pgamma_nogenes=pgamma_noout_nocor-pgamma_noout_norand) %>%
  mutate(big_dataset=as.factor(if_else(n_samples>n_samples.big, "big data", "small data"))) %>%
  mutate(bias_type=as.factor(if_else(quality_bias<=bias.high.cutoff, "low bias", "high bias"))) %>% 
  filter(n_degs>=dataset_n_deg_cutoff)

diff_dunn.cor   <- round(cor(dtmp$diff_dunn, dtmp$quality_bias), 2)
diff_wbrat.cor  <- round(cor(dtmp$diff_wbrat, dtmp$quality_bias), 2)
diff_pgamma.cor <- round(cor(dtmp$diff_pgamma, dtmp$quality_bias), 2)

diff_dunn.smallData.cor   <- round(cor(dtmp[dtmp$big_dataset=="small data",]$diff_dunn, dtmp[dtmp$big_dataset=="small data",]$quality_bias), 2)
diff_wbrat.smallData.cor  <- round(cor(dtmp[dtmp$big_dataset=="small data",]$diff_wbrat, dtmp[dtmp$big_dataset=="small data",]$quality_bias), 2)
diff_pgamma.smallData.cor <- round(cor(dtmp[dtmp$big_dataset=="small data",]$diff_pgamma, dtmp[dtmp$big_dataset=="small data",]$quality_bias), 2)

diff_dunn.bigData.cor   <- round(cor(dtmp[dtmp$big_dataset=="big data",]$diff_dunn, dtmp[dtmp$big_dataset=="big data",]$quality_bias), 2)
diff_wbrat.bigData.cor  <- round(cor(dtmp[dtmp$big_dataset=="big data",]$diff_wbrat, dtmp[dtmp$big_dataset=="big data",]$quality_bias), 2)
diff_pgamma.bigData.cor <- round(cor(dtmp[dtmp$big_dataset=="big data",]$diff_pgamma, dtmp[dtmp$big_dataset=="big data",]$quality_bias), 2)

diff_dunn_nogenes.cor  <- round(cor(dtmp$diff_dunn_nogenes, dtmp$quality_bias), 2)
diff_wbrat_nogenes.cor <- round(cor(dtmp$diff_wbrat_nogenes, dtmp$quality_bias), 2)
diff_pgamma_nogenes.cor <- round(cor(dtmp$diff_pgamma_nogenes, dtmp$quality_bias), 2)

diff_dunn_nogenes.smallData.cor  <- round(cor(dtmp[dtmp$big_dataset=="small data",]$diff_dunn_nogenes, dtmp[dtmp$big_dataset=="small data",]$quality_bias), 2)
diff_wbrat_nogenes.smallData.cor <- round(cor(dtmp[dtmp$big_dataset=="small data",]$diff_wbrat_nogenes, dtmp[dtmp$big_dataset=="small data",]$quality_bias), 2)
diff_pgamma_nogenes.smallData.cor <- round(cor(dtmp[dtmp$big_dataset=="small data",]$diff_pgamma_nogenes, dtmp[dtmp$big_dataset=="small data",]$quality_bias), 2)

diff_dunn_nogenes.bigData.cor  <- round(cor(dtmp[dtmp$big_dataset=="big data",]$diff_dunn_nogenes, dtmp[dtmp$big_dataset=="big data",]$quality_bias), 2)
diff_wbrat_nogenes.bigData.cor <- round(cor(dtmp[dtmp$big_dataset=="big data",]$diff_wbrat_nogenes, dtmp[dtmp$big_dataset=="big data",]$quality_bias), 2)
diff_pgamma_nogenes.bigData.cor <- round(cor(dtmp[dtmp$big_dataset=="big data",]$diff_pgamma_nogenes, dtmp[dtmp$big_dataset=="big data",]$quality_bias), 2)


# Dunn Index Plot
# ------------------------------------------------------------------------------
path <- file.path(pca.dir.path, "pca.dunn.pdf")
ymax <- round(max(dtmp$dunn_noout - dtmp$dunn_all) * 1.8, 1)
ymin <- round(min(dtmp$dunn_noout - dtmp$dunn_all) * 1.8, 1)
ymax <- if_else(ymax<0, ymax-0.25, ymax)
ymin <- if_else(ymin<0, ymin-0.25, ymin)
plot <- dtmp %>% 
  # filter(n_degs>dataset_n_deg_cutoff) %>%
  ggplot(aes(x=reorder(dataset, diff_dunn), y=diff_dunn)) +
#  geom_bar(stat='identity', fill="#295D8A") +
  geom_bar(aes(fill=quality_bias), stat='identity') +
  ylim(c(ymin, ymax)) +
  geom_text(aes(label=paste(round(dunn_noout, 2), "-", round(dunn_all, 2))), hjust=ifelse(dtmp$diff_dunn>=0, -0.25, 1.25), size=6) +
  labs(title=paste("PCA Clustering Evaluation Boost"),
       subtitle=paste0("datasets selection: diff genes\u2265", dataset_n_deg_cutoff,
         "\nbias-diff correlation=", diff_dunn.cor),
       x="Datasets",
       y="Difference in Dunn Index (no outliers - all samples)") + 
  coord_flip() +
  theme_minimal(base_size = 20)

ggsave(filename=path, plot=plot, width=12, height=15)

# wbrat Index Plot
# ------------------------------------------------------------------------------
path <- file.path(pca.dir.path, "pca.wbrat.pdf")
ymax <- round(max(dtmp$wbrat_noout - dtmp$wbrat_all) * 1.8, 1)
ymin <- round(min(dtmp$wbrat_noout - dtmp$wbrat_all) * 1.8, 1)
ymax <- if_else(ymax<0, ymax-0.25, ymax)
ymin <- if_else(ymin<0, ymin-0.25, ymin)
plot <- dtmp %>% 
  ggplot(aes(x=reorder(dataset, diff_wbrat), y=diff_wbrat)) +
#  geom_bar(stat='identity', fill="#295D8A") +
  geom_bar(aes(fill=quality_bias), stat='identity') +
  ylim(c(ymin, ymax)) +
#  geom_text(aes(label=paste(round(wbrat_noout, 2), "-", round(wbrat_all, 2))), hjust=-0.25, size=7) +
  geom_text(aes(label=paste(round(wbrat_noout, 2), "-", round(wbrat_all, 2))), hjust=ifelse(dtmp$diff_wbrat>=0, -0.25, 1.25), size=6) +
  labs(title=paste("PCA Clustering Evaluation Boost"),
       subtitle=paste0("datasets selection: diff genes\u2265", dataset_n_deg_cutoff,
                       "\nbias-diff correlation=", diff_wbrat.cor),
       x="Datasets",
       y="Difference in wbrat Index (no outliers - all samples)") + 
  coord_flip() +
  theme_minimal(base_size = 20)
ggsave(filename=path, plot=plot, width=12, height=15)

# Pearson's Gamma Index Plot
# ------------------------------------------------------------------------------
path <- file.path(pca.dir.path, "pca.pgamma.pdf")
ymax <- round(max(dtmp$pgamma_noout - dtmp$pgamma_all) * 1.8, 3)
ymin <- round(min(dtmp$pgamma_noout - dtmp$pgamma_all) * 1.8, 3)
ymax <- if_else(ymax<0, ymax-0.25, ymax)
ymin <- if_else(ymin<0, ymin-0.25, ymin)
plot <- dtmp %>% 
  ggplot(aes(x=reorder(dataset, diff_pgamma), y=diff_pgamma)) +
#  geom_bar(stat='identity', fill="#295D8A") +
  geom_bar(aes(fill=quality_bias), stat='identity') +
  ylim(c(ymin, ymax)) +
  geom_text(aes(label=paste(round(pgamma_noout, 2), "-", round(pgamma_all, 2))), hjust=ifelse(dtmp$diff_pgamma>=0, -0.25, 1.25), size=6) +
  labs(title=paste("PCA Clustering Evaluation Boost (Normalized Gamma)"),
       subtitle=paste0("datasets selection: diff genes\u2265", dataset_n_deg_cutoff,
                       "\nbias-diff correlation=", diff_pgamma.cor),
       x="Datasets",
       y="Difference in Normalized Gamma (dataset without outliers vs with all samples)") + 
  coord_flip() +
  theme_minimal(base_size = 20)
ggsave(filename=path, plot=plot, width=15, height=15)

# Difference with all indices
# ------------------------------------------------------------------------------
path <- file.path(pca.dir.path, "pca.indices.pdf")
ymax <- round(max(dtmp$diff_dunn, dtmp$diff_wbrat, dtmp$diff_pgamma) * 1.8, 3)
ymin <- round(min(dtmp$diff_dunn, dtmp$diff_wbrat, dtmp$diff_pgamma) * 1.8, 3)
ymax <- if_else(ymax<0, ymax-0.25, ymax)
ymin <- if_else(ymin<0, ymin-0.25, ymin)

dtmp2 <-dtmp %>%
  select(dataset, diff_dunn, diff_wbrat, diff_pgamma, quality_bias, big_dataset) %>%
  gather(diff_dunn, diff_wbrat, diff_pgamma, key = "feature", value = "value") %>%
  mutate(feature=as.factor(feature)) 

plot <- dtmp2 %>%
  ggplot(aes(x = reorder(dataset, value), y = value)) +
  geom_bar(aes(fill=quality_bias), stat='identity') +
  ylim(c(ymin, ymax)) +
  geom_text(aes(label=paste(round(value, 2))), hjust=ifelse(dtmp2$value>=0, -0.25, 1.25), size=6) +
  labs(title=paste("PCA Clustering Evaluation Boost"),
       subtitle=paste0(
         "datasets selection: diff genes\u2265", dataset_n_deg_cutoff,
         "\nbias-diff cor: ", "Dunn=", diff_dunn.cor, " | wbrat=", diff_wbrat.cor, " | Gamma=", diff_pgamma.cor,
           "\nbias-diff cor (small data): ", "Dunn=", diff_dunn.smallData.cor, " | wbrat=", diff_wbrat.smallData.cor, " | Gamma=", diff_pgamma.smallData.cor,
           "\nbias-diff cor (big data): ", "Dunn=", diff_dunn.bigData.cor, " | wbrat=", diff_wbrat.bigData.cor, " | Gamma=", diff_pgamma.bigData.cor),
       x="Datasets",
       y="Difference in Dunn Index or Normalized Gamma Index") + 
  coord_flip() +
  facet_wrap( ~ feature, scales="free", drop=TRUE) +
  theme_minimal(base_size = 20)
ggsave(filename=path, plot=plot, width=20, height=20)


# # Difference with all indices - SCATTER PLOT
# # ------------------------------------------------------------------------------
# path <- file.path(out.dir.path, "pca.indices.plow.scatter.pdf")
# 
# dtmp2 <-dtmp %>%
#   filter(deg>dataset_n_deg_cutoff) %>%
#   filter(n_samples>samples.cutoff) %>%
#   # select(dataset, diff_dunn, diff_wbrat, diff_pgamma, quality_bias, big_dataset) %>%
#   gather(diff_dunn, diff_wbrat, diff_pgamma, key = "feature", value = "value") %>%
#   mutate(feature=as.factor(feature)) 
# plot <- dtmp2 %>%
#   ggplot(aes(x = quality_bias , y = value)) +
#   geom_point() +
#   geom_smooth(method="glm", formula=y~x)+
#   stat_poly_eq(formula = my.formula, 
#                aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
#                size = rel(4),
#                parse = TRUE) +  
#   facet_wrap(~big_dataset*feature)+
#   theme_minimal(base_size = 20)
# ggsave(filename=path, plot=plot, width=10, height=8)
# 
# 
# 
# # Difference with all indices - SCATTER PLOT
# # ------------------------------------------------------------------------------
# path <- file.path(out.dir.path, "pca.indices.plow.scatter.pdf")
# 
# dtmp2 <-dtmp %>%
#   filter(deg>dataset_n_deg_cutoff) %>%
#   filter(n_samples>samples.cutoff) %>%
#   rename(`Dunn index`=diff_dunn) %>% 
#   rename(`wbrat index`=diff_wbrat) %>% 
#   rename(`Gamma index`=diff_pgamma) %>% 
#   # select(dataset, diff_dunn, diff_wbrat, diff_pgamma, quality_bias, big_dataset) %>%
#   gather(`Dunn index`, `wbrat index`, `Gamma index`, key = "feature", value = "value") %>% 
#   mutate(feature=as.factor(feature)) %>% 
#   mutate(sign=as.factor(if_else(value>=0.1, "positive", if_else(value>=-0.1, "near-zero", "negative"))))
# plot <- dtmp2 %>%
#   ggplot(aes(x = quality_bias , y = value, color=sign)) +
#   geom_point() +
#   ylim(c(-1, 1))+
#   geom_hline(yintercept=0, linetype="dashed", color = "red") +
#   geom_smooth(aes(x = quality_bias , y = value), 
#               method="glm", 
#               formula = my.formula,
#               inherit.aes=F)+
#   stat_poly_eq(mapping=aes(x = quality_bias , y = value, label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
#                size = rel(4),
#                parse = TRUE, 
#                formula = my.formula,
#                inherit.aes=F) +
#   facet_wrap(~big_dataset*feature)+
#   facet_wrap(~feature)+
#   theme_minimal(base_size = 20) +
#   labs(title=paste("PCA Clustering Evaluation Boost"),
#        subtitle="With vs without outlier Samples",
#        x="Dataset bias",
#        y="Index difference\nwith - without outlier samples removal") 
# 
# ggsave(filename=path, plot=plot, width=12, height=6)
# 

#===============================================================================
# Dunn Index Plot When Romoving Genes
# ------------------------------------------------------------------------------
path <- file.path(pca.dir.path, "pca.genes.removal.dunn.pdf")
ymax <- round(max(dtmp$diff_dunn_nogenes) * 1.8, 1)
ymin <- round(min(dtmp$diff_dunn_nogenes) * 1.8, 1)
ymax <- if_else(ymax<0, ymax-0.25, ymax)
ymin <- if_else(ymin<0, ymin-0.25, ymin)
plot <- dtmp %>% 
  ggplot(aes(x=reorder(dataset, diff_dunn_nogenes), y=diff_dunn_nogenes)) +
  geom_bar(aes(fill=quality_bias), stat='identity') +
  ylim(c(ymin, ymax)) +
#  geom_text(aes(label=paste(round(dunn_noout_nocor, 2), "-", round(dunn_noout_norand, 2))), hjust=-0.25, size=7) +
  geom_text(aes(label=paste(round(dunn_noout_nocor, 2), "-", round(dunn_noout_norand, 2))), hjust=ifelse(dtmp$diff_dunn_nogenes>=0, -0.25, 1.25), size=6) +
  labs(title=paste("PCA Clustering Evaluation Boost"),
       subtitle=paste0("datasets selection: diff genes\u2265", dataset_n_deg_cutoff,
                       "\nbias-diff correlation=", diff_dunn_nogenes.cor),
       x="Datasets",
       y="Difference in Dunn Index without outliers (correlated - random genes removal)") + 
  coord_flip() +
  theme_minimal(base_size = 20)
ggsave(filename=path, plot=plot, width=12, height=15)

# wbrat Index Plot When Romoving Genes
# ------------------------------------------------------------------------------
path <- file.path(pca.dir.path, "pca.genes.removal.wbrat.pdf")
ymax <- round(max(dtmp$diff_wbrat_nogenes) * 1.8, 1)
ymin <- round(min(dtmp$diff_wbrat_nogenes) * 1.8, 1)
ymax <- if_else(ymax<0, ymax-0.25, ymax+0.25)
ymin <- if_else(ymin<0, ymin-0.25, ymin)
plot <- dtmp %>% 
  ggplot(aes(x=reorder(dataset, diff_wbrat_nogenes), y=diff_wbrat_nogenes)) +
  geom_bar(aes(fill=quality_bias), stat='identity') +
  ylim(c(ymin, ymax)) +
#  geom_text(aes(label=paste(round(dunn_noout_nocor, 2), "-", round(dunn_noout_norand, 2))), hjust=-0.25, size=7) +
  geom_text(aes(label=paste(round(wbrat_noout_nocor, 2), "-", round(wbrat_noout_norand, 2))), hjust=ifelse(dtmp$diff_wbrat_nogenes>=0, -0.25, 1.25), size=6) +
  labs(title=paste("PCA Clustering Evaluation Boost"),
       subtitle=paste0("datasets selection: diff genes\u2265", dataset_n_deg_cutoff,
                       "\nbias-diff correlation=", diff_wbrat_nogenes.cor),
       x="Datasets",
       y="Difference in wbrat Index without outliers (correlated - random genes removal)") + 
  coord_flip() +
  theme_minimal(base_size = 20)
ggsave(filename=path, plot=plot, width=12, height=15)


# Pearson's Gamma Index Plot When Romoving Genes
# ------------------------------------------------------------------------------
path <- file.path(pca.dir.path, "pca.genes.removal.gamma.pdf")
ymax <- round(max(dtmp$diff_pgamma_nogenes) * 1.8, 1)
ymin <- round(min(dtmp$diff_pgamma_nogenes) * 1.8, 1)
plot <- dtmp %>% 
  ggplot(aes(x=reorder(dataset, diff_pgamma_nogenes), y=diff_pgamma_nogenes)) +
#  geom_bar(stat='identity', fill="#295D8A") +
  geom_bar(aes(fill=quality_bias), stat='identity') +
  ylim(c(ymin, ymax)) +
  geom_text(aes(label=paste(round(pgamma_noout_nocor, 2), "-", round(pgamma_noout_norand, 2))), hjust=ifelse(dtmp$diff_pgamma_nogenes>=0, -0.25, 1.25), size=6) +
  labs(title=paste("PCA Clustering Evaluation Boost"),
       subtitle=paste0("datasets selection: diff genes\u2265", dataset_n_deg_cutoff,
                       "\nbias-diff correlation=", diff_pgamma_nogenes.cor),
       x="Datasets",
       y="Difference in Gamma Index without outliers (correlated - random genes removal)") + 
  coord_flip() +
  theme_minimal(base_size = 20)
ggsave(filename=path, plot=plot, width=15, height=15)


# Difference with all indices When Romoving Genes
# ------------------------------------------------------------------------------
#path <- file.path(out.dir.path, "pca.genes.removal.indices.pdf")
#plot <- d %>% 
#  select(dataset, dunn_noout_norand, dunn_noout_nocor, wbrat_noout_norand, wbrat_noout_nocor, pgamma_noout_norand, pgamma_noout_nocor, p_group_cor) %>%
#  rename(quality_bias=p_group_cor) %>%
#  mutate(diff_dunn = dunn_noout_nocor - dunn_noout_norand) %>%
#  mutate(diff_wbrat = wbrat_noout_nocor - wbrat_noout_norand) %>%
#  mutate(diff_pgamma = pgamma_noout_nocor - pgamma_noout_norand) %>%
#  select(dataset, diff_dunn, diff_wbrat, diff_pgamma, quality_bias) %>%
#  gather(diff_dunn, diff_wbrat, diff_pgamma, key = "feature", value = "value") %>%
#  ggplot(aes(x = reorder(dataset, value), y = value)) +
##  geom_bar(stat='identity', fill="#295D8A") +
#  geom_bar(aes(fill=quality_bias), stat='identity') +
#  geom_text(aes(label=paste(round(value, 2))), hjust=-0.25, size=7) +
#  labs(title=paste("PCA Clustering Evaluation Boost"),
#       x="Datasets",
#       y="Difference in Dunn Index or Normalized Gamma Index When Romoving Genes") + 
#  coord_flip() +
#  facet_wrap( ~ feature) +
#  theme_minimal(base_size = 20)
#ggsave(filename=path, plot=plot, width=20, height=15)

path <- file.path(pca.dir.path, "pca.genes.removal.indices.pdf")
ymax <- round(max(dtmp$diff_dunn_nogenes, dtmp$diff_wbrat_nogenes, dtmp$diff_pgamma_nogenes) * 1.8, 3)
ymin <- round(min(dtmp$diff_dunn_nogenes, dtmp$diff_wbrat_nogenes, dtmp$diff_pgamma_nogenes) * 1.8, 3)
ymax <- if_else(ymax<0, ymax-0.25, ymax)
ymin <- if_else(ymin<0, ymin-0.25, ymin)

dtmp2 <- dtmp %>%
#  select(dataset, diff_dunn_nogenes, diff_wbrat_nogenes, diff_pgamma_nogenes, quality_bias, big_dataset) %>%
  gather(diff_dunn_nogenes, diff_wbrat_nogenes, diff_pgamma_nogenes, key = "feature", value = "value")

plot <- dtmp2 %>%
  ggplot(aes(x = reorder(dataset, value), y = value)) +
  geom_bar(aes(fill=quality_bias), stat='identity') +
  ylim(c(ymin, ymax)) +
  geom_text(aes(label=paste(round(value, 2))), hjust=ifelse(dtmp2$value>=0, -0.25, 1.25), size=6) +
  labs(title=paste("PCA Clustering Evaluation Boost"),
       subtitle=paste0(
         "datasets selection: diff genes\u2265", dataset_n_deg_cutoff,
         "\ncorrelation(bias vs quality): ", "Dunn=", diff_dunn_nogenes.cor, " | wbrat=", diff_wbrat_nogenes.cor, " | Gamma=", diff_pgamma_nogenes.cor,
         "\nbias-diff cor (small data): ", "Dunn=", diff_dunn_nogenes.smallData.cor, " | wbrat=", diff_wbrat_nogenes.smallData.cor, " | Gamma=", diff_pgamma_nogenes.smallData.cor,
         "\nbias-diff cor (big data): ", "Dunn=", diff_dunn_nogenes.bigData.cor, " | wbrat=", diff_wbrat_nogenes.bigData.cor, " | Gamma=", diff_pgamma_nogenes.bigData.cor),
       x="Datasets",
       y="Difference in Dunn Index or Normalized Gamma Index When Romoving Genes") + 
  coord_flip() +
#  facet_wrap( ~ big_dataset*bias_type*feature, scales="free_x") +
  facet_wrap( ~ feature, scales="free", drop=TRUE) +
  theme_minimal(base_size = 20)
ggsave(filename=path, plot=plot, width=20, height=25)



# # Difference with all indices - SCATTER PLOT
# # ------------------------------------------------------------------------------
# path <- file.path(out.dir.path, "pca.genes.removal.indices.plow.scatter.pdf")
# 
# dtmp2 <-dtmp %>%
#   filter(deg>dataset_n_deg_cutoff) %>%
#   filter(n_samples>samples.cutoff) %>%
#   rename(`Dunn index`=diff_dunn_nogenes) %>% 
#   rename(`wbrat index`=diff_wbrat_nogenes) %>% 
#   rename(`Gamma index`=diff_pgamma_nogenes) %>% 
#   # select(dataset, diff_dunn, diff_wbrat, diff_pgamma, quality_bias, big_dataset) %>%
#   gather(`Dunn index`, `wbrat index`, `Gamma index`, key = "feature", value = "value") %>% 
#   mutate(feature=as.factor(feature)) %>% 
#   mutate(sign=as.factor(if_else(value>=0.1, "positive", if_else(value>=-0.1, "near-zero", "negative"))))
# plot <- dtmp2 %>%
#   ggplot(aes(x = quality_bias , y = value, color=sign)) +
#   geom_point() +
#   ylim(c(-1, 1))+
#   geom_hline(yintercept=0, linetype="dashed", color = "red") +
#   geom_smooth(aes(x = quality_bias , y = value), 
#               method="glm", 
#               formula = my.formula,
#               inherit.aes=F)+
#   stat_poly_eq(mapping=aes(x = quality_bias , y = value, label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
#                size = rel(4),
#                parse = TRUE, 
#                formula = my.formula,
#                inherit.aes=F) +
#   facet_wrap(~big_dataset*feature)+
#   facet_wrap(~feature)+
#   theme_minimal(base_size = 20) +
#   labs(title=paste("PCA Clustering Evaluation Boost"),
#        subtitle="Biased vs Random Genes Removal | Outlier Samples Removed",
#        x="Dataset bias",
#        y="Index difference\nbiased - random genes removal") 
#   
# ggsave(filename=path, plot=plot, width=12, height=6)

# # Entropy Index Plot
# # ------------------------------------------------------------------------------
# path <- file.path(out.dir.path, "pca.entropy.pdf")
# ymax <- round(max(d$entropy_noout - d$entropy_all) * 1.5, 3)
# ymin <- round(min(d$entropy_noout - d$entropy_all) * 1.5, 3)
# plot <- d %>% 
#   select(dataset, entropy_all, entropy_noout) %>% 
#   mutate(diff_entropy=entropy_noout-entropy_all) %>% 
#   ggplot(aes(x=reorder(dataset, diff_entropy), y=diff_entropy)) +
#   geom_bar(stat='identity', fill="#295D8A") +
#   ylim(c(ymin, ymax)) +
#   geom_text(aes(label=paste(round(entropy_noout, 2), "-", round(entropy_all, 2))), hjust=-0.25) +
#   labs(y="Difference in Entropy (dataset without outliers vs with all samples)",
#        x="",
#        title=paste("PCA Clustering Evaluation Boost (Entropy)")) + 
#   coord_flip() +
#   theme_minimal() +
#   ggsave(path, width=15.2, height=12.8)
# 
# # Difference with all indices
# # ------------------------------------------------------------------------------
# path <- file.path(out.dir.path, "pca.indices.pdf")
# plot <- d %>% 
#   select(dataset, dunn_all, dunn_noout, pgamma_all, 
#          pgamma_noout, entropy_all, entropy_noout ) %>%
#   mutate(diff_dunn = dunn_noout - dunn_all) %>%
#   mutate(diff_pgamma = pgamma_noout - pgamma_all) %>%
#   mutate(diff_entropy = entropy_noout - entropy_all) %>%
#   select(dataset, diff_dunn, diff_pgamma, diff_entropy) %>%
#   gather(diff_dunn, diff_pgamma, diff_entropy, key = "feature", value = "value") %>%
#   ggplot(aes(x = reorder(dataset, value), y = value)) +
#   geom_bar(stat = "identity", fill="#295D8A") +
#   coord_flip() +
#   facet_wrap( ~ feature) +
#   theme_minimal() +
#   ggsave(path, width=15.2, height=12.8)

# Difference with all indices - SCATTER PLOT
# ------------------------------------------------------------------------------
samples.cutoff2 <- samples.cutoff
# bias.vline <- 0.1
dtmp2 <-dtmp %>%
  filter(n_samples>=samples.cutoff2) %>%
  gather(diff_dunn, diff_wbrat, diff_pgamma, 
         diff_dunn_nogenes, diff_wbrat_nogenes, diff_pgamma_nogenes, 
         key = "feature", value = "value") %>%
#  filter(feature!="diff_wbrat") %>% 
#  filter(feature!="diff_wbrat_nogenes") %>% 
  mutate(analysis=feature) %>%
  mutate(analysis=if_else(grepl("nogenes",analysis), "no biased genes", "no outliers")) %>%
  mutate(analysis=factor(analysis, levels = c("no outliers", "no biased genes"))) %>% 
  mutate(feature=if_else(startsWith(feature, "diff_pgamma"), "Gamma index", feature)) %>% 
  mutate(feature=if_else(startsWith(feature, "diff_wbrat"), "wbrat index", feature)) %>% 
  mutate(feature=if_else(startsWith(feature, "diff_dunn"), "Dunn index", feature)) %>% 
  mutate(feature=as.factor(feature)) %>%
  mutate(sign=if_else(value>=near.zero.cutoff, "positive", if_else(value>(-near.zero.cutoff), "near-zero", "negative"))) %>% 
  mutate(sign=factor(sign, levels=c("negative", "positive", "near-zero"))) %>% 
  mutate(data=as.factor(if_else(n_samples>n_samples.big, "big data", if_else(n_samples>=n_samples.medium, "medium data", "small data"))))
  
path <- file.path(out.dir.path, "pca.indices.plow.scatter.pdf")
plot <- dtmp2 %>%
  ggplot(aes(x = quality_bias , y = value, color=sign)) +
  geom_point() +
  ylim(c(-1, 1))+
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  # geom_vline(xintercept=bias.vline <- 0.20, linetype="dashed", color = "red") +
  geom_smooth(aes(x = quality_bias , y = value), 
              method="glm", 
              formula = my.formula,
              inherit.aes=F)+
  stat_poly_eq(mapping=aes(x = quality_bias , y = value, label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               size = rel(2),
               parse = TRUE, 
               formula = my.formula,
               inherit.aes=F) +
  facet_grid(rows=vars(analysis), cols=vars(feature))+
  # theme_minimal(base_size = 10) +
  theme_minimal() +
  theme(legend.position="bottom")+
  labs(title=paste("PCA Clustering Evaluation Boost"),
       subtitle=paste0(
         "With vs without biased genes or outlier samples",
         "\nDatasets selection: diff genes\u2265", dataset_n_deg_cutoff,
         " and samples\u2265", samples.cutoff2,
         "\nSign: ", -near.zero.cutoff, "<near-zero<", near.zero.cutoff
       ),
       x="Dataset bias",
       y="Index difference") 
plot
ggsave(filename=path, plot=plot, width=5, height=5)

path <- file.path(out.dir.path, "pca.indices.plow.scatter.extraFacets.pdf")
plot <- plot+
  # facet_grid(rows=vars(analysis, big_dataset), cols=vars(feature))+
  facet_grid(cols=vars(analysis, data), rows=vars(feature))+
  labs(subtitle=paste0(
    "With vs without biased genes or outlier samples",
    "\nDatasets selection: diff genes\u2265", dataset_n_deg_cutoff,
    " and samples\u2265", samples.cutoff2,
    "\nSign: ", -near.zero.cutoff, "<near-zero<", near.zero.cutoff,
    " | Data  size: ", "small<", n_samples.medium, " & big>", n_samples.big
  ))+
  # theme_minimal(base_size = 11)+
  theme(legend.position="bottom")
ggsave(filename=path, plot=plot, width=10, height=5)



# Difference with all indices - PROPORTIONS OF POSITIVES
# ------------------------------------------------------------------------------
path <- file.path(out.dir.path, "pca.indices.positives.pdf")
plot <- dtmp2 %>%
  ggplot(aes(x=bias_type, fill=sign)) +
    geom_bar() +
    geom_text(stat='count', aes(label=..count..),position = position_stack(vjust = 0.5)) +
    facet_grid(rows=vars(analysis), cols=vars(feature))+
    theme_minimal() +
  theme(legend.position="bottom")+
  labs(title=paste("PCA Clustering Evaluation Boost"),
       subtitle=paste0(
         "Sign of clustering indices difference (removal vs control)",
         "\nDatasets selection: diff genes\u2265", dataset_n_deg_cutoff,
         " and samples\u2265", samples.cutoff2,
         "\nDataset bias: low \u2264 ", bias.high.cutoff, "< high",
         " | Sign: ", -near.zero.cutoff, "<near-zero<", near.zero.cutoff
       ),
       x="",
       y="Datasets"
  ) 
ggsave(filename=path, plot=plot, width=5, height=5)

path <- file.path(out.dir.path, "pca.indices.positives.extraFacets.pdf")
plot <- plot+
  # facet_grid(rows=vars(analysis, big_dataset), cols=vars(feature))+
  facet_grid(cols=vars(analysis, data), rows=vars(feature))+
  labs(subtitle=paste0(
    "Sign of clustering indices difference (removal vs control)",
    "\nDatasets selection: diff genes\u2265", dataset_n_deg_cutoff,
    " and samples\u2265", samples.cutoff2,
    "\nDataset bias: low \u2264 ", bias.high.cutoff, "< high",
    " | Sign: ", -near.zero.cutoff, "<near-zero<", near.zero.cutoff,
    " | Data size: ", "small<", n_samples.medium, " & big>", n_samples.big
    
  )) 
ggsave(filename=path, plot=plot, width=10, height=5)



################################################################################
## Subset Analysis
################################################################################
datasets.table <- read_csv(datasets.path, show_col_types=F) %>% 
  rename(mesh_terms=`mesh_terms (pipe-separated list)`)
dataids.list <- (datasets.table %>% filter(select==1 & batches==0))$GEO_Series
# take only subsets
dataids.list <- dataids.list[grepl("sub", dataids.list)]

datasets.stats0 <- tibble()
for(DATAID in dataids.list){
  file <- file.path("output", "main", "qualityVsExp", DATAID, paste0(DATAID, ".cor.tsv"))
  if (file.exists(file)){
    datasets.stats0 <- bind_rows(list(datasets.stats0, read_tsv(file, show_col_types=F)))
  }
}

datasets.stats <- datasets.stats0 %>% 
  # left_join(deg.stats.table, by="dataset") %>% 
  mutate(n_degs=if_else(is.na(n_degs), 0, n_degs)) %>%
  mutate(n_degs_fdr_only=if_else(is.na(n_degs_fdr_only), 0, n_degs_fdr_only)) %>%
  filter(n_degs>0) %>% 
  # mutate(selected=!(n_degs<dataset_n_deg_cutoff | p_group_cor>=dataset_p_group_cor_cutoff))
  mutate(selected=n_degs>=dataset_n_deg_cutoff & p_group_cor<dataset_p_group_cor_cutoff) %>% 
  # mutate(selected_krusk=n_degs>=dataset_n_deg_cutoff & p_group_kruskal>0.05) %>% 
  left_join(datasets.table %>% rename(dataset=GEO_Series) , by = "dataset")

# DEGs vs BIAS
# ------------------------------------------------------------------------------
# fomrula
my.formula <- y ~ x
poly.formula <- y ~ poly(x, 3, raw=T)

# path <- file.path(out.dir.path, "deg.plow.homogen.pdf")
# plot <- datasets.stats %>%
#   filter(n_degs>=dataset_n_deg_cutoff) %>%
#   filter(n_samples>=samples.cutoff) %>%
#   filter(`Samples homogeneity`==TRUE) %>%
#   select(dataset, p_group_cor, selected, n_degs, n_samples) %>% 
#   mutate(data=as.factor(if_else(n_samples>n_samples.big, "big data", if_else(n_samples>=n_samples.medium, "medium data", "small data")))) %>% 
#   rename(samples=n_samples) %>% 
#   ggplot(aes(y=n_degs, x=p_group_cor, size=samples, color=data)) +
#   geom_point() +
#   geom_smooth(method="glm", formula=my.formula) +
#   stat_poly_eq(formula = my.formula, 
#                aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
#                size = rel(6),
#                parse = TRUE) +  
#   labs(title=paste("Differential Genes vs Dataset Bias - `Samples homogeneity`"), 
#        subtitle = paste0("Datasets selection: diff genes\u2265", dataset_n_deg_cutoff,
#                         " and samples\u2265", samples.cutoff,
#                         "\nData size: ", "small<", n_samples.medium, " | big>", n_samples.big, " samples"
#        ),
#        y="Differential genes", 
#        x="Dataset Bias") +
#   theme_minimal(base_size = 20)
# ggsave(filename=path, plot=plot, width=12, height=10)

# path <- file.path(out.dir.path, "deg.plow.subsets.poly.pdf")
# plot <- datasets.stats %>%
#   mutate(simi=ifelse(str_detect(dataset, "simi"), "simi", "diff")) %>%
#   # filter(n_degs>=dataset_n_deg_cutoff) %>%
#   filter(n_samples>=samples.cutoff) %>%
#   select(dataset, p_group_cor, selected, n_degs, n_samples, simi) %>% 
#   mutate(data=as.factor(if_else(n_samples>n_samples.big, "big data", if_else(n_samples>=n_samples.medium, "medium data", "small data")))) %>% 
#   rename(samples=n_samples) %>% 
#   ggplot(aes(y=n_degs, x=p_group_cor, color=simi)) +
#   geom_point() +
#   geom_smooth(method="glm", formula=poly.formula) +
#   stat_poly_eq(formula = poly.formula, 
#                aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
#                size = rel(6),
#                parse = TRUE) +  
#   labs(title=paste("Differential Genes vs Dataset Bias - Regress All"), 
#        subtitle = paste0("Datasets selection: diff genes\u2265", dataset_n_deg_cutoff,
#                         " and samples\u2265", samples.cutoff,
#                         "\nData size: ", "small<", n_samples.medium, " | big>", n_samples.big, " samples"
#        ),
#        y="Differential genes", 
#        x="Dataset Bias") +
#   theme_minimal(base_size = 20)
# ggsave(filename=path, plot=plot, width=12, height=10)


path <- file.path(out.dir.path, "deg.plow.subsets.pdf")
plot <- datasets.stats %>%
  mutate(simi=ifelse(str_detect(dataset, "simi"), "simi", "diff")) %>%
  # filter(n_degs>=dataset_n_deg_cutoff) %>%
  filter(n_samples>=samples.cutoff) %>%
  select(dataset, p_group_cor, selected, n_degs, n_samples, simi) %>% 
  mutate(data=as.factor(if_else(n_samples>n_samples.big, "big data", if_else(n_samples>=n_samples.medium, "medium data", "small data")))) %>% 
  rename(samples=n_samples) %>% 
  ggplot(aes(y=n_degs, x=p_group_cor, color=simi)) +
  geom_point() +
  # geom_smooth(method="glm", formula=my.formula) +
  # stat_poly_eq(formula = my.formula, 
  #              aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
  #              size = rel(6),
  #              parse = TRUE) +  
  labs(title=paste("Differential Genes vs Quality Imbalance - Regress All"), 
       subtitle = paste0("Datasets selection: diff genes\u2265", dataset_n_deg_cutoff,
                        " and samples\u2265", samples.cutoff,
                        "\nData size: ", "small<", n_samples.medium, " | big>", n_samples.big, " samples"
       ),
       y="Differential genes", 
       x="Quality Imbalance") +
  theme_minimal(base_size = 20)
ggsave(filename=path, plot=plot, width=12, height=10)

path <- file.path(out.dir.path, "deg.plow.subsets.color.pdf")
plot <- datasets.stats %>%
  filter(n_degs>=dataset_n_deg_cutoff) %>%
  filter(n_samples>=samples.cutoff) %>%
  select(dataset, p_group_cor, selected, n_degs, n_samples) %>% 
  #mutate(data=as.factor(if_else(n_samples>n_samples.big, "big data", if_else(n_samples>=n_samples.medium, "medium data", "small data")))) %>% 
  mutate(data_origin = str_extract(dataset, "GSE[0-9]*")) %>%
  rename(samples=n_samples) %>% 
  ggplot(aes(y=n_degs, x=p_group_cor, color=data_origin)) +
  geom_point() +
  labs(title=paste("Differential Genes vs Quality Imbalance"), 
       subtitle = paste0("Datasets selection: diff genes\u2265", dataset_n_deg_cutoff,
                        " and samples\u2265", samples.cutoff,
                        "\nData size: ", "small<", n_samples.medium, " | big>", n_samples.big, " samples"
       ),
       y="Differential genes", 
       x="Quality Imbalance") +
  theme_minimal(base_size = 20)
ggsave(filename=path, plot=plot, width=12, height=10)

path <- file.path(out.dir.path, "deg.plow.subsets.regr.color.pdf")
plot <- datasets.stats %>%
  filter(n_degs>=dataset_n_deg_cutoff) %>%
  filter(n_samples>=samples.cutoff) %>%
  select(dataset, p_group_cor, selected, n_degs, n_samples) %>% 
  #mutate(data=as.factor(if_else(n_samples>n_samples.big, "big data", if_else(n_samples>=n_samples.medium, "medium data", "small data")))) %>% 
  mutate(data_origin = str_extract(dataset, "GSE[0-9]*")) %>%
  rename(samples=n_samples) %>% 
  ggplot(aes(y=n_degs, x=p_group_cor, color=data_origin)) +
  geom_point() +
  geom_smooth(method="glm", formula=my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               size = rel(6),
               parse = TRUE) +  
  labs(
      #  title=paste("Differential Genes vs Quality Imbalance"), 
      #  subtitle = paste0("Datasets selection: diff genes\u2265", dataset_n_deg_cutoff,
      #                   " and samples\u2265", samples.cutoff,
      #                   "\nData size: ", "small<", n_samples.medium, " | big>", n_samples.big, " samples"
      #  ),
       y="Differential genes", 
       x="Quality Imbalance") +
  theme_minimal(base_size = 20)
ggsave(filename=path, plot=plot, width=12, height=10)

path <- file.path(out.dir.path, "figure-2-deg.plow.subsets.regr.facet.pdf")
plot <- datasets.stats %>%
  # filter(n_degs>=dataset_n_deg_cutoff) %>%
  # filter(n_samples>=samples.cutoff) %>%
  select(dataset, p_group_cor, selected, n_degs, n_samples) %>% 
  #mutate(data=as.factor(if_else(n_samples>n_samples.big, "big data", if_else(n_samples>=n_samples.medium, "medium data", "small data")))) %>% 
  mutate(data_origin = str_extract(dataset, "GSE[0-9]*")) %>%
  rename(samples=n_samples) %>% 
  ggplot(aes(y=n_degs, x=p_group_cor)) +
  geom_point(color="#295D8A") +
  geom_smooth(method="glm", formula=my.formula, color="#295D8A") +
  stat_poly_eq(formula = my.formula, 
               aes(label = ..eq.label..), 
               size = rel(6),
               parse = TRUE,
               color="#295D8A") +
  stat_poly_eq(formula = my.formula, 
               aes(label = ..rr.label..), 
               size = rel(6),
               parse = TRUE,
               vjust= 2.2,
               color="#295D8A") +  
  labs(
       # title=paste("Differential Genes vs Quality Imbalance"), 
       # subtitle = expression('Datasets selection: Subsets stratified by P'[low]),
       y="Differential genes", 
       x="Quality Imbalance") +
  facet_wrap(~ data_origin, scales="free_y") +
  theme_minimal(base_size = 22) +
  theme(aspect.ratio = 1)
ggsave(filename=path, plot=plot, width=15, height=7)
