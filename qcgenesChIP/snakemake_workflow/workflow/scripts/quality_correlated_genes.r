#!/usr/bin/env Rscript
# ==============================================================================
# Collect and analyse genes correlated with quality in different datasets 
# Output directory: output/qualityCorGenes/
# ==============================================================================

# ==============================================================================
#TODO

# Examples of fgsea at https://stephenturner.github.io/deseq-to-fgsea/

# LIBRARIES AND SYSTEM OPTIONS
# ------------------------------------------------------------------------------
suppressMessages(library(fgsea, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(yaml, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
options(readr.num_columns = 0)

# GLOBAL VARIABLES  
# ------------------------------------------------------------------------------
config.path <- file.path(".", "config", "config.yaml")
config.data <- yaml.load_file(config.path)
msigdb.kegg.path <- config.data$MSIGDB_PATH
msigdb.posi.path <- config.data$MSIGDB_POSITIONAL_PATH
sra.path    <- file.path(".", "config", "metadata", "Mega_SraRunTable.csv")
datasets.path <- file.path(".", "config", "metadata", "Datasets.csv")
out.dir.path <- file.path("output", "qualityCorGenes")
gsea.table.cutoff <- config.data$GSEA_TABLE_PADJ_CUTOFF
gsea.plot.cutoff <- config.data$GSEA_PLOT_PADJ_CUTOFF2
dataset_p_group_cor_cutoff <- config.data$MAX_DATASET_P_VS_GROUP_CORRELATION
dataset_n_deg_cutoff <- config.data$MIN_DATASET_N_DEG

dir.create(out.dir.path, showWarnings=FALSE, recursive=TRUE)

# GLOBAL DATA
# ------------------------------------------------------------------------------
sra.table <- read_csv(sra.path) %>% 
  rename(sample=Run)
dataids.list <- unique(sra.table$GEO_Series)

datasets.table <- read_csv(datasets.path)

# datasets deg
deg.stats <- c()
for(DATAID in dataids.list){
  # print(paste0("quality_correlated_gemes.r: dataset ", DATAID))
  cotrl <- (datasets.table %>% filter(GEO_Series==DATAID))$Control
  treat <- (datasets.table %>% filter(GEO_Series==DATAID))$Treat
  file <- file.path("output", "RASflowResults", DATAID, "trans", "dea", "DEA", 
                    "gene-level", paste0("deg_", cotrl, "_", treat, ".tsv") )
  if (file.exists(file)){
    deg.stats <- append(deg.stats, length(readLines(file)))
  }
  else{
    deg.stats <- append(deg.stats, -1)
  }
}
deg.stats.table <- as_tibble(cbind(dataset=dataids.list, deg=deg.stats)) %>% 
  mutate(deg=as.numeric(deg)-1)


# datasets stats
datasets.stats <- data_frame()
for(DATAID in dataids.list){
  file <- file.path("output", "qualityVsExp", DATAID, paste0(DATAID, ".cor.tsv"))
  if (file.exists(file)){
    datasets.stats <- bind_rows(list(datasets.stats, read_tsv(file)))
  }
}
datasets.stats <- datasets.stats %>% 
  left_join(deg.stats.table, by="dataset") %>% 
  mutate(selected=!(deg<dataset_n_deg_cutoff | p_group_cor>dataset_p_group_cor_cutoff))

path <- file.path(out.dir.path, paste0("datasets.stats.tsv"))
write_tsv(datasets.stats, path)  



dataids.list <- (datasets.stats %>% filter(selected==TRUE))$dataset

# POSITIVELY CORRELATED GENES
  # ==============================================================================
d <- c()
for(DATAID in dataids.list){
  # print(paste0("quality_correlated_gemes.r: dataset ", DATAID))
  file <- file.path("output", "qualityVsExp", DATAID, paste0(DATAID, ".cor_pos_genes.tsv"))
  if (file.exists(file)){
    d <- append(d, read_tsv(file)$SYMBOL)
  }
}
pos_genes <- tibble(genes=d) %>% 
  group_by(genes) %>% 
  summarise(n=n()) %>% 
  filter(n>=2 & genes!="NA") %>% 
  arrange(desc(n))
out.file1 <- file.path(out.dir.path, "pos.cor.genes.tsv")
write_tsv(pos_genes, out.file1)

path <- file.path(out.dir.path, paste0("pos.cor.genes.png"))
a <- pos_genes %>% 
  rename(gene=genes, datasets=n) %>% 
  top_n(20,datasets) %>% 
  arrange(datasets) %>% 
  ggplot(aes(x=reorder(gene, datasets), y=datasets)) +
  geom_bar(stat="identity") +
  labs(title=paste("Top positively correlated genes in ", length(dataids.list), "datasets"), x="Gene") +
  coord_flip() +
  theme_minimal() +
  ggsave(path, device='png', width=15.2, height=12.8)

# NEGATIVELY CORRELATED GENES
# ==============================================================================
d <- c()
for(DATAID in dataids.list){
  file <- file.path("output", "qualityVsExp", DATAID, paste0(DATAID, ".cor_neg_genes.tsv"))
  if (file.exists(file)){
    d <- append(d, read_tsv(file)$SYMBOL)
  }
}
neg_genes <- tibble(genes=d) %>% 
  group_by(genes) %>% 
  summarise(n=n()) %>% 
  filter(n>=2 & genes!="NA") %>% 
  arrange(desc(n))
out.file2 <- file.path(out.dir.path, "neg.cor.genes.tsv")
write_tsv(neg_genes, out.file2)

path <- file.path(out.dir.path, paste0("neg.cor.genes.png"))
b <- neg_genes %>% 
  rename(gene=genes, datasets=n) %>% 
  top_n(20,datasets) %>% 
  arrange(datasets) %>% 
  ggplot(aes(x=reorder(gene, datasets), y=datasets)) +
  geom_bar(stat="identity") +
  labs(title=paste("Top negatively correlated genes in ", length(dataids.list), "datasets"), x="Gene") +
  coord_flip() +
  theme_minimal() +
  ggsave(path, device='png', width=15.2, height=12.8)

# GENE SET ENRICHMENT ANALYSIS
# ==============================================================================
# ------------------------------------------------------------------------------
# FUNCTION
# ------------------------------------------------------------------------------
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
  path1 <- file.path(out.dir, paste0(out.file.prefix, ".gsea.plot.png"))
  path2 <- file.path(out.dir, paste0(out.file.prefix, ".gsea.genes.annot.tsv"))
  
  ranks <- deframe(gene_rank_table)
  pathways <- gmtPathways(genesets)
  fgseaRes <- fgsea(pathways=pathways, stats=ranks, nperm=1000)
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
    filter(padj<table.cutoff)

  if(dim(fgseaResTidy)[1]>0){
    
    write_tsv(fgseaResTidy, path0)
    
    ggplot(fgseaResTidy %>% top_n(top_n_to_plot, desc(abs(NES))), aes(reorder(pathway, NES), NES)) +
      geom_col(aes(fill=padj<plot.cutoff)) +
      coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score",
           title=paste("Top pathways from GSEA")) + 
  theme_minimal() +
      ggsave(path1, width=15.2, height=12.8)
    
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
# ------------------------------------------------------------------------------

gsea_analysis(pos_genes, "pos.kegg", out.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.kegg.path)
gsea_analysis(neg_genes, "neg.kegg", out.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.kegg.path)
gsea_analysis(pos_genes, "pos.posi", out.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.posi.path)
gsea_analysis(neg_genes, "neg.posi", out.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.posi.path)


# AGGREGATED STATISTICS
# ==============================================================================

d <- datasets.stats %>% filter(selected==TRUE)

path <- file.path(out.dir.path, "cor.size.png")
plot <- d %>% select(dataset, p_bases_cor, p_bases_cor_ctrl, 
                    p_bases_cor_treat, p_bytes_cor, p_bytes_cor_ctrl, 
                    p_bytes_cor_treat) %>% 
  gather(p_bases_cor, p_bases_cor_ctrl, 
       p_bases_cor_treat, p_bytes_cor, p_bytes_cor_ctrl, 
       p_bytes_cor_treat, key="Feature", value="value") %>% 
  ggplot(aes(x=reorder(dataset, value), y=value))+
  geom_bar(stat='identity')+
  facet_wrap(~Feature) +
  coord_flip() +
  labs(y="Correlation coefficient (Low-Quality Probability vs Total Bases or Total Bytes)",
      title=paste("Low-Quality Probability vs Dataset Size")) + 
  theme_minimal() +
  ggsave(path, width=15.2, height=12.8)

path <- file.path(out.dir.path, "cor.group.png")
plot <- d %>% select(dataset, p_group_cor) %>% 
  ggplot(aes(x=reorder(dataset, p_group_cor), y=p_group_cor))+
  geom_bar(stat='identity')+
  coord_flip() +
  labs(y="Correlation coefficient (Low-Quality Probability vs Sample Group)",
       title=paste("Low-Quality Probability vs Sample Group")) +
  theme_minimal() +
  ggsave(path, width=15.2, height=12.8)
  
path <- file.path(out.dir.path, "pca.dunn.png")
plot <- d %>% select(dataset, dunn_all, dunn_noout) %>%  
  mutate(diff_dunn=dunn_noout-dunn_all) %>%  
  ggplot(aes(x=reorder(dataset, diff_dunn), y=diff_dunn))+
  geom_bar(stat='identity')+
  labs(y="Difference in Dunn Index (dataset without outliers vs with all samples)",
       x="",
       title=paste("PCA Clustering Evaluation Boost (Dunn Index)")) + 
  coord_flip() +
  theme_minimal() +
  ggsave(path, width=15.2, height=12.8)

path <- file.path(out.dir.path, "pca.pgamma.png")
plot <- d %>% select(dataset, pgamma_all, pgamma_noout) %>% 
  mutate(diff_pgamma=pgamma_noout-pgamma_all) %>% 
  ggplot(aes(x=reorder(dataset, diff_pgamma), y=diff_pgamma))+
  geom_bar(stat='identity')+
  labs(y="Difference in Normalized Gamma (dataset without outliers vs with all samples)",
       x="",
       title=paste("PCA Clustering Evaluation Boost (Normalized Gamma)")) + 
  coord_flip() +
  theme_minimal() +
  ggsave(path, width=15.2, height=12.8)

  # # Difference with Dunn and difference with pGamma
  # d %>% select(dataset, dunn_all, dunn_noout, pgamma_all, pgamma_noout) %>% 
  #   mutate(diff_dunn=dunn_noout-dunn_all) %>% 
  #   mutate(diff_pgamma=pgamma_noout-pgamma_all) %>% 
  #   select(dataset, diff_dunn, diff_pgamma) %>% 
  #   gather(diff_dunn, diff_pgamma, key="feature", value="value") %>% 
  #   ggplot(aes(x=reorder(dataset, value), y=value))+
  #   geom_bar(stat = "identity")+
  #   coord_flip() +
  #   facet_wrap(~feature) +
  #   theme_minimal()
  # 