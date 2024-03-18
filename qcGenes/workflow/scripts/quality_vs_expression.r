# ==============================================================================
# ANALZE 1 RNA-SEQ DATASET
# ==============================================================================

# ______________________________________________________________________________
# LIBRARIES AND SYSTEM OPTIONS
# ______________________________________________________________________________
# install.packages(c("FactoMineR", "factoextra", "fpc"))
# BiocManager::install(c("fgsea", "fpc"), ask = F, update = FALSE)
#suppressMessages(library(EnsDb.Hsapiens.v86, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(yaml, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(FactoMineR, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(factoextra, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(ggrepel, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(fpc, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(fgsea, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(rstatix, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
options(readr.num_columns = 0) # to mute readr functions
# library(ImpactEffectsize)


# ______________________________________________________________________________
# ARGUMENTS
# ______________________________________________________________________________
args = commandArgs(trailingOnly=TRUE)
if (length(args)<1 | length(args)>6) {
  # stop("Please provide GEO_Series ID as command line parameter, and optionally a subdirectory name", call.=FALSE)
  stop("Wrong number of parameters", call.=FALSE)
}

# test.dataid <- "GSE155800"
# args <- c(
#   test.dataid,
#   "main",
#   file.path("output", "main", "scores", paste(test.dataid, "scores", "txt", sep=".")),
#   file.path("output", "main", "gene_expression", test.dataid, paste(test.dataid, "samples", "rlg", "tsv", "gz", sep=".")),
#   file.path("output", "main", "gene_expression", test.dataid, paste(test.dataid, "diff", "genes", "tsv", "gz", sep=".")),
#   file.path("output", "main", "qualityVsExp", test.dataid)
# )

dataid <- args[1]
subdir <- args[2]
scores.path <- args[3]
exp.data.paths <- args[4]
dea.path <- args[5]
out.dir.path <- args[6]
dir.create(out.dir.path, showWarnings=FALSE, recursive=TRUE)

if(is.na(subdir)){
  subdir <- ""
}


# ______________________________________________________________________________
# GLOBAL VARIABLES
# ______________________________________________________________________________

config.path                   <- file.path(".", "config", "config.yaml")
config                        <- yaml.load_file(config.path)

QI_METHOD                    <- config$QI_METHOD
QI_TEST_WITH_PERMUTATIONS    <- config$QI_TEST_WITH_PERMUTATIONS
N_PERMUTATIONS               <- config$N_PERMUTATIONS

gs2d.path                     <- config$GS2D_PATH
datasets.metadata.path        <- config$DATASETS_FILE
sra.path                      <- config$SAMPLES_FILE
pval_floor                    <- 10e-10  # lowest p-value to display boxplots
outlier_sample_cutoff_p       <- 0.01     # min low-quality p-value to define an outlier sample
outlier_sample_cutoff_rank    <- 2       # max sample rank to consider outlier
outlier_sample_cutoff_n_group <- 2       # min number of samples to keep per group
outlier_sample_min_n          <- 2       # min number of samples to remove
outlier_sample_min_n_group    <- 2       # min number of samples to remove per group
outlier_sample_max_n_group    <- 2       # max number of samples to remove per group
cutoff_quant                  <- 0.25    # e.g. 0.05 to remove 5% lowest and also 5% highest correlatied genes
cutoff_absol                  <- 0.4     # e.g. 0.7 to remove genes with correlation < -0.7 and also >0.7

msigdb.kegg.path              <- config$MSIGDB_CURATED_PATH
msigdb.posi.path              <- config$MSIGDB_POSITIONAL_PATH
msigdb.gs2d.path              <- config$GS2D_MSIGDB_PATH
gsea.table.cutoff             <- config$GSEA_TABLE_PADJ_CUTOFF
gsea.plot.cutoff              <- config$GSEA_PLOT_PADJ_CUTOFF2

gsea.input_genes           <- config$GSEA_INPUT_GENES
gsea.pathway.perc.cutoff   <- config$GSEA_PATHWAY_REPRESENTATION_PERCENTAGE_CUTOFF
gsea.input.perc.cutoff     <- config$GSEA_INPUTGENES_REPRESENTATION_PERCENTAGE_CUTOFF
gsea.table.cutoff          <- config$GSEA_TABLE_PADJ_CUTOFF
gsea.plot.cutoff           <- config$GSEA_PLOT_PADJ_CUTOFF2

msigdb.hallmark.path   <- config$MSIGDB_HALLMARK_PATH
msigdb.posi.path       <- config$MSIGDB_POSITIONAL_PATH
msigdb.curated.path    <- config$MSIGDB_CURATED_PATH
msigdb.curated.cp.path <- config$MSIGDB_CURATED_CP_PATH
msigdb.regulatory.path <- config$MSIGDB_REGULATORY_PATH
msigdb.cells.path      <- config$MSIGDB_CELLS_PATH
msigdb.gs2d.path       <- config$GS2D_MSIGDB_PATH

fgsea.dir.path <- file.path(out.dir.path, "fgsea")
dir.create(fgsea.dir.path, recursive = T, showWarnings = F)


# ______________________________________________________________________________
# GENE IDS
# ______________________________________________________________________________
gene_ids_file <- config$ENSID_PATH
ensembl_ids <- read_tsv(gene_ids_file, show_col_types = FALSE)

full_gene_ids_file <- config$ENSIDFULLVER_PATH
ensembl_ids_full <- read_tsv(full_gene_ids_file, show_col_types = FALSE)


# DESEQ2 FILES
# ______________________________________________________________________________
# DIFFERENTIAL GENES TABLE
dea.table <- read_tsv(dea.path, show_col_types = FALSE) %>% 
  filter(!is.na(padj))
n.degs <- sum(dea.table$padj<=0.05 & abs(dea.table$log2FoldChange)>=1)
n.degs.FdrOnly <- sum(dea.table$padj<=0.05)

# EXPRESSION DATA TABLE
data0 <- read_tsv(exp.data.paths, show_col_types = FALSE) %>% 
  filter(!is.na(geneid))

# FILTERING
data0$var <- apply(data0 %>% dplyr::select(-geneid), 1, var)
var.min.cutoff <- quantile(data0$var[data0$var>0], 0.05)[1]
my.annot <- ensembl_ids_full %>% 
  rename(geneid=gene_id_ver) %>% 
  select(geneid, gene_name) %>% 
  distinct()
data0 <- data0 %>%
  left_join(my.annot, by="geneid") %>% 
  filter(var>var.min.cutoff) %>%
  group_by(gene_name) %>%
  arrange(desc(var)) %>%
  filter(row_number()==1) %>% 
  ungroup() %>% 
  select(-gene_name)


deseq2.sample.names <- names(data0)[!names(data0) %in% c("geneid", "var")]


# GLOBAL DATA
# ______________________________________________________________________________
gs2d.table <- read_tsv(gs2d.path, show_col_types = FALSE) # "GeneSet 2 Disease" data

# vector of disease MeSH terms
disease.meshs <- as.character(
  read_csv(datasets.metadata.path, show_col_types = FALSE) %>% 
    dplyr::filter(GEO_Series==dataid) %>% 
    dplyr::rename(mesh_terms=`mesh_terms (pipe-separated list)`) %>% 
    dplyr::select(mesh_terms))
disease.meshs <- unlist(strsplit(disease.meshs, "|", fixed = TRUE)) 

sra.table <- read_csv(sra.path, show_col_types = FALSE) %>% 
  rename(sample=Run) %>% 
  filter(GEO_Series==dataid) %>% 
  filter(Selected==1) %>% 
  filter(sample %in% deseq2.sample.names) # FILTER HERE TO FIT DESEQ2 SAMPLES => OK

datasets.table <- read_csv(datasets.metadata.path, show_col_types = FALSE)


# DATASET-SPECIFIC DATA
# ______________________________________________________________________________
metadata.table  <- sra.table

dataset.metadata <-  as.list(datasets.table %>% dplyr::filter(GEO_Series==dataid))
control.name <- dataset.metadata$Control
treats.name <- dataset.metadata$Treat

# QUALITY SCORES
scores.table <- read_tsv(scores.path, col_names=c("id", "p", "desc"), show_col_types = FALSE) %>% 
  dplyr::select(-desc) %>% 
  left_join(metadata.table, by=c("id" = "sample")) %>% 
  filter(!is.na(group)) %>%   # FILTER HERE TO FIT DESEQ2 SAMPLES => OK
  mutate(row_ungrouped = row_number(desc(p))) %>% 
  group_by(group) %>% 
  mutate(n = n()) %>% 
  mutate(rank = dense_rank(desc(p))) %>% 
  mutate(row = row_number(desc(p))) %>% 
  mutate(outlier = ifelse(p>outlier_sample_cutoff_p & rank<=outlier_sample_cutoff_rank & n-row>=outlier_sample_cutoff_n_group, 1, 0)) %>%
  mutate(outlier = ifelse(row_ungrouped<=outlier_sample_min_n, 1, outlier)) %>%
  mutate(outlier = ifelse(row<=outlier_sample_min_n_group, 1, outlier)) %>%
  mutate(outlier = ifelse(row>outlier_sample_max_n_group, 0, outlier)) %>%
  mutate(outlier = ifelse(n<=3 & row>1, 0, outlier)) %>%
  mutate(group_factor=as.factor(group)) %>% 
  mutate(group_numeric=ifelse(group == control.name, 1, 2)) %>%
  # mutate(group_numeric=as.numeric(as.factor(group_factor))) %>%
  select(-group_factor) %>%
  ungroup()
# scores.table$group_numeric


# ______________________________________________________________________________
# ______________________________________________________________________________
# QI INDEX
# ______________________________________________________________________________
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
    sqrt( (gmd(x)^2 + gmd(x)^2) /  2 )
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
# ctdiff.permut.test(scores.table$p, scores.table$group_numeric, myfunc=mean)
# ctdiff.permut.test(scores.table$p, scores.table$group_numeric, myfunc=median)

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

# ______________________________________________________________________________
# ______________________________________________________________________________

dataset.p.group.cor <- -1
dataset.p.group.test <- -1
group_1_values <- scores.table[scores.table$group_numeric==1,]$p
group_2_values <- scores.table[scores.table$group_numeric==2,]$p

if(QI_METHOD %in% c("pearson", "spearman")){
  dataset.p.group.cor <- abs(cor(scores.table$p, scores.table$group_numeric, method=QI_METHOD))
  dataset.p.group.test <- cor.test(scores.table$p, scores.table$group_numeric, method=QI_METHOD, alternative = "two.sided", exact = F)$p.value
} else if(QI_METHOD=="cohen"){
  dataset.p.group.cor <- abs(cohens_d( scores.table %>% arrange(group_numeric, subject), p~group_numeric )$effsize)
  dataset.p.group.test <- t.test(group_1_values, group_2_values, alternative = "two.sided")$p.value
} else if(QI_METHOD=="ctdiff"){
  # Lötsch J, Ultsch A. A non-parametric effect-size measure capturing changes in central tendency and data distribution shape. PLoS One. 2020 Sep 24;15(9):e0239623. doi: 10.1371/journal.pone.0239623. PMID: 32970758; PMCID: PMC7514071.
  dataset.p.group.cor <- abs(ctdiff(group_1_values, group_2_values))
  dataset.p.group.test <- wilcox.test(group_1_values, group_2_values, alternative = "two.sided", exact = F)$p.value
} else if(QI_METHOD=="permut"){
  dataset.p.group.test <- permut.test.myfunc(scores.table$p, scores.table$group_numeric, myfunc=median)
  dataset.p.group.cor <- 1 - dataset.p.group.test
}

# else if(QI_METHOD=="impact"){
#   # Lötsch J, Ultsch A. A non-parametric effect-size measure capturing changes in central tendency and data distribution shape. PLoS One. 2020 Sep 24;15(9):e0239623. doi: 10.1371/journal.pone.0239623. PMID: 32970758; PMCID: PMC7514071.
#   dataset.p.group.cor <- Impact(scores.table$p, as.factor(scores.table$group_numeric))$Impact
#   dataset.p.group.test <- wilcox.test(group_1_values, group_2_values, alternative = "two.sided", exact = F)$p.value
# }

# ______________________________________________________________________________
# ______________________________________________________________________________
# QI SIGNIFICANCE TEST
# ______________________________________________________________________________
# ______________________________________________________________________________

if(QI_TEST_WITH_PERMUTATIONS){
  if(QI_METHOD == "ctdiff"){
    dataset.p.group.test <- ctdiff.permut.test(scores.table$p, scores.table$group_numeric, permutations=N_PERMUTATIONS)    
  } else if(QI_METHOD %in% c("pearson", "cohen")){
    dataset.p.group.test <- permut.test.myfunc(scores.table$p, scores.table$group_numeric, permutations=N_PERMUTATIONS, myfunc=mean)
  } else{
    dataset.p.group.test <- permut.test.myfunc(scores.table$p, scores.table$group_numeric, permutations=N_PERMUTATIONS, myfunc=median)
  }
}


# The variable QI_METHOD will be used below only with the correlation function cor()
QI_METHOD <- ifelse(QI_METHOD %in% c("pearson", "cohen", "ctdiff"), "pearson", "spearman")
# QI_METHOD <- "pearson"


# ______________________________________________________________________________
# ______________________________________________________________________________

# dataset.p.group.cor <- abs(cor(scores.table$p, scores.table$group_numeric, method=QI_METHOD))
dataset.p.group.cor.batched <- (cluster.stats(dist(scores.table$p), scores.table$group_numeric)$pearsongamma+1)/2
dataset.p.bases.cor <- cor(scores.table$p, scores.table$Bases, method=QI_METHOD)
dataset.p.bases.ctrol.cor <- cor(scores.table[scores.table$group_numeric==0, ]$p, scores.table[scores.table$group_numeric==0, ]$Bases, method=QI_METHOD)
dataset.p.bases.treat.cor <- cor(scores.table[scores.table$group_numeric==1, ]$p, scores.table[scores.table$group_numeric==1, ]$Bases, method=QI_METHOD)
dataset.p.bytes.cor <- cor(scores.table$p, scores.table$Bytes, method=QI_METHOD)
dataset.p.bytes.ctrol.cor <- cor(scores.table[scores.table$group_numeric==0, ]$p, scores.table[scores.table$group_numeric==0, ]$Bytes, method=QI_METHOD)
dataset.p.bytes.treat.cor <- cor(scores.table[scores.table$group_numeric==1, ]$p, scores.table[scores.table$group_numeric==1, ]$Bytes, method=QI_METHOD)

# OUTLIERS AND COUNTS
outliers.table     <- scores.table %>% dplyr::filter(outlier==1) %>% dplyr::select(id) # 1-col table of outlier dataset ids
n.outliers         <- length(outliers.table$id)           # number of outlier datasets
scores.noout.table <- scores.table %>% dplyr::filter(outlier==0) # scores' table with no outlier datasets



# FIRST DATA VERSION
d0 <- data0 %>%
   dplyr::select(geneid, scores.table$id)
groups.labels0 <- scores.table$group
group.name0 <- unique(groups.labels0)


# EXPRESSION DATA WITHOUT OUTLIERS
d <- data0 %>%
  dplyr::select(geneid, scores.noout.table$id)
groups.labels <- scores.noout.table$group
group.name <- unique(groups.labels)

# CORRELATIONS AND DISTANCES
suppressWarnings(
  #  cor_coef <- apply(d[,scores.noout.table$id], 1, 
  cor_coef <- apply(data0[,scores.table$id], 1, 
                    function(x){ 
                      cor(x, 
                          scores.table$p, 
                          use="pairwise.complete.obs", 
                          method=QI_METHOD
                      ) 
                    } )
)
suppressWarnings(
  cor_coef_noout <- apply(data0[,scores.noout.table$id], 1, 
                          function(x){ 
                            cor(x, 
                                scores.noout.table$p, 
                                use="pairwise.complete.obs", 
                                method=QI_METHOD
                            ) 
                          } )
)
#d <- d %>% mutate(cor=cor_coef)
d$cor <- cor_coef

# SCORES PLOTS
# ______________________________________________________________________________
scores.bar.plot <- scores.table %>% 
  ggplot(aes(x=reorder(id, p), y=p, fill=group)) + 
  geom_bar(stat="identity") + 
  labs(title=paste0("Low Quality Probabilities in ", dataid),
       subtitle=paste0("Correlation coefficient of P_low vs Group: ", round(dataset.p.group.cor, 4)),
       x="Samples",
       y="Low-Quality Probability") +
  coord_flip() +
  theme_minimal()

scores_bytes.scatter.plot <- scores.table %>%  
  ggplot(aes(x=p, y=Bytes, color=group)) + 
  geom_point(size = 8) + 
  geom_smooth(method='lm', formula='y ~ x') +
  labs(title=paste0("Quality vs Bytes in ", dataid),
       x="Low-Quality Probability") +
  theme_minimal()

scores_bases.scatter.plot <- scores.table %>% 
  ggplot(aes(x=p, y=Bases, color=group)) + 
  geom_point(size = 8) + 
  geom_smooth(method='lm', formula='y ~ x') +
  labs(title=paste0("Quality vs Bases in ", dataid),
       x="Low-Quality Probability") +
  theme_minimal()

# capture.output(suppressWarnings(suppressMessages(
  ggexport(
    plotlist = list(
      scores.bar.plot, 
      scores_bytes.scatter.plot, 
      scores_bases.scatter.plot), 
    nrow = 1, ncol = 3, width=6000, height=3000, res=300, verbose=F,
    filename = file.path(out.dir.path, paste0(dataid, ".scores.png")))
# )), file='/dev/null')

  
# PCA PLOTS
# ______________________________________________________________________________
palette <- c("#999999", "#377EB8")
if(length(group.name)==3){
  palette <- c("#999999", "#377EB8", "#000000")
}

# Function returning plots from a PCA {FactoMineR} object
# res.pca (PCA): object returned by funtion PCA() of FactoMineR package
# plot.title (String): PCA plot title
# control.name (String): name of control group (must fit with a group name in scores.table)
# scores.table (data.frame): describes samples (column id for sample id String, p for low-quality probability between 0 and 1, group for sample group name String, outlier for identifying outliers by 1 or otherwise 0)
# label_only_outliers (boolean): TRUE to label only outliers, FALSE (default) to label all samples
# RETURNS A LIST OF 4 GGPLOTS AND A LIST
# ggplot 1: scree plot
# ggplot 2: top variables contributions in dimension 1 plot
# ggplot 3: top variables contributions in dimension 2 plot
# ggplot 4: PCA plot
# list of PCA validation metrics (list) as returned by function cluster.stats() from package fpc
PCA_plots <- function(res.pca, plot.title, control.name, scores.table, label_only_outliers=FALSE, design_bias="", degs=-1){
  
  var <- get_pca_var(res.pca)
  # Diagnostic plots
  scree.plot0 <- fviz_eig(res.pca, addlabels = TRUE, rotate=T, subtitle=plot.title)
  contrib1.plot0 <- fviz_contrib(res.pca, choice = "var", axes = 1, top = 20, rotate=T, subtitle=plot.title)
  contrib2.plot0 <- fviz_contrib(res.pca, choice = "var", axes = 2, top = 20, rotate=T, subtitle=plot.title)
  # Integration of the scores table
  res.pca.enriched.table <-  as_tibble(res.pca$ind$coord, rownames="id") %>% 
      left_join(scores.table, by="id") %>% 
      mutate(outlier=ifelse(outlier==1, TRUE, FALSE)) %>% 
      mutate(log_p=round(-log10(1-p), digits = 3)) %>% 
      mutate(labels=ifelse(outlier==TRUE, paste(id, p, sep="\n"), "")) %>% 
      column_to_rownames(var="id")

  res.pca.enriched.dist <- dist(res.pca.enriched.table[, c("Dim.1", "Dim.2")])
#  clustering.labels <- ifelse(res.pca.enriched.table$group==control.name, 1, 2)
  clustering.labels <- res.pca.enriched.table$group_numeric
  pca.valid <- cluster.stats(res.pca.enriched.dist, clustering.labels)
  # pca plot
  pca.plot0.data <- as_tibble(res.pca$ind$coord, rownames="id") %>% 
    left_join(scores.table, by="id") %>% 
    mutate(outlier=ifelse(outlier==1, TRUE, FALSE)) %>% 
    mutate(log_p=round(-log10(1-p), digits = 3)) %>% 
    mutate(labels=paste(id, p, sep="\n"))
  
  if(label_only_outliers==TRUE){
    pca.plot0.data <- pca.plot0.data %>% 
      mutate(labels=ifelse(outlier==TRUE, paste(id, p, sep="\n"), ""))
  }

  if(!is.na(sum(scores.table$batch))){
    pca.plot0.data <- pca.plot0.data %>% 
      mutate(labels=paste(labels, batch, sep="\n"), "")
  }
  pca.plot0 <- pca.plot0.data %>% 
    ggplot(aes(x = Dim.1, y = Dim.2))+
    geom_point(aes(color = group, size = log_p, shape = outlier)) +
    geom_text_repel(aes(label=labels)) +
    labs(title=plot.title, 
         subtitle = paste0("DesignBias=", round(design_bias,2),
                           ", Gamma=", round(pca.valid$pearsongamma, 2),
                           ", Dunn1=", round(pca.valid$dunn, 2),
                           ", WbRatio=", round(pca.valid$wb.ratio, 2),
                           if_else(degs>0, paste0(", DEGs=", degs), "")) ) + 
#         subtitle = paste0("Dunn index: ", round(pca.valid$dunn, 2)) ) + 
    guides(size=FALSE) +
    theme_minimal() + 
    theme(panel.background = element_rect(fill = "white", colour = "grey50"))
  
  pca.valid0 <- pca.valid
  
  return(list(scree.plot0, contrib1.plot0, contrib2.plot0, pca.plot0, pca.valid0))
}


# PCA PLOT FOR ALL GENES AND ALL SAMPLES
# ______________________________________________________________________________
norm.table <- as.data.frame(d0 %>% column_to_rownames("geneid"))
#pca.table <- log2(norm.table[which(apply(norm.table, 1, var) != 0), ]+0.00000001)
pca.table <- norm.table[which(apply(norm.table, 1, var) != 0), ] # For the RLOG data
pca.table <- pca.table %>% 
  rownames_to_column(var="gene_id_ver") %>% 
  left_join(ensembl_ids_full %>% select(gene_id_ver, gene_name) %>% distinct(), by="gene_id_ver") %>% 
  mutate(rowname=paste(gene_name, gene_id_ver, sep = "_")) %>% 
  column_to_rownames() %>%
  select(-gene_id_ver, -gene_name)

pca.table.symbols2geneid <- ensembl_ids_full %>%
  rename(GENEID=gene_id_ver, SYMBOL=gene_name) %>%
  select(SYMBOL, GENEID) %>%
  unique()


# rownames(pca.table) <- paste(pca.table.symbols, pca.table.ids, sep="_")
pca.table <- t(pca.table)
res.pca <- PCA(pca.table, scale.unit = TRUE, ncp = 5, graph = FALSE)
plots_and_valid <- PCA_plots(res.pca, paste0("PCA (all genes & all samples): ", dataid), control.name, scores.table, TRUE, design_bias=dataset.p.group.cor, degs=n.degs)
scree.plot0    <- plots_and_valid[[1]]
contrib1.plot0 <- plots_and_valid[[2]]
contrib2.plot0 <- plots_and_valid[[3]]
pca.plot0      <- plots_and_valid[[4]]
pca.valid0     <- plots_and_valid[[5]]


# PCA PLOT FOR ALL GENES BUT NO OUTLIER SAMPLES
# ______________________________________________________________________________
if(ncol(d)>5){
  norm.table <- as.data.frame(d %>% dplyr::select(-cor) %>% column_to_rownames("geneid"))
#  pca.table <- log2(norm.table[which(apply(norm.table, 1, var) != 0), ]+0.00000001)
  pca.table <- norm.table[which(apply(norm.table, 1, var) != 0), ]
  pca.table <- pca.table %>% 
    rownames_to_column(var="gene_id_ver") %>% 
    left_join(ensembl_ids_full %>% select(gene_id_ver, gene_name) %>% distinct(), by="gene_id_ver") %>% 
    mutate(rowname=paste(gene_name, gene_id_ver, sep = "_")) %>% 
    column_to_rownames() %>%
    select(-gene_id_ver, -gene_name)
  
  # pca.table.ids <- rownames(pca.table)
  # #pca.table.symbols <- ensembldb::select(EnsDb.Hsapiens.v86, 
  # #                                       keys=pca.table.ids, 
  # #                                       keytype = "GENEID", 
  # #                                       columns = c("SYMBOL","GENEID"))
  # pca.table.symbols <- merge(data.frame(GENEID=pca.table.ids), 
  #                            pca.table.symbols2geneid, 
  #                            all.x=TRUE)$SYMBOL
  # rownames(pca.table) <- paste(pca.table.symbols, pca.table.ids, sep="_")
  pca.table <- t(pca.table)
  res.pca <- PCA(pca.table, scale.unit = TRUE, ncp = 5, graph = FALSE)
  plots_and_valid <- PCA_plots(res.pca, paste0("PCA (all genes, ", n.outliers, " outlier(s) removed):"), control.name, scores.table, design_bias=dataset.p.group.cor)
  scree.plot1    <- plots_and_valid[[1]]
  contrib1.plot1 <- plots_and_valid[[2]]
  contrib2.plot1 <- plots_and_valid[[3]]
  pca.plot1      <- plots_and_valid[[4]]
  pca.valid1     <- plots_and_valid[[5]]
}else{
  scree.plot1    <- ggplot()
  contrib1.plot1 <- ggplot()
  contrib2.plot1 <- ggplot()
  pca.plot1      <- ggplot()
  pca.valid1     <- ggplot()
}

# PCA PLOT OF NON CORRELATED GENES
# ______________________________________________________________________________
dtmp <- d %>% dplyr::select(-cor)
dtmp$cor <- cor_coef_noout
cutoffs.c1 <- quantile(dtmp$cor, probs=c(cutoff_quant, 1-cutoff_quant), na.rm=T)
cutoffs.c2 <- c(-cutoff_absol, cutoff_absol)
cutoffs <- c( cutoffs.c1[1], cutoffs.c1[2] )
#cutoffs <- c( cutoffs.c2[1], cutoffs.c2[2] )
#cutoffs <- c( max(cutoffs.c1[1], cutoffs.c2[1]), min(cutoffs.c1[2], cutoffs.c2[2]) )
cutoffs <- c( min(cutoffs.c1[1], cutoffs.c2[1]), max(cutoffs.c1[2], cutoffs.c2[2]) )

if(ncol(d)>5){
  norm.table0 <- as.data.frame(dtmp %>% dplyr::select(-cor) %>% column_to_rownames("geneid"))
  zerovar.genes <- apply(norm.table0, 1, var) == 0
  n.zerovar.genes <- sum(zerovar.genes)

  dcut <- dtmp[!zerovar.genes, ] %>% dplyr::filter(cor>=cutoffs[1] & cor<=cutoffs[2])
  n.pos.cor.genes <- dim(dtmp[!zerovar.genes, ] %>% dplyr::filter(cor<=cutoffs[1]))[1]
  n.neg.cor.genes <- dim(dtmp[!zerovar.genes, ] %>% dplyr::filter(cor>=cutoffs[2]))[1]
  norm.table <- as.data.frame(dcut %>% dplyr::select(-cor) %>% column_to_rownames("geneid"))
  n.filtered.genes <- nrow(dcut)
  n.filtered.out.genes <- nrow(dtmp[!zerovar.genes, ]) - nrow(dcut)
 
#  pca.table <- log2(norm.table[which(apply(norm.table, 1, var) != 0), ]+0.00000001)
  pca.table <- norm.table[which(apply(norm.table, 1, var) != 0), ]
  pca.table <- pca.table %>% 
    rownames_to_column(var="gene_id_ver") %>% 
    left_join(ensembl_ids_full %>% select(gene_id_ver, gene_name) %>% distinct(), by="gene_id_ver") %>% 
    mutate(rowname=paste(gene_name, gene_id_ver, sep = "_")) %>% 
    column_to_rownames() %>%
    select(-gene_id_ver, -gene_name)
  # pca.table.ids <- rownames(pca.table)
  # pca.table.symbols <- merge(data.frame(GENEID=pca.table.ids), 
  #                            pca.table.symbols2geneid, 
  #                           all.x=TRUE)$SYMBOL
  # rownames(pca.table) <- paste(pca.table.symbols, pca.table.ids, sep="_")
  pca.table <- t(pca.table)
  res.pca <- PCA(pca.table, scale.unit = TRUE, ncp = 5, graph = FALSE)
  plot.title <- paste0("PCA (", n.filtered.out.genes, " (",n.pos.cor.genes, "+", n.neg.cor.genes, ") genes and ", n.outliers, " outlier(s) removed): ", dataid)
  plots_and_valid <- PCA_plots(res.pca, plot.title, control.name, scores.table, design_bias=dataset.p.group.cor)
  scree.plot2    <- plots_and_valid[[1]]
  contrib1.plot2 <- plots_and_valid[[2]]
  contrib2.plot2 <- plots_and_valid[[3]]
  pca.plot2      <- plots_and_valid[[4]]
  pca.valid2     <- plots_and_valid[[5]]
}else{
  scree.plot2    <- ggplot()
  contrib1.plot2 <- ggplot()
  contrib2.plot2 <- ggplot()
  pca.plot2      <- ggplot()
  pca.valid2     <- ggplot()
}
rm(dtmp)


# PCA PLOT USING RANDOM GENE SELECTION
# ______________________________________________________________________________
if(ncol(d)>5){
  random_ids <- sample(d[!zerovar.genes, ]$geneid, n.filtered.genes)
  norm.table <- as.data.frame(d %>% 
      filter(geneid %in% random_ids) %>% 
      column_to_rownames("geneid") %>%
      dplyr::select(-cor)
  )
#  pca.table <- log2(norm.table[which(apply(norm.table, 1, var) != 0), ]+0.00000001)
  pca.table <- norm.table[which(apply(norm.table, 1, var) != 0), ]
  pca.table <- pca.table %>% 
    rownames_to_column(var="gene_id_ver") %>% 
    left_join(ensembl_ids_full %>% select(gene_id_ver, gene_name) %>% distinct(), by="gene_id_ver") %>% 
    mutate(rowname=paste(gene_name, gene_id_ver, sep = "_")) %>% 
    column_to_rownames() %>%
    select(-gene_id_ver, -gene_name)
  # #pca.table <- pca.table %>% 
  # #  rownames_to_column() %>% 
  # #  sample_n(dim(pca.table)[1] - n.filtered.genes) %>% 
  # #  column_to_rownames()
  # pca.table.ids <- rownames(pca.table)
  # #pca.table.symbols <- ensembldb::select(EnsDb.Hsapiens.v86, 
  # #                                       keys=pca.table.ids, 
  # #                                       keytype = "GENEID", 
  # #                                       columns = c("SYMBOL","GENEID"))
  # pca.table.symbols <- merge(data.frame(GENEID=pca.table.ids), 
  #                           pca.table.symbols2geneid, 
  #                           all.x=TRUE)$SYMBOL
  # rownames(pca.table) <- paste(pca.table.symbols, pca.table.ids, sep="_")
  pca.table <- t(pca.table)
  res.pca <- PCA(pca.table, scale.unit = TRUE, ncp = 5, graph = FALSE)
  plot.title <- paste0("PCA (", n.filtered.out.genes, " random genes and ", n.outliers, " outlier(s) removed): ", dataid)
  plots_and_valid <- PCA_plots(res.pca, plot.title, control.name, scores.table, design_bias=dataset.p.group.cor)
  scree.plot3    <- plots_and_valid[[1]]
  contrib1.plot3 <- plots_and_valid[[2]]
  contrib2.plot3 <- plots_and_valid[[3]]
  pca.plot3      <- plots_and_valid[[4]]
  pca.valid3     <- plots_and_valid[[5]]
}else{
  scree.plot3    <- ggplot()
  contrib1.plot3 <- ggplot()
  contrib2.plot3 <- ggplot()
  pca.plot3      <- ggplot()
  pca.valid3     <- ggplot()
}


# MEAN AND SD OF PCA CLUSTERING QUALITY USING RANDOM GENE SELECTION
# ______________________________________________________________________________
if(ncol(d)>5){

dunn1.rand.values <- c()
wbrat.rand.values <- c()
gamma.rand.values <- c()

for(i in 1:3){
  random_ids <- sample(d[!zerovar.genes, ]$geneid, n.filtered.genes)
  norm.table <- as.data.frame(d %>% 
      filter(geneid %in% random_ids) %>% 
      column_to_rownames("geneid") %>%
      dplyr::select(-cor)
  )
  pca.table <- norm.table[which(apply(norm.table, 1, var) != 0), ]
  pca.table <- pca.table %>% 
    rownames_to_column(var="gene_id_ver") %>% 
    left_join(ensembl_ids_full %>% select(gene_id_ver, gene_name) %>% distinct(), by="gene_id_ver") %>% 
    mutate(rowname=paste(gene_name, gene_id_ver, sep = "_")) %>% 
    column_to_rownames() %>%
    select(-gene_id_ver, -gene_name)
  # pca.table.ids <- rownames(pca.table)
  # pca.table.symbols <- merge(data.frame(GENEID=pca.table.ids), 
  #                           pca.table.symbols2geneid, 
  #                           all.x=TRUE)$SYMBOL
  # rownames(pca.table) <- paste(pca.table.symbols, pca.table.ids, sep="_")
  pca.table <- t(pca.table)
  res.pca <- PCA(pca.table, scale.unit = TRUE, ncp = 5, graph = FALSE)
  plot.title <- paste0("PCA (", n.filtered.out.genes, " random genes and ", n.outliers, " outlier(s) removed): ", dataid)
  plots_and_valid <- PCA_plots(res.pca, plot.title, control.name, scores.table, design_bias=dataset.p.group.cor)
  pca.valid.random     <- plots_and_valid[[5]]
  dunn1.rand.values <- c(dunn1.rand.values, pca.valid.random$dunn)
  wbrat.rand.values <- c(wbrat.rand.values, pca.valid.random$wb.ratio)
  gamma.rand.values <- c(gamma.rand.values, pca.valid.random$pearsongamma)
}
rand.clustering.table <- data.frame(dunn1.rand.values, wbrat.rand.values, gamma.rand.values)
write_tsv(rand.clustering.table, file.path(out.dir.path, paste0(dataid, '.random.clust.quality.tsv')))

}



# SAVE PCA PLOTS AND GENE LISTS
# ______________________________________________________________________________
# capture.output(suppressWarnings(suppressMessages(
  ggexport(plotlist = list(scree.plot0, scree.plot1, scree.plot2, scree.plot3, 
                           contrib1.plot0, contrib1.plot1, contrib1.plot2, contrib1.plot3, 
                           contrib2.plot0, contrib2.plot1, contrib2.plot2, contrib2.plot3), 
           nrow = 3, ncol = 4, width=7200, height=7200, res=300, verbose=F,
           filename = file.path(out.dir.path, paste0(dataid, '.pca_diag.png'))
           )
# )),  file='/dev/null')

# capture.output(suppressWarnings(suppressMessages(
  ggexport(plotlist = list(pca.plot0, pca.plot1), 
           nrow = 1, ncol = 2, width=3600, height=1800, res=300, verbose=FALSE,
           filename=file.path(out.dir.path, paste0(dataid, '.pca.png'))
  )
# )),  file='/dev/null')

# capture.output(suppressWarnings(suppressMessages(
  ggexport(plotlist = list(pca.plot2, pca.plot3), 
           nrow = 1, ncol = 2, width=4800, height=2400, res=300, verbose=FALSE,
           filename=file.path(out.dir.path, paste0(dataid, '.pca.genesout.png'))
  )
# )),  file='/dev/null')

dcut2 <- d %>% dplyr::filter(cor<=cutoffs[1] | cor>=cutoffs[2])
if(length(dcut2$geneid>0)){
#  dcut2.symbols <- ensembldb::select(EnsDb.Hsapiens.v86, keys=dcut2$geneid, keytype="GENEID", columns=c("SYMBOL","GENEID"))
  dcut2 <- dcut2 %>% left_join(pca.table.symbols2geneid, by=c("geneid" = "GENEID")) %>% 
    dplyr::select(geneid, SYMBOL, cor, setdiff(colnames(dcut2), c("geneid", "SYMBOL", "cor")) ) %>% 
    arrange(desc(cor))
}
write_tsv(dcut2 %>% dplyr::filter(cor>0), file.path(out.dir.path, paste0(dataid, '.cor_pos_genes.tsv')))
write_tsv(dcut2 %>% dplyr::filter(cor<0), file.path(out.dir.path, paste0(dataid, '.cor_neg_genes.tsv')))
d0tmp <- d0
d0tmp$cor <- cor_coef
write_tsv(d0tmp, file.path(out.dir.path, paste0(dataid, '.cor_genes.tsv')))
rm(d0tmp)



# CORRELATION VS DISEASE VS DEG
# ______________________________________________________________________________
# yaml.file <- yaml.load_file(file.path(projects.dir, dataid, "configs", "config_main.yaml"))
#dea.table <- read.table(dea.path, header = TRUE, row.names = 1)
#dea.table <- as_tibble(dea.table, rownames="geneid")
          
# Merge tables and bin correlations
tags <- c("[-1|-0.8)","[-0.8|-0.6)", "[-0.6|-0.4)", "[-0.4|-0.2)", "[-0.2|0)", 
          "[0|0.2)","[0.2|0.4)", "[0.4|0.6)","[0.6|0.8)", "[0.8|1)")

# dcor %>% print(n=20)

dcor <- d %>% 
  dplyr::select(geneid, cor) %>%
  dplyr::filter(!is.na(cor)) %>% 
  left_join(pca.table.symbols2geneid, by=c("geneid" = "GENEID")) %>% 
  rename(symbol=SYMBOL) %>% 
  left_join(gs2d.table %>% 
              dplyr::filter(name %in% disease.meshs) %>% 
              group_by(symbol) %>% 
              top_n(1, desc(fdr)),
            by="symbol") %>% 
  dplyr::select(geneid, entrez_id=gene_id, symbol, cor, 
         gs2d_fc=fold_change, gs2d_fdr=fdr, disease=name) %>% 
  left_join(dea.table, by="geneid") %>% 
  mutate(cor_bin = case_when(
    cor <  -0.8              ~ tags[1],
    cor >= -0.8 & cor < -0.6 ~ tags[2],
    cor >= -0.6 & cor < -0.4 ~ tags[3],
    cor >= -0.4 & cor < -0.2 ~ tags[4],
    cor >= -0.2 & cor <    0 ~ tags[5],
    cor >=    0 & cor <  0.2 ~ tags[6],
    cor >=  0.2 & cor <  0.4 ~ tags[7],
    cor >=  0.4 & cor <  0.6 ~ tags[8],
    cor >=  0.6 & cor <  0.8 ~ tags[9],
    cor >=  0.8              ~ tags[10]
    )) %>% 
  mutate(cor_bin=factor(cor_bin, levels = tags, ordered = FALSE)) %>% 
  mutate(disease_gene=ifelse(gs2d_fdr<=0.05 & gs2d_fc>=2, TRUE, FALSE)) %>%
  mutate(disease_gene=ifelse(is.na(disease_gene), FALSE, disease_gene)) %>% 
  # mutate(diff_gene=padj<=0.05 & abs(log2FoldChange)>=1 )
  mutate(diff_gene=if_else(padj<=0.05, TRUE, FALSE)) %>% 
  mutate(diff_gene=ifelse(is.na(diff_gene), FALSE, diff_gene))
  

# sort(dcor$padj, decreasing = T)[1:10]
# sum(is.na(dcor$geneid))
# sum(is.na(dcor$padj))
# sum(!is.na(dcor$padj))
# 
# dcor %>% 
#   filter(is.na(padj))

# sum(dcor$disease_gene)
# sum(dcor$diff_gene)

dcor$gs2d_fdr[dcor$gs2d_fdr==0] <-rnorm(1, min(dcor$gs2d_fdr[dcor$gs2d_fdr>0], na.rm=T), min(dcor$gs2d_fdr[dcor$gs2d_fdr>0], na.rm=T)/10)
dcor$padj[dcor$padj==0] <-rnorm(1, min(dcor$padj[dcor$padj>0], na.rm=T), min(dcor$padj[dcor$padj>0], na.rm=T)/10)

# Correlation coefficients between abs(gene-quality probability correlations) and FDR (either differential gene expression or disease literature association)

cor.q.fdr.allgenes <- NA
cor.q.fdr.pos.allgenes <- NA
cor.q.fdr.neg.allgenes <- NA
cor.q.fdr.difgenes <- NA
cor.q.fdr.pos.difgenes <- NA
cor.q.fdr.neg.difgenes <- NA
cor.q.fdr.disgenes <- NA
cor.q.fdr.pos.disgenes <- NA
cor.q.fdr.neg.disgenes <- NA

# Differential statistics on all genes
if(dim(dcor)[1]>0){
  cor.q.fdr.allgenes <- cor(abs(dcor$cor), dcor$padj, method=QI_METHOD, use="pairwise.complete.obs")
}
if(dim(dcor %>% dplyr::filter(cor>0))[1]>0){
  cor.q.fdr.pos.allgenes <- cor(abs((dcor %>% dplyr::filter(cor>0))$cor), (dcor %>% dplyr::filter(cor>0))$padj, method=QI_METHOD, use="pairwise.complete.obs")
}
if(dim(dcor %>% dplyr::filter(cor<0))[1]>0){
  cor.q.fdr.neg.allgenes <- cor(abs((dcor %>% dplyr::filter(cor<0))$cor), (dcor %>% dplyr::filter(cor<0))$padj, method=QI_METHOD, use="pairwise.complete.obs")
}

# Differential genes
if(dim(dcor %>% dplyr::filter(diff_gene==TRUE))[1]>0){
  cor.q.fdr.difgenes <- cor(abs((dcor %>% dplyr::filter(diff_gene==TRUE))$cor), -log10((dcor %>% dplyr::filter(diff_gene==TRUE))$padj), method=QI_METHOD, use="pairwise.complete.obs")
}
if(dim(dcor %>% dplyr::filter(cor>0 & diff_gene==TRUE))[1]>0){
  cor.q.fdr.pos.difgenes <- cor(abs((dcor %>% dplyr::filter(cor>0 & diff_gene==TRUE))$cor), -log10((dcor %>% dplyr::filter(cor>0 & diff_gene==TRUE))$padj), method=QI_METHOD, use="pairwise.complete.obs")
}
if(dim(dcor %>% dplyr::filter(cor<0 & diff_gene==TRUE))[1]>0){
  cor.q.fdr.neg.difgenes <- cor(abs((dcor %>% dplyr::filter(cor<0 & diff_gene==TRUE))$cor), -log10((dcor %>% dplyr::filter(cor<0 & diff_gene==TRUE))$padj), method=QI_METHOD, use="pairwise.complete.obs")
}

# Disease genes
if(dim(dcor %>% dplyr::filter(disease_gene==TRUE))[1]>0){
  cor.q.fdr.disgenes <- cor(abs((dcor %>% dplyr::filter(disease_gene==TRUE))$cor), (dcor %>% dplyr::filter(disease_gene==TRUE))$gs2d_fdr, method=QI_METHOD, use="pairwise.complete.obs")
}
if(dim(dcor %>% dplyr::filter(cor>0 & disease_gene==TRUE))[1]>0){
  cor.q.fdr.pos.disgenes <- cor(abs((dcor %>% dplyr::filter(cor>0 & disease_gene==TRUE))$cor), (dcor %>% dplyr::filter(cor>0 & disease_gene==TRUE))$gs2d_fdr, method=QI_METHOD, use="pairwise.complete.obs")
}
if(dim(dcor %>% dplyr::filter(cor<0 & disease_gene==TRUE))[1]>0){
  cor.q.fdr.neg.disgenes <- cor(abs((dcor %>% dplyr::filter(cor<0 & disease_gene==TRUE))$cor), (dcor %>% dplyr::filter(cor<0 & disease_gene==TRUE))$gs2d_fdr, method=QI_METHOD, use="pairwise.complete.obs")
}


# QUALITY VS ALL GENES
# ______________________________________________________________________________
# General histogram
cor.genes.plot <- dcor %>% 
  group_by(cor_bin) %>% 
  summarise(n=n()) %>% 
  ggplot(aes(x=cor_bin, y=n)) +
  geom_bar(stat="identity", colour="black", fill="white") +
  labs(title=paste0("Quality vs All Genes in ", dataid),
       subtitle=paste0("Correlation of P_low and sample groups: ", round(dataset.p.group.cor, 4)),
       x="gene-quality correlation", 
       y="Genes") +
  geom_text(aes(label=n), position=position_dodge(width=0.9), vjust=-0.25) +
  theme_minimal()   
 
# # cor vs dea fc
# cor.fc.plot <- dcor %>% 
#   filter(!is.na(log2FoldChange)) %>% 
#   ggplot(aes(cor_bin, log2FoldChange))+
#   geom_boxplot() + 
#   labs(title=paste0("Quality vs All Genes (Differential Expression Fold Change) in ", dataid), x="Correlation Coefficients", y="Log2 Fold Change") +
#   theme_minimal()   
# plot.y.max <- max(ggplot_build(cor.fc.plot)$data[[1]]$ymax)*1.5
# plot.y.min <- min(ggplot_build(cor.fc.plot)$data[[1]]$ymin)*1.5
# cor.fc.plot <- cor.fc.plot + coord_cartesian(ylim = c(plot.y.min,plot.y.max))

cor.fc.plot <- dcor %>% 
  filter(!is.na(log2FoldChange)) %>% 
  mutate(CorSign=ifelse(cor>=0,"POS", "NEG")) %>% 
  ggplot(aes(cor, log2FoldChange, fill=CorSign))+
  geom_point() + 
  stat_cor(label.x.npc = "left", label.y.npc = "top") +
  geom_smooth(method='lm', formula='y ~ x') +
  labs(title=paste0("Quality vs All Genes (Differential Expression Fold Change) in ", dataid), 
       x="gene-quality correlation", 
       y="Log2 Fold Change") +
  theme_minimal()   


# cor vs dea pval
# cor.q.fdr.allgenes <- cor(abs(dcor$cor), dcor$padj, method=QI_METHOD, use="pairwise.complete.obs")
# cor.pval.plot <- dcor %>% 
#   filter(!is.na(padj)) %>% 
#   # mutate(padj=ifelse(padj<pval_floor, pval_floor, padj)) %>% 
#   # top_frac(0.1, desc(padj)) %>%
#   ggplot(aes(cor_bin, -log10(padj))) +
#   geom_boxplot() + 
#   # labs(title=paste0("Quality vs All Genes (Best 10% Differential Expression FDR) in ", dataid), 
#   labs(title=paste0("Quality vs All Genes (Differential Expression FDR) in ", dataid), 
#        subtitle = paste("Correlation -log10(FDR) vs abs(gene-quality correlation):", round(cor.q.fdr.allgenes, 2)),
#        x="Correlation Coefficients", 
#        y="-Log10 of False Discovery Rate (FDR)") +
#   theme_minimal()   
# plot.y.max <- max(ggplot_build(cor.pval.plot)$data[[1]]$ymax)*1.5
# plot.y.min <- 0
# cor.pval.plot <- cor.pval.plot + coord_cartesian(ylim = c(plot.y.min,plot.y.max))

cor.points.plot <- dcor %>% 
  # dplyr::filter(disease_gene==TRUE) %>%
  # mutate(padj=ifelse(padj<pval_floor, pval_floor, padj)) %>%
  select(cor, padj) %>% 
  mutate(CorSign=ifelse(cor>=0,"POS", "NEG")) %>% 
  # ggplot(aes(abs(cor), -log10(padj))) +
  ggplot(aes(cor, -log10(padj), fill=CorSign)) +
  geom_point() + 
  stat_cor(label.x.npc = "left", label.y.npc = "top") +
  # stat_regline_equation(label.x = 0, label.y = max(-log10(dcor$padj)-2, na.rm=T)) +
  geom_smooth(method='lm', formula='y ~ x') +
  labs(title=paste0("Quality vs All Genes' FDR in ", dataid), 
       subtitle = paste("Correlation of -log10(FDR) vs |gene-quality correlation|:", round(cor.q.fdr.allgenes, 2)),
       x="gene-quality correlation", 
       y="-Log10 of False Discovery Rate (FDR)") +
  theme_minimal()  


# QUALITY VS DISEASE GENES
# ______________________________________________________________________________
# cor vs disease genes
cor.dis.plot <- dcor %>% 
  group_by(cor_bin) %>% 
  summarise(n=sum(disease_gene, na.rm=T)) %>% 
  ggplot(aes(x=cor_bin, y=n)) +
  geom_bar(stat="identity", colour="black", fill="white") +
  labs(title=paste0("Quality vs Disease Genes in ", dataid),
       subtitle=paste0("Disease genes: ", round(sum(dcor$disease_gene, na.rm=T), 2)), 
       x="gene-quality correlation", 
       y="Disease Genes") +
  geom_text(aes(label=n), position=position_dodge(width=0.9), vjust=-0.25) +
  theme_minimal()   

cor.dis.perc.plot <- dcor %>% 
  group_by(cor_bin) %>% 
  summarise(n=sum(disease_gene, na.rm=T), total=n(), perc=n/total, label=paste(n, total, sep = "/")) %>% 
  ggplot(aes(x=cor_bin, y=perc)) +
  geom_bar(stat="identity", colour="black", fill="white") +
  labs(title=paste0("Quality vs Disease Genes (%) in ", dataid), 
       x="gene-quality correlation", 
       y="Percentage of Disease Genes") +
  geom_text(aes(label=label), position=position_dodge(width=0.9), vjust=-0.25) +
  theme_minimal()   

# cor.q.fdr.disgenes <- cor(abs((dcor %>% dplyr::filter(disease_gene==TRUE))$cor), (dcor %>% dplyr::filter(disease_gene==TRUE))$gs2d_fdr, method=QI_METHOD, use="pairwise.complete.obs")
# missing_categories <- setdiff(tags, unique((dcor %>% dplyr::filter(disease_gene==TRUE))$cor_bin))
# missing_categ_rows <- tibble(cor_bin=missing_categories, gs2d_fdr=1, disease_gene=TRUE)
# dcor2 <- rbind(dcor %>% select(cor_bin, gs2d_fdr, disease_gene), missing_categ_rows)
# cor.dis.box.plot <- dcor2 %>%
#   filter(disease_gene==TRUE) %>% 
#   ggplot(aes(cor_bin, -log10(gs2d_fdr))) +
#   geom_boxplot() +
#   labs(title=paste0("Quality vs Disease Genes (Disease Genes FDR) in ", dataid), 
#        subtitle = paste("Correlation -log10(FDR) vs abs(gene-quality correlation):", round(cor.q.fdr.disgenes, 2)),
#        x="Correlation Coefficients", 
#        y="-Log10 of False Discovery Rate (FDR)") +
#   theme_minimal()
# suppressWarnings(
#   plot.y.max <- max(ggplot_build(cor.dis.box.plot)$data[[1]]$ymax)*1.5
# )
# plot.y.min <- 0
# cor.dis.box.plot <- cor.dis.box.plot + coord_cartesian(ylim = c(plot.y.min,plot.y.max))

cor.dis.points.plot <- dcor %>% 
  dplyr::filter(disease_gene==TRUE) %>%
  # mutate(gs2d_fdr=ifelse(gs2d_fdr<pval_floor, pval_floor, gs2d_fdr)) %>%
  select(cor, gs2d_fdr) %>% 
  mutate(CorSign=ifelse(cor>=0,"POS", "NEG")) %>% 
  # ggplot(aes(abs(cor), -log10(padj))) +
  ggplot(aes(cor, -log10(gs2d_fdr), fill=CorSign)) +
  geom_point() + 
  stat_cor(label.x.npc = "left", label.y.npc = "top") +
  # stat_regline_equation(label.x = 0.5, label.y = formula.y.pos) +
  geom_smooth(method='lm', formula='y ~ x') +
  labs(title=paste0("Quality vs Disease Genes' FDR in ", dataid), 
       subtitle = paste("Correlation of -log10(FDR) vs |gene-quality correlation|:", round(cor.q.fdr.disgenes, 2)),
       x="gene-quality correlation", 
       y="-Log10 of False Discovery Rate (FDR)") +
  theme_minimal()  


# QUALITY VS DIFFERENTIALLY EXPRESSED GENES
# ______________________________________________________________________________
# cor vs dea genes
cor.deg.plot <- dcor %>%
  dplyr::filter(!is.na(padj)) %>% 
  group_by(cor_bin) %>%
  summarise(n=sum(diff_gene)) %>%
  ggplot(aes(x=cor_bin, y=n)) +
  geom_bar(stat="identity", colour="black", fill="white") +
  labs(title=paste0("Quality vs Differential Genes in ", dataid),
       subtitle=paste0("Differential genes: ", round(sum(dcor$diff_gene, na.rm=T), 2)), 
       x="gene-quality correlation", 
       y="Differential Genes") +
  geom_text(aes(label=n), position=position_dodge(width=0.9), vjust=-0.25) +
  theme_minimal()

cor.deg.perc.plot <- dcor %>% 
  dplyr::filter(!is.na(padj)) %>% 
  group_by(cor_bin) %>% 
  summarise(n=sum(diff_gene), total=n(), perc=n/total, label=paste(n, total, sep = "/")) %>% 
  ggplot(aes(x=cor_bin, y=perc)) +
  geom_bar(stat="identity", colour="black", fill="white") +
  labs(title=paste0("Quality vs Differential Genes (%) in ", dataid), 
       x="gene-quality correlation", 
       y="Percentage of Differential Genes") +
  geom_text(aes(label=label), position=position_dodge(width=0.9), vjust=-0.25) +
  theme_minimal()   

# cor.q.fdr.difgenes <- cor(abs((dcor %>% dplyr::filter(diff_gene==TRUE))$cor), -log10((dcor %>% dplyr::filter(diff_gene==TRUE))$padj), method=QI_METHOD, use="pairwise.complete.obs")
# missing_categories <- setdiff(tags, unique((dcor %>% dplyr::filter(diff_gene==TRUE))$cor_bin))
# missing_categ_rows <- tibble(cor_bin=missing_categories, padj=1, diff_gene=TRUE)
# dcor2 <- rbind(dcor %>% select(cor_bin, padj, diff_gene), missing_categ_rows)
# cor.deg.box.plot <- dcor2 %>% 
#   dplyr::filter(diff_gene==TRUE) %>%
#   # mutate(padj=ifelse(padj<pval_floor, pval_floor, padj)) %>% 
#   select(cor_bin, padj) %>% 
#   ggplot(aes(cor_bin, -log10(padj))) +
#   geom_boxplot() + 
#   labs(title=paste0("Quality vs Differential Genes' FDR in ", dataid), 
#        subtitle = paste("Correlation -log10(FDR) vs abs(gene-quality correlation):", round(cor.q.fdr.difgenes, 2)),
#        x="Correlation Coefficients", 
#        y="-Log10 of False Discovery Rate (FDR)") +
#   theme_minimal()  
# suppressWarnings(
#   plot.y.max <- max(ggplot_build(cor.deg.box.plot)$data[[1]]$ymax)*1.5
# )
# plot.y.min <- 0
# cor.deg.box.plot <- cor.deg.box.plot + coord_cartesian(ylim = c(plot.y.min,plot.y.max))

cor.deg.points.plot <- dcor %>% 
  dplyr::filter(diff_gene==TRUE) %>%
  # mutate(padj=ifelse(padj<pval_floor, pval_floor, padj)) %>%
  select(cor, padj) %>% 
  mutate(CorSign=ifelse(cor>=0,"POS", "NEG")) %>% 
  # ggplot(aes(abs(cor), -log10(padj))) +
  ggplot(aes(cor, -log10(padj), fill=CorSign)) +
  geom_point() + 
  stat_cor(label.x.npc = "left", label.y.npc = "top") +
  # stat_regline_equation(label.x = 0, label.y = max(-log10(dcor$padj)-2, na.rm=T)) +
  geom_smooth(method='lm', formula='y ~ x') +
  labs(title=paste0("Quality vs Differential Genes' FDR in ", dataid), 
       subtitle = paste("Correlation of -log10(FDR) vs |gene-quality correlation|:", round(cor.q.fdr.difgenes, 2)),
       x="gene-quality correlation", 
       y="-Log10 of False Discovery Rate (FDR)") +
  theme_minimal()  

# capture.output(suppressWarnings(suppressMessages(
  ggexport(plotlist = list(cor.genes.plot, cor.dis.plot,      cor.deg.plot, 
                           cor.fc.plot,    cor.dis.perc.plot, cor.deg.perc.plot, 
                           cor.points.plot,  cor.dis.points.plot,  cor.deg.points.plot), 
           # cor.pval.plot,  cor.dis.box.plot,  cor.deg.box.plot), 
           nrow = 3, ncol = 3, width=6200, height=4800, res=300, verbose=F,
           filename = file.path(out.dir.path, paste0(dataid, '.cor.png'))
           )
# )), file='/dev/null')


# QUALITY VS BATCH
# ______________________________________________________________________________

# cor.p.batch <- NA
pval.p.batch <- NA
n.batches <- 0

if(!is.na(sum(scores.table$batch)) & length(unique(scores.table$batch))>1){
  # cor.p.batch <- cor(scores.table$p, scores.table$batch)
  n.batches <- length(unique(scores.table$batch))
  pval.p.batch <- kruskal.test(p ~ batch, data = scores.table)$p.value
  cor.p.batch.plot <- scores.table %>% 
    mutate(batch=as.factor(batch)) %>% 
    ggplot(aes(x=p, y=batch)) + 
    geom_boxplot() +
    geom_jitter(aes(color=group, shape=batch), size=8, position=position_dodge(0.5)) + 
    labs(title=paste0("Quality vs Batches in ", dataid),
         subtitle = paste("Kruskal-Wallis rank sum test's p-value:", formatC(pval.p.batch, format = "e", digits = 2)),
         x="Low-quality probability"
    ) +
    coord_flip() +
    theme_minimal()
  
  batch.exp.samples.plot <- d0 %>% gather(key="id", value="value", -geneid) %>% 
    left_join(scores.table %>% select(id, p, group, batch, outlier, group_numeric), by="id") %>% 
#    filter(value>50) %>%
#    mutate(value=log2(value)) %>% 
    mutate(id=as.factor(id)) %>% 
    mutate(batch=as.factor(batch)) %>%
    ggplot(aes(x=value, y=reorder(id, value), color=batch)) + 
    geom_boxplot() +
    facet_grid(batch ~ ., scales="free_y", space="free_y") +
    labs(title=paste0("Genes Expression in Samples and Batches from ", dataid),
         x="Log2 Gene Expression", y="Samples"
    ) +
    theme_minimal()

  batch.exp.plot <- d0 %>% gather(key="id", value="value", -geneid) %>% 
    left_join(scores.table %>% select(id, p, group, batch, outlier, group_numeric), by="id") %>% 
#    filter(value>50) %>%
#    mutate(value=log2(value)) %>% 
    mutate(id=as.factor(id)) %>% 
    mutate(batch=as.factor(batch)) %>% 
    ggplot(aes(x=value, y=batch, color=batch)) + 
    geom_boxplot() +
    labs(title=paste0("Gene Expression in Batches from ", dataid),
         x="Log2 Gene Expression", y="Batch"
    ) +
    coord_flip() +
    theme_minimal()

  
  capture.output(suppressWarnings(suppressMessages(
    ggexport(cor.p.batch.plot, nrow = 1, ncol = 1, width=2400, height=2400, 
             res=300, verbose=F, 
             filename = file.path(out.dir.path, paste0(dataid, '.batch.png'))
    )
  )), file='/dev/null')
  
  capture.output(suppressWarnings(suppressMessages(
    ggexport(
      ggarrange(
        ggarrange(cor.p.batch.plot, batch.exp.plot, nrow = 2, labels = c("A", "B")),
        batch.exp.samples.plot,
        ncol = 2, 
        labels = c("", "C")                                        
      ) 
      , nrow = 1, ncol = 1, width=4800, height=4800, 
      res=300, verbose=F, 
      filename = file.path(out.dir.path, paste0(dataid, '.diag.batch.png'))
    )
  )), file='/dev/null')
  
#summary(d0 %>% gather(key="id", value, -geneid) %>% 
#          left_join(scores.table %>% select(id, p, group, batch, outlier, group_numeric), by="id"))
}

# TABLE OF STATISTICS
# ______________________________________________________________________________
n.samples <- nrow(scores.table)
if(ncol(d)>5){ dataset.stats.table <- data.frame(dataset=dataid,
                                 p_group_cor=dataset.p.group.cor,
                                 p_group_test=dataset.p.group.test, # !!!!!!!!!!!!!!!!!!!!!!!
                                 p_group_cor_batched=dataset.p.group.cor.batched,
                                 p_bases_cor=dataset.p.bases.cor, 
                                 p_bases_cor_ctrl=dataset.p.bases.ctrol.cor,
                                 p_bases_cor_treat=dataset.p.bases.treat.cor,
                                 p_bytes_cor=dataset.p.bytes.cor, 
                                 p_bytes_cor_ctrl=dataset.p.bytes.ctrol.cor,
                                 p_bytes_cor_treat=dataset.p.bytes.treat.cor,
                                 dunn_all=pca.valid0$dunn,
                                 dunn_noout=pca.valid1$dunn,
                                 dunn_noout_nocor=pca.valid2$dunn,
                                 dunn_noout_norand=pca.valid3$dunn,
                                 wbrat_all=pca.valid0$wb.ratio,
                                 wbrat_noout=pca.valid1$wb.ratio,
                                 wbrat_noout_nocor=pca.valid2$wb.ratio,
                                 wbrat_noout_norand=pca.valid3$wb.ratio,
                                 pgamma_all=pca.valid0$pearsongamma,
                                 pgamma_noout=pca.valid1$pearsongamma,
                                 pgamma_noout_nocor=pca.valid2$pearsongamma,
                                 pgamma_noout_norand=pca.valid3$pearsongamma,
                                 entropy_all=pca.valid0$entropy,
                                 entropy_noout=pca.valid1$entropy,
                                 entropy_noout_nocor=pca.valid2$entropy,
                                 entropy_noout_norand=pca.valid3$entropy,
                                 cor_q_fdr_allgenes=cor.q.fdr.allgenes,
                                 cor_q_fdr_pos_allgenes=cor.q.fdr.pos.allgenes,
                                 cor_q_fdr_neg_allgenes=cor.q.fdr.neg.allgenes,
                                 cor_q_fdr_disgenes=cor.q.fdr.disgenes,
                                 cor_q_fdr_pos_disgenes=cor.q.fdr.pos.disgenes,
                                 cor_q_fdr_neg_disgenes=cor.q.fdr.neg.disgenes,
                                 cor_q_fdr_difgenes=cor.q.fdr.difgenes,
                                 cor_q_fdr_pos_difgenes=cor.q.fdr.pos.difgenes,
                                 cor_q_fdr_neg_difgenes=cor.q.fdr.neg.difgenes,
                                 pval_p_batch=pval.p.batch,
                                 n_batches=n.batches,
                                 n_samples=n.samples,
                                 n_degs=n.degs,
                                 n_degs_fdr_only=n.degs.FdrOnly)
}else{dataset.stats.table <- data.frame(dataset=dataid,
                                    p_group_cor=dataset.p.group.cor, 
                                    p_group_test=dataset.p.group.test, # !!!!! => OK
                                    p_group_cor_batched=dataset.p.group.cor.batched,
                                    p_bases_cor=dataset.p.bases.cor, 
                                    p_bases_cor_ctrl=dataset.p.bases.ctrol.cor,
                                    p_bases_cor_treat=dataset.p.bases.treat.cor,
                                    p_bytes_cor=dataset.p.bytes.cor, 
                                    p_bytes_cor_ctrl=dataset.p.bytes.ctrol.cor,
                                    p_bytes_cor_treat=dataset.p.bytes.treat.cor,
                                    dunn_all=pca.valid0$dunn,
                                    dunn_noout=NA,
                                    dunn_noout_nocor=NA,
                                    dunn_noout_norand=NA,
                                    wbrat_all=pca.valid0$wb.ratio,
                                    wbrat_noout=NA,
                                    wbrat_noout_nocor=NA,
                                    wbrat_noout_norand=NA,
                                    pgamma_all=pca.valid0$pearsongamma,
                                    pgamma_noout=NA,
                                    pgamma_noout_nocor=NA,
                                    pgamma_noout_norand=NA,
                                    entropy_all=pca.valid0$entropy,
                                    entropy_noout=NA,
                                    entropy_noout_nocor=NA,
                                    entropy_noout_norand=NA,
                                    cor_q_fdr_allgenes=cor.q.fdr.allgenes,
                                    cor_q_fdr_pos_allgenes=cor.q.fdr.pos.allgenes,
                                    cor_q_fdr_neg_allgenes=cor.q.fdr.neg.allgenes,
                                    cor_q_fdr_disgenes=cor.q.fdr.disgenes,
                                    cor_q_fdr_pos_disgenes=cor.q.fdr.pos.disgenes,
                                    cor_q_fdr_neg_disgenes=cor.q.fdr.neg.disgenes,
                                    cor_q_fdr_difgenes=cor.q.fdr.difgenes,
                                    cor_q_fdr_pos_difgenes=cor.q.fdr.pos.difgenes,
                                    cor_q_fdr_neg_difgenes=cor.q.fdr.neg.difgenes,
                                    pval_p_batch=pval.p.batch,
                                    n_batches=n.batches,
                                    n_samples=n.samples,
                                    n_degs=n.degs,
                                    n_degs_fdr_only=n.degs.FdrOnly)
}

write_delim(dataset.stats.table, file = file.path(out.dir.path, paste0(dataid, ".cor.tsv")), delim="\t")


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
  # out.file.prefix="pos.kegg"
  # out.dir=out.dir.path
  # top_n_to_plot=40
  # table.cutoff=1
  # plot.cutoff=0.15
  # genesets = msigdb.kegg.path

  path0 <- file.path(out.dir, paste0(out.file.prefix, ".gsea.sets.tsv"))
  path1 <- file.path(out.dir, paste0(out.file.prefix, ".gsea.plot.png"))
  path2 <- file.path(out.dir, paste0(out.file.prefix, ".gsea.genes.annot.tsv"))
  
#  ranks <- deframe(gene_rank_table)
  ranks <- deframe(gene_rank_table[!duplicated(gene_rank_table$symbol),])
  pathways <- gmtPathways(genesets)
  # fgseaRes <- fgsea(pathways=pathways, stats=ranks, nperm=1000) # Examples of fgsea at https://stephenturner.github.io/deseq-to-fgsea/
  fgseaRes <- fgseaSimple(pathways=pathways, stats=ranks, nperm=1000, scoreType = "pos") # Examples of fgsea at https://stephenturner.github.io/deseq-to-fgsea/
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
    filter(padj<table.cutoff) %>% 
    mutate(significant=padj<plot.cutoff)

  if(dim(fgseaResTidy)[1]>0){
    
    write_tsv(fgseaResTidy, path0)
    
    gsea.plot <- ggplot(fgseaResTidy %>% top_n(top_n_to_plot, desc(abs(NES))), aes(reorder(pathway, NES), NES)) +
#      geom_col(aes(fill=significant)) +
      geom_col(aes(fill=-log10(padj))) +
      coord_flip() +
      labs(x="Top pathways", y="Normalized Enrichment Score",
           title=paste("Gene Set Enrichment Analysis")) + 
      theme_minimal(base_size = 11)
    ggsave(filename=path1, plot=gsea.plot, width=15.2, height=12.8)
    
    fgseaGenesAnnot <-  pathways %>% 
      enframe("pathway", "SYMBOL") %>% 
      unnest(cols = c(SYMBOL)) %>% 
      inner_join(gene_rank_table, by=c("SYMBOL"="symbol"))  %>% 
      filter(pathway %in% fgseaResTidy$pathway)
    write_tsv(fgseaGenesAnnot, path2)
    } else{
      print("quality_correlated_genes.r: no correlated genes passed the cutoff")
      file.create(path0)
      file.create(path1)
      file.create(path2)
    }
}

pos_genes <- dcor %>% 
  ungroup() %>% 
  filter(padj<0.05) %>% 
  select(symbol, log2FoldChange) %>% 
  filter(log2FoldChange>0) %>% 
  group_by(symbol) %>% 
  slice_max(order_by=log2FoldChange, n=1, with_ties = FALSE) %>% 
  ungroup() %>% 
  arrange(desc(log2FoldChange))

neg_genes <- dcor %>% 
  ungroup() %>% 
  filter(padj<0.05) %>% 
  select(symbol, log2FoldChange) %>% 
  filter(log2FoldChange<0) %>% 
  mutate(log2FoldChange=abs(log2FoldChange)) %>% 
  group_by(symbol) %>% 
  slice_max(order_by=log2FoldChange, n=1, with_ties = FALSE) %>% 
  ungroup() %>% 
  arrange(desc(log2FoldChange))
posneg_genes <- bind_rows(pos_genes, neg_genes) %>% arrange(desc(log2FoldChange))


# gsea.table.cutoff <- 0.5
if(nrow(pos_genes)>50){gsea_analysis(pos_genes, "posLFC.kegg", out.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.kegg.path)}
if(nrow(neg_genes)>50){gsea_analysis(neg_genes, "negLFC.kegg", out.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.kegg.path)}
if(nrow(pos_genes)>50){gsea_analysis(pos_genes, "posLFC.posi", out.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.posi.path)}
if(nrow(neg_genes)>50){gsea_analysis(neg_genes, "negLFC.posi", out.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.posi.path)}
if(nrow(pos_genes)>50){gsea_analysis(pos_genes, "posLFC.gs2d", out.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.gs2d.path)}
if(nrow(neg_genes)>50){gsea_analysis(neg_genes, "negLFC.gs2d", out.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.gs2d.path)}

if(nrow(posneg_genes)>50){gsea_analysis(posneg_genes, "posnegLFC.kegg", out.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.kegg.path)}
if(nrow(posneg_genes)>50){gsea_analysis(posneg_genes, "posnegLFC.posi", out.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.posi.path)}
if(nrow(posneg_genes)>50){gsea_analysis(posneg_genes, "posnegLFC.gs2d", out.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.gs2d.path)}



# ______________________________________________________________________________
# ORA FUNCTION
# ______________________________________________________________________________
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
  path1 <- file.path(out.dir, paste0(out.file.prefix, ".gsea.plot.png"))
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
  
  gsea.plot <- ggplot() +
    labs(x="Top pathways", y="Fold change",
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
      labs(x="Top pathways", y="Fold change",
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

# ______________________________________________________________________________
# ORA ANALYSIS
# ______________________________________________________________________________

ora.gene.universe <- unique(dcor[!is.na(dcor$baseMean),]$symbol)

pos_genes <- dcor %>% 
  filter(!is.na(baseMean)) %>% 
  ungroup() %>% 
  filter(padj<0.05) %>% 
  select(symbol, log2FoldChange) %>% 
  filter(log2FoldChange>0) %>% 
  mutate(log2FoldChange=abs(log2FoldChange)) %>% 
  group_by(symbol) %>% 
  slice_max(order_by=log2FoldChange, n=1, with_ties = FALSE) %>% 
  ungroup() %>% 
  arrange(desc(log2FoldChange))

neg_genes <- dcor %>% 
  filter(!is.na(baseMean)) %>% 
  ungroup() %>% 
  filter(padj<0.05) %>% 
  select(symbol, log2FoldChange) %>% 
  filter(log2FoldChange<0) %>% 
  mutate(log2FoldChange=abs(log2FoldChange)) %>% 
  group_by(symbol) %>% 
  slice_max(order_by=log2FoldChange, n=1, with_ties = FALSE) %>% 
  ungroup() %>% 
  arrange(desc(log2FoldChange))
posneg_genes <- bind_rows(pos_genes, neg_genes) %>% arrange(desc(log2FoldChange))


# fgsea.dir.path

top_pos_genes <- pos_genes$symbol[1:min(gsea.input_genes, nrow(pos_genes))]
top_neg_genes <- neg_genes$symbol[1:min(gsea.input_genes, nrow(neg_genes))]
gsea.pos.noBias.hall.plot <- ora_analysis("Hallmark pathways enrichment", "Positively Regulated Genes", top_pos_genes, ora.gene.universe, "fgsea.pos.hall", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.hallmark.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)
gsea.pos.noBias.posi.plot <- ora_analysis("Chromosomes enrichment", "Positively Regulated Genes", top_pos_genes, ora.gene.universe, "fgsea.pos.posi", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.posi.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff, create.default.files=F)
gsea.pos.noBias.cure.plot <- ora_analysis("Curated pathways enrichment", "Positively Regulated Genes", top_pos_genes, ora.gene.universe, "fgsea.pos.cure", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.curated.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff, create.default.files=T)
gsea.pos.noBias.path.plot <- ora_analysis("Canonical pathways enrichment", "Positively Regulated Genes", top_pos_genes, ora.gene.universe, "fgsea.pos.path", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.curated.cp.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff, create.default.files=F)
gsea.pos.noBias.regu.plot <- ora_analysis("Regulatory target genes enrichment", "Positively Regulated Genes", top_pos_genes, ora.gene.universe, "fgsea.pos.regu", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.regulatory.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)
gsea.pos.noBias.cell.plot <- ora_analysis("Cell type signature genes enrichment", "Positively Regulated Genes", top_pos_genes, ora.gene.universe, "fgsea.pos.cell", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.cells.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)
gsea.pos.noBias.gs2d.plot <- ora_analysis("Disease genes enrichment", "Positively Regulated Genes", top_pos_genes, ora.gene.universe, "fgsea.pos.gs2d", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.gs2d.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)

gsea.neg.noBias.hall.plot <- ora_analysis("Hallmark pathways enrichment", "Negatively Regulated Genes", top_neg_genes, ora.gene.universe, "fgsea.neg.hall", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.hallmark.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)
gsea.neg.noBias.posi.plot <- ora_analysis("Chromosomes enrichment", "Negatively Regulated Genes", top_neg_genes, ora.gene.universe, "fgsea.neg.posi", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.posi.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff, create.default.files=F)
gsea.neg.noBias.cure.plot <- ora_analysis("Curated pathways enrichment", "Negatively Regulated Genes", top_neg_genes, ora.gene.universe, "fgsea.neg.cure", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.curated.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff, create.default.files=T)
gsea.neg.noBias.path.plot <- ora_analysis("Canonical pathways enrichment", "Negatively Regulated Genes", top_neg_genes, ora.gene.universe, "fgsea.neg.path", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.curated.cp.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff, create.default.files=F)
gsea.neg.noBias.regu.plot <- ora_analysis("Regulatory target genes enrichment", "Negatively Regulated Genes", top_neg_genes, ora.gene.universe, "fgsea.neg.regu", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.regulatory.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)
gsea.neg.noBias.cell.plot <- ora_analysis("Cell type signature genes enrichment", "Negatively Regulated Genes", top_neg_genes, ora.gene.universe, "fgsea.neg.cell", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.cells.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)
gsea.neg.noBias.gs2d.plot <- ora_analysis("Disease genes enrichment", "Negatively Regulated Genes", top_neg_genes, ora.gene.universe, "fgsea.neg.gs2d", fgsea.dir.path, 40, gsea.table.cutoff, gsea.plot.cutoff, msigdb.gs2d.path, gsea.pathway.perc.cutoff, gsea.input.perc.cutoff)

