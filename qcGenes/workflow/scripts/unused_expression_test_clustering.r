suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE) )
suppressMessages(library(yaml, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE) )
suppressMessages(library(DESeq2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE) )
suppressMessages(library(biomaRt, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE) )
suppressMessages(library(tximport, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE) )
suppressMessages(library(GenomicFeatures, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(writexl, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(fgsea, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(data.table, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(org.Hs.eg.db, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))


# BiocManager::install(c("DESeq2",
# "biomaRt",
# "tximport",
# "GenomicFeatures",
# "org.Hs.eg.db"
# ))
# install.packages("writexl")

# ARGUMENTS
# ------------------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
if (length(args)<2 | length(args)>2) {
  stop("Please provide GEO_Series ID and subdir as command line parameter, and optionally a subdirectory name", call.=FALSE)
}
dataid <- "GSE133039sub_diff_tu_low_0"
subdir <- "main"

# dataid = "GSE113255"
# subdir = "main"

out.path <- file.path("review", "expr_clustering_stuff", dataid)
quant.path <- file.path("data", "output", dataid, "salmon")
cor.genes.path <- file.path("output", "main", "qualityVsExp", dataid, paste0(dataid,".cor_pos_genes.tsv"))
dir.create(out.path, recursive=T)


### global variables from config
# ------------------------------------------------------------------------------
# ENSEMBL <- config$ENSEMBL
config <- yaml.load_file('config/config.yaml')
TRANSCRIPTOME <- config$TRANS_PATH
GENOME <- config$GENOM_PATH
ANNOTATION <- config$ANNOT_PATH
BIOMART_ENS_IDS <- config$ENSID_PATH
ensembl_dataset <- "hsapiens_gene_ensembl"
fdr.thr <- 0.05  # threshold of FDR/adjusted P-value for significantlly differentially expressed genes
# num.control <- length(controls)  # number of comparisons that the user wants to do
# num.treat <- length(treats)  # should equals to num.control
msigdb.kegg.path   <- config$MSIGDB_CURATED_PATH
msigdb.posi.path   <- config$MSIGDB_POSITIONAL_PATH
msigdb.gs2d.path   <- config$GS2D_MSIGDB_PATH
gsea.table.cutoff   <- config$GSEA_TABLE_PADJ_CUTOFF
gsea.plot.cutoff   <- config$GSEA_PLOT_PADJ_CUTOFF2

# metadata 
# ------------------------------------------------------------------------------
# read mega file
sra.path    <- file.path(".", "config", "metadata", "Mega_SraRunTable.csv")
datasets.path <- file.path(".", "config", "metadata", "Datasets.csv")
meta.data <- read_csv(sra.path, show_col_types=F)  %>% filter(GEO_Series == dataid)
meta.data <- meta.data %>% filter(Selected == 1) %>% dplyr::select(Run, group, Age, gender)

dataset.meta <- read_csv(datasets.path, show_col_types=F) %>% filter(GEO_Series == dataid)
control <- dataset.meta$Control
treat <- dataset.meta$Treat
# meta.data$group <- as.factor(meta.data$group)
# meta.data$batch <- as.factor(meta.data$batch)
# meta.data$gender <- as.factor(meta.data$gender)
# meta.data$Age <- as.factor(meta.data$Age)

# functions
# ------------------------------------------------------------------------------
# function from rasflow, works fine :)
convert.id2symbol <- function(gene.id) {
  gene.symbol <- gene.id  # initialize the gene symbol with the gene id
  
  # it may happen that no symbol can be found for any id. In that case, "queryMany" will throw an error
  # so "try" is used here to take care of that error
  try({
    gene.symbol.all <- queryMany(gene.id, scopes = 'ensembl.gene', fields = 'symbol')
    
    h <- hash()
    for (i in 1:nrow(gene.symbol.all)) {
      query <- gene.symbol.all$query[i]
      symbol <- gene.symbol.all$symbol[i]
      if (has.key(query, h)) {  # if there's duplicate for the same query
        h[[query]] <- paste(hash::values(h, keys = query), symbol, sep = ', ')
      } else {
        if (is.na(symbol)) {  # if there's no hit for the query, keep the original id
          h[[query]] <- query
        } else {
          h[[query]] <- symbol
        }
      }
    }
    
    for (i in c(1:length(gene.symbol))) {
      gene.symbol[i] <- h[[gene.id[i]]]
    }
  })
  
  return(gene.symbol)
}

# builöd noVersion file to go for gene level analysis (like in rasflow)
# remove the version in the transcript ID
remove_version <- function(files) {  # input files (file names with directory) are output from Salmon
  for (i in c(1:length(files))) {
      quant.file <- files[i]
      quant.table <- read.table(quant.file, header = TRUE, stringsAsFactors = FALSE)
      trans.id.version <- quant.table$Name
      trans.id <- rep('ID', length(trans.id.version))
      for (j in c(1:length(trans.id.version))) {
      trans.id[j] <- strsplit(trans.id.version[j], ".", fixed = TRUE)[[1]][1]
      }
      quant.table$Name <- trans.id
      write.table(quant.table, file.path(quant.path, samples[i], "quant_noVersion.sf"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  return(trans.id)
}


# # load salmons tpms and merge files.
# # ------------------------------------------------------------------------------
# names.vector <- unique(meta.data$Run)
# tpm.file.names.vector <- paste0(names.vector, "_tpm.tsv")
# tpm.file.paths.vector <- file.path(quant.path, tpm.file.names.vector)
# table.file.name <- "all.samples.tpm.tsv"
# table.file.name.xls <- "all.samples.tpm.xlsx"

# # loop over all and merge
# data.tpm <- read.delim(tpm.file.paths.vector[1], header=F) %>% as_tibble()
# colnames(data.tpm) <- c("trans_id_ver", names.vector[1])

# for(i in 2:length(names.vector)){
#   f <- tpm.file.paths.vector[i]
#   temp_table <- read.delim(f, header=F) %>% as_tibble()
#   colnames(temp_table) <- c("trans_id_ver", names.vector[i])
#   data.tpm <- data.tpm %>% left_join(temp_table, by = "trans_id_ver")
# }
# data.tpm <- data.tpm %>% 
#   mutate(trans_id_ver=as.character(trans_id_ver)) %>% 
# #  mutate(total = select(., sample1:sample18) %>% rowSums) %>% 
#   mutate(total = dplyr::select(., -trans_id_ver) %>% rowSums) %>% 
#   filter(total>0) %>% 
#   dplyr::select(-total)

# biomart <- as_tibble(read.delim(ANNOTATION, header=F,  comment.char = "#", stringsAsFactors=F)) %>% 
#   filter(V3=="transcript") %>% 
#   mutate(gene_id=str_extract(V9, regex("(?<=gene_id )(.+?)(?=;)", dotall = TRUE))) %>% 
#   mutate(gene_ver=str_extract(V9, regex("(?<=gene_version )(.+?)(?=;)", dotall = TRUE))) %>% 
#   mutate(gene_id_ver=paste0(gene_id, ".", gene_ver)) %>% 
#   mutate(gene_name=str_extract(V9, regex("(?<=gene_name )(.+?)(?=;)", dotall = TRUE))) %>% 
#   mutate(trans_id=str_extract(V9, regex("(?<=transcript_id )(.+?)(?=;)", dotall = TRUE))) %>% 
#   mutate(trans_ver=str_extract(V9, regex("(?<=transcript_version )(.+?)(?=;)", dotall = TRUE))) %>% 
#   mutate(trans_id_ver=paste0(trans_id, ".", trans_ver)) %>% 
#   mutate(trans_name=str_extract(V9, regex("(?<=transcript_name )(.+?)(?=;)", dotall = TRUE))) %>% 
#   dplyr::select(trans_id_ver, trans_name, gene_id_ver, gene_name)

# data.tpm %>% 
#   left_join(biomart, by="trans_id_ver") %>% 
#   dplyr::select(trans_id_ver, trans_name, gene_id_ver, gene_name, everything()) %>% 
#   write_tsv(file.path(out.path, table.file.name))

# data.tpm %>% 
#   left_join(biomart, by="trans_id_ver") %>% 
#   dplyr::select(trans_id_ver, trans_name, gene_id_ver, gene_name, everything()) %>% 
#   writexl::write_xlsx(file.path(out.path, table.file.name.xls))

# Deseq2, load data from quants  perform differental analysis
# ------------------------------------------------------------------------------
# build tx2gene from gtf:
annot.file  <- config$ANNOT_PATH
txdb <- makeTxDbFromGFF(file=annot.file)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- biomaRt::select(txdb, k, "GENEID", "TXNAME")

### load salomns quants as deseq2 object.
### the original quant files from Salmon
samples.all <- meta.data$Run
group.all <- meta.data$group

# get complete counts
samples <- factor(samples.all) 
group <- factor(group.all) 
group <- relevel(group, ref = control)

### import quantification as txi
# because ensemble is used 
files <- file.path(quant.path, samples, "quant.sf")
names(files) <- samples

trans.id <- remove_version(files)
files.noVersion <- file.path(quant.path, samples, "quant_noVersion.sf")
names(files.noVersion) <- samples

# load data from salom, gene level
txi_all <- tximport(files.noVersion, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "no")
colData = data.frame(samples, group)
design <- model.matrix(~group)
dds <- DESeqDataSetFromTximport(txi_all, colData = colData, design = design)
dds <- estimateSizeFactors(dds)
# gene abundance
counts_deseq2 <- counts(dds, normalized=FALSE)
counts_deseq2 <- as_tibble(counts_deseq2, rownames="ensembl_id")
write.table(counts_deseq2, file.path(out.path, "all_samples_gene_abundance.tsv"), sep = "\t", row.names=F)

# gene norm 
normalized_counts_deseq2 <- counts(dds, normalized=TRUE)
normalized_counts_deseq2 <- as_tibble(normalized_counts_deseq2, rownames="ensembl_id")
write.table(normalized_counts_deseq2, file.path(out.path, "all_samples_gene_norm.tsv"), sep = "\t", row.names=F)
# gene norm + rlog
normalized_counts_deseq2 <- counts(dds, normalized=TRUE)
normalized_counts_deseq2 <- rlog(round(normalized_counts_deseq2), blind=T)
normalized_counts_deseq2 <- as_tibble(normalized_counts_deseq2, rownames="ensembl_id")
### TODO: check if it actually makes sense to compute rlog from normalized counts!


# Deseq2,  perform differental analysis
# ------------------------------------------------------------------------------
## perform DEA
txi_all <- tximport(files.noVersion, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "no")
dds <- DESeqDataSetFromTximport(txi_all, colData = colData, design = design)
# Filtering
filter.need = T
if (filter.need) {
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
}
dds <- estimateSizeFactors(dds, type="poscounts")
dds <- DESeq(dds)

## export the results
res.dea <- results(dds)
res.dea <- res.dea[complete.cases(res.dea), ]  # remove any rows with NA

dea <- as.data.frame(res.dea)
dea <- dea[order(dea$padj, -abs(dea$log2FoldChange), decreasing = FALSE), ]  # sort the table: ascending of FDR then descending of absolute valued of logFC

deg <- dea[dea$padj < 0.05, ]
if (nrow(deg) > 1) {
  deg <- deg[order(deg$padj, -abs(deg$log2FoldChange), decreasing = FALSE), ]  # sort the table: ascending of FDR then descending of absolute valued of logFC
}
# save the DEA result and DEGs to files
write.table(dea, paste(out.path, '/dea_', control, '_', treat, '.tsv', sep = ''), row.names = T, col.names = NA, quote = FALSE, sep = '\t')
write.table(deg, paste(out.path, '/deg_', control, '_', treat, '.tsv', sep = ''), row.names = T, col.names = NA, quote = FALSE, sep = '\t')

cor.genes <- read_delim(cor.genes.path, delim = "\t") %>% dplyr::select(geneid, SYMBOL, cor)

length(intersect(cor.genes$geneid, rownames(deg)))

# hierarchical clustering for degs, check if correlating genes fall into same clusters.

# cluster
counts <- counts(dds)
counts <- counts[rownames(deg),]

hc_deg <- counts %>% 
    scale() %>%
  dist(method = "euclidian" ) %>% 
  hclust(method = "ward.D") 

pdf(file.path(out.path, "tree.pdf"), width=7, height=7)
plot(hc_deg, labels=F)
dev.off()

clusters_hc <-  cutree(hc_deg, k=4)

cor.genes <- cor.genes %>% mutate(clust = enframe(clusters_hc[cor.genes$geneid], name = "ID", value = "NewValue")$NewValue)

write_csv(cor.genes, file.path(out.path, "cor.genes.clust.csv"))
# plot clusters
# inp <- results_paired_th17_comp[results_paired_th17_comp$clust == 2,]
inp <- results_paired_th17_comp

# shape
clust_list <- ifelse(
  inp$clust == 1, 1, ifelse(
    inp$clust == 2, 2, ifelse(
      inp$clust == 3, 3, 4
      )
    )
  )

names(clust_list)[clust_list == 1] <- '1'
names(clust_list)[clust_list == 2] <- '2'
names(clust_list)[clust_list == 3] <- '3'
# names(clust_list)[clust_list == 4] <- '4'
# color
keyvals <- ifelse(
  inp$logFC > logFC_iact & inp$adj.P.Val < p_val_thresh , '#9C0E0F', # th17
  'darkgrey')
keyvals[is.na(keyvals)] <- 'darkgrey'
names(keyvals)[keyvals == '#9C0E0F'] <- 'Th17_interactor'
names(keyvals)[keyvals == 'darkgrey'] <- 'ns'



EnhancedVolcano(inp,
                lab = rownames(inp),
                x = "logFC", 
                y = "adj.P.Val",
                selectLab = rownames(inp)[which(names(keyvals) %in% c('Th17_interactor'))],
                colCustom = keyvals,
                # outcomment line below for complete volcano
                xlim = c(0, max(inp[["logFC"]], na.rm = TRUE) +1.5),
                pCutoff = p_val_thresh,
                FCcutoff = logFC_iact,
                pointSize = 1.5,
                labSize = 4.0,
                shapeCustom = clust_list,
                legendLabels = c("NS", expression(Log[2] ~ FC), "p-value", "Both"),
                title = "Th17 limma results",
                subtitle = "imputed data (min + knn)") 
ggsave(file = file.path(outdir,'/paired_volcanoplot_th17_clusters.pdf'), width=2000, height=2000, units = "px", dpi=300)
