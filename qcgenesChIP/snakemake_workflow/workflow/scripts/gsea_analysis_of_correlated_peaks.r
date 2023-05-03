suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(ChIPseeker, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(org.Hs.eg.db, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene, warn.conflicts = FALSE, quietly = TRUE, 
                         verbose = FALSE))
suppressMessages(library(fgsea, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

pathways.hallmark <- gmtPathways("/mnt/scratch1/projects/qcGenes.git/data/ref/msigdb/c2.all.v7.4.symbols.gmt")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
CORR_THRESHOLD = 0.3

argv <- commandArgs(trailingOnly=TRUE)

# GSEA Analysis Function ------------------------------------------------------------------------

correlated_pathways <- function(threshold, label){
  
  granges <- list()
  
  if(threshold > 0){
    for(i in 1:length(argv)){
      granges[[i]] <- read_csv(file = argv[i]) %>% 
        mutate(end = genomic_range + 499) %>% 
        dplyr::select(chromosome, start = genomic_range, end, cor, p, bin) %>% 
        filter(cor > threshold) %>% 
        makeGRangesFromDataFrame()
    }
  } else {
    for(i in 1:length(argv)){
      granges[[i]] <- read_csv(file = argv[i]) %>% 
        mutate(end = genomic_range + 499) %>% 
        dplyr::select(chromosome, start = genomic_range, end, cor, p, bin) %>% 
        filter(cor < threshold) %>% 
        makeGRangesFromDataFrame()
    }
  }

  full_range <- sort.GenomicRanges(unique(do.call(c, granges)))
  
  overlaps <- vector(mode="integer", length=length(full_range))
  
  for(i in 1:length(granges)){
    overlaps <- overlaps + countOverlaps(full_range, granges[[i]])
  }
  
  genes_across_datasets <- 
    as_tibble(annotatePeak(full_range, tssRegion=c(-6000, 6000), 
                           TxDb=txdb, annoDb = "org.Hs.eg.db")) %>% 
    mutate(dataset_count = overlaps) %>% 
    group_by(SYMBOL) %>% 
    summarize(dataset_corr_count = max(dataset_count))
  
  fgsea_stats <- genes_across_datasets[["dataset_corr_count"]]
  names(fgsea_stats) <- genes_across_datasets[["SYMBOL"]]
  print("fgsea_stats: ")
  print(fgsea_stats)
  fgseaRes <- as_tibble(fgsea::fgsea(pathways = pathways.hallmark, stats = fgsea_stats, minSize = 15, 
                           maxSize = 500, scoreType = "pos")) %>% 
    dplyr::select(-leadingEdge)
  
  write_csv(fgseaRes, path = paste0("output/", label, "_correlated_pathways.csv"))
  write_csv(genes_across_datasets, path = paste0("output/", label, "_correlated_genes.csv"))
  
  ggplot(fgseaRes %>% filter(padj < 0.05) %>% arrange(padj) %>% slice_head(n=25)) + aes(reorder(pathway, NES), NES) + 
    geom_col(aes(fill=padj)) + 
    labs(x = "Pathway",  y = "Normalized Enrichment Score (NES)", 
        title = paste0("GSEA of ", label, " correlated Genes across Datasets")) + 
    coord_flip() + theme_minimal()
  ggsave(paste0("output/", label, "_correlated_pathways.png"), width=15.2, height=12.8)

}

# Positively Correlated -----------------------------------------------------------------------

correlated_pathways(+CORR_THRESHOLD, "positively")

# Negatively Correlated -----------------------------------------------------------------------

correlated_pathways(-CORR_THRESHOLD, "negatively")
