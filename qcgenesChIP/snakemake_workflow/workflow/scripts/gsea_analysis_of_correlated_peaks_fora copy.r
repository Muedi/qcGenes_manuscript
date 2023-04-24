suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(ChIPseeker, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(org.Hs.eg.db, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene, warn.conflicts = FALSE, quietly = TRUE, 
                         verbose = FALSE))
suppressMessages(library(fgsea, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(yaml, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

config.path <- file.path(".", "config", "config.yaml")
config.data <- yaml.load_file(config.path)
pathways.curated <- gmtPathways(config.data$MSIGDB_PATH) #  "/mnt/scratch1/projects/qcGenes.git/data/ref/msigdb/c2.all.v7.4.symbols.gmt")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
CORR_THRESHOLD <- 0.2
GSEA.INPUT.GENES = config.data$GSEA_INPUT_GENES
gsea.pathway.perc.cutoff <- config.data$GSEA_PATHWAY_REPRESENTATION_PERCENTAGE_CUTOFF
gsea.input.perc.cutoff <- config.data$GSEA_INPUTGENES_REPRESENTATION_PERCENTAGE_CUTOFF
gsea.table.cutoff <- config.data$GSEA_TABLE_PADJ_CUTOFF
gsea.plot.cutoff <- config.data$GSEA_PLOT_PADJ_CUTOFF2
top_n_to_plot<-40
argv <- commandArgs(trailingOnly=TRUE)

# GSEA Analysis Function ------------------------------------------------------------------------

correlated_pathways <- function(threshold, label){
  
  granges <- list()
  
  if(threshold > 0){
    for(i in 1:length(argv)){
      granges[[i]] <- read_csv(file = argv[i]) %>% 
        mutate(end = genomic_range + 999) %>% 
        dplyr::select(chromosome, start = genomic_range, end, cor, p, bin) %>% 
        filter(cor > threshold) %>% 
        makeGRangesFromDataFrame()
    }
  } else {
    for(i in 1:length(argv)){
      granges[[i]] <- read_csv(file = argv[i]) %>% 
        mutate(end = genomic_range + 999) %>% 
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
    summarize(dataset_corr_count = max(dataset_count)) %>%
    arrange(desc(dataset_corr_count))
  
  fgsea_stats <- genes_across_datasets[["dataset_corr_count"]]
  names(fgsea_stats) <- genes_across_datasets[["SYMBOL"]]
  
  fgseaRes <- as_tibble(fgsea::fgsea(pathways = pathways.curated, stats = fgsea_stats, minSize = 15, 
                           maxSize = 500, scoreType = "pos")) %>% 
    dplyr::select(-leadingEdge)
  print(fgseaRes)
  #ora instead of fgsea  
  universe <- genes_across_datasets %>% pull(SYMBOL)
  genes <- universe[1:GSEA.INPUT.GENES]
  
  pathways.universe <- unique(unlist(pathways.curated))
  pathways.universe <- pathways.universe[ pathways.universe %in% universe ]
  pathways_tibble <- enframe(pathways.curated) %>% 
    unnest_longer(value)
  pathways_tibble <- pathways_tibble[ pathways_tibble$value %in% universe, ]
  pathways_tibble <- pathways_tibble %>%
    group_by(name) %>%
    filter(n()>=5)
  pathways2 <- deframe(pathways_tibble %>% nest(data=c(value)))
  pathways2 <- lapply(pathways2, function(x) x$value )


  foraRes <- fora(pathways2, genes, universe)
  foraResTidy <- foraRes %>%
    as_tibble() %>%
    mutate(input_genes=length(genes)) %>%
    mutate(universe_genes=length(universe)) %>%
    filter(padj<gsea.table.cutoff) %>%
    filter((overlap/input_genes)>gsea.input.perc.cutoff) %>%
    filter((overlap/size)>gsea.pathway.perc.cutoff) %>%
    mutate(significant=padj<gsea.plot.cutoff) %>%
    mutate(fold_change= round((overlap/input_genes) / (size/universe_genes), 2)) %>%
    arrange(desc(fold_change))

  write_csv(foraResTidy, path = paste0("output/", label, "_correlated_pathways.csv"))
  write_csv(genes_across_datasets, path = paste0("output/", label, "_correlated_genes.csv"))
  print(foraRes)


  # ggplot(fgseaRes %>% filter(padj < 0.05)) + aes(reorder(pathway, NES), NES) + 
  #   geom_col(aes(fill=padj)) + 
  #   labs(x = "Pathway",  y = "Normalized Enrichment Score (NES)", 
  #       title = paste0("GSEA of ", label, " correlated Genes across Datasets")) + 
  #   coord_flip() + theme_minimal()
  # ggsave(paste0("output/", label, "_correlated_pathways.png"), width=15.2, height=12.8)

  gsea.plot <- foraResTidy %>% 
      filter(padj<gsea.plot.cutoff) %>%
      slice_max(fold_change, n=top_n_to_plot) %>%
      # ggplot(aes(reorder(pathway, fold_change), fold_change)) +
      mutate(pathway=str_to_lower(str_wrap(str_replace_all(str_trunc(pathway, 60), "_", " "), 60))) %>% 
      ggplot(aes(reorder(pathway, fold_change), fold_change)) +
      geom_col(aes(fill=-log10(padj))) +
      coord_flip() +
      labs(x="Top pathways", y="Fold change",
           title=paste("Curated pathways enrichment"),
           subtitle=paste0(label, " correlated genes")
           ) +
      # theme_minimal(base_size = 10)
    theme_minimal()
    ggsave(filename=paste0("output/", label, "_correlated_pathways.png"), plot=gsea.plot, width=8, height=8)

}

# Positively Correlated -----------------------------------------------------------------------

correlated_pathways(+CORR_THRESHOLD, "positively")

# Negatively Correlated -----------------------------------------------------------------------

correlated_pathways(-CORR_THRESHOLD, "negatively")
