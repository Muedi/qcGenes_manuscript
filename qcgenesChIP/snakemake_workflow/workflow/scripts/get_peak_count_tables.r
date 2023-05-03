suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(ChIPseeker, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(org.Hs.eg.db, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene, warn.conflicts = FALSE, quietly = TRUE, 
                         verbose = FALSE))
suppressMessages(library(fgsea, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

# # for testing
# setwd(file.path(getwd(), "qcgenesChIP/snakemake_workflow"))
# folder <- "output/macs2/GSE107734/"
# bed_files <- grep("\\peaks.bed$",list.files(folder), value = T)
# bed_files <- file.path(folder, bed_files)
# output <- c("output/counts/GSE107734/peak_counts.csv", 
#             "output/counts/GSE107734/peak_counts_normalized.csv")


pathways.hallmark <- gmtPathways("/mnt/scratch1/projects/qcGenes.git/data/ref/msigdb/c2.all.v7.4.symbols.gmt")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
CORR_THRESHOLD = 0.3

gene_lengths <- data.frame(gene_id=genes(txdb)$gene_id, width=width(genes(txdb)))
symbols <- mapIds(org.Hs.eg.db,
                  keys=gene_lengths$gene_id,
                  column="SYMBOL",
                  keytype="ENTREZID")
gene_lengths$SYMBOL <- symbols 



argv <- commandArgs(trailingOnly=TRUE)
# Filter the arguments to keep only bed files
bed_files <- grep("\\.bed$", argv, value = TRUE)
scores  <- grep("\\.txt$", argv, value = TRUE)
output  <- grep("\\.csv$", argv, value = TRUE)
print(output)

# Annotate each bed file and combine them into a single dataframe
annotated_peaks <- lapply(bed_files, function(bed_file) {
  peakdf <- read_delim(bed_file, col_names=F)
  
  names(peakdf) <- c("chr",
                    "start",
                    "end",
                    "length",
                    "abs_summit",
                    "pileup",
                    "-log10(pvalue)",
                    "fold_enrichment",
                    "-log10(qvalue)",
                    "name")


  # print(head(peakdf))
  peaks <- makeGRangesFromDataFrame(peakdf, keep.extra.columns=T)

  annot <- annotatePeak(peaks, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
  annot <- as_tibble(annot)
  return(annot)
})
# counted peaks
names <- str_extract_all(bed_files, "SRR\\d+")
names(annotated_peaks) <- names
# loop over list 
annot2 <- list()
for (i in 1:length(names)) {
   annot2[[names[[i]]]] <- annotated_peaks[[names[[i]]]] %>% count(SYMBOL)
}
summarized_peaks <- annot2 %>% purrr::reduce(full_join, by = "SYMBOL") %>% set_names("SYMBOL", names(annot2))
print(paste0("Output1: ", output[[1]]))
write.csv(summarized_peaks, output[[1]], row.names = FALSE)
print("written.")


summarized_peaks_NORM <- summarized_peaks %>% 
  left_join(gene_lengths %>% dplyr::select(SYMBOL, width), by="SYMBOL") %>%
  mutate_if(is.numeric, funs( . / width)) %>% dplyr::select(-width)
print(paste0("Output1: ", output[[2]]))
write.csv(summarized_peaks_NORM, output[[2]], row.names = FALSE)
print("written.")

# TODO: 
# compute correlations!
# add reads in peaks (only need to filter it out of the annot df, no hasssle! :)
# the current counting is faulty, so I'll fall back on counted peaks atm
  # names(peakdf) <- c("chr",
  #                   "start",
  #                   "end",
  #                   "name",
  #                   "score",
  #                   "strand",
  #                   "thickstart",
  #                   "thickend",
  #                   "itemRGB",
  #                   "blockCount",
  #                   "blockSizes",
  #                   "blockStarts")