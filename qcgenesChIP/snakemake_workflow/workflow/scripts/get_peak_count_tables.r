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
                    "name",
                    "score",
                    "strand",
                    "thickstart",
                    "thickend",
                    "itemRGB",
                    "blockCount")
  # print(head(peakdf))

  peaks <- makeGRangesFromDataFrame(peakdf, keep.extra.columns=T)

  annot <- annotatePeak(peaks, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
  annot <- as_tibble(annot)

  return(annot)
})

# for (i in 1:length(bed_files)) {
#   peakdf <- read_delim(bed_files[[i]], col_names=F)
#   names(peakdf) <- c("chr",
#                     "start",
#                     "end",
#                     "name",
#                     "score",
#                     "strand",
#                     "thickstart",
#                     "thickend",
#                     "itemRGB",
#                     "blockCount")
#   peaks <- makeGRangesFromDataFrame(peakdf, keep.extra.columns=T)

# }

# Group peaks by gene and summarize the results
# to get the number of od peaks in genes over all datasets
# use number of datasets to normalize this number
combined_peaks <- do.call(rbind, annotated_peaks)
summarized_peaks <- combined_peaks %>% count(SYMBOL) %>% mutate(n_norm=n/length(bed_files))

# Save summarized peaks to a file
write.csv(summarized_peaks, output, row.names = FALSE)
