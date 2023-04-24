suppressMessages(library(ChIPseeker, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

args <- commandArgs(trailingOnly = TRUE)
bed_file = args[1]
assembly = args[2]
out_file = args[3]

txdb = NULL
if (assembly == 'hg38') {
  suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
	txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
}
if (assembly == 'mm10') {
  suppressMessages(library(TxDb.Mmusculus.UCSC.mm10.knownGene, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
	txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
}
if (assembly == 'dm6') {
  suppressMessages(library(TxDb.Dmelanogaster.UCSC.dm6.ensGene, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
	txdb = TxDb.Dmelanogaster.UCSC.dm6.ensGene
}
if (assembly == 'ce10') {
  suppressMessages(library(BSgenome.Celegans.UCSC.ce10, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
	txdb = BSgenome.Celegans.UCSC.ce10
}

readsAnno = annotatePeak(bed_file, tssRegion=c(-1000, 150), TxDb=txdb, verbose=FALSE)
percentages = data.frame(show(readsAnno))

write.table(percentages,file=out_file, sep="\t")


