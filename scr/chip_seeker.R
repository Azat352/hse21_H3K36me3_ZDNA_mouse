#!/usr/bin/R
source('lib.R')

###
.libPaths(c("~/libr", .libPaths()))
#if (!requireNamespace("BiocManager", quietly = TRUE))
#	install.packages("BiocManager")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
BiocManager::install("ChIPseeker")

library(ChIPseeker)
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)

###

#NAME <- 'H3K4me3_A549.intersect_with_DeepZ'
#NAME <- 'mouseZ-DNA2'
NAME <- 'H3K36me3_MEL.ENCFF390YLK.filtered'
# NAME <- 'H3K36me3_MEL.ENCFF503JQE.filtered'
BED_FN <- paste0(DATA_DIR, NAME, '.bed')

###

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

peakAnno <- annotatePeak(BED_FN, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

pdf(paste0(OUT_DIR, 'chip_seeker.', NAME, '.plotAnnoPie.pdf'))
png(paste0(OUT_DIR, 'chip_seeker.', NAME, '.plotAnnoPie.png'))
plotAnnoPie(peakAnno)
dev.off()

peak <- readPeakFile(BED_FN)
pdf(paste0(OUT_DIR, 'chip_seeker.', NAME, '.covplot.pdf'))
covplot(peak, weightCol="V5")
dev.off()

