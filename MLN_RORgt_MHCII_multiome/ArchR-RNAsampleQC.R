setwd("~/0-workspace/CCR7_DC/MLN_RORgt_MHCII_multiome/")

suppressPackageStartupMessages(library(ArchR))
library(GenomicRanges)
library(Biostrings)
library(Matrix)
suppressPackageStartupMessages(library(dplyr))
library(parallel)
library(BSgenome.Mmusculus.UCSC.mm10)
addArchRThreads(threads = 16)
addArchRGenome("mm10")
set.seed(1)

setwd("~/0-workspace/CCR7_DC/MLN_RORgt_MHCII_multiome/ArchR/")

archr.obj <- loadArchRProject("ArchROutput/")

# Restrict to 10,468 cells including DC and TC that pass RNA QC ####
rna.md <- read.csv("../Seurat/results/meta-data.csv", row.names = 1)
rownames(rna.md) <- paste0(archr.obj$Sample[[1]], "#", rownames(rna.md))
idx <- rownames(rna.md)

pref.archr <- "ArchROutput-RNAsampleQC/"
archr.obj <- subsetArchRProject(
  ArchRProj = archr.obj, cells = idx, outputDirectory = pref.archr, force = T
)
# numberOfCells(1): 10468
# medianTSS(1): 21.6165
# medianFrags(1): 13811.5

archr.obj <- addIterativeLSI(archr.obj, iterations = 5, varFeatures = 50000, scaleDims = T, force = T) %>%
  addUMAP(name = "UMAP", nNeighbors = 30, force = T)

### Save res ####
write.csv(getEmbedding(archr.obj), file = paste0(pref.archr, "Embeddings/UMAP.csv"), quote = F)
write.csv(getReducedDims(archr.obj), file = paste0(pref.archr, "IterativeLSI/LSI.csv"), quote = F)
write.csv(archr.obj@cellColData, file = paste0(pref.archr, "cellColData.csv"), quote = F)

archr.obj <- saveArchRProject(archr.obj)
