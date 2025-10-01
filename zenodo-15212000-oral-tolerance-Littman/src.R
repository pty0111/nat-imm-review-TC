# Littman
suppressPackageStartupMessages({
  library(curl)
  library(Seurat)
  library(ArchR)
  library(Rmagic)
  library(anndata)
  library(Matrix)
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(rasterpdf)
  library(viridis)
  library(future)
})
reticulate::py_require("magic-impute")
reticulate::py_require("anndata")

set.seed(1)
plan("multicore", workers = 16)
options(future.globals.maxSize = Inf)

source("utils.R")

pal <- list(
  Clusters = c(
    "#00bf00", "#489de8", "#d40663", "#f8c72f", "#077315",
    "#785cd4", "#e67109", "#0eefff", "#f081e6", "#260691",
    "#49709c", "#9e7d3f", "#bd537a", "#4e225c", "#f202ed",
    "#fec55f", "#062e0b", "#9c6fa8", "#078d94", "#5c1a1a",
    "#827c68", "#aebeff", "#9c2903", "#ffc5af", "#4f5715",
    "#0249f0", "#f43525", "#0077ff", "#7f227e", "#dfddff",
    "#7e85d7", "#fff64f", "#5fed0e", "#543018", "#f31220"
  ),
  Clusters_long = c(
    "#000000", "#FAD09F", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
    "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
    "#5A0007", "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
    "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
    "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
    "#372101", "#FFFF00", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
    "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
    "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
    "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
    "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
    "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
    "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
    "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")
)
pal$RabiClusters <- pal$Clusters_long[1:length(unique(sro$RabiClusters))]
names(pal$RabiClusters) <- sort(unique(sro$RabiClusters))
pal$RabiClusters['na'] <- "lightgray"

dir.create("Seurat")
pref.p.sro <- "plots/Seurat/paper/"
dir.create(pref.p.sro, recursive = T)

sro3 <- readRDS("data/mouse/OUTPUT/annotatedPlus7kb.rds") # main data fig. 2a
write.csv(sro3@meta.data, "Seurat/paper/meta-data.csv", quote = F)

# 1. Re-process DeltaPlus7kb data ####
pref.sro <- 'Seurat/merged/'; pref.p.sro <- 'plots/Seurat/merged/'
dir.create(pref.sro); dir.create(pref.p.sro)

data.folder <- "data/2024-04-26/cellranger/"
samples <- c(ctrl = "count-7KB-CTRL/", mut = "count-7KB-MUT")

sro.list <- lapply(1:length(samples), function(i){
  counts <- Read10X_h5(paste0(data.folder, samples[i], "/outs/filtered_feature_bc_matrix.h5"), 
                       use.names = TRUE, unique.features = TRUE)
  sro <- CreateSeuratObject(
    project = names(samples)[i],
    counts = counts
  )
  new.cell.names <- paste0(rownames(sro@meta.data), "_", i)
  sro <- RenameCells(sro, new.names = new.cell.names)
  return(list(sro = sro))
})
sro <- merge(sro.list[[1]]$sro, y = sapply(sro.list,"[[",1)[2:length(samples)], project = "Prdm16")
sro$sample <- sro$orig.ident
sro <- PercentageFeatureSet(sro, pattern = "^mt-", col.name = "MtFrac_RNA")

thr <- data.frame(
  nf.min = 500, nf.max = 10000, mp.max = 5)
cell.discard <- 
  sro$nFeature_RNA < thr$nf.min | sro$nFeature_RNA > thr$nf.max |
  sro$MtFrac_RNA > thr$mp.max
table(cell.discard)
# cell.discard
# FALSE  TRUE 
# 23091   702 

Idents(sro) <- sro$sample
dir.create("plots/Seurat/merged/QC", recursive = T)
pdf("plots/Seurat/merged/QC/violin-sample-RNA-QC.pdf", width = 8, height = 6)
plot.all.QC(sro, ident = "sample", thr = thr)
dev.off()

# sro <- sro[, !cell.discard] # not run here
write.csv(sro@meta.data, file = "Seurat/initial-meta-data.csv")
saveRDS(sro, file = "Seurat/initial-SRO.rds")

# 2. RNA analysis ####
pref.sro <- 'Seurat/merged/'; pref.p.sro <- 'plots/Seurat/merged/'
sro <- readRDS(paste0('Seurat/initial-SRO.rds'))
sro <- sro[, !cell.discard]

s.genes <- rownames(sro)[toupper(rownames(sro)) %in% toupper(cc.genes.updated.2019$s.genes)]
g2m.genes <- rownames(sro)[toupper(rownames(sro)) %in% toupper(cc.genes.updated.2019$g2m.genes)]
sro <- CellCycleScoring(sro, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
sro$Phase <- gsub("G2M", "G2/M", sro$Phase)

res.list <- seq(0.6, 2, 0.2)
cell.count <- rowSums(sro@assays$RNA@counts > 0)
sro <- NormalizeData(sro[cell.count > 1, ]) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 5000)
sro <- ScaleData(sro, features = rownames(sro)) %>%
  RunPCA(features = rownames(sro), npcs = 50) %>%
  FindNeighbors(dims = 1:30, k.param = 30) %>%
  FindClusters(resolution = res.list) %>%
  RunUMAP(dims = 1:30, n.neighbors = 30, metric = "cosine", min.dist = 0.4, spread = 1)

paper.md <- read.csv("Seurat/paper/meta-data.csv", row.names = 1)
sro$RabiClusters <- paper.md[rownames(sro@meta.data), ]$RabiClusters
sro$RabiClusters[is.na(sro$RabiClusters)] <- "na"

## save res ####
write.csv(sro@reductions$pca@cell.embeddings, file = paste0(pref.sro, "PCA.csv"), quote = F)
write.csv(sro@reductions$umap@cell.embeddings, file = paste0(pref.sro, "UMAP.csv"), quote = F)
write.csv(sro@meta.data, file = paste0(pref.sro, "meta-data.csv"), quote = F)
saveRDS(sro, file = paste0(pref.sro, "SRO.rds"))

## cellxgene ####
unimp.expr <- t(as.matrix(sro@assays$RNA$data))
unimp.expr.dgC <- as(unimp.expr, "sparseMatrix") 
adata <- AnnData(
  X = unimp.expr.dgC,
  obs = sro@meta.data,
  obsm = list(
    X_umap = sro@reductions$umap@cell.embeddings
  )
)
write_h5ad(adata, paste0(pref.sro, "unimputed-expr.h5ad"))