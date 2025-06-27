# Colonna
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

reticulate::py_require("anndata")
set.seed(1)
plan("multicore", workers = 16)
options(future.globals.maxSize = Inf)

source("utils.R")

pal <- list(
  Clusters = c(
    "#5fed0e", "#489de8", "#d40663", "#d1c50f", "#077315",
    "#e67109", "#785cd4", "#260691", "#9e7d3f", "#bd537a",
    "#49709c", "#aebeff", "#9c2903", "#9c6fa8", "#827c68",
    "#062e0b", "#1ee3c5"
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
dir.create('plots')

# 1. Load data ####
pref.sro <- 'Seurat/'; pref.p.sro <- 'plots/Seurat/'
dir.create(pref.sro); dir.create(pref.p.sro)

data.folder <- "data/"
samples <- c(hCD2_neg = "GSE289270_hCD2_neg", hCD2_pos = "GSE289270_hCD2_pos")

sro.list <- lapply(1:length(samples), function(i){
  conf.table <- read.csv(paste0("data/", samples[i], "_assignment_confidence_table.csv.gz"), row.names = 1)
  rownames(conf.table) <- conf.table$Barcode
  counts <- readMM(paste0("data/", samples[i], "_matrix.mtx.gz"))
  barcodes <- readLines(paste0("data/", samples[i], "_barcodes.tsv.gz"))
  feature.info <- read.table(paste0("data/", samples[i], "_features.tsv.gz"), 
                             sep = "\t", header = F,
                             stringsAsFactors = F, col.names = c("ID", "symbol", "type"))
  feature.info$symbol.unique <- make.unique(feature.info$symbol)
  rownames(feature.info) <- feature.info$symbol.unique
  rownames(counts) <- feature.info$symbol.unique; colnames(counts) <- barcodes
  
  counts <- counts[, conf.table$Barcode]
  sro <- CreateSeuratObject(
    project = names(samples)[i],
    counts = counts, meta.data =  conf.table
  )
  return(list(sro = sro))
})
sro <- merge(sro.list[[1]]$sro, y = sapply(sro.list,"[[",1)[2:length(samples)], 
             add.cell.ids = names(samples), project = "Colonna")
sro$sample <- sro$orig.ident
sro$tissue <- 'mLN'
sro$age <- '15days'
sro$Assignment <- stringr::str_split_i(sro$Assignment, "_", 1)
sro$Assignment <- ifelse(grepl("Total", sro$Assignment), yes = stringr::str_split_i(sro$Assignment, "-", 2),
                         no = sro$Assignment)

feature.info <- read.table(paste0("data/", samples[1], "_features.tsv.gz"), 
                           sep = "\t", header = F,
                           stringsAsFactors = F, col.names = c("ID", "symbol", "type"))
feature.info$symbol.unique <- make.unique(feature.info$symbol)
rownames(feature.info) <- feature.info$symbol.unique

sro@assays$RNA <- AddMetaData(sro@assays$RNA, metadata = feature.info[rownames(sro), ])

sro <- PercentageFeatureSet(sro, pattern = "^mt-", col.name = "MtFrac_RNA")

thr <- data.frame(nc.min = 1000, nf.min = 2500, nf.max = 8000, 
                  mp.max = 5)
cell.discard <- sro$nCount_RNA < thr$nc.min | 
  sro$nFeature_RNA < thr$nf.min | sro$nFeature_RNA > thr$nf.max |
  sro$MtFrac_RNA > thr$mp.max | sro$Assignment %in% c("Blank", "Multiplet", "Unassigned")
table(cell.discard)
# cell.discard
# FALSE  TRUE 
# 9816  4588 

Idents(sro) <- sro$sample
dir.create("plots/Seurat/QC", recursive = T)
pdf("plots/Seurat/QC/violin-4I-sample-RNA-QC.pdf", width = 8, height = 6)
plot.all.QC(sro, ident = "sample", thr = thr)
dev.off()

pdf("plots/Seurat/QC/violin-4I-HTO-RNA-QC.pdf", width = 15, height = 6)
plot.all.QC(sro, ident = "Assignment", thr = thr)
dev.off()

# sro <- sro[, !cell.discard] # not run here
write.csv(sro@meta.data, file = "Seurat/initial-meta-data-4I.csv")
saveRDS(sro, file = "Seurat/initial-SRO-4I.rds")

# 2. RNA analysis ####
pref.sro <- 'Seurat/merged-4I/'; pref.p.sro <- 'plots/Seurat/merged-4I/'
dir.create(pref.sro); dir.create(pref.p.sro)
sro <- readRDS(paste0('Seurat/initial-SRO-4I.rds'))
sro <- sro[, !cell.discard]

s.genes <- rownames(sro)[toupper(rownames(sro)) %in% toupper(cc.genes.updated.2019$s.genes)]
g2m.genes <- rownames(sro)[toupper(rownames(sro)) %in% toupper(cc.genes.updated.2019$g2m.genes)]
sro <- CellCycleScoring(sro, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
sro$Phase <- gsub("G2M", "G2/M", sro$Phase)

res.list <- seq(0.1, 1, 0.1)
cell.count <- rowSums(sro@assays$RNA@counts > 0)
sro <- NormalizeData(sro[cell.count > 1, ]) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 5000)
sro <- ScaleData(sro, features = rownames(sro)) %>%
  RunPCA(features = rownames(sro), npcs = 50) %>%
  FindNeighbors(dims = 1:30, k.param = 30) %>%
  FindClusters(resolution = res.list) %>%
  RunUMAP(dims = 1:30, n.neighbors = 30, metric = "cosine", min.dist = 0.4, spread = 1)


# res 0.4
# res0.4.to.anno <- c(7 = 'TC I', 5 = 'TC II', 9 = 'TC III', 10 = 'TC IV')
res0.4.to.anno <- c(7 = 'RORgt DC I', 5 = 'RORgt DC II', 9 = 'RORgt DC III', 10 = 'RORgt DC IV')

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
