Sys.setenv(RETICULATE_PYTHON = "/data1/lesliec/tyler/utils/miniforge3/envs/multiome/bin/python")
setwd("~/0-workspace/CCR7_DC/oral-tolerance-Gardner/")

suppressPackageStartupMessages({
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
saveRDS(pal, file = "plots/palette.rds")
pal <- readRDS(file = "plots/palette.rds")

adata.e <- read_h5ad("data/GSE285182_early_life_int.h5ad")
adata.a <- read_h5ad("data/GSE273746_named_int_td_adata_OG_FINAL.h5ad")

write.csv(adata.e$obs, "data/GSE285182_adata_obs.csv", quote = F)
write.csv(adata.a$obs, "data/GSE273746_adata_obs.csv", quote = F)

# 1. Load data ####

pref.sro <- 'Seurat/'; pref.p.sro <- 'plots/Seurat/'
dir.create(pref.sro); dir.create(pref.p.sro)

data.folder <- "data/"
samples <- c(adult_LN = "GSE273746_RAW-Adult-RLT/GSM8436371_tomato_ln", 
             adult_mLN = "GSE273746_RAW-Adult-RLT/GSM8436372_tomato_mln", 
             adult_spleen = "GSE273746_RAW-Adult-RLT/GSM8436373_tomato_sp", 
             early_LN = "GSE285182_RAW-Early-Life-RLT/GSM8697645_early_life_ln",
             early_mLN = "GSE285182_RAW-Early-Life-RLT/GSM8697646_early_life_mln",
             early_spleen = "GSE285182_RAW-Early-Life-RLT/GSM8697647_early_life_sp"
            )

sro.list <- lapply(1:length(samples), function(i){
  counts <- Read10X_h5(paste0(data.folder, samples[i], "_filtered_feature_bc_matrix.h5"))
  sro <- CreateSeuratObject(
    project = names(samples)[i],
    counts = counts
  )
  return(list(sro = sro))
})

sro.a <- merge(sro.list[[1]]$sro, y = sapply(sro.list,"[[",1)[2:3], 
             add.cell.ids = names(samples)[1:3], project = "Gardner.adult")

sro.e <- merge(sro.list[[4]]$sro, y = sapply(sro.list,"[[",1)[5:6], 
             add.cell.ids = names(samples)[4:6], project = "Gardner.early")

sro.a$sample <- sro.a$orig.ident
sro.a$tissue <- stringr::str_split_i(sro.a$sample, "_", i=2)
sro.a$age <- stringr::str_split_i(sro.a$sample, "_", i=1)
sro.a <- PercentageFeatureSet(sro.a, pattern = "^mt-", col.name = "MtFrac_RNA")

sro.e$sample <- sro.e$orig.ident
sro.e$tissue <- stringr::str_split_i(sro.e$sample, "_", i=2)
sro.e$age <- stringr::str_split_i(sro.e$sample, "_", i=1)
sro.e <- PercentageFeatureSet(sro.e, pattern = "^mt-", col.name = "MtFrac_RNA")

thr.a <- data.frame(nc.max = 60000, nf.max = 8000, 
                  mp.max = 5)

thr.e <- data.frame(nc.max = 35000, nf.max = 5500, 
                  mp.max = 5)

17975+13482+27642

14682+11428+22408

cell.discard.a <- sro.a$nCount_RNA > thr.a$nc.max |
  sro.a$nFeature_RNA > thr.a$nf.max |
  sro.a$MtFrac_RNA > thr.a$mp.max
table(cell.discard.a)

cell.discard.e <- sro.e$nCount_RNA > thr.e$nc.max |
  sro.e$nFeature_RNA > thr.e$nf.max |
  sro.e$MtFrac_RNA > thr.e$mp.max
table(cell.discard.e)

# cell.discard.a
# FALSE  TRUE 
# 38409  1676 
# cell.discard.e
# FALSE  TRUE 
# 52351  6748 

Idents(sro.a) <- sro.a$sample
dir.create("plots/Seurat/QC", recursive = T)
pdf("plots/Seurat/QC/violin-adult-sample-RNA-QC.pdf", width = 8, height = 6)
plot.all.QC(sro.a, ident = "sample", thr = thr.a)
dev.off()

Idents(sro.e) <- sro.e$sample
dir.create("plots/Seurat/QC", recursive = T)
pdf("plots/Seurat/QC/violin-early-sample-RNA-QC.pdf", width = 8, height = 6)
plot.all.QC(sro.e, ident = "sample", thr = thr.e)
dev.off()

# sro <- sro[, !cell.discard] # not run here
write.csv(sro.a@meta.data, file = "Seurat/adult-initial-meta-data.csv")
saveRDS(sro.a, file = "Seurat/adult-initial-SRO.rds")

write.csv(sro.e@meta.data, file = "Seurat/early-initial-meta-data.csv")
saveRDS(sro.e, file = "Seurat/early-initial-SRO.rds")

pref.sro.a <- 'Seurat/adult/'; pref.p.sro.a <- 'plots/Seurat/adult/'
dir.create(pref.sro.a); dir.create(pref.p.sro.a)
pref.sro.e <- 'Seurat/early/'; pref.p.sro.e <- 'plots/Seurat/early/'
dir.create(pref.sro.e); dir.create(pref.p.sro.e)

# # 2. RNA analysis ####
# sro.a <- readRDS(paste0('Seurat/adult-initial-SRO.rds'))
# sro.a <- sro.a[, !cell.discard.a]

# sro.e <- readRDS(paste0('Seurat/early-initial-SRO.rds'))
# sro.e <- sro.e[, !cell.discard.e]

obs.a <- read.csv("data/GSE273746_adata_obs.csv", row.names = 1)
obs.e <- read.csv("data/GSE285182_adata_obs.csv", row.names = 1)

sro.a$barcode <- stringr::str_split_i(rownames(sro.a@meta.data), "_", 3)
sro.a$paper.annotations <- obs.a[sro.a$barcode, ]$final_cell_type

table(sro.a$barcode %in% rownames(obs.a))

sro.e$barcode <- stringr::str_split_i(rownames(sro.e@meta.data), "_", 3)
sro.e$paper.clusters <- obs.e[sro.e$barcode, ]$leiden

table(sro.e$barcode %in% rownames(obs.e))

table(sro.e$paper.clusters, useNA='ifany')

# sro.input <- sro.a
# pref.sro.input <- pref.sro.a
# pref.p.sro.input <- pref.p.sro.a

sro.input <- sro.e
pref.sro.input <- pref.sro.e
pref.p.sro.input <- pref.p.sro.e

s.genes <- rownames(sro.input)[toupper(rownames(sro.input)) %in% toupper(cc.genes.updated.2019$s.genes)]
g2m.genes <- rownames(sro.input)[toupper(rownames(sro.input)) %in% toupper(cc.genes.updated.2019$g2m.genes)]
sro.input <- CellCycleScoring(sro.input, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
sro.input$Phase <- gsub("G2M", "G2/M", sro.input$Phase)

res.list <- seq(0.2, 2, 0.2)
cell.count <- rowSums(sro.input@assays$RNA@counts > 0)
sro.input <- NormalizeData(sro.input[cell.count > 1, ]) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 5000)
sro.input <- ScaleData(sro.input, features = rownames(sro.input)) %>%
  RunPCA(features = rownames(sro.input), npcs = 50) %>%
  FindNeighbors(dims = 1:30, k.param = 30) %>%
  FindClusters(resolution = res.list) %>%
  RunUMAP(dims = 1:30, n.neighbors = 30, metric = "cosine", min.dist = 0.4, spread = 1)

## save res ####
write.csv(sro.input@reductions$pca@cell.embeddings, file = paste0(pref.sro.input, "PCA.csv"), quote = F)
write.csv(sro.input@reductions$umap@cell.embeddings, file = paste0(pref.sro.input, "UMAP.csv"), quote = F)
write.csv(sro.input@meta.data, file = paste0(pref.sro.input, "meta-data.csv"), quote = F)
saveRDS(sro.input, file = paste0(pref.sro.input, "SRO.rds"))

## cellxgene ####
unimp.expr <- t(as.matrix(sro.input@assays$RNA$data))
unimp.expr.dgC <- as(unimp.expr, "sparseMatrix") 
adata <- AnnData(
  X = unimp.expr.dgC,
  obs = sro.input@meta.data,
  obsm = list(
    X_umap = sro.input@reductions$umap@cell.embeddings
  )
)
write_h5ad(adata, paste0(pref.sro.input, "unimputed-expr.h5ad"))


# 2. RNA analysis ####
sro.a <- readRDS(paste0('Seurat/adult/SRO.rds'))

sro.input <- sro.a
pref.sro.input <- pref.sro.a
pref.p.sro.input <- pref.p.sro.a

pdf(paste0(pref.p.sro.input, "QC/violin-paper.annotations-RNA-QC.pdf"), width = 16, height = 6)
plot.all.QC(sro.input, ident = "paper.annotations")
dev.off()

## plots ####
dir.create(paste0(pref.p.sro.input, 'QC/'), recursive = T)

for (res in res.list){
  colname <- paste0("RNA_snn_res.", res)
  pl <- plot.all.QC(sro.input, ident = colname, col = pal$Clusters_long)
  pdf(paste0(pref.p.sro.input, "QC/violin-QC-", colname, ".pdf"), width = 12, height = 6)
  for (i in 1:length(pl)){
    print(pl[[i]])
  }
  dev.off()
}

pdf(paste0(pref.p.sro.input, "QC/violin-sample-RNA-QC.pdf"), width = 8, height = 6)
plot.all.QC(sro.input, ident = "sample")
dev.off()

pdf(paste0(pref.p.sro.input, "QC/UMAP-QC.pdf"), width = 12, height = 10)
plot.continuous.value(sro.input, idx = rownames(sro.input@meta.data), vis = sro.input@reductions$umap@cell.embeddings,
                      val = sro.input$nCount_RNA, val.name='nCount_RNA', point.size=1)
plot.continuous.value(sro.input, idx = rownames(sro.input@meta.data), vis = sro.input@reductions$umap@cell.embeddings,
                      val = sro.input$nFeature_RNA, val.name='nFeature_RNA', point.size=1)
plot.continuous.value(sro.input, idx = rownames(sro.input@meta.data), vis = sro.input@reductions$umap@cell.embeddings,
                      val = sro.input$MtFrac_RNA, val.name='MtFrac_RNA', point.size=1)
plot.continuous.value(sro.input, idx = rownames(sro.input@meta.data), vis = sro.input@reductions$umap@cell.embeddings,
                      val = sro.input$S.Score, val.name='S.Score', point.size=1)
plot.continuous.value(sro.input, idx = rownames(sro.input@meta.data), vis = sro.input@reductions$umap@cell.embeddings,
                      val = sro.input$G2M.Score, val.name='G2M.Score', point.size=1)
dev.off()

pdf(paste0(pref.p.sro.input, "UMAP-Clusters.pdf"), width = 12, height = 10)
for (res in res.list){
  colname <- paste0("RNA_snn_res.", res)
  print(
    plot.clusters(sro.input, groups = sro.input@meta.data[[colname]],
                  clusters.col = colname, col = pal$Clusters_long,
                  label.size = 5, point.size = 0.5,
                  pref.C = T)
  )
}
dev.off()

pdf(paste0(pref.p.sro.input, "UMAP-sample.pdf"), width = 12, height = 10)
plot.clusters(sro.input, groups = sro.input$sample, clusters.col = "sample",
              col = pal$Clusters, label.size = 5, labels = F, point.size = 0.5,
              label.pad = 1, pref.C = F)
dev.off()

pdf(paste0(pref.p.sro.input, "UMAP-age.pdf"), width = 12, height = 10)
plot.clusters(sro.input, groups = sro.input$age, clusters.col = "age",
              col = pal$Clusters, label.size = 5, labels = F, point.size = 0.5,
              label.pad = 1, pref.C = F)
dev.off()

pdf(paste0(pref.p.sro.input, "UMAP-tissue.pdf"), width = 12, height = 10)
plot.clusters(sro.input, groups = sro.input$tissue, clusters.col = "tissue",
              col = pal$Clusters, label.size = 5, labels = F, point.size = 0.5,
              label.pad = 1, pref.C = F)
dev.off()


# for adult
# pdf(paste0(pref.p.sro.input, "UMAP-paper.annotations.pdf"), width = 12, height = 10)
# plot.clusters(sro.input, groups = sro.input$paper.annotations, clusters.col = "paper.annotations",
#               col = pal$Clusters, label.size = 5, labels = T, point.size = 0.5,
#               label.pad = 1, pref.C = F)
# dev.off()


# for child
pdf(paste0(pref.p.sro.input, "UMAP-paper.clusters.pdf"), width = 12, height = 10)
plot.clusters(sro.input, groups = sro.input$paper.clusters, clusters.col = "paper.clusters",
              col = pal$Clusters, label.size = 5, labels = T, point.size = 0.5,
              label.pad = 1, pref.C = F)
dev.off()


FindSubCluster.custom <- function(sro, sro.subset, master.res, 
                                  newcolname.pref,
                                  subres.list = seq(0.1, 0.5, 0.1)
                                  ){
  cell.count <- rowSums(sro.subset@assays$RNA@counts > 0)
  sro.subset <- NormalizeData(sro.subset[cell.count > 1, ]) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 5000)
  
  sro.subset <- ScaleData(sro.subset, features = rownames(sro.subset)) %>%
    RunPCA(features = rownames(sro.subset), npcs = 50) %>%
    FindNeighbors(dims = 1:30, k.param = 30) %>%
    FindClusters(resolution = subres.list) %>%
    RunUMAP(dims = 1:30, n.neighbors = 30, metric = "cosine", min.dist = 0.4, spread = 1)
  
  for (res in subres.list){
    colname.in.subset <- paste0("RNA_snn_res.", res)
    newcolname <- paste0(newcolname.pref, res)
    subclusters <- ifelse(test = is.na(sro.subset@meta.data[colnames(sro), colname.in.subset]),
                          yes = sro.subset@meta.data[colnames(sro), colname.in.subset],
                          no = paste0("sub", sro.subset@meta.data[colnames(sro), colname.in.subset])
    )
    sro@meta.data[[newcolname]] <- coalesce(subclusters, sro@meta.data[[master.res]])
  }
  return(sro)
}


pref.sro.a <- 'Seurat/adult/'; pref.p.sro.a <- 'plots/Seurat/adult/'
pref.sro.e <- 'Seurat/early/'; pref.p.sro.e <- 'plots/Seurat/early/'

sro.e <- readRDS(paste0(pref.sro.e, 'SRO.rds'))

sro.input <- sro.e
pref.sro.input <- pref.sro.e
pref.p.sro.input <- pref.p.sro.e

subres.list <- seq(0.1, 0.5, 0.1)
sro.subset <- subset(sro, RNA_snn_res.0.2 == 7)
sro.input <- FindSubCluster.custom(sro.input, sro.subset, 'RNA_snn_res.0.2', 'res.0.2_C7_subres.', subres.list)

raster_pdf(paste0(pref.p.sro.input, "UMAP-RNA_snn_res.0.2_subC7.pdf"), width = 12, height = 10, res = 300)
for (subres in subres.list){
    group.name <- paste0("res.0.2_C7_subres.", subres)
    print(plot.clusters(sro.input, groups = sro.input@meta.data[[group.name]], clusters.col = group.name,
              col = pal$Clusters_long, label.size = 5, labels = T,
              point.size = 1, label.pad = 1, pref.C = T))
}

dev.off()

table(sro.input$`res.0.2_C7_subres.0.1`)

cl.to.anno <- c(
    '1' = 'ILC3', 
    '10' = 'ILC3',
    'sub0' = 'eTAC I',
    'sub1' = 'eTAC III',
    'sub2' = 'eTAC II',
    'sub3' = 'Proliferating eTAC'
)

head(sro$res.0.2_C7_subres.0.1)

sro.input$paper.annot <- ifelse(
    sro.input$`res.0.2_C7_subres.0.1` %in% names(cl.to.anno),
    yes = cl.to.anno[sro.input$`res.0.2_C7_subres.0.1`],
    no = 'na'
)

raster_pdf(paste0(pref.p.sro.input, "UMAP-paper.annot.pdf"), width = 12, height = 10, res = 300)
group.name <- "paper.annot"
print(plot.clusters(sro.input, groups = sro.input@meta.data[[group.name]], clusters.col = group.name,
          col = pal$Clusters, label.size = 5, labels = T,
          point.size = 1, label.pad = 1, pref.C = F))
dev.off()

## save res ####
write.csv(sro.input@meta.data, file = paste0(pref.sro.input, "meta-data.csv"), quote = F)
saveRDS(sro.input, file = paste0(pref.sro.input, "SRO.rds"))



cl.to.annot <- c(
  "0" = "LTi",
  "1" = "LTi",
  "2" = "LTi",
  "3" = "LTi",
  "4" = "LTi",
  "5" = "LTi-like ILC",
  "6" = "R-eTAC1",
  "7" = "LTi",
  "9" = "R-eTAC2",
  "10" = "LTi",
  "11" = "R-cDC2",
  "12" = "R-eTAC3",
  "14" = "R-cDC1",
  "15" = "LTi",
  "16" = "R-mDC",
  "17" = "LTi",
  "18" = "LTi"
)


sro.input$Cluster.annot <- sro.input$paper.annot

sro.input$paper.annot <- as.character(cl.to.annot[as.character(sro.input$paper.clusters)])

# for child
pdf(paste0(pref.p.sro.input, "UMAP-paper.clusters.pdf"), width = 12, height = 10)
plot.clusters(sro.input, groups = sro.input$paper.clusters, clusters.col = "paper.clusters",
              col = pal$Clusters, label.size = 5, labels = T, point.size = 0.5,
              label.pad = 1, pref.C = F)
dev.off()


# for child
pdf(paste0(pref.p.sro.input, "UMAP-Cluster.annot.pdf"), width = 12, height = 10)
plot.clusters(sro.input, groups = sro.input$Cluster.annot, clusters.col = "Cluster.annot",
              col = pal$Clusters, label.size = 5, labels = T, point.size = 0.5,
              label.pad = 1, pref.C = F)
dev.off()


# for child
pdf(paste0(pref.p.sro.input, "UMAP-paper.annot.pdf"), width = 12, height = 10)
plot.clusters(sro.input, groups = sro.input$paper.annot, clusters.col = "paper.annot",
              col = pal$Clusters, label.size = 5, labels = T, point.size = 0.5,
              label.pad = 1, pref.C = F)
dev.off()


## save res ####
write.csv(sro.input@meta.data, file = paste0(pref.sro.input, "meta-data.csv"), quote = F)
saveRDS(sro.input, file = paste0(pref.sro.input, "SRO.rds"))


