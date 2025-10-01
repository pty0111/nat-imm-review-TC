# integrate all datasets
Sys.setenv(RETICULATE_PYTHON = "/data1/lesliec/tyler/utils/miniforge3/envs/multiome/bin/python")
setwd("~/0-workspace/nat-imm-review-TC/")

suppressPackageStartupMessages({
  library(curl)
  library(Seurat)
  library(anndata)
  library(Rmagic)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
  library(ComplexHeatmap)
  library(circlize)
  library(rasterpdf)
  library(future)
  library(pals)
})
reticulate::py_require("anndata")

set.seed(1)
options(future.globals.maxSize = Inf)
plan(strategy = "multicore", workers = 16)

pal <- readRDS("plots/palette.rds")
source("utils.R")

# ############################################################################ #
# 1. Load SROs ####
# ############################################################################ #
## MLN_RORgt ####
sro1 <- readRDS("../GSE205066-MLN_RORgt_MHCII_multiome-Brown/Seurat/results/SRO.rds")
sro1$MtFrac_RNA <- sro1$percent.MT

sro1@assays$RNA@meta.features$symbol.unique <- make.unique(sro1@assays$RNA@meta.features$Symbol)
new.genes <- toupper(sro1@assays$RNA@meta.features$symbol.unique)
sro1 <- RenameGenesSeurat(sro1, new.genes)
sro1$sample <- sro1$orig.ident

sro1 <- sro1 %>% 
  subset(Clusters %in% 1:16) %>%
  NormalizeData() %>% FindVariableFeatures(nfeatures = 5000) %>%
  ScaleData() %>% RunPCA(npcs = 50)
sro1$orig.ident <- "MLN_RORgt_MHCII_multiome"
sro1$Cluster.prev <- paste0("mto C", sro1$Clusters)
sro1$Cluster.annot <- sro1$Clusters.annot
sro1$paper.annot <- sro1$Cluster.annot

Idents(sro1) <- sro1$Cluster.prev

## Kedmi ####
sro2 <- readRDS("../GSE200148-kedmi-2022/results/SRO.rds")
sro2$MtFrac_RNA <- sro2$percent.MT
sro2$sample <- sro2$Sample

sro2 <- sro2 %>% 
  subset(RNA_snn_res.0.6 %in% c(14, 6, 12, 9)) %>%
  NormalizeData() %>% FindVariableFeatures(nfeatures = 5000) %>%
  ScaleData() %>% RunPCA(npcs = 50)

sro2$RNA_snn_res.0.6 <- as.character(sro2$RNA_snn_res.0.6)
cl2annot.kedmi <- c('14' = 'ILC3', '6' = 'JC', '12' = 'JC', '9' = 'JC')
sro2$Cluster.annot <- as.character(cl2annot.kedmi[sro2$RNA_snn_res.0.6])

sro2$orig.ident <- "Kedmi"
sro2$Cluster.prev <- paste0("Kedmi C", sro2$RNA_snn_res.0.6)
sro2$paper.annot <- sro2$Cluster.annot

Idents(sro2) <- sro2$Cluster.prev

## Wang ####
sro3 <- readRDS("../GSE176282-wang-2021/results/SRO.rds")
sro3$MtFrac_RNA <- sro3$percent.MT
sro3$sample <- sro3$Sample

sro3 <- sro3 %>% 
  subset(RNA_snn_res.0.6 %in% c(3, 12, 9, 13)) %>%
  NormalizeData() %>% FindVariableFeatures(nfeatures = 5000) %>%
  ScaleData() %>% RunPCA(npcs = 50)

sro3$RNA_snn_res.0.6 <- as.character(sro3$RNA_snn_res.0.6)
cl2annot.wang <- c('3' = 'TC I', '12' = 'TC II, III, IV', '9' = 'Ki67 TC', '13' = 'LTi')
sro3$Cluster.annot <- as.character(cl2annot.wang[sro3$RNA_snn_res.0.6])

sro3$orig.ident <- "Wang"
sro3$Cluster.prev <- paste0("Wang C", sro3$RNA_snn_res.0.6)
sro3$paper.annot <- sro3$Cluster.annot

Idents(sro3) <- sro3$Cluster.prev

## Lyu ####
sro4 <- readRDS("../GSE184175-lyu-2022/results/SRO.rds")
s.genes <- rownames(sro4)[toupper(rownames(sro4)) %in% toupper(cc.genes.updated.2019$s.genes)]
g2m.genes <- rownames(sro4)[toupper(rownames(sro4)) %in% toupper(cc.genes.updated.2019$g2m.genes)]
sro4 <- CellCycleScoring(sro4, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
sro4$Phase <- gsub("G2M", "G2/M", sro4$Phase)
sro4$MtFrac_RNA <- sro4$percent.MT
sro4$sample <- sro4$Sample

sro4 <- sro4 %>% 
  subset(Clusters %in% c(2, 6, 11)) %>%
  NormalizeData() %>% FindVariableFeatures(nfeatures = 5000) %>%
  ScaleData() %>% RunPCA(npcs = 50)
sro4$orig.ident <- "Lyu"
sro4$Cluster.prev <- paste0("Lyu C", sro4$Clusters)
sro4$Cluster.annot <- sro4$annotations
sro4$paper.annot <- sro4$Cluster.annot

Idents(sro4) <- sro4$Cluster.prev

## TC_all_LN ####
sro5 <- readRDS("../GSE294005-TC_all_LN-Brown/results/SRO.rds")
sro5$MtFrac_RNA <- sro5$percent.MT
sro5 <- sro5 %>% subset(Clusters2 %in% c(1:7)) %>%
  NormalizeData() %>% FindVariableFeatures(nfeatures = 5000) %>%
  ScaleData() %>% RunPCA(npcs = 50)
sro5$orig.ident <- "TC_all_LN"
sro5$Cluster.prev <- paste0("TC_all_LN C", sro5$Clusters)
sro5$Cluster.annot <- sro5$Annotation
sro5$paper.annot <- sro5$Cluster.annot

Idents(sro5) <- sro5$Cluster.prev

## oral tolerance Colonna ####
sro6 <- readRDS("../GSE289268-oral-tolerance-Colonna/Seurat/merged-4I/SRO.rds")
sro6 <- RenameGenesSeurat(sro6, toupper(rownames(sro6)))

cl.annot <- c('5' = 'TC II', '7' = 'TC I', '9' = 'TC III', '10' = 'TC IV', '0' = 'ILC3')
sro6$annotations <- as.character(cl.annot[as.character(sro6$RNA_snn_res.0.4)])
sro6$annotations[is.na(sro6$annotations)] <- 'na'

sro6 <- sro6 %>% subset(RNA_snn_res.0.4 %in% c(5, 7, 9, 10, 0))
sro6 <- sro6 %>% 
  NormalizeData() %>% FindVariableFeatures(nfeatures = 5000) %>%
  ScaleData() %>% RunPCA(npcs = 50)
sro6$orig.ident <- "Colonna"
sro6$Cluster.prev <- paste0("Colonna C", sro6$RNA_snn_res.0.4)
sro6$Cluster.annot <- sro6$annotations
sro6$paper.annot <- sro6$Cluster.annot

Idents(sro6) <- sro6$Cluster.prev

## oral tolerance Littman ####
sro7 <- readRDS("../zenodo-15212000-oral-tolerance-Littman-oral-tolerance-Littman/Seurat/merged/SRO.rds")
sro7 <- RenameGenesSeurat(sro7, toupper(rownames(sro7)))
sro7$RabiClusters[is.na(sro7$RabiClusters)] <- 'na'

cl.annot <- c('14' = 'TCs', '1' = 'ILC3', '5' = 'ILC3')
sro7$annotations <- as.character(cl.annot[as.character(sro7$RNA_snn_res.0.6)])
sro7$annotations[is.na(sro7$annotations)] <- 'na'

sro7 <- sro7 %>% subset(RNA_snn_res.0.6 %in% c(14, 1))
sro7 <- sro7 %>% 
  NormalizeData() %>% FindVariableFeatures(nfeatures = 5000) %>%
  ScaleData() %>% RunPCA(npcs = 50)
sro7$orig.ident <- "Littman"
sro7$Cluster.prev <- paste0("Littman C", sro7$RNA_snn_res.0.6)
sro7$Cluster.annot <- sro7$annotations
sro7$paper.annot <- sro7$RabiClusters

Idents(sro7) <- sro7$Cluster.prev

## oral tolerance Gardner ####
### adult
sro8 <- readRDS("../GSE273746-GSE285182-oral-tolerance-Gardner-oral-tolerance-Gardner/Seurat/adult/SRO.rds")
sro8 <- RenameGenesSeurat(sro8, toupper(rownames(sro8)))
sro8$paper.annotations[is.na(sro8$paper.annotations)] <- 'na'

cl.annot <- c('1' = 'ILC3', '3' = 'TCs', '4' = 'ILC3')
sro8$annotations <- as.character(cl.annot[as.character(sro8$RNA_snn_res.0.2)])
sro8$annotations[is.na(sro8$annotations)] <- 'na'

sro8 <- sro8 %>% subset(RNA_snn_res.0.2 %in% c(1, 3, 4))
sro8 <- sro8 %>% 
  NormalizeData() %>% FindVariableFeatures(nfeatures = 5000) %>%
  ScaleData() %>% RunPCA(npcs = 50)
sro8$orig.ident <- "Gardner.A"
sro8$Cluster.prev <- paste0("Gardner.A C", sro8$RNA_snn_res.0.2)
sro8$Cluster.annot <- sro8$annotations
sro8$paper.annot <- sro8$paper.annotations

Idents(sro8) <- sro8$Cluster.prev

### early life
sro9 <- readRDS("../GSE273746-GSE285182-oral-tolerance-Gardner-oral-tolerance-Gardner/Seurat/early/SRO.rds")
sro9 <- RenameGenesSeurat(sro9, toupper(rownames(sro9)))

cl.annot <- c('1' = 'ILC3', '7' = 'TCs', '10' = 'ILC3')
sro9$annotations <- as.character(cl.annot[as.character(sro9$RNA_snn_res.0.2)])
sro9$annotations[is.na(sro9$annotations)] <- 'na'

sro9 <- sro9 %>% subset(RNA_snn_res.0.2 %in% c(1, 7))
table(sro9$annotations)
sro9 <- sro9 %>% 
  NormalizeData() %>% FindVariableFeatures(nfeatures = 5000) %>%
  ScaleData() %>% RunPCA(npcs = 50)
sro9$orig.ident <- "Gardner.E"
sro9$Cluster.prev <- paste0("Gardner.E C", sro9$RNA_snn_res.0.2)
sro9$Cluster.annot <- sro9$annotations
sro9$paper.annot <- sro9$paper.annot

Idents(sro9) <- sro9$Cluster.prev

## clean up Meta data ####
md.columns <- c('orig.ident', 'sample', 'nCount_RNA', 'nFeature_RNA', 'MtFrac_RNA',
                'S.Score', 'G2M.Score', 'Phase', 'Cluster.prev', 'Cluster.annot', 'paper.annot')
sro1@meta.data <- sro1@meta.data[md.columns]
sro2@meta.data <- sro2@meta.data[md.columns]
sro3@meta.data <- sro3@meta.data[md.columns]
sro4@meta.data <- sro4@meta.data[md.columns]
sro5@meta.data <- sro5@meta.data[md.columns]
sro6@meta.data <- sro6@meta.data[md.columns]
sro7@meta.data <- sro7@meta.data[md.columns]
sro8@meta.data <- sro8@meta.data[md.columns]
sro9@meta.data <- sro9@meta.data[md.columns]

# ############################################################################ #
# 2. Integrate ####
# ############################################################################ #
## Prep integration ####
sro.list <- list(sro1, sro2, sro3, sro4, sro5, sro6, sro7, sro8, sro9)

genes <- Reduce(intersect, list(rownames(sro1), rownames(sro2), rownames(sro3),
                                rownames(sro4), rownames(sro5), rownames(sro6),
                                rownames(sro7), rownames(sro8), rownames(sro9)))
genes.info <- sro1@assays$RNA@meta.features[genes, ]

pref.i <- "integrate-ILC3-TC-v5/"; dir.create(pref.i)
pref <- paste0(pref.i, "cca/"); dir.create(pref)

write.csv(genes.info, file = paste0(pref.i, "gene-info.csv"))

## Run integration ####
intg.anchors <- FindIntegrationAnchors(
  object.list = sro.list,
  normalization.method = "LogNormalize",
  anchor.features = 5000,
  reduction = reduction
)
sro.i <- IntegrateData(
  anchorset = intg.anchors,
  normalization.method = "LogNormalize",
  features.to.integrate = rownames(genes.info)
)

sro.i <- ScaleData(sro.i, features = rownames(sro.i)) %>%
  RunPCA(features = rownames(sro.i), npcs = 50) %>%
  FindNeighbors(dims = 1:30, k.param = 30) %>%
  RunUMAP(dims = 1:30, n.neighbors = 30, metric = "cosine", min.dist = 0.4, spread = 1)

sro.i$MtFrac_RNA_quant <- cut(sro.i$MtFrac_RNA,
                              breaks=c(-Inf, summary(sro.i$MtFrac_RNA)['1st Qu.'], summary(sro.i$MtFrac_RNA)['Median'],
                                       summary(sro.i$MtFrac_RNA)['3rd Qu.'], Inf),
                              labels=c("0~25%","25~50%","50~75%", "75~100%"))

## Final annotations ####
kedmi.wang.cl.to.annot <- c(
    'Kedmi C14' = 'LTi-like ILC3',
    'Kedmi C12' = 'Ki67+ JC',
    'Kedmi C6' = 'JC1',
    'Kedmi C9' = 'JC2',
    'Wang C13' = 'LTi-like ILC3',
    'Wang C9' = 'Ki67+ JC',
    'Wang C3' = 'JC1',
    'Wang C12' = 'JC2'
)

sro.i$paper.annot <- ifelse(
    test = sro.i$Cluster.prev %in% names(kedmi.wang.cl.to.annot),
    yes = as.character(kedmi.wang.cl.to.annot[sro.i$Cluster.prev]),
    no = sro.i$paper.annot
)

relabel.annot <- c(
    ILC = "ILC3", ILC1 = "ILC3", ILC2 = "ILC3", ILC3 = "ILC3", ILC3p = "ILC3",
    JC1 = "JC 1", JC2 = "JC 2", "Ki67+ JC" = "Ki67+", "Ki67+ TC" = "Ki67+ TC",
    LTi = "ILC3", "LTi Variation 1" = "ILC3", "LTi Variation 2" = "ILC3",
    "LTi Variation 3" = "ILC3", "LTi Variation 4" = "ILC3", "LTi Variation 5" = "ILC3",
    "LTi Variation 6" = "ILC3", "LTi Variation 7" = "ILC3", "LTi Variation 8" = "ILC3",
    "LTi-like ILC" = "ILC3", "LTi-like ILC3" = "ILC3", "LTi-like ILC3s" = "ILC3",
    "LTi-like R-eTAC" = "LTi-like eTAC",
    Mig_DC_1 = "other", "NCR+ ILC3" = "ILC3", Nrg1_Pos = "Nrg1+ / TC I",
    "Pro. ILC" = "ILC3", "Pro. R-eTAC" = "Ki67+ eTAC",
    "Proliferating NCR+ ILC3" = "ILC3", "Proliferating TC" = "Ki67+ TC",
    "R-cDC1" = "other", "R-cDC2" = "other",
    "R-eTAC1" = "R-eTAC1", "R-eTAC2" = "R-eTAC2", "R-eTAC3" = "R-eTAC3",
    "R-mDC" = "other",
    "RORgt DC I" = "RORgt DC I", "RORgt DC II" = "RORgt DC II",
    "RORgt DC III" = "RORgt DC III", "RORgt DC IV" = "RORgt DC IV",
    "T cell zone macs" = "other",
    "TC 1" = "TC I", "TC 2" = "TC II", "TC 3" = "TC III", "TC 4" = "TC IV",
    "TC I" = "TC I", "TC II" = "TC II", "TC III" = "TC III", "TC IV" = "TC IV",
    "Tingible body macs" = "other", "Tolerogenic DC" = "TolDC", cDC2A = "other",
    "eTACs I" = "eTAC I", "eTACs II" = "eTAC II", "early/transitional TC" = "early TC",
    "na" = "other"
)

sro.i$annotations <- as.character(relabel.annot[sro.i@meta.data$paper.annot])
sro.i$annotations <- factor(sro.i$annotations, levels = names(pal$annotations))

## Save results ####
write.csv(sro.i@reductions$pca@cell.embeddings, file = paste0(pref, "PCA.csv"))
write.csv(sro.i@reductions$umap@cell.embeddings, file = paste0(pref, "UMAP.csv"))
write.csv(sro.i@meta.data, file = paste0(pref, "meta-data.csv"))
saveRDS(sro.i, paste0(pref, "SRO.rds"))

## h5ad ####
unimp.expr <- t(as.matrix(sro.i@assays$RNA$data))
unimp.expr.dgC <- as(unimp.expr, "sparseMatrix") 
adata <- AnnData(
  X = unimp.expr.dgC,
  obs = sro.i@meta.data,
  obsm = list(
    X_umap = sro.i@reductions$umap@cell.embeddings
  )
)
write_h5ad(adata, paste0(pref, "unimputed-expr.h5ad"))

# 3. plots ####
pref.i <- "integrate-ILC3-TC-v5/"; dir.create(pref.i)
reduction <- "cca"; pref <- paste0(pref.i, reduction, "/"); dir.create(pref)
pref.p <- "plots/"; dir.create(pref.p, recursive = T)

sro.i <- readRDS(paste0(pref, "SRO.rds"))

## QC ####
# PCA QC
raster_pdf(paste0(pref.p, 'pca-QC.pdf'), width = 15, height = 12, res = 150)
plot.continuous.value(sro.i, idx = rownames(sro.i@meta.data), vis = sro.i@reductions$pca@cell.embeddings,
                      val = sro.i$nCount_RNA, val.name='nCount_RNA', point.size=1)
plot.continuous.value(sro.i, idx = rownames(sro.i@meta.data), vis = sro.i@reductions$pca@cell.embeddings,
                      val = sro.i$nFeature_RNA, val.name='nFeature_RNA', point.size=1)
plot.continuous.value(sro.i, idx = rownames(sro.i@meta.data), vis = sro.i@reductions$pca@cell.embeddings,
                      val = sro.i$MtFrac_RNA, val.name='MtFrac_RNA', point.size=1)
plot.continuous.value(sro.i, idx = rownames(sro.i@meta.data), vis = sro.i@reductions$pca@cell.embeddings,
                      val = sro.i$S.Score, val.name='S.Score', point.size=1)
plot.continuous.value(sro.i, idx = rownames(sro.i@meta.data), vis = sro.i@reductions$pca@cell.embeddings,
                      val = sro.i$G2M.Score, val.name='G2M.score', point.size=1)
plot.clusters(sro.i, groups = sro.i$MtFrac_RNA_quant, vis = sro.i@reductions$pca@cell.embeddings,
              axis.titles = c("PC1", "PC2"),
              clusters.col = "MtFrac_RNA_quant", labels = F,
              label.size = 5,
              label.pad = 1, pref.C = F)
plot.clusters(sro.i, groups = sro.i$Phase, vis = sro.i@reductions$pca@cell.embeddings,
              axis.titles = c("PC1", "PC2"),
              clusters.col = "Phase", labels = F,
              label.size = 5,
              label.pad = 1, pref.C = F)
dev.off()

raster_pdf(paste0(pref.p, "UMAP-QC.pdf"), width = 15, height = 12, res = 150)
plot.continuous.value(sro.i, idx = colnames(sro.i), vis=sro.i@reductions$umap@cell.embeddings,
                      val = sro.i$nCount_RNA, val.name='nCount_RNA', point.size=1)
plot.continuous.value(sro.i, idx = colnames(sro.i), vis=sro.i@reductions$umap@cell.embeddings,
                      val = sro.i$nFeature_RNA, val.name='nFeature_RNA', point.size=1)
plot.continuous.value(sro.i, idx = colnames(sro.i), vis=sro.i@reductions$umap@cell.embeddings,
                      val = sro.i$MtFrac_RNA, val.name='MtFrac_RNA', point.size=1)
plot.continuous.value(sro.i, idx = colnames(sro.i), vis=sro.i@reductions$umap@cell.embeddings,
                      val = sro.i$S.Score, val.name='S.Score', point.size=1)
plot.continuous.value(sro.i, idx = colnames(sro.i), vis=sro.i@reductions$umap@cell.embeddings,
                      val = sro.i$G2M.Score, val.name='G2M.Score', point.size=1)
plot.clusters(sro.i, groups = sro.i$MtFrac_RNA_quant, vis = sro.i@reductions$umap@cell.embeddings,
              clusters.col = "MtFrac_RNA_quant", labels = F,
              label.size = 5,
              label.pad = 1, pref.C = F)
plot.clusters(sro.i, groups = sro.i$Phase, vis = sro.i@reductions$umap@cell.embeddings,
              clusters.col = "Phase", labels = F,
              label.size = 5,
              label.pad = 1, pref.C = F)
dev.off()

## group ####
sources <- c(
    "MLN_RORgt_MHCII_multiome", "TC_all_LN", "Colonna",
    "Kedmi", "Wang", "Littman", 
    "Lyu", "Gardner.A", "Gardner.E"
)

titles <- c(
    "Akagbosu et al. 2022 - 2w", "Cabric et al. 2025 - 2w", "Rodrigues et al. 2025 - 2w",
    "Kedmi et al. 2022 - 4w-adult", "Wang et al. 2021 - 4w", "Fu et al. 2025 - 3w",
    "Lyu et al. 2022 - adult", "Sun et al. 2025 - adult", "Sun et al. 2025 - 2w"
)

## annotations
group.name <- "annotations"
pl <- list()
for (i in 1:length(sources)){
    source.curr <- sources[i]
    title <- titles[i]
    pal.curr <- if (source.curr == "Lyu") pal$lyu else pal[[group.name]]
    p.curr <- plot.clusters.highlight.one(
        SRO = sro.i,
        idx1 = colnames(sro.i), idx2 = sro.i$orig.ident == source.curr,
        groups = get.named.vector.sro(sro.i, group.name),
        clusters.col = source.curr,
        labels = F, label.size = 5,
        col = pal.curr, pref.C = F) + guides(colour = guide_legend(title = NULL,
                                                                   label.theme = element_text(size = 14), 
                                                                   override.aes = list(size=4), ncol = 1)) +
        theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.line = element_blank(),
              plot.title = element_text(size = 21, hjust = 0.5, vjust = -2.5), legend.text.align = 0,
              legend.position = 'right', legend.box.margin = margin(l = -30, t = 0, b = 0, r = 0)) + ggtitle(title)
    
    if (source.curr == "Colonna"){
        p.curr <- p.curr + scale_color_manual(
            values = c("RORgt DC I" = stepped2()[16], "RORgt DC II" = stepped2()[15], "RORgt DC III" = stepped2()[14], "RORgt DC IV" = stepped2()[13], 
                       'ILC3' = "#14a38e"),
            labels = c(expression("ROR" * gamma * "t DC I"), expression("ROR" * gamma * "t DC II"), 
                       expression("ROR" * gamma * "t DC III"), expression("ROR" * gamma * "t DC IV"),
                       'ILC3')
        )
    }
    
    
    pl[[i]] <- p.curr
}

ncols <- 3
nrows <- ceiling(length(pl)/ncols)

raster_pdf(file = paste0(pref.p, "UMAP-annotations.pdf"), width = 7.5*ncols, height = 5*nrows, res = 300)
plot_grid(plotlist = pl, ncol = ncols, align = 'hv')
dev.off()

## sample
group.name <- 'sample'

pl <- list()
for (i in 1:length(sources)){
    source.curr <- sources[i]
    title <- titles[i]
    pal.curr <- pal[[group.name]]
    p.curr <- plot.clusters.highlight.one(
        SRO = sro.i,
        idx1 = colnames(sro.i), idx2 = sro.i$orig.ident == source.curr,
        groups = get.named.vector.sro(sro.i, group.name),
        clusters.col = source.curr,
        labels = F, label.size = 5,
        col = pal.curr, pref.C = F) + guides(colour = guide_legend(title = NULL,
                                                                   label.theme = element_text(size = 14), 
                                                                   override.aes = list(size=4), ncol = 1)) +
        theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.line = element_blank(),
              plot.title = element_text(size = 21, hjust = 0.5, vjust = -2.5), legend.text.align = 0,
              legend.position = 'right', legend.box.margin = margin(l = -30, t = 0, b = 0, r = 0)) + ggtitle(title)
    
    pl[[i]] <- p.curr
}

ncols <- 3
nrows <- ceiling(length(pl)/ncols)

raster_pdf(file = paste0(pref.p, "UMAP-sample.pdf"), width = 8.5*ncols, height = 5*nrows, res = 300)
plot_grid(plotlist = pl, ncol = ncols, align = 'hv')
dev.off()

## data source
pl <- list()
pl[[1]] <- plot.clusters(sro.i, groups = sro.i$orig.ident, clusters.col = 'source', pref.C = F,
                    col = pal$source, point.size = 1, point.alpha = 0.5, labels = F, label.size = 7,
                    label.pad = 1) + guides(colour = guide_legend(override.aes = list(size=5))) +
    theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.line = element_blank())

for (i in 1:length(sources)){
    source.curr <- sources[i]
    title <- titles[i]
    p.curr <- plot.clusters.highlight.one(
        SRO = sro.i,
        idx1 = colnames(sro.i), idx2 = sro.i$orig.ident == source.curr,
        groups = get.named.vector.sro(sro.i, "orig.ident"),
        clusters.col = source.curr,
        labels = F, label.size = 5,
        col = pal$source, pref.C = F) + guides(colour = guide_legend(title = NULL,
                                                                   label.theme = element_text(size = 14), 
                                                                   override.aes = list(size=4), ncol = 1)) +
        theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.line = element_blank(),
              plot.title = element_text(size = 21, hjust = 0.5, vjust = -2.5), legend.text.align = 0,
              legend.position = 'right', legend.box.margin = margin(l = -30, t = 0, b = 0, r = 0)) + ggtitle(title)
    
    pl[[i+1]] <- p.curr
}

ncols <- 3
nrows <- ceiling(length(pl)/ncols)
raster_pdf(file = paste0(pref.p, "UMAP-source.pdf"), width = 9*ncols, height = 6*nrows, res = 300)
plot_grid(plotlist = pl, ncol = ncols, align = 'hv')
dev.off()
