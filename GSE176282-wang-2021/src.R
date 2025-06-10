setwd("~/0-workspace/CCR7_DC/GSE176282-wang-2021/")
Sys.setenv(RETICULATE_PYTHON = "/data1/lesliec/tyler/utils/miniforge3/envs/multiome/bin/python")

library(Seurat)
library(Rphenograph)
library(Rmagic)
library(Matrix)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

set.seed(1)
options(future.globals.maxSize = Inf)

# pal <- list(
#   R.clusters = c('R1' = '#077315', 'R2' = '#d1c50f', 'R3' = '#2e2bed',
#                  'R4' = '#e60e0e', 'R5' = 'pink'
#   ),
#   Cluster.annot = c(
#     'TC I' = "#2e2bed", 'TC II, III, IV' = "#e60e0e", 'Ki67 TC' = "#d1c50f", 'LTi' = "#077315", 'B cell' = "#f081e6",
#     'Other' = 'darkgray', 'low QC' = 'lightgray')
# )
# pal$Clusters <- pal$Cluster.annot[c('TC I', 'Other', 'Other', 'TC II, III, IV', 'Other', 'Other', 'Other', 'Other', 'Other', 'B cell', 'Ki67 TC',
#                                     'Other', 'Other', 'low QC', 'Other', 'Other', 'Other', 'Other', 'Other', 'LTi', 'Other')]
# cl.col = c(
#   "#1ee3c5", "#7c2da6",
#   "#de9309", "#489de8", "#5fed0e",
#   "#bd537a", "#2f026e", "#1f758c", "#9c6fa8", "#6b6580",
#   "#9e7d3f", "#49709c", "#ccf3ff", "#5c1a1a", "#d1d18a",
#   "#c9a6a1", "#827c68", "#b54800", "#79a695"
# )
# cc <- 1
# for (i in 1:length(pal$Clusters)){
#   if (pal$Clusters[i] == 'darkgray'){
#     pal$Clusters[i] <- cl.col[cc]
#     cc <- cc+1
#   }
# }
# names(pal$Clusters) <- 1:length(pal$Clusters)
# saveRDS(pal, file = "plots/palette.rds")
pal <- readRDS("plots/palette.rds")

plot.QC.violin <- function(sro, ident, feature, yintercept, br, Log10 = F, pal = NULL){
  if(feature == "nCount_RNA"){ name <- "number of transcripts"
  } else if(feature == "nFeature_RNA"){ name <- "number of detected genes"
  } else if(feature == "percent.MT") name <- "percentage of mitochondrial transcripts"
  ft <- sro@meta.data[, feature]
  if(Log10){
    ft <- log10(ft)
    if(!is.null(yintercept)) yintercept <- log10(yintercept)
  }
  gp <- ggplot(data.frame(id = sro@active.ident, ft = ft),
               aes(x = id, y = ft, color = id, fill = id)) +
    geom_violin(show.legend = F) +
    geom_hline(yintercept = yintercept, linetype = 2, color = "violetred4") +
    scale_y_continuous(breaks = seq(0, max(ft), by = br)) +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) +
    labs(x = ident, y = "", title = ifelse(Log10, paste0("log10(", name, ")"), name))
  if(!is.null(pal)) gp <- gp + scale_color_manual(values = pal) + scale_fill_manual(values = pal)
  return(gp)
}

plot.all.QC <- function(sro, ident, thr = NULL, col = NULL){
  list(
    plot.QC.violin(sro, ident = ident, feature = "nCount_RNA", Log10 = T, yintercept = c(thr$nc.min, thr$nc.max), br = 0.2, pal = col),
    plot.QC.violin(sro, ident = ident, feature = "nFeature_RNA", yintercept = c(thr$nf.min, thr$nf.max), br = 500, pal = col),
    plot.QC.violin(sro, ident = ident, feature = "percent.MT", yintercept = thr$mp.max, br = 2, pal = col)
  )
}

run.PhenoGraph <- function(sro, npcs, k){
  sro@misc[["PhenoGraph"]] <- Rphenograph(sro@reductions$pca@cell.embeddings[, 1:npcs], k = k)
  adj <- as_adjacency_matrix(sro@misc$PhenoGraph[[1]], sparse = F)
  rownames(adj) <- colnames(sro); colnames(adj) <- colnames(sro)
  sro@graphs[["PhenoGraph"]] <- as.Graph(adj)
  sro$Clusters <- sro@misc$PhenoGraph[[2]]$membership
  sro$Clusters <- factor(sro$Clusters, levels = min(sro$Clusters):max(sro$Clusters))
  return(sro)
}

plot.groups <- function(
    sro = NULL, idx = NULL, nudge_y = 0.5,
    vis = sro@reductions$umap@cell.embeddings,
    clusters = sro$Clusters, clusters.label = clusters,
    cl.name = "Clusters", col = pal[[cl.name]],
    pref.C = T, labels = T, line = F, point.size = 1
){
  gp <- ggplot() + labs(color = cl.name) + theme_classic() +
    theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
  if(!line) gp <- gp + theme(axis.line = element_blank())
  if(is.null(idx)) idx <- rep(T, nrow(vis))
  gp <- gp + geom_point(aes(x = vis[idx, 1], y = vis[idx, 2], color = clusters[idx]),
                        size = point.size, alpha = 0.7)
  if(!is.null(col)) gp <- gp + scale_color_manual(values = col)
  if(labels){
    C <- levels(clusters)
    CL <- paste0(ifelse(pref.C, "C", ""), levels(clusters.label))
    vis.cl <- data.frame(matrix(NA, nrow = length(C), ncol = 2, dimnames = list(C, NULL)))
    for(c in C){
      vis.cl[c, 1] <- median(vis[idx & clusters == c, 1], na.rm = T)
      vis.cl[c, 2] <- median(vis[idx & clusters == c, 2], na.rm = T)
    }
    gp + geom_label_repel(aes(x = vis.cl[, 1], vis.cl[, 2], color = C, label = CL),
                          size = 5, show.legend = F, nudge_y = nudge_y) + guides(colour = guide_legend(override.aes = list(size=5)))
  } else gp + guides(colour = guide_legend(override.aes = list(size=5)))
}

plot.continuous.value <- function(
    sro = NULL, idx = NULL,
    vis = sro@reductions$umap@cell.embeddings,
    val, val.name,
    scale.color=scale_color_distiller(palette = "Spectral"), point.size = 1
){
  ggplot() +
    geom_point(
      mapping = aes(x = vis[idx, 1], y = vis[idx, 2], color = val[idx]),
      size = point.size, alpha = 0.8
    ) +
    scale.color + theme_classic() + labs(color = val.name)+
    theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
}

break.down.bar.plot <- function(md, group1, group2, ...){
  df <- md[, c(group1, group2)] %>% table() %>%
    as.data.frame() %>% subset(Freq > 0)
  colnames(df)[1:2] <- c("group1", "group2")
  sum.count <- table(md[, group1])
  df$Perc <- round(df$Freq / as.numeric(sum.count[df$group1]), 3) * 100
  df <- df[order(df$Perc, decreasing = F), ]
  ggplot(df, aes(x = group1, y = Perc)) +
    geom_bar(aes(fill = group2), stat = "identity", position = "stack") +
    geom_label(aes(label = paste0(Perc, "%"), color = group2), show.legend = F,
               stat = "identity", position = "stack", ...) +
    scale_color_manual(values = pal[[group2]]) +
    scale_fill_manual(values = pal[[group2]]) +
    scale_y_continuous(breaks = seq(0, 100, by = 10)) +
    labs(x = group1, y = "number of cells", fill = group2) + theme_bw() +
    theme(panel.grid.major = element_blank())
}

# create sro ####
gene.mtx <- Read10X_h5('GSM5362108_RNA_GFP_filtered_feature_bc_matrix.h5')
rownames(gene.mtx) <- toupper(rownames(gene.mtx))
sro <- CreateSeuratObject(counts = gene.mtx, project = "wang")
sro$Sample <- sro$orig.ident
Idents(sro) <- sro$Sample
# sro <- sro[, rownames(md)]
sro <- PercentageFeatureSet(sro, pattern = "^MT-", col.name = "percent.MT")

thr <- data.frame(nc.min = 1000, nc.max = 30000, nf.min = 600, nf.max = 5000, 
                  mp.max = 10)
cell.discard <- sro$nCount_RNA < thr$nc.min | sro$nCount_RNA > thr$nc.max |
  sro$nFeature_RNA < thr$nf.min | sro$nFeature_RNA > thr$nf.max |
  sro$percent.MT > thr$mp.max

pdf("plots/QC/sample-QC.pdf", width = 10, height = 10)
plot.all.QC(sro, ident = "Sample", thr = thr)
dev.off()

sro <- sro[, !cell.discard]
# 10,198 out of 11,697 cells in our analysis

cell.count <- rowSums(sro@assays$RNA@counts > 0)
sro <- NormalizeData(sro[cell.count > 1, ]) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 5000)
sro <- ScaleData(sro, features = rownames(sro)) %>%
  RunPCA(features = rownames(sro), npcs = 50) %>%
  run.PhenoGraph(npcs = 30, k = 30) %>%
  RunUMAP(dims = 1:30, n.neighbors = 30, metric = "cosine", min.dist = 0.4, spread = 1)
Idents(sro) <- sro$Clusters

pdf("plots/QC/cluster-QC.pdf", width = 10, height = 10)
plot.all.QC(sro, ident = "Clusters")
dev.off()

s.genes <- intersect(toupper(cc.genes.updated.2019$s.genes), rownames(sro))
g2m.genes <- intersect(toupper(cc.genes.updated.2019$g2m.genes), rownames(sro))
sro <- CellCycleScoring(sro, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
sro$Phase <- gsub("G2M", "G2/M", sro$Phase)

# add annotations ####
sro <- readRDS("results/SRO.rds")

annot <- openxlsx::read.xlsx(
  xlsxFile = "wang-annotation.xlsx",
  rowNames = F, colNames = T, check.names = F
)
annot$Annotation[is.na(annot$Annotation)] <- "Other"
sro$Cluster.annot <- annot[sro$Clusters, "Annotation"]
sro$Cluster.annot <- factor(sro$Cluster.annot, levels = c("TC I", "TC II, III, IV", "Ki67 TC", 
                                                          "LTi", "B cell", "low QC", "Other"))

# save files ####
write.csv(sro@meta.data, file = "results/meta-data.csv")
write.csv(sro@reductions$pca@cell.embeddings, file = "results/PCA.csv")
write.csv(sro@reductions$umap@cell.embeddings, file = "results/UMAP.csv")
# writeMM(sro@assays$RNA@data, file = "results/unimputed-expr.mtx")
write.csv(sro@assays$RNA@data, file = "results/unimputed-expr.csv")
saveRDS(sro, file = "results/SRO.rds")

sro.imp <- magic(sro)
# saveRDS(sro.imp, file = "results/sro.imp.rds")
# saveRDS(sro.imp@assays$MAGIC_RNA@data, file = "results/imputed-expr.rds")
write.csv(sro.imp@assays$MAGIC_RNA@data, file = "results/imputed-expr.csv")

sro.imp <- ScaleData(sro.imp, features = rownames(sro.imp), assay = "MAGIC_RNA")
saveRDS(sro.imp@assays$MAGIC_RNA@scale.data, file = "results/scale-data.rds")

# ############################################################################ #
# Plot clusters ####
# ############################################################################ #
sro <- readRDS("results/SRO.rds")

pdf("plots/QC/cluster-QC-UMAP.pdf", width = 15, height = 12)
plot.groups(sro)
plot.continuous.value(sro, idx = colnames(sro),
                      val = sro$nCount_RNA, val.name='nCount_RNA', point.size=1)
plot.continuous.value(sro, idx = colnames(sro),
                      val = sro$nFeature_RNA, val.name='nFeature_RNA', point.size=1)
plot.continuous.value(sro, idx = colnames(sro),
                      val = sro$percent.MT, val.name='percent.MT', point.size=1)
dev.off()

pdf("plots/clusters.pdf", width = 15, height = 12)
plot.groups(sro)
dev.off()

# Wang
# R1  C20 LTi
# R2  C11 Ki
# R3  C1 TC I
# R4  C4 mixed TC2/3/4
wang.cl.to.newcl <- c("20"="R1",
                      "11"="R2",
                      "1"="R3",
                      "4"="R4")

sro$R.clusters <- wang.cl.to.newcl[as.character(sro$Clusters)]
sro$R.clusters <- factor(sro$R.clusters, levels = c("R1", "R2", "R3", "R4"))
pdf("plots/R-clusters.pdf", width = 15, height = 12)
plot.groups(sro, clusters = sro$R.clusters, cl.name = 'Cluster', pref.C = F, label = T, col = pal[['R.clusters']])
dev.off()

idx <- sro$Clusters %in% c(1, 4, 11, 20)
pdf("plots/R-clusters-RORgt+cells.pdf", width = 15, height = 12)
plot.groups(sro, clusters = sro$R.clusters, idx = idx, cl.name = 'Cluster', pref.C = F, label = T, col = pal[['R.clusters']])
dev.off()

pdf("plots/annotations.pdf", width = 15, height = 12)
plot.groups(sro, clusters = sro$Cluster.annot, cl.name = 'Cluster.annot', pref.C = F, col = pal[['Cluster.annot']])
dev.off()

# RORgt+ cells only 
idx <- sro$Clusters %in% c(1, 4, 11, 20)
annot.to.include <- as.character(unique(sro@meta.data[sro@meta.data$Clusters %in% c(1,4,11,20),'Cluster.annot']))
ggsave(
  filename = "plots/annotations-RORgt+cells.pdf", width = 10, height = 14,
  plot = plot.groups(
    sro = sro, pref.C = F, labels = T, idx = idx, 
    clusters = sro$Cluster.annot, cl.name = "Annotation", col = pal$Cluster.annot[annot.to.include]
  )
)

# ############################################################################ #
# Redo UMAP and plot clusters ####
# ############################################################################ #
sro.subset <- subset(sro, subset = Clusters %in% c(1, 4, 11, 20)) %>%
  run.PhenoGraph(npcs = 30, k = 30) %>%
  RunUMAP(dims = 1:30, n.neighbors = 30, metric = "cosine", min.dist = 0.4, spread = 0.7)

write.csv(sro.subset@reductions$umap@cell.embeddings, file = "results/UMAP-subset.csv")

pdf("plots/R-clusters-RORgt+cells.pdf", width = 15, height = 12)
plot.groups(sro.subset, clusters = sro.subset$R.clusters, cl.name = 'Cluster', pref.C = F, label = T, col = pal[['R.clusters']])
dev.off()

# RORgt+ cells only 
annot.to.include <- as.character(unique(sro@meta.data[sro@meta.data$Clusters %in% c(1,4,11,20),'Cluster.annot']))
ggsave(
  filename = "plots/annotations-RORgt+cells.pdf", width = 15, height = 12,
  plot = plot.groups(
    sro = sro.subset, pref.C = F, labels = T, 
    clusters = sro.subset$Cluster.annot, cl.name = "Annotation", col = pal$Cluster.annot[annot.to.include]
  )
)
# ############################################################################ #
# for module scores ####
sro <- readRDS("results/SRO.rds")
## TC I ~ IV ####
TC.module <- c('Mki67', 'Gal', 'Nrxn1', 'Aire', 'Kif21a', 'Pigr', 'Col17a1', 
               'Hk2', 'LTb', 'Dnase1l3', 'Ahcyl2', 'Nlrc5', 'Itgb8', 'Ccl22', 'Ccl5', 'Il2ra')
TC_signatures <- read.csv('../MLN_RORgt_MHCII_multiome/Seurat/results/markers/C2-5_top130.csv')
TC_I <- toupper(TC_signatures[TC_signatures$celltype == 'TC I',]$gene_name)
TC_II <- toupper(TC_signatures[TC_signatures$celltype == 'TC II',]$gene_name)
TC_III <- toupper(TC_signatures[TC_signatures$celltype == 'TC III',]$gene_name)
TC_IV <- toupper(TC_signatures[TC_signatures$celltype == 'TC IV',]$gene_name)

TC.modules <- list(TC_I, TC_II, 
                   TC_III, TC_IV)
sro <- AddModuleScore(
  object = sro,
  features = TC.modules,
  name = 'TC'
)

## ILC3, JC1, JC2, JC3 ####
JC1 <- toupper(c('Slc7a10', 'Dcaf12l2', 'Olig1', 'Gal', 'Atp1b1', 'Dsg1b', 'Ttyh1', 'Tbx5', 'Cnr1', 'Ank', 'Fam81a', 'B3galt1', 'Ube2e2', 'Syt1', 'Zfand6'))
JC2 <- toupper(c('Egfl6', 'Tnni1', '1110008L16Rik', 'Cep112', 'Asic1', 'Ly9', 'Fabp1', 'Col17a1', 'Pgam2', 
                 'Poc1a', 'Clic3', 'Prdm16', 'Ppp2r2c', 'Gstt2'))
JC3 <- toupper(c('Gm26917', 'Ptbp2', 'Zc3h7a', 'Lcor', 'Nfat5', 'Smg1', 'Cep350', 'Mdm4', 'Chuk', 
                 'Mapk8ip3', 'Prpf39', 'Eml5', 'Phip', 'Rnf111', 'Trpm7'))
ILC3 <- toupper(c('Cxcr6', 'Clnk', 'Fam184b', 'Klrb1b', 'Klrb1f', 'Chad', 'Apol7e', 'Ncr1', 
                  'Il22', 'Arg1', 'Il2rb', 'Dgat1', 'Il18rap', 'Gzmb', 'Ccdc184'))
JC.modules <- list(JC1, JC2, 
                   JC3)
sro <- AddModuleScore(
  object = sro,
  features = JC.modules,
  name = 'JC'
)

sro <- AddModuleScore(
  object = sro,
  features = list(ILC3),
  name = 'ILC3'
)


# ############################################################################ #
# Re-cluster ####
# ############################################################################ #
sro <- readRDS("results/SRO.rds")
reslist <- seq(0.2, 2, 0.2)
sro <- sro %>%
    FindNeighbors(dims = 1:30, k.param = 30) %>%
    FindClusters(resolution = reslist)

pref <- "results/"; dir.create(pref)
write.csv(sro@meta.data, file = paste0(pref, "meta-data.csv"), quote = F)
saveRDS(sro, file = paste0(pref, "SRO.rds"))

## plot ####
pal$Clusters2 = c(
    "#00bf00", "#489de8", "#d40663", "#f8c72f", "#077315",
    "#785cd4", "#e67109", "#0eefff", "#f081e6", "#260691",
    "#49709c", "#9e7d3f", "#bd537a", "#4e225c", "#f202ed",
    "#fec55f", "#062e0b", "#9c6fa8", "#078d94", "#5c1a1a",
    "#827c68", "#aebeff", "#9c2903", "#ffc5af", "#4f5715",
    "#0249f0", "#f43525", "#0077ff", "#7f227e", "#dfddff",
    "#7e85d7", "#fff64f", "#5fed0e", "#543018", "#f31220"
)

pref.p <- 'plots/'
pdf(paste0(pref.p, "UMAP-Clusters.pdf"), width = 15, height = 12)
for (res in reslist){
    colname <- paste0("RNA_snn_res.", res)
    print(
        plot.groups(sro, clusters = sro@meta.data[[colname]],
                    cl.name = colname, col = pal$Clusters2,
                    point.size = 0.5,
                    pref.C = T)
    )
}
dev.off()
