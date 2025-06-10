setwd("~/0-workspace/CCR7_DC/GSE184175-lyu-2022/")

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
# pal ####
# pal <- list(
#   Clusters = c(
#     "#e60e0e", "#077315", "#1ee3c5", "#f081e6", "#7c2da6",
#     "#de9309", "#d1c50f", "#489de8", "#2e2bed", "#5fed0e",
#     "#bd537a", "#2f026e", "#1f758c", "#9c6fa8", "#6b6580",
#     "#9e7d3f", "#49709c", "#ccf3ff", "#5c1a1a", "#d1d18a",
#     "#c9a6a1", "#827c68", "#b54800", "#79a695"
#   ),
#   R.clusters = c('R1' = '#077315', 'R2' = '#2e2bed', 'R3' = '#e60e0e',
#                  'R4' = '#d1c50f', 'R5' = 'pink'
#   ),
#   Cluster.annot.from.paper = c(
#     'eTACs I' = "#2e2bed", 'eTACs II' ="#e60e0e",
#     'Treg I' = '#b2df8a', 'Treg II' = '#1f758c', 'Treg III' = '#5fed0e', 
#     'Th17' = '#1ee3c5', 'ILC2s' = '#6a3d9a',
#     'LTi-like ILC3s' = '#077315', 'T-bet+ ILC3s I' = '#ccf3ff', 'T-bet+ ILC3s II' = '#489de8', 
#     'gamma-delta T cells I' = '#ff7f00', 'gamma-delta T cells II' = '#fdbf6f', 
#     'B cells'="#f081e6",
#     'Doublets' = 'gray'
#   ),
#   annotations = c(
#     "#e60e0e", "#077315", "#1ee3c5", "#f081e6", "#7c2da6",
#     "#de9309", "#d1c50f", "#489de8", "#2e2bed", "#5fed0e",
#     "#bd537a", "#2f026e", "#1f758c", "#9c6fa8", "#6b6580",
#     "#9e7d3f", "#49709c", "#ccf3ff", "#5c1a1a", "#d1d18a",
#     "#c9a6a1", "#827c68", "#b54800", "#79a695"
#   )
#   # Cluster.annot.from.paper = c(
#   #   'eTAC I' = "#e60e0e", 'eTAC II' ="#c994c7", 
#   #   'LTi cell' ="#d1c50f", 'pTreg' ="#ccf3ff", 'gd T cell' ="#2e2bed",
#   #   'Th17'="#489de8", 'NCR1+ ILC3'="#077315", 'ILC2'="#5fed0e", 'Proliferating'="#de9309", 
#   #   'B cell'="#1ee3c5", 'Other'='darkgray'
#   # )
#   # Cluster.annot = c(
#   #   'TC I' = "#e60e0e", 'TC II/III/IV' ="#f081e6", 'Transitional ILC3-TC' ="#7c2da6", 
#   #   'LTi cell' ="#d1c50f", 'pTreg' ="#ccf3ff", 'gd T cell' ="#2e2bed",
#   #   'Th17'="#489de8", 'NCR1+ ILC3'="#077315", 'ILC2'="#5fed0e", 'Proliferating'="#de9309", 
#   #   'B cell'="#1ee3c5",
#   #   'doublet'='lightgray', 'Other'='darkgray'
#   # ),
# )
# pal$Clusters <- pal$Cluster.annot.from.paper[c('Treg I', 'Th17', 'LTi-like ILC3s', 'T-bet+ ILC3s I', 'gamma-delta T cells I', 'Treg II',
#                                'eTACs I', 'B cells', 'Doublets', 'gamma-delta T cells II', 'ILC2s', 'eTACs II', 'T-bet+ ILC3s II', "Treg III")]
# names(pal$Clusters) <- 0:(length(pal$Clusters)-1)
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

add.custom.umap <- function(sro, umap, slot.name, key, 
                            global = T, assay = "RNA"){
  sro[[slot.name]] <- CreateDimReducObject(embeddings = as.matrix(umap), 
                                           key = key, global = global, 
                                           assay = assay)
  return(sro)
}
# create sro ####
gene.mtx <- Read10X_h5('GSM5579608_filtered_feature_bc_matrix.h5')
rownames(gene.mtx) <- toupper(rownames(gene.mtx))
head(gene.mtx)
md <- read.csv("GSM5579608_cell_metadata.tsv.gz", sep='\t') # original md fromt the paper
rownames(md) <- md$barcode
colnames(md)[4] <- 'percent.MT'
gene.mtx <- gene.mtx[, md$barcode]
dim(gene.mtx)
# 13,547 cells

sro <- CreateSeuratObject(counts = gene.mtx, project = "lyu", meta.data = md)
sro$Sample <- sro$orig.ident
Idents(sro) <- sro$Sample
# sro <- PercentageFeatureSet(sro, pattern = "^MT-", col.name = "percent.MT")
thr <- data.frame(nf.max = 5000, nf.min = 600, mp.max = 10)
cell.discard <- sro$nFeature_RNA > thr$nf.max |
  sro$nFeature_RNA < thr$nf.min |
  sro$percent.MT > thr$mp.max
# cell.discard = 0
pdf("plots/QC/sample-QC.pdf", width = 10, height = 10)
plot.all.QC(sro, ident = "Sample", thr = thr)
dev.off()

sro <- sro[, !cell.discard]

cell.count <- rowSums(sro@assays$RNA@counts > 0)
sro <- NormalizeData(sro[cell.count > 1, ]) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 5000)
sro <- ScaleData(sro, features = rownames(sro))
sro$Clusters <- factor(sro$seurat_clusters)
Idents(sro) <- factor(sro$seurat_clusters)
pdf("plots/QC/cluster-QC.pdf", width = 10, height = 10)
plot.all.QC(sro, ident = "seurat_clusters")
dev.off()

# add paper's umap
umap.df <- data.frame(UMAP_1 = md[rownames(sro@meta.data),]$UMAP_1,
                      UMAP_2 = md[rownames(sro@meta.data),]$UMAP_2,
                      row.names = rownames(sro@meta.data))
sro <- add.custom.umap(sro, 
                       umap.df, 
                       'umap', 'umap')

# s.genes <- intersect(toupper(cc.genes.updated.2019$s.genes), rownames(sro))
# g2m.genes <- intersect(toupper(cc.genes.updated.2019$g2m.genes), rownames(sro))
# sro <- CellCycleScoring(sro, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
# sro$Phase <- gsub("G2M", "G2/M", sro$Phase)

# ############################################################################ #
# annotations from paper ####
# ############################################################################ #
length(intersect(names(pal$Cluster.annot.from.paper), cl.to.paper.annotation)) == length(pal$Cluster.annot.from.paper)
cl.to.paper.annotation <- c('0' = 'Treg I',
                            '1' = 'Th17', 
                            '2' = 'LTi-like ILC3s', 
                            '3' = 'T-bet+ ILC3s I', 
                            '4' = 'gamma-delta T cells I', 
                            '5' = 'Treg II', 
                            '6' = 'eTACs I',
                            '7' = 'B cells',
                            '8' = 'Doublets',
                            '9' = 'gamma-delta T cells II',
                            '10' = 'ILC2s',
                            '11' = 'eTACs II',
                            '12' = 'T-bet+ ILC3s II',
                            '13' = 'Treg III')

sro$annotations <- cl.to.paper.annotation[as.character(sro$seurat_clusters)]
sro$annotations <- factor(sro$annotations, levels = names(pal$Cluster.annot.from.paper))

# save #####
saveRDS(sro, file = "results/SRO.rds")

write.csv(sro@meta.data, file = "results/meta-data.csv")
# write.csv(sro@reductions$pca@cell.embeddings, file = "results/PCA.csv")
write.csv(sro@reductions$umap@cell.embeddings, file = "results/UMAP.csv")
saveRDS(as.matrix(sro@assays$RNA@data), file = "results/unimputed-expr.rds")
# writeMM(sro@assays$RNA@data, file = "results/unimputed-expr.mtx")
# write.csv(sro@assays$RNA@data, file = "results/unimputed-expr.csv")

sro.imp <- magic(sro)
# saveRDS(sro.imp, file = "results/sro.imp.rds")
saveRDS(sro.imp@assays$MAGIC_RNA@data, file = "results/imputed-expr.rds")
# write.csv(sro.imp@assays$MAGIC_RNA@data, file = "results/imputed-expr.csv")

# sro.imp <- ScaleData(sro.imp, features = rownames(sro.imp), assay = "MAGIC_RNA")
# saveRDS(sro.imp@assays$MAGIC_RNA@scale.data, file = "results/scale-data.rds")

# ############################################################################ #
# Plot clusters ####
# ############################################################################ #
sro <- readRDS("results/SRO.rds")

pdf("plots/QC/cluster-QC-UMAP.pdf", width = 15, height = 12)
plot.groups(
  sro = sro, pref.C = T, labels = T, vis = sro@reductions$umap@cell.embeddings,
  clusters = sro$seurat_clusters, cl.name = "seurat_clusters", col = pal$Clusters
)
plot.continuous.value(sro, idx = rownames(sro@meta.data),
                      val = sro$nCount_RNA, val.name='nCount_RNA', point.size=1)
plot.continuous.value(sro, idx = rownames(sro@meta.data),
                      val = sro$nFeature_RNA, val.name='nFeature_RNA', point.size=1)
plot.continuous.value(sro, idx = rownames(sro@meta.data),
                      val = sro$percent.MT, val.name='percent.MT', point.size=1)
dev.off()

ggsave(
  filename = "plots/clusters.pdf", width = 15, height = 12,
  plot = plot.groups(
    sro = sro, pref.C = T, labels = T, vis = sro@reductions$umap@cell.embeddings,
    clusters = sro$seurat_clusters, cl.name = "seurat_clusters", col = pal$Clusters
  )
)

###
# Lyu
# R1 C2 ILC3
# R2 C6 eTAC I
# R3 C11 eTAC II
lyu.cl.to.newcl <- c("2"="R1", 
                     "6"="R2",
                     "11"="R3"
)

sro$R.clusters <- lyu.cl.to.newcl[as.character(sro$Clusters)]
sro$R.clusters <- factor(sro$R.clusters, levels = c("R1", "R2", "R3"))
pdf("plots/R-clusters.pdf", width = 15, height = 12)
plot.groups(sro, clusters = sro$R.clusters, cl.name = 'Cluster', pref.C = F, label = T, col = pal[['R.clusters']])
dev.off()

idx <- sro$Clusters %in% c(2, 6, 11)
pdf("plots/R-clusters-RORgt+cells.pdf", width = 15, height = 12)
plot.groups(sro, clusters = sro$R.clusters, idx = idx, cl.name = 'Cluster', pref.C = F, label = T, col = pal[['R.clusters']])
dev.off()

ggsave(
  filename = "plots/annotations.pdf", width = 15, height = 12,
  plot = plot.groups(
    sro = sro, pref.C = F, labels = T,
    clusters = sro$annotations, cl.name = "Annotation", col = pal$Cluster.annot.from.paper
  )
)

## RORgt+ cells only ####
idx <- sro$Clusters %in% c(2, 6, 11)
annot.to.include <- as.character(unique(sro@meta.data[sro@meta.data$Clusters %in% c(2, 6, 11),'annotations']))
ggsave(
  filename = "plots/annotations-RORgt+cells.pdf", width = 10, height = 14,
  plot = plot.groups(
    sro = sro, pref.C = F, labels = T, idx = idx, 
    clusters = sro$annotations, cl.name = "Annotation", col = pal$Cluster.annot.from.paper[annot.to.include]
  )
)

# ############################################################################ #
# Redo UMAP and plot clusters ####
# ############################################################################ #
sro.subset <- subset(sro, subset = Clusters %in% c(2, 6, 11)) %>%
  ScaleData(features = rownames(sro)) %>%
  RunPCA(features = rownames(sro), npcs = 50) %>%
  run.PhenoGraph(npcs = 30, k = 30) %>%
  RunUMAP(dims = 1:30, n.neighbors = 30, metric = "cosine", min.dist = 0.4, spread = 0.7)

write.csv(sro.subset@reductions$umap@cell.embeddings, file = "results/UMAP-subset.csv")

pdf("plots/R-clusters-RORgt+cells.pdf", width = 15, height = 12)
plot.groups(sro.subset, clusters = sro.subset$R.clusters, cl.name = 'Cluster', pref.C = F, label = T, col = pal[['R.clusters']])
dev.off()

# RORgt+ cells only 
annot.to.include <- as.character(unique(sro@meta.data[sro@meta.data$Clusters %in% c(2, 6, 11), 'annotations']))
ggsave(
  filename = "plots/annotations-RORgt+cells.pdf", width = 15, height = 12,
  plot = plot.groups(
    sro = sro.subset, pref.C = F, labels = T, 
    clusters = sro.subset$annotations, cl.name = "Annotation", col = pal$Cluster.annot.from.paper[annot.to.include]
  )
)

# ############################################################################ #
# Dot plot ####
# ############################################################################ #
sro <- readRDS("results/SRO.rds")
sro.subset <- subset(sro, subset = Clusters %in% c(2, 6, 11))

genes <- c("CXCR6", "RORA", "AIRE", "GAL", 'COL17A1', 'H2-K1', 
           'DNASE1L3', 'NLRC5', 
           'ITGB8', 'CCL22')
Idents(sro) <- sro$annotations
Idents(sro.subset) <- sro.subset$annotations

ggsave(
  filename = "plots/dot-all-cells.pdf", width = 15, height = 12,
  plot = DotPlot(sro, assay = "RNA", group.by = "annotations",
                 features = genes, cols = c("blue", "red"), dot.scale = 7) +
    scale_x_discrete(breaks = genes, labels = genes) +
    theme(axis.text.x.bottom = element_text(angle = 90, vjust = 0, size = 11)) + labs(x = "", y = "Annotation")
)


ggsave(
  filename = "plots/dot-RORgt+cells.pdf", width = 15, height = 12,
  plot = DotPlot(sro.subset, assay = "RNA", group.by = "annotations",
                 features = genes, cols = c("blue", "red"), dot.scale = 7) +
    scale_x_discrete(breaks = genes, labels = genes) +
    theme(axis.text.x.bottom = element_text(angle = 90, vjust = 0, size = 11)) + labs(x = "", y = "Annotation")
)

