# Kedmi
suppressPackageStartupMessages({
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
})
set.seed(1)
options(future.globals.maxSize = Inf)
# pal <- list(
#   Cluster.annot = c(
#     'JC' = "#e60e0e", 'ILC3' = "#077315", 'B cell' = "#f081e6", 'T cell' = "#ff7f00", 'cDC1' = "#7c2da6",
#     'Other' = 'gray'),
#   R.clusters = c('R1' = '#077315', 'R2' = '#d1c50f', 'R3' = '#2e2bed',
#                  'R4' = '#e60e0e', 'R5' = 'pink'
#   )
# )
# pal$Clusters <- pal$Cluster.annot[c('JC', 'Other', 'Other', 'Other', 'Other', 'Other', 'Other', 'Other', 'Other', 'JC', 'Other', 'JC', 'T cell', 'Other',
#                                     'cDC1', 'JC', 'B cell', 'ILC3', 'Other')]
# cl.col <- c(
#   "#1ee3c5", "#de9309", "#5fed0e", "#489de8", "#2e2bed",
#   "#bd537a", "#2f026e", "#1f758c", "#9c6fa8", "#6b6580",
#   "#9e7d3f", "#49709c", "#ccf3ff", "#5c1a1a", "#d1d18a",
#   "#c9a6a1", "#827c68", "#b54800", "#79a695"
# )
# cc <- 1
# for (i in 1:length(pal$Clusters)){
#   if (pal$Clusters[i] == 'gray'){
#     pal$Clusters[i] <- cl.col[cc]
#     cc <- cc+1
#   }
# }
# names(pal$Clusters) <- 1:length(pal$Clusters)
# pal$Clusters[1] <- pal$R.clusters['R3']
# pal$Clusters[10] <- pal$R.clusters['R5']
# pal$Clusters[12] <- pal$R.clusters['R4']
# pal$Clusters[16] <- pal$R.clusters['R2']
# pal$Clusters[18] <- pal$R.clusters['R1']
# saveRDS(pal, file = "plots/palette.rds")
pal <- readRDS("plots/palette.rds")
pal$Clusters2 = c(
    "#00bf00", "#489de8", "#d40663", "#f8c72f", "#077315",
    "#785cd4", "#e67109", "#0eefff", "#f081e6", "#260691",
    "#49709c", "#9e7d3f", "#bd537a", "#4e225c", "#f202ed",
    "#fec55f", "#062e0b", "#9c6fa8", "#078d94", "#5c1a1a",
    "#827c68", "#aebeff", "#9c2903", "#ffc5af", "#4f5715",
    "#0249f0", "#f43525", "#0077ff", "#7f227e", "#dfddff",
    "#7e85d7", "#fff64f", "#5fed0e", "#543018", "#f31220"
)

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
gene.mtx <- Read10X_h5('GSM6012929_RNA_GFP_Pos_5p.h5')
rownames(gene.mtx) <- toupper(rownames(gene.mtx))

# Samples were merged and filtered, removing genes with no counts, and 
# retaining cells with 600 to 5,000 genes and 1,000 to 30,000 counts, 
# leaving a total of 35,797 cells and 22,979 genes in the dataset.

sro <- CreateSeuratObject(counts = gene.mtx, project = "kedmi")
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
# 10,554 out of 11,772 cells in our analysis

cell.count <- rowSums(sro@assays$RNA@counts > 0)
sro <- NormalizeData(sro[cell.count > 1, ]) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 5000)
sro <- ScaleData(sro, features = rownames(sro)) %>%
  RunPCA(features = rownames(sro), npcs = 50) %>%
  run.PhenoGraph(npcs = 30, k = 30) %>%
  RunUMAP(dims = 1:30, n.neighbors = 30, metric = "cosine", min.dist = 0.4, spread = 1)
Idents(sro) <- sro$Clusters

Idents(sro) <- sro$Clusters
pdf("plots/QC/cluster-QC.pdf", width = 10, height = 10)
plot.all.QC(sro, ident = "Clusters", col=pal$Clusters)
dev.off()


sro.r <- subset(sro, R.clusters %in% c("R1", "R2", "R3", "R4", "R5"))
Idents(sro.r) <- sro.r$R.clusters
pdf("plots/QC/R.cluster-QC.pdf", width = 10, height = 10)
plot.all.QC(sro.r, ident = "R.clusters", col=pal$R.clusters)
dev.off()


s.genes <- intersect(toupper(cc.genes.updated.2019$s.genes), rownames(sro))
g2m.genes <- intersect(toupper(cc.genes.updated.2019$g2m.genes), rownames(sro))
sro <- CellCycleScoring(sro, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
sro$Phase <- gsub("G2M", "G2/M", sro$Phase)

# add annotations ####
sro <- readRDS("results/SRO.rds")

annot <- openxlsx::read.xlsx(
  xlsxFile = "kedmi-annotation.xlsx",
  rowNames = F, colNames = T, check.names = F
)
annot$X2[is.na(annot$X2)] <- "Other"
sro$Cluster.annot <- annot[sro$Clusters, "X2"]
sro$Cluster.annot <- factor(sro$Cluster.annot, levels = c("JC", "ILC3", "B cell", "T cell", "cDC1", "Other"))

# save files ####
write.csv(sro@meta.data, file = "results/meta-data.csv")
write.csv(sro@reductions$pca@cell.embeddings, file = "results/PCA.csv")
write.csv(sro@reductions$umap@cell.embeddings, file = "results/UMAP.csv")
# writeMM(sro@assays$RNA@data, file = "results/unimputed-expr.mtx")
write.csv(sro@assays$RNA@data, file = "results/unimputed-expr.csv")
saveRDS(sro, file = "results/SRO.rds")

sro.imp <- magic(sro)
# saveRDS(sro.imp, file = "results/sro.imp.rds")
saveRDS(sro.imp@assays$MAGIC_RNA@data, file = "results/imputed-expr.rds")
write.csv(sro.imp@assays$MAGIC_RNA@data, file = "results/imputed-expr.csv")

# sro.imp <- ScaleData(sro.imp, features = rownames(sro.imp), assay = "MAGIC_RNA")
# saveRDS(sro.imp@assays$MAGIC_RNA@scale.data, file = "results/scale-data.rds")

# ############################################################################ #
# Re-clusters ####
# ############################################################################ #
sro <- readRDS("results/SRO.rds")
reslist <- seq(0.2, 2, 0.2)
sro <- sro %>%
    FindNeighbors(dims = 1:30, k.param = 30) %>%
    FindClusters(resolution = reslist)

pref <- "results/"; dir.create(pref)
write.csv(sro@meta.data, file = paste0(pref, "meta-data.csv"), quote = F)
saveRDS(sro, file = paste0(pref, "SRO.rds"))
