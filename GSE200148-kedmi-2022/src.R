setwd("~/0-workspace/CCR7_DC/GSE200148-kedmi-2022/")

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

###
# Kedmi
# R1 C18 ILC3
# R2 C16 Ki
# R3 C1 TC I
# R4 12 mixed TC2/3/4
# R5 C10 lowQC
kedmi.cl.to.newcl <- c("18"="R1",
                       "16"="R2",
                       "1"="R3",
                       "12"="R4",
                       "10"="R5"
)

sro$R.clusters <- kedmi.cl.to.newcl[as.character(sro$Clusters)]
sro$R.clusters <- factor(sro$R.clusters, levels = c("R1", "R2", "R3", "R4", "R5"))
pdf("plots/R-clusters.pdf", width = 15, height = 12)
plot.groups(sro, clusters = sro$R.clusters, cl.name = 'Cluster', pref.C = F, label = T, col = pal[['R.clusters']])
dev.off()

idx <- sro$Clusters %in% c(1, 10, 12, 16, 18)
pdf("plots/R-clusters-RORgt+cells.pdf", width = 15, height = 12)
plot.groups(sro, clusters = sro$R.clusters, idx = idx, cl.name = 'Cluster', pref.C = F, label = T, col = pal[['R.clusters']])
dev.off()

##
pdf("plots/annotations.pdf", width = 15, height = 12)
plot.groups(sro, clusters = sro$Cluster.annot, cl.name = 'Cluster.annot', pref.C = F, col = pal[['Cluster.annot']])
dev.off()

# RORgt+ cells only 
idx <- sro$Clusters %in% c(1, 10, 12, 16, 18)
annot.to.include <- as.character(unique(sro@meta.data[sro@meta.data$Clusters %in% c(1, 10, 12, 16, 18),'Cluster.annot']))
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
sro.subset <- subset(sro, subset = Clusters %in% c(1, 10, 12, 16, 18)) %>%
  run.PhenoGraph(npcs = 30, k = 30) %>%
  RunUMAP(dims = 1:30, n.neighbors = 30, metric = "cosine", min.dist = 0.4, spread = 0.7)

write.csv(sro.subset@reductions$umap@cell.embeddings, file = "results/UMAP-subset.csv")

pdf("plots/R-clusters-RORgt+cells.pdf", width = 15, height = 12)
plot.groups(sro.subset, clusters = sro.subset$R.clusters, cl.name = 'Cluster', pref.C = F, label = T, col = pal[['R.clusters']])
dev.off()

# RORgt+ cells only 
annot.to.include <- as.character(unique(sro@meta.data[sro@meta.data$Clusters %in% c(1, 10, 12, 16, 18), 'Cluster.annot']))
ggsave(
  filename = "plots/annotations-RORgt+cells.pdf", width = 15, height = 12,
  plot = plot.groups(
    sro = sro.subset, pref.C = F, labels = T, 
    clusters = sro.subset$Cluster.annot, cl.name = "Annotation", col = pal$Cluster.annot[annot.to.include]
  )
)

# ############################################################################ #
# MAST DEG ####
# ############################################################################ #
sro <- readRDS("results/SRO.rds")
imp.expr <- readRDS("results/imputed-expr.rds")
sro@assays$RNA@meta.features$gene_name <- rownames(sro@assays$RNA@meta.features)
fn <- "C10-vs-rest"
# MAST.DE(fn, sro, group.by = "Clusters", ident.1 = 10, ident.2 = c(1:9, 11:19))
pdf(width = 20, height = 50, file = paste0("plots/markers/", fn, ".pdf"))
DE.heatmap(
  sro, expr = imp.expr, fn = fn, pairwise = T,
  cells = rownames(subset(sro@meta.data, Clusters %in% c(10, c(1:9, 11:19)))),
  fc.thr = 2.5
)
dev.off()
####
fn <- "R5-vs-R1,2,3,4"
MAST.DE(fn, sro, group.by = "R.clusters", ident.1 = 'R5', ident.2 = c('R1','R2','R3','R4'))
pdf(width = 20, height = 50, file = paste0("plots/markers/", fn, ".pdf"))
DE.heatmap(
  sro, expr = imp.expr, fn = fn, pairwise = T,
  cells = rownames(subset(sro@meta.data, R.clusters %in% c('R1','R2','R3','R4','R5'))),
  fc.thr = 2.5, split = 'R.clusters'
)
dev.off()

pdf(width = 20, height = 30, file = paste0("plots/markers/", fn, "-top50.pdf"))
DE.heatmap(
  sro, expr = imp.expr, fn = fn, pairwise = T,
  cells = rownames(subset(sro@meta.data, R.clusters %in% c('R1','R2','R3','R4','R5'))),
  fc.thr = 2.5, split = 'R.clusters', n=50
)
dev.off()

MAST.DE <- function(fn, sro, ...){
  m <- suppressWarnings(FindMarkers(object = sro, test.use = "MAST", logfc.threshold = 0, ...))
  m$gene_name <- sro@assays$RNA@meta.features[rownames(m), "gene_name"]
  write.csv(m, file = paste0("results/markers/", fn, ".csv"), row.names = T, quote = F)
}

select.markers <- function(fn, pairwise = F, fc.thr = 1.5, apv.thr = 0.01, n = Inf){
  markers <- read.csv(file = paste0("results/markers/", fn, ".csv"), header = T, stringsAsFactors = F) %>%
    subset(!grepl("(^MT-)|(^RPS)|(^RPL)|(^MRPL)|(^MRPS)", toupper(gene_name)))
  if(pairwise){
    markers <- subset(markers, abs(avg_log2FC) > log2(fc.thr) & p_val_adj < apv.thr)
    markers$cluster <- ifelse(markers$avg_log2FC > 0, "up", "down")
    if(!is.infinite(n)) markers <- markers %>% group_by(cluster) %>% top_n(n, abs(avg_log2FC))
    markers <- markers[order(markers$avg_log2FC), ]
  } else{
    markers <- subset(markers, avg_log2FC > log2(fc.thr) & p_val_adj < apv.thr)
    if(any(duplicated(markers$gene_name))) markers <- markers %>% group_by(gene_name) %>% top_n(1, avg_log2FC)
    if(!is.infinite(n)) markers <- markers %>% group_by(cluster) %>% top_n(n, abs(avg_log2FC))
    markers <- markers[order(markers$cluster, markers$avg_log2FC), ]
  }
  return(markers)
}

DE.heatmap <- function(
    sro = NULL, expr = sro@assays$RNA@data, md = sro@meta.data,
    cells = rownames(md), split = "Clusters",
    cluster_columns = T, cluster_column_slices = T, ...
){
  markers <- select.markers(...)
  ca <- columnAnnotation(
    df = sro@meta.data[cells, c("R.clusters", "Cluster.annot")],
    col = list(R.clusters = pal$R.clusters, Cluster.annot = pal$Cluster.annot)
  )
  col.ramp <- colorRamp2(c(-10, -5, seq(-3, 3, 1), 5, 10),
                         c("#173c68", "#173c68", "#1c5785", "#166db0", "#a6d1e3", "#ffffff",
                           "#f09173", "#d46052", "#b32339", "#690822", "#690822"))
  Heatmap(
    matrix = t(scale(t(expr[markers$gene_name, cells]))), name = "scaled\nimputed\nexpression", col = col.ramp,
    show_column_names = F, column_split = md[cells, split], top_annotation = ca,
    cluster_rows = F, row_names_side = "right",
    cluster_columns = cluster_columns, cluster_column_slices = cluster_column_slices,
    row_names_gp = gpar(fontface = "italic", fontsize = 7), use_raster = T
  )
}

# ############################################################################ #
# Violin plot ####
# ############################################################################ #
sro <- readRDS("results/SRO.rds")
expr <- sro@assays$RNA@data
expr <- t(expr)
md <- sro@meta.data
md <- md[!(is.na(md$R.clusters)),]
df <- cbind(md, expr[rownames(md), ])

plot.violin <- function(
    vis, mapping,
    groups, title, legend.title, jitter = F, pal = NULL, y.axis.title = 'Expression level',
    et = element_text(size = 15)
){
  eb <- element_blank()
  gp <- ggplot(mapping = mapping) + theme_classic() +
    geom_violin(aes(fill = .data[[groups]]), vis, show.legend = F) + 
    labs(title=title, y = y.axis.title, x = "") +
    theme(axis.text = et,
          axis.title = et, title = et,
          legend.title = et, legend.text = et) + guides(fill=guide_legend(title=legend.title))
  if (jitter){
    gp <- gp + geom_jitter(aes(fill = .data[[groups]]), vis, shape=16, position=position_jitter(0.2))
  }
  if (!(is.null(pal))){
    gp <- gp + scale_fill_manual(values = pal)
  }
  return(gp)
}

ggsave(
  filename = "plots/violin-Rclusters-Malat1.pdf", width = 10, height = 10,
  plot = plot.violin(df, mapping=aes(x=R.clusters, y=MALAT1), groups='R.clusters', pal=pal$R.clusters,
                     title='MALAT1', legend.title = 'Cluster')
)



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
