# #####
setwd("~/CCR7_DC/GSE200148-kedmi-2022/")

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

set.seed(1)
options(future.globals.maxSize = Inf)

pal <- readRDS("plots/palette.rds")
pal[['Cluster.annot']] <- c(
  'JC' = "#e60e0e", 'ILC3' = "#2e2bed",
  'TC I' ='#f081e6',
  'TC II, III, IV' = '#7c2da6',
  'Ki67 TC' = '#5c1a1a',
  'LTi' = '#2e2bed'
  )
pal$Sample <- c("#e60e0e", "blue")

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
      vis.cl[c, 1] <- median(vis[idx & clusters == c, 1])
      vis.cl[c, 2] <- median(vis[idx & clusters == c, 2])
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

MAST.DE <- function(fn, sro, ...){
  m <- suppressWarnings(FindMarkers(object = sro, test.use = "MAST", logfc.threshold = 0, ...))
  m <- cbind(m, sro@assays$RNA@meta.features[rownames(m), 7:9])
  write.csv(m, file = paste0("results/markers/", fn, ".csv"), row.names = T, quote = F)
  return(m)
}

MAST.DE.multiple <- function(fn, sro, idents = sro@active.ident, ...){
  Idents(sro) <- idents
  m <- FindAllMarkers(object = sro, test.use = "MAST", logfc.threshold = 0, ...)
  m <- cbind(m, sro@assays$RNA@meta.features[m$gene, 7:9])
  write.csv(m, file = paste0("results/markers/", fn, ".csv"), row.names = F, quote = F)
  return(m)
}

select.markers <- function(fn, pairwise = F, fc.thr = 1.5, apv.thr = 0.01, n = Inf){
  markers <- read.csv(file = paste0("results/markers/", fn, ".csv"), header = T, stringsAsFactors = F) %>%
    subset(!grepl("(^MT-)|(^RPL)|(^RPS)|(^MRPL)", toupper(cellranger.gene_name)) &
             !grepl("(^MT-)|(^RPL)|(^RPS)|(^MRPL)", toupper(EnsDb.gene_name)) &
             !grepl("(^MT-)|(^RPL)|(^RPS)|(^MRPL)", toupper(seqc.gene_name)))
  if(pairwise){
    colnames(markers)[1] <- "gene"
    markers <- subset(markers, abs(avg_log2FC) > log2(fc.thr) & p_val_adj < apv.thr)
    markers$cluster <- ifelse(markers$avg_log2FC > 0, "up", "down")
    if(!is.infinite(n)) markers <- markers %>% group_by(cluster) %>% top_n(n, avg_log2FC)
    markers <- markers[order(markers$avg_log2FC), ]
  } else{
    markers <- subset(markers, avg_log2FC > log2(fc.thr) & p_val_adj < apv.thr)
    if(any(duplicated(markers$gene))) markers <- markers %>% group_by(gene) %>% top_n(1, avg_log2FC)
    if(!is.infinite(n)) markers <- markers %>% group_by(cluster) %>% top_n(n, avg_log2FC)
    markers <- markers[order(markers$cluster, markers$avg_log2FC), ]
  }
  return(markers)
}

DE.heatmap <- function(
    expr = sro@assays$RNA@data, cells = rep(T, nrow(md)),
    split = md$Clusters, md = sro@meta.data, markers,
    ...
){
  markers$gene_name <- ifelse(
    test = is.na(markers$cellranger.gene_name),
    no = markers$cellranger.gene_name,
    yes = ifelse(
      test = is.na(markers$EnsDb.gene_name),
      no = markers$EnsDb.gene_name,
      yes = stringr::str_to_title(markers$seqc.gene_name)
    )
  )
  
  col.ramp <- colorRamp2(c(-10, -5, seq(-3, 3, 1), 5, 10),
                         c("#173c68", "#173c68", "#1c5785", "#166db0", "#a6d1e3", "#ffffff",
                           "#f09173", "#d46052", "#b32339", "#690822", "#690822"))
  Heatmap(
    matrix = t(scale(t(expr[markers$gene, cells]))),
    name = "scaled\nimputed\nexpression", col = col.ramp,
    top_annotation = columnAnnotation(df = md[cells, c(5:7, 9:11)], col = pal),
    column_split = split[cells], column_title_gp = gpar(fontsize = 10),
    column_names_centered = T, show_column_names = F,
    cluster_columns = T, cluster_column_slices = F,
    row_names_side = "right", row_labels = markers$gene_name,
    row_names_gp = gpar(fontface = "italic", fontsize = 5), cluster_rows = F,
    use_raster = T, ...
  )
}

# ########################################################################### #
# anchor-based integration ####
# ########################################################################### #
# kedmi
sro1 <- readRDS("results/SRO.rds") %>% 
  subset(Clusters %in% c(1, 10, 12, 16, 18)) %>%
  NormalizeData() %>% FindVariableFeatures(nfeatures = 5000) %>%
  ScaleData() %>% RunPCA(npcs = 50)
# wang
sro2 <- readRDS("../GSE176282-wang-2021/results/SRO.rds") %>% 
  subset(Clusters %in% c(1, 4, 11, 20)) %>%
  NormalizeData() %>% FindVariableFeatures(nfeatures = 5000) %>%
  ScaleData() %>% RunPCA(npcs = 50)
# normalize and identify variable features for each dataset independently
sro.list <- list(sro1, sro2)

# select features that are repeatedly variable across datasets for integration
integration.features <- SelectIntegrationFeatures(object.list = sro.list, nfeatures = 5000)
anchors <- FindIntegrationAnchors(object.list = sro.list, anchor.features = integration.features)
sro.combined <- IntegrateData(anchorset = anchors)
sro.combined$Cluster.prev <- sro.combined$Clusters

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(sro.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
sro.combined <- ScaleData(sro.combined, features = rownames(sro.combined)) %>%
  RunPCA(features = rownames(sro.combined), npcs = 50) %>%
  run.PhenoGraph(npcs = 30, k = 30) %>%
  RunUMAP(dims = 1:30, n.neighbors = 30, metric = "cosine", min.dist = 0.4, spread = 1)

Idents(sro.combined) <- sro.combined$Clusters

# save files ####
DefaultAssay(sro.combined) <- "RNA"
write.csv(sro.combined@meta.data, file = "integrated-with-wang/results/meta-data.csv")
write.csv(sro.combined@reductions$umap@cell.embeddings, file = "integrated-with-wang/results/UMAP.csv")
saveRDS(sro.combined, file = "integrated-with-wang/results/SRO.rds")
write.csv(sro.combined@assays$RNA@data, file = "integrated-with-wang/results/unimputed-expr.csv")
saveRDS(as.matrix(sro.combined@assays$RNA@data), file = "integrated-with-wang/results/unimputed-expr.rds")

sro.combined <- magic(sro.combined)
write.csv(sro.combined@assays$MAGIC_RNA@data, file = "integrated-with-wang/results/imputed-expr.csv")
saveRDS(sro.combined@assays$MAGIC_RNA@data, file = "integrated-with-wang/results/imputed-expr.rds")

# plots ####
pdf("integrated-with-wang/plots/QC/cluster-QC.pdf", width = 15, height = 10)
Idents(sro.combined) <- sro.combined$Clusters
plot.all.QC(sro.combined, ident = "Clusters", col = pal$Clusters)
dev.off()

# cell.count <- rowSums(sro.combined@assays$RNA@counts > 0)
# sro.combined <- NormalizeData(sro.combined[cell.count > 1, ]) %>%
#   FindVariableFeatures(selection.method = "vst", nfeatures = 5000)

sro.combined <- readRDS("integrated-with-wang/results/SRO.rds")

# UMAP plots #####
pdf("integrated-with-wang/plots/clusters.pdf", width = 15, height = 12)
plot.groups(sro.combined)
dev.off()

ggsave(
  filename = "integrated-with-wang/plots/Sample.pdf", width = 15, height = 12,
  plot = plot.groups(
    sro = sro.combined, pref.C = F, labels = F,
    clusters = sro.combined$Sample, cl.name = "Sample"
  )
)

# integrated plots
kedmi.idx <- sro.combined$orig.ident == 'kedmi'
wang.idx <- sro.combined$orig.ident == 'wang'

df <- cbind(sro.combined@reductions$umap@cell.embeddings, sro.combined@meta.data)
eb <- element_blank()
gp <- ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  theme_classic() + theme(axis.ticks = eb, axis.text = eb, axis.title = eb)

# annotations
pdf("integrated-with-wang/plots/annotations.pdf", width = 15, height = 12)
# kedmi
cluster_labels <- df[kedmi.idx, ] %>% group_by(Cluster.annot) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)
gp + geom_point(data = df[!kedmi.idx, ], color = "grey", size = 1, alpha = 0.7) +
  geom_point(aes(color = Cluster.annot), df[kedmi.idx, ], size = 2, alpha = 0.9) +
  geom_label_repel(data = cluster_labels, 
                   aes(x = UMAP_1, y= UMAP_2, label = Cluster.annot, color=Cluster.annot), 
                   size=5, show.legend = F) +
  scale_color_manual(values = pal$Cluster.annot) + labs(color = "Kedmi.annot")
# wang
cluster_labels <- df[wang.idx, ] %>% group_by(Cluster.annot) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)
gp + geom_point(data = df[!wang.idx, ], color = "grey", size = 1, alpha = 0.7) +
  geom_point(aes(color = Cluster.annot), df[wang.idx, ], size = 2, alpha = 0.9) +
  geom_label_repel(data = cluster_labels, 
                   aes(x = UMAP_1, y= UMAP_2, label = Cluster.annot, color=Cluster.annot), 
                   size=5, show.legend = F) +
  scale_color_manual(values = pal$Cluster.annot) + labs(color = "Wang.annot")
dev.off()

# clusters
pdf("integrated-with-wang/plots/clusters.pdf", width = 15, height = 12)
plot.groups(sro.combined)
# kedmi
cluster_labels <- df[kedmi.idx, ] %>% group_by(Cluster.prev) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)
gp + geom_point(data = df[!kedmi.idx, ], color = "grey", size = 1, alpha = 0.7) +
  geom_point(aes(color = Cluster.prev), df[kedmi.idx, ], size = 2, alpha = 0.9) +
  geom_label_repel(data = cluster_labels, 
                   aes(x = UMAP_1, y= UMAP_2, label = Cluster.prev, color=Cluster.prev), 
                   size=5, show.legend = F) +
  scale_color_manual(values = pal$Clusters) + labs(color = "Kedmi.Clusters")
# wang
cluster_labels <- df[wang.idx, ] %>% group_by(Cluster.prev) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)
gp + geom_point(data = df[!wang.idx, ], color = "grey", size = 1, alpha = 0.7) +
  geom_point(aes(color = Cluster.prev), df[wang.idx, ], size = 2, alpha = 0.9) +
  geom_label_repel(data = cluster_labels, 
                   aes(x = UMAP_1, y= UMAP_2, label = Cluster.prev, color=Cluster.prev), 
                   size=5, show.legend = F) +
  scale_color_manual(values = pal$Clusters) + labs(color = "Wang.Clusters")
dev.off()

pdf("integrated-with-wang/plots/QC/cluster-QC-UMAP.pdf", width = 15, height = 12)
plot.groups(sro.combined)
plot.continuous.value(sro.combined, idx = colnames(sro.combined), vis=sro.combined@reductions$umap@cell.embeddings,
                      val = sro.combined$nCount_RNA, val.name='nCount_RNA', point.size=1)
plot.continuous.value(sro.combined, idx = colnames(sro.combined), vis=sro.combined@reductions$umap@cell.embeddings,
                      val = sro.combined$nFeature_RNA, val.name='nFeature_RNA', point.size=1)
plot.continuous.value(sro.combined, idx = colnames(sro.combined), vis=sro.combined@reductions$umap@cell.embeddings,
                      val = sro.combined$percent.MT, val.name='percent.MT', point.size=1)
dev.off()


# Breakdown bar ####
ggsave(
  plot = break.down.bar.plot(md = sro.combined@meta.data, group1 = "Clusters", group2 = "Sample",
                             size = 2, label.padding = unit(0.1, "lines")),
  filename = "integrated-with-wang/plots/breakdown-barplots/Samples-in-Clusters.pdf", width = 10, height = 10
)

