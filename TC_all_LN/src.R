# TC all LN ####
setwd("~/0-workspace/CCR7_DC/TC_all_LN/")

library(Seurat)
library(Rmagic)
library(EnsDb.Mmusculus.v79)
library(anndata)
library(Matrix)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(dplyr)

set.seed(1)
options(future.globals.maxSize = Inf)

# pal ####
# adding a new colour for the transitional as we haven't had that before in previous publications - something very pale is fine!
# pal <- list(
#   sample = c(
#     'colonic' = "#bf0000", 'hepatic' = "#ccf3ff",
#     'mediastinal' = "#ff8400", "para-aortic" = "#d1c50f",
#     "Peyers patches" = "#2e2bed" ,"salivary gland" = "#ff19a0", "skin" = "#11d100","small intestine" = '#00b0cf'
#   ),
#   Clusters = c(
#     "#e60e0e", "#077315", "#1ee3c5", "#f081e6", "#7c2da6",
#     "#de9309", "#d1c50f", "#489de8", "#2e2bed", "#5fed0e",
#     "#bd537a", "#2f026e", "#1f758c", "#9c6fa8", "#6b6580",
#     "#9e7d3f", "#49709c", "#ccf3ff", "#5c1a1a", "#d1d18a",
#     "#c9a6a1", "#827c68", "#b54800", "#79a695"
#   ),
#   hash_id = c('B0301' = "#3322A6", 'B0302'="#92a1f7", 'B0303' = "#ff8400",
#               'B0304'="#f51b8b", 'B0305' ="#e080b5", "B0306" =  "#1ee3c5",
#               "B0307" = "#5fed0e", "B0308" = "#4f5715",
#               'Negative' = '#d1c50f'
#   ),
#   # Broad.Annotation = c(
#   #   `Ki67+ TC` = "#e60e0e",
#   #   `TC I` = "#fc12f5",
#   #   `TC II` = "#BC23FF",
#   #   `TC III` = "#72418F",
#   #   `TC IV` = "#2e2bed",
#   #   `NA` = "lightgrey"
#   # ),
# 
#   # Annotation = c(
#   #   `Proliferating Ki67+ TC I` = "#e051bc",
#   #   `TC I` = "#D790FF",
#   #   `TC II` = "#BC23FF",
#   #   `TC III` = "#72418F",
#   #   `TC IV` = "#3A2465",
#   #   `Ex TC IV` = "#7c7191",
#   #   `ILC3 - TC` = "#7c7feb",
#   #   `Proliferating Ki67+ ILC3` = "#0976de",
#   #   `Ncr1+ ILC3` = "#003a85",
#   #   `LTi` = "#14a38e"
#   # ),
#   Annotation = c(
#     `early/transitional TC` = '#ccf3ff',
#     `Ki67+ TC` = "#e051bc",
#     `TC I` = "#D790FF",
#     `TC II` = "#BC23FF",
#     `TC III` = "#72418F",
#     `TC IV` = "#3A2465",
#     `DC contaminant` = "#e8c91c",
#     `ILC3` = "#de9309",
#     `low QC` = "#827c68"
#   )
# )
# annotations <- c('TC III', "TC II", "TC I", "early/transitional TC", "TC IV",
#                  "ILC3", "Ki67+ TC",
#                  "DC contaminant", "low QC")
# annotations <- c("Ki67+ TC", "early/transitional TC", 'TC I', "TC II", "TC III", "TC IV",
#                  "ILC3", "DC contaminant", "low QC")
# pal$Clusters <- pal$Annotation[annotations]
# names(pal$Clusters) <- 1:length(pal$Clusters)
# pal$gut <- c("gut" = "blue", "non-gut" = "red")
# saveRDS(pal, file = "plots/palette.rds")
pal <- readRDS(file = "plots/palette.rds")

plot.QC.violin <- function(sro, ident, feature, yintercept = NULL, br, Log10 = F, pal = NULL){
  if(feature == "nCount_RNA"){ name <- "number of transcripts"
  } else if(feature == "nFeature_RNA"){ name <- "number of detected genes"
  } else if(feature == "percent.MT"){ name <- "percentage of mitochondrial transcripts"
  } else if(feature == "S.Score"){ name <- "S-phase score"
  } else if(feature == "G2M.Score") name <- "G2/M-phase score"
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
    theme_classic() +
    labs(x = ident, y = "", title = ifelse(Log10, paste0("log10(", name, ")"), name))
  if(!is.null(pal))
    gp <- gp + scale_color_manual(values = pal) + scale_fill_manual(values = pal)
  return(gp)
}

plot.all.QC <- function(sro, ident, thr = NULL, col = NULL, CC = F){
  L <- labs(color = ident, x = "number of transcripts", y = "number of detected genes")
  gp1 <- FeatureScatter(sro, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme_classic() + L
  if(!is.null(col)){
    names(col) <- levels(sro@meta.data[, ident])
    gp1 <- gp1 + scale_color_manual(values = col)
  }
  gp2 <- ggplot(sro@meta.data, aes(x = nCount_RNA, y = nFeature_RNA)) + geom_density_2d_filled() + theme_classic() + L
  pl <- list(
    plot.QC.violin(sro, ident = ident, feature = "nCount_RNA", Log10 = T, yintercept = c(thr$nc.min, thr$nc.max), br = 0.2, pal = col),
    plot.QC.violin(sro, ident = ident, feature = "nFeature_RNA", yintercept = c(thr$nf.min, thr$nf.max), br = 500, pal = col),
    plot.QC.violin(sro, ident = ident, feature = "percent.MT", yintercept = thr$mp.max, br = 2, pal = col)
  )
  if(is.null(thr)){
    pl[[4]] <- gp1
    pl[[5]] <- gp2
  } else {
    hl <- geom_hline(yintercept = thr$nf.min, linetype = 2, color = "violetred4")
    vl <- geom_vline(xintercept = thr$nc.min, linetype = 2, color = "violetred4")
    pl[[4]] <- gp1 + hl + vl
    pl[[5]] <- gp2 + hl + vl
  }
  if (CC) {
    pl <- c(pl, list(
      plot.QC.violin(sro, ident = ident, feature = "S.Score", br = 0.2, pal = col),
      plot.QC.violin(sro, ident = ident, feature = "G2M.Score", br = 0.2, pal = col)
    ))
  }
  return(pl)
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
    gp + geom_label(aes(x = vis.cl[, 1], vis.cl[, 2], color = C, label = CL), color = 'black',
                    size = 5, show.legend = F, nudge_y = nudge_y) + guides(colour = guide_legend(override.aes = list(size=5)))
      # geom_label(aes(x = vis.cl[, 1], vis.cl[, 2], color = C, label = CL), color = 'black', label.size = NA,
      #            size = 5, show.legend = F, nudge_y = nudge_y) + guides(colour = guide_legend(override.aes = list(size=5)))
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

# Basic processing ####
dir.create("plots/QC", recursive=T)
dir.create("results", recursive=T)

## Create sro ####
counts <- Read10X_h5('data/filtered_feature_bc_matrix.h5') # gene by cell matrix
rownames(counts) <- toupper(rownames(counts))

dmx.cl <- read_h5ad("data/adataFinal_CB-2623_TC_all_LN.h5ad")
rownames(dmx.cl[["var"]]) <- gsub("-.+$", "", rownames(dmx.cl[["var"]]))
rownames(dmx.cl[["obs"]]) <- paste0(dmx.cl[["obs"]]$barcode_sequence, "-1")
obs <- dmx.cl[["obs"]][colnames(counts), ]
obs$sample <- dmx.cl[["var"]][obs$hash_id, "feature_name"]

sro <- CreateSeuratObject(
  project = "TC_all_LN",
  counts = counts, meta.data = obs
)
Idents(sro) <- sro$sample

## sample QC ####
sro <- PercentageFeatureSet(sro, pattern = "^MT.", col.name = "percent.MT")
thr <- data.frame(nc.min = 1500, nf.min = 500, mp.max = 5) # filter by number of transcripts, genes, and mt
cell.discard <- sro$nCount_RNA < thr$nc.min |
  sro$nFeature_RNA < thr$nf.min |
  sro$percent.MT > thr$mp.max |
  sro$hash_id %in% c("Negative", "Doublet") | is.na(sro$hash_id)

Idents(sro) <- sro$sample
pdf("plots/QC/sample-QC.pdf", width = 10, height = 10)
plot.all.QC(sro, ident = "Sample", thr = thr)
dev.off()

Idents(sro) <- sro$hash_id
pdf("plots/QC/hashID-QC.pdf", width = 10, height = 10)
plot.all.QC(sro, ident = "hashID", thr = thr)
dev.off()

## cluster and umap ####
sro <- sro[, !cell.discard]

cell.count <- rowSums(sro@assays$RNA@counts > 0)
sro <- NormalizeData(sro[cell.count > 1, ]) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 5000)
sro <- ScaleData(sro, features = rownames(sro)) %>%
  RunPCA(features = rownames(sro), npcs = 50) %>%
  FindNeighbors(dims = 1:30, k.param = 30) %>%
  FindClusters(resolution = c(0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0)) %>%
  RunUMAP(dims = 1:30, n.neighbors = 30, metric = "cosine", min.dist = 0.4, spread = 1)

sro$Clusters <- factor(as.integer(sro$RNA_snn_res.0.6))
Idents(sro) <- sro$Clusters

s.genes <- rownames(sro)[toupper(rownames(sro)) %in% toupper(cc.genes.updated.2019$s.genes)]
g2m.genes <- rownames(sro)[toupper(rownames(sro)) %in% toupper(cc.genes.updated.2019$g2m.genes)]
sro <- CellCycleScoring(sro, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
sro$Phase <- gsub("G2M", "G2/M", sro$Phase)

## cluster QC ####
Idents(sro) <- sro$Clusters
pdf("plots/QC/cluster-QC.pdf", width = 10, height = 10)
plot.all.QC(sro, ident = "Clusters", col = pal$Clusters)
dev.off()

pdf("plots/QC/cluster-QC-UMAP.pdf", width = 15, height = 12)
plot.continuous.value(sro, idx = rownames(sro@meta.data), vis = sro@reductions$umap@cell.embeddings,
                      val = sro$nCount_RNA, val.name='nCount_RNA', point.size=1)
plot.continuous.value(sro, idx = rownames(sro@meta.data), vis = sro@reductions$umap@cell.embeddings,
                      val = sro$nFeature_RNA, val.name='nFeature_RNA', point.size=1)
plot.continuous.value(sro, idx = rownames(sro@meta.data), vis = sro@reductions$umap@cell.embeddings,
                      val = sro$percent.MT, val.name='percent.MT', point.size=1)
plot.continuous.value(sro, idx = rownames(sro@meta.data), vis = sro@reductions$umap@cell.embeddings,
                      val = sro$S.Score, val.name='S.Score', point.size=1)
plot.continuous.value(sro, idx = rownames(sro@meta.data), vis = sro@reductions$umap@cell.embeddings,
                      val = sro$G2M.Score, val.name='G2M.Score', point.size=1)
dev.off()
# ############################################################################ #
# add annotations ####
sro$sample <- factor(sro$sample, 
                     levels = c("colonic", "small intestine", "hepatic",
                                "Peyers patches", "mediastinal", "skin",
                                "salivary gland", "para-aortic"))
# 0 = TC III,
# 1 = TC II,
# 3 = early/transitional TC,
# 4 = TC IV,
# 5 = ILC3
# 6 = Ki67+ TC, and
# cluster 7 (DC contaminant),
# cluster 8 (low QC),
# annot <- openxlsx::read.xlsx(
#   xlsxFile = "Foxn1creTap63 annotations.xlsx",
#   rowNames = T, colNames = T, check.names = F
# )
# annot$Cluster.annotation <- annot$Cluster.annotation %>%
#   factor(levels = c('Aire+', 'IRF8hi', 'Tuft', 'Klk enriched', 'CCL21+ mTEC',
#     'Ki67+ TAA', 'Lung goblet (top), Enterohepatic bottom',
#     'Keratinocyte', 'Neurendocrine', 'M cell', 'Ciliated', 'cTEC',
#     'low QC', 'low QC?', 'Tuft ?low QC'))
annotations <- c('TC III', "TC II", "TC I", "early/transitional TC", "TC IV", 
                 "ILC3", "Ki67+ TC", 
                 "DC contaminant", "low QC")
sro$Annotation <- annotations[sro$Clusters] %>% 
  factor(levels = c("Ki67+ TC", "early/transitional TC", "TC I", "TC II", 'TC III', "TC IV", "ILC3", 
                    "DC contaminant", "low QC"))

## re-order clusters
cl.ord <- c(7, 4, 3, 2, 1, 5, 6, 8, 9)
sro$Clusters2 <- factor(match(sro$Clusters, cl.ord))

## Save results ####
saveRDS(sro, file = "results/SRO.rds")
write.csv(sro@meta.data, file = "results/meta-data.csv")
write.csv(sro@reductions$pca@cell.embeddings, file = "results/PCA.csv")
write.csv(sro@reductions$umap@cell.embeddings, file = "results/UMAP.csv")
saveRDS(as.matrix(sro@assays$RNA@data), file = "results/unimputed-expr.rds")

sro.imp <- magic(sro)
saveRDS(sro.imp@assays$MAGIC_RNA@data, file = "results/imputed-expr.rds")

# Plot ####
sro <- readRDS("results/SRO.rds")
imp.expr <- readRDS("results/imputed-expr.rds")

## UMAP ####
ggsave(
  filename = "plots/Sample-UMAP.pdf", width = 15, height = 12,
  plot = plot.groups(
    sro = sro, pref.C = F, labels = T, vis = sro@reductions$umap@cell.embeddings,
    clusters = sro$sample, cl.name = "Sample", col = pal$sample
  )
)

ggsave(
  filename = "plots/hash_id-UMAP.pdf", width = 15, height = 12,
  plot = plot.groups(
    sro = sro, pref.C = F, labels = T, vis = sro@reductions$umap@cell.embeddings,
    clusters = sro$hash_id, cl.name = "hash_id", col = pal$hash_id
  )
)

# one that has the ILC3 included and another that has no ILC3?
pdf("plots/Clusters-UMAP.pdf", width = 15, height = 12)
plot.groups(
  sro = sro, pref.C = T, labels = T, vis = sro@reductions$umap@cell.embeddings,
  clusters = sro$Clusters2, cl.name = "Clusters", col = pal$Clusters
)
plot.groups(
  sro = sro, pref.C = T, labels = T, vis = sro@reductions$umap@cell.embeddings,
  idx = !(sro$Annotation %in% c("DC contaminant", "low QC")),
  clusters = sro$Clusters2, cl.name = "Clusters", col = pal$Clusters
)
plot.groups(
  sro = sro, pref.C = T, labels = T, vis = sro@reductions$umap@cell.embeddings,
  idx = !(sro$Annotation %in% c("ILC3", "DC contaminant", "low QC")),
  clusters = sro$Clusters2, cl.name = "Clusters", col = pal$Clusters
)
dev.off()

# ggsave(
#   filename = "plots/Clusters-UMAP.pdf", width = 15, height = 12,
#   plot = plot.groups(
#     sro = sro, pref.C = T, labels = T, vis = sro@reductions$umap@cell.embeddings,
#     clusters = sro$Clusters, cl.name = "Clusters", col = pal$Clusters2
#   )
# )

# one that has the ILC3 included and another that has no ILC3?
pdf("plots/Annotations-UMAP.pdf", width = 15, height = 12)
plot.groups(
  sro = sro, pref.C = F, labels = T, vis = sro@reductions$umap@cell.embeddings,
  clusters = sro$Annotation, cl.name = "Annotation", col = pal$Annotation
)
plot.groups(
  sro = sro, pref.C = F, labels = T, vis = sro@reductions$umap@cell.embeddings,
  idx = !(sro$Annotation %in% c("DC contaminant", "low QC")),
  clusters = sro$Annotation, cl.name = "Annotation", col = pal$Annotation
)
plot.groups(
  sro = sro, pref.C = F, labels = T, vis = sro@reductions$umap@cell.embeddings,
  idx = !(sro$Annotation %in% c("ILC3", "DC contaminant", "low QC")),
  clusters = sro$Annotation, cl.name = "Annotation", col = pal$Annotation
)
dev.off()
# ggsave(
#   filename = "plots/Annotations-UMAP.pdf", width = 15, height = 12,
#   plot = plot.groups(
#     sro = sro, pref.C = F, labels = T, vis = sro@reductions$umap@cell.embeddings,
#     clusters = sro$Annotation, cl.name = "Annotation", col = pal$Annotation
#   )
# )

## Aire expr
plot_grid(
  plot.continuous.value(
    sro, vis = sro@reductions$umap@cell.embeddings,
    idx = rownames(sro@meta.data), point.size = 0.5, 
    val = sro@assays$RNA@data["AIRE", ], val.name = "AIRE"
  ) + ggtitle("log-normalized gene expression"),
  plot.continuous.value(
    sro, vis = sro@reductions$umap@cell.embeddings,
    idx = rownames(sro@meta.data), point.size = 0.5, 
    val = imp.expr["AIRE", ], val.name = "AIRE"
  ) + ggtitle("imputed, log-normalized gene expression")
) %>% ggsave(filename = "plots/AIRE-expr.pdf", width = 12, height = 5)

plot_grid(
  plot.continuous.value(
    sro, vis = sro@reductions$umap@cell.embeddings,
    idx = rownames(sro@meta.data), point.size = 0.5, 
    val = sro@assays$RNA@data["AIRE", ], val.name = "RALDH2"
  ) + ggtitle("log-normalized gene expression"),
  plot.continuous.value(
    sro, vis = sro@reductions$umap@cell.embeddings,
    idx = rownames(sro@meta.data), point.size = 0.5, 
    val = imp.expr["AIRE", ], val.name = "RALDH2"
  ) + ggtitle("imputed, log-normalized gene expression")
) %>% ggsave(filename = "plots/RALDH2-expr.pdf", width = 12, height = 5)

## Breakdown ####
# so for the bar plots we can then show the proportions of TC clusters (not including ILC3) across the samples, 
# and we can group them: colonic, small intestine, hepatic, Peyer's patches, mediastinal, skin, salivary gland, para-aortic
# from left to right
# Ki67+, transitional, I, II, III, IV

pdf("plots/breakdown-barplots-1.pdf", width = 10, height = 8)
break.down.bar.plot(sro@meta.data[!(sro@meta.data$Annotation %in% c("ILC3", "DC contaminant", "low QC")),], 
                    group1 = "sample", group2 = "Annotation")
dev.off()
pdf("plots/breakdown-barplots-2.pdf", width = 15, height = 15)
break.down.bar.plot(sro@meta.data[!(sro@meta.data$Annotation %in% c("ILC3", "DC contaminant", "low QC")),], group1 = "Annotation", group2 = "sample")
break.down.bar.plot(sro@meta.data[!(sro@meta.data$Annotation %in% c("ILC3", "DC contaminant", "low QC")),], group1 = "sample", group2 = "Annotation")
dev.off()

table(sro$sample)

# ############################################################################ #
# Redo UMAP and plot clusters ####
sro.subset <- subset(sro, Annotation %in% c("Ki67+ TC", "early/transitional TC", "TC I", "TC II", 'TC III', "TC IV")) %>%
  # FindNeighbors(dims = 1:30, k.param = 30) %>%
  # FindClusters(resolution = c(0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0)) %>%
  RunUMAP(dims = 1:30, n.neighbors = 30, metric = "cosine", min.dist = 0.4, spread = 0.7)

## re-order clusters in sro.subset
sro.subset <- readRDS("results/SRO_subset.rds")
sro.subset$Clusters2 <- sro@meta.data[colnames(sro.subset), "Clusters2"]
sro.subset$Clusters2 <- droplevels(sro.subset$Clusters2)

saveRDS(sro.subset, file = "results/SRO_subset.rds")
write.csv(sro.subset@reductions$umap@cell.embeddings, file = "results/UMAP-subset.csv")

pdf("plots/Annotations-UMAP-subset.pdf", width = 15, height = 12)
plot.groups(
  sro = sro.subset, pref.C = F, labels = T, vis = sro.subset@reductions$umap@cell.embeddings,
  idx = rep(T, ncol(sro.subset)),
  clusters = sro.subset$Annotation, cl.name = "Annotation", col = pal$Annotation
)
dev.off()

pdf("plots/Clusters-UMAP-subset.pdf", width = 15, height = 12)
plot.groups(
  sro = sro.subset, pref.C = T, labels = T, vis = sro.subset@reductions$umap@cell.embeddings,
  idx = rep(T, ncol(sro.subset)),
  clusters = sro.subset$Clusters2, cl.name = "Clusters", col = pal$Clusters
)
dev.off()

## Expr overlay ####
sro.subset <- readRDS("results/SRO_subset.rds")
imp.expr <- readRDS("results/imputed-expr.rds")
geneList <- c('CD36', 'MGLL', "PRDM16")
geneList <- c("RALDH2")
geneList <- c("ALDH1A2"); fn <- "RALDH2"
"RALDH2" %in% rownames(imp.expr)
geneList <- geneList[geneList %in% rownames(imp.expr)]

if (length(geneList) >= 5){
  n.col <- 5} else{
    n.col <- length(geneList)}
marker.p <- lapply(geneList, function(gene){
  p <- plot.continuous.value(
    sro.subset, vis = sro.subset@reductions$umap@cell.embeddings,
    idx = rownames(sro.subset@meta.data), point.size = 0.5,
    val = imp.expr[gene, ], val.name = 'imputed\nexpression'
  ) + ggtitle(gene)
  return(p)
})

marker.p.unimputed <- lapply(geneList, function(gene){
  p <- plot.continuous.value(
    sro.subset, vis = sro.subset@reductions$umap@cell.embeddings,
    idx = rownames(sro.subset@meta.data), point.size = 0.5,
    val = sro.subset@assays$RNA$data[gene, ], val.name = 'expression'
  ) + ggtitle(gene)
  return(p)
})

plot_grid(plotlist = marker.p, ncol = n.col) %>% 
  ggsave(filename = paste0("plots/gene-expr/UMAP-gene-expr-imputed-", fn, ".pdf"), 
         width = n.col*6, height = ceiling(length(marker.p)/n.col)*4)
plot_grid(plotlist = marker.p.unimputed, ncol = n.col) %>% 
  ggsave(filename = paste0("plots/gene-expr/UMAP-gene-expr-unimputed-", fn, ".pdf"), 
         width = n.col*6, height = ceiling(length(marker.p.unimputed)/n.col)*4)

## Violin ####
# a violin plot for these genes grouping LNs into gut vs non gut (para-aortic, mediastinal, skin) 
# as I wonder whether the lipid genes might be an adaptation to the gut environment? 
# Although confounded by the fact that tc IV has most of these genes and they are low in non gut
table(sro.subset$sample)
VlnPlot(sro.subset, features = geneList, group.by = "sample") %>% 
  ggsave(filename = paste0("plots/gene-expr/violin-CD36_MGLL_PRDM16.pdf"),
         width = 18, height = 4)
sro.subset$gut <- ifelse(test = sro.subset$sample %in% c("para-aortic", "mediastinal", "salivary gland", "skin"),
                         yes = "non-gut",
                         no = "gut")

VlnPlot(sro.subset, features = geneList, group.by = "gut", cols = pal$gut) %>% 
  ggsave(filename = paste0("plots/gene-expr/violin-CD36_MGLL_PRDM16-gut.pdf"),
         width = 12, height = 4)
VlnPlot(sro.subset, features = geneList, group.by = "sample", cols = pal$sample) %>% 
  ggsave(filename = paste0("plots/gene-expr/violin-CD36_MGLL_PRDM16-sample.pdf"),
         width = 18, height = 4)
VlnPlot(sro.subset, features = geneList, group.by = "Annotation", cols = pal$Annotation) %>% 
  ggsave(filename = paste0("plots/gene-expr/violin-CD36_MGLL_PRDM16-Annotation.pdf"),
         width = 18, height = 4)
# ############################################################################ #
# fgsea ####
library(fgsea)
library(msigdbr)
sro.subset <- readRDS("results/SRO_subset.rds")
annot <- unique(sro.subset@meta.data[, c("Clusters", "Annotation")])
degs <- read.csv("results/markers/TC-clusters.csv")
degs$rank <- degs$avg_log2FC
degs$feature <- degs$gene
degs$group <- annot[match(degs$cluster, annot$Clusters), "Annotation"]
degs$group <- droplevels(degs$group)

genesets.hallmark = msigdbr(species = "Mus musculus", category = "H")
genesets.hallmark$gene_symbol <- toupper(genesets.hallmark$gene_symbol)
pathways.hallmark = split(x = genesets.hallmark$gene_symbol, f = genesets.hallmark$gs_name)
genesets.reactome = msigdbr(species = "Mus musculus", subcategory = 'CP:REACTOME')
genesets.reactome$gene_symbol <- toupper(genesets.reactome$gene_symbol)
pathways.reactome = split(x = genesets.reactome$gene_symbol, f = genesets.reactome$gs_name)

genesets.go.bp = msigdbr(species = "Mus musculus", subcategory = 'GO:BP')
genesets.go.bp$gene_symbol <- toupper(genesets.go.bp$gene_symbol)
pathways.go.bp = split(x = genesets.go.bp$gene_symbol, f = genesets.go.bp$gs_name)


pathway.list <- list(
  "Hallmark" = pathways.hallmark,
  'Reactome' = pathways.reactome,
  'GO-Biological_Process' = pathways.go.bp)

pathway.list <- list('GO-Biological_Process' = pathways.go.bp)
  # "C2_CGP" = pathways.c2.cgp,
  # 'Reactome' = pathways.reactome,
  # 'KEGG' = pathways.kegg,
  # 'GO-Biological_Process' = pathways.go.bp,
  # 'GO-Cellular_Component' = pathways.go.cc,
  # 'GO-Molecular_Function' = pathways.go.mf

pref.p <- "plots/fgsea/avg_log2FC/"; pref <- "results/fgsea/avg_log2FC/"
group.name <- "TC-clusters"

for (grp in unique(degs$group)){
  for (j in 1:length(pathway.list)){
    curr.pathway <- pathway.list[[j]]
    curr.pathway.name <- names(pathway.list)[j]
    print(paste(grp, curr.pathway.name, sep=", "))
    
    marker.genes <- degs %>%
      dplyr::filter(group == as.character(grp)) %>%
      arrange(desc(rank)) %>% 
      dplyr::select(feature, rank)

    cluster.name <- paste0(as.character(grp))
    ranks <- tibble::deframe(marker.genes)
    # test
    fgseaRes <- fgsea(pathways=curr.pathway, stats=ranks, nproc=8)
    fgseaResTidy <- fgseaRes %>%
      as_tibble() %>%
      arrange(padj)
    
    # res.table.pathway <- curr.pathway %>%
    #   tibble::enframe("pathway", "feature") %>%
    #   tidyr::unnest(cols = c("feature")) %>%
    #   inner_join(marker.genes, by="feature")
    
    topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=20), pathway]
    topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=20), pathway]
    # topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    topPathways <- topPathwaysUp
    waterfall.top20 <- plotGseaTable(curr.pathway[topPathways], ranks, fgseaRes, 
                                     gseaParam=0.5)
    collapsedPathways <- collapsePathways(fgseaRes[ES > 0][order(pval)][padj < 0.05],
                                          curr.pathway, ranks)
    mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
      order(-NES), pathway]
    waterfall.main <- plotGseaTable(curr.pathway[mainPathways], ranks, fgseaRes,
                                    gseaParam = 0.5)
    
    bar.p <- ggplot(fgseaResTidy[fgseaResTidy$pathway %in% topPathways,], aes(reorder(pathway, NES), NES)) +
      geom_col(aes(fill=padj<0.05)) +
      coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score",
           title=curr.pathway.name) + scale_fill_manual(values = c("TRUE" = '#00BFC4', "FALSE" = '#F8766D')) +
      theme_minimal()
    bar.collapsed.p <- ggplot(fgseaResTidy[fgseaResTidy$pathway %in% mainPathways,], aes(reorder(pathway, NES), NES)) +
      geom_col(aes(fill=padj<0.05)) +
      coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score",
           title=curr.pathway.name) + scale_fill_manual(values = c("TRUE" = '#00BFC4', "FALSE" = '#F8766D')) +
      theme_minimal()
    
    dir.create(paste0(pref.p, group.name, '/', curr.pathway.name), recursive=T, showWarnings = F)
    dir.create(paste0(pref, group.name, '/', curr.pathway.name), recursive=T, showWarnings = F)
    
    ggsave(plot = waterfall.top20,
           file = paste0(pref.p, group.name, '/', curr.pathway.name, '/', gsub("/", "_", cluster.name), '-waterfall-top20.pdf'),
           width = 16, height = 16)
    # ggsave(plot = waterfall.main,
    #        file = paste0(pref.p, group.name, '/', curr.pathway.name, '/', gsub("/", "_", cluster.name), '-waterfall-collapsed_main.pdf'),
    #        width = 16, height = 10)
    
    ggsave(plot = bar.p, 
           file = paste0(pref.p, group.name, '/', curr.pathway.name, '/', gsub("/", "_", cluster.name), '-bar-top20.pdf'),
           width = 16, height = 6)
    # ggsave(plot = bar.collapsed.p,
    #        file = paste0(pref.p, group.name, '/', curr.pathway.name, '/', gsub("/", "_", cluster.name), '-bar-collapsed_main.pdf'),
    #        width = 12, height = 6)
    # 
    data.table::fwrite(fgseaResTidy, 
                       file=paste0(pref, group.name, '/', curr.pathway.name, '/', gsub("/", "_", cluster.name), '.tsv'),
                       sep="\t", sep2=c("", ",", ""))
    # data.table::fwrite(res.table.pathway, 
    #                    file=paste0(pref, group.name, '/', curr.pathway.name, '/', cluster.name, '-res.table.tsv'),
    #                    sep="\t", sep2=c("", ",", ""))
  }
}

## search and save lipid or metabolism
pathway.list <- list("Hallmark", 'Reactome', 'GO-Biological_Process')
padj.threshold <- 0.1

# per group
# for (grp in unique(degs$group)){
#   for (curr.pathway.name in pathway.list){
#     group <- gsub("/", "_", grp)
#     res <- data.table::fread(file = paste0(pref, group.name, "/", curr.pathway.name, "/", group, ".tsv"), 
#                              sep = "\t", sep2 = c("", ",", "")) %>% subset(ES > 0)
#     res <- res[res$padj < padj.threshold,]
#     res.lm <- res[grep("LIPID|METABOLISM", res$pathway),]
#     data.table::fwrite(res.lm, 
#                        file=paste0(pref, group.name, '/', curr.pathway.name, '/', gsub("/", "_", group), '-lipid_or_metabolism-padj-', padj.threshold, '.tsv'),
#                        sep="\t", sep2=c("", ",", ""))
#   }
# }

for (curr.pathway.name in pathway.list){
  res2 <- data.frame()
  for (grp in unique(degs$group)){
    group <- gsub("/", "_", grp)
    res <- data.table::fread(file = paste0(pref, group.name, "/", curr.pathway.name, "/", group, ".tsv"), 
                             sep = "\t", sep2 = c("", ",", "")) %>% subset(ES > 0)
    res <- res[res$padj < padj.threshold,]
    res.lm <- res[grep("LIPID|METABOLISM", res$pathway),]
    if (nrow(res.lm) > 0){
      res.lm <- cbind('group' = grp, res.lm)
      res2 <- rbind(res2, res.lm)
    }
  }
  data.table::fwrite(res2,
                     file=paste0(pref, group.name, '/', curr.pathway.name, '/', 'lipid_or_metabolism-padj-', padj.threshold, '.tsv'),
                     sep="\t", sep2=c("", ",", ""))
}

# all
for (curr.pathway.name in pathway.list){
  res2 <- data.frame()
  for (grp in unique(degs$group)){
    group <- gsub("/", "_", grp)
    res <- data.table::fread(file = paste0(pref, group.name, "/", curr.pathway.name, "/", group, ".tsv"), 
                             sep = "\t", sep2 = c("", ",", "")) %>% subset(ES > 0)
    res.lm <- res[res$padj < padj.threshold,]
    if (nrow(res.lm) > 0){
      res.lm <- cbind('group' = grp, res.lm)
      res2 <- rbind(res2, res.lm)
    }
  }
  data.table::fwrite(res2,
                     file=paste0(pref, group.name, '/', curr.pathway.name, '/', 'padj-', padj.threshold, '.tsv'),
                     sep="\t", sep2=c("", ",", ""))
}


# TC IV lipid related GO:BP TERMS
DE.heatmap <- function(
    sro = NULL, expr = sro@assays$RNA@data, md = sro@meta.data,
    cells = rownames(md), split = "Clusters", rowsplit = NULL,
    cluster_columns = T, cluster_column_slices = T, markers, cluster_row_slices = F,
    legend.name = "scaled\nimputed\nexpression", 
    split.columns = c("Sample", "Clusters", 'celltype'),
    ca.col = list(Sample = pal$sample, Clusters = pal$Clusters, celltype=pal$celltype),
    ...
){
  ca <- columnAnnotation(
    df = md[cells, split.columns],
    col = ca.col
  )
  col.ramp <- colorRamp2(c(-10, -5, seq(-3, 3, 1), 5, 10),
                         c("#173c68", "#173c68", "#1c5785", "#166db0", "#a6d1e3", "#ffffff",
                           "#f09173", "#d46052", "#b32339", "#690822", "#690822"))
  Heatmap(
    matrix = t(scale(t(expr[markers, cells]))), name = legend.name, col = col.ramp,
    show_column_names = F, column_split = md[cells, split], top_annotation = ca,
    cluster_rows = T, row_names_side = "right", 
    row_labels = markers, row_split = rowsplit, cluster_row_slices = cluster_row_slices, 
    cluster_columns = cluster_columns, cluster_column_slices = cluster_column_slices,
    row_title_rot = 0, column_title_rot = 45,
    row_names_gp = gpar(fontface = "italic", fontsize = 12), use_raster = T
  )
}

group <- "TC IV"
curr.pathway.name <- "GO-Biological_Process"
res <- data.table::fread(file = paste0(pref, group.name, "/", curr.pathway.name, "/", group, ".tsv"),
                         sep = "\t", sep2 = c("", ",", "")) %>% subset(ES > 0)
res2 <- res[res$pathway %in% c("GOBP_STEROL_TRANSPORT", "GOBP_INTRACELLULAR_LIPID_TRANSPORT")]
# data.table::fwrite(res2,
#                    file=paste0(pref, group.name, '/', curr.pathway.name, '/', 'TC IV-GOBP_terms.tsv'),
#                    sep="\t", sep2=c("", ",", ""))

genes <- unique(unname(unlist(sapply(res2$leadingEdge, FUN = stringr::str_split, pattern = ","))))
sro.subset@meta.data$Annotation <- droplevels(sro.subset@meta.data$Annotation)
pdf(width = 16, height = 8, file = paste0("plots/fgsea/avg_log2FC/TC-clusters/GO-Biological_Process/heatmap-GOBP_terms.pdf"))
DE.heatmap(sro.subset, expr = imp.expr,
           markers = genes,
           split = 'Annotation',
           split.columns = c("sample", "Annotation"), cluster_column_slices = F,
           ca.col = list(sample = pal$sample, "Annotation" = pal$Annotation))
dev.off()

# markers from a GMT file
loadGMT <- function(gmtFile) {
  lines <- readLines(gmtFile)  # Read lines from the GMT file
  geneSetList <- lapply(lines, function(line) {
    parts <- strsplit(line, "\t")[[1]]  # Split each line by tab
    df <- data.frame(name = parts[1], description = parts[2], genes = parts[-c(1,2)])
    return(df)
  })
  geneSetDF <- do.call(rbind, geneSetList)
  return(geneSetDF)
}
geneSets <- loadGMT("20240214_genesets_pruned.gmt")

# get a df from gobp results
gobp.dfs <- apply(res2, MARGIN = 1, FUN = function(x){
  df <- data.frame(name = x['pathway'],
                   description = x['pathway'],
                   genes = stringr::str_split(x['leadingEdge'], pattern = ',')[[1]])
  return(df)
})
gobp.df <- do.call(rbind, gobp.dfs)

## all lipid genes ####
all.genes <- rbind(gobp.df, geneSets)
if (any(duplicated(all.genes$genes))) all.genes <- all.genes %>% group_by(genes) %>% top_n(1, name)
all.genes <- all.genes[all.genes$genes %in% rownames(imp.expr),]

sro.subset@meta.data$Annotation <- droplevels(sro.subset@meta.data$Annotation)

unique(sro.subset$Annotation)[1]
idx <- rownames(sro.subset@meta.data[sro.subset@meta.data$Annotation %in% c("TC I", "TC II", "TC III", "TC IV"), ])
pdf(width = 16, height = 300, file = paste0("plots/heatmaps/0-heatmap-GOBP-all_lipid_metabolism-TC_only.pdf"))
DE.heatmap(sro.subset, expr = imp.expr, cells = idx, 
           markers = all.genes$genes,
           split = 'Annotation', cluster_row_slices = T,
           rowsplit = all.genes$name,
           split.columns = c("sample", "Annotation"), cluster_column_slices = F,
           ca.col = list(sample = pal$sample, "Annotation" = pal$Annotation))
dev.off()

## each gene set separately ####
nGenes <- geneSets %>% group_by(name) %>% summarise(count = n())
summary(nGenes$count)
lapply(unique(geneSets$name), function(pathway){
  genes.df <- rbind(gobp.df, geneSets[geneSets$name == pathway,])
  genes.df <- genes.df[genes.df$genes %in% rownames(imp.expr),]
  if (any(duplicated(genes.df$genes))) genes.df <- genes.df %>% group_by(genes) %>% top_n(1, name)
  fig.h <- round(nrow(genes.df)*0.5)
  pdf(width = 16, height = fig.h, file = paste0("plots/heatmaps/heatmap-GOBP-", pathway, ".pdf"))
  print(DE.heatmap(sro.subset, expr = imp.expr,
             markers = genes.df$genes,
             split = 'Annotation', cluster_row_slices = T,
             rowsplit = genes.df$name,
             split.columns = c("sample", "Annotation"), cluster_column_slices = F,
             ca.col = list(sample = pal$sample, "Annotation" = pal$Annotation)))
  dev.off()
})

# Cytokine signatures ####
cytokines.df <- openxlsx::read.xlsx(
  xlsxFile = "Cui-Nature-2023-ExFig10.xlsx",
  rowNames = F, colNames = T, check.names = F
)
cytokines <- toupper(cytokines.df$Receptor.Gene)
cytokines[!cytokines %in% rownames(sro.subset)]
cytokines <- cytokines[cytokines %in% rownames(sro.subset)]
length(cytokines)
sro.subset <- readRDS("results/SRO_subset.rds")
imp.expr <- readRDS("results/imputed-expr.rds")
dir.create("plots/cytokines")


make.chunks <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
chunks <- make.chunks(cytokines, 4)

## vln plot ####
lapply(1:length(chunks), FUN = function(i){
  gene.list <- chunks[[i]]
  n.col <- 3
  n.row <- ceiling(length(gene.list)/n.col)
  VlnPlot(object = sro.subset, features = gene.list, group.by = 'Annotation', assay = "RNA", cols = pal$Annotation, ncol = n.col, pt.size = 0.1) %>%
    ggsave(filename = paste0("plots/cytokines/violin-cytokines-", i, ".pdf"), 
           width = 6*n.col, height = 4*n.row)
})

## Dot plots ####
geneList <- rev(cytokines) ; fn <- "cytokines"
group.name <- "Annotation"
ggsave(
  filename = paste0("plots/cytokines/dot-", fn, ".pdf"), width = 8, height = 20,
  plot = DotPlot(sro.subset, assay = "RNA", group.by = group.name,
                 features = geneList, cols = c("blue", "red"), dot.scale = 7) +
    scale_x_discrete(breaks = geneList, labels = geneList) +
    theme(axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 11)) + labs(x = "", y = group.name) + coord_flip()
)

## umap overlay ####
lapply(1:length(chunks), FUN = function(i){
  gene.list <- chunks[[i]]
  n.col <- 5
  n.row <- ceiling(length(gene.list)/n.col)
  
  umap.unimputed <- lapply(gene.list, function(gene){
    p <- plot.continuous.value(
      sro.subset, vis = sro.subset@reductions$umap@cell.embeddings,
      idx = rownames(sro.subset@meta.data), point.size = 0.5,
      val = sro.subset@assays$RNA[gene, ], val.name = 'expression'
    ) + ggtitle(gene)
    return(p)
  })
  umap.imputed <- lapply(gene.list, function(gene){
    p <- plot.continuous.value(
      sro.subset, vis = sro.subset@reductions$umap@cell.embeddings,
      idx = rownames(sro.subset@meta.data), point.size = 0.5,
      val = imp.expr[gene, ], val.name = 'imputed\nexpression'
    ) + ggtitle(gene)
    return(p)
  })
  plot_grid(plotlist = umap.unimputed, ncol = n.col) %>%
    ggsave(filename = paste0("plots/cytokines/umap-", fn, "-unimputed-", i, ".pdf"),
           width = n.col*6, height = n.row*4)
  plot_grid(plotlist = umap.imputed, ncol = n.col) %>%
    ggsave(filename = paste0("plots/cytokines/umap-", fn, "-imputed-", i, ".pdf"),
           width = n.col*6, height = n.row*4)
})

## heatmap ####

sro.subset@meta.data$Annotation <- droplevels(sro.subset@meta.data$Annotation)
pdf(width = 16, height = 14, file = paste0("plots/cytokines/heatmap-cytokines.pdf"))
DE.heatmap(sro.subset, expr = imp.expr,
           markers = cytokines,
           split = 'Annotation',
           cluster_rows = T,
           split.columns = c("sample", "Annotation"), cluster_column_slices = F,
           ca.col = list(sample = pal$sample, "Annotation" = pal$Annotation))
dev.off()

ceiling(length(cytokines) * 0.12)+3


# gene umap ####
"RALDH2" %in% rownames(sro.subset)
umap.unimputed <- lapply(c("RALDH2"), function(gene){
  p <- plot.continuous.value(
    sro.subset, vis = sro.subset@reductions$umap@cell.embeddings,
    idx = rownames(sro.subset@meta.data), point.size = 0.5,
    val = sro.subset@assays$RNA$data[gene, ], val.name = 'imputed\nexpression'
  ) + ggtitle(gene)
  return(p)
})
umap.imputed <- lapply(c("RALDH2"), function(gene){
  p <- plot.continuous.value(
    sro.subset, vis = sro.subset@reductions$umap@cell.embeddings,
    idx = rownames(sro.subset@meta.data), point.size = 0.5,
    val = imp.expr[gene, ], val.name = 'imputed\nexpression'
  ) + ggtitle(gene)
  return(p)
})
plot_grid(plotlist = umap.unimputed, ncol = n.col) %>%
  ggsave(filename = paste0("plots/cytokines/umap-", fn, "-unimputed-", i, ".pdf"),
         width = n.col*6, height = n.row*4)
plot_grid(plotlist = umap.imputed, ncol = n.col) %>%
  ggsave(filename = paste0("plots/gene-expr/umap-", fn, "-imputed-", i, ".pdf"),
         width = n.col*6, height = n.row*4)

# ############################################################################ #
# Violin plots of Glutamate family of receptors ####
# ############################################################################ #
RenameGenesSeurat <- function(obj, newnames) { # Replace gene names in different slots of a Seurat object. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA
  
  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]] <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]] <- newnames
    if (length(RNA@scale.data)) rownames( RNA@scale.data) <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}
sro <- readRDS("results/SRO.rds")
sro.subset <- subset(sro, Annotation == "TC I")

cell.count <- rowSums(sro.subset@assays$RNA@counts > 0)
dir.create("plots/gene-expr/glutamate")

ampa.receptors <- sort(rownames(sro.subset)[grepl("GRIA", rownames(sro.subset))])
# ampa.receptors <- names(cell.count[ampa.receptors][cell.count[ampa.receptors] > 3])
nmda.receptors <- sort(rownames(sro.subset)[grepl("GRIN", rownames(sro.subset))])
nmda.receptors <- names(cell.count[nmda.receptors][cell.count[nmda.receptors] > 3])

group.name <- "Annotation"
gene.list <- ampa.receptors
fn <- paste0("violin-", group.name, "-AMPA-all4")
n.col <- 4
n.row <- ceiling(length(gene.list)/n.col)

VlnPlot(object = sro.subset, features = gene.list, group.by = group.name, assay = "RNA", slot = 'data', 
        cols = pal$Annotation, ncol = n.col, pt.size = 0.1) %>%
  ggsave(filename = paste0("plots/gene-expr/glutamate/", fn, ".pdf"), 
         width = 4*n.col, height = 4*n.row)

## new violin plot format ####
# https://github.com/ycl6/StackedVlnPlot
prep.violin.df <- function(sro, features, group.by, idents = NULL){
  # Subset data.frame
  sro.df <- FetchData(sro, features, slot = "data")
  # Add cell ID and identity classes
  sro.df$Cell <- rownames(sro.df)
  sro.df$Idents <- sro@meta.data[, group.by]
  if (!is.null(idents)){
    sro.df <- sro.df[sro.df$Idents %in% idents,]
  }
  # Use melt to change data.frame format
  sro.df <- reshape2::melt(sro.df, id.vars = c("Cell","Idents"), 
                           measure.vars = features, variable.name = "Feat", value.name = "Expr")
  return(sro.df)
}
plot.violin <- function(df, title = ""){
  # Features on x-axis
  p <- ggplot(df, aes(factor(Feat), Expr, fill = Feat)) +
    geom_violin(scale = "width", adjust = 1, trim = TRUE) +
    scale_y_continuous(position="right") +
    # scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
    #   c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
    facet_grid(rows = vars(Idents), scales = "free", switch = "y") +
    theme_cowplot(font_size = 12) +
    theme(legend.position = "none", panel.spacing = unit(0, "lines"),
          plot.title = element_text(hjust = 0.5),
          panel.background = element_rect(fill = NA, color = "black"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          strip.text.y.left = element_text(angle = 0),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggtitle(title) + xlab("Feature") + ylab("Expression Level")
  return(p)
}

vln.df <- prep.violin.df(sro.subset, gene.list, "Annotation", c("TC I"))
table(vln.df$Feat)
gria2 <- vln.df[vln.df$Feat == "GRIA1",]
gria2$Feat <- "GRIA2"
gria2$Expr <- 0
vln.df <- rbind(vln.df, gria2)
vln.df$Feat <- factor(vln.df$Feat, levels = c("GRIA1", "GRIA2", "GRIA3", "GRIA4"))
c <- plot.violin(vln.df)
ggsave(filename = paste0("plots/gene-expr/glutamate/", fn, "-format2.pdf"), c,
       width = 8, height = 4)
# ggsave(filename = paste0("plots/gene-expr/glutamate/", fn, "-format2.pdf"), c,
#        width = 8, height = 4)

# THRA/B ####
sro.subset <- readRDS(file = "results/SRO_subset.rds")
sro.subset$gut <- ifelse(test = sro.subset$sample %in% c("para-aortic", "mediastinal", "salivary gland", "skin"),
                         yes = "non-gut",
                         no = "gut")
fn <- "THRA_THRB_PRDM16"
Idents(sro.subset) <- sro.subset$gut
VlnPlot(object = sro.subset, features = c("THRA", "THRB", "PRDM16"), group.by = "gut", assay = "RNA", slot = 'data', 
        cols = pal$gut, ncol = 1, pt.size = 0.1) %>%
  ggsave(filename = paste0("plots/gene-expr/glutamate/violin-", fn, ".pdf"), 
         width = 4*1, height = 4*3)

VlnPlot(object = sro.subset, features = c("THRA", "THRB", "PRDM16"), group.by = "Annotation", split.by = "gut", assay = "RNA", slot = 'data', 
        cols = pal$gut, ncol = 1, pt.size = 0.1) %>%
  ggsave(filename = paste0("plots/gene-expr/glutamate/violin-splitBy-gut-", fn, ".pdf"), 
         width = 8, height = 4*3)


markers <- data.frame(gene = c("THRA", "THRB", "PRDM16"), gene_name = c("THRA", "THRB", "PRDM16"))
pdf(width = 20, height = 8, file = paste0("plots/gene-expr/glutamate/heatmap-", fn, ".pdf"))
print(DE.heatmap(sro.subset, expr = imp.expr, fn = fn, markers = markers,
                 split = 'gut', column_title_rot = 45, row_title_rot = 0,
                 split.columns = c("sample", "Annotation"), cluster_column_slices = F,
                 ca.col = list(sample = pal$sample, "Annotation" = pal$Annotation)))
dev.off()



p1 <- VlnPlot(sro.subset, "THRA", cols = pal$gut, pt.size = 0.1, split.by = "gut", group.by = "Annotation") +theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))  + geom_boxplot(position=position_dodge(1), alpha = 0)
p2 <- VlnPlot(sro.subset, "THRB", cols = pal$gut, pt.size = 0.1, split.by = "gut", group.by = "Annotation")+theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + geom_boxplot(position=position_dodge(1), alpha = 0)
p3 <- VlnPlot(sro.subset, "PRDM16", cols = pal$gut, pt.size = 0.1, split.by = "gut", group.by = "Annotation")+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))  + geom_boxplot(position=position_dodge(1), alpha = 0)
pdf(width = 10, height = 12, file = paste0("plots/gene-expr/glutamate/violin-splitBy-gut-groupBy-Annot-", fn, ".pdf"))
plot_grid(p1, p2, p3, ncol = 1)
dev.off()
sro.subset@meta.data$Annotation <- droplevels(sro.subset@meta.data$Annotation)
groups <- levels(sro.subset@meta.data$Annotation)
idx <- rownames(sro.subset@meta.data[sro.subset@meta.data$Annotation == groups[4],])
res <- FoldChange(sro.subset[,idx], ident.1 = "gut", ident.2 = "non-gut", group.by = "gut", 
           subset.ident = NULL, assay = "RNA", slot = "data")
res[c("THRA", "THRB", "PRDM16"),]

# dot plot ####
sro.subset <- readRDS("results/SRO_subset.rds")
# DR3 (Tnfrsf25), OX40L (Tnfsf4), and Bhlhe40
geneList <- c("TNFRSF25", "TNFSF4", "BHLHE40"); 
fn <- "TNFRSF25_TNFSF4_BHLHE40"
group.name <- "Annotation"
ggsave(
  filename = paste0("plots/gene-expr/dot-", fn, ".pdf"), width = 8, height = 5,
  plot = DotPlot(sro.subset, assay = "RNA", group.by = group.name,
                 features = geneList, cols = c("blue", "red"), dot.scale = 6) +
    scale_x_discrete(breaks = geneList, labels = geneList) +
    theme(axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 11)) + labs(x = "", y = group.name) + coord_flip()
)

# GEO submission ####
# Analysis results from Seurat including clusters, annotations, and UMAP

sro.subset <- readRDS("results/SRO_subset.rds")
DimPlot(sro.subset, group.by = "Annotation", label = T, label.box = T)+DimPlot(sro.subset, group.by = "Clusters", label = T, label.box = T)
rna.md <- sro.subset@meta.data[c("sample", "hash_id", "nCount_RNA", "nFeature_RNA", "percent.MT", "S.Score", "G2M.Score", "Phase", "Clusters2", "Annotation")]
colnames(rna.md)[9] <- 'Clusters'

rna.md$UMAP_1 <- sro.subset@reductions$umap@cell.embeddings[, 1]
rna.md$UMAP_2 <- sro.subset@reductions$umap@cell.embeddings[, 2]

write.table(rna.md, "results/results.tsv", sep = "\t", quote = F)
table(sro.subset$Clusters2)
