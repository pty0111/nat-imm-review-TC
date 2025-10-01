# TC all LN ####
suppressPackageStartupMessages({
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
})

set.seed(1)
options(future.globals.maxSize = Inf)

# pal ####
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

# ############################################################################ #
# Redo UMAP ####
sro.subset <- subset(sro, Annotation %in% c("Ki67+ TC", "early/transitional TC", "TC I", "TC II", 'TC III', "TC IV")) %>%
  RunUMAP(dims = 1:30, n.neighbors = 30, metric = "cosine", min.dist = 0.4, spread = 0.7)

## re-order clusters in sro.subset
sro.subset <- readRDS("results/SRO_subset.rds")
sro.subset$Clusters2 <- sro@meta.data[colnames(sro.subset), "Clusters2"]
sro.subset$Clusters2 <- droplevels(sro.subset$Clusters2)

saveRDS(sro.subset, file = "results/SRO_subset.rds")
write.csv(sro.subset@reductions$umap@cell.embeddings, file = "results/UMAP-subset.csv")
