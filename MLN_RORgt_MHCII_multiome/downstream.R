# -------------------------------------------------------------------------
setwd("~/0-workspace/CCR7_DC/MLN_RORgt_MHCII_multiome/")

library(ArchR) # v1.0.1
library(Seurat) # v4.0.4
library(Signac) # v1.4.0.9006
library(Rphenograph) # v0.99.1
library(TFBSTools)
library(motifmatchr)
library(parallel)
library(gprofiler2)
library(Rmagic)
library(Matrix)
library(future)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ComplexHeatmap)
library(circlize)

set.seed(1)
options(future.globals.maxSize = Inf)


# -------------------------------------------------------------------------
setwd("ArchR/")
addArchRThreads(threads = 32)
addArchRGenome("mm10")

atac.sample <- "3228_Sample_CB-1288_MLN_RORgt_MHCII_multiome_ATAC_IGO_12437_C_19"
arrow.files <- createArrowFiles(
  inputFiles = "../cellranger-arc/outs/atac_fragments.tsv.gz", sampleNames = atac.sample,
  minFrags = 10 ^ 3.5, maxFrags = 10 ^ 4.5, excludeChr = c("chrM", "chrY"), addGeneScoreMat = F
)
doubScores <- addDoubletScores(arrow.files)
archr.obj <- ArchRProject(arrow.files) %>% filterDoublets() # 1599 doublets removed!
# numberOfCells: 11050
# medianTSS: 21.5545
# medianFrags: 13676.5

saveArchRProject(archr.obj)
write.csv(archr.obj@cellColData, file = "ArchROutput/cellColData_before_filtering_by_scRNA-seq.csv")


# -------------------------------------------------------------------------
setwd("Seurat/")

run.PhenoGraph <- function(sro, npcs, k){
  sro@misc[["PhenoGraph"]] <- Rphenograph(sro@reductions$pca@cell.embeddings[, 1:npcs], k = k)
  adj <- as_adjacency_matrix(sro@misc$PhenoGraph[[1]], sparse = F)
  rownames(adj) <- colnames(sro); colnames(adj) <- colnames(sro)
  sro@graphs[["PhenoGraph"]] <- as.Graph(adj)
  sro$Clusters <- sro@misc$PhenoGraph[[2]]$membership
  sro$Clusters <- factor(sro$Clusters, levels = min(sro$Clusters):max(sro$Clusters))
  return(sro)
}
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
    theme_classic() +
    labs(x = ident, y = "", title = ifelse(Log10, paste0("log10(", name, ")"), name))
  if(!is.null(pal)) gp <- gp + scale_color_manual(values = pal) + scale_fill_manual(values = pal)
  return(gp)
}
plot.all.QC <- function(sro, ident, thr = NULL, col = NULL){
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
  return(pl)
}

counts <- readMM("../cellranger-arc/outs/filtered_feature_bc_matrix/matrix.mtx.gz")
barcodes <- readLines("../cellranger-arc/outs/filtered_feature_bc_matrix/barcodes.tsv.gz")
feature.info <- read.table("../cellranger-arc/outs/filtered_feature_bc_matrix/features.tsv.gz", sep = "\t", header = F,
                           stringsAsFactors = F, col.names = c("ID", "symbol", "type", "chr", "start", "end"))[, 1:3]
rownames(counts) <- feature.info$ID; colnames(counts) <- barcodes

gi <- feature.info$type == "Gene Expression"
gene.info <- feature.info[gi, ]
counts <- counts[gi, ]

seurat.obj <- CreateSeuratObject(counts = counts, project = "MLN_RORgt_MHCII_multiome")
seurat.obj$percent.MT <- PercentageFeatureSet(seurat.obj, features = subset(gene.info, grepl("^MT-", toupper(symbol)))$ID)
seurat.obj@assays$RNA <- AddMetaData(seurat.obj@assays$RNA, metadata = gene.info$symbol, col.name = "Symbol")

cells.keep <- read.csv("../ArchR/ArchROutput/cellColData_before_filtering_by_scRNA-seq.csv", row.names = 1, header = T) %>%
  rownames() %>% gsub(pattern = "^.+#", replacement = "")
thr <- data.frame(nc.min = 1000, nc.max = 50000, nf.min = 500, nf.max = 6000, mp.max = 15)
cell.discard <- seurat.obj$nCount_RNA < thr$nc.min | seurat.obj$nCount_RNA > thr$nc.max |
  seurat.obj$nFeature_RNA < thr$nf.min | seurat.obj$nFeature_RNA > thr$nf.max |
  seurat.obj$percent.MT > thr$mp.max | !(colnames(seurat.obj) %in% cells.keep)
# 5720 discarded in total (10468 left).
# 778 discarded regardless of scATAC-seq filtering.
# 503 discarded in addition to scATAC-seq filtering.
pdf("plots/QC/sample-QC.pdf", width = 20, height = 20)
plot.all.QC(seurat.obj, ident = "Sample", thr = thr)
dev.off()

seurat.obj <- seurat.obj[, !cell.discard]
cell.count <- rowSums(seurat.obj@assays$RNA@counts > 0)
seurat.obj <- NormalizeData(seurat.obj[cell.count > 1, ]) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 5000)
seurat.obj <- ScaleData(seurat.obj, features = rownames(seurat.obj)) %>% RunPCA(features = rownames(seurat.obj)) %>%
  run.PhenoGraph(npcs = 30, k = 30) %>% RunUMAP(graph = "PhenoGraph")

seurat.obj.sub <- subset(seurat.obj, Clusters %in% 1:16) %>%
  run.PhenoGraph(npcs = 30, k = 30) %>% RunUMAP(graph = "PhenoGraph", spread = 0.7, min.dist = 0.4)
write.csv(seurat.obj.sub@reductions$umap@cell.embeddings, file = "results/UMAP_C1-16.csv")

s.genes <- rownames(subset(seurat.obj@assays$RNA@meta.features, toupper(Symbol) %in% cc.genes.updated.2019$s.genes))
g2m.genes <- rownames(subset(seurat.obj@assays$RNA@meta.features, toupper(Symbol) %in% cc.genes.updated.2019$g2m.genes))
seurat.obj <- CellCycleScoring(seurat.obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
seurat.obj$Phase <- gsub("G2M", "G2/M", seurat.obj$Phase)

rna.annot <- c(
  `Proliferating TC` = 19, # Early TCs
  `TC ` = c(17, # similar to Aire+ mTECs
            2 , # similar to one of the post Aire Amit populations in mTEC III that expressed PIGR/Claudin7
            18, # Aire hi but shares DC genes similar to TC IV
            11), # RORgt+ DC
  `ILC3p` = 1,
  `Proliferating NCR+ ILC3` = 15,
  `NCR+ ILC3` = 21,
  `LTi Variation ` = c(3, 6, 8, 9, 10, 13, 14, 16),
  `Neuronal/Microglial Cell` = 12, # Ptprc-
  `pDC` = 5 , # Low RORc
  `DC` = 7 , # Low RORc
  `Macrophage ` = c(4, # Low RORc
                    20) # Low RORc
)
seurat.obj$Clusters.annot <- factor(names(rna.annot)[match(seurat.obj$Clusters, rna.annot)], levels = names(rna.annot))
rename <- order(as.integer(rna.annot))
seurat.obj$Clusters <- factor(rename[seurat.obj$Clusters], levels = 1:21)
Idents(seurat.obj) <- seurat.obj$Clusters

seurat.obj$Clusters2 <- ifelse(seurat.obj$Clusters %in% 9:16, "9-16", seurat.obj$Clusters) %>%
  factor(levels = c(1:8, "9-16", 17:21))
Idents(seurat.obj) <- seurat.obj$Clusters2
pdf("~/Desktop/cluster-QC.pdf", width = 20, height = 20)
plot.all.QC(
  sro = seurat.obj, ident = "Clusters2",
  col = c(rna.pal2[1:8], `9-16` = as.character(rna.pal[11]), rna.pal[17:21])
)
dev.off()

saveRDS(seurat.obj, file = "results/SRO.rds")
write.csv(seurat.obj@meta.data, file = "results/meta-data.csv")
write.csv(seurat.obj@reductions$pca@cell.embeddings, file = "results/PCA.csv")
write.csv(seurat.obj@reductions$umap@cell.embeddings, file = "results/UMAP.csv")

seurat.obj.imp <- magic(seurat.obj)
saveRDS(seurat.obj.imp@assays$MAGIC_RNA@data, file = "results/imputed-expr.rds")

seurat.obj <- readRDS(file = "Seurat/results/SRO.rds")
saveRDS(as.matrix(seurat.obj@assays$RNA@data), file = "Seurat/results/unimputed-expr.rds")
# -------------------------------------------------------------------------
setwd("Seurat/")

MAST.DE <- function(fn, sro, ...){
  m <- suppressWarnings(FindMarkers(object = sro, test.use = "MAST", logfc.threshold = 0, ...))
  m$gene_name <- sro@assays$RNA@meta.features[rownames(m), "Symbol"]
  write.csv(m, file = paste0("results/markers/", fn, ".csv"), row.names = T, quote = F)
}
MAST.DE.multiple <- function(fn, sro, idents = sro@active.ident, ...){
  Idents(sro) <- idents
  m <- FindAllMarkers(object = sro, test.use = "MAST", logfc.threshold = 0, ...)
  m$gene_name <- sro@assays$RNA@meta.features[m$gene, "Symbol"]
  write.csv(m, file = paste0("results/markers/", fn, ".csv"), row.names = F, quote = F)
  return(m)
}
select.markers <- function(fn, pairwise = F, fc.thr = 1.5, apv.thr = 0.01, n = Inf){
  markers <- read.csv(file = paste0("results/markers/", fn, ".csv"), header = T, stringsAsFactors = F) %>%
    subset(!grepl("(^MT-)|(^RP)", toupper(gene_name)))
  if(pairwise){
    colnames(markers)[1] <- "gene"
    markers <- subset(markers, abs(avg_log2FC) > log2(fc.thr) & p_val_adj < apv.thr)
    markers$cluster <- ifelse(markers$avg_log2FC > 0, "up", "down")
    if(!is.infinite(n)) markers <- markers %>% group_by(cluster) %>% top_n(n, abs(avg_log2FC))
    markers <- markers[order(markers$avg_log2FC), ]
  } else{
    markers <- subset(markers, avg_log2FC > log2(fc.thr) & p_val_adj < apv.thr)
    if(any(duplicated(markers$gene))) markers <- markers %>% group_by(gene) %>% top_n(1, avg_log2FC)
    if(!is.infinite(n)) markers <- markers %>% group_by(cluster) %>% top_n(n, abs(avg_log2FC))
    markers <- markers[order(markers$cluster, markers$avg_log2FC), ]
  }
  return(markers)
}
DE.heatmap <- function(sro = NULL, expr = sro@assays$RNA@data, cluster_columns = T,
                       split = md$Clusters, annot = md$Clusters.annot, annot.pal = rna.pal[1:16],
                       md = sro@meta.data, cells = rep(T, nrow(md)), ...){
  markers <- select.markers(...)
  ca <- columnAnnotation(col = list(Cluster = annot.pal), Cluster = annot[cells])
  col.ramp <- colorRamp2(c(-10, -5, seq(-3, 3, 1), 5, 10),
                         c("#173c68", "#173c68", "#1c5785", "#166db0", "#a6d1e3", "#ffffff",
                           "#f09173", "#d46052", "#b32339", "#690822", "#690822"))
  Heatmap(t(scale(t(expr[markers$gene, cells]))), name = "scaled\nimputed\nexpression", col = col.ramp,
          top_annotation = ca, show_column_names = F, column_split = split[cells],
          cluster_rows = F, row_names_side = "right", row_labels = markers$gene_name,
          cluster_columns = cluster_columns, cluster_column_slices = F,
          row_names_gp = gpar(fontface = "italic", fontsize = 5), use_raster = T)
}

plan(strategy = "multicore", workers = 8)
MAST.DE.multiple(fn = "C1-21", sro = seurat.obj)
MAST.DE.multiple(fn = "C1-16", sro = subset(seurat.obj, Clusters %in% 1:16))
MAST.DE.multiple(fn = "C1-5",  sro = subset(seurat.obj, Clusters %in% 1:5))
MAST.DE.multiple(fn = "C2-5",  sro = subset(seurat.obj, Clusters %in% 2:5))
MAST.DE.multiple(fn = "C6-8",  sro = subset(seurat.obj, Clusters %in% 6:8))
MAST.DE.multiple(fn = "C6-16", sro = subset(seurat.obj, Clusters %in% 6:16))
MAST.DE.multiple(fn = "C9-16", sro = subset(seurat.obj, Clusters %in% 9:16))
MAST.DE(fn = "C2_vs_C4", sro = seurat.obj, ident.1 = 2, ident.2 = 4)
MAST.DE(fn = "C5_vs_C1-4", sro = seurat.obj, ident.1 = 5, ident.2 = 1:4)
MAST.DE(fn = "C1-5_vs_C8-16", sro = seurat.obj, ident.1 = 1:5, ident.2 = 8:16)

Idents(seurat.obj) <- ifelse(seurat.obj$Clusters %in% 2:5, "TC 1-4", as.character(seurat.obj$Clusters.annot)) %>%
  gsub(pattern = "LTi Variation \\d+", replacement = "LTi")
MAST.DE.multiple(fn = "groups", sro = subset(seurat.obj, Clusters %in% c(2:5, 8:16)))

Idents(seurat.obj) <- gsub("LTi Variation \\d+", "LTi", as.character(seurat.obj$Clusters.annot))
MAST.DE.multiple(fn = "groups_individualTCs", sro = subset(seurat.obj, Clusters %in% c(2:5, 8:16)))
MAST.DE.multiple(fn = "C1-16_groupedLTis", sro = subset(seurat.obj, Clusters %in% 1:16))

imp.expr <- readRDS("results/imputed-expr.rds")
md <- read.csv("results/meta-data.csv", header = T, row.names = 1)
md$Clusters <- factor(md$Clusters, levels = c(1:8, 10, 12, 9, 13, 11, 14, 15, 16))
md$Clusters.annot <- factor(md$Clusters.annot, levels = names(rna.annot))

fn <- "groups"; n <- 50
ci <- which(md$Clusters %in% c(2:5, 8:16)); ci <- ci[order(md$Clusters[ci])]
pdf(width = 20, height = 20, file = paste0("plots/markers/", fn, "_top", n, ".pdf"))
DE.heatmap(fn = fn, n = n, md = md, expr = imp.expr, cells = ci, cluster_columns = F)
dev.off()

fn <- "groups_individualTCs"; n <- 50
ci <- which(md$Clusters %in% c(2:5, 8:16)); ci <- ci[order(md$Clusters[ci])]
pdf(width = 30, height = 25, file = paste0("plots/markers/", fn, "_top", n, ".pdf"))
DE.heatmap(fn = fn, n = n, md = md, expr = imp.expr, cells = ci, cluster_columns = F)
dev.off()

fn <- "C1-21"; n <- 20
pdf(width = 35, height = 30, file = paste0("plots/markers/", fn, "_top", n, ".pdf"))
DE.heatmap(fn = fn, n = n, md = md, expr = imp.expr, annot.pal = rna.pal)
dev.off()

fn <- "C1-16"; n <- 30
pdf(width = 30, height = 30, file = paste0("plots/markers/", fn, "_top", n, ".pdf"))
DE.heatmap(fn = fn, n = n, md = md, expr = imp.expr, cells = md$Clusters %in% 1:16)
dev.off()

fn <- "C1-16_groupedLTis"; n <- 50
pdf(width = 30, height = 40, file = paste0("plots/markers/", fn, "_top", n, ".pdf"))
DE.heatmap(fn = fn, n = n, md = md, expr = imp.expr, cells = md$Clusters %in% 1:16,
           split = ifelse(md$Clusters %in% 9:16, "9-16", md$Clusters),
           annot = gsub("LTi Variation \\d+", "LTi", as.character(md$Clusters.annot)),
           annot.pal = c(rna.pal[1:8], LTi = as.character(rna.pal[11])))
dev.off()

fn <- "C1-5"; n <- 50
pdf(width = 15, height = 20, file = paste0("plots/markers/", fn, "_top", n, ".pdf"))
DE.heatmap(fn = fn, n = n, md = md, expr = imp.expr, cells = md$Clusters %in% 1:5)
dev.off()

fn <- "C2-5"; n <- 130
pdf(width = 15, height = 35, file = paste0("plots/markers/", fn, "_top", n, ".pdf"))
DE.heatmap(fn = fn, n = n, md = md, expr = imp.expr, cells = md$Clusters %in% 2:5)
dev.off()

fn <- "C6-8"; n <- 50
pdf(width = 10, height = 15, file = paste0("plots/markers/", fn, "_top", n, ".pdf"))
DE.heatmap(fn = fn, n = n, md = md, expr = imp.expr, cells = md$Clusters %in% 6:8)
dev.off()

fn <- "C6-16"; n <- 50
pdf(width = 20, height = 20, file = paste0("plots/markers/", fn, "_top", n, ".pdf"))
DE.heatmap(fn = fn, n = n, md = md, expr = imp.expr, cells = md$Clusters %in% 6:16)
dev.off()

fn <- "C9-16"; n <- 50
pdf(width = 20, height = 10, file = paste0("plots/markers/", fn, "_top", n, ".pdf"))
DE.heatmap(fn = fn, n = n, md = md, expr = imp.expr, cells = md$Clusters %in% 9:16)
dev.off()

fn <- "C2_vs_C4"; n <- 100
pdf(width = 10, height = 20, file = paste0("plots/markers/", fn, "_top", n, ".pdf"))
DE.heatmap(fn = fn, pairwise = T, n = n, md = md, expr = imp.expr, cells = md$Clusters %in% c(2, 4))
dev.off()

fn <- "C5_vs_C1-4"; n <- 50
pdf(width = 15, height = 10, file = paste0("plots/markers/", fn, "_top", n, ".pdf"))
DE.heatmap(fn = fn, pairwise = T, n = n, md = md, expr = imp.expr, cells = md$Clusters %in% 1:5)
dev.off()

fn <- "C1-5_vs_C8-16"; n <- 50
pdf(width = 20, height = 15, file = paste0("plots/markers/", fn, "_top", n, ".pdf"))
DE.heatmap(fn = fn, pairwise = T, n = n, md = md, expr = imp.expr, cells = md$Clusters %in% c(1:5, 8:16))
dev.off()


# -------------------------------------------------------------------------
setwd("ArchR/")

rna.md <- read.csv("../Seurat/results/meta-data.csv", row.names = 1, header = T) %>% subset(Clusters %in% 1:16)
cells.keep <- paste(atac.sample, rownames(rna.md), sep = "#") %>% intersect(archr.obj$cellNames)
archr.obj <- subsetArchRProject(archr.obj, cells = cells.keep, outputDirectory = "ArchROutput_SeuratC1-16")
# numberOfCells: 10145
# medianTSS: 21.648
# medianFrags: 13885

archr.obj <- addIterativeLSI(archr.obj, iterations = 5, varFeatures = 100000, scaleDims = T) %>%
  addClusters(k.param = 30, resolution = 1.2) %>% addUMAP(name = "UMAP", nNeighbors = 30)
archr.obj$Clusters.scRNAseq <- paste0("C", rna.md[gsub("^.+#", "", archr.obj$cellNames), "Clusters"])

plotPDF(
  plotEmbedding(archr.obj, name = "nFrags", embedding = "UMAP"),
  plotEmbedding(archr.obj, name = "PromoterRatio", embedding = "UMAP"),
  plotEmbedding(archr.obj, name = "DoubletScore", embedding = "UMAP"),
  name = "after_filtering_by_scRNA-seq.pdf", ArchRProj = archr.obj, width = 15, height = 12
)

saveArchRProject(archr.obj)
write.csv(archr.obj@cellColData, file = "ArchROutput_SeuratC1-16/cellColData.csv")
write.csv(getEmbedding(archr.obj), file = "ArchROutput_SeuratC1-16/Embeddings/UMAP.csv")
write.csv(getReducedDims(archr.obj), file = "ArchROutput_SeuratC1-16/IterativeLSI/LSI.csv")


# -------------------------------------------------------------------------
setwd("ArchR/")

groups <- c("C1", rep("C2-4", 3), "C5,7", "C6,8-13", "C5,7", rep("C6,8-13", 6))
archr.obj$Groups <- groups[as.integer(gsub("C", "", archr.obj$Clusters))]

arrow.file <- "3228_Sample_CB-1288_MLN_RORgt_MHCII_multiome_ATAC_IGO_12437_C_19.arrow"
group.nFrags <- sapply(unique(archr.obj$Groups), function(g){
  cat("Extracting fragments for", g, "\n")
  group.cells <- subset(archr.obj@cellColData, Groups == g) %>% rownames()
  group.frags <- getFragmentsFromArrow(ArrowFile = arrow.file, cellNames = group.cells)
  cbind(as.character(seqnames(group.frags)), start(group.frags), end(group.frags)) %>%
    write.table(file = paste0("ArchROutput_SeuratC1-16/GroupFragments/", g, ".bed"),
                row.names = F, col.names = F, quote = F, sep = "\t")
  return(length(group.frags))
})

peaks <- readRDS("ArchROutput_SeuratC1-16/PeakCalls/union-peaks.rds")
archr.obj <- addPeakSet(archr.obj, peakSet = peaks) %>% addPeakMatrix()
# 176952 peaks with maximum length 1158 and median length 206

saveArchRProject(archr.obj)
getMatrixFromProject(archr.obj, useMatrix = "PeakMatrix", binarize = F) %>%
  saveRDS("ArchROutput_SeuratC1-16/PeakCalls/counts.rds")

tc.cells <- archr.obj$cellNames[archr.obj$Clusters.scRNAseq %in% paste0("C", 1:5)]
tc.peaks <- readRDS(file = "ArchROutput_SeuratC1-16/PeakCalls/union-peaks-TCs.rds")
tc.archr.obj <- subsetArchRProject(
  ArchRProj = archr.obj, cells = tc.cells,
  outputDirectory = "ArchROutput_TCs/"
) %>% addPeakSet(peakSet = tc.peaks, force = T) %>% addPeakMatrix(force = T)
saveArchRProject(tc.archr.obj)


# -------------------------------------------------------------------------
setwd("ArchR/")

archr.obj$Groups <- archr.obj$Clusters %>%
  ifelse(test = archr.obj$Clusters == "C7", yes = "C6") %>%
  ifelse(test = archr.obj$Clusters %in% paste0("C", c(6, 8:13)), yes = "C7-13")

archr.obj$Clusters2 <- ifelse(archr.obj$Clusters == "C7", "C6",
                              ifelse(archr.obj$Clusters == "C6", "C7",
                                     archr.obj$Clusters))
cl.pal2 <- rna.pal2[c(5, 4, 2, 3, 6, 8, 11:17)]
names(cl.pal2) <- paste0("C", 1:length(cl.pal2))
useClusters2 <- c("C3", "C4", "C2", "C1", paste0("C", 5:13))

group.pal <- rna.pal2[c(5, 4, 2, 3, 6, 8, 11)]
names(group.pal) <- paste0("C", c(1:6, "7-13"))
useGroups <- c("C3", "C4", "C2", "C1", "C5", "C6", "C7-13")
group.pal <- group.pal[useGroups]

group.bw <- getGroupBW(
  archr.obj, groupBy = "Groups", normMethod = "ReadsInTSS",
  tileSize = 100, maxCells = 10000, ceiling = 4
)

Ciita.range <- readRDS("ArchROutput_SeuratC1-16/Plots/tracks/Ciita.rds")
plotPDF(plotBrowserTrack(archr.obj, region = Ciita.range,
                         groupBy = "Groups", pal = group.pal, useGroups = useGroups),
        tileSize = 100, ylim = c(0, 1), sizes = c(10, 1, 4),
        ArchRProj = archr.obj, name = "Tracks-Ciita",
        width = 15, height = 15, baseSize = 15, facetbaseSize = 15)
plotPDF(plotBrowserTrack(archr.obj, region = Ciita.range,
                         groupBy = "Clusters2", pal = cl.pal2, useGroups = useClusters2),
        tileSize = 100, ylim = c(0, 1), sizes = c(10, 1, 4),
        ArchRProj = archr.obj, name = "Tracks-Ciita",
        width = 15, height = 20, baseSize = 15, facetbaseSize = 15)

plotPDF(plotBrowserTrack(archr.obj, geneSymbol = c("Rora"), upstream = 500000, downstream = 2000000,
                         groupBy = "Groups", pal = group.pal, useGroups = useGroups),
        tileSize = 100, ylim = c(0, 1), sizes = c(10, 1, 4),
        ArchRProj = archr.obj, name = "Tracks-Rora",
        width = 15, height = 15, baseSize = 15, facetbaseSize = 15)

plotPDF(plotBrowserTrack(archr.obj, geneSymbol = c("Cxcr6", "Itgax"), upstream = 50000, downstream = 50000,
                         groupBy = "Groups", pal = group.pal, useGroups = useGroups),
        tileSize = 100, ylim = c(0, 1), sizes = c(10, 1, 4),
        ArchRProj = archr.obj, name = "Tracks-Cxcr6,Itgax",
        width = 15, height = 15, baseSize = 15, facetbaseSize = 15)

plotPDF(plotBrowserTrack(archr.obj, geneSymbol = c("Rag1", "Rag2", "Clec9a"), upstream = 50000, downstream = 50000,
                         groupBy = "Groups", pal = group.pal, useGroups = useGroups),
        tileSize = 100, ylim = c(0, 1), sizes = c(10, 1, 4),
        ArchRProj = archr.obj, name = "Tracks-Rag1,Rag2,Clec9a",
        width = 15, height = 15, baseSize = 15, facetbaseSize = 15)

plotPDF(plotBrowserTrack(archr.obj, geneSymbol = c("Rorc", "Bcl11a"),
                         upstream = 200000, downstream = 200000,
                         groupBy = "Groups", pal = group.pal, useGroups = useGroups),
        tileSize = 100, ylim = c(0, 1), sizes = c(10, 1, 4),
        ArchRProj = archr.obj, name = "Tracks-Rorc,Bcl11a",
        width = 15, height = 15, baseSize = 15, facetbaseSize = 15)

plotPDF(plotBrowserTrack(archr.obj, geneSymbol = c("Chad", "Aire", "Irf8", "Spib"),
                         upstream = 50000, downstream = 50000,
                         groupBy = "Groups", pal = group.pal, useGroups = useGroups),
        tileSize = 100, ylim = c(0, 1), sizes = c(10, 1, 4),
        ArchRProj = archr.obj, name = "Tracks-Chad,Aire,Irf8,Spib",
        width = 15, height = 15, baseSize = 15, facetbaseSize = 15)

plotPDF(plotBrowserTrack(archr.obj, geneSymbol = c("Irf8"),
                         upstream = 100000, downstream = 100000,
                         groupBy = "Groups", pal = group.pal, useGroups = useGroups),
        tileSize = 100, ylim = c(0, 1), sizes = c(10, 1, 4),
        ArchRProj = archr.obj, name = "Irf8",
        width = 15, height = 15, baseSize = 15, facetbaseSize = 15)

plotPDF(plotBrowserTrack(archr.obj, geneSymbol = "Itgb8", upstream = 200000, downstream = 200000,
                         groupBy = "Groups", pal = group.pal, useGroups = useGroups),
        tileSize = 100, ylim = c(0, 1), sizes = c(10, 1, 4),
        ArchRProj = archr.obj, name = "Tracks-Itgb8",
        width = 15, height = 15, baseSize = 15, facetbaseSize = 15)

plotPDF(
  plotBrowserTrack(
    archr.obj, geneSymbol = "Rorc", upstream = 100000, downstream = 100000,
    groupBy = "Groups", pal = group.pal, useGroups = useGroups
  ), tileSize = 100, ylim = c(0, 1), sizes = c(10, 1, 4),
  ArchRProj = archr.obj, name = "Rorc",
  width = 15, height = 15, baseSize = 15, facetbaseSize = 15
)

plotPDF(
  plotBrowserTrack(
    archr.obj, geneSymbol = c("Adam23"), upstream = 500000, downstream = 500000,
    groupBy = "Groups", pal = group.pal, useGroups = useGroups
  ), tileSize = 100, ylim = c(0, 1), sizes = c(10, 1, 2, 2),
  ArchRProj = archr.obj, name = "Tracks-Adam23",
  width = 15, height = 15, baseSize = 15, facetbaseSize = 15
)

## TP tracks ####
plotPDF(plotBrowserTrack(archr.obj, geneSymbol = c("Pigr"),
                         upstream = 100000, downstream = 100000,
                         groupBy = "Groups", pal = group.pal, useGroups = useGroups,
                         tileSize = 100, ylim = c(0, 1), sizes = c(10, 1, 2, 2)),
        ArchRProj = archr.obj, name = "Tracks-Groups-Pigr",
        width = 8, height = 6, baseSize = 15, facetbaseSize = 15, addDOC = T)

marker.md <- data.frame(getGenes(ArchRProj = archr.obj, symbols = c("Cadm2")))
rownames(marker.md) <- marker.md$symbol
gene <- "Cadm2"
direction <- marker.md[gene, "strand"]
width <- marker.md[gene, "width"]
window.size <- 200000
if (direction == "+"){
  upstream <- window.size
  downstream <- width+window.size
} else{
  upstream <- width+window.size
  downstream <- window.size
}

plotPDF(plotBrowserTrack(archr.obj, geneSymbol = c("Cadm2"),
                         upstream = upstream, downstream = downstream,
                         groupBy = "Groups", pal = group.pal, useGroups = useGroups,
                         tileSize = 100, ylim = c(0, 1), sizes = c(10, 1, 2, 2)),
        ArchRProj = archr.obj, name = "Tracks-Groups-Cadm2",
        width = 20, height = 10, baseSize = 15, facetbaseSize = 15, addDOC = T)

marker.md <- data.frame(getGenes(ArchRProj = archr.obj, symbols = c("Rorc")))
rownames(marker.md) <- marker.md$symbol
gene <- "Rorc"
direction <- marker.md[gene, "strand"]
width <- marker.md[gene, "width"]
window.size <- 1000
if (direction == "+"){
  upstream <- window.size
  downstream <- width+window.size
} else{
  upstream <- width+window.size
  downstream <- window.size
}

width <- 500
starts <- c(marker.md$start+6000-width/2, marker.md$start+7000-width/2)
ends <- c(marker.md$start+6000+width/2, marker.md$start+7000+width/2)
enhancer.range <- GRanges(seqnames = c("chr3", "chr3"), 
       IRanges(starts, ends),
       strand = "+")
features.to.plot <- GRangesList(`peaks` = getPeakSet(archr.obj), `enhancer` = enhancer.range)

plotPDF(
  plotBrowserTrack(
    archr.obj, geneSymbol = c("Rorc"),
    upstream = upstream, downstream = downstream,
    features = features.to.plot,
    groupBy = "Groups", pal = group.pal, useGroups = useGroups,
    tileSize = 100, ylim = c(0, 1), sizes = c(10, 1, 2, 2)), 
  ArchRProj = archr.obj, name = "Tracks-Groups-Rorc",
  width = 12, height = 8, baseSize = 15, facetbaseSize = 15, addDOC = F)

## ####

g.list <- data.frame(readxl::read_xlsx("../Dotplot gene lists.xlsx"), check.names = F)
colnames(g.list) <- gsub("\\/", "-", colnames(g.list))
dummy <- lapply(colnames(g.list), function(list.name){
  gene.list <- na.omit(g.list[, list.name])
  cat(gene.list, "\n", sep = "\t")
  plotPDF(plotBrowserTrack(archr.obj, geneSymbol = gene.list, upstream = 200000, downstream = 200000,
                           groupBy = "Groups", pal = group.pal),
          tileSize = 100, ylim = c(0, 1), sizes = c(10, 1, 4),
          ArchRProj = archr.obj, name = "Tracks-Groups",
          width = 15, height = 15, baseSize = 15, facetbaseSize = 15)
})


# -------------------------------------------------------------------------
md.atac <- read.csv("ArchR/ArchROutput_SeuratC1-16/cellColData_before_filtering_by_scRNA-seq.csv",
                    row.names = 1, header = T, stringsAsFactors = F)
cl <- read.csv("ArchR/ArchROutput_SeuratC1-16/cellColData.csv", row.names = 1, header = T)[rownames(md.atac), "Clusters"]
cl <- gsub("X", "C6", gsub("C6", "C7", gsub("C7", "X", cl)))
md.atac <- data.frame(
  Barcode = gsub("^.+#", "", rownames(md.atac)),
  md.atac[, -c(1, 7)],
  Discarded = ifelse(is.na(cl), "Yes", "No"),
  Cluster = as.integer(gsub("C", "", cl)),
  read.csv("ArchR/ArchROutput_SeuratC1-16/Embeddings/UMAP.csv", header = T, row.names = 1,
           col.names = c("id", "UMAP_1", "UMAP_2"))[rownames(md.atac), ],
  read.csv("ArchR/ArchROutput_SeuratC1-16/IterativeLSI/LSI.csv", header = T, row.names = 1,
           col.names = c("id", paste0("LSI_", 1:30)))[rownames(md.atac), ]
)
write.table(md.atac, file = "MS/GEO/ATAC_ArchR-results.tsv", sep = "\t", quote = F, col.names = T, row.names = F)


# -------------------------------------------------------------------------
get.pwm <- function(motif.id, pseudocount = 0.008){
  filename <- paste0("Mus_musculus_2023_09_25_4-00_pm/pwms_all_motifs/", motif.id, ".txt")
  if(file.exists(filename)){
    pfm <- fread(filename)
    if(nrow(pfm) > 0){
      pfm <- pfm[, c('A', 'C', 'G', 'T')] %>% as.matrix() %>% t()
      pwm <- apply(pfm + pseudocount, 2, function(x){ x / sum(x) })
      pwm <- PWMatrix(profileMatrix = log(pwm * 4), ID = motif.id)
      return(pwm)
    }
  } else return(NULL)
}

loci <- "Tcf7"
tf <- "Runx3"
name <- paste0(loci, "_loci-", tf, "_motifs")

ext <- 100000
gene.loci <- ArchR::geneAnnoMm10$genes[which(ArchR::geneAnnoMm10$genes$symbol == loci)] %>%
  resize(width = ext * 2, fix = "center")
peaks.tc <- readRDS("ArchR/ArchROutput_SeuratC1-16/PeakCalls/union-peaks-TCs.rds") %>%
  subsetByOverlaps(ranges = gene.loci)

motif.info <- fread("Mus_musculus_2023_09_25_4-00_pm/TF_Information_all_motifs_plus.txt")
motif.db <- subset(motif.info, TF_Name == tf)
motif.pwms <- sapply(unique(motif.db$Motif_ID), get.pwm, simplify = F)
motif.pwms <- do.call(PWMatrixList, motif.pwms[!sapply(motif.pwms, is.null)])

mm10 <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
motif.match <- matchMotifs(
  subject = peaks.tc, genome = mm10,
  pwms = motif.pwms, out = "positions"
)
motif.match <- lapply(names(motif.match), function(m){
  x <- motif.match[[m]]
  if(length(x) == 0) return(x)
  x$Motif_ID <- m
  x$Peak <- NA
  return(x)
}) %>% do.call(what = c)
hits <- findOverlaps(query = motif.match, subject = peaks.tc)
motif.match$Peak[queryHits(hits)] <- as.character(peaks.tc)[subjectHits(hits)]
write.csv(
  x = data.frame(sort(motif.match)), row.names = F,
  file = paste0("ArchR/ArchROutput_SeuratC1-16/Plots/", name, ".csv")
)

archr.obj <- loadArchRProject("ArchR/ArchROutput_SeuratC1-16/")
archr.obj$Groups <- archr.obj$Clusters %>%
  ifelse(test = archr.obj$Clusters == "C7", yes = "C6") %>%
  ifelse(test = archr.obj$Clusters %in% paste0("C", c(6, 8:13)), yes = "C7-13")

load("palette.RData")
group.pal <- rna.pal2[c(5, 4, 2, 3, 6, 8, 11)]
names(group.pal) <- paste0("C", c(1:6, "7-13"))
useGroups <- c("C3", "C4", "C2", "C1", "C5", "C6", "C7-13")
group.pal <- group.pal[useGroups]

plotPDF(
  plotBrowserTrack(
    ArchRProj = archr.obj, region = gene.loci,
    features = GRangesList(`TC peaks` = peaks.tc, `Runx3 motif matches` = motif.match),
    groupBy = "Groups", pal = group.pal, useGroups = useGroups
  ), tileSize = 100, ylim = c(0, 1), sizes = c(8, 2, 2, 2),
  ArchRProj = archr.obj, name = name,
  width = 30, height = 12, baseSize = 15, facetbaseSize = 15
)

# TP tracks ####
addArchRThreads(threads = 32)
addArchRGenome("mm10")

archr.obj <- loadArchRProject("ArchR/ArchROutput_SeuratC1-16/")
archr.obj$Groups <- archr.obj$Clusters %>%
  ifelse(test = archr.obj$Clusters == "C7", yes = "C6") %>%
  ifelse(test = archr.obj$Clusters %in% paste0("C", c(6, 8:13)), yes = "C7-13")

load("palette.RData")
group.pal <- rna.pal2[c(5, 4, 2, 3, 6, 8, 11)]
names(group.pal) <- paste0("C", c(1:6, "7-13"))
useGroups <- c("C3", "C4", "C2", "C1", "C5", "C6", "C7-13")
group.pal <- group.pal[useGroups]

## motif ####
## Runx3 match in +7kb RORgt enhancer (match using RORgt gene body)
myAnnotation <- readRDS("ArchR/Rorc_annotation.RDS")
annot.df <- data.frame(myAnnotation$genes)
tss.start <- annot.df[annot.df$tx_external_name == "Rorc-204", ]$start
width <- 500
starts <- c(tss.start+7000-width/2, tss.start+11000-width/2)
ends <- c(tss.start+7000+width/2, tss.start+11000+width/2)
enhancer.range <- GRanges(seqnames = c("chr3", "chr3"), 
                          IRanges(starts, ends),
                          strand = "+")
gene <- "Rorct"
tf <- "Runx3"
name <- paste0("Tracks-Groups-", gene, "-7kb-locus-", tf, "-motifs")
gene.loc <- geneAnnoMm10$genes[which(geneAnnoMm10$genes$symbol == "Rorc")]
window <- 1000
ext <- round(gene.loc@ranges@width/2) + window
gene.loc <- gene.loc %>% resize(width = ext * 2, fix = "center")
peaks.tc <- readRDS("ArchR/ArchROutput_SeuratC1-16/PeakCalls/union-peaks-TCs.rds") %>%
  subsetByOverlaps(ranges = gene.loc)
# peaks.tc <- getPeakSet(archr.obj) %>%
#   subsetByOverlaps(ranges = gene.loc)

motif.info <- readRDS("Mus_musculus_2022_01_14_6-40_pm/CisBP-TF_Info.rds")
motif.db <- subset(motif.info, TF_Name == tf)
motif.pwms <- readRDS("Mus_musculus_2022_01_14_6-40_pm/CisBP-PWMs.rds")
motif.pwms <- motif.pwms[motif.db$Motif_ID]

mm10 <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
motif.match <- matchMotifs(
  subject = enhancer.range, genome = mm10,
  pwms = motif.pwms, out = "positions"
)
motif.match <- lapply(names(motif.match), function(m){
  x <- motif.match[[m]]
  if(length(x) == 0) return(x)
  x$Motif_ID <- m
  x$Peak <- NA
  return(x)
}) %>% do.call(what = c)
hits <- findOverlaps(query = motif.match, subject = peaks.tc)
motif.match$Peak[queryHits(hits)] <- as.character(peaks.tc)[subjectHits(hits)]
write.csv(
  x = data.frame(sort(motif.match)), row.names = F,
  file = paste0("ArchR/ArchROutput_SeuratC1-16/Plots/", name, ".csv")
)

## add cbfb2 peaks
read.narrowPeak <- function(fn, format = 'idr'){
  if (format == 'idr'){
    col.names = c("seqnames", "start", 'end', 'name', 'score', 'strand', 'signalValue',
                  'pval', 'qval', 'summit', 'localIDR', 'globalIDR',
                  'rep1_start', 'rep1_end', 'rep1_signalValue', 'rep1_summit',
                  'rep2_start', 'rep2_end', 'rep2_signalValue', 'rep2_summit')
  } else if (format == 'narrowPeak') {
    col.names = c("seqnames", "start", 'end', 'name', 'score', 'strand', 'signalValue',
                  'pval', 'qval', 'summit')
  } else if (format == "bed"){
    col.names = c("seqnames", "start", 'end', 'name', 'score', 'strand')
  }
  peaks <- read.table(fn, header = F, col.names = col.names)
  return(peaks)
}

# annot is chipseeker output
create.GR <- function(annot, input.is.0.indexed = F){
  mapper <- setNames(c("-", "+", "*", "*"), c("-", "+", ".", "*"))
  if (input.is.0.indexed){
    annot$start <- annot$start + 1
  }
  gr <- GRanges(seqnames = annot$seqnames, 
                IRanges(annot$start, annot$end),
                strand = mapper[annot$strand],
                annot[,!(colnames(annot) %in% 
                           c("seqnames", "ranges", "strand", "seqlevels", 
                             "seqlengths", "isCircular", "start", "end", 
                             "width", "element"))])
  return(gr)
}

cbfb2.peaks <- read.narrowPeak("../GSE90794-Cbfb2-ChIP-seq/alignment-results/alignment-results_MACS2_vsInput_peaks.filt.narrowPeak.gz", format = "narrowPeak")
cbfb2.peaks.gr <- create.GR(cbfb2.peaks, input.is.0.indexed = T)
cbfb2.peaks.gr <- cbfb2.peaks.gr %>% subsetByOverlaps(ranges = gene.loc)

##

features.to.show <- GRangesList(
  "TC peaks" = peaks.tc, "Runx3 motifs" = motif.match, 
  "enhancers" = enhancer.range, "Cbfb2 peaks" = cbfb2.peaks.gr)
plotPDF(
  plotBrowserTrack(
    ArchRProj = archr.obj, region = gene.loc,
    features = features.to.show,
    geneAnnotation = myAnnotation,
    groupBy = "Groups", pal = group.pal, useGroups = useGroups,
    tileSize = 100, ylim = c(0, 1), sizes = c(10, 1, 2, 2)),
  ArchRProj = archr.obj, name = name,
  width = 12, height = 8, baseSize = 15, facetbaseSize = 15, addDOC = F
)

### Match to peaks in Rorc gene 
myAnnotation <- readRDS("ArchR/Rorc_annotation.RDS")
annot.df <- data.frame(myAnnotation$genes)
tss.start <- annot.df[annot.df$tx_external_name == "Rorc-204", ]$start
width <- 500
starts <- c(tss.start+7000-width/2)
ends <- c(tss.start+7000+width/2)
enhancer.range <- GRanges(seqnames = c("chr3", "chr3"), 
                          IRanges(starts, ends),
                          strand = "+")
gene <- "Rorct"
tf <- "Runx3"
name <- paste0("Tracks-Groups-", gene, "-locus-", tf, "-motifs")
gene.loc <- geneAnnoMm10$genes[which(geneAnnoMm10$genes$symbol == "Rorc")]
window <- 1000
ext <- round(gene.loc@ranges@width/2) + window
gene.loc <- gene.loc %>% resize(width = ext * 2, fix = "center")
peaks.tc <- readRDS("ArchR/ArchROutput_SeuratC1-16/PeakCalls/union-peaks-TCs.rds") %>%
  subsetByOverlaps(ranges = gene.loc)
# peaks.tc <- getPeakSet(archr.obj) %>%
#   subsetByOverlaps(ranges = gene.loc)

motif.info <- readRDS("Mus_musculus_2022_01_14_6-40_pm/CisBP-TF_Info.rds")
motif.db <- subset(motif.info, TF_Name == tf)
motif.pwms <- readRDS("Mus_musculus_2022_01_14_6-40_pm/CisBP-PWMs.rds")
motif.pwms <- motif.pwms[motif.db$Motif_ID]

mm10 <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
motif.match <- matchMotifs(
  subject = peaks.tc, genome = mm10,
  pwms = motif.pwms, out = "positions"
)
motif.match <- lapply(names(motif.match), function(m){
  x <- motif.match[[m]]
  if(length(x) == 0) return(x)
  x$Motif_ID <- m
  x$Peak <- NA
  return(x)
}) %>% do.call(what = c)
hits <- findOverlaps(query = motif.match, subject = peaks.tc)
motif.match$Peak[queryHits(hits)] <- as.character(peaks.tc)[subjectHits(hits)]
write.csv(
  x = data.frame(sort(motif.match)), row.names = F,
  file = paste0("ArchR/ArchROutput_SeuratC1-16/Plots/", name, ".csv")
)

features.to.show <- GRangesList(
  "TC peaks" = peaks.tc, "Runx3 motifs" = motif.match, 
  "enhancers" = enhancer.range, "Cbfb2 peaks" = cbfb2.peaks.gr)
plotPDF(
  plotBrowserTrack(
    ArchRProj = archr.obj, region = gene.loc,
    features = features.to.show,
    geneAnnotation = myAnnotation,
    groupBy = "Groups", pal = group.pal, useGroups = useGroups,
    tileSize = 100, ylim = c(0, 1), sizes = c(10, 1, 2, 2)),
  ArchRProj = archr.obj, name = name,
  width = 12, height = 8, baseSize = 15, facetbaseSize = 15, addDOC = F
)

### Match to peaks in Prdm16 gene
load("palette.RData")
group.pal <- rna.pal2[c(5, 4, 2, 3, 6, 8, 11)]
names(group.pal) <- paste0("C", c(1:6, "7-13"))
useGroups <- c("C3", "C4", "C2", "C1", "C5", "C6", "C7-13")
group.pal <- group.pal[useGroups]

gene <- "Prdm16"
tf1 <- "Thra"
tf2 <- "Thrb"
name <- paste0("Tracks-Groups-", gene, "-locus-", tf1, "_", tf2, "-motifs")
name.1 <- paste0("Tracks-Groups-", gene, "-locus-", tf1, "-motifs")
name.2 <- paste0("Tracks-Groups-", gene, "-locus-", tf2, "-motifs")
gene.loc <- geneAnnoMm10$genes[which(geneAnnoMm10$genes$symbol == gene)]

start <- 154614056-3000
end <- 154633518+5000
prdm16.promoter <- GRanges(seqnames = c("chr4"), 
                           IRanges(start, end),
                           strand = "-")
window <- 100000
ext <- round(gene.loc@ranges@width/2) + window
gene.loc <- gene.loc %>% resize(width = ext * 2, fix = "center")

peaks.all <- getPeakSet(archr.obj) %>%
  subsetByOverlaps(ranges = gene.loc)
peaks.tc <- readRDS("ArchR/ArchROutput_SeuratC1-16/PeakCalls/union-peaks-TCs.rds") %>%
  subsetByOverlaps(ranges = gene.loc)

tf.info <- readRDS("chromVAR/All/results/top-tf-info.rds")
tf.info2 <- readRDS("chromVAR/TC/results/top-tf-info.rds")
thra.motif <- tf.info2[tf.info2$Symbol == tf1,]$top.motifs
thrb.motif <- tf.info[tf.info$Symbol == tf2,]$top.motifs
# motif.info <- readRDS("Mus_musculus_2022_01_14_6-40_pm/CisBP-TF_Info.rds")
# motif.db <- subset(motif.info, TF_Name == tf)
motif.pwms <- readRDS("Mus_musculus_2022_01_14_6-40_pm/CisBP-PWMs.rds")
motif.pwms.1 <- motif.pwms[thra.motif]
motif.pwms.2 <- motif.pwms[thrb.motif]

mm10 <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
motif.match.1 <- matchMotifs(
  subject = peaks.tc, genome = mm10,
  pwms = motif.pwms.1, out = "positions"
)
motif.match.2 <- matchMotifs(
  subject = peaks.tc, genome = mm10,
  pwms = motif.pwms.2, out = "positions"
)
motif.match.1 <- lapply(names(motif.match.1), function(m){
  x <- motif.match.1[[m]]
  if(length(x) == 0) return(x)
  x$Motif_ID <- m
  x$Peak <- NA
  return(x)
}) %>% do.call(what = c)
hits <- findOverlaps(query = motif.match.1, subject = peaks.tc)
motif.match.1$Peak[queryHits(hits)] <- as.character(peaks.tc)[subjectHits(hits)]
write.csv(
  x = data.frame(sort(motif.match.1)), row.names = F,
  file = paste0("ArchR/ArchROutput_SeuratC1-16/Plots/", name.1, ".csv")
)

motif.match.2 <- lapply(names(motif.match.2), function(m){
  x <- motif.match.2[[m]]
  if(length(x) == 0) return(x)
  x$Motif_ID <- m
  x$Peak <- NA
  return(x)
}) %>% do.call(what = c)
hits <- findOverlaps(query = motif.match.2, subject = peaks.tc)
motif.match.2$Peak[queryHits(hits)] <- as.character(peaks.tc)[subjectHits(hits)]
write.csv(
  x = data.frame(sort(motif.match.2)), row.names = F,
  file = paste0("ArchR/ArchROutput_SeuratC1-16/Plots/", name.2, ".csv")
)

features.to.show <- GRangesList(
  "TC peaks" = peaks.tc, 
  "Thra motifs" = motif.match.1, "Thrb motifs" = motif.match.2)
plotPDF(
  plotBrowserTrack(
    ArchRProj = archr.obj, region = gene.loc,
    features = features.to.show,
    groupBy = "Groups", pal = group.pal, useGroups = useGroups,
    tileSize = 100, ylim = c(0, 1), sizes = c(10, 1, 2, 2)),
  ArchRProj = archr.obj, name = name,
  width = 12, height = 8, baseSize = 15, facetbaseSize = 15, addDOC = F
)

name.short <- paste0("Tracks-Groups-", gene, "-promoter-", tf1, "_", tf2, "-motifs")
plotPDF(
  plotBrowserTrack(
    ArchRProj = archr.obj, region = prdm16.promoter,
    features = features.to.show,
    groupBy = "Groups", pal = group.pal, useGroups = useGroups,
    tileSize = 100, ylim = c(0, 1), sizes = c(10, 1, 2, 2)),
  ArchRProj = archr.obj, name = name.short,
  width = 12, height = 8, baseSize = 15, facetbaseSize = 15, addDOC = F
)

## Ichiro Runx3/Cbfb2 paper ####
gene.annot <- getGeneAnnotation(archr.obj)
gene.annot$genes <- gene.annot$genes[gene.annot$genes$symbol %in% c("Rorc")]

myAnnotation <- readRDS("ArchR/Rorc_annotation.RDS")
annot.df <- data.frame(myAnnotation$genes)
tss.start <- annot.df[annot.df$tx_external_name == "Rorc-204", ]$start
width <- 500

starts <- c(tss.start+7000-width/2, tss.start+11000-width/2)
ends <- c(tss.start+7000+width/2, tss.start+11000+width/2)
enhancer.range <- GRanges(seqnames = c("chr3", "chr3"), 
                          IRanges(starts, ends),
                          strand = "+")
gene <- "Rorct"

gene.loc <- geneAnnoMm10$genes[which(geneAnnoMm10$genes$symbol == "Rorc")]
window <- 1000
ext <- round(gene.loc@ranges@width/2) + window
gene.loc <- gene.loc %>% resize(width = ext * 2, fix = "center")
peaks.all <- getPeakSet(archr.obj) %>%
  subsetByOverlaps(ranges = gene.loc)
peaks.tc <- readRDS("ArchR/ArchROutput_SeuratC1-16/PeakCalls/union-peaks-TCs.rds") %>%
  subsetByOverlaps(ranges = gene.loc)

cbfb2.peaks <- read.narrowPeak("../GSE90794-Cbfb2-ChIP-seq/alignment-results/alignment-results_MACS2_vsInput_peaks.filt.narrowPeak.gz", format = "narrowPeak")
cbfb2.peaks.gr <- create.GR(cbfb2.peaks, input.is.0.indexed = T)
cbfb2.peaks.gr <- cbfb2.peaks.gr %>% subsetByOverlaps(ranges = gene.loc)

useGroups <- c("C3", "C4", "C2", "C1", "C7-13")

name <- paste0("Tracks-Groups-", gene, "-locus")
features.to.show <- GRangesList(
  "TC peaks" = peaks.tc, #'All peaks' = peaks.all,
  "enhancers" = enhancer.range, "Cbfb2 peaks" = cbfb2.peaks.gr)
plotPDF(
  plotBrowserTrack(
    ArchRProj = archr.obj, region = gene.loc,
    features = features.to.show,
    geneAnnotation = gene.annot,
    groupBy = "Groups", 
    pal = group.pal, useGroups = useGroups,
    tileSize = 100, ylim = c(0, 1), sizes = c(10, 1, 2, 1)),
  ArchRProj = archr.obj, name = name,
  width = 12, height = 8, baseSize = 15, facetbaseSize = 15, addDOC = F
)

name <- paste0("Tracks-Groups-", gene, "-locus-wo-enhancers")
features.to.show2 <- GRangesList(
  "TC peaks" = peaks.tc, 
  "Cbfb2 peaks" = cbfb2.peaks.gr)
plotPDF(
  plotBrowserTrack(
    ArchRProj = archr.obj, region = gene.loc,
    features = features.to.show2,
    geneAnnotation = gene.annot,
    groupBy = "Groups", 
    pal = group.pal, useGroups = useGroups,
    tileSize = 100, ylim = c(0, 1), sizes = c(10, 1, 2, 1)),
  ArchRProj = archr.obj, name = name,
  width = 12, height = 8, baseSize = 15, facetbaseSize = 15, addDOC = F
)


plot.clusters <- function(
    SRO = NULL, idx = NULL,
    vis = SRO@reductions$umap@cell.embeddings, 
    axis.titles = c("UMAP1", "UMAP2"),
    groups = SRO$Clusters, clusters.col = 'Clusters', 
    col = pal[[clusters.col]], new.levels = NULL,
    pref.C = T, labels = T, line = T,
    point.size = 1.5, point.alpha = 0.7,
    label.size = point.size * 2, label.pad = 1, border.size = 0.5, 
    legend.size = 3, show_legend = T
){
  if(class(groups) == "factor"){ 
    clusters <- groups
  } else {
    clusters <- factor(groups, levels = sort(unique(groups)))
  }
  gp <- ggplot() + theme_classic() +
    labs(x = axis.titles[1], y = axis.titles[2], color = clusters.col) +
    theme(axis.ticks = element_blank(), axis.text = element_blank())
  if(!line) gp <- gp + theme(axis.line = element_blank())
  if(!is.null(idx)){
    vis <- vis[idx, ]
    if(is.null(new.levels)){
      new.levels <- as.character(sort(unique(clusters[idx])))
    }
    clusters <- factor(as.character(clusters[idx]), levels = new.levels)
    if(!is.null(col)) col <- col[new.levels]
  }
  gp <- gp + geom_point(
    mapping = aes(x = vis[, 1], y = vis[, 2], color = clusters),
    size = point.size, alpha = point.alpha, show.legend = show_legend
  )
  if(!is.null(col)) gp <- gp + scale_color_manual(values = col)
  if(labels){
    C <- levels(clusters)
    CL <- paste0(ifelse(pref.C, "C", ""), levels(clusters))
    vis.cl <- data.frame(matrix(NA, nrow = length(C), ncol = 2, dimnames = list(C, NULL)))
    for(c in C){
      vis.cl[c, 1] <- median(vis[clusters == c, 1], na.rm = T)
      vis.cl[c, 2] <- median(vis[clusters == c, 2], na.rm = T)
    }
    gp + geom_label_repel(
      aes(x = vis.cl[, 1], vis.cl[, 2], color = C, label = CL),
      color = 'black', seed = 1,
      label.size = border.size,
      size = label.size, label.padding = unit(label.pad, "mm"), show.legend = F
    ) +
      geom_label_repel(
        aes(x = vis.cl[, 1], vis.cl[, 2], color = C, label = CL),
        size = label.size, seed = 1,
        label.size = NA, label.padding = unit(label.pad, "mm"), show.legend = F
      )+ guides(colour = guide_legend(override.aes = list(size = legend.size)))
  } else gp + guides(colour = guide_legend(override.aes = list(size = legend.size)))
}

library(ggrepel)

umap <- read.csv("~/0-workspace/CCR7_DC/MLN_RORgt_MHCII_multiome/Seurat/results/UMAP_C1-16-TConly.csv", row.names = 1)
idx <- grepl("TC", sro$Clusters.annot)


p <- plot.clusters(
    SRO = sro, idx = rownames(umap), vis = umap,
    groups = sro$Clusters.annot, clusters.col = 'Annotations', 
    col = rna.pal, pref.C = F, labels = T, line = T,
    point.size = 1.5, point.alpha = 0.7,
    label.size = 4, label.pad = 1, border.size = 0.5, 
    legend.size = 3, show_legend = T)
p.c <- plot.clusters(
  SRO = sro, idx = rownames(umap),vis = umap,
  groups = sro$Clusters, clusters.col = 'Clusters', 
  col = rna.pal2, pref.C = T, labels = T, line = T,
  point.size = 1.5, point.alpha = 0.7,
  label.size = 4, label.pad = 1, border.size = 0.5, 
  legend.size = 3, show_legend = T)
pdf("Seurat/plots/UMAP-annotations.pdf", width = 8, height = 6)
print(p)
print(p.c)
dev.off()

sro.tc <- sro[, idx]
sro.tc$Clusters.annot <- droplevels(sro.tc$Clusters.annot)
sro.tc$Clusters <- droplevels(sro.tc$Clusters)
sro.tc <- sro.tc %>% RunUMAP(dims = 1:30, n.neighbors = 30, metric = "cosine", min.dist = 0.4, spread = 1)
p3 <- plot.clusters(
  SRO = sro.tc, 
  groups = sro.tc$Clusters.annot, clusters.col = 'Annotations', 
  col = rna.pal, pref.C = F, labels = T, line = T,
  point.size = 1.5, point.alpha = 0.7,
  label.size = 4, label.pad = 1, border.size = 0.5, 
  legend.size = 3, show_legend = T)
p4<- plot.clusters(
  SRO = sro.tc, 
  groups = sro.tc$Clusters, clusters.col = 'Clusters', 
  col = rna.pal2, pref.C = T, labels = T, line = T,
  point.size = 1.5, point.alpha = 0.7,
  label.size = 4, label.pad = 1, border.size = 0.5, 
  legend.size = 3, show_legend = T)
pdf("Seurat/plots/UMAP-annotations2.pdf", width = 8, height = 6)
print(p3)
print(p4)
dev.off()


# Subcluster TC I ####
## TC I contains TC I and transitional
addArchRThreads(threads = 32)
addArchRGenome("mm10")

## SRO ####
sro <- readRDS("Seurat/results/SRO.rds")
idx <- sro$Clusters.annot == "TC 1"
sro.subset <- sro[, idx]
subres.list <- seq(0.1, 1, 0.1)
sro.subset <- ScaleData(sro.subset, features = rownames(sro.subset)) %>%
  RunPCA(features = rownames(sro.subset), npcs = 50) %>%
  FindNeighbors(dims = 1:30, k.param = 30) %>%
  FindClusters(resolution = subres.list)

table(sro.subset$RNA_snn_res.0.8)

master.res <- "Clusters"
for (res in subres.list){
  subset.colname <- paste0("RNA_snn_res.", res)
  newcolname <- paste0("Clusters_subres.", res)
  tmp <- sro.subset@meta.data[colnames(sro), ][[subset.colname]]
  subclusters <- ifelse(
    test = is.na(tmp),
    yes = tmp,
    no = paste0("sub", tmp)
  )
  sro@meta.data[[newcolname]] <- coalesce(subclusters, sro@meta.data[[master.res]])
}
table(sro$Clusters_subres.0.8)

p <- plot.clusters(
  SRO = sro, 
  groups = sro$Clusters_subres.0.8, clusters.col = 'Clusters_subres.0.8', 
  col = pal$Clusters, pref.C = F, labels = T, line = T,
  point.size = 1.5, point.alpha = 0.7,
  label.size = 4, label.pad = 1, border.size = 0.5, 
  legend.size = 3, show_legend = T)
pdf("Seurat/plots/UMAP-Clusters_subres.0.8.pdf", width = 8, height = 6)
print(p)
dev.off()

module.scores <- read.csv("Seurat/results/module_scores.csv", row.names = 1)
module.scores <- module.scores[c(11:16)]
sro@meta.data <- cbind(sro@meta.data, module.scores)
saveRDS(sro, "Seurat/results/SRO.rds")
write.csv(sro@meta.data, "Seurat/results/meta-data.csv")

pdf("Seurat/plots/vln-transitional-tc1-markers.pdf", width = 12, height = 6)
VlnPlot(object = sro, features = 'transitional.markers1', group.by = "ATAC.Clusters.subC3")
VlnPlot(object = sro, features = 'tc1.markers1', group.by = "ATAC.Clusters.subC3")
VlnPlot(object = sro, features = 'transitional.markers1', group.by = "Clusters_subres.0.8")
VlnPlot(object = sro, features = 'tc1.markers1', group.by = "Clusters_subres.0.8")
dev.off()

## ArchR ####
archr.obj <- loadArchRProject("ArchR/ArchROutput_SeuratC1-16/")
idx <- archr.obj$cellNames[archr.obj$Clusters == "C3"]
archr.obj.c3 <- subsetArchRProject(
  ArchRProj = archr.obj, cells = idx,
  outputDirectory = "ArchR/ArchROutput_C3",
  dropCells = TRUE, logFile = NULL, threads = getArchRThreads(), force = FALSE
)

archr.obj.c3 <- addIterativeLSI(archr.obj.c3, iterations = 5, varFeatures = 100000, scaleDims = T, force = T)
res.list <- seq(0.1, 1, 0.1)
for (res in res.list){
  archr.obj.c3 <- archr.obj.c3 %>% 
    addClusters(
      k.param = 30, resolution = res, maxClusters = 50,
      force = T, name = paste0("Clusters.res.", res)
    )
}

md <- data.frame(archr.obj@cellColData)
md$Clusters.subC3 <- data.frame(archr.obj.c3@cellColData)[archr.obj$cellNames, ]$Clusters.res.0.2
md$Clusters.subC3 <- ifelse(is.na(md$Clusters.subC3), md$Clusters.subC3, paste0("sub", md$Clusters.subC3))
md$Clusters.subC3 <- coalesce(md$Clusters.subC3, md$Clusters)
archr.obj$Clusters.subC3 <- md$Clusters.subC3
rownames(md) <- gsub("^.+#", "", rownames(md))

sro$ATAC.Clusters.subC3 <- md[rownames(sro@meta.data), ]$Clusters.subC3
table(sro$ATAC.Clusters.subC3)

pal <- list(
  Clusters = c(
    "#00bf00", "#489de8", "#d40663", "#f8c72f", "#077315",
    "#785cd4", "#e67109", "#0eefff", "#f081e6", "#260691",
    "#49709c", "#9e7d3f", "#bd537a", "#4e225c", "#f202ed",
    "#fec55f", "#062e0b", "#9c6fa8", "#078d94", "#5c1a1a",
    "#827c68", "#aebeff", "#9c2903", "#ffc5af", "#4f5715",
    "#0249f0", "#f43525", "#0077ff", "#7f227e", "#dfddff",
    "#7e85d7", "#fff64f", "#5fed0e", "#543018", "#f31220"
  )
)

p <- plot.clusters(
  SRO = sro, 
  groups = sro$ATAC.Clusters.subC3, clusters.col = 'ATAC.Clusters.subC3', 
  col = pal$Clusters, pref.C = F, labels = T, line = T,
  point.size = 1.5, point.alpha = 0.7,
  label.size = 4, label.pad = 1, border.size = 0.5, 
  legend.size = 3, show_legend = T)
pdf("Seurat/plots/UMAP-ATAC.Clusters.subC3.pdf", width = 8, height = 6)
print(p)
dev.off()

saveArchRProject(archr.obj)

## Re-group cells and run chromVAR ####
archr.obj <- loadArchRProject("ArchR/ArchROutput_SeuratC1-16/")
setwd("ArchR/")
saveArchRProject(ArchRProj = archr.obj, outputDirectory = "ArchROutput_SeuratC1-16-subC", load = FALSE)
archr.obj <- loadArchRProject("ArchROutput_SeuratC1-16-subC/")

table(archr.obj$Clusters.subC3)

table(archr.obj$Groups)
table(archr.obj$Clusters)

archr.obj$Clusters2 <- ifelse(archr.obj$Clusters == "C7", "C6",
                              ifelse(archr.obj$Clusters == "C6", "C7",
                                     archr.obj$Clusters))

archr.obj$Groups.new <- archr.obj$Clusters.subC3 %>%
  ifelse(test = archr.obj$Clusters.subC3 == "C7", yes = "C6") %>%
  ifelse(test = archr.obj$Clusters.subC3 %in% paste0("C", c(6, 8:13)), yes = "C7-13") %>%
  ifelse(test = archr.obj$Clusters.subC3 %in% paste0("subC", c(1, 2, 4, 5)), yes = "subC1245")

table(archr.obj$Groups.new)

arrow.file <- "ArchROutput_SeuratC1-16-subC/ArrowFiles/3228_Sample_CB-1288_MLN_RORgt_MHCII_multiome_ATAC_IGO_12437_C_19.arrow"
group.nFrags <- sapply(unique(archr.obj$Groups.new), function(g){
  cat("Extracting fragments for", g, "\n")
  group.cells <- subset(archr.obj@cellColData, Groups.new == g) %>% rownames()
  group.frags <- getFragmentsFromArrow(ArrowFile = arrow.file, cellNames = group.cells)
  cbind(as.character(seqnames(group.frags)), start(group.frags), end(group.frags)) %>%
    write.table(file = paste0("ArchROutput_SeuratC1-16-subC/GroupFragments/", g, ".bed"),
                row.names = F, col.names = F, quote = F, sep = "\t")
  return(length(group.frags))
})
## call peaks first
peaks <- readRDS("ArchROutput_SeuratC1-16-subC/PeakCalls/union-peaks.rds")
archr.obj <- addPeakSet(archr.obj, peakSet = peaks, force = T) %>% addPeakMatrix()
# Groups: 176952 peaks with maximum length 1158 and median length 206
# Groups.new: 172973 peaks with maximum length 1129 and median length 206

saveArchRProject(archr.obj)
getMatrixFromProject(archr.obj, useMatrix = "PeakMatrix", binarize = F) %>%
  saveRDS("ArchROutput_SeuratC1-16-subC/PeakCalls/counts.rds")


table(archr.obj$Clusters.scRNAseq)
setwd("ArchR")
tc.cells <- archr.obj$cellNames[archr.obj$Clusters.scRNAseq %in% paste0("C", 1:5)]
tc.peaks <- readRDS(file = "ArchROutput_SeuratC1-16-subC/PeakCalls/union-peaks-TCs.rds")

tc.archr.obj <- subsetArchRProject(
  ArchRProj = archr.obj, cells = tc.cells,
  outputDirectory = "ArchROutput_TCs/"
) %>% addPeakSet(peakSet = tc.peaks, force = T) %>% addPeakMatrix(force = T)
saveArchRProject(tc.archr.obj)
getMatrixFromProject(tc.archr.obj, useMatrix = "PeakMatrix", binarize = F) %>%
  saveRDS("ArchROutput_TCs/PeakCalls/counts.rds")


# export BW files ####
archr.obj <- loadArchRProject("ArchROutput_SeuratC1-16/", showLogo = F)
table(archr.obj$Groups)
cl.to.annot <- c("C1" = "TC IV", "C2" = "TC III", "C3" = "TC I", "C4" = "TC II", 
                 "C5" = "Proliferating NCR+ ILC3", "C6" = "NRC+ ILC3", "C7-13" = "LTi")
archr.obj@cellColData$ATAC.annotations <- cl.to.annot[archr.obj$Groups]
table(archr.obj@cellColData$ATAC.annotations)
getGroupBW(
  ArchRProj = archr.obj,
  groupBy = "ATAC.annotations",
  normMethod = "ReadsInTSS",
  tileSize = 100,
  maxCells = 1000,
  ceiling = 4
)

rtracklayer::export.bed(getPeakSet(archr.obj), "ArchROutput_SeuratC1-16/union_peaks.bed")

tc.peaks <- readRDS(file = "ArchROutput_SeuratC1-16-subC/PeakCalls/union-peaks-TCs.rds")
rtracklayer::export.bed(tc.peaks, "ArchROutput_SeuratC1-16/tc_peaks.bed")

# dot plot ####
RenameGenesSeurat <- function(obj, newnames) { # Replace gene names in different slots of a Seurat object. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA
  
  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]] <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]] <- newnames
    if (length(RNA@scale.data)) rownames(RNA@scale.data) <- newnames
    if (length(RNA@meta.features)) {
      RNA@meta.features$prev.names <- rownames(RNA@meta.features)
      rownames(RNA@meta.features) <- newnames
    }
    if (length(RNA@var.features)){
      RNA@var.features <- rownames(RNA@meta.features[match(RNA@var.features, RNA@meta.features$prev.names), ])
    }
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}
sro <- readRDS("Seurat/results/SRO.rds")
sro@assays$RNA@meta.features$Symbol.unique <- make.unique(sro@assays$RNA@meta.features$Symbol)
sro <- RenameGenesSeurat(sro, sro@assays$RNA@meta.features$Symbol.unique)

sro.subset <- subset(sro, Clusters %in% 1:16)
load("palette.RData")

rna.annot <- c(names(rna.annot)[1:8], rep("LTi", 8))
cl <- read.csv("Seurat/results/meta-data.csv", header = 1, row.names = 1) %>% subset(Clusters %in% 1:16)
cl$Clusters <- ifelse(cl$Clusters %in% 9:16, "9-16", as.character(cl$Clusters)) %>% factor(levels = c(as.character(1:8), "9-16"))
cl$Clusters.annot <- factor(gsub("LTi Variation \\d+$", "LTi", cl$Clusters.annot), levels = rna.annot[1:9])
cl <- cl[order(cl$Clusters, decreasing = T), ]

sro.subset$Clusters.annot <- cl[rownames(sro.subset@meta.data), ]$Clusters.annot

# DR3 (Tnfrsf25), OX40L (Tnfsf4), and Bhlhe40
geneList <- c("Tnfrsf25", "Tnfsf4", "Bhlhe40"); 
geneList %in% rownames(sro.subset)
fn <- "Tnfrsf25_Tnfsf4_Bhlhe40"
group.name <- "Clusters.annot"

ggsave(
  filename = paste0("Seurat/plots/gene-expr/dot-", fn, ".pdf"), width = 8, height = 5,
  plot = DotPlot(sro.subset, assay = "RNA", group.by = group.name,
                 features = geneList, cols = c("blue", "red"), dot.scale = 6) +
    scale_x_discrete(breaks = geneList, labels = geneList) +
    theme(axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 11)) + labs(x = "", y = group.name) + coord_flip()
)


# violin ####
# is it possible, for the RORgt MHCII multiome dataset to have a violin plot of 
# the LTi grouped cluster, the NCR1+ cluster and the TC I, II, III, IV (not Ki67) cluster 
# for Prdm16 expression (unimputed) using the same colors we used in the original UMAP?
sro <- readRDS("Seurat/results/SRO.rds")
sro@assays$RNA@meta.features$Symbol.unique <- make.unique(sro@assays$RNA@meta.features$Symbol)
sro <- RenameGenesSeurat(sro, sro@assays$RNA@meta.features$Symbol.unique)

sro.subset <- subset(sro, Clusters %in% 1:16)
load("palette.RData")
rna.pal3 <- c(rna.pal[1:8], LTi = as.character(rna.pal[11]))

rna.annot <- c(names(rna.annot)[1:8], rep("LTi", 8))
cl <- read.csv("Seurat/results/meta-data.csv", header = 1, row.names = 1) %>% subset(Clusters %in% 1:16)
cl$Clusters <- ifelse(cl$Clusters %in% 9:16, "9-16", as.character(cl$Clusters)) %>% factor(levels = c(as.character(1:8), "9-16"))
cl$Clusters.annot <- factor(gsub("LTi Variation \\d+$", "LTi", cl$Clusters.annot), levels = rna.annot[1:9])
cl <- cl[order(cl$Clusters, decreasing = T), ]

sro.subset$Clusters.annot <- cl[rownames(sro.subset@meta.data), ]$Clusters.annot
sro.subset2 <- subset(sro.subset, Clusters.annot %in% c(
  "TC 1", "TC 2", "TC 3", "TC 4", "NCR+ ILC3", 'LTi'
))
sro.subset2$Clusters.annot <- droplevels(sro.subset2$Clusters.annot)

v.p <- VlnPlot(sro.subset2, features = "Prdm16", pt.size = 0, group.by = "Clusters.annot", 
        cols = rna.pal3)
ggsave("Seurat/plots/vln-Prdm16.pdf", v.p, width = 8, height = 4)
