# MLN
suppressPackageStartupMessages({
  library(ArchR)
  library(Seurat)
  library(Signac)
  library(Rphenograph)
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
})
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
