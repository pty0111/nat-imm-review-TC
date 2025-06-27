# integrate all datasets
Sys.setenv(RETICULATE_PYTHON = "/data1/lesliec/tyler/utils/miniforge3/envs/multiome/bin/python")
setwd("~/0-workspace/CCR7_DC/oral-tolerance-integrate/")

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
    library(anndata)
    library(rasterpdf)
})
reticulate::py_require("anndata")

set.seed(1)
options(future.globals.maxSize = Inf)

source("utils.R")


load("../MLN_RORgt_MHCII_multiome/palette.RData")
pal <- list(annotation = rna.pal)
pal$Clusters = c(
    "#00bf00", "#489de8", "#d40663", "#f8c72f", "#077315",
    "#785cd4", "#e67109", "#0eefff", "#f081e6", "#260691",
    "#49709c", "#9e7d3f", "#bd537a", "#4e225c", "#f202ed",
    "#fec55f", "#062e0b", "#9c6fa8", "#078d94", "#5c1a1a",
    "#827c68", "#aebeff", "#9c2903", "#ffc5af", "#4f5715",
    "#0249f0", "#f43525", "#0077ff", "#7f227e", "#dfddff",
    "#7e85d7", "#fff64f", "#5fed0e", "#543018", "#f31220"
)
pal$Clusters_long = c(
    "#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
    "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
    "#5A0007", "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
    "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
    "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
    "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
    "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
    "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
    "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
    "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
    "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
    "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
    "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")


# 3. plots ####
pref.i <- "integrate-ILC3-TC-v5/"; dir.create(pref.i)
reduction <- "cca"; pref <- paste0(pref.i, reduction, "/"); dir.create(pref)
pref.p <- paste0("plots/", pref); dir.create(pref.p, recursive = T)
sro.i <- readRDS(paste0(pref, "SRO.rds"))

group.name <- 'sample'
groups <- sort(unique(sro.i@meta.data[[group.name]]))
groups

sources <- c(
    "MLN_RORgt_MHCII_multiome", "TC_all_LN", "Colonna",
    "Kedmi", "Lyu", "Wang",
    "Littman", "Gardner.A", "Gardner.E"
)

pl <- list()
for (i in 1:length(sources)){
    source.curr <- sources[i]
    p.curr <- plot.clusters.highlight.one(
        SRO = sro.i,
        idx1 = colnames(sro.i), idx2 = sro.i$orig.ident == source.curr,
        groups = get.named.vector.sro(sro.i, group.name),
        clusters.col = source.curr,
        labels = F, label.size = 5,
        col = pal[[group.name]], pref.C = F) + guides(colour = guide_legend(override.aes = list(size=5), ncol = 1))
    pl[[i]] <- p.curr
}

ncols <- 3
nrows <- ceiling(length(pl)/ncols)

pdf(file = paste0(pref.p, "UMAP-", group.name, ".pdf"), width = 10*ncols, height = 6*nrows)
plot_grid(plotlist = pl, ncol = ncols, align = 'hv')
dev.off()

res0.4.to.anno <- c('Colonna C7' = 'RORgt DC I', 'Colonna C5' = 'RORgt DC II', 'Colonna C9' = 'RORgt DC III', 'Colonna C10' = 'RORgt DC IV', 'Colonna C0' = 'ILC3')

sro.i@meta.data[sro.i@meta.data$orig.ident == 'Colonna',]$paper.annot <- res0.4.to.anno[sro.i@meta.data[sro.i@meta.data$orig.ident == 'Colonna',]$Cluster.prev]

group.name <- 'Cluster.prev'
groups <- sort(unique(sro.i@meta.data[[group.name]]))
groups

sources <- c(
    "MLN_RORgt_MHCII_multiome", "TC_all_LN", "Colonna",
    "Kedmi", "Lyu", "Wang",
    "Littman", "Gardner.A", "Gardner.E"
)

pl <- list()
for (i in 1:length(sources)){
    source.curr <- sources[i]
    p.curr <- plot.clusters.highlight.one(
        SRO = sro.i,
        idx1 = colnames(sro.i), idx2 = sro.i$orig.ident == source.curr,
        groups = get.named.vector.sro(sro.i, group.name),
        clusters.col = source.curr,
        labels = F, label.size = 5,
        col = pal[[group.name]], pref.C = F) + guides(colour = guide_legend(override.aes = list(size=5), ncol = 1))
    pl[[i]] <- p.curr
}

ncols <- 3
nrows <- ceiling(length(pl)/ncols)

raster_pdf(file = paste0(pref.p, "UMAP-", group.name, ".pdf"), width = 10*ncols, height = 6*nrows, res = 150)
plot_grid(plotlist = pl, ncol = ncols, align = 'hv')
dev.off()



# write.csv(sro.i@meta.data, file = paste0(pref, "meta-data.csv"))
# saveRDS(sro.i, paste0(pref, "SRO.rds"))

gardner.e.md <- read.csv("../oral-tolerance-Gardner/Seurat/early/meta-data.csv", row.names = 1)
gardner.e.md <- gardner.e.md %>% subset(RNA_snn_res.0.2 %in% c(1, 7, 10))
rownames(gardner.e.md) <- paste0(rownames(gardner.e.md), "_9")
gardner.e.md$paper.annot[is.na(gardner.e.md$paper.annot)] <- 'na'

table(gardner.e.md$paper.annot)

sro.i@meta.data[sro.i@meta.data$orig.ident == 'Gardner.E',]$paper.annot <- gardner.e.md[rownames(sro.i@meta.data[sro.i@meta.data$orig.ident == 'Gardner.E',]), ]$paper.annot

kedmi.wang.cl.to.annot <- c(
    'Kedmi C14' = 'LTi-like ILC3',
    'Kedmi C12' = 'Ki67+',
    'Kedmi C6' = 'JC1',
    'Kedmi C9' = 'JC2',
    'Wang C13' = 'LTi-like ILC3',
    'Wang C9' = 'Ki67+',
    'Wang C3' = 'JC1',
    'Wang C12' = 'JC2'
)

sro.i$paper.annot <- ifelse(
    test = sro.i$Cluster.prev %in% names(kedmi.wang.cl.to.annot),
    yes = as.character(kedmi.wang.cl.to.annot[sro.i$Cluster.prev]),
    no = sro.i$paper.annot
)

table(sro.i$paper.annot)[table(sro.i$paper.annot) < 20]

annot1.to.annot3 <- c(
    "ILC3" = "ILC3",
    "ILC3p" = "ILC3",
    "Ki67+ TC" = "Ki67+ TC",
    "LTi" = "LTi", 
    "LTi Variation 1" = "LTi Variation 1", "LTi Variation 2" = "LTi Variation 2",
    "LTi Variation 3" = "LTi Variation 3", "LTi Variation 4" = "LTi Variation 4",
    "LTi Variation 5" = "LTi Variation 5", "LTi Variation 6" = "LTi Variation 6",
    "LTi Variation 7" = "LTi Variation 7", "LTi Variation 8" = "LTi Variation 8",
    "LTi-like ILC" = "ILC", "LTi-like ILC3s" = "ILC3",
    'LTi-like ILC3' = 'ILC3',
    "LTi-like R-eTAC" = "LTi-like eTAC",
    "NCR+ ILC3" = "ILC3",
    "Nrg1_Pos" = "Nrg1_Pos",
    "Pro. ILC" = "ILC",
    "Pro. R-eTAC" = "Ki67+ eTAC",
    "Proliferating NCR+ ILC3" = "ILC3",
    "Proliferating TC" = "Ki67+ TC",
    "R-eTAC1" = "eTAC I",
    "R-eTAC2" = "eTAC II", 
    "R-eTAC3" = "eTAC III",
    "TC 1" = "TC I", "TC 2" = "TC II", "TC 3" = "TC III", "TC 4" = "TC IV",
    "TC I" = "TC I", "TC II" = "TC II", "TC III" = "TC III", "TC IV" = "TC IV",
    "eTACs I" = "eTAC I",
    "eTACs II" = "eTAC II",
    'Proliferating eTAC' = 'Ki67+ eTAC',
    "ILC" = "ILC", "ILC1" = "ILC1", "ILC2" = "ILC2",
    "na" = "na",
    "Tingible body macs" = "other",
    "T cell zone macs" = "other",
    "Mig_DC_1" = "other", 
    "R-cDC1" = "other", 
    "R-cDC2" = "other",
    "cDC2A" = "other",
    "R-mDC" = "other"
)


unique(sro.i$paper.annot)[unique(sro.i$paper.annot) %ni% names(annot1.to.annot3)]

names(annot1.to.annot3)[names(annot1.to.annot3) %ni% unique(sro.i$paper.annot)]

sro.i$paper.annot3 <- ifelse(
    test = sro.i$paper.annot %in% names(annot1.to.annot3),
    yes = as.character(annot1.to.annot3[sro.i$paper.annot]),
    no = sro.i$paper.annot
)
group.name <- 'paper.annot3'
groups <- sort(unique(sro.i@meta.data[[group.name]]))
groups[groups %ni% names(pal$paper.annot3)]

pal$paper.annot3 <-c(
    "Proliferating TC" = "#e051bc",
    "Ki67 TC" = "#e051bc",
    "Ki67+ TC" = "#e051bc",
    "Ki67+" = "#e051bc",
    "Ki67+ eTAC" = "#e051bc",
    "early/transitional TC" = "#ccf3ff",
    "TC I" = "#D790FF",
    "TC II" = "#BC23FF",
    "TC III" = "#72418F",
    "TC IV" = "#3A2465",
    "TC II, III, IV" = "#72418F",
    "TCs" = "#BC23FF",
    "RORgt DC I" = "#D790FF",
    "RORgt DC II" = "#BC23FF",
    "RORgt DC III" = "#72418F",
    "RORgt DC IV" = "#3A2465",
    "LTi-like eTAC" = "#ccf3ff",
    "Proliferating eTAC" = "#e051bc",
    "RORgt+_DC_like" = '#D790FF',
    "RORgt+_eTAC" = '#72418F',
    "JC1" = "#D790FF",
    "JC2" = "#BC23FF",
    "eTAC I" = "#D790FF",
    "eTAC II" = "#BC23FF",
    "eTAC III" = "#72418F",
    "JC"= "#BC23FF",
    "ILC1" ='#17e81b',
    "ILC2" ='#77d9ae',
    "ILC3p" = "#7daedb",
    "Proliferating NCR+ ILC3" = "#0976de",
    "ILC3" = "#14a38e",
    "NCR+ ILC3" = "#003a85",
    "LTi-like ILC3s" = "#14a38e",
    "LTi-like ILC3" = "#14a38e",
    "LTi Variation 1" = "#17e81b",
    "LTi Variation 2" = "#77d9ae",
    "LTi Variation 3" = "#14a38e",
    "LTi Variation 4" = "#05a15c",
    "LTi Variation 5" = "#4e7506",
    "LTi Variation 6" = "#0d6930",
    "LTi Variation 7" = "#405933",
    "LTi Variation 8" = "#0d3b0e",
    'LTi' = "#14a38e",
    "ILC" ='#14a38e',
    "Nrg1_Pos" = '#D790FF',
    "Tolerogenic DC" = '#72418F',
    "mDC" = '#FAD09F',
    "pDC" = '#e67109',
    "tDC" = '#ffc5af',
    "migDC" = 'orange',
    "cDC1" = 'orange',
    "cDC2" = '#e60e0e',
    "cDC2A" = '#e60e0e',
    "cDC2_A" = '#e60e0e',
    "cDC2_B" = 'pink',
    "Tingible body macs" =  '#0077ff',
    "T cell zone macs" = '#489de8',
    "na" = 'black',
    'other' = 'yellow'
    # 'other' = '#636363'
)

sro.i$paper.annot3 <- factor(sro.i$paper.annot3, levels = names(pal$paper.annot3))
sources <- c(
    "MLN_RORgt_MHCII_multiome", "TC_all_LN", "Colonna",
    "Kedmi", "Lyu", "Wang",
    "Littman", "Gardner.A", "Gardner.E"
)

pl <- list()
for (i in 1:length(sources)){
    source.curr <- sources[i]
    p.curr <- plot.clusters.highlight.one(
        SRO = sro.i,
        idx1 = colnames(sro.i), idx2 = sro.i$orig.ident == source.curr,
        groups = get.named.vector.sro(sro.i, group.name),
        clusters.col = source.curr,
        labels = F, label.size = 5,
        col = pal[[group.name]], pref.C = F) + guides(colour = guide_legend(override.aes = list(size=5), ncol = 1))
    pl[[i]] <- p.curr
}

ncols <- 3
nrows <- ceiling(length(pl)/ncols)

raster_pdf(file = paste0(pref.p, "UMAP-", group.name, ".pdf"), width = 10*ncols, height = 6*nrows, res = 150)
plot_grid(plotlist = pl, ncol = ncols, align = 'hv')
dev.off()


