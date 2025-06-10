Sys.setenv(RETICULATE_PYTHON = "/data1/lesliec/tyler/utils/miniforge3/envs/multiome/bin/python")
setwd("~/0-workspace/CCR7_DC/oral-tolerance-Littman/")

suppressPackageStartupMessages({
  library(Seurat) # v4.4
  library(Matrix)
  library(future)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(ComplexHeatmap)
  library(circlize)
  library(ggrepel)
})

set.seed(1)
options(future.globals.maxSize = Inf)

pal <- readRDS("plots/palette.rds")

source("utils.R")

# ############################################################################ #
# MAST ####
# ############################################################################ #
pref.sro <- "Seurat/merged/"; pref.p.sro <- "plots/Seurat/merged/"
dir.create(paste0(pref.sro, "markers"))
dir.create(paste0(pref.p.sro, "markers"))

sro <- readRDS(paste0(pref.sro, "SRO.rds"))
imp.expr <- readRDS(paste0(pref.sro, "imputed-expr.rds"))

# ############################################################################ #
# 1. all one-vs-rest for option2 with C15 divided into 2 subC ####

group.name <- "RNA_snn_res.1"
fn <- paste0(group.name, "-all")
n <- 30

sro.input <- sro

# res <- MAST.DE.multiple(fn, sro.input, idents = group.name, output.dir = paste0(pref.sro, 'markers/'))
markers <- select.markers(fn, output.dir = paste0(pref.sro, 'markers/'), n=n, pairwise = F)
# write.csv(markers, paste0(pref.sro, 'markers/', fn, '-top-', n, '-markers.csv'), row.names = T, quote = F)

all.markers <- select.markers(fn, output.dir = paste0(pref.sro, 'markers/'), pairwise = F)
# write.csv(all.markers, paste0(pref.sro, 'markers/', fn, '-sig-markers.csv'), row.names = T, quote = F)

head(sro.input@meta.data)
pal[[group.name]] <- pal$Clusters[1:length(unique(sro.input@meta.data[[group.name]]))]
names(pal[[group.name]]) <- sort(unique(sro.input@meta.data[[group.name]]))

pal[['sample']] <- pal$Clusters[1:length(unique(sro.input@meta.data[['sample']]))]
names(pal[['sample']]) <- sort(unique(sro.input@meta.data[['sample']]))

ca.col.list <- list(pal[[group.name]], pal$RabiClusters, pal$sample)
names(ca.col.list) <- c(group.name, "RabiClusters", 'sample')

h <- nrow(markers)*0.2+2
pdf(width = 40, height = h, file = paste0(pref.p.sro, "markers/", fn, "-top-", n, ".pdf"))
DE.heatmap(sro.input, expr = imp.expr, legend.name = "scaled\nimputed\nexpression",
           cells = colnames(sro.input), markers = markers$gene, md = sro.input@meta.data,
           split = group.name, cluster_column_slices = T, column_title_rot = 45,
           rowsplit = markers$cluster, cluster_rows = T, cluster_row_slices = T,
           split.columns = names(ca.col.list),
           ca.col = ca.col.list)
dev.off()
