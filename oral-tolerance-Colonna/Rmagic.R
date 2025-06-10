setwd("~/0-workspace/CCR7_DC/oral-tolerance-Colonna/")

library(Seurat)
library(Rmagic)
library(reticulate)
py_require("magic-impute")

# sro <- readRDS("Seurat/merged-5A/SRO.rds")
# sro.imp <- magic(sro)
# saveRDS(sro.imp@assays$MAGIC_RNA@data, file = "Seurat/merged-5A/imputed-expr.rds")

sro <- readRDS("Seurat/merged-4I/SRO.rds")
sro.imp <- magic(sro)
saveRDS(sro.imp@assays$MAGIC_RNA@data, file = "Seurat/merged-4I/imputed-expr.rds")
