setwd("~/0-workspace/CCR7_DC/oral-tolerance-Littman/")

library(Seurat)
library(Rmagic)
library(reticulate)
py_require("magic-impute")

sro <- readRDS("Seurat/merged/SRO.rds")

sro.imp <- magic(sro)
saveRDS(sro.imp@assays$MAGIC_RNA@data, file = "Seurat/merged/imputed-expr.rds")
