setwd("~/0-workspace/CCR7_DC/oral-tolerance-Gardner/")

library(Seurat)
library(Rmagic)
library(reticulate)
py_require("magic-impute")

sro.a <- readRDS("Seurat/adult/SRO.rds")

sro.imp.a <- magic(sro.a)
saveRDS(sro.imp.a@assays$MAGIC_RNA@data, file = "Seurat/adult/imputed-expr.rds")


sro.e <- readRDS("Seurat/early/SRO.rds")

sro.imp.e <- magic(sro.e)
saveRDS(sro.imp.e@assays$MAGIC_RNA@data, file = "Seurat/early/imputed-expr.rds")
