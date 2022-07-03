#' Integrative clustering.R
# Reference-1: https://satijalab.org/seurat/
# Reference-2: https://github.com/immunogenomics/harmony

#' @title Integrative clustering
#' @description Clustering and Dimension reduction (Seurat) after batch correction with harmony
#' @param orig.seurat nHVGs dims harmony.dims batch.term
#' @return harmony.seurat
#' @details TBA
#' @name seurat2harmony
#' @export


library(scater)
library(Seurat)
library(harmony)
library(scater)
library(SingleCellExperiment)

seurat2harmony <- function(orig.seurat,nHVGs,dims,harmony.dims,batch.term){

  orig.seurat=orig.seurat %>% NormalizeData() %>% FindVariableFeatures(nfeatures=nHVGs) %>%
    ScaleData() %>% RunPCA(npcs=dims) %>%
    FindNeighbors(dims=1:dims) %>% FindClusters() %>%
    RunUMAP(dims = 1:dims)

  harmony.seurat=orig.seurat
  harmony.seurat=RunHarmony(harmony.seurat,batch.term)
  harmony.seurat <- harmony.seurat %>%
    RunUMAP(reduction = "harmony", dims = 1:harmony.dims) %>%
    RunTSNE(reduction = "harmony", dims = 1:harmony.dims) %>%
    FindNeighbors(reduction = "harmony", dims = 1:harmony.dims) %>%
    FindClusters() %>% identity()

  return(harmony.seurat)
}
