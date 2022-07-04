#' Projection of external dataset

# Based source code offered by  follow publication
# Jaitin, Diego Adhemar, et al. "Lipid-associated macrophages control metabolic homeostasis in a Trem2-dependent manner." Cell 178.3 (2019): 686-698.

#' @title Projecting external dataset (transcriptome data) to my own dataset
#' @description Projection the public/external dataset using KernelKNN (with pearson correlation option) method to your own dataset
#' @param seurat.obj RNAseq.query_dataset
#' @return proj.pca.location: Projected coordinates of public/external dataset on the dimension reduction plot from your own dataset
#' @details RNAseq.query_dataset: external dataset, seurat.obj: my own dataset
#' @name projection2seurat
#' @export

library(Seurat)
library(scater)
library(scRNAseq)
library(dplyr)
library(harmony)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(magrittr)
library(Matrix)
library(Matrix.utils)
library(DESeq2)
library(tidyverse)
library(KernelKnn)

# RNAseq.query_dataset
# seurat.obj

projection2seurat <- function(seurat.obj, RNAseq.query_dataset){

  proj.geneset=intersect(VariableFeatures(seurat.obj), rownames(RNAseq.query_dataset))

  scRNAseq.reference_dataset=seurat.obj@assays$RNA@data[proj.geneset,]
  RNAseq.query_dataset=RNAseq.query_dataset[proj.geneset,]

  set.seed(12345)
  indexesN <- knn.index.dist(data=t(scRNAseq.reference_dataset), TEST_data = t(RNAseq.query_dataset), k = 20, method = "pearson_correlation")
  iN2 <- indexesN$test_knn_idx
  rownames(iN2) <- colnames(B)
  iN3 <- apply(iN2, 2, function(x) colnames(A)[x])
  idxx <- apply(iN3, 2, function(x) seurat.obj@reductions$pca@cell.embeddings[x,1])
  idxy <- apply(iN3, 2, function(x) seurat.obj@reductions$pca@cell.embeddings[x,2])
  prjx <- rowMeans(idxx)
  prjy <- rowMeans(idxy)
  names(prjx) <- rownames(iN2)
  names(prjy) <- rownames(iN2)

  proj.pca.location=data.frame(PC_1=prjx, PC_2=prjy, sampleID=rownames(iN2))

  return(proj.pca.location)

}
