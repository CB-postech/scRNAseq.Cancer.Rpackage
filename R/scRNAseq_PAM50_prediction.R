#' PAM50 prediction

# Based on genefu package: https://www.bioconductor.org/packages/release/bioc/html/genefu.html
# Deena M.A. Gendoo, Natchar Ratanasirigulchai, Markus S. Schroeder, Laia Pare, Joel S Parker, Aleix Prat, Benjamin Haibe-Kains (2022). genefu: Computation of Gene Expression-Based Signatures in Breast Cancer. R package version 2.28.0, http://www.pmgenomics.ca/bhklab/software/genefu.

#' @title PAM50 prediction with Seurat object
#' @description PAM50 classification with single-cell RNA-seq data
#' @param seurat.obj
#' @return PAM50xClusters
#' @details TBA
#' @name seurat2PAM50
#' @export

library(Seurat)
library(scater)
library(scRNAseq)
library(rlist)
library(genefu)

seurat2PAM50 <-function(seurat.obj) {

  lognorm.counts=as.matrix(seurat.obj@assays$RNA@data)
  lognorm.counts=lognorm.counts[rowSums(lognorm.counts)>0,]

  annot=data.frame(Gene.Symbol=rownames(lognorm.counts))
  pam50_res.RNA=molecular.subtyping(sbt.model=c("pam50"), t(lognorm.counts), annot)
  seurat.obj$pam50.prediction=pam50_res.RNA$subtype[colnames(seurat.obj)]
  PAM50xClusters=table(seurat.obj$pam50.prediction, seurat.obj$seurat_clusters)

  return (PAM50xClusters)
}
