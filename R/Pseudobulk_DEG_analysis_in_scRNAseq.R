#' Pseudobulk differentially expressed gene (DEG) analysis
# Reference "https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html"

#' @title Pseudobulk DEG analysis with Seurat object
#' @description Aggregation of single-cell RNA-seq data by batch (or specimen)
#' @param seurat.obj metadata
#' @return res_tbl # Data frame of DEG results using DESeq2 package
#' @details TBA
#' @name seurat2pseudobulk_DEA
#' @export

library(Seurat)
library(scater)
library(scRNAseq)
library(rlist)
library(DESeq2)
library(dplyr)
library(tidyverse)

seurat2pseudobulk_DEA <- function(seurat.obj, metadata) {

  pb_counts=seurat.obj@assays$RNA@counts
  pb_counts=pb_counts[rowSums(pb_counts)>0,]
  df=as.data.frame(t(as.matrix(pb_counts)))
  df$donors=seurat.obj$specimenID[rownames(df)]

  pb_counts.aggr=aggregate(. ~ donors, data=df, FUN=sum)
  pb_counts.aggr=pb_counts.aggr %>% remove_rownames %>% column_to_rownames(var="donors")
  pb_counts.aggr=t(pb_counts.aggr)

  # metadata = data.frame(row.names = colnames(pb_counts.aggr), donors = colnames(pb_counts.aggr))
  # metadata$tumor.subtype=factor(metadata$tumor.subtype) # depend on the study design
  # metadata$donors=factor(metadata$donors)

  library(DESeq2)
  dds <- DESeqDataSetFromMatrix(pb_counts.aggr, colData = metadata, design = ~ tumor.subtype)
  rld <- rlog(dds)

  ### Run DESeq2 differential expression analysis
  dds <- DESeq(dds)
  contrast <- c("tumor.subtype", "Group_1", "Group_2")
  res <- results(dds, contrast = contrast)
  res_tbl <- res %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble() %>% na.omit()

  return (res_tbl)
}
