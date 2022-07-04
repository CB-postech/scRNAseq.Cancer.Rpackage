#' Pipeline for regulon analysis

# Reference: Garcia-Alonso L, Holland C, Ibrahim M, Turei D, Saez-Rodriguez J (2019). “Benchmark and integration of resources for the estimation of human transcription factor activities.” Genome Research. doi: 10.1101/gr.240663.118.
# Reference: https://bioconductor.org/packages/release/data/experiment/html/dorothea.html

#' @title TF-Regulon analysis using DoRothEA package
#' @description Estimation of gene regulatory network (TF ~ target genes)
#' @param seurat.obj regulon
#' @return summarized_viper_scores_df
#' @details
#' @name seurat2dorothea
#' @export

library(Seurat)
library(dorothea)
library(dplyr)
library(tibble)
library(pheatmap)
library(tidyr)
library(viper)

### User should pre-define the regulon to the dataset-specific what you have
# dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
# regulon <- dorothea_regulon_human %>% dplyr::filter(confidence %in% c("A","B","C")) # user can change the confidence level in dorothea

seurat2dorothea <- function(seurat.obj, regulon){

  seurat.obj <- run_viper(seurat.obj, regulon)
  DefaultAssay(object = seurat.obj) <- "dorothea"
  seurat.obj <- ScaleData(seurat.obj)

  seurat.obj$celltypexsubtype=paste(seurat.obj$celltype, seurat.obj$tumor.subtype, sep = "_")

  viper_scores_df <- GetAssayData(seurat.obj, slot = "scale.data", assay = "dorothea") %>% data.frame(check.names = F) %>% t()
  CellsClusters <- data.frame(cell = colnames(seuratset), cell_groups = seuratset[["celltypexsubtype"]][[1]], check.names = F)

  viper_scores_clusters <- viper_scores_df %>% data.frame() %>% rownames_to_column("cell") %>% gather(tf, activity, -cell) %>% inner_join(CellsClusters)
  summarized_viper_scores <- viper_scores_clusters %>% group_by(tf, cell_groups) %>% summarise(avg = mean(activity), std = sd(activity))
  highly_variable_tfs <- summarized_viper_scores %>% group_by(tf) %>% mutate(var = var(avg))  %>% ungroup() %>% distinct(tf)

  summarized_viper_scores_df <- summarized_viper_scores %>% semi_join(highly_variable_tfs, by = "tf") %>%
    dplyr::select(-std) %>% spread(tf, avg) %>% data.frame(row.names = 1, check.names = FALSE)

  return (summarized_viper_scores_df)
}
