# scRNAseq.Cancer.Rpackage

Workflow for single-cell RNA-seq analysis (especially focused on breast cancer dataset)

## Installation

> devtools::install_github('CB-postech/scRNAseq.Cancer.Rpackage')

## Pre-requisite R packages you should install before use this package
```
Seurat (>= 4.0.0), scater, scRNAseq, genefu, harmony, dplyr, KernelKnn, tidyverse, DESeq2, rlist, data.table    
```

## Scripts

#### scRNA Analysis Scripts
*Integrative_clustering.R*   
-make seurat object from raw umi counts    
-run harmony for batch correction and re-clustering    

#### Projection of external transcriptome data to dimension reduction plot (e.g., UMAP or PCA) processed with your own scRNA-seq data
*Projection_external_transcriptome_dataset.R*    
-make external transcriptome data to pseudobulk form    
-projection using KernelKNN method

#### Pseudobulk DEG ananlysis 
*Pseudobulk_DEG_analysis_in_scRNAseq.R*    
-make scRNA-seq data to pseudobulk form (aggregate by batch or donor)    
-using DESeq2 to retrieve differentially expressed genes from user-defined contrast design    

#### PAM50 classification in scRNA-seq to reveal molecular subtypes of breast cancer cells
*scRNAseq_PAM50_prediction.R*    
-prepare normalized counts matrix
-using genefu package to predict molecular subtypes of breast cancer cells

#### Cell-of-Origin prediction in scRNA-seq 
*seurat2CellofOrigin.R*    
-prepare normalized counts matrix    
-using scHCL package to predict cell-of-origins of your own single-cell rna-seq data    
-refine the results from scHCL by considering the tissue origin where you resected    

#### Estimation of TF-regulon activity
*seurat2dorothea.R*    
-user should pre-define confidence level and species to set regulon    
-from seurat object, calculate the activity by using DoRothEA package    

