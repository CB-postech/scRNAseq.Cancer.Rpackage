# scRNAseq.Cancer.Rpackage

Workflow for single-cell RNA-seq analysis (especially focused on breast cancer dataset)

## Installation

devtools::install_github('CB-postech/scRNAseq.Cancer.Rpackage')


## Scripts

#### scRNA Analysis Scripts
*Integrative_clustering.R*   
-make seurat object from raw umi counts    
-run harmony for batch correction and re-clustering    

#### Projection of external transcriptome data to dimension reduction plot (e.g., UMAP or PCA) processed with your own scRNA-seq data
*Projection_external_transcriptome_dataset.R*    
-make external transcriptome data to pseudobulk form    
-projection using KernelKNN method

