# devtools::install_github('CB-postech/scRNAseq.Cancer.Rpackage')
library(scRNAseq.Cancer.Rpackage)

### Initiual Clustering followed by (batch correction and re-clustering)

# orig.seurat: Seurat object created from single-cell rna-seq UMI counts (this code is focused on 10x genomics platform.)
# nHVGs: number of highly variable genes which is used for initial clustering
# dims: number of dimensions from PCA which is used for initial clustering
# harmony.dims: number of dimensions of harmony corrected dims to re-clustering
# batch.term: this attributes should included in the metadata of orig.seurat object
harmony.seurat = seurat2harmony(orig.seurat, nHVGs=2000, dims=30, harmony.dims=20,batch.term=batch.term)
  
### Projection public(external) dataset    
# bulkRNAseq_data: rows: genes, columns: samples    
proj.result = projection2seurat(harmony.seurat, bulkRNAseq_data)

### Pseudobulk differentially expressed gene (DEG) analysis per cell types
target.celltype <- "Cancer.cell"
Idents(harmony.seurat) <- "celltype.attribute" # celltype.attribute indicates column name in meta.data of harmony.seurat contains cell.type label information
celltype_specific.seurat = subset(harmony.seurat, cells = target.celltype)
DEG_results = seurat2pseudobulk_DEA(celltype_specific.seurat, metadata) # metadata for DESeq2 analysis (row.names = colnames of aggregated counts matrix from scRNA-seq data, contrast.group = contrast.group)

