library(scRNAseq.Cancer.Rpackage)

# orig.seurat: Seurat object created from single-cell rna-seq UMI counts (this code is focused on 10x genomics platform.)
# nHVGs: number of highly variable genes which is used for initial clustering
# dims: number of dimensions from PCA which is used for initial clustering
# harmony.dims: number of dimensions of harmony corrected dims to re-clustering
# batch.term: this attributes should included in the metadata of orig.seurat object

harmony.seurat=seurat2harmony(orig.seurat, nHVGs=2000, dims=30, harmony.dims=20,batch.term=batch.term)
