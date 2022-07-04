#' Pipeline for predict cell-of-origin

# Reference: https://github.com/ggjlab/scHCL

#' @title cell-of-origin prediction and summary
#' @description cell-of-origin using scHCL package
#' @param seurat.obj scHCL.coef.threshold scHCL_prediction_geneset Target.tissue
#' @return scHCL_summary
#' @details
#' @name seurat2CellofOrigin
#' @export


library(Seurat)
library(scater)
library(dplyr)
library(scHCL)

seurat2CellofOrigin <- function(scdata, scHCL.coef.threshold, scHCL_prediction_geneset, Target.tissue, numbers_plot=3){

  scdata = seurat.obj@assays$RNA@counts
  scdata = scdata[rowSums(scdata)>0,]
  scHCL_prediction_geneset = intersect(rownames(scHCL::ref.expr), rownames(scdata))

  scdata = scdata[scHCL_prediction_geneset,]
  ref.expr = ref.expr[scHCL_prediction_geneset,]

  tst.expr<-scdata
  tst.expr<-as.matrix(t(t(tst.expr)/colSums(tst.expr))*100000) # Normalization (scale factor: 10^5)
  tst.expr<-log(tst.expr+1)
  cors <- cor(log(ref.expr+1),tst.expr)

  cors_index <- apply(cors,2,gettissue,numbers_plot)
  cors_index <- sort(unique(as.integer(cors_index)))
  scblast.result <- apply(cors,2,function(x) rownames(cors)[which.max(x)])
  cors_in = cors[cors_index,]
  colnames(cors_in)<-colnames(scdata)
  cors_out = reshape2::melt(cors_in)[c(2,1,3)]
  colnames(cors_out)<- c("Cell","Cell type","Score")
  cors_out <- as.data.frame(cors_out %>% group_by(Cell) %>% top_n(n=numbers_plot,wt=Score))

  result <- list()
  cors[which(is.na(cors))]<-0
  result[["cors_matrix"]] <- cors
  result[['top_cors']]<-numbers_plot
  result[['scHCL']]<-scblast.result
  result[['scHCL_probility']]<-cors_out

  # Revise the results from scHCL
  # When your target tissue is not ranked at first position, but positioned at 2nd or 3rd, then replace 1st ranked value to your target tissue

  assigned_celltype=c()
  coef.threshold=scHCL.coef.threshold
  scHCL_probility_matrix=result$scHCL_probility
  scHCL_probility_matrix$Score[scHCL_probility_matrix$Score <= coef.threshold]=NA
  scHCL_probility_matrix=na.omit(scHCL_probility_matrix)

  for (cellbarcode in colnames(seurat.obj)){

    if (!(cellbarcode %in% as.vector(scHCL_probility_matrix$Cell))){
      assigned_celltype_tmp <- "Unassigned"
    } else {
      barcode_scHCL.res=scHCL_probility_matrix[scHCL_probility_matrix$Cell==cellbarcode,]
      if (grepl("Breast",barcode_scHCL.res[,"Cell type"])){
        assigned_celltype_tmp=as.vector(barcode_scHCL.res[which.max(barcode_scHCL.res[grepl(Target.tissue,barcode_scHCL.res[,"Cell type"]),"Score"]), "Cell type"])
      } else {
        assigned_celltype_tmp=as.vector(barcode_scHCL.res[which.max(barcode_scHCL.res$Score),"Cell type"])
      }
    }
    assigned_celltype=c(assigned_celltype,assigned_celltype_tmp)
  }

  scHCL_summary=data.frame(row.names = colnames(mpe.cancercell), assigned_celltype=assigned_celltype)

  return(scHCL_summary)
}
