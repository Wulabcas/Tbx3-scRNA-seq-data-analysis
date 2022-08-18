  library(Seurat)
  library(patchwork)
  library(dplyr)
library(DoubletFinder)

rm(list = ls())
Tbx3_merge <- readRDS("/public-supool/home/zhchen/XS/3_cluster/Tbx3_merge_0911.rds")

if(T){
  library(DoubletFinder)
  sweep.res.list_kidney <- c()
  sweep.stats_kidney <- c()
  bcmvn_kidney <- c()
  mpK <- c()
  annotations <- c() ## 
  homotypic.prop <- c()
  nExp_poi <- c()
  nExp_poi.adj <- c()
  seurat_filterDouble <- c()
}

Idents(Tbx3_merge) <- Tbx3_merge$RNA_snn_res.0.05
Seurat_object <- Tbx3_merge ## 
rm(Tbx3_merge)
if(T){
  sweep.res.list_kidney <- paramSweep_v3(Seurat_object, PCs = 1:40, sct = T)
  sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
  bcmvn_kidney <- find.pK(sweep.stats_kidney)
  mpK<-as.numeric(as.vector(bcmvn_kidney$pK[which.max(bcmvn_kidney$BCmetric)]))
  
  #  nExp
  annotations <- Idents(Seurat_object) ## 
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.05*ncol(Seurat_object@assays$RNA@data)) ## Doublet ,BGI 0.05
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  #  Doublet
  seurat_filterDouble <- doubletFinder_v3(Seurat_object, PCs = 1:40, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  
}
#### 
saveRDS(seurat_filterDouble,file = "seurat_filterDouble_0911.rds")
rm(list = ls())

