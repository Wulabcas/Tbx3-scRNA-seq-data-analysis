setwd("/public-supool/home/zhchen/XS/3_cluster") 
options(future.globals.maxSize= 891289600000)
rm(list = ls())
if(T){
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(stringr)
  library(ggplot2)
  library(tidyverse)
} ##load package

Tbx3_merge <- readRDS("/public-supool/home/zhchen/XS/2_preproccess/Tbx3_merge_0910.rds")


##QC1 
 pdf(file=paste0("./Tbx3_merge_0911", "_QC1", ".pdf"),width = 18, height = 7) 
  p <- VlnPlot(Tbx3_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.01) 
  print(p)
  dev.off()
  ## QC2
  plot1 <- FeatureScatter(Tbx3_merge, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(Tbx3_merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  pdf(file=paste0("./Tbx3_merge", "_QC2", ".pdf"),width = 18, height = 7) 
  print(plot1 + plot2)
  dev.off()



Allmarkers <- list()
p <- c()
for (i in c(0.05,0.1,0.3,0.5)) {
    Tbx3_merge <- RunPCA(Tbx3_merge,features = VariableFeatures(Tbx3_merge),npcs = 50) %>%
    FindNeighbors(reduction = "pca",dims = 1:50,nn.eps = 0.5) %>%
    FindClusters(resolution = i, n.start = 10,random.seed = 12580) %>% 
    RunUMAP(dims = 1:10,seed.use = 12580) %>%
    RunTSNE(dims = 1:10,seed.use = 12580)
  ######## plot ######
  #### tsne 10，7
  #### umap
  pdf(file=paste0("./Tbx3_merge_0911_",i, "_umap", ".pdf"),width = 10, height = 7)
  p <- DimPlot(Tbx3_merge, reduction = "umap",label =T)
  print(p)
  dev.off() 
  ##### tsne ###
  pdf(file=paste0("./Tbx3_merge_0911_",i, "_tsne", ".pdf"),width = 10, height = 7)
  p <- DimPlot(Tbx3_merge, reduction = "tsne",label =T)
  print(p)
  dev.off()
  ### marker
}

saveRDS(Tbx3_merge,file = "./Tbx3_merge_0911.rds")
rm(list = ls())
