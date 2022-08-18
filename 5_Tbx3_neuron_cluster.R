library(Seurat)
library(patchwork)
library(dplyr)

options(future.globals.maxSize= 891289600000)
rm(list = ls())
load("/public-supool/home/zhchen/Lab_share/public_source/cellcycle_gene/mouse_cell_circle_gene.Rdata")
setwd("/public-supool/home/zhchen/XS/5_subset_neuron")
Tbx3_neuron <- readRDS("/public-supool/home/zhchen/XS/4_rm_doublet/Tbx3_neuron.rds")


# text <- Tbx3_neuron[, sample(colnames(Tbx3_neuron), size =500, replace=F)]


if(T){
	Tbx3_neuron <- NormalizeData(Tbx3_neuron,normalization.method = "LogNormalize", scale.factor = 10000) %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) %>%
    ScaleData(vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Tbx3_neuron))
	Tbx3_neuron <- RunPCA(Tbx3_neuron,features = VariableFeatures(Tbx3_neuron),npcs = 50) %>%
    FindNeighbors(reduction = "pca",dims = 1:50,nn.eps = 0.5) %>%
    FindClusters(resolution = 0.01, n.start = 10,random.seed = 12580) %>%  ## 可选
    RunUMAP(dims = 1:10,seed.use = 12580) %>%
    RunTSNE(dims = 1:10,seed.use = 12580)
}
### 先保存
saveRDS(Tbx3_neuron,file ="./Tbx3_neuron_final.rds")

### 聚类
resulation <- c(0.1,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5)
Marker <- c()
p <- c()
setwd("/public-supool/home/zhchen/XS/5_subset_neuron/marker")
for(i in 1:length(resulation)){
  Tbx3_neuron <- FindNeighbors(Tbx3_neuron, reduction = "pca", dims = 1:30, nn.eps = 0.5) %>% 
  FindClusters(resolution =resulation[i] , n.start = 10,random.seed = 12580)
  Marker <- FindAllMarkers(Tbx3_neuron, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.25)
  write.csv(Marker,file = paste0( "Allmarker_re",resulation[i], ".csv"))
  #### tsne 单张图就 10，7;双排为18;
  pdf(file=paste0("tsne_split_stim_re",resulation[i], ".pdf"),width = 18, height = 7)
  p <- DimPlot(Tbx3_neuron, reduction = "tsne",label=TRUE,split.by = "group")
  print(p)
  dev.off()
  pdf(file=paste0("tsne_group_stim_re",resulation[i], ".pdf"),width = 10, height = 7)
  p <- DimPlot(Tbx3_neuron, reduction = "tsne",label=TRUE,group = "group")
  print(p)
  dev.off()
  pdf(file=paste0("tsne_re",resulation[i], ".pdf"),width = 10, height = 7)
  p <- DimPlot(Tbx3_neuron, reduction = "tsne",label=TRUE)
  print(p)
  dev.off()
  
  ### UMAP ####
  
  pdf(file=paste0("umap_split_state_re",resulation[i], ".pdf"),width = 18, height = 7)
  p <- DimPlot(Tbx3_neuron, reduction = "umap",label=TRUE,split.by = "group")
  print(p)
  dev.off()
  pdf(file=paste0("umap_group_state_re",resulation[i], ".pdf"),width = 10, height = 7)
  p <- DimPlot(Tbx3_neuron, reduction = "umap",label=TRUE,group = "group")
  print(p)
  dev.off()
  pdf(file=paste0("umap_re",resulation[i], ".pdf"),width = 10, height = 7)
  p <- DimPlot(Tbx3_neuron, reduction = "umap",label=TRUE)
  print(p)
  dev.off()
}


### 再保存
setwd("/public-supool/home/zhchen/XS/5_subset_neuron")
saveRDS(Tbx3_neuron,file ="./Tbx3_neuron_final.rds")
rm(list = ls())
