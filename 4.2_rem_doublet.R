library(DoubletFinder)
library(Seurat)
library(patchwork)
library(dplyr)
rm(list = ls())
load("/public-supool/home/zhchen/Lab_share/public_source/cellcycle_gene/mouse_cell_circle_gene.Rdata")


text <- readRDS("/public-supool/home/zhchen/XS/4_rm_doublet/seurat_filterDouble_0911.rds")


doublet_cell_id <- colnames(text)[which(text@meta.data$DF.classifications_0.25_0.26_4068 == "Doublet")]
Tbx3_remo_double <- subset(text,cells = doublet_cell_id,invert = TRUE) ## subset single cell 
rm(text)

### 重新pipline 测试ok
if(T){
  Tbx3_remo_double <- NormalizeData(Tbx3_remo_double,normalization.method = "LogNormalize", scale.factor = 10000) %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) %>%
    ScaleData(vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Tbx3_remo_double))
  Tbx3_remo_double <- RunPCA(Tbx3_remo_double,features = VariableFeatures(Tbx3_remo_double),npcs = 50) %>%
    FindNeighbors(reduction = "pca",dims = 1:50,nn.eps = 0.5) %>%
    FindClusters(resolution = 0.01, n.start = 10,random.seed = 12580) %>%  ## 
    RunUMAP(dims = 1:10,seed.use = 12580) %>%
    RunTSNE(dims = 1:10,seed.use = 12580)
}

saveRDS(Tbx3_remo_double,file="./Tbx3_remo_double_0911.rds")

### Find cluster
resulation <- c(0.01,0.02,0.03,0.05,0.1)
Marker <- c()
p <- c()
for(i in 1:length(resulation)){
  Tbx3_remo_double <- FindNeighbors(Tbx3_remo_double, reduction = "pca", dims = 1:30, nn.eps = 0.5) %>% 
    FindClusters(resolution =resulation[i] , n.start = 10,random.seed = 12580) %>%
    RunUMAP(reduction = "pca", dims = 1:30,seed.use = 12580) %>%
    RunTSNE(reduction = "pca", dims = 1:30,seed.use = 12580)
  Marker <- FindAllMarkers(Tbx3_remo_double, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.25)
  write.csv(Marker,file = paste0( "Allmarker_re_0911",resulation[i], ".csv"))
  #### tsne 单张图就 10，7;双排为18，1
  pdf(file=paste0("tsne_split_stim_re_0911_",resulation[i], ".pdf"),width = 18, height = 7)
  p <- DimPlot(Tbx3_remo_double, reduction = "tsne",label=TRUE,split.by = "group")
  print(p)
  dev.off()
  pdf(file=paste0("tsne_group_stim_re_0911_",resulation[i], ".pdf"),width = 10, height = 7)
  p <- DimPlot(Tbx3_remo_double, reduction = "tsne",label=TRUE,group = "group")
  print(p)
  dev.off()
  pdf(file=paste0("tsne_re_0911_",resulation[i], ".pdf"),width = 10, height = 7)
  p <- DimPlot(Tbx3_remo_double, reduction = "tsne",label=TRUE)
  print(p)
  dev.off()
  
  ### UMAP ####
  
  pdf(file=paste0("umap_split_state_re_0911_",resulation[i], ".pdf"),width = 18, height = 7)
  p <- DimPlot(Tbx3_remo_double, reduction = "umap",label=TRUE,split.by = "group")
  print(p)
  dev.off()
  pdf(file=paste0("umap_group_state_re_0911_",resulation[i], ".pdf"),width = 10, height = 7)
  p <- DimPlot(Tbx3_remo_double, reduction = "umap",label=TRUE,group = "group")
  print(p)
  dev.off()
  pdf(file=paste0("umap_re_0911_",resulation[i], ".pdf"),width = 10, height = 7)
  p <- DimPlot(Tbx3_remo_double, reduction = "umap",label=TRUE)
  print(p)
  dev.off()
  
  
}

saveRDS(Tbx3_remo_double,file ="./Tbx3_remo_double_0911.rds")
rm(list = ls())
