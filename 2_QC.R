
########
options(future.globals.maxSize= 891289600000)
rm(list = ls())
if(T){
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(stringr)
  library(ggplot2)
  library(tidyverse)
} ##加载包

setwd("/public-supool/home/zhchen/XS/2_preproccess")
load("/public-supool/home/zhchen/Lab_share/public_source/cellcycle_gene/mouse_cell_circle_gene.Rdata")

### 走流程
Tbx3_merge <- readRDS("/public-supool/home/zhchen/XS/merge/Tbx3.rds")
Tbx3_merge <- PercentageFeatureSet(Tbx3_merge,pattern = "^mt-", col.name = "percent.mt")%>% subset(subset = nFeature_RNA > 800 & nFeature_RNA < 8000 & percent.mt < 10) %>% NormalizeData(verbose = FALSE) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 3000)
Tbx3_merge <- CellCycleScoring(Tbx3_merge, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Tbx3_merge <- ScaleData(Tbx3_merge,features = rownames(Tbx3_merge))
### 看是否需要去除cell cycle
Tbx3_merge <- RunPCA(Tbx3_merge, features = c(s.genes, g2m.genes))
pdf(file=paste0("Tbx3_merge_cc2","_pca", ".pdf"),width = 10, height = 7)
p <- DimPlot(Tbx3_merge, reduction = "pca")
print(p)
dev.off()

rm(list = ls())

