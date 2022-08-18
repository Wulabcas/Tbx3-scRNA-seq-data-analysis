setwd("/public-supool/home/zhchen/XS/merge")
rm(list = ls())
library(Seurat)
library(patchwork)
library(dplyr)

##### load data############
path_way <- c("/public-supool/home/zhchen/XS/data/Mouse_Hyp_male_Ctrl_0520_12_20-1",
             "/public-supool/home/zhchen/XS/data/Mouse_Hyp_male_Ctrl_0520_12_20-2",
             "/public-supool/home/zhchen/XS/data/Mouse_Hyp_male_Ctrl_0520_12_20-3",
             "/public-supool/home/zhchen/XS/data/Mouse_Hyp_female_Ctrl_0520_6-1",
             "/public-supool/home/zhchen/XS/data/Mouse_Hyp_female_Ctrl_0520_6-2",
             "/public-supool/home/zhchen/XS/data/Mouse_Hyp_male_CKO_0520_15_13-1",
             "/public-supool/home/zhchen/XS/data/Mouse_Hyp_male_CKO_0520_15_13-2",
             "/public-supool/home/zhchen/XS/data/Mouse_Hyp_male_CKO_0520_15_13-3",
             "/public-supool/home/zhchen/XS/data/Mouse_Hyp_female_CKO_0520_4-1",
             "/public-supool/home/zhchen/XS/data/Mouse_Hyp_female_CKO_0520_4-2") ## path

samble <- c("C1_male","C2_male","C3_male","C4_fem","C5_fem","K1_male","K2_male","K3_male","K4_fem","K5_fem") ## sample names
condition <- c("WT","WT","WT","WT","WT","CKO","CKO","CKO","CKO","CKO")
gender <- c("male","male","male","fem","fem","male","male","male","fem","fem")
data_xs <- list() ## 
for (i in 1:length(samble)){
	data_xs[[i]] <- Read10X(data.dir = path_way[i],gene.column = 1) %>% CreateSeuratObject(project = paste0(samble[i], "_", "data"), min.cells = 3, min.features = 200) 
    names(data_xs)[i] <- samble[i] ## 
	data_xs[[i]]$group <- paste0(condition[i])
	data_xs[[i]]$gender <- paste0(gender[i])
}


### wt 和 ko 合并
WT <- merge(data_xs[[1]],y = c(data_xs[[2]],data_xs[[3]],data_xs[[4]],data_xs[[5]]), add.cell.ids = c("C1_male","C2_male","C3_male","C4_fem","C5_fem"), project = "Contral")
KO <- merge(data_xs[[6]],y = c(data_xs[[7]],data_xs[[8]],data_xs[[9]],data_xs[[10]]),   add.cell.ids = c("K1_male","K2_male","K3_male","K4_fem","K5_fem"), project = "Mut")
data_xs_merge <- list(WT,KO)  ## 

Tbx3 <- merge(WT,y=KO,add.cell.ids=c("WT","CKO"),project = "Tbx3")

saveRDS(Tbx3,file ="./Tbx3.rds")

rm(list = ls())
