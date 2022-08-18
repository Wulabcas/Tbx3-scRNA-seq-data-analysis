
#### correlation
wt_10x <- as.data.frame(tbx3_10x_WT_ave$RNA)
colnames(wt_10x) <- c("x3","x12","x13","x14","x16","x18")
wt_10x$gene <- rownames(wt_10x)
head(wt_10x)


wt_BGI <- as.data.frame(tbx3_BGI_WT_ave$RNA)
colnames(wt_BGI) <- rev(c("d9","d14","d17","d41","d42","d46","d53","d58","d60","d65","d66","d67"))
wt_BGI$gene <- rownames(wt_BGI)
head(wt_BGI)

cor_data <- merge(wt_10x,wt_BGI,by="gene")
marker <- read.csv("D:/personal user/CZH/Analysis/For XS/出图_20210927/BGI_WT_MARKER.csv")
top100 <- marker %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)

top100 <- top100$gene
top100  <- top100 [!duplicated(top100 )]
rownames(cor_data) <- cor_data$gene
cor_data2 <- cor_data[top100,]
cor_data2 <- na.omit(cor_data2)
cor_data2 <- cor_data2[,-1]

cor_caculate <- cor(log2(cor_data2+1),method = 'pearson')
cor_caculate <- cor_caculate[-c(1:6),]
cor_caculate <- cor_caculate[,c(1:6)]

pheatmap::pheatmap(cor_caculate,border="white")


### Sankey Plot
library(ggalluvial)


head(cor_caculate)
cor_data3 <- melt(cor_caculate)

cor_data3

cor_data4 <- cor_data3[c(9,23,27,48,43,46,44,42,52,65),]
cor_data4$factor <- factor(cor_data4$value)
## factor
ggplot(data = cor_data4,
       aes(axis1 = Var1, axis2 = Var2,
           y = value)) +
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_x_discrete(limits = c("Var1", "Var2"), expand = c(.2, .05)) +
  xlab("Demographic") +
  geom_alluvium(aes(fill = factor)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_manual(name = "Values",values = c("#FFEBEE","#FFCDD2","#EF9A9A","#E57373","#EF5350","#F44336","#E53935","#D32F2F","#C62828","#B71C1C"))+
  theme_minimal() +
  ggtitle("Cor between 10x and BGI data")


ggplot(data = cor_data4,
       aes(axis1 = Var1, axis2 = Var2,
           y = value)) +
  scale_x_discrete(limits = c("Var1", "Var2"), expand = c(.2, .05)) +
  xlab("Demographic") +
  geom_alluvium(aes(fill = factor)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() +
  ggtitle("Cor between 10x and BGI data")
  

cell_id <- c(9,14,17,41,46,60,58,65,53,42,65,53,42)
DefaultAssay(tbx3) <- "RNA"
marker <- c()
p <- c()


#### Feature plot
for (i in 1:length(cell_id)) {
  marker <- read.delim(paste0("D:/personal user/CZH/Analysis/For XS/For_XS_analysis/featurePlot/",cell_id[i],".txt"),sep = ",",header = F)
  for (j in 1:length(marker)) {
    p <- FeaturePlot(tbx3,features = marker[j],label = T)
    png(paste0("D:/personal user/CZH/Analysis/For XS/For_XS_analysis/featurePlot/",cell_id[i],"/",marker[j],".png"))
    print(p)
    dev.off()
  }


}



#### Violin plot
for (i in 1:length(cell_id)) {
  marker <- read.delim(paste0("D:/personal user/CZH/Analysis/For XS/For_XS_analysis/featurePlot/",cell_id[i],".txt"),sep = ",",header = F)
  for (j in 1:length(marker)) {
    p <- VlnPlot(tbx3,features = marker[j],slot = "data",idents = paste0(cell_id2) )
    png(paste0("D:/personal user/CZH/Analysis/For XS/For_XS_analysis/featurePlot/",cell_id[i],"/",marker[j],"_vin_.png"))
    print(p)
    dev.off()
  }
  
  
}






