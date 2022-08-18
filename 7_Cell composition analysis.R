
if(T){
  library(tidyverse)
  library(magrittr)
  library(conos)
  library(cacoa)
  library(Matrix)
  library(ggpubr)
  library(cowplot)
  library(pagoda2)
  library(reshape2)
  library(RColorBrewer)
  library(tidyr)
  library(reshape2)
}


load("E:/wu_lab/For XS/20220709/Tbx3_metadata.Rdata")
a <- as.matrix(table(metadata$nchip,metadata$cell))
data.frame(a)
a <- as.data.frame(a)
a <- data.frame(spread(a, Var2, Freq))
write.csv(a,file = "D:/annaconda_install/Lib/site-packages/sccoda/datasets/text.csv")

colnames(a) <- c("Var1",paste0("d",seq(1:67)))
a$chip <- rownames(a)
a$group <- c(rep("WT",3),rep("CKO",3))
write.csv(a[,c("Var1","group")],file = "D:/annaconda_install/Lib/site-packages/sccoda/datasets/text1.csv")

metadata$cell <- str_replace(metadata$cell,"N","d")

tbx3_cnts <- table(metadata$nchip,metadata$cell)
tbx3_group <- c(FALSE,FALSE,FALSE,TRUE,TRUE,TRUE)
names(tbx3_group) <- rownames(tbx3_cnts)
tbx3_sample_group <- c(rep("K",3),rep("C",3))
names(tbx3_sample_group) <- rownames(tbx3_cnts)

tbx3_freqs <- tbx3_cnts %>% {. / rowSums(.)}
groups.name <- c( 'K','C')


res <- cacoa:::runCoda(tbx3_cnts, tbx3_group, n.seed=239)
dfs <- cacoa:::estimateCdaSpace(tbx3_cnts, tbx3_group)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
color74 = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

sg_pal <- c(C="red", K="gray") # 
palette <- color74 # 

theme_text_rot <- theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
sample.names <- rownames(tbx3_cnts)


gg_boxes <- list(freqs=tbx3_freqs, counts=tbx3_cnts) %>% lapply(function(mat) {
  df <- melt(t(mat)) %>% set_colnames(c('type', 'sample.name', 'value'))
  df[['group']] <- str_remove(df$sample.name,"\\d+")

  ggplot(df, aes(x=type, y=value, fill=group)) +
    geom_boxplot(outlier.shape = NA) +
    # geom_jitter(position=position_jitterdodge(jitter.width=0.2), alpha=1, size = 0.1) +
    scale_fill_manual(values=sg_pal) +
    stat_compare_means(aes(label=cacoa:::pvalueToCode(..p.adj.., ns.symbol="")),
                       label.x=1.5, label.y=max(df$value), size=2.5) +
    theme(legend.title=element_blank(), axis.title.x=element_blank()) +
    theme_text_rot +
    cacoa:::theme_legend_position(c(1, 1))
})

gg_boxes$counts %<>% {. + ylab("Counts")}
gg_boxes$freqs %<>% {. + ylab("Proportions")}

gg_surface <- ggplot(dfs$red, aes(x=S1, y=S2)) +
  geom_abline(slope=-5, intercept=0.1) +
  geom_point(aes(colour=tbx3_sample_group)) +
  geom_hline(yintercept=0, linetype="dashed", size=0.2) +
  geom_vline(xintercept=0, linetype="dashed", size=0.2) +
  labs(colour="Condition", x="CDA-1", y="CDA-2") +
  scale_color_manual(values=sg_pal) +
  theme(panel.grid=element_blank(), legend.position="none")
gg_surface


gg_coda <- res %$% cacoa:::plotCellLoadings(
  loadings, padj, palette=palette, jitter.alpha=0.0,
  ref.level=groups.name[2], target.level=groups.name[1],
  ref.load.level=res$ref.load.level, annotation.x=1.0
) + theme(
  panel.border=element_rect(size=0.1, fill="transparent"),
  panel.grid.minor.x=element_blank()
)
gg_coda # 
gg_coda + theme_bw()
plot_grid(
  gg_boxes$counts + ylab("Counts"),
  gg_boxes$freqs + ylab("Proportions"),
  gg_surface + theme(plot.margin=margin(t=10)),
  gg_coda,
  nrow=1
)
gg_boxes$counts + theme_bw()
gg_boxes$freqs + theme_bw()
tree_theme <- theme(
  legend.key.height=unit(10, "pt"), legend.key.width=unit(14, "pt"),
  legend.position="bottom", plot.margin=margin(),
  axis.text.y=element_text(hjust=1, vjust=0.5, margin=margin()), axis.text.x=element_blank(),
  axis.ticks=element_blank()
)
# gg_toy_tree <- cacoa:::plotContrastTree(
#   tbx3_cnts, tbx3_group, ref.level=groups.name[2], target.level=groups.name[1], plot.theme=NULL,
#   adjust.pvalues=TRUE, loadings.mean=rowMeans(res$loadings), palette=sg_pal
# ) + coord_flip() + tree_theme + theme(legend.margin=margin(l=10, t=-30)) +
#   guides(color=guide_legend(direction="vertical", title="Condition"))
#
# gg_toy_tree
dev.off()
p1 <- gg_coda + theme_bw()+theme(legend.position = "none")+scale_x_continuous(limits = c(-0.5,0.5))
p1




library(tidyverse)
library(ggpubr)
library(ggsci)
library(introdataviz)

df <- melt(t(tbx3_cnts)) %>% set_colnames(c('type', 'sample.name', 'value'))
df[['group']] <- str_remove(df$sample.name,"\\d+")

df <- df[df$type %in% c("N9","N14","N17","N41","N42","N46","N53","N58","N60","N65","N66","N67","N18","N36","N61","N64"),]

ggplot(df,aes(x = type,y = value,fill = group)) +
  # split violin
  geom_split_violin(alpha = .5, trim = F,color = NA,width = 1) +
  # mean point
  stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
  # errorbar
  stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
               size = 0.3,
               position = position_dodge(0.2)) +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90,color = 'black',hjust = 1),
        legend.position = 'top') +
  # scale_fill_brewer(palette = 'Set1') +
  scale_fill_jco(name = '') +
  ylim(0,1000) +
  stat_compare_means(aes(group=group),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "NS")),label = "p.signif",
                     label.y = 7,size = 5)


# devtools/remotes won't install Suggested packages from Bioconductor
BiocManager::install(c("CellBench", "BiocStyle"))
remotes::install_github("Oshlack/speckle",
                        dependencies = "Suggest")
library(speckle)
library(limma)
library(ggplot2)


head(metadata)

tbx3_met <- metadata[,c("nchip","group","RNA_snn_res.4","cell")]
tbx3_met$cell <- factor(tbx3_met$cell,levels = paste0("N",seq(1,67)))
tbx3_met <- tbx3_met[order(tbx3_met$cell),]
propeller(clusters =tbx3_met$nchip, sample =  tbx3_met$RNA_snn_res.4,
          group = tbx3_met$nchip)



p_result <- propeller(clusters =tbx3_met$cell, sample = tbx3_met$nchip,
          group = tbx3_met$group)

# Plot cell type proportions
tbx3_met <- tbx3_met[order(tbx3_met$cell),]
plotCellTypeProps(clusters=tbx3_met$RNA_snn_res.4, sample=tbx3_met$nchip)+theme_bw()
plotCellTypeProps(clusters=tbx3_met$nchip, sample= tbx3_met$RNA_snn_res.4)+theme_bw()
p <- plotCellTypeProps(clusters=tbx3_met$group, sample= tbx3_met$cell)+theme_bw()
p



library(ggpubr)
p_result$group <- "T"
p_result[p_result$FDR > 0.05,]$group <- "F"
p_result$BaselineProp.clusters <- factor(p_result$BaselineProp.clusters,levels = paste0("N",seq(1:67)))
p_result <- p_result[order(p_result$BaselineProp.clusters),]
p_result_plot <- ggdotchart(p_result, x = "BaselineProp.clusters", y = "FDR",
                dot.size = 5,
                color = "group",     
                palette = c("#00AFBB", "#E7B800"), 
                sorting = "none",
                add = "segments",                             
                ggtheme = theme_pubr(),                        
                xlab=""
)
p_result_plot


############# scConda ---------------------------------------------------


library(ggpubr)


scConda_result <- read.csv("E:/wu_lab/For XS/20220709/text_result_fdr0.001_effect_df.csv")
scConda_result$Cell.Type <- factor(scConda_result$Cell.Type,levels = paste0("N",seq(1:67)))
p1_sort <- ggbarplot(scConda_result, x = "Cell.Type", y = "log2.fold.change",
          fill = "Final.Parameter.1",           
          color = "white",            
          palette = "jco",            
          sort.val = "desc",          
          sort.by.groups = FALSE,     
          x.text.angle = 90,          
          ylab = "expect-log2.fold.change",
          legend.title = "Model result",
          rotate = TRUE,
          ggtheme = theme_minimal()
)
p1_sort # 8 4

p2_unsort <- ggbarplot(scConda_result, x = "Cell.Type", y = "log2.fold.change",
                     fill = "Final.Parameter.1",           
                     color = "white",            
                     palette = "jco",            
                     sort.by.groups = FALSE,     
                     x.text.angle = 90,          
                     ylab = "expect-log2.fold.change",
                     legend.title = "Model result",
                     rotate = TRUE,
                     ggtheme = theme_minimal()
)
p2_unsort # 8 6

p_result$BaselineProp.clusters <- factor(p_result$BaselineProp.clusters,levels = paste0("N",seq(1:67)))
p_result <- p_result[order(p_result$BaselineProp.clusters),]
p_result$group <- "TRUE"
p_result[p_result$FDR > 0.05,"group"] <- "FALSE"

p_result$FDR_log <- -log2(p_result$FDR+1)

p1_sort_method2 <- ggbarplot(p_result, x = "BaselineProp.clusters", y = "FDR_log",
                     fill = "group",           
                     color = "white",            
                     palette = "jco",           
                     sort.val = "desc",          
                     sort.by.groups = FALSE,     
                     x.text.angle = 90,          
                     ylab = "expect-log2.fold.change",
                     legend.title = "Model result",
                     rotate = TRUE,
                     ggtheme = theme_minimal()
)
p1_sort_method2 # 10 3
ggdotchart(p_result, x = "BaselineProp.clusters", y = "FDR",
           color = "group",                                
           palette = c("#00AFBB", "#E7B800"), # Custom color palette
           add = "segments",                             
           rotate = TRUE,                                
           group = "group",                                
           dot.size = 3,                                        
           font.label = list(color = "white", size = 9,
                             vjust = 0.5),               
           ggtheme = theme_pubr()                       
)

########  combine 3 sofeware
p_result_combine <- data.frame(cell = p_result$BaselineProp.clusters,
                               scConda_result = scConda_result$Final.Parameter.1,
                               p_result = p_result$group,
                               gg_coda = c(FALSE,rep(TRUE,3),FALSE,rep(TRUE,4),FALSE,FALSE,FALSE,rep(TRUE,9),FALSE,FALSE,TRUE,FALSE,rep(TRUE,12),FALSE,rep(TRUE,4),FALSE,TRUE,FALSE,TRUE,FALSE,
                                           rep(TRUE,3),FALSE,rep(TRUE,5),FALSE,TRUE,FALSE,rep(TRUE,8))


)

p_result_combine[p_result_combine==TRUE] <- as.numeric(1)
p_result_combine[p_result_combine==FALSE] <- as.numeric(0)
p_result_combine$cell <- paste0("d",seq(1,67))
rownames(p_result_combine) <- p_result_combine$cell
p_result_combine <- p_result_combine[,-1]
p_result_combine$p_result <- as.numeric(p_result_combine$p_result)
pheatmap::pheatmap(p_result_combine,cluster_rows=F,
                   cluster_cols=F,border_color = "white",
                   color = c("#424242","#FFEB3B")) # colorRampPalette(c("#E0E0E0","#FFEB3B"))(256)
pheatmap::pheatmap(p_result_combine,cluster_rows=F,
                   cluster_cols=F,border_color = "white",
                   color = c("#000000","#DE2D26")) # colorRampPalette(c("#E0E0E0","#FFEB3B"))(256)


## 22 3
pheatmap::pheatmap(p_result_combine,cluster_rows=F,
                   cluster_cols=F,border_color = "white",
                   color = c("#64B5F6","#E57373")) # colorRampPalette(c("#E0E0E0","#FFEB3B"))(256)










