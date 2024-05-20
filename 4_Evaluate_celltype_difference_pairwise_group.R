#
##
###

library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(distdimscr)
library(ggpubr)

load("/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/3_integrateData_with_cellType_annotation.RData")

# sample info
# MB.combined$sample = mapvalues(MB.combined$sample,"1","MB757")
# MB.combined$sample = mapvalues(MB.combined$sample,"2","MB758")
# MB.combined$sample = mapvalues(MB.combined$sample,"3","MB759")
# MB.combined$sample = mapvalues(MB.combined$sample,"4","MB760")
MB.combined$sample = mapvalues(MB.combined$sample,"5","MB761")
MB.combined$sample = mapvalues(MB.combined$sample,"6","MB762")
MB.combined$sample = mapvalues(MB.combined$sample,"7","MB763")
MB.combined$sample = mapvalues(MB.combined$sample,"8","MB764")
MB.combined$sample = mapvalues(MB.combined$sample,"9","MB765")
MB.combined$sample = mapvalues(MB.combined$sample,"10","MB766")
MB.combined$sample = mapvalues(MB.combined$sample,"11","MB767")
MB.combined$sample = mapvalues(MB.combined$sample,"12","MB768")
MB.combined$sample = mapvalues(MB.combined$sample,"13","MB785")
MB.combined$sample = mapvalues(MB.combined$sample,"14","MB786")
MB.combined$sample = mapvalues(MB.combined$sample,"15","MB787")
MB.combined$sample = mapvalues(MB.combined$sample,"16","MB788")

# group info
# MB.combined$group = mapvalues(MB.combined$group,"group1","injury_acetate")
MB.combined$group = mapvalues(MB.combined$group,"group2","injury_SCFA")
MB.combined$group = mapvalues(MB.combined$group,"group3","injury_control")
MB.combined$group = mapvalues(MB.combined$group,"group4","sham_control")

# count each type cells
a = table(MB.combined$sample,MB.combined$celltype)
bar_input = as.data.frame(a)
colnames(bar_input) = c("sample","celltype","Freq")

col = c(RColorBrewer::brewer.pal(9,'Set1'),
        RColorBrewer::brewer.pal(6,'Pastel1'),
        RColorBrewer::brewer.pal(10,'Set3'))
bar_input %>%
  ggplot(aes(x = sample, y = Freq, fill = celltype))+
  geom_bar(stat="identity", position = 'fill')+
  scale_fill_manual(values= col)
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/4_celltype_frequency_barplot.pdf",width = 10,height = 8)

#
b = as.data.frame(matrix(a[1:length(a)],nrow = length(row.names(a)),ncol = length(colnames(a))))
row.names(b) = row.names(a)
colnames(b) = colnames(a)
b = prop.table(a,1)
b = as.data.frame(b)
colnames(b) = c("sample","celltype","Freq")
b$group = if_else(b$sample %in% c("MB761","MB762","MB763","MB764"),"3_injury_SCFA",
                  if_else(b$sample %in% c("MB765","MB766","MB767","MB768"),"2_injury_control","1_sham_control"))

dodge <- position_dodge(width = 0.8)

b%>%ggplot(aes(celltype,Freq))+geom_boxplot(aes(fill=group),scale = "width",position = dodge,width = 0.7)+
  # geom_boxplot(aes(fill=group),position = dodge,outlier.size = 0,notch = F,width = 0.3)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "right"
  )+
  ggtitle("celltype frequency")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/4_celltype_frequency_group_comparison.pdf",height = 6,width = 8)

## chi-squre plot for Freq of each group
a = table(MB.combined$group,MB.combined$celltype)
id = colnames(a)
a = as.vector(a)
dim(a) = c(3,length(unique(MB.combined$celltype)))
a = data.frame(a)
rownames(a) = c("injury_control","injury_SCFA","sham_control")
colnames(a) = id

chi_p = data.frame()
for(i in 1:ncol(a)){
  ic_is = matrix(c(a[1,i],a[2,i],rowSums(a)[1]-a[1,i],rowSums(a)[2]-a[2,i]),ncol = 2,nrow = 2)
  chi_R = chisq.test(ic_is)
  chi_p[3,i] = chi_R$p.value
  
  ic_sc = matrix(c(a[1,i],a[3,i],rowSums(a)[1]-a[1,i],rowSums(a)[3]-a[3,i]),ncol = 2,nrow = 2)
  chi_R = chisq.test(ic_sc)
  chi_p[1,i] = chi_R$p.value
  
  is_sc = matrix(c(a[2,i],a[3,i],rowSums(a)[2]-a[2,i],rowSums(a)[3]-a[3,i]),ncol = 2,nrow = 2)
  chi_R = chisq.test(is_sc)
  chi_p[2,i] = chi_R$p.value
  
}
rownames(chi_p) = c("injury_control_vs_sham_control",
                    "injury_SCFA_vs_sham_control",
                    "injury_SCFA_vs_injury_control")
colnames(chi_p) = colnames(a)
write.csv(chi_p,file = "/ix/gkohanbash/SimonsLab/R_Share/table/No_acetate_group/4_combine_celltype_chi_square_p_value_matrix.csv")

chi_p[chi_p < 1e-4] = 4
chi_p[chi_p > 0.05 & chi_p <= 1 | chi_p == 'NaN'] = 0
chi_p[chi_p > 0.01 & chi_p < 0.05] = 1
chi_p[chi_p > 0.001 & chi_p < 0.01] = 2
chi_p[chi_p > 0.0001 & chi_p < 0.001] = 3

color = RColorBrewer::brewer.pal(9,"Reds")[c(6,7,9)]
pheatmap::pheatmap(chi_p,cellwidth = 10,cellheight = 10,
                   cluster_rows = F,cluster_cols = F,
                   color = color)
# save /ix/gkohanbash/SimonsLab/R_Share/table/No_acetate_group/4_combined_celltype_chi_square_p_value_heatmap.pdf 4*6 

# Bhattacharyya distance injury vs sham
MB.combined@active.ident = as.factor(MB.combined$group)
MB.Bhattacharyya = subset(MB.combined,idents = c("injury_control","sham_control"))
comparison = "group"

## Metadata
overall.metadata <- MB.Bhattacharyya@meta.data
format(object.size(overall.metadata),unit="MB")

## UMAP
overall.umap <- Embeddings(MB.Bhattacharyya,reduction="umap")
format(object.size(overall.umap),unit="MB")

## PCA embeddings
overall.pca <- Embeddings(MB.Bhattacharyya,reduction="pca")[,1:30]
format(object.size(overall.pca),unit="MB")

# check out the same cell type between conditions
overall.data <- cbind(overall.umap,overall.metadata)

# Check out cell numbers in each sample
knitr::kable(table(overall.data[,comparison],overall.data$celltype))
type_number = as.data.frame(table(overall.data[,comparison],overall.data$celltype))
var = as.character(unique(type_number$Var1))
type_number = type_number[type_number$Freq >= 500,]
type_number = type_number[duplicated(type_number$Var2),]
condition = as.character(unique(type_number$Var2))

bhatt.total = data.frame()
for (j in 1:length(condition)) {
  
  condition1 <- rownames(overall.data)[overall.data$celltype==condition[j] & overall.data[,comparison]==var[1]]
  condition2 <- rownames(overall.data)[overall.data$celltype==condition[j] & overall.data[,comparison]==var[2]]
  
  condition1.pca <- overall.pca[condition1,]
  condition2.pca <- overall.pca[condition2,]
  
  bhatt.dist <- bhatt.dist.rand <- vector("logical",length=100)
  set.seed("0222")
  
  for (i in 1:100) {
    bhatt.dist[[i]] <- dim_dist(embed_mat_x=condition1.pca,embed_mat_y=condition2.pca,dims_use=1:30,num_cells_sample=500,distance_metric="bhatt_dist",random_sample=FALSE)
    bhatt.dist.rand[[i]] <- dim_dist(embed_mat_x=condition1.pca,embed_mat_y=condition2.pca,dims_use=1:30,num_cells_sample=500,distance_metric="bhatt_dist",random_sample=TRUE)
  }
  
    # Combine the results and plot
  bhatt.dist <- data.frame(cells.distance=bhatt.dist,comparison="real")
  bhatt.dist.rand <- data.frame(cells.distance=bhatt.dist.rand,comparison="random")
  
  bhatt.res <- rbind(bhatt.dist,bhatt.dist.rand)
  bhatt.res$celltype = condition[j]
  bhatt.total = rbind(bhatt.total,bhatt.res)
}
fill = numeric()
logfc = numeric()
for (i in 1:length(condition)) {
  logfc[i] = round(mean(bhatt.total$cells.distance[which(bhatt.total$comparison == "real" & bhatt.total$celltype == condition[i])]) / mean(bhatt.total$cells.distance[which(bhatt.total$comparison == "random" & bhatt.total$celltype == condition[i])]),digits = 2)
  fill[2*i-1] = mean(bhatt.total$cells.distance[which(bhatt.total$comparison == "random" & bhatt.total$celltype == condition[i])])
  fill[2*i] = mean(bhatt.total$cells.distance[which(bhatt.total$comparison == "real" & bhatt.total$celltype == condition[i])])
}

label = data.frame(celltype = condition,cells.distance = rep(max(bhatt.total$cells.distance)-0.2,length(condition)), fc = logfc)
label$fc = paste0(label$fc,"-fold")

bhatt.total$comparison = mapvalues(bhatt.total$comparison,"real","injury control vs sham control")
bhatt.is = bhatt.total
ggbarplot(bhatt.total,x="celltype",y="cells.distance",color = "comparison",position = position_dodge(0.8),add = "mean_se",fill = "comparison") +
  theme_bw() +
  xlab("Comparison type") +
  ylab("Bhattacharrya distance")+
  stat_compare_means(aes(group = bhatt.total$comparison),label.y = c(max(bhatt.total$cells.distance)-0.1),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),label = "p.signif",method = "wilcox")+
  geom_text(data = label, aes(label = fc))
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/4_Bhattacharrya_distance_injury_control_vs_sham.pdf",width = 12,height = 5)

# Bhattacharyya distance injury_acetate vs injury
library(distdimscr)
library(ggplot2)
library(Seurat)
library(ggpubr)

MB.Bhattacharyya = subset(MB.combined,idents = c("injury_control","injury_SCFA"))
comparison = "group"

## Metadata
overall.metadata <- MB.Bhattacharyya@meta.data
format(object.size(overall.metadata),unit="MB")

## UMAP
overall.umap <- Embeddings(MB.Bhattacharyya,reduction="umap")
format(object.size(overall.umap),unit="MB")

## PCA embeddings
overall.pca <- Embeddings(MB.Bhattacharyya,reduction="pca")[,1:30]
format(object.size(overall.pca),unit="MB")

# check out the same cell type between conditions
overall.data <- cbind(overall.umap,overall.metadata)

# Check out cell numbers in each sample
knitr::kable(table(overall.data[,comparison],overall.data$celltype))
type_number = as.data.frame(table(overall.data[,comparison],overall.data$celltype))
var = as.character(unique(type_number$Var1))
type_number = type_number[type_number$Freq >= 500,]
type_number = type_number[duplicated(type_number$Var2),]
condition = as.character(unique(type_number$Var2))

bhatt.total = data.frame()
for (j in 1:length(condition)) {
  
  condition1 <- rownames(overall.data)[overall.data$celltype==condition[j] & overall.data[,comparison]==var[1]]
  condition2 <- rownames(overall.data)[overall.data$celltype==condition[j] & overall.data[,comparison]==var[2]]
  
  condition1.pca <- overall.pca[condition1,]
  condition2.pca <- overall.pca[condition2,]
  
  bhatt.dist <- bhatt.dist.rand <- vector("logical",length=100)
  set.seed("0222")
  
  for (i in 1:100) {
    bhatt.dist[[i]] <- dim_dist(embed_mat_x=condition1.pca,embed_mat_y=condition2.pca,dims_use=1:30,num_cells_sample=500,distance_metric="bhatt_dist",random_sample=FALSE)
    bhatt.dist.rand[[i]] <- dim_dist(embed_mat_x=condition1.pca,embed_mat_y=condition2.pca,dims_use=1:30,num_cells_sample=500,distance_metric="bhatt_dist",random_sample=TRUE)
  }
  
  # Combine the results and plot
  bhatt.dist <- data.frame(cells.distance=bhatt.dist,comparison="real")
  bhatt.dist.rand <- data.frame(cells.distance=bhatt.dist.rand,comparison="random")
  
  bhatt.res <- rbind(bhatt.dist,bhatt.dist.rand)
  bhatt.res$celltype = condition[j]
  bhatt.total = rbind(bhatt.total,bhatt.res)
}
fill = numeric()
logfc = numeric()
for (i in 1:length(condition)) {
  logfc[i] = round(mean(bhatt.total$cells.distance[which(bhatt.total$comparison == "real" & bhatt.total$celltype == condition[i])]) / mean(bhatt.total$cells.distance[which(bhatt.total$comparison == "random" & bhatt.total$celltype == condition[i])]),digits = 2)
  fill[2*i-1] = mean(bhatt.total$cells.distance[which(bhatt.total$comparison == "random" & bhatt.total$celltype == condition[i])])
  fill[2*i] = mean(bhatt.total$cells.distance[which(bhatt.total$comparison == "real" & bhatt.total$celltype == condition[i])])
}

label = data.frame(celltype = condition,cells.distance = rep(max(bhatt.total$cells.distance)-0.2,length(condition)), fc = logfc)
label$fc = paste0(label$fc,"-fold")

bhatt.total$comparison = mapvalues(bhatt.total$comparison,"real","injury SCFA vs injury control")
bhatt.SCFAi = bhatt.total
ggbarplot(bhatt.total,x="celltype",y="cells.distance",color = "comparison",position = position_dodge(0.8),add = "mean_se",fill = "comparison") +
  theme_bw() +
  xlab("Comparison type") +
  ylab("Bhattacharrya distance")+
  stat_compare_means(aes(group = bhatt.total$comparison),label.y = c(max(bhatt.total$cells.distance)-0.1),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),label = "p.signif",method = "wilcox")+
  geom_text(data = label, aes(label = fc))
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/4_Bhattacharrya_distance_injury_SCFA_vs_injury_control.pdf",width = 12,height = 5)

# Bhattacharyya distance injury butyrate vs acetate
MB.Bhattacharyya = subset(MB.combined,idents = c("injury_SCFA","sham_control"))
comparison = "group"

## Metadata
overall.metadata <- MB.Bhattacharyya@meta.data
format(object.size(overall.metadata),unit="MB")

## UMAP
overall.umap <- Embeddings(MB.Bhattacharyya,reduction="umap")
format(object.size(overall.umap),unit="MB")

## PCA embeddings
overall.pca <- Embeddings(MB.Bhattacharyya,reduction="pca")[,1:30]
format(object.size(overall.pca),unit="MB")

# check out the same cell type between conditions
overall.data <- cbind(overall.umap,overall.metadata)

# Check out cell numbers in each sample
knitr::kable(table(overall.data[,comparison],overall.data$celltype))
type_number = as.data.frame(table(overall.data[,comparison],overall.data$celltype))
var = as.character(unique(type_number$Var1))
type_number = type_number[type_number$Freq >= 500,]
type_number = type_number[duplicated(type_number$Var2),]
condition = as.character(unique(type_number$Var2))

bhatt.total = data.frame()
for (j in 1:length(condition)) {
  
  condition1 <- rownames(overall.data)[overall.data$celltype==condition[j] & overall.data[,comparison]==var[1]]
  condition2 <- rownames(overall.data)[overall.data$celltype==condition[j] & overall.data[,comparison]==var[2]]
  
  condition1.pca <- overall.pca[condition1,]
  condition2.pca <- overall.pca[condition2,]
  
  bhatt.dist <- bhatt.dist.rand <- vector("logical",length=100)
  set.seed("0222")
  
  for (i in 1:100) {
    bhatt.dist[[i]] <- dim_dist(embed_mat_x=condition1.pca,embed_mat_y=condition2.pca,dims_use=1:30,num_cells_sample=500,distance_metric="bhatt_dist",random_sample=FALSE)
    bhatt.dist.rand[[i]] <- dim_dist(embed_mat_x=condition1.pca,embed_mat_y=condition2.pca,dims_use=1:30,num_cells_sample=500,distance_metric="bhatt_dist",random_sample=TRUE)
  }
  
  # Combine the results and plot
  bhatt.dist <- data.frame(cells.distance=bhatt.dist,comparison="real")
  bhatt.dist.rand <- data.frame(cells.distance=bhatt.dist.rand,comparison="random")
  
  bhatt.res <- rbind(bhatt.dist,bhatt.dist.rand)
  bhatt.res$celltype = condition[j]
  bhatt.total = rbind(bhatt.total,bhatt.res)
}
fill = numeric()
logfc = numeric()
for (i in 1:length(condition)) {
  logfc[i] = round(mean(bhatt.total$cells.distance[which(bhatt.total$comparison == "real" & bhatt.total$celltype == condition[i])]) / mean(bhatt.total$cells.distance[which(bhatt.total$comparison == "random" & bhatt.total$celltype == condition[i])]),digits = 2)
  fill[2*i-1] = mean(bhatt.total$cells.distance[which(bhatt.total$comparison == "random" & bhatt.total$celltype == condition[i])])
  fill[2*i] = mean(bhatt.total$cells.distance[which(bhatt.total$comparison == "real" & bhatt.total$celltype == condition[i])])
}

label = data.frame(celltype = condition,cells.distance = rep(max(bhatt.total$cells.distance)-0.2,length(condition)), fc = logfc)
label$fc = paste0(label$fc,"-fold")

bhatt.total$comparison = mapvalues(bhatt.total$comparison,"real","injury SCFA vs sham control")
ggbarplot(bhatt.total,x="celltype",y="cells.distance",color = "comparison",position = position_dodge(0.8),add = "mean_se",fill = "comparison") +
  theme_bw() +
  xlab("Comparison type") +
  ylab("Bhattacharrya distance")+
  stat_compare_means(aes(group = bhatt.total$comparison),label.y = c(max(bhatt.total$cells.distance)-0.1),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),label = "p.signif",method = "wilcox")+
  geom_text(data = label, aes(label = fc))
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/4_Bhattacharrya_distance_injury_SCFA_vs_sham_control.pdf",width = 10,height = 5)

# integrate 3 comparison
bhatt.total = rbind(bhatt.SCFAi,bhatt.total,bhatt.is)
ggbarplot(bhatt.total,x="celltype",y="cells.distance",color = "comparison",position = position_dodge(0.8),add = "mean_se",fill = "comparison") +
  theme_bw() +
  xlab("Comparison type") +
  ylab("Bhattacharrya distance")+
  stat_compare_means(aes(group = bhatt.total$comparison),label.y = c(max(bhatt.total$cells.distance)-0.1),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),label = "p.signif",method = "wilcox")+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/4_Bhattacharrya_distance_all_group_pairwise_comparison.pdf",width = 10,height = 5)

# overview 
DefaultAssay(MB.combined) = "SCT"
save(MB.combined,file = "/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/4_integrateData_with_group_info.RData")

MB.combined$group = mapvalues(MB.combined$group,"sham_control","1_sham_control")
MB.combined$group = mapvalues(MB.combined$group,"injury_control","2_injury_control")
MB.combined$group = mapvalues(MB.combined$group,"injury_SCFA","3_injury_SCFA")

# boxplot of nGene & nUMI in different celltype and groups
vln.df=data.frame(nUMI = MB.combined$nCount_SCT,
                  nGene = MB.combined$nFeature_SCT)
vln.df$celltype = MB.combined$celltype
vln.df$group = MB.combined$group

dodge <- position_dodge(width = 1)
vln.df%>%ggplot(aes(celltype,nUMI))+geom_violin(aes(fill=group),scale = "width",position = dodge)+
  geom_boxplot(aes(fill=group),position = dodge,outlier.size = 0,notch = T,width = 0.3)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "right"
  )+
  ggtitle("nUMI")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/4_group_celltype_nUMI_distribution.pdf",width = 10,height = 5)

vln.df%>%ggplot(aes(celltype,nGene))+geom_violin(aes(fill=group),scale = "width",position = dodge)+
  geom_boxplot(aes(fill=group),position = dodge,outlier.size = 0,notch = T,width = 0.3)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "right"
  )+
  ggtitle("nGene")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/4_group_celltype_nGene_distribution.pdf",width = 10,height = 5)
