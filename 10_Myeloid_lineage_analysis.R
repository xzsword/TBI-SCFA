#
##
###


# Myeloid_lineage
library(Seurat)
suppressMessages(library(clusterProfiler))
library(enrichplot)
library(ggplot2)
library(AUCell)
library(msigdbr) 
library(dplyr)
library(ballgown)
library(plyr)
library(scWGCNA)
library(Matrix)
library(tidyverse)
library(WGCNA)
library(SeuratObject)
library(corrplot)
library(condiments)
library(slingshot)
library(tidyr)
library(RColorBrewer)
library(viridis)
library(tradeSeq)
library(CytoTRACE)
library(scales)
library(ggpubr)
nCores = 16

load("/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/4_integrateData_with_group_info.RData")

# Myeloid lineage
MB.combined@active.ident = as.factor(MB.combined$celltype)

# get Myeloid cells
MB.Myeloid = MB.combined[,MB.combined$celltype %in% c("Myeloid_cell")]
# MB.Myeloid$operation = if_else(MB.Myeloid$batch =="batch1","op1","op2")

MB.Myeloid$group = mapvalues(MB.Myeloid$group,"sham_control","1_sham_control")
MB.Myeloid$group = mapvalues(MB.Myeloid$group,"injury_control","2_injury_control")
MB.Myeloid$group = mapvalues(MB.Myeloid$group,"injury_SCFA","3_injury_SCFA")

MB.Myeloid = SCTransform(MB.Myeloid,
                           variable.features.n = 3000,
                           ncells = 5000,
                           vars.to.regress = c("percent.mt"))

vln.df=data.frame(nUMI = MB.Myeloid$nCount_SCT,
                  nGene = MB.Myeloid$nFeature_SCT)
vln.df$celltype = MB.Myeloid$celltype
vln.df$group = MB.Myeloid$group

dodge <- position_dodge(width = 1)
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

# recluster Myeloid
MB.Myeloid = RunPCA(MB.Myeloid,npcs = 50, assay = "SCT")
ElbowPlot(MB.Myeloid, ndims=50)
pc.num=1:12

MB.Myeloid = RunUMAP(MB.Myeloid , reduction = "pca", dims = pc.num)
MB.Myeloid = FindNeighbors(MB.Myeloid, verbose = FALSE, reduction = "pca",dims = pc.num)
MB.Myeloid = FindClusters(MB.Myeloid,resolution = 0.5)

p1 = UMAPPlot(MB.Myeloid,group.by = "group")
p2 = UMAPPlot(MB.Myeloid,group.by = "seurat_clusters",label = T) + NoLegend()
cowplot::plot_grid(p1,p2)
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_recluster_UMAP.pdf",width = 12,height = 5)

subtype = MB.Myeloid$seurat_clusters

subtype = mapvalues(subtype,"0","MG_1")
subtype = mapvalues(subtype,"1","MG_2")
subtype = mapvalues(subtype,"2","MG_3")
subtype = mapvalues(subtype,"3","MG_4")
subtype = mapvalues(subtype,"4","MG_5")
subtype = mapvalues(subtype,"5","MG_6")
subtype = mapvalues(subtype,"6","MAC")
subtype = mapvalues(subtype,"7","MG_7")
subtype = mapvalues(subtype,"8","MG_8")
subtype = mapvalues(subtype,"9","MG_9")
subtype = mapvalues(subtype,"10","MG_10")

MB.Myeloid$subtype = subtype

# subtype identification
# seurat object Marker gene shown in violin plot
gene = c("Gapdh","Cd14","Cx3cr1","Fcrls", # Macrophage/microglia
         "Tmem119", # MG
         "Tnf", # immediate-early genes
         "Top2a", # proliferative
         "Plac8", # MNC
         "Pf4","Mrc1","Cd68","Cd163" # MAC
)

# subtype is the feature for x axis

library(reshape2)
vln.df=as.data.frame(MB.Myeloid[["SCT"]]@data[gene,])
vln.df$gene=rownames(vln.df)
vln.df=melt(vln.df,id="gene")
colnames(vln.df)[c(2,3)]=c("CB","exp")

MB.Myeloid$CB = colnames(MB.Myeloid)
anno=MB.Myeloid@meta.data[,c("CB","subtype")]
vln.df=inner_join(vln.df,anno,by="CB")
vln.df$gene=factor(vln.df$gene,levels = gene) #为了控制画图的基因顺序

vln.df%>%ggplot(aes(subtype,exp))+geom_violin(aes(fill=gene),scale = "width")+
  facet_grid(vln.df$gene~.,scales = "free_y")+
  scale_fill_brewer(palette = "Set3",direction = 1)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "none"
  )
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_subtype_marker_expression_violin.pdf",height = 8,width = 6)

UMAPPlot(MB.Myeloid,group.by = "subtype",label = T) + NoLegend()
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_subtype_UMAP.pdf",width = 5,height = 5)

## identify M1/M2 macrophage
# get M1 and M2 marker from msigdb C2
load("/ix/gkohanbash/SimonsLab/R_Share/processing/5_KEGG_AUC_genelist.RData")
marker = h.sets[grep("MACROPHAGE_M1_VS_M2",names(h.sets))]
names(marker) = c("M2","M1")

# use AUCell to score the cells
# AUC
# input -- genes/features as rows and cells as columns.
exprMatrix <- MB.Myeloid@assays$SCT@counts # counts
exprMatrix <- as.matrix(exprMatrix)
cellsUMAP = MB.Myeloid@reductions$umap@cell.embeddings

# ranking cells by gene expression high -> low
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=nCores, plotStats=F)

#Calculate enrichment for the gene signatures (AUC)
cells_AUC <- AUCell_calcAUC(marker, cells_rankings,
                            aucMaxRank=nrow(cells_rankings)*0.05)

save(cells_AUC, file = "/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/10_Myeloid_macrophage_M1_M2_AUC.RData")

M.marker.score = as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
M.marker.score$final = apply(M.marker.score,1,max)
M_type = if_else(M.marker.score$M2 < M.marker.score$final, "M1","M2")

MB.Myeloid$M_type = M_type
UMAPPlot(MB.Myeloid,group.by="M_type",split.by="group")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_M_type_UMAP.pdf",width = 12,height = 5)

# compare the M2/M1 Freq in each group
boxinput = data.frame()
for(i in 1:length(unique(MB.Myeloid$sample))){
  for(j in 1:length(unique(MB.Myeloid$subtype))){
a = table(MB.Myeloid$M_type[MB.Myeloid$sample == unique(MB.Myeloid$sample)[i] & MB.Myeloid$subtype == unique(MB.Myeloid$subtype)[j]])
a = as.vector(a)
boxinput[length(unique(MB.Myeloid$subtype))*(i-1)+j,1] = a[1]
boxinput[length(unique(MB.Myeloid$subtype))*(i-1)+j,2] = a[2]
boxinput[length(unique(MB.Myeloid$subtype))*(i-1)+j,3] = unique(MB.Myeloid$sample)[i]
boxinput[length(unique(MB.Myeloid$subtype))*(i-1)+j,4] = unique(MB.Myeloid$subtype)[j]   
  }
}
colnames(boxinput) = c("M1","M2","sample","subtype")
boxinput$group = if_else(boxinput$sample %in% unique(MB.Myeloid$sample)[1:4],"3_injury_SCFA",
                         if_else(boxinput$sample %in% unique(MB.Myeloid$sample)[5:8],"2_injury_control","1_sham_control"))

boxinput = na.omit(boxinput)

# calculate M2/M1 Freq in each sample
boxinput$`M2/M1` = boxinput$M2/boxinput$M1

vln.df=data.frame(nUMI = MB.combined$nCount_SCT,
                  nGene = MB.combined$nFeature_SCT)
vln.df$celltype = MB.combined$celltype
vln.df$group = MB.combined$group

dodge <- position_dodge(width = 0.7)
boxinput%>%ggplot(aes(subtype,`M2/M1`))+
  geom_boxplot(aes(fill=group),position = dodge,outlier.size = 0,notch = F,width = 0.5)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "right"
  )+
  ggtitle("M2/M1")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_M_type_Freq_boxplot.pdf",width = 12,height = 5)

## M1 type to M2 type repolarization in injury SCFA group
M.marker.score$group = MB.Myeloid$group
M.marker.score$type = M_type
M1 = M.marker.score %>% filter(type == "M1")

dataset = data.frame(M2.score = M1$M2,
                     group = M1$group)
dataset = na.omit(dataset)

group=levels(factor(dataset$group))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

dodge <- position_dodge(width = 1)

dataset%>%ggplot(aes(group,M2.score))+geom_violin(aes(fill=group),scale = "width",position = dodge)+
  # facet_grid(dataset$condition~.,scales = "free_y")+
  geom_boxplot(aes(fill=group),position = dodge,outlier.size = 0,notch = T,width = 0.3)+
  # scale_fill_brewer(palette = "Set3",direction = 1)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1,size = 13),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "right",
    axis.text.y.left = element_text(size = 15),
    legend.text = element_text(size = 13)
  )+
  stat_compare_means(comparisons = my_comparisons,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),label = "p.signif")+
  ggtitle("M2.score in M1 Macrophage/Microglia")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_M2_score_in_M1_Cells_group_violin.pdf",width = 5,height = 5)

# MG M1/M2 validation
boxinput = data.frame(Arg1 = as.numeric(GetAssayData(MB.Myeloid["Arg1",])),
                      subtype = MB.Myeloid$subtype,
                      group = factor(MB.Myeloid$group,level = c("sham_control","injury_control","injury_SCFA")))

boxinput = boxinput[boxinput$subtype != "MAC",]

dodge <- position_dodge(width = 0.7)
boxinput%>%ggplot(aes(subtype,Arg1))+
  geom_boxplot(aes(fill=group),position = dodge,outlier.size = 0,notch = F,width = 0.5)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "right"
  )+
  stat_compare_means(aes(group = group),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),label = "p.signif",method = "anova")+
  ggtitle("Arg1")

exprMatrix <- MB.Myeloid@assays$SCT@counts# counts
exprMatrix <- as.matrix(exprMatrix)
cellsUMAP = MB.Myeloid@reductions$umap@cell.embeddings
logMat <- log2(exprMatrix+2)

# GO
## build GO ref geneset -- gene name, not entrenz 
load("/ix/gkohanbash/SimonsLab/R_Share/processing/5_KEGG_AUC_genelist.RData")
d = data.frame(pathway = names(h.sets))

## Neurodegenerative disease
geneSetName = c("WP_TNF_ALPHA_SIGNALING_PATHWAY","WP_IL1_SIGNALING_PATHWAY","WP_IL6_SIGNALING_PATHWAY", # M1 like MG
                "WP_IL4_SIGNALING_PATHWAY","WP_IL10_ANTIINFLAMMATORY_SIGNALING_PATHWAY") # M2 like MG

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)

for (i in 1:length(geneSetName)) {
  boxinput = data.frame(Score = as.numeric(meanByGs[i,]),
                        subtype = MB.Myeloid$subtype,
                        group = factor(MB.Myeloid$group,level = c("sham_control","injury_control","injury_SCFA")))
  
  boxinput = boxinput[boxinput$subtype != "MAC",]
  
  # subtype separate
  dodge <- position_dodge(width = 0.7)
  boxinput%>%ggplot(aes(subtype,Score))+
    geom_boxplot(aes(fill=group),position = dodge,outlier.size = 0,notch = F,width = 0.5)+
    scale_x_discrete("")+scale_y_continuous("")+
    theme_bw()+
    theme(
      axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
      panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
      legend.position = "right"
    )+
    stat_compare_means(aes(group = group),
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("***", "**", "*", "ns")),label = "p.signif",method = "anova")+
    ggtitle(paste0("MG_",geneSetName[i]))
  ggsave(paste0("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_MG_subtype_",geneSetName[i],".pdf"))
  
  # celltype integrate - violin
  group=levels(factor(boxinput$group))
  comp=combn(group, 2)
  my_comparisons=list()
  for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}
  
  dodge <- position_dodge(width = 1)
  
  boxinput%>%ggplot(aes(group,Score))+geom_violin(aes(fill=group),scale = "width",position = dodge)+
    # facet_grid(dataset$condition~.,scales = "free_y")+
    geom_boxplot(aes(fill=group),position = dodge,outlier.size = 0,notch = T,width = 0.3)+
    # scale_fill_brewer(palette = "Set3",direction = 1)+
    scale_x_discrete("")+scale_y_continuous("")+
    theme_bw()+
    theme(
      axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1,size = 13),
      panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
      legend.position = "right",
      axis.text.y.left = element_text(size = 15),
      legend.text = element_text(size = 13)
    )+
    stat_compare_means(comparisons = my_comparisons,
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("***", "**", "*", "ns")),label = "p.signif")+
    ggtitle(paste0("MG_",geneSetName[i]))
  ggsave(paste0("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_MG_",geneSetName[i],".pdf"),height = 5,width = 5)
}


# compare the DEGs number among groups
MB.Myeloid@active.ident = as.factor(MB.Myeloid$group)
markers =  FindAllMarkers(MB.Myeloid,
                          only.pos = T,
                          test.use = "MAST",
                          logfc.threshold = 0.25)

library(UpSetR)
UpSetR::upset(fromList(list(sham_control = markers$gene[markers$cluster == "sham_control"], 
                            injury_control = markers$gene[markers$cluster == "injury_control"],
                            injury_SCFA = markers$gene[markers$cluster == "injury_SCFA"])))
# save 	10_Myeloid_group_DEGs_counts.pdf 6*8

# get SCFA diet group specific DEGs
dup = markers[duplicated(markers$gene),]
SCFA_DEG = markers[markers$cluster %in% c("injury_SCFA"),]
SCFA_DEG = setdiff(SCFA_DEG$gene,dup$gene)

go <- enrichGO(gene = SCFA_DEG, OrgDb = "org.Mm.eg.db", ont="BP",qvalueCutoff = 0.05,keyType = "SYMBOL",pvalueCutoff = 0.05)
go@result = go@result[-grep("viral",go@result$Description),]
go@result = go@result[-grep("virus",go@result$Description),]

edox2 <- pairwise_termsim(go)
treeplot(edox2, hclust_method = "average")
ggsave(paste0("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_SCFA_group_DEG_GO_treeplot.pdf"),width = 20,height = 8)

# get control diet group specific DEGs
dup = markers[duplicated(markers$gene),]
SCFA_DEG = markers[markers$cluster %in% c("injury_control"),]
SCFA_DEG = setdiff(SCFA_DEG$gene,dup$gene)

go <- enrichGO(gene = SCFA_DEG, OrgDb = "org.Mm.eg.db", ont="BP",qvalueCutoff = 0.05,keyType = "SYMBOL",pvalueCutoff = 0.05)

edox2 <- pairwise_termsim(go)
treeplot(edox2, hclust_method = "average")
ggsave(paste0("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_injury_control_group_DEG_GO_treeplot.pdf"),width = 20,height = 8)

# find markers for each group
MB.Myeloid@active.ident = as.factor(MB.Myeloid$subtype)
markers =  FindAllMarkers(MB.Myeloid,
                          only.pos = T,
                          test.use = "MAST",
                          logfc.threshold = 0.25)

top5 <- markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
DoHeatmap(MB.Myeloid, features = top5$gene) + NoLegend()
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_subtype_marker_heatmap.pdf",width = 15,height = 8)
save(markers, file = "/ix/gkohanbash/SimonsLab/R_Share/processing/10_Myeloid_subtype_markers.RData")

# cell cycle

## mgi - hsa gene symbol transformation
require("biomaRt")
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

## get the cell cycle marker
### marker from: Lundberg A, Lindström LS, Parker JS, Löverli E, Perou CM, Bergh J, Tobin NP. A pan-cancer analysis of the frequency of DNA alterations across cell cycle activity levels. Oncogene. 2020 Aug;39(32):5430-5440. doi: 10.1038/s41388-020-1367-4. Epub 2020 Jun 24. PMID: 32581248; PMCID: PMC8159764.

g1.features = read.table("/ix/gkohanbash/SimonsLab/R_Share/marker/cell_cycle/G1.txt")
g1s.features = read.table("/ix/gkohanbash/SimonsLab/R_Share/marker/cell_cycle/G1-S.txt")
s.features = read.table("/ix/gkohanbash/SimonsLab/R_Share/marker/cell_cycle/S.txt")
g2.features = read.table("/ix/gkohanbash/SimonsLab/R_Share/marker/cell_cycle/G2.txt")
g2m.features = read.table("/ix/gkohanbash/SimonsLab/R_Share/marker/cell_cycle/G2-M.txt")
m.features = read.table("/ix/gkohanbash/SimonsLab/R_Share/marker/cell_cycle/M.txt")

g1.features = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                     values = g1.features$V1, 
                     mart = mouse, 
                     attributesL = c("hgnc_symbol"), 
                     martL = human, uniqueRows=T)
g1.features = g1.features$MGI.symbol

g1s.features = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                      values = g1s.features$V1, 
                      mart = mouse, 
                      attributesL = c("hgnc_symbol"), 
                      martL = human, uniqueRows=T)
g1s.features = g1s.features$MGI.symbol

s.features = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                    values = s.features$V1, 
                    mart = mouse, 
                    attributesL = c("hgnc_symbol"), 
                    martL = human, uniqueRows=T)
s.features = s.features$MGI.symbol

g2.features = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                     values = g2.features$V1, 
                     mart = mouse, 
                     attributesL = c("hgnc_symbol"), 
                     martL = human, uniqueRows=T)
g2.features = g2.features$MGI.symbol

g2m.features = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                      values = g2m.features$V1, 
                      mart = mouse, 
                      attributesL = c("hgnc_symbol"), 
                      martL = human, uniqueRows=T)
g2m.features = g2m.features$MGI.symbol

m.features = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                    values = m.features$V1, 
                    mart = mouse, 
                    attributesL = c("hgnc_symbol"), 
                    martL = human, uniqueRows=T)
m.features = m.features$MGI.symbol

## calculate the cycle score
name <- "Cell.Cycle"
features <- list(G1.Score = g1.features, 
                 G1S.Score = g1s.features,
                 S.Score = s.features,
                 G2.Score = g2.features,
                 G2M.Score = g2m.features,
                 M.Score = m.features)

object.cc <- AddModuleScore(object = MB.Myeloid, features = features, 
                            name = name, ctrl = 100)
cc.columns <- grep(pattern = name, x = colnames(x = object.cc[[]]), 
                   value = TRUE)
cc.scores <- object.cc[[cc.columns]]
# rm(object.cc)
CheckGC()
assignments <- apply(cc.scores, 1,function(scores, first = "G1",second = "G1S", third = "S",fourth = "G2", fifth = "G2M",sixth = "M")
{c(first, second,third,fourth,fifth,sixth)[which(x = scores == max(scores))]
})

cc.scores <- merge(x = cc.scores, y = data.frame(assignments),by = 0)
colnames(x = cc.scores) <- c("rownames", "G1.Score", "G1S.Score","S.Score", "G2.Score", "G2M.Score", "M.Score", 
                             "Phase")
rownames(x = cc.scores) <- cc.scores$rownames
cc.scores <- cc.scores[, c("G1.Score", "G1S.Score","S.Score", "G2.Score", "G2M.Score", "M.Score","Phase")]
MB.Myeloid[[colnames(x = cc.scores)]] <- cc.scores

UMAPPlot(MB.Myeloid,group.by = "Phase")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cell_cycle_phase_UMAP.pdf",width = 5,height = 4)

save(MB.Myeloid, file = "/ix/gkohanbash/SimonsLab/R_Share/processing/10_Myeloid_Seurat_object.RData")

a = table(MB.Myeloid$Phase,MB.Myeloid$subtype)
bar_input = as.data.frame(a)
colnames(bar_input) = c("phase","subtype","Freq")

col = c(RColorBrewer::brewer.pal(9,'Set1'),
        RColorBrewer::brewer.pal(6,'Pastel1'),
        RColorBrewer::brewer.pal(10,'Set3'))
bar_input %>%
  ggplot(aes(x = subtype, y = Freq, fill = phase))+
  geom_bar(stat="identity", position = 'fill',width = 0.5)+
  scale_fill_manual(values= col)
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cell_cycle_phase_subtype_frequency_barplot.pdf",width = 10,height = 8)

a = table(MB.Myeloid$Phase,MB.Myeloid$group)
bar_input = as.data.frame(a)
colnames(bar_input) = c("phase","group","Freq")

col = c(RColorBrewer::brewer.pal(9,'Set1'),
        RColorBrewer::brewer.pal(6,'Pastel1'),
        RColorBrewer::brewer.pal(10,'Set3'))
bar_input %>%
  ggplot(aes(x = group, y = Freq, fill = phase))+
  geom_bar(stat="identity", position = 'fill',width = 0.5)+
  scale_fill_manual(values= col)
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cell_cycle_phase_group_frequency_barplot.pdf",width = 10,height = 8)


## pairwise chi-square test in cluster frequency for different groups
a = as.vector(a)
dim(a) = c(length(unique(MB.Myeloid$Phase)),3)
a = data.frame(t(a))
rownames(a) = c("injury_control","injury_SCFA","sham_control")
colnames(a) = c("G1","G1S","G2","G2M","M","S")

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
write.csv(chi_p,file = "/ix/gkohanbash/SimonsLab/R_Share/table/No_acetate_group/10_Myeloid_cell_cycle_phase_chi_square_p_value_matrix.csv")

chi_p[chi_p < 1e-4] = 4
chi_p[chi_p > 0.05 & chi_p <= 1 | chi_p == 'NaN'] = 0
chi_p[chi_p > 0.01 & chi_p < 0.05] = 1
chi_p[chi_p > 0.001 & chi_p < 0.01] = 2
chi_p[chi_p > 0.0001 & chi_p < 0.001] = 3

color = RColorBrewer::brewer.pal(9,"Reds")[c(1,4,6,7,9)]
pheatmap::pheatmap(chi_p,cellwidth = 10,cellheight = 10,
                   cluster_rows = F,cluster_cols = F,
                   color = color)
# save as "/ix/gkohanbash/SimonsLab/R_Share/fig/10_Myeloid_group_pairwise_group_phase_frequency_chi_square.pdf" 4*6

# identification of neurodegenerative disease-associated microglia (DAM)
load("/ix/gkohanbash/SimonsLab/R_Share/processing/10_Myeloid_Seurat_object.RData")
load("/ix/gkohanbash/SimonsLab/R_Share/processing/10_Myeloid_subtype_markers.RData")

# get homestatic, early stage, late stage MG specific marker from ref
# cite: Mathys H, Adaikkan C, Gao F, Young JZ, Manet E, Hemberg M, De Jager PL, Ransohoff RM, Regev A, Tsai LH. Temporal Tracking of Microglia Activation in Neurodegeneration at Single-Cell Resolution. Cell Rep. 2017 Oct 10;21(2):366-380. doi: 10.1016/j.celrep.2017.09.039. PMID: 29020624; PMCID: PMC5642107.
# LATE-STAGE Here = DAM in  Keren-Shaul H, Spinrad A, Weiner A, Matcovitch-Natan O, Dvir-Szternfeld R, Ulland TK, David E, Baruch K, Lara-Astaiso D, Toth B, Itzkovitz S, Colonna M, Schwartz M, Amit I. A Unique Microglia Type Associated with Restricting Development of Alzheimer's Disease. Cell. 2017 Jun 15;169(7):1276-1290.e17. doi: 10.1016/j.cell.2017.05.018. Epub 2017 Jun 8. PMID: 28602351.

a1 = read.csv("/ix/gkohanbash/SimonsLab/R_Share/marker/Neurodegenerative_MG/Cell_Rep/Homeostatic_vs_Early_stage_1.csv")
a2 = read.csv("/ix/gkohanbash/SimonsLab/R_Share/marker/Neurodegenerative_MG/Cell_Rep/Homeostatic_vs_Early_stage_2.csv")
a3 = read.csv("/ix/gkohanbash/SimonsLab/R_Share/marker/Neurodegenerative_MG/Cell_Rep/Homeostatic_vs_Late_stage.csv")
a4 = read.csv("/ix/gkohanbash/SimonsLab/R_Share/marker/Neurodegenerative_MG/Cell_Rep/Early_stage_s_Late_stage.csv")

a1 = a1[abs(a1$mle) > 2,]
a2 = a2[abs(a2$mle) > 2,]
a3 = a3[abs(a3$mle) > 2,]
a4 = a4[abs(a4$mle) > 2,]

a1 = a1[order(a1$mle,decreasing = T),]
a2 = a2[order(a2$mle,decreasing = T),]
a3 = a3[order(a3$mle,decreasing = T),]

Early = intersect(setdiff(union(a1$X[a1$mle < 0],
                                a2$X[a2$mle < 0]),a3$X[a3$mle < 0]),a4$X[a4$mle > 0])
Late = intersect(a3$X[a3$mle < 0],a4$X[a4$mle < 0])

k=150
Homeostatic = intersect(a1$X[1:k],a3$X[1:k])

# input -- genes/features as rows and cells as columns.
exprMatrix <- MB.Myeloid@assays$SCT@counts # counts
exprMatrix <- as.matrix(exprMatrix)
cellsUMAP = MB.Myeloid@reductions$umap@cell.embeddings

# GO
## build GO ref geneset -- gene name, not entrenz 
ref = list(Homeostatic = Homeostatic,
           Early_stage = Early,
           Late_stage = Late)

# ranking cells by gene expression high -> low
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=nCores, plotStats=F)

#Calculate enrichment for the gene signatures (AUC)
cells_AUC <- AUCell_calcAUC(ref, cells_rankings,
                            aucMaxRank=nrow(cells_rankings)*0.05)#The percentage to take into account can be modified with the argument aucMaxRank
#by default only the top 5% of the genes in the ranking are used
#the AUC estimates the proportion of genes in the gene-set that are highly expressed in each cell.
#Cells expressing many genes from the gene-set will have higher AUC values than cells expressing fewer

#Determine the cells with the given gene signatures or active gene sets - find threshold
set.seed(123)
save(cells_AUC, file = "/ix/gkohanbash/SimonsLab/R_Share/processing/10_Myeloid_DAM_AUC.RData")

Myeloid_DAM_AUC = cells_AUC@assays@data@listData[["AUC"]]

input = data.frame()

group = unique(MB.Myeloid$group)
subtype =  c("MG_1","MG_2","MG_3","MG_4","MG_5","MG_6","MG_7","MG_8","MG_9","MG_10")

for (i in 1:length(group)) {
  for (j in 1:length(subtype)) {
    cells = colnames(MB.Myeloid)[MB.Myeloid$group == group[i] & MB.Myeloid$subtype == subtype[j]]
    
    input[3*i - 2,j] = mean(Myeloid_DAM_AUC[1,cells])
    input[3*i - 1,j] = mean(Myeloid_DAM_AUC[2,cells])
    input[3*i ,j] = mean(Myeloid_DAM_AUC[3,cells])
  }
}

rownames(input) = c("Injury_SCFA_Homeostatic","Injury_SCFA_DAM_early","Injury_SCFA_DAM_late",
                    "Injury_control_Homeostatic","Injury_control_DAM_early","Injury_control_DAM_late",
                    "Sham_control_Homeostatic","Sham_control_DAM_early","Sham_control_DAM_late")
colnames(input) = subtype

# for young sample
load("/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/10_MG_Ruchi_data_matching_Seurat.RData")
Myeloid$subtype = Myeloid$predicted.cellType

exprMatrix <- Myeloid@assays$RNA@counts # counts
exprMatrix <- as.matrix(exprMatrix)
cellsUMAP = Myeloid@reductions$umap@cell.embeddings

# GO
## build GO ref geneset -- gene name, not entrenz 
ref = list(Homeostatic = Homeostatic,
           Early_stage = Early,
           Late_stage = Late)

# ranking cells by gene expression high -> low
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=nCores, plotStats=F)

#Calculate enrichment for the gene signatures (AUC)
cells_AUC <- AUCell_calcAUC(ref, cells_rankings,
                            aucMaxRank=nrow(cells_rankings)*0.05)#The percentage to take into account can be modified with the argument aucMaxRank

#Determine the cells with the given gene signatures or active gene sets - find threshold
set.seed(123)
save(cells_AUC, file = "/ix/gkohanbash/SimonsLab/R_Share/processing/10_Myeloid_DAM_AUC_of_Ruchi_condition_data.RData")

Myeloid_DAM_AUC_R = cells_AUC@assays@data@listData[["AUC"]]

input_r = data.frame()

subtype =  c("MG_1","MG_2","MG_3","MG_4","MG_5","MG_6","MG_7","MG_8","MG_9","MG_10")

for (j in 1:length(subtype)) {
  cells = colnames(Myeloid)[Myeloid$condition == "Naive" & Myeloid$subtype == subtype[j]]
  
  input_r[1,j] = mean(Myeloid_DAM_AUC_R[1,cells])
  input_r[2,j] = mean(Myeloid_DAM_AUC_R[2,cells])
  input_r[3,j] = mean(Myeloid_DAM_AUC_R[3,cells])
}

rownames(input_r) = c("Young_Normal_Homeostatic","Young_Normal_DAM_early","Young_Normal_DAM_late")
colnames(input_r) = subtype

input = rbind(input,input_r)

# p value comparison - Mann-Whitney U Test
p_value = data.frame()

for(i in 1:length(subtype)) {
  cell1 = colnames(MB.Myeloid)[MB.Myeloid$subtype == subtype[i] & MB.Myeloid$group == "sham_control"]
  cell2 = colnames(Myeloid)[Myeloid$subtype == subtype[i] & Myeloid$condition == "Naive"]
  
  p1 = wilcox.test(Myeloid_DAM_AUC[1,cell1],Myeloid_DAM_AUC_R[1,cell2])
  p2 = wilcox.test(Myeloid_DAM_AUC[2,cell1],Myeloid_DAM_AUC_R[2,cell2])
  p3 = wilcox.test(Myeloid_DAM_AUC[3,cell1],Myeloid_DAM_AUC_R[3,cell2])
  
  p_value[7,i] = p1$p.value
  p_value[8,i] = p2$p.value
  p_value[9,i] = p3$p.value
  
  cell3 = colnames(MB.Myeloid)[MB.Myeloid$subtype == subtype[i] & MB.Myeloid$group == "injury_control"]
  cell4 = colnames(MB.Myeloid)[MB.Myeloid$subtype == subtype[i] & MB.Myeloid$group == "injury_SCFA"]
  
  p1 = wilcox.test(Myeloid_DAM_AUC[1,cell3],Myeloid_DAM_AUC[1,cell4])
  p2 = wilcox.test(Myeloid_DAM_AUC[2,cell3],Myeloid_DAM_AUC[2,cell4])
  p3 = wilcox.test(Myeloid_DAM_AUC[3,cell3],Myeloid_DAM_AUC[3,cell4])
  
  p_value[1,i] = p1$p.value
  p_value[2,i] = p2$p.value
  p_value[3,i] = p3$p.value
}

p_value[10:12,] = NA

for (r in 1:nrow(p_value)) {
  for (c in 1:ncol(p_value)) {
    if (is.na(p_value[r,c])) {
      p_value[r,c] = ""
    }else{
      
      if(r < 4) {
        p_value[r,c] = if_else(as.numeric(p_value[r,c]) < 0.05,"*","")
      }
      else{
        p_value[r,c] = if_else(as.numeric(p_value[r,c])  < 0.05,"#","")
      }}
  }
}

# annotation
# left
g = data.frame(sample = rownames(input),
               group = c(rep("Injury_SCFA",3),rep("Injury_control",3),rep("Sham_control",3),rep("Normal_young",3)))
group = g[,2]
names(group) = g$sample

color.use <- c(hue_pal()(3),"white")
names(color.use) = c("Sham_control","Injury_control","Injury_SCFA","Normal_young")

# colors = colors[unique(left_anno)]
left_annotation <- HeatmapAnnotation(group = group, col = list(group = color.use), 
                                     which = "row", show_legend = T, show_annotation_name = T, 
                                     simple_anno_size = grid::unit(0.4, "cm"))

# legend
lgd_sig = Legend(pch = "*", type = "points", labels = "sig in Injury_SCFA vs Injury_control")
lgd_sig_1 = Legend(pch = "#", type = "points",labels = "sig in Sham_control vs Young_Normal")

# heatmap
library(ComplexHeatmap)

color.heatmap.use = colorRampPalette(colors = c("deepskyblue3", "white", "tomato3"))(100)

if (min(input, na.rm = T) == max(input, na.rm = T)) {
  legend.break <- max(input, na.rm = T)
}else {
  legend.break <- c(round(min(input, na.rm = T), digits = 1), 
                    round(max(input, na.rm = T), digits = 1))
}

ht = Heatmap(input, col = color.heatmap.use, row_split = factor(g$group,levels = c("Normal_young","Sham_control","Injury_control","Injury_SCFA")),
             name = "Mean AUC", left_annotation = left_annotation,#bottom_annotation = col_annotation, top_annotation = top_annotation,
             cluster_rows = F, 
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.text(p_value[i, j], x, y, gp = gpar(fontsize = 12))},
             cluster_columns = F, row_names_side = "left", 
             row_names_rot = 0, row_names_gp = gpar(fontsize = 8), 
             column_names_gp = gpar(fontsize = 8), width = unit(10,"cm"), 
             height = unit(8, "cm"), column_title = "Neurodegenerative MG signature score", 
             column_title_gp = gpar(fontsize = 10), 
             column_names_rot = 90, heatmap_legend_param = list(title_gp = gpar(fontsize = 8,fontface = "plain"), title_position = "leftcenter-rot", 
                                                                border = NA, at = legend.break, legend_height = unit(20,"mm"), labels_gp = gpar(fontsize = 8), grid_width = unit(2,"mm")))

draw(ht,  annotation_legend_list = list(lgd_sig,lgd_sig_1))
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/long-term_TBI_consquence/10_Myeloid_DAM_signature_AUC_comparison.pdf 10*20  


# validate
# cite: Keren-Shaul H, Spinrad A, Weiner A, Matcovitch-Natan O, Dvir-Szternfeld R, Ulland TK, David E, Baruch K, Lara-Astaiso D, Toth B, Itzkovitz S, Colonna M, Schwartz M, Amit I. A Unique Microglia Type Associated with Restricting Development of Alzheimer's Disease. Cell. 2017 Jun 15;169(7):1276-1290.e17. doi: 10.1016/j.cell.2017.05.018. Epub 2017 Jun 8. PMID: 28602351.
MB.MG = MB.Myeloid[,MB.Myeloid$subtype != "MAC"]

MB.MG@active.ident = as.factor(MB.MG$subtype)
markers =  FindAllMarkers(MB.MG,
                          only.pos = T,
                          test.use = "MAST",
                          logfc.threshold = 0.25)


DAM = readxl::read_xlsx("/ix/gkohanbash/SimonsLab/R_Share/marker/Neurodegenerative_MG/mmc2.xlsx")

DAM$diff = log2((DAM$`Microglia3  (average UMI count)`+1)/(1+DAM$`Microglia1 (average UMI count)`))

# filter marker
DAM = DAM[-grep("mir",DAM$...1),]
DAM = DAM[order(DAM$diff,decreasing = T),]

Homeostatic = DAM$...1[DAM$`up/down` == -1]

k = 50
subtype_ref = data.frame(Homeostatic = Homeostatic[1:k], # only has 27 genes
                         DAM = DAM$...1[1:k])

# subtype_ref = data.frame(Homeostatic = c("Hexb","Cst3","Cx3cr1","Ctsd","Csf1r","Ctss","Sparc","Tmsb4x","P2ry12","C1qa","C1qb",rep(NA,1)),
#                          Stage1_DAM = c("Tyrobp","Ctsb","Ctsd","Apoe","B2m","Fth1","Lyz2",rep(NA,5)),
#                          Stage2_DAM = c("Trem2","Axl","Cst7","Ctsl","Lpl","Cd9","Csf1","Ccl6","Itgax","Clec7a","Lilrb4","Timp2"))

# count cluster markers and ref markers overlap and evaluate by Fisher exact test
load( "/ix/gkohanbash/SimonsLab/R_Share/processing/10_Myeloid_subtype_markers.RData")
markers = markers[markers$cluster != "MAC",]

overlap_count = data.frame()
overlap_fdr = data.frame()
overlap_p = vector(mode = "numeric")

cluster = unique(markers$cluster)

for (i in 1:ncol(subtype_ref)) {
  for (j in 1:length(cluster)) {
    ref = na.omit(subtype_ref[,i])
    test = markers$gene[which(markers$cluster == cluster[j])]
    o = length(intersect(ref,test))    # the number of marker genes that in ref gene list 
    k = length(ref)            # ref gene count
    n = nrow(MB.Myeloid)      # total gene count
    m = length(test)           # the number of input DEGs
    overlap_p = c(overlap_p,1-phyper(o-1,k,n-k,m))  # Fisher exact test
    overlap_count[j,i] = o
  }
}

rownames(overlap_count) = cluster#paste0("cluster_",cluster)
colnames(overlap_count) = colnames(subtype_ref)

overlap_fdr = p.adjust(overlap_p, method = "BH", n = length(overlap_p))
# dim(overlap_fdr) = c(29,8)
# rownames(overlap_fdr) = paste0("cluster_",cluster)
# colnames(overlap_fdr) = colnames(cell_marker_ref)

overlap_count$cluster = rownames(overlap_count)
hm <- overlap_count %>% gather(Celltype, value, 1:ncol(subtype_ref))
overlap_fdr = if_else(overlap_fdr > 0.05,0,
                      if_else(overlap_fdr > 0.01, 1,
                              if_else(overlap_fdr > 0.001, 2,
                                      if_else(overlap_fdr > 1e-04,3,
                                              if_else(overlap_fdr > 1e-05,4,
                                                      if_else(overlap_fdr > 1e-6,5,6))))))

# buuild confusion matrix input
hm$fdr = overlap_fdr

# plot the confusion matrix
ggplot(hm, aes(x=cluster, y=Celltype, fill=fdr)) +
  geom_tile() + theme_bw() + coord_equal() +
  scale_fill_distiller(palette="Greens", direction=1) +
  #guides(fill=F) + # removing legend for `fill`
  labs(title ="DEGs overlap") + # using a title instead
  geom_text(aes(label=value), color="black") # printing values
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_Neurodegenerative_associated_MG_identification_by_markers.pdf")
# color = FDR value
# number in cell = number of cluster marker gene overlaped with ref cell markers
# FDR: 0 = fdr -> (0.05,1]
#      1 = fdr -> (0.01,0.05]
#      2 = fdr -> (1e-03,1e-02]
#      3 = fdr -> (1e-04,1e-03]
#      4 = fdr -> (1e-05,1e-04]
#      5 = fdr -> (1e-06,1e-05]
#      6 = fdr -> [0,1e-06]

# compare DAM MG activation between groups
# cite: Mathys H, Adaikkan C, Gao F, Young JZ, Manet E, Hemberg M, De Jager PL, Ransohoff RM, Regev A, Tsai LH. Temporal Tracking of Microglia Activation in Neurodegeneration at Single-Cell Resolution. Cell Rep. 2017 Oct 10;21(2):366-380. doi: 10.1016/j.celrep.2017.09.039. PMID: 29020624; PMCID: PMC5642107.
exprMatrix <- MB.Myeloid@assays$SCT@counts # counts
exprMatrix <- as.matrix(exprMatrix)
cellsUMAP = MB.Myeloid@reductions$umap@cell.embeddings
logMat <- log2(exprMatrix+2)

# GO
## build GO ref geneset -- gene name, not entrenz 
load("/ix/gkohanbash/SimonsLab/R_Share/processing/5_GO_AUC_genelist.RData")

## Neurodegenerative disease
geneSetName = c("GOBP_MICROGLIAL_CELL_ACTIVATION","GOBP_POSITIVE_REGULATION_OF_MICROGLIAL_CELL_ACTIVATION")

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)

boxinput = data.frame(Score = as.numeric(meanByGs[1,]),
                      subtype = MB.Myeloid$subtype,
                      group = factor(MB.Myeloid$group,level = c("sham_control","injury_control","injury_SCFA")))

boxinput = boxinput[boxinput$subtype != "MAC",]

dodge <- position_dodge(width = 0.7)
boxinput%>%ggplot(aes(subtype,Score))+
  geom_boxplot(aes(fill=group),position = dodge,outlier.size = 0,notch = F,width = 0.5)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "right"
  )+
  stat_compare_means(aes(group = group),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),label = "p.signif",method = "anova")+
  ggtitle("Microglia Activation")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/long-term_TBI_consquence/10_MG_activation_split_by_subtype_ANOVA.pdf")

# SEPARATE FOR MG SUBTYPE: late_stage_DAM: MG_5, MG_6, MG_7 and injury_control
vln.df = boxinput[boxinput$subtype == "MG_5",]

group=levels(factor(vln.df$group))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

dodge <- position_dodge(width = 1)
vln.df%>%ggplot(aes(group,Score))+geom_violin(aes(fill=group),scale = "width",position = dodge)+
  geom_boxplot(aes(fill=group),position = dodge,outlier.size = 0,notch = T,width = 0.3)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "right"
  )+
  stat_compare_means(comparisons = my_comparisons,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),label = "p.signif")+
  ggtitle("MG_5_GOBP_MICROGLIAL_ACTIVATION")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/long-term_TBI_consquence/10_MG_5_activation_score_comparison.pdf",width = 5,height = 5)

vln.df = boxinput[boxinput$subtype == "MG_6",]

group=levels(factor(vln.df$group))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

dodge <- position_dodge(width = 1)
vln.df%>%ggplot(aes(group,Score))+geom_violin(aes(fill=group),scale = "width",position = dodge)+
  geom_boxplot(aes(fill=group),position = dodge,outlier.size = 0,notch = T,width = 0.3)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "right"
  )+
  stat_compare_means(comparisons = my_comparisons,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),label = "p.signif")+
  ggtitle("MG_6_GOBP_MICROGLIAL_ACTIVATION")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/long-term_TBI_consquence/10_MG_6_activation_score_comparison.pdf",width = 5,height = 5)

vln.df = boxinput[boxinput$subtype == "MG_7",]

group=levels(factor(vln.df$group))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

dodge <- position_dodge(width = 1)
vln.df%>%ggplot(aes(group,Score))+geom_violin(aes(fill=group),scale = "width",position = dodge)+
  geom_boxplot(aes(fill=group),position = dodge,outlier.size = 0,notch = T,width = 0.3)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "right"
  )+
  stat_compare_means(comparisons = my_comparisons,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),label = "p.signif")+
  ggtitle("MG_7_GOBP_MICROGLIAL_ACTIVATION")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/long-term_TBI_consquence/10_MG_7_activation_score_comparison.pdf",width = 5,height = 5)

vln.df = boxinput[boxinput$subtype == "MG_1",]

group=levels(factor(vln.df$group))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

dodge <- position_dodge(width = 1)
vln.df%>%ggplot(aes(group,Score))+geom_violin(aes(fill=group),scale = "width",position = dodge)+
  geom_boxplot(aes(fill=group),position = dodge,outlier.size = 0,notch = T,width = 0.3)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "right"
  )+
  stat_compare_means(comparisons = my_comparisons,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),label = "p.signif")+
  ggtitle("MG_1_GOBP_MICROGLIAL_ACTIVATION")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/long-term_TBI_consquence/10_MG_1_activation_score_comparison.pdf",width = 5,height = 5)

vln.df = boxinput[boxinput$subtype == "MG_2",]

group=levels(factor(vln.df$group))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

dodge <- position_dodge(width = 1)
vln.df%>%ggplot(aes(group,Score))+geom_violin(aes(fill=group),scale = "width",position = dodge)+
  geom_boxplot(aes(fill=group),position = dodge,outlier.size = 0,notch = T,width = 0.3)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "right"
  )+
  stat_compare_means(comparisons = my_comparisons,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),label = "p.signif")+
  ggtitle("MG_2_GOBP_MICROGLIAL_ACTIVATION")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/long-term_TBI_consquence/10_MG_2_activation_score_comparison.pdf",width = 5,height = 5)

vln.df = boxinput[boxinput$subtype == "MG_3",]

group=levels(factor(vln.df$group))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

dodge <- position_dodge(width = 1)
vln.df%>%ggplot(aes(group,Score))+geom_violin(aes(fill=group),scale = "width",position = dodge)+
  geom_boxplot(aes(fill=group),position = dodge,outlier.size = 0,notch = T,width = 0.3)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "right"
  )+
  stat_compare_means(comparisons = my_comparisons,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),label = "p.signif")+
  ggtitle("MG_3_GOBP_MICROGLIAL_ACTIVATION")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/long-term_TBI_consquence/10_MG_3_activation_score_comparison.pdf",width = 5,height = 5)

# axon growth inhibitory molecules secretion
gene = c("Acan","Bcan","Ncan","Vcan","Xylt1","Tnr","Epha4","Ephb2","Efnb3","Ntn1","Rgma","Rgmb","Slit1","Slit2","Slitrk1","Ptprz1","Sema3a",
         "Sema3f","Plxna1","Plxnb1","Nrp1","Unc5b","Dcc","Neo1","Robo1","Robo2","Robo3","Draxin")

# Reactive astrocytes are the main cellular component of the glial scar. -- ASC_9 & ASC_10
DotPlot(MB.Myeloid, features = gene, group.by = "subtype")+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_subtype_axon_growth_inhibitory_molecules_dotplot.pdf",height = 3,width = 9)

DotPlot(MB.Myeloid[,MB.Myeloid$subtype == "MAC"], features = gene, group.by = "group")+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_MAC_axon_growth_inhibitory_molecules_dotplot.pdf",height = 3,width = 9)

# axon permissive molecules secretion
gene = c("Cspg4","Cspg5","Tnc","Sdc1","Sdc2","Sdc3","Sdc4","Bdnf","Lif","Igf1","Fgf2","Tgfa","Lama1","Lama2",
         "Lama4","Lama5","Lamb1","Lamc1","Col4a1","Fn1","Hspg2","Gpc1","Gpc3","Gpc5","Dcn","Lgals1","Ncam1","Matn2")
# [1] "Ntf3" "Gdnf" "Cntf" not shown in the dataset

# Reactive astrocytes are the main cellular component of the glial scar. -- ASC_9 & ASC_10
DotPlot(MB.Myeloid, features = gene, group.by = "subtype")+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_subtype_axon_permissive_molecules_dotplot.pdf",height = 3,width = 9)

DotPlot(MB.Myeloid[,MB.Myeloid$subtype == "MAC"], features = gene, group.by = "group")+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_MAC_axon_permissive_molecules_dotplot.pdf",height = 3,width = 9)

# pro-inflammatory cytokine & chemokine in Myeloid cell
gene = c("Il1a","Il1b","Il6","Tnf","Ccl2")

# MG_5 & MG_7 are two main cell subtype secret cytokines
DotPlot(MB.Myeloid, features = gene, group.by = "subtype")+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_subtype_pro_inflammatory_molecules_dotplot.pdf",height = 3,width = 9)

DotPlot(MB.Myeloid[,MB.Myeloid$subtype == "MG_5"], features = gene, group.by = "group")+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_MG_5_pro_inflammatory_molecules_dotplot.pdf",height = 3,width = 9)

DotPlot(MB.Myeloid[,MB.Myeloid$subtype == "MG_7"], features = gene, group.by = "group")+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_MG_7_pro_inflammatory_molecules_dotplot.pdf",height = 3,width = 9)

DotPlot(MB.Myeloid[,MB.Myeloid$subtype == "MAC"], features = gene, group.by = "group")+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_MAC_pro_inflammatory_molecules_dotplot.pdf",height = 3,width = 9)

# scWGCNA to show different module expression in groups
k = 100    # the number of cells used to create the pseudocell/metacells
seurat_object = MB.Myeloid

## select the feature that you want to use to build pseudo-metacells
seurat_object$metacell_group = MB.Myeloid$subtype
celltype = as.data.frame(table(seurat_object$metacell_group))    # change cell type
cell = character()
for (i in unique(celltype$Var1)){
  if(min(celltype[celltype$Var1 == i,]$Freq) > k*1.2){
    cell = c(cell,i)
  }else{
    cell = cell
  }
}
seurat_object = seurat_object[,seurat_object$metacell_group %in% cell]                # change cell type
genes.keep = VariableFeatures(seurat_object)

# construct metacells
seurat_list = list()
for (group in unique(seurat_object$metacell_group)) {
  print(group)
  cur_seurat = seurat_object[,seurat_object$metacell_group == group]
  cur_seurat = cur_seurat[genes.keep,]
  cur_metacell_seurat = scWGCNA::construct_metacells(
    cur_seurat, name = group,
    k = k, reduction = "umap",
    assay = "RNA", slot = "counts"                                        # change assay
  )
  cur_metacell_seurat$group = as.character(group)          # change condition
  seurat_list[[group]] = cur_metacell_seurat
}

metacell_seurat = merge(seurat_list[[1]],seurat_list[2:length(seurat_list)])
save(metacell_seurat,
     file = "/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/10_Myeloid_gene_expression_count_metacell.RData")

# visualize metacell
metacell_seurat = SCTransform(metacell_seurat,
                              ncells = 10000)
metacell_seurat = RunPCA(metacell_seurat, features = rownames(metacell_seurat))
metacell_seurat = RunUMAP(metacell_seurat , reduction = "pca", dims = 1:20)

DimPlot(metacell_seurat, group.by = "group", reduction = "umap",label = T) + NoLegend()   # change cell type
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_gene_metacell_UMAP.pdf",width = 5,height = 5)

#scWGCNA
### additional variable 
dataset = as.data.frame(metacell_seurat@assays$SCT@scale.data)
metadata = metacell_seurat@meta.data
networkType = "signed" # c("signed","unsigned")
feature = metadata$group # 进行对比的feature, the target feature, i.e metadata$gender
corType = "pearson"  # "pearson"  or  "bicor" 

###
maxPOutliers = ifelse(corType == "pearson", 1, 0.05) # 对二元变量，如样本性状信息计算相关性时， 或基因表达严重依赖于疾病状态时，需设置下面参数

# 关联样品性状的二元变量时，设置
robustY = ifelse(corType == "pearson", T, F)

## use all 3000 variable genes as input
mydata = t(dataset)

#filter sample and gene
gsg = goodSamplesGenes(mydata, verbose = 3)
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(mydata)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(mydata)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  mydata = mydata[gsg$goodSamples, gsg$goodGenes]
}

#find outlier
sampleTree = hclust(dist(mydata), method = "average")
save(sampleTree, file = "/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/10_Myeloid_scWGCNA_sampleTree.RData")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree,
     main = "Sample clustering to detect outliers",
     sub="",
     xlab="", 
     cex.lab = 1.5,
     cex.axis = 1.5,
     cex.main = 2)

#remove outlier
h = (quantile(sampleTree[["height"]],probs = seq(0,1,0.1))[10]-quantile(sampleTree[["height"]],probs = seq(0,1,0.1))[2])*18 + quantile(sampleTree[["height"]],probs = seq(0,1,0.1))[10] 
# 超出 Q3 的 3 倍四分位间距（IQR）视为离群值

plot(sampleTree,
     main = "Sample clustering to detect outliers",
     sub="", 
     xlab="",
     cex.lab = 1.5, 
     cex.axis = 1.5, 
     cex.main = 2) +
  abline(h = h, col = "red") 
clust = cutreeStatic(sampleTree, cutHeight = h, minSize = 10) 
keepSamples = (clust==1)
mydata = mydata[keepSamples, ]
nGenes = ncol(mydata)
nSamples = nrow(mydata)
dim(mydata)

#metadata
metadata = metadata[keepSamples,]

#soft threshold power
#finding power R^2>0.9
powers = c(1:30)
sft = pickSoftThreshold(mydata, 
                        powerVector = powers,
                        networkType = networkType)

save(sft, file = "/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/10_Myeloid_scWGCNA_softThreshold.RData")

par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",
     type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#select soft power
power = sft$powerEstimate

if (is.na(power)){
  power = ifelse(nSamples<20, ifelse(networkType == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(networkType == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(networkType == "unsigned", 7, 14),
                               ifelse(networkType == "unsigned", 6, 12))       
                 )
  )
}

save(mydata, metadata, file = "/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/10_Myeloid_scWGCNA_input.RData")
#
cor = WGCNA::cor
net = blockwiseModules(mydata, power = power,
                       TOMType = "unsigned",  # ideally for single cell in scWGCNA pipeline, TOMType = "unsigned"
                       networkType = "signed", # ideally for single cell in scWGCNA pipeline, networkType = "signed"
                       corType = corType, #bicor for binary variate or outliers
                       mergeCutHeight = 0.15,
                       maxPOutliers = maxPOutliers,
                       minModuleSize = 20,
                       reassignThreshold = 0,
                       numericLabels = TRUE, 
                       pamStage = TRUE,
                       pamRespectsDendro = T,
                       saveTOMs = TRUE,
                       nThreads = 16)
save(net, file = "/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/10_Myeloid_scWGCNA_net.RData")

#dendrogram
mergedColors = labels2colors(net$colors)

plotDendroAndColors(net$dendrograms[[1]],
                    mergedColors[net$blockGenes[[1]]],
                    "Merged Dynamic Hybrid Tree Cut",
                    dendroLabels = FALSE, 
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

#annotate color bar for dendrogram
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
dynamicColors = mergedColors[net$blockGenes[[1]]]

################## 观察WGCNA的各个模块不同feature中具体表达 #############
library(caret)
e <- dummyVars(~as.factor(group),metadata)                                      # alter stage to the target feature ~as.factor(target feature)
f <- predict(e,metadata)

group.selected = model.matrix(~ as.factor(f[,1])) #生成哑变量
group.selected = as.data.frame(group.selected)
for(i in 1:length(unique(metadata$group))-1){
  group.selected[,i+2] = model.matrix(~ as.factor(f[,i+1]))[,2] #生成哑变量
}

trait.exp.cor.mat = cor(group.selected,as.data.frame(mydata)[,net[["blockGenes"]][[1]]],use = "p")
trait.exp.cor.mat[is.na(trait.exp.cor.mat)] = 0
library(fifer)
trait.exp.cor.mat.color = data.frame(dynamicColors)
for(i in 1:length(unique(metadata$group))){
  trait.exp.cor.mat.color[,i+1] = number.to.colors(trait.exp.cor.mat[i+1,], colors = c("green", "red"), num = 100)
}

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_scWGCNA_module.pdf",width = 10, height = 5)
plotDendroAndColors(geneTree, 
                    trait.exp.cor.mat.color,
                    c("Gene Module",sort(unique(metadata$group))),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#calculation module preservation
library(caret)
inTraining <- createDataPartition(metadata$group, p = 0.75, list = FALSE)
train<- mydata[inTraining,]
test<-mydata[-inTraining,]
setLabels = c("Train", "Test");
multiExpr = list(Train = list(data = train), Test = list(data = test));
multiColor = list(Train = moduleColors)

#2. preservation分析
mp = modulePreservation(multiExpr, multiColor,
                        referenceNetworks = 1,
                        nPermutations = 200, #重复计算数，正式分析官方推荐200 
                        randomSeed = 1,
                        quickCor = 0,
                        verbose = 3)
save(mp,file = "/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/10_Myeloid_mp-10%.RData") 


ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);
print(cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
            signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )
#Z.pres大于10，代表strong preserved，好的module
#大于2小于10代表weak preserved
#小于2代表not preserved，不好的module

modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1]

#这几个module的Z低于10
row.names(statsZ[statsZ$Zsummary.pres<10,])  #"cyan"       "lightgreen";  10% 都大于0 
# [1] "purple"

#去掉Z<10的module
plotMods = !(modColors %in% row.names(statsZ[statsZ$Zsummary.pres<10,]))
text = modColors[plotMods];
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])

#3.preservation的可视化
mains = c("Preservation Median rank", "Preservation Zsummary");

sizeGrWindow(10, 5);
pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_scWGCNA_module_preservation.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2){
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2){
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
  # For Zsummary, add threshold lines
  if (p==2){
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  }
}
dev.off()

#output modules
dir.create('/ix/gkohanbash/SimonsLab/R_Share/table/No_acetate_group/scWGCNA/Myeloid_scWGCNA_module_gene_list')
setwd('/ix/gkohanbash/SimonsLab/R_Share/table/No_acetate_group/scWGCNA/Myeloid_scWGCNA_module_gene_list')
mydata = as.data.frame(mydata)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),mydata[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"csv",sep = "."),quote = F)
}

######    output hub genes 
hub = chooseTopHubInEachModule(mydata, 
                               dynamicColors, 
                               power = power, 
                               type= networkType)

write.csv(hub, file = "/ix/gkohanbash/SimonsLab/R_Share/table/No_acetate_group/scWGCNA/Myeloid_scWGCNA_module_gene_list/module_hub_genes.csv",row.names = T,col.names = T)

# GO cluster markers for modules
text = text[! text %in% c("gold","grey")]
for(i in 1:length(text)){
  gene <- colnames(mydata)[moduleColors==text[i]]
  go <- enrichGO(gene = gene, OrgDb = "org.Mm.eg.db", ont="BP",qvalueCutoff = 0.05,keyType = "SYMBOL")
  
  p3 = cnetplot(go,layout = "circle")
  p4 = goplot(go,showCategory = 10)
  cowplot::plot_grid(p3,p4,rel_widths = c(1.5,2))
  ggsave(paste0("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_GO_Myeloid_module_",text[i],".pdf"),width = 16,height = 8)
}

View(go@result)
go@result = go@result[-grep("glio",go@result$Description),]

# use selected WGCNA module enriched item  to show the module function
load("/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/Extra_scWGCNA_pathway_group_by_GK.RData")
tab = tab.m
tab = as.data.frame(tab)
rownames(tab) = tab[,1]
tab = tab[,-1]

left_anno = tab[,ncol(tab)]
tab = tab[,-ncol(tab)]

library(ComplexHeatmap)

color.heatmap = "BuGn"
colnames(tab)
color.use <- c("black","blue","brown","green", "magenta","pink","purple","red","turquoise","yellow")
color.heatmap.use = colorRampPalette(colors = c("white", "red"))(100)

df <- data.frame(group = colnames(tab))
rownames(df) <- colnames(tab)
names(color.use) <- colnames(tab)
col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
                                    which = "column", show_legend = FALSE, show_annotation_name = FALSE, 
                                    simple_anno_size = grid::unit(0.2, "cm"))

df <- data.frame(group = left_anno)
rownames(df) <- rownames(tab)

dat = df
dat$module = c(rep(1,4),
               rep(2,7),
               rep(3,2),
               rep(4,5),
               rep(5,4),
               rep(6,2),
               rep(7,5),
               rep(8,4),
               rep(9,6),
               rep(10,5)) # the WGCNA module rank. rep(1,5) is black module has 5 pathway, 1 = black, 5 = 5 pathways
dat = dat[order(dat$module,dat$group),]

tab = tab[rownames(dat),]
df = as.data.frame(df[rownames(dat),])
colnames(df) = "group"
rownames(df) = rownames(tab)

# colors = brewer.pal(9,"Set1")
# names(colors) <- unique(left_anno)
colors = colors[unique(left_anno)]
left_annotation <- HeatmapAnnotation(df = df, col = list(group = colors), 
                                     which = "row", show_legend = T, show_annotation_name = T, 
                                     simple_anno_size = grid::unit(0.4, "cm"))

tab[is.na(tab)] = 0
if (min(tab, na.rm = T) == max(tab, na.rm = T)) {
  legend.break <- max(tab, na.rm = T)
}else {
  legend.break <- c(round(min(tab, na.rm = T), digits = 1), 
                    round(max(tab, na.rm = T), digits = 1))
}

Heatmap(tab, col = color.heatmap.use, 
        name = "WGCNA module GO significance", bottom_annotation = col_annotation, left_annotation = left_annotation,
        cluster_rows = F, 
        cluster_columns = F, row_names_side = "left", 
        row_names_rot = 0, row_names_gp = gpar(fontsize = 8), 
        column_names_gp = gpar(fontsize = 8), width = unit(10,"cm"), 
        height = unit(20, "cm"), column_title = "WGCNA gene module enriched GO pathway", 
        column_title_gp = gpar(fontsize = 10), 
        column_names_rot = 90, heatmap_legend_param = list(title_gp = gpar(fontsize = 8,fontface = "plain"), title_position = "leftcenter-rot", 
                                                           border = NA, at = legend.break, legend_height = unit(20,"mm"), labels_gp = gpar(fontsize = 8), grid_width = unit(2,"mm")))

# save 10_Myeloid_scWGCNA_module_pathway_heatmap.pdf    12*14
# p value significance
# 0: p = (0.05,1]
# 1: p = (0.01,0.05]
# 2: p = (0.001,0.01]
# 3: p = (0.0001,0.001]
# 4: p = (0.00001,0.0001]


# cluster group proportion
a = table(MB.Myeloid$group,MB.Myeloid$subtype)
bar_input = as.data.frame(a)
colnames(bar_input) = c("group","subtype","Freq")

col = c(RColorBrewer::brewer.pal(9,'Set1'),
        RColorBrewer::brewer.pal(6,'Pastel1'),
        RColorBrewer::brewer.pal(10,'Set3'))
bar_input %>%
  ggplot(aes(x = subtype, y = Freq, fill = group ))+
  geom_bar(stat="identity", position = 'fill')+
  scale_fill_manual(values= col)
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_subtype_frequency_barplot.pdf",width = 12,height = 8)

bar_input %>%
  ggplot(aes(x = group, y = Freq, fill = subtype))+
  geom_bar(stat="identity", position = 'fill',width = 0.5)+
  scale_fill_manual(values= col)
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_group_frequency_barplot.pdf",width = 10,height = 8)

## pairwise chi-square test in cluster frequency for different groups
id = colnames(a)
a = as.vector(a)
dim(a) = c(3,length(unique(MB.Myeloid$subtype)))
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
write.csv(chi_p,file = "/ix/gkohanbash/SimonsLab/R_Share/table/No_acetate_group/10_Myeloid_seurat_cluster_chi_square_p_value_matrix.csv")

chi_p[chi_p < 1e-4] = 4
chi_p[chi_p > 0.05 & chi_p <= 1 | chi_p == 'NaN'] = 0
chi_p[chi_p > 0.01 & chi_p < 0.05] = 1
chi_p[chi_p > 0.001 & chi_p < 0.01] = 2
chi_p[chi_p > 0.0001 & chi_p < 0.001] = 3

color = RColorBrewer::brewer.pal(9,"Reds")[c(1,4,6,7,9)]
pheatmap::pheatmap(chi_p,cellwidth = 10,cellheight = 10,
                   cluster_rows = F,cluster_cols = F,
                   color = color)
# save as "/ix/gkohanbash/SimonsLab/R_Share/fig/10_Myeloid_group_pairwise_cluster_frequency_chi_square.pdf" 4*6

# AUC
# input -- genes/features as rows and cells as columns.
exprMatrix <- MB.Myeloid@assays$SCT@counts # counts
exprMatrix <- as.matrix(exprMatrix)
cellsUMAP = MB.Myeloid@reductions$umap@cell.embeddings

# GO
## build GO ref geneset -- gene name, not entrenz 
load("/ix/gkohanbash/SimonsLab/R_Share/processing/5_GO_AUC_genelist.RData")

# ranking cells by gene expression high -> low
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=nCores, plotStats=F)

#Calculate enrichment for the gene signatures (AUC)
cells_AUC <- AUCell_calcAUC(h.sets, cells_rankings,
                            aucMaxRank=nrow(cells_rankings)*0.05)#The percentage to take into account can be modified with the argument aucMaxRank
#by default only the top 5% of the genes in the ranking are used
#the AUC estimates the proportion of genes in the gene-set that are highly expressed in each cell.
#Cells expressing many genes from the gene-set will have higher AUC values than cells expressing fewer

#Determine the cells with the given gene signatures or active gene sets - find threshold
set.seed(123)
save(cells_AUC, file = "/ix/gkohanbash/SimonsLab/R_Share/processing/10_Myeloid_GO_AUC.RData")

## select significant activity-different GO item between injury and sham group (containing SCFA feeding group)
Myeloid_GO_AUC = cells_AUC@assays@data@listData[["AUC"]]
Myeloid_GO_AUC = Myeloid_GO_AUC[,colnames(MB.Myeloid)]

## cluster based on GO AUC by NMF to find the metabolism similarity among seurat clusters
AUC.seurat =  CreateSeuratObject(counts = Myeloid_GO_AUC, min.cells = 3,meta.data = MB.Myeloid@meta.data)
AUC.seurat@active.ident = as.factor(MB.Myeloid$group)

# find cluster AUC similarity by average spearman coef
spearman = cor(Myeloid_GO_AUC,method = "spearman")
spearman = as.data.frame(spearman)

spearman_mean = data.frame()
for (i in 1:length(unique(MB.Myeloid$subtype))) {
  for (j in 1:length(unique(MB.Myeloid$subtype))) {
    a = spearman[colnames(MB.Myeloid)[MB.Myeloid$subtype == unique(MB.Myeloid$subtype)[i]],colnames(MB.Myeloid)[MB.Myeloid$subtype == unique(MB.Myeloid$subtype)[j]]]
    spearman_mean[i,j] = mean(as.matrix(a))
  }
}
rownames(spearman_mean) = unique(MB.Myeloid$subtype)
colnames(spearman_mean) = rownames(spearman_mean)

## visualized by corrplot
corrplot::corrplot(cor=as.matrix(spearman_mean),order="hclust",
                   method = "square",
                   type = "full",
                   tl.cex=0.8, pch=F, tl.srt = 45,
                   hclust.method = "complete",
                   insig = "label_sig",
                   sig.level=c(.001, .01, .05), pech.cex=0.003, pch.col="black",
                   cl.lim = c(-0.3, 0.3), is.corr = F,
                   col=colorRampPalette(c("deepskyblue3", "white", "tomato3"))(10),
                   tl.col="black")
# ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/10_Myeloid_subtype_GO_AUC_spearman_mean_heatmap.pdf",width = 7,height = 7)

save(spearman_mean,file = "/ix/gkohanbash/SimonsLab/R_Share/processing/10_Myeloid_subtype_GO_AUC_spearman_mean.RData")

# find seurat cluster marker patheay
AUC.seurat@active.ident = as.factor(MB.Myeloid$subtype)
AUC.markers = FindAllMarkers(AUC.seurat,
                             only.pos = F,
                             test.use = "wilcox",
                             logfc.threshold = 0)

save(AUC.markers,file = "/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/10_Myeloid_subtype_GO_AUC_marker.RData")

AUC.markers = AUC.markers[AUC.markers$p_val_adj < 0.05,]

# GOBP-POSITIVE-REGULATION-OF-TRANSCRIPTION-FROM-RNA-POLYMERASE-II-PROMOTER-IN-RESPONSE-TO-ENDOPLASMIC-RETICULUM-STRESS in MG_1
geneSetName <- rownames(cells_AUC)[grep("GOBP_POSITIVE_REGULATION_OF_TRANSCRIPTION_FROM_RNA_POLYMERASE_II_PROMOTER_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_POSITIVE_REGULATION_OF_TRANSCRIPTION_FROM_RNA_POLYMERASE_II_PROMOTER_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# GOBP-RESPONSE-TO-GRANULOCYTE-MACROPHAGE-COLONY-STIMULATING-FACTOR
geneSetName <- rownames(cells_AUC)[grep("GOBP_RESPONSE_TO_GRANULOCYTE_MACROPHAGE_COLONY_STIMULATING_FACTOR" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_RESPONSE_TO_GRANULOCYTE_MACROPHAGE_COLONY_STIMULATING_FACTOR_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# GOBP-INTEGRATED-STRESS-RESPONSE-SIGNALING
geneSetName <- rownames(cells_AUC)[grep("GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# GOBP-POSITIVE-REGULATION-OF-TRANSCRIPTION-FROM-RNA-POLYMERASE-II-PROMOTER-IN-RESPONSE-TO-STRESS
geneSetName <- rownames(cells_AUC)[grep("GOBP_POSITIVE_REGULATION_OF_TRANSCRIPTION_FROM_RNA_POLYMERASE_II_PROMOTER_IN_RESPONSE_TO_STRESS" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_POSITIVE_REGULATION_OF_TRANSCRIPTION_FROM_RNA_POLYMERASE_II_PROMOTER_IN_RESPONSE_TO_STRESS_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# GOBP-POSITIVE-REGULATION-OF-INTRINSIC-APOPTOTIC-SIGNALING-PATHWAY-BY-P53-CLASS-MEDIATOR
geneSetName <- rownames(cells_AUC)[grep("GOBP_POSITIVE_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_BY_P53_CLASS_MEDIATOR" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_POSITIVE_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_BY_P53_CLASS_MEDIATOR_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# GOBP-ESTABLISHMENT-OF-PROTEIN-LOCALIZATION-TO-ENDOPLASMIC-RETICULUM
geneSetName <- rownames(cells_AUC)[grep("GOBP_ESTABLISHMENT_OF_PROTEIN_LOCALIZATION_TO_ENDOPLASMIC_RETICULUM" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_ESTABLISHMENT_OF_PROTEIN_LOCALIZATION_TO_ENDOPLASMIC_RETICULUM_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# GOBP-POSITIVE-REGULATION-OF-INTRINSIC-APOPTOTIC-SIGNALING-PATHWAY-IN-RESPONSE-TO-DNA-DAMAGE
geneSetName <- rownames(cells_AUC)[grep("GOBP_POSITIVE_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_DNA_DAMAGE" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_POSITIVE_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_DNA_DAMAGE_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# GOBP-POSITIVE-REGULATION-OF-METALLOENDOPEPTIDASE-ACTIVITY
geneSetName <- rownames(cells_AUC)[grep("GOBP_POSITIVE_REGULATION_OF_METALLOENDOPEPTIDASE_ACTIVITY" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_POSITIVE_REGULATION_OF_METALLOENDOPEPTIDASE_ACTIVITY_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# GOBP-L-ASPARTATE-TRANSMEMBRANE-TRANSPORT
geneSetName <- rownames(cells_AUC)[grep("GOBP_L_ASPARTATE_TRANSMEMBRANE_TRANSPORT" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_L_ASPARTATE_TRANSMEMBRANE_TRANSPORT_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# GOBP-DETOXIFICATION-OF-INORGANIC-COMPOUND
geneSetName <- rownames(cells_AUC)[grep("GOBP_DETOXIFICATION_OF_INORGANIC_COMPOUND" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_DETOXIFICATION_OF_INORGANIC_COMPOUND_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# GOBP-MITOCHONDRIAL-ELECTRON-TRANSPORT-CYTOCHROME-C-TO-OXYGEN
geneSetName <- rownames(cells_AUC)[grep("GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_CYTOCHROME_C_TO_OXYGEN" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_CYTOCHROME_C_TO_OXYGEN_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# GOBP-MACROPHAGE-FUSION
geneSetName <- rownames(cells_AUC)[grep("GOBP_MACROPHAGE_FUSION" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_MACROPHAGE_FUSION_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# GOBP-MACROPHAGE-INFLAMMATORY-PROTEIN-1-ALPHA-PRODUCTION
geneSetName <- rownames(cells_AUC)[grep("GOBP_MACROPHAGE_INFLAMMATORY_PROTEIN_1_ALPHA_PRODUCTION" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_MACROPHAGE_INFLAMMATORY_PROTEIN_1_ALPHA_PRODUCTION_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# GOBP-REGULATION-OF-ANTIGEN-PROCESSING-AND-PRESENTATION-OF-PEPTIDE-OR-POLYSACCHARIDE-ANTIGEN-VIA-MHC-CLASS-II
geneSetName <- rownames(cells_AUC)[grep("GOBP_REGULATION_OF_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_OR_POLYSACCHARIDE_ANTIGEN_VIA_MHC_CLASS_II" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_REGULATION_OF_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_OR_POLYSACCHARIDE_ANTIGEN_VIA_MHC_CLASS_II_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# GOBP-POSITIVE-REGULATION-OF-MICROGLIAL-CELL-ACTIVATION
geneSetName <- rownames(cells_AUC)[grep("GOBP_POSITIVE_REGULATION_OF_MICROGLIAL_CELL_ACTIVATION" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_POSITIVE_REGULATION_OF_MICROGLIAL_CELL_ACTIVATION_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)

dev.off()

# GOBP-MICROGLIAL-CELL-MIGRATION
geneSetName <- rownames(cells_AUC)[grep("GOBP_MICROGLIAL_CELL_MIGRATION" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_MICROGLIAL_CELL_MIGRATION_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)

dev.off()

# GOBP-POSITIVE-REGULATION-OF-NEUROINFLAMMATORY-RESPONSE
geneSetName <- rownames(cells_AUC)[grep("GOBP_POSITIVE_REGULATION_OF_NEUROINFLAMMATORY_RESPONSE" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_POSITIVE_REGULATION_OF_NEUROINFLAMMATORY_RESPONSE_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)

dev.off()

# GOBP-POSITIVE-REGULATION-OF-CALCIUM-MEDIATED-SIGNALING
geneSetName <- rownames(cells_AUC)[grep("GOBP_POSITIVE_REGULATION_OF_CALCIUM_MEDIATED_SIGNALING" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_POSITIVE_REGULATION_OF_CALCIUM_MEDIATED_SIGNALING_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)

dev.off()

# GOBP-POSITIVE-REGULATION-OF-AMYLOID-BETA-CLEARANCE
geneSetName <- rownames(cells_AUC)[grep("GOBP_POSITIVE_REGULATION_OF_AMYLOID_BETA_CLEARANCE" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_POSITIVE_REGULATION_OF_AMYLOID_BETA_CLEARANCE_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)

dev.off()

# GOBP-ANTIGEN-PROCESSING-AND-PRESENTATION-OF-PEPTIDE-ANTIGEN-VIA-MHC-CLASS-IB
geneSetName <- rownames(cells_AUC)[grep("GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_IB" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_IB_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)

dev.off()

# GOBP-POSITIVE-REGULATION-OF-COMPLEMENT-ACTIVATION
geneSetName <- rownames(cells_AUC)[grep("GOBP_POSITIVE_REGULATION_OF_COMPLEMENT_ACTIVATION" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_POSITIVE_REGULATION_OF_COMPLEMENT_ACTIVATION_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)

dev.off()

# GOBP-MONOCYTE-AGGREGATION
geneSetName <- rownames(cells_AUC)[grep("GOBP_MONOCYTE_AGGREGATION" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_MONOCYTE_AGGREGATION_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# GOBP-POSITIVE-REGULATION-OF-CELL-ADHESION-MOLECULE-PRODUCTION
geneSetName <- rownames(cells_AUC)[grep("GOBP_POSITIVE_REGULATION_OF_CELL_ADHESION_MOLECULE_PRODUCTION" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_POSITIVE_REGULATION_OF_CELL_ADHESION_MOLECULE_PRODUCTION_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# GOBP-POSITIVE-REGULATION-OF-FEVER-GENERATION
geneSetName <- rownames(cells_AUC)[grep("GOBP_POSITIVE_REGULATION_OF_FEVER_GENERATION" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_POSITIVE_REGULATION_OF_FEVER_GENERATION_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# GOBP-NEGATIVE-REGULATION-OF-LEUKOCYTE-ADHESION-TO-VASCULAR-ENDOTHELIAL-CELL
geneSetName <- rownames(cells_AUC)[grep("GOBP_NEGATIVE_REGULATION_OF_LEUKOCYTE_ADHESION_TO_VASCULAR_ENDOTHELIAL_CELL" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_NEGATIVE_REGULATION_OF_LEUKOCYTE_ADHESION_TO_VASCULAR_ENDOTHELIAL_CELL_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# GOBP-LEUKOCYTE-CHEMOTAXIS-INVOLVED-IN-INFLAMMATORY-RESPONSE
geneSetName <- rownames(cells_AUC)[grep("GOBP_LEUKOCYTE_CHEMOTAXIS_INVOLVED_IN_INFLAMMATORY_RESPONSE" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_LEUKOCYTE_CHEMOTAXIS_INVOLVED_IN_INFLAMMATORY_RESPONSE_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# GOBP-MEIOTIC-CHROMOSOME-CONDENSATION
geneSetName <- rownames(cells_AUC)[grep("GOBP_MEIOTIC_CHROMOSOME_CONDENSATION" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_MEIOTIC_CHROMOSOME_CONDENSATION_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# GOBP-SPINDLE-ASSEMBLY-INVOLVED-IN-MEIOSIS
geneSetName <- rownames(cells_AUC)[grep("GOBP_SPINDLE_ASSEMBLY_INVOLVED_IN_MEIOSIS" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_SPINDLE_ASSEMBLY_INVOLVED_IN_MEIOSIS_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# GOBP-POSITIVE-REGULATION-OF-PROSTAGLANDIN-BIOSYNTHETIC-PROCESS
geneSetName <- rownames(cells_AUC)[grep("GOBP_POSITIVE_REGULATION_OF_PROSTAGLANDIN_BIOSYNTHETIC_PROCESS" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_POSITIVE_REGULATION_OF_PROSTAGLANDIN_BIOSYNTHETIC_PROCESS_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# GOBP-POSITIVE-REGULATION-OF-UNSATURATED-FATTY-ACID-BIOSYNTHETIC-PROCESS
geneSetName <- rownames(cells_AUC)[grep("GOBP_POSITIVE_REGULATION_OF_UNSATURATED_FATTY_ACID_BIOSYNTHETIC_PROCESS" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_POSITIVE_REGULATION_OF_UNSATURATED_FATTY_ACID_BIOSYNTHETIC_PROCESS_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# GOBP-HYPOXIA-INDUCIBLE-FACTOR-1ALPHA-SIGNALING-PATHWAY
geneSetName <- rownames(cells_AUC)[grep("GOBP_HYPOXIA_INDUCIBLE_FACTOR_1ALPHA_SIGNALING_PATHWAY" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_HYPOXIA_INDUCIBLE_FACTOR_1ALPHA_SIGNALING_PATHWAY_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# GOBP-POSITIVE-REGULATION-OF-MACROPHAGE-CYTOKINE-PRODUCTION
geneSetName <- rownames(cells_AUC)[grep("GOBP_POSITIVE_REGULATION_OF_MACROPHAGE_CYTOKINE_PRODUCTION" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_POSITIVE_REGULATION_OF_MACROPHAGE_CYTOKINE_PRODUCTION_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# GOBP-POSITIVE-REGULATION-OF-MONOCYTE-DIFFERENTIATION
geneSetName <- rownames(cells_AUC)[grep("GOBP_POSITIVE_REGULATION_OF_MONOCYTE_DIFFERENTIATION" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_GOBP_POSITIVE_REGULATION_OF_MONOCYTE_DIFFERENTIATION_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# AUC for metabolism pathway
# build KEGG ref geneset -- gene name, not entrenz 
load("/ix/gkohanbash/SimonsLab/R_Share/processing/5_KEGG_AUC_genelist.RData")

#ranking cells by gene expression high -> low
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=nCores, plotStats=F)

#Calculate enrichment for the gene signatures (AUC)
cells_AUC <- AUCell_calcAUC(h.sets, cells_rankings,
                            aucMaxRank=nrow(cells_rankings)*0.05)#The percentage to take into account can be modified with the argument aucMaxRank
#by default only the top 5% of the genes in the ranking are used
#the AUC estimates the proportion of genes in the gene-set that are highly expressed in each cell.
#Cells expressing many genes from the gene-set will have higher AUC values than cells expressing fewer

#Determine the cells with the given gene signatures or active gene sets - find threshold
set.seed(123)
save(cells_AUC, file = "/ix/gkohanbash/SimonsLab/R_Share/processing/10_Myeloid_KEGG_AUC.RData")

## select significant activity-different GO item between injury and sham group (containing SCFA feeding group)
Myeloid_KEGG_AUC = as.data.frame(cells_AUC@assays@data@listData[["AUC"]])
Myeloid_KEGG_AUC = Myeloid_KEGG_AUC[,colnames(MB.Myeloid)]

## cluster based on KEGG AUC by NMF to find the metabolism similarity among seurat clusters
AUC.seurat =  CreateSeuratObject(counts = Myeloid_KEGG_AUC, min.cells = 3,meta.data = MB.Myeloid@meta.data)
AUC.seurat@active.ident = as.factor(MB.Myeloid$group)

# find cluster AUC similarity by average spearman coef
spearman = cor(Myeloid_KEGG_AUC,method = "spearman")
spearman = as.data.frame(spearman)

spearman_mean = data.frame()
for (i in 1:length(unique(MB.Myeloid$subtype))) {
  for (j in 1:length(unique(MB.Myeloid$subtype))) {
    a = spearman[colnames(MB.Myeloid)[MB.Myeloid$subtype == unique(MB.Myeloid$subtype)[i]],colnames(MB.Myeloid)[MB.Myeloid$subtype == unique(MB.Myeloid$subtype)[j]]]
    spearman_mean[i,j] = mean(as.matrix(a))
  }
}
rownames(spearman_mean) = unique(MB.Myeloid$subtype)
colnames(spearman_mean) = rownames(spearman_mean)

## visualized by corrplot
corrplot::corrplot(cor=as.matrix(spearman_mean),order="hclust",
                   method = "square",
                   type = "full",
                   tl.cex=0.8, pch=F, tl.srt = 45,
                   hclust.method = "complete",
                   insig = "label_sig",
                   sig.level=c(.001, .01, .05), pech.cex=0.003, pch.col="black",
                   cl.lim = c(-0.3, 0.3), is.corr = F,
                   col=colorRampPalette(c("deepskyblue3", "white", "tomato3"))(10),
                   tl.col="black")
# ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/10_Myeloid_subtype_KEGG_AUC_spearman_mean_heatmap.pdf",width = 7,height = 7)

save(spearman_mean,file = "/ix/gkohanbash/SimonsLab/R_Share/processing/10_Myeloid_subtype_KEGG_AUC_spearman_mean.RData")

# find seurat cluster marker patheay
AUC.seurat@active.ident = as.factor(MB.Myeloid$subtype)
AUC.markers = FindAllMarkers(AUC.seurat,
                             only.pos = F,
                             test.use = "wilcox",
                             logfc.threshold = 0)

save(AUC.markers,file = "/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/10_Myeloid_subtype_KEGG_AUC_marker.RData")

AUC.markers = AUC.markers[AUC.markers$p_val_adj < 0.05,]

# RAY_ALZHEIMERS_DISEASE
geneSetName <- rownames(cells_AUC)[grep("RAY_ALZHEIMERS_DISEASE" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_RAY_ALZHEIMERS_DISEASE_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# KEGG-RIBOSOME
geneSetName <- rownames(cells_AUC)[grep("KEGG_RIBOSOME" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_KEGG_RIBOSOME_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# REACTOME-EUKARYOTIC-TRANSLATION-ELONGATION
geneSetName <- rownames(cells_AUC)[grep("REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# REACTOME-REGULATION-OF-EXPRESSION-OF-SLITS-AND-ROBOS
geneSetName <- rownames(cells_AUC)[grep("REACTOME_REGULATION_OF_EXPRESSION_OF_SLITS_AND_ROBOS" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_REACTOME_REGULATION_OF_EXPRESSION_OF_SLITS_AND_ROBOS_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# BIOCARTA-FOSB-PATHWAY
geneSetName <- rownames(cells_AUC)[grep("BIOCARTA_FOSB_PATHWAY" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_BIOCARTA_FOSB_PATHWAY_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# WP-NEUROINFLAMMATION
geneSetName <- rownames(cells_AUC)[grep("WP_NEUROINFLAMMATION" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_WP_NEUROINFLAMMATION_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# REACTOME-METALLOTHIONEINS-BIND-METALS
geneSetName <- rownames(cells_AUC)[grep("REACTOME_METALLOTHIONEINS_BIND_METALS" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_REACTOME_METALLOTHIONEINS_BIND_METALS_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# BIOCARTA-CCR5-PATHWAY
geneSetName <- rownames(cells_AUC)[grep("BIOCARTA_CCR5_PATHWAY" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_BIOCARTA_CCR5_PATHWAY_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# RAY-ALZHEIMERS-DISEASE
geneSetName <- rownames(cells_AUC)[grep("RAY_ALZHEIMERS_DISEASE" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_RAY_ALZHEIMERS_DISEASE_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# REACTOME-INTERLEUKIN-10-SIGNALING
geneSetName <- rownames(cells_AUC)[grep("REACTOME_INTERLEUKIN_10_SIGNALING" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_REACTOME_INTERLEUKIN_10_SIGNALING_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# PHONG-TNF-TARGETS-UP
geneSetName <- rownames(cells_AUC)[grep("PHONG_TNF_TARGETS_UP" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_PHONG_TNF_TARGETS_UP_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# REACTOME-CHEMOKINE-RECEPTORS-BIND-CHEMOKINES
geneSetName <- rownames(cells_AUC)[grep("REACTOME_CHEMOKINE_RECEPTORS_BIND_CHEMOKINES" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_REACTOME_CHEMOKINE_RECEPTORS_BIND_CHEMOKINES_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# WP-MACROPHAGE-MARKERS
geneSetName <- rownames(cells_AUC)[grep("WP_MACROPHAGE_MARKERS" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_WP_MACROPHAGE_MARKERS_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# REACTOME-CLEC7A-INFLAMMASOME-PATHWAY
geneSetName <- rownames(cells_AUC)[grep("REACTOME_CLEC7A_INFLAMMASOME_PATHWAY" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_REACTOME_CLEC7A_INFLAMMASOME_PATHWAY_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# REACTOME-INTERLEUKIN-1-PROCESSING
geneSetName <- rownames(cells_AUC)[grep("REACTOME_INTERLEUKIN_1_PROCESSING" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_REACTOME_INTERLEUKIN_1_PROCESSING_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# REACTOME-CONDENSATION-OF-PROMETAPHASE-CHROMOSOMES
geneSetName <- rownames(cells_AUC)[grep("REACTOME_CONDENSATION_OF_PROMETAPHASE_CHROMOSOMES" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_REACTOME_CONDENSATION_OF_PROMETAPHASE_CHROMOSOMES_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# EGUCHI-CELL-CYCLE-RB1-TARGETS
geneSetName <- rownames(cells_AUC)[grep("EGUCHI_CELL_CYCLE_RB1_TARGETS" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_EGUCHI_CELL_CYCLE_RB1_TARGETS_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()


# REACTOME-FIBRONECTIN-MATRIX-FORMATION
geneSetName <- rownames(cells_AUC)[grep("REACTOME_FIBRONECTIN_MATRIX_FORMATION" , rownames(cells_AUC))]

# plot pathway activity by signature gene mean value
nBreaks = 10
logMat <- log2(exprMatrix+2)

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))

gset = geneSetName[1]

cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(as.numeric(meanByGs[gset,]), breaks=nBreaks, right=FALSE,include.lowest=TRUE))], names(meanByGs[gset,]))

pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cluster_REACTOME_FIBRONECTIN_MATRIX_FORMATION_UMAP.pdf",width = 5,height = 5)
plot(cellsUMAP, main=gset, axes=FALSE, xlab="", ylab="",
     sub="Pathway Activity",
     col=cellColor[rownames(cellsUMAP)], pch=16)
dev.off()

# pseudotime - cytoTRACE
input = as.matrix(MB.Myeloid@assays$SCT@counts)
results <- CytoTRACE(input,ncores = 16)

save(results,file = "/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/10_Myeloid_cytoTRACE.RData")
## visualization
dir.create(paste0("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/cytoTRACE/Myeloid"))
dir.create(paste0("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/cytoTRACE/Myeloid/subtype"))

pheno = as.character(MB.Myeloid$group)
names(pheno) = colnames(MB.Myeloid)
plotCytoTRACE(results,
              phenotype = pheno, 
              emb = MB.Myeloid@reductions[["umap"]]@cell.embeddings,
              otherName = "group",
              # enableFast = FALSE,
              outputDir = paste0("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/cytoTRACE/Myeloid/"))

pheno = as.character(MB.Myeloid$subtype)
names(pheno) = colnames(MB.Myeloid)
plotCytoTRACE(results,
              phenotype = pheno, 
              emb = MB.Myeloid@reductions[["umap"]]@cell.embeddings,
              otherName = "subtype",
              # enableFast = FALSE,
              outputDir = paste0("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/cytoTRACE/Myeloid/subtype/"))

# boxplot of cytoTRACEpseudotime in different subtypes and groups
# calculate p value
input = data.frame(Pseudotime = results$CytoTRACE,
                   group = MB.Myeloid$group,
                   celltype = MB.Myeloid$celltype)
dataset = input

group=levels(factor(dataset$group))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

dodge <- position_dodge(width = 1)

dataset%>%ggplot(aes(group,Pseudotime))+geom_violin(aes(fill=group),scale = "width",position = dodge)+
  # facet_grid(dataset$condition~.,scales = "free_y")+
  geom_boxplot(aes(fill=group),position = dodge,outlier.size = 0,notch = T,width = 0.1)+
  scale_fill_brewer(palette = "Set3",direction = 1)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1,size = 13),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "right",
    axis.text.y.left = element_text(size = 15),
    legend.text = element_text(size = 13)
  )+
  stat_compare_means(comparisons = my_comparisons,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),label = "p.signif")+
  ggtitle("cytoTRACE Pseudotime Myeloid")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_cytoTRACE_pseudotime_distribution.pdf",width = 10,height = 5)

# scVelo
## load scvelo module
library(reticulate)
library(SeuratDisk)
reticulate::use_python("/ihome/gkohanbash/zxiong/.local/share/r-miniconda/envs/r-reticulate/bin/python",required = T)

mat = import("matplotlib")
mat$use("Agg")

scv <- import("scvelo")
scv$set_figure_params('scvelo')  # for beautified visualization

# input generation
load("/ix/gkohanbash/SimonsLab/R_Share/processing/merged_splice.RData")
Myeloid = merged[,colnames(MB.Myeloid)]
Myeloid$subtype = as.character(MB.Myeloid$subtype)
Myeloid[["umap"]] = MB.Myeloid@reductions[["umap"]]
Myeloid@reductions[["umap"]]@assay.used = "RNA"

SaveH5Seurat(Myeloid, filename = "/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/10_Myeloid_velocity.h5Seurat")
Convert("/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/10_Myeloid_velocity.h5Seurat", dest = "h5ad")

# load input data
adata <- scv$read("/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/10_Myeloid_velocity.h5ad")
adata

# scv$pl$proportions(adata,save = "/ix/gkohanbash/SimonsLab/R_Share/fig/10_Myeloid_velocyto_pctdemo.svg")

# normalized the input -- log(x + 1)
scv$pp$filter_and_normalize(adata, min_shared_counts = as.integer(20), n_top_genes = as.integer(2000))
scv$pp$moments(adata, n_pcs = as.integer(20),n_neighbors = as.integer(30))
scv$tl$recover_dynamics(adata) ## model
scv$tl$velocity(adata,
                mode='dynamical'
                ,diff_kinetics = T                                             
                # to consider the situation that different cluster's transcription velocity stable state is different
) 
scv$tl$velocity_graph(adata)
scv$pl$velocity_embedding_stream(adata, basis = "umap",
                                 #save = "/ix/gkohanbash/SimonsLab/R_Share/fig/10_Myeloid_velocyto_demo.svg", 
                                 color = "subtype")
# save 10_Myeloid_velocyto_differentiation 6*8

# condiments - show cluster 4,5,6 are at different lineages by slingshot
dir.create(paste0("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/condiments/Myeloid"))

set.seed(2071)
theme_set(theme_classic())

MB.Myeloid@active.ident = as.factor(MB.Myeloid$subtype)
input = MB.Myeloid[,-grep("MAC",MB.Myeloid$subtype)]
input = input[,input$subtype != "MG_10"]
sce = as.SingleCellExperiment(input,assay = "SCT")

# input 
df = as.data.frame(input@reductions[["umap"]]@cell.embeddings)
df$condition = input$group
df$cluster = input$subtype

# imbalance score to see if condition distribution merge evently 
scores <- imbalance_score(Object = df %>% 
                            dplyr::select(UMAP_1, UMAP_2) %>% 
                            as.matrix(),
                          conditions = df$condition)
df$scores <- scores$scores
df$scaled_scores <- scores$scaled_scores 

ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = scaled_scores)) +
  geom_point(alpha = .3) +
  scale_color_viridis_c(option = "C")+
  labs(col = "Score")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/condiments/Myeloid/imbalance_score.pdf")

# save imbalance_score.pdf
# the higher score represent the part of cells more come from one single condition

# Step 1 -  topology test -- is the trajactory of all condition can fit one common trajactory? H0 = can fit one common trajactory
sce <- slingshot(sce, reducedDim = "UMAP",
                 clusterLabels = sce$subtype,
                 # start.clus = "MG_2",
                 reweight = F, reassign = F)
SlingshotDataSet(sce)
# lineages: 2 
# Lineage1: MG_4  MG_8  MG_2  MG_1  MG_7  MG_6  MG_9  MG_3  
# Lineage2: MG_4  MG_8  MG_2  MG_1  MG_7  MG_5  

sce$condition =  as.character(sce$group)

sdss <- slingshot_conditions(SlingshotDataSet(sce), sce$condition, approx_points = FALSE,
                             extend = "n", reweight = F, reassign = F)

curves <- bind_rows(lapply(sdss, slingCurves, as.df = TRUE),
                    .id = "condition")

ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = condition)) +
  geom_point(size = .7, alpha = .3) +
  scale_color_brewer(palette = "Accent") +
  labs(col = "Lineage") +
  geom_path(data = curves %>% arrange(Lineage, Order),
            aes(group = interaction(Lineage, condition)), size = 1.5) +
  theme(legend.position = c(.15, .35),
        legend.background = element_blank())

ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/condiments/Myeloid/slingshot_lineage_UMAP.pdf")
save(sce,file = "/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/10_Myeloid_slingshot.RData")

## plot
df <-  slingPseudotime(sce) %>% as.data.frame() 
df$condition <- sce$condition
df <- df %>% 
  pivot_longer(-condition, names_to = "Lineage",
               values_to = "pst") %>%
  filter(!is.na(pst))

ggplot(df, aes(x = pst)) +
  geom_density(alpha = .8, aes(fill = condition), col = "transparent") +
  geom_density(aes(col = condition), fill = "transparent", size = 1.5) +
  labs(x = "Pseudotime", fill = "condition") +
  facet_wrap(~Lineage, scales = "free_x") +
  guides(col = "none", fill = guide_legend(
    override.aes = list(size = 1.5)
  )) +
  scale_fill_brewer(palette = "Accent") +
  scale_color_brewer(palette = "Accent") +
  theme(legend.position = c(.8, .8), legend.direction = "horizontal")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/condiments/Myeloid/slingshot_lineage_progression_UMAP.pdf")

# save 	lineage_progression_diff.pdf

# lineage 1-4 is the most different lineage between control diet and SCFA diet
input$lineage_1 = sce$slingPseudotime_1

dataset = data.frame(Pseudotime = input$lineage_1,
                     group = input$group)
dataset = na.omit(dataset)

group=levels(factor(dataset$group))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

dodge <- position_dodge(width = 1)

dataset%>%ggplot(aes(group,Pseudotime))+geom_violin(aes(fill=group),scale = "width",position = dodge)+
  # facet_grid(dataset$condition~.,scales = "free_y")+
  geom_boxplot(aes(fill=group),position = dodge,outlier.size = 0,notch = T,width = 0.1)+
  # scale_fill_brewer(palette = "Set3",direction = 1)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1,size = 13),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "right",
    axis.text.y.left = element_text(size = 15),
    legend.text = element_text(size = 13)
  )+
  stat_compare_means(comparisons = my_comparisons,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),label = "p.signif")+
  ggtitle("Slingshot Pseudotime - lineage 1")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_Slingshot_Pseudotime_lineage_1_group_violin.pdf",width = 10,height = 6)

## due to lineage 2 is similar between injury_control and injury_SCFA
# Lamian XDE for control/SCFA diet in lineage 2
options(warn=-1)
suppressMessages(library(Lamian))

input = MB.Myeloid[,MB.Myeloid$subtype != "MAC"]
input = input[,input$subtype != "MG_10"]

input@active.ident = as.factor(input$group)
input$lineage_1 = sce$slingPseudotime_2
input = subset(input,idents  = c("injury_SCFA","injury_control"))
input$diet = if_else(input$group == "injury_control","control","SCFA")

DefaultAssay(input) = "SCT"
# input = NormalizeData(input,normalization.method = "LogNormalize",scale.factor = 10000)
# input <- ScaleData(input, features = rownames(input))

# input <- FindVariableFeatures(object = input,
#                               mean.function = ExpMean,
#                               dispersion.function = LogVMR,
#                               x.low.cutoff = 0.0125,
#                               x.high.cutoff = 3,
#                               y.cutoff = 0.5,
#                               nfeatures = 3000)
# DefaultAssay(input) = "RNA"

input = input[,!is.na(input$lineage_1)]
# for specific lineage
# XDE test input for Lamian - data preparation
set.seed(12345)
dataset = as.data.frame(GetAssayData(input,slot = "data"))#[VariableFeatures(input),]
# Values are library-size-normalized log-transformed gene expression matrix 
colnames(dataset) = gsub("-","_",colnames(dataset))

# pseudotime = read.table("/bgfs/gkohanbash/TBI_Shared_files/cytoTRACE/Microglia/condition/umap/CytoTRACE_plot_table.txt")
pseudotime = na.omit(input$lineage_1)
pseudotime = rank(pseudotime)    # pseudotime input should be the pseudotime rank: 1, 2, 3, 4, etc.
pseudotime[which.max(pseudotime)] = ceiling(max(pseudotime))
names(pseudotime) = colnames(dataset)

cellanno = data.frame(Cell = colnames(dataset), 
                      Sample = input$sample) # different samples/patients which cells from

design = data.frame(intercept = rep(1,length(unique(input$sample))),
                    group = c(1,1,1,1,0,0,0,0))
# SCFA = 1; control = 0
rownames(design) = unique(input$sample)

Res <- lamian.test(expr = as.matrix(dataset), # matrix: Values are library-size-normalized log-transformed gene expression matrix 
                   cellanno = cellanno, 
                   pseudotime = pseudotime, 
                   design = design, 
                   test.type = 'variable', 
                   permuiter = 100)

## get differential dynamic genes statistics 
stat <- Res$statistics
head(stat)

stat <- stat[order(stat[,1], -stat[,3]), ]
## identify XDE genes with FDR.overall < 0.05 cutoff
diffgene <- rownames(stat[stat[, grep('^fdr.*overall$', colnames(stat))] < 0.05,])
str(diffgene)

## population fit
Res$populationFit <- getPopulationFit(Res, gene = diffgene, type = 'variable')

## clustering
Res$covariateGroupDiff <- getCovariateGroupDiff(testobj = Res, gene = diffgene)

DDGType <- getDDGType(Res)
cluster_gene <- function(testobj, 
                         gene,
                         k,
                         k.auto = FALSE,
                         type = 'Time',
                         method = 'kmeans', 
                         scale.difference = F){
  if (toupper(type) == 'TIME'){ 
    if ('populationFit' %in% names(testobj)) {
      fit <- testobj$populationFit
    } else {
      fit <- getPopulationFit(testobj, gene = gene, type = type)
    }
  } else if (toupper(type) == 'VARIABLE'){
    if ('covariateGroupDiff' %in% names(testobj)){
      fit <- testobj$covariateGroupDiff
    } else{
      fit <- getCovariateGroupDiff(testobj = testobj, gene = gene)  
    }
  }
  if (scale.difference){
    mat.scale <- scalematrix(fit[gene, ,drop=F])
  } else {
    max <- apply(abs(fit[gene, ,drop=F]), 1, max)
    mat.scale <- fit[gene, ,drop=F]/max
  }
  
  if (method == 'kmeans'){
    set.seed(12345)
    # 
    if (k.auto){
      clu <- mykmeans(mat.scale, maxclunum = 20)$cluster
    } else {
      clu <- kmeans(mat.scale, k, iter.max = 1000)$cluster
    }
  } else if (method == 'hierarchical') {
    clu <- cutree(hclust(dist(mat.scale)), k = k)
  } else if (method == 'louvain'){
    graph = scran::buildSNNGraph(mat.scale, transposed=T,k=k,d=NA)
    res = igraph::cluster_louvain(graph)$membership
    if (max(res) <= k){
      hclu <- hclust(dist(mat.scale))
      clu <- cutree(hclu,k)
    } else {
      cc <- aggregate(mat.scale, list(res), mean)
      cc <- as.matrix(cc[,-1])
      hclu <- hclust(dist(cc))
      clu <- cutree(hclu,k)
      clu <- clu[res]      
    }
    names(clu) = row.names(mat.scale)
  } else if (method == 'GMM'){
    samplen = 2e2
    colnames(mat.scale) = paste0('cell', seq(1, ncol(mat.scale)))
    set.seed(12345)
    sampid = sample(1:ncol(mat.scale), samplen)
    if (nrow(mat.scale) > samplen){
      res <- mclust::Mclust(data = mat.scale[, sampid], G = k, modelNames = 'EII', verbose = FALSE)
    } else {
      res <- mclust::Mclust(data = mat.scale[, sampid], G = k, modelNames = 'VII', verbose = FALSE)
    }
    clu <- apply(res$z, 1, which.max)
  }
  
  # order clusters by genes' earliest max expr position
  v <- sapply(unique(clu), function(i){
    ap <- which(colMeans(mat.scale[names(clu)[clu==i], -ncol(mat.scale), drop=FALSE]) * colMeans(mat.scale[names(clu)[clu==i], -1, drop = FALSE]) < 0)
    if(length(ap) == 0){
      1
    } else{
      ap[which.min(abs(ap-ncol(mat.scale)/2))]  
    }
  })
  names(v) <- unique(clu)
  
  corv <- apply(mat.scale,1,cor,1:ncol(mat.scale))
  corv <- tapply(corv,list(clu),mean)
  corv <- corv[names(v)]
  # self study
  v[corv < 0] <- ncol(mat.scale)-v[corv < 0]
  if (toupper(type) == 'VARIABLE'){
    v <- v * (2*(corv > 0)-1)
  }
  
  trans <- cbind(as.numeric(names(sort(v))),1:length(v))
  n <- names(clu)
  clu <- trans[match(clu,trans[,1]),2]
  names(clu) <- n
  
  if (toupper(type) == 'VARIABLE'){
    clu2 <- paste0(clu, ';',rowMeans(fit[names(clu), , drop=F]) > 0)  
  } else {
    clu2 <- paste0(clu, ';TRUE')
  }
  uclu2 <- sort(unique(clu2))
  clu2 <- match(clu2,uclu2)
  names(clu2) <- n
  return(clu2)  
}
Res$cluster <- cluster_gene(Res, gene = diffgene, type = "variable", 
                            k = 3, scale.difference = F, 
                            method = "kmeans")

# Res$cluster <- clusterGene_edit(Res, gene = diffgene, type = 'variable', k.auto = T)
table(Res$cluster)

colnames(Res$populationFit[[1]]) <- colnames(Res$populationFit[[2]]) <- colnames(Res$expr) 
names(Res[["cluster"]]) = diffgene
# Res[["populationFit"]]=Res[["populationFit"]][c("group_1","group_2","group_3","group_4")]
# order the order of groups

save(Res,file = "/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/10_Myeloid_Lamian_RNA_lineage_2.RData")

library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- apply(data, 1, sd)
  (data - cm) / csd
}

plotXDEHm_edit = function (testobj, showRowName = FALSE, cellWidthTotal = 250, 
                           cellHeightTotal = 400, showCluster = FALSE, colann = NULL, 
                           rowann = NULL, annotation_colors = NULL, subsampleCell = TRUE, 
                           numSubsampleCell = 1000, sep = NA, break.0 = TRUE) 
{
  fit <- testobj$populationFit
  if ("DDGType" %in% names(testobj)) {
    DDGType <- testobj$DDGType
  }
  else {
    DDGType <- getDDGType(testobj)
  }
  DDGType <- DDGType[rownames(testobj$covariateGroupDiff)]
  if (subsampleCell) {
    id <- round(seq(1, ncol(fit[[1]]), length.out = numSubsampleCell))
    for (i in seq_len(length(fit))) {
      fit[[i]] <- fit[[i]][, id]
    }
    if (sum(DDGType == "meanSig", na.rm = TRUE) > 0) {
      meanid <- which(DDGType == "meanSig")
      max <- apply(abs(testobj$covariateGroupDiff[-meanid, 
                                                  id, drop = FALSE]), 1, max)
      FitDiff.scale1 <- testobj$covariateGroupDiff[-meanid, 
                                                   id, drop = FALSE]/max
      FitDiff.scale <- rbind(FitDiff.scale1, testobj$covariateGroupDiff[meanid, 
                                                                        id, drop = FALSE])
      FitDiff.scale <- FitDiff.scale[rownames(testobj$covariateGroupDiff), 
                                     , drop = FALSE]
    }
    else {
      max <- apply(abs(testobj$covariateGroupDiff[, id, 
                                                  drop = FALSE]), 1, max)
      FitDiff.scale <- testobj$covariateGroupDiff[, id, 
                                                  drop = FALSE]/max
    }
    FitDiff.sd <- scalematrix(testobj$covariateGroupDiff[, 
                                                         id, drop = FALSE])
    colnames(FitDiff.scale) <- paste0("FitDiff:cell", seq(1, 
                                                          ncol(FitDiff.scale)))
    testobj$pseudotime <- sort(sample(testobj$pseudotime, 
                                      numSubsampleCell))
    print("subsample done!")
  }
  else {
    if (sum(DDGType == "meanSig", na.rm = TRUE) > 0) {
      meanid <- which(DDGType == "meanSig")
      max <- apply(abs(testobj$covariateGroupDiff[-meanid, 
                                                  , drop = FALSE]), 1, max)
      FitDiff.scale1 <- testobj$covariateGroupDiff[-meanid, 
                                                   , drop = FALSE]/max
      FitDiff.scale <- rbind(FitDiff.scale1, testobj$covariateGroupDiff[meanid, 
                                                                        , drop = FALSE])
      FitDiff.scale <- FitDiff.scale[rownames(testobj$covariateGroupDiff), 
                                     , drop = FALSE]
    }
    else {
      max <- apply(abs(testobj$covariateGroupDiff), 1, 
                   max)
      FitDiff.scale <- testobj$covariateGroupDiff/max
    }
    colnames(FitDiff.scale) <- paste0("FitDiff:cell", seq(1, 
                                                          ncol(FitDiff.scale)))
    FitDiff.sd <- scalematrix(testobj$covariateGroupDiff)
  }
  max <- apply(abs(testobj$covariateGroupDiff), 1, max)
  alluniformdiff <- testobj$covariateGroupDiff/max
  oridata <- testobj$covariateGroupDiff
  fit.bak = fit
  clu <- testobj$cluster
  rownames(testobj$cellanno) <- testobj$cellanno[, 1]
  testobj$cellanno <- testobj$cellanno[names(testobj$pseudotime), 
  ]
  if ("expr.ori" %in% names(testobj)) {
    testobj$expr <- testobj$expr.ori[, names(testobj$pseudotime)]
  }
  else {
    testobj$expr <- testobj$expr[, names(testobj$pseudotime)]
  }
  if ("DDGType" %in% names(testobj)) {
    DDGType <- testobj$DDGType
  }
  else {
    DDGType <- getDDGType(testobj)
  }
  fit.scale <- do.call(cbind, fit)
  fit.scale <- fit.scale[names(testobj$cluster), ]
  fit.scale <- scalematrix(fit.scale)
  colnames(fit.scale) <- paste0(rep(names(fit), each = ncol(fit.scale)/length(fit)), 
                                ";cell", seq(1, ncol(fit.scale)))
  changepoint <- sapply(names(clu), function(i) {
    ap <- which(FitDiff.sd[i, -ncol(FitDiff.sd)] * FitDiff.sd[i, 
                                                              -1] < 0)
    ap[which.min(abs(ap - ncol(FitDiff.sd)/2))]
  })
  res <- data.frame(clu = clu, cor = sapply(names(clu), function(i) cor(FitDiff.scale[i, 
                                                                                      seq(1, ncol(FitDiff.scale))], seq(1, ncol(FitDiff.scale)))), 
                    changepoint = changepoint, DDGType = DDGType[names(clu)])
  res1 <- res[res$DDGType == "trendSig", ]
  res2 <- res[res$DDGType == "bothSig", ]
  res3 <- res[res$DDGType == "other", ]
  res4 <- res[res$DDGType == "meanSig", ]
  o1 <- rownames(res1)[order(res1$clu, res1$cor > 0, res1$changepoint)]
  pn <- rowMeans(alluniformdiff[rownames(res2), , drop = FALSE])
  o2 <- rownames(res2)[order(pn > 0, res2$clu, res2$cor > 
                               0, res2$changepoint)]
  o3 <- rownames(res3)[order(res3$clu, res3$cor > 0, res3$changepoint)]
  o4 <- rownames(res4)[order(res4$clu)]
  res <- res[c(o1, o2, o4), ]
  rle <- rle(paste0(res$clu, res$DDGType))$lengths
  clu[rownames(res)] <- rep(seq_len(length(rle)), rle)
  fit.scale <- fit.scale[rownames(res), ]
  FitDiff.scale <- FitDiff.scale[rownames(res), ]
  cellanno <- testobj$cellanno
  expr = testobj$expr
  expr <- expr[, names(testobj$pseudotime)]
  tmp <- lapply(names(fit), function(i) {
    expr[rownames(fit.scale), colnames(expr) %in% cellanno[cellanno[, 
                                                                    2] %in% rownames(testobj$design)[testobj$design[, 
                                                                                                                    2] == sub(".*_", "", i)], 1]]
  })
  expr.scale <- do.call(cbind, tmp)
  expr.scale <- scalematrix(expr.scale)
  expr.scale <- expr.scale[rownames(fit.scale), ]
  if (is.null(colann)) {
    colann <- data.frame(pseudotime = testobj$pseudotime[colnames(expr.scale)], 
                         group = as.character(testobj$design[cellanno[match(colnames(expr.scale), 
                                                                            cellanno[, 1]), 2], 2]), expression = "Original", 
                         stringsAsFactors = FALSE)
    col.group = colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(length(unique(colann$group)) + 
                                                                       1)
    names(col.group) = c("NA", unique(colann$group))
  }
  rownames(colann) = colnames(expr.scale)
  col.expression = brewer.pal(n = 8, name = "Pastel2")[seq_len(3)]
  names(col.expression) = c("Original", "ModelFitted", "ModeledGroupDiff")
  col.pseudotime = colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(length(unique(colann$pseudotime)))
  names(col.pseudotime) = unique(colann$pseudotime)
  if (is.null(rowann)) {
    rowann = data.frame(cluster = factor(as.character(clu), 
                                         levels = as.character(seq_len(max(clu)))), DDGType = factor(as.character(DDGType[names(clu)]), 
                                                                                                     levels = as.character(unique(res[, 4]))), stringsAsFactors = FALSE)
    rownames(rowann) = names(clu)
  }
  rowann <- rowann[rownames(fit.scale), , drop = FALSE]
  rowann[, "DDGType"] <- factor(as.character(rowann[, "DDGType"]), 
                                levels = c("trendSig", "meanSig", "bothSig", "nonDDG", 
                                           "other"))
  if (length(unique(clu)) < 8) {
    col.clu = brewer.pal(8, "Set1")[seq_len(length(unique(clu)))]
  }
  else {
    col.clu = colorRampPalette(brewer.pal(8, "Set1"))(length(unique(clu)))
  }
  col.clu <- sample(col.clu)
  names(col.clu) = levels(rowann$clu)
  if (is.null(colann) | is.null(annotation_colors)) {
    col.meanDiff = c("blue", "red")
    names(col.meanDiff) <- c("Positive", "Negative")
    col.DDGType = brewer.pal(8, "Set3")[seq_len(3)]
    names(col.DDGType) = c("trendSig", "bothSig", "meanSig")
    annotation_colors = list(pseudotime = col.pseudotime, 
                             group = col.group, expression = col.expression, 
                             cluster = col.clu, DDGType = col.DDGType)
  }
  col.gs <- c("pink", "skyblue")
  names(col.gs) <- c("No", "Yes")
  col.limmaPb <- c("pink", "skyblue")
  names(col.limmaPb) <- c("nonDiff", "Diff")
  annotation_colors[["gs"]] <- col.gs
  annotation_colors[["limmaPb"]] <- col.limmaPb
  col.signalType <- brewer.pal(8, "Set3")[seq_len(3)]
  names(col.signalType) <- c("trend only", "mean only", "both")
  annotation_colors[["signalType"]] <- col.signalType
  cpl = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  if (break.0) {
    cpl <- c(cpl[seq_len(40)], cpl[60:100])
  }
  plist <- list()
  if (!is.na(sep)) {
    rownames(expr.scale) <- sub(sep, "", rownames(expr.scale))
    rownames(rowann) <- sub(sep, ":.*", rownames(rowann))
    rownames(oridata) <- sub(sep, "", rownames(oridata))
  }
  p1data <- expr.scale
  p1data[p1data > quantile(as.vector(p1data), 0.95, na.rm = TRUE)] <- quantile(as.vector(p1data), 
                                                                               0.95, na.rm = TRUE)
  p1data[p1data < quantile(as.vector(p1data), 0.05, na.rm = TRUE)] <- quantile(as.vector(p1data), 
                                                                               0.05, na.rm = TRUE)
  col_fun = colorRamp2(seq(-max(abs(p1data)), max(abs(p1data)), 
                           length.out = 50), colorRampPalette(c("blue3", "skyblue", 
                                                                "white", "pink", "red3"))(50))
  pt_col_fun = colorRamp2(seq(1, max(colann$pseudotime)), 
                          colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(max(colann$pseudotime)))
  annotation_colors$pseudotime <- pt_col_fun
  ht1 <- Heatmap(p1data, cluster_rows = FALSE, cluster_columns = FALSE, 
                 show_row_names = showRowName, show_column_names = FALSE, 
                 col = col_fun, heatmap_legend_param = list(legend_direction = "horizontal"), 
                 top_annotation = HeatmapAnnotation(df = colann, col = annotation_colors, 
                                                    show_annotation_name = FALSE), left_annotation = rowAnnotation(df = rowann, 
                                                                                                                   col = annotation_colors), width = 2)
  colann.fit1 <- data.frame(pseudotime = rep(seq_len(ncol(fit[[1]])), 
                                             length(fit)), group = gsub(sub("_.*", "_", names(fit)[1]), 
                                                                        "", sub(";.*", "", colnames(fit.scale))), expression = "ModelFitted", 
                            stringsAsFactors = FALSE)
  colann.fit2 <- data.frame(pseudotime = seq(1, ncol(FitDiff.scale)), 
                            group = "NA", expression = "ModeledGroupDiff", stringsAsFactors = FALSE)
  colann.fit <- rbind(colann.fit1, colann.fit2)
  col.pseudotime = colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(length(unique(colann.fit$pseudotime)))
  names(col.pseudotime) = unique(colann.fit$pseudotime)
  annotation_colors$pseudotime <- col.pseudotime
  fit.scale <- cbind(fit.scale, FitDiff.scale)
  rownames(colann.fit) = colnames(fit.scale)
  if (!is.na(sep)) {
    rownames(fit.scale) <- sub(sep, "", rownames(fit.scale))
  }
  p2data <- fit.scale[, rownames(colann.fit)[colann.fit[, 
                                                        "expression"] != "ModeledGroupDiff"]]
  p2data[p2data > quantile(as.vector(p2data), 0.99, na.rm = TRUE)] <- quantile(as.vector(p2data), 
                                                                               0.99, na.rm = TRUE)
  p2data[p2data < quantile(as.vector(p2data), 0.01, na.rm = TRUE)] <- quantile(as.vector(p2data), 
                                                                               0.01, na.rm = TRUE)
  col_fun = colorRamp2(seq(-max(abs(p2data)), max(abs(p2data)), 
                           length.out = 50), colorRampPalette(c("blue3", "skyblue", 
                                                                "white", "pink", "red3"))(50))
  pt_col_fun = colorRamp2(seq(1, length(unique(colann.fit$pseudotime))), 
                          colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(length(unique(colann.fit$pseudotime))))
  annotation_colors$pseudotime <- pt_col_fun
  ht2 <- Heatmap(p2data, cluster_rows = FALSE, cluster_columns = FALSE, 
                 show_row_names = showRowName, show_column_names = FALSE, 
                 heatmap_legend_param = list(legend_direction = "horizontal"), 
                 col = col_fun, top_annotation = HeatmapAnnotation(df = colann.fit[colnames(p2data), 
                 ], col = annotation_colors, show_annotation_name = FALSE), 
                 width = 2)
  p3data <- fit.scale[, rownames(colann.fit)[colann.fit[, 
                                                        "expression"] == "ModeledGroupDiff"]]
  p3data[p3data > quantile(as.vector(p3data), 0.99, na.rm = TRUE)] <- quantile(as.vector(p3data), 
                                                                               0.99, na.rm = TRUE)
  p3data[p3data < quantile(as.vector(p3data), 0.01, na.rm = TRUE)] <- quantile(as.vector(p3data), 
                                                                               0.01, na.rm = TRUE)
  col_fun = colorRamp2(seq(-max(abs(p3data)), max(abs(p3data)), 
                           length.out = 50), colorRampPalette(c("blue3", "skyblue", 
                                                                "white", "pink", "red3"))(50))
  ht3 <- Heatmap(p3data, cluster_rows = FALSE, cluster_columns = FALSE, 
                 show_row_names = showRowName, show_column_names = FALSE, 
                 heatmap_legend_param = list(legend_direction = "horizontal"), 
                 col = col_fun, top_annotation = HeatmapAnnotation(df = colann.fit[colnames(p3data), 
                 ], col = annotation_colors, show_annotation_name = FALSE), 
                 width = 1)
  p4data <- fit.scale[, rownames(colann.fit)[colann.fit[, 
                                                        "expression"] == "ModeledGroupDiff"]]
  p4data <- (p4data - rowMeans(p4data))/apply(p4data, 1, sd)
  p4data[p4data > quantile(as.vector(p4data), 0.99, na.rm = TRUE)] <- quantile(as.vector(p4data), 
                                                                               0.99, na.rm = TRUE)
  p4data[p4data < quantile(as.vector(p4data), 0.01, na.rm = TRUE)] <- quantile(as.vector(p4data), 
                                                                               0.01, na.rm = TRUE)
  # p4data[rownames(p4data) %in% sub(":.*", "", rownames(rowann)[rowann$DDGType == 
  #                                                                "meanSig"]), ] <- 0
  col_fun = colorRamp2(seq(-max(abs(p4data)), max(abs(p4data)), 
                           length.out = 50), colorRampPalette(c("blue3", "skyblue", 
                                                                "white", "pink", "red3"))(50))
  ht4 <- Heatmap(p4data, cluster_rows = FALSE, cluster_columns = FALSE, 
                 show_row_names = showRowName, show_column_names = FALSE, 
                 heatmap_legend_param = list(legend_direction = "horizontal"), 
                 col = col_fun, top_annotation = HeatmapAnnotation(df = colann.fit[colnames(p4data), 
                 ], col = annotation_colors, show_annotation_name = FALSE), 
                 width = 1)
  rownames(alluniformdiff) <- sub(":.*", "", rownames(alluniformdiff))
  p5data <- rowMeans(alluniformdiff[rownames(fit.scale), , 
                                    drop = FALSE]) %*% matrix(1, nrow = 1, ncol = ncol(p4data))
  rownames(p5data) <- rownames(p4data)
  p5data[rownames(p5data) %in% sub(":.*", "", rownames(rowann)[rowann$DDGType == 
                                                                 "trendSig"]), ] <- 0
  col_fun = colorRamp2(seq(-max(abs(p5data)), max(abs(p5data)), 
                           length.out = 50), colorRampPalette(c("blue3", "skyblue", 
                                                                "white", "pink", "red3"))(50))
  ht5 <- Heatmap(p5data, cluster_rows = FALSE, cluster_columns = FALSE, 
                 show_row_names = showRowName, show_column_names = FALSE, 
                 heatmap_legend_param = list(legend_direction = "horizontal"), 
                 col = col_fun, width = 1)
  htlist <- ht1 + ht2 + ht4 + ht5
  draw(htlist, merge_legend = FALSE, annotation_legend_side = "right", 
       heatmap_legend_side = "bottom")
}


pdf("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_lineage_2_injury_control_SCFA_diet_Lamian_diffgene.pdf",width = 15, height = 10)
plotXDEHm_edit(Res, cellWidthTotal = 100, cellHeightTotal = 100, subsampleCell = F, sep = ':.*',showRowName = T)
# save plot 10*15 Microglia_lineage_1/2/3 
dev.off()


a = plotXDEHm_edit(Res, cellWidthTotal = 100, cellHeightTotal = 100, subsampleCell = F, sep = ':.*',showRowName = T)

genes = environment(a@layout[["graphic_fun_list"]][[1]])[["ht"]]@row_names_param[["labels"]]
cluster = environment(a@layout[["graphic_fun_list"]][[1]])[["ht_main"]]@left_annotation@anno_list[["cluster"]]@fun@var_env[["value"]]

Gene_cluster = data.frame(Gene = genes,
                          cluster = cluster)
write.csv(Gene_cluster,file = "/ix/gkohanbash/SimonsLab/R_Share/table/No_acetate_group/10_Myeloid_lineage_2_Lamian_gene.csv",row.names = F,col.names = T,quote = F)

gene <- Gene_cluster$Gene[Gene_cluster$cluster %in% c(1,2,3,5)]
go <- enrichGO(gene = gene, OrgDb = "org.Mm.eg.db", ont="BP",qvalueCutoff = 0.05,keyType = "SYMBOL")

p3 = cnetplot(go,layout = "circle")
p4 = goplot(go,showCategory = 5)
cowplot::plot_grid(p3,p4,rel_widths = c(1.5,2))
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_lineage_2_injury_control_SCFA_diet_Lamian_diffgene_cluster_1_2_3_5_GO.pdf",width = 16,height = 8)

# for neurodegenerative disease
gene = read.csv("/ix/gkohanbash/SimonsLab/R_Share/table/No_acetate_group/10_Myeloid_lineage_2_Lamian_gene.csv")
gene = gene$Gene[gene$cluster %in% c(3,5)]

trans <- bitr(gene, fromType="SYMBOL",
              toType="ENTREZID", OrgDb="org.Mm.eg.db")
kk <- enrichKEGG(gene= trans$ENTREZID,
                 organism= 'mmu',
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 1)

kk@result$p.adjust = kk@result$pvalue

edox2 <- pairwise_termsim(kk)
treeplot(edox2, hclust_method = "average")

kk@result = kk@result[c(1,9,11,13,14),]
dotplot(kk)
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_lineage_2_injury_control_SCFA_diet_Lamian_diffgene_cluster_3_5_KEGG_Neurodegenerative_disease_enrich.pdf",
       height = 6,width = 15)

# SCENIC
library(Seurat)
library(SCENIC)
library(pkgmaker)
library(foreach)
library(pheatmap)
library(scales)

dir.create("/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/SCENIC/Myeloid")

setwd("/ix/gkohanbash/SimonsLab/R_Share/processing/SCENIC/Myeloid")

intermediate_file_path = "/ix/gkohanbash/SimonsLab/R_Share/processing/SCENIC/Myeloid/"
intermediate_file_path_new = "/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/SCENIC/Myeloid/"
db_dir = "/ix/gkohanbash/Zujuan_dhr11/TBI/ref/SCENIC/cisTarget_databases"
organism = "mgi"

# use SCT adjusted counts as input
input_matrix = GetAssayData(MB.Myeloid,slot = "counts") #data.frame (#row = gene, column = sample, rownames(input) = gene_symbol, colnames(input) = sample_ID)
input_matrix = as.data.frame(input_matrix)
anno = MB.Myeloid@meta.data
cellInfo <- as.character(MB.Myeloid$subtype) # interset metadata that will be grouped in the SCENIC analysis
col = hue_pal()(16)
cluster = col[1:length(unique(MB.Myeloid$subtype))]
names(cluster) = unique(MB.Myeloid$subtype)
colVars <- list(cluster=cluster)

thread = 16
dir.create(paste0(intermediate_file_path_new,"fig"))

#1.analysis

#input data reading
exprMat <- as.matrix(input_matrix)

#create cell info
cellInfo <- data.frame(cellInfo) #such as cell seurat_clusters, cell type, nGene, nUMI
colnames(cellInfo) <- c("cluster")
rownames(cellInfo)<- rownames(anno)
cbind(table(cellInfo$cluster))

dir.create(paste0(intermediate_file_path,"int"))
saveRDS(cellInfo, file=paste0(intermediate_file_path,"int/cellInfo.Rds"))

#cell cluster color
colVars$cluster <- colVars$cluster[intersect(names(colVars$cluster), cellInfo$cluster)]

#initiation step i
scenicOptions <- initializeScenic(org = organism, 
                                  dbDir = db_dir, 
                                  nCores = thread)#create SCENIC input

scenicOptions@inputDatasetInfo$cellInfo <- cellInfo
scenicOptions@inputDatasetInfo$colVars <- colVars

saveRDS(scenicOptions, file=paste0(intermediate_file_path,"int/scenicOptions.Rds"))

genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)

exprMat_filtered <- exprMat[genesKept, ]

runCorrelation(exprMat_filtered, scenicOptions)
#Spearman correlation between the TF and the potential target

runGenie3(exprMat_filtered, scenicOptions
          # ,resumePreviousRun = T  # Reload partial results from a previous (interrupted) run and resume execution
) #耗CPU

save(scenicOptions,file = paste0(intermediate_file_path,'Genie3.RData'))

scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings

runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions) # Toy run settings
scenicOptions@settings[["nCores"]] = 8
runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)
runSCENIC_4_aucell_binarize(scenicOptions)

save(scenicOptions,file = paste0(intermediate_file_path,'runSCENIC.RData'))

# heatmap for regulon
load(paste0(intermediate_file_path,'runSCENIC.RData'))
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonAUC <- regulonAUC[,colnames(regulonAUC) %in% rownames(cellInfo)]
cellInfo$sample <- rownames(cellInfo)
cellInfo <- cellInfo[cellInfo$sample %in% colnames(regulonAUC),]

regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$cluster),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
regulonActivity_byCellType_Scaled = regulonActivity_byCellType_Scaled[-grep("extend",rownames(regulonActivity_byCellType_Scaled)),]

pheatmap::pheatmap(regulonActivity_byCellType_Scaled, 
                   #fontsize_row=3, 
                   color=colorRampPalette(c("blue","white","red"))(100), 
                   breaks=seq(-3, 3, length.out = 100),
                   treeheight_row=10, 
                   treeheight_col=10,
                   border_color=NA)
# save /ix/gkohanbash/SimonsLab/R_Share/fig/10_Myeloid_subtype_Regulon_activity_heatmap 11*8

# downstream analysis of SCENIC
library(scFunctions)
library(launcheR)
library(SCENIC)
library(AUCell)

# for groups
metadata_sub = MB.Myeloid@meta.data
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
binaryRegulonActivity = as.data.frame(binaryRegulonActivity)
binaryRegulonActivity = binaryRegulonActivity[,colnames(MB.Myeloid)]
binaryRegulonActivity = as.matrix(binaryRegulonActivity)

regulons <- loadInt(scenicOptions, "regulons")
Regulon <- regulonAUC@assays@data@listData[["AUC"]]
Regulon <- as.data.frame(Regulon)
Regulon = Regulon[,colnames(MB.Myeloid)]

binary_regulons_trans <- as.matrix(binaryRegulonActivity)

# calculate regulon specificity score (RSS) based on Jensen-Shannon divergence
# need metadata has a column named cell_type
rrs_df = calculate_rrs(metadata_sub,
                       binary_regulons = binary_regulons_trans,
                       cell_type_column = "group") # the name of the type which you want to calculate the RSS for

# visualization
plot_rrs_ranking(rrs_df,
                 cell_type = "all",    # the cell type name
                 ggrepel_force = 1,
                 ggrepel_point_padding = 0.2,
                 top_genes = 30,
                 plot_extended = FALSE) # whether you'd like to plot the extend regulon which is low confidence regulon-relation or contains genes based on motif prediction
# save  /ix/gkohanbash/SimonsLab/R_Share/fig/10_Myeloid_SCENIC_regulon_RSS_rank_for_groups 6*12

# visualize the distribution of RSS on all cell types
library(ggridges)

rrs_df_nona <- subset(rrs_df,RSS > 0)
ggplot(rrs_df_nona,aes(RSS,cell_type, fill = cell_type)) +
  geom_density_ridges(scale = 5, alpha = 0.75) +
  geom_vline(xintercept = 0.1) +
  theme(legend.position = "none")

# heatmap of the specific regulon on each cell type
rrs_df_wide <- rrs_df %>%
  spread(cell_type,RSS)
rownames(rrs_df_wide) <- rrs_df_wide$regulon 
rrs_df_wide <- rrs_df_wide[,2:ncol(rrs_df_wide)]

## Subset all regulons that don't have at least an RSS of 0.7 for one cell type
a = as.matrix(rrs_df_wide)
a = as.vector(a)
b = data.frame(RRS = a,
               group = c(rep("injury_control",nrow(rrs_df_wide)),
                         rep("injury_SCFA",nrow(rrs_df_wide)),
                         rep("sham_control",nrow(rrs_df_wide))),
               regulon = rep(rownames(rrs_df_wide),ncol(rrs_df_wide)))
b = b[-grep("extend",b$regulon),]
top5 <- b %>% group_by(group) %>% top_n(5, RRS)

rrs_df_wide_specific <- rrs_df_wide[unique(top5$regulon),] # filter the RSS > 0.4 as the specific threshold for specific regulon

library(heatmaply)
heatmaply(rrs_df_wide_specific)
pheatmap(rrs_df_wide_specific,color = viridis(n=256, alpha = 1, begin = 0, end = 1, option = "viridis"),border_color =NA,
         angle_col = 0)

# save 10_Myeloid_group_specific_regulon_heatmap

# # calculating the connection specificity index (CSI) for all regulons

# The CSI is a major of connectedness between the different regulons. 
# Regulons that share high CSI likely are co-regulating downstream genes and are together responsible for cell function. 

calculate_csi = function (regulonAUC, calc_extended = FALSE, verbose = FALSE) 
{
  compare_pcc <- function(vector_of_pcc, pcc) {
    pcc_larger <- length(vector_of_pcc[vector_of_pcc > pcc])
    if (pcc_larger == length(vector_of_pcc)) {
      return(0)
    }
    else {
      return(length(vector_of_pcc))
    }
  }
  calc_csi <- function(reg, reg2, pearson_cor) {
    test_cor <- pearson_cor[reg, reg2]
    total_n <- ncol(pearson_cor)
    pearson_cor_sub <- subset(pearson_cor, rownames(pearson_cor) == 
                                reg | rownames(pearson_cor) == reg2)
    sums <- apply(pearson_cor_sub, MARGIN = 2, FUN = compare_pcc, 
                  pcc = test_cor)
    fraction_lower <- length(sums[sums == nrow(pearson_cor_sub)])/total_n
    return(fraction_lower)
  }
  regulonAUC_sub <- regulonAUC
  if (calc_extended == TRUE) {
    regulonAUC_sub <- subset(regulonAUC_sub, grepl("extended", 
                                                   rownames(regulonAUC_sub)))
  }
  else if (calc_extended == FALSE) {
    regulonAUC_sub <- subset(regulonAUC_sub, !grepl("extended", 
                                                    rownames(regulonAUC_sub)))
  }
  regulonAUC_sub <- t(regulonAUC_sub)
  pearson_cor <- cor(regulonAUC_sub)
  pearson_cor_df <- as.data.frame(pearson_cor)
  pearson_cor_df$regulon_1 <- rownames(pearson_cor_df)
  pearson_cor_long <- pearson_cor_df %>% gather(regulon_2, 
                                                pcc, -regulon_1) %>% mutate(regulon_pair = paste(regulon_1, 
                                                                                                 regulon_2, sep = "_"))
  regulon_names <- unique(colnames(pearson_cor))
  num_of_calculations <- length(regulon_names) * length(regulon_names)
  csi_regulons <- data.frame(matrix(nrow = num_of_calculations, 
                                    ncol = 3))
  colnames(csi_regulons) <- c("regulon_1", "regulon_2", "CSI")
  num_regulons <- length(regulon_names)
  f <- 0
  for (reg in regulon_names) {
    if (verbose == TRUE) {
      print(reg)
    }
    for (reg2 in regulon_names) {
      f <- f + 1
      fraction_lower <- calc_csi(reg, reg2, pearson_cor)
      csi_regulons[f, ] <- c(reg, reg2, fraction_lower)
    }
  }
  csi_regulons$CSI <- as.numeric(csi_regulons$CSI)
  return(csi_regulons)
}

regulons_csi <- calculate_csi(Regulon,
                              calc_extended = FALSE)

# plot heatmap for regulon module
library(pheatmap)

k = 6 # the cluster numbers

plot_csi_modules = function (csi_df, nclust = 10, font_size_regulons = 6) 
{
  csi_test_mat <- csi_df %>% spread(regulon_2, CSI)
  future_rownames <- csi_test_mat$regulon_1
  csi_test_mat <- as.matrix(csi_test_mat[, 2:ncol(csi_test_mat)])
  rownames(csi_test_mat) <- future_rownames
  pheatmap(csi_test_mat, show_colnames = FALSE, color = viridis(n = 10), 
           cutree_cols = nclust, cutree_rows = nclust, fontsize_row = font_size_regulons, 
           cluster_cols = TRUE, cluster_rows = TRUE, treeheight_row = nclust, 
           treeheight_col = 10, clustering_distance_rows = "euclidean", 
           clustering_distance_cols = "euclidean", widt = 2000, 
           height = 3200)
}

# plot heatmap based on paired CSI between regulon pairs
plot_csi_modules(regulons_csi,
                 nclust = k,   # the cluster number
                 font_size_regulons = 8)
# save 10_Myeloid_Regulon_module_group 9*9

# get the CSI module regulon dataframe
csi_csi_wide <- regulons_csi %>%
  spread(regulon_2,CSI)

future_rownames <- csi_csi_wide$regulon_1
csi_csi_wide <- as.matrix(csi_csi_wide[,2:ncol(csi_csi_wide)])
rownames(csi_csi_wide) <- future_rownames

regulons_hclust <- hclust(dist(csi_csi_wide,method = "euclidean"))

clusters <- cutree(regulons_hclust,k = k)

# adjust cluster number to keep the same with the heatmap order
clusters = mapvalues(clusters, "1","Module_1")
clusters = mapvalues(clusters, "2","Module_3")
clusters = mapvalues(clusters, "3","Module_4")
clusters = mapvalues(clusters, "4","Module_5")
clusters = mapvalues(clusters, "5","Module_6")
clusters = mapvalues(clusters, "6","Module_2")

clusters_df <- data.frame("regulon" = names(clusters),
                          "csi_cluster" = clusters)

# enrich GO pathway for TF target gene of each module to annotate the module function
## Database: TRRUST v2
## cite: Han H, Cho JW, Lee S, Yun A, Kim H, Bae D, Yang S, Kim CY, Lee M, Kim E, Lee S, Kang B, Jeong D, Kim Y, Jeon HN, Jung H, Nam S, Chung M, Kim JH, Lee I. TRRUST v2: an expanded reference database of human and mouse transcriptional regulatory interactions. Nucleic Acids Res. 2018 Jan 4;46(D1):D380-D386. doi: 10.1093/nar/gkx1013. PMID: 29087512; PMCID: PMC5753191.
tf_target = read.delim("/ix/gkohanbash/Zujuan_dhr11/refs/Murine_TF_target_TRRUST/trrust_rawdata.mouse.tsv",header = F)
colnames(tf_target) = c("TF","Target","Direction","ID")

for (i in unique(clusters_df$csi_cluster)) {
  tf = clusters_df$regulon[clusters_df$csi_cluster == i]
  tf = gsub(" \\(.*","",tf)
  
  target = tf_target[tf_target$TF %in% tf,]
  
  activation = target$Target[target$Direction == "Activation"]
  go.a <- enrichGO(gene = activation, OrgDb = "org.Mm.eg.db", ont="BP",qvalueCutoff = 0.05,keyType = "SYMBOL",pvalueCutoff = 0.05)
  edox2 <- pairwise_termsim(go.a)
  treeplot(edox2, hclust_method = "average")
  ggsave(paste0("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_TF_",i,"_activated_target_gene_GO_treeplot.pdf"),width = 20,height = 8)
  
  repression = target$Target[target$Direction == "Repression"]
  go.r <- enrichGO(gene = repression, OrgDb = "org.Mm.eg.db", ont="BP",qvalueCutoff = 0.05,keyType = "SYMBOL",pvalueCutoff = 0.05)
  edox2 <- pairwise_termsim(go.r)
  treeplot(edox2, hclust_method = "average")
  ggsave(paste0("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/10_Myeloid_TF_",i,"_repressed_target_gene_GO_treeplot.pdf"),width = 20,height = 8)
}

# find module activate/repress cellchat outgoing signal
library(CellChat)
CellChatDB <- CellChatDB.mouse
ligand = data.frame(ligand = CellChatDB[["interaction"]][["ligand"]],
                    pathway = CellChatDB[["interaction"]][["pathway_name"]])


tf_target = read.delim("/ix/gkohanbash/Zujuan_dhr11/refs/Murine_TF_target_TRRUST/trrust_rawdata.mouse.tsv",header = F)
colnames(tf_target) = c("TF","ligand","Direction","ID")

tf = tf_target[tf_target$ligand %in% ligand$ligand,]
tf = left_join(tf,ligand,by = "ligand")
tf = tf[!duplicated(tf),]

cdf = clusters_df$regulon[clusters_df$csi_cluster == "Module_6"]
cdf = gsub(" \\(.*","",cdf)

p = tf[tf$TF %in% cdf,]
signal_pathway = table(p$pathway[p$Direction == "Activation"])
# Activation Repression
b = names(signal_pathway)
a = as.matrix(signal_pathway)

col = c(RColorBrewer::brewer.pal(12,name = "Paired"),RColorBrewer::brewer.pal(9,name = "Set1"),'black')
wordcloud::wordcloud(words = b, freq = a,
                     colors = col, random.color = TRUE, rot.per = 0.3, min.freq = 2)

# save 10_Myeloid_TF_Module_6_Activate_downstream_cellchat_outgoing_signal_pathway 5*7

signal_pathway = table(p$pathway[p$Direction == "Repression"])
# Activation Repression
b = names(signal_pathway)
a = as.matrix(signal_pathway)

col = c(RColorBrewer::brewer.pal(12,name = "Paired"),RColorBrewer::brewer.pal(9,name = "Set1"),'black')
wordcloud::wordcloud(words = b, freq = a,
                     colors = col, random.color = TRUE, rot.per = 0.3, min.freq = 2)
# save 10_Myeloid_TF_Module_6_Repress_downstream_cellchat_outgoing_signal_pathway 4*6

# global statistic in CSI module
# Check how many regulons are in each cluster
clusters_df_stats <- clusters_df %>%
  group_by(csi_cluster) %>%
  mutate("regulon" = as.character(regulon)) %>%
  tally()

ggplot(clusters_df_stats,aes(as.factor(csi_cluster),n,fill=as.factor(csi_cluster))) +
  geom_bar(color= "black",stat="identity") +
  theme(legend.position="none") +
  scale_fill_brewer(palette = "Set3") +
  labs(x = "HC clusters",
       y = "# Regulons")

## Check average regulon size per cluster
clusters_df_regsizes <- clusters_df %>%
  separate(regulon, into = c("regulon_name","regulon_size"), sep=" ") %>%
  mutate("regulon_size" = gsub("\\(","",regulon_size)) %>%
  mutate("regulon_size" = gsub("\\g)","",regulon_size)) %>%
  mutate("regulon_size" = as.numeric(regulon_size))

ggplot(clusters_df_regsizes,aes(log10(regulon_size),as.factor(csi_cluster),fill=as.factor(csi_cluster))) + 
  geom_density_ridges() +
  scale_fill_brewer(palette = "Set3") +
  theme(legend.position = "none")
## Picking joint bandwidth of 0.213

## Plot correlation between number of regulons and regulon size
library(ggrepel)

clusters_df_regsizes_summary <- clusters_df_regsizes %>%
  group_by(csi_cluster) %>%
  summarise("mean_regulon_size" = mean(regulon_size),
            "median_regulon_size" = median(regulon_size),
            "sd_regulon_size" = sd(regulon_size))

clusters_meta <-  full_join(clusters_df_stats,clusters_df_regsizes_summary,by="csi_cluster")

ggplot(clusters_meta,aes(n,log10(median_regulon_size),label=as.factor(csi_cluster),fill=as.factor(csi_cluster))) +
  geom_point(color = "black",pch = 21) +
  geom_label_repel() +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set3")

# heatmap of CSI module activity in each cell type
calc_csi_module_activity = function (clusters_df, regulonAUC, metadata, cell_type_column) 
{
  metadata$cell_type <- metadata[, cell_type_column]
  cell_types <- unique(metadata$cell_type)
  regulons <- unique(clusters_df$regulon)
  regulonAUC_sub <- regulonAUC
  regulonAUC_sub <- regulonAUC_sub[regulons, ]
  csi_activity_matrix_list <- list()
  csi_cluster_activity <- data.frame(csi_cluster = c(), mean_activity = c(), 
                                     cell_type = c())
  cell_type_counter <- 0
  regulon_counter <- for (ct in cell_types) {
    cell_type_counter <- cell_type_counter + 1
    cell_type_aucs <- rowMeans(regulonAUC_sub[, rownames(subset(metadata, 
                                                                cell_type == ct))])
    cell_type_aucs_df <- data.frame(regulon = names(cell_type_aucs), 
                                    activtiy = cell_type_aucs, cell_type = ct)
    csi_activity_matrix_list[[ct]] <- cell_type_aucs_df
  }
  for (ct in names(csi_activity_matrix_list)) {
    for (cluster in unique(clusters_df$csi_cluster)) {
      csi_regulon <- subset(clusters_df, csi_cluster == 
                              cluster)
      csi_regulon_activtiy <- subset(csi_activity_matrix_list[[ct]], 
                                     regulon %in% csi_regulon$regulon)
      csi_activtiy_mean <- mean(csi_regulon_activtiy$activtiy)
      this_cluster_ct_activity <- data.frame(csi_cluster = cluster, 
                                             mean_activity = csi_activtiy_mean, cell_type = ct)
      csi_cluster_activity <- rbind(csi_cluster_activity, 
                                    this_cluster_ct_activity)
    }
  }
  csi_cluster_activity[is.na(csi_cluster_activity)] <- 0
  csi_cluster_activity_wide <- csi_cluster_activity %>% spread(cell_type, 
                                                               mean_activity)
  rownames(csi_cluster_activity_wide) <- csi_cluster_activity_wide$csi_cluster
  csi_cluster_activity_wide <- as.matrix(csi_cluster_activity_wide[2:ncol(csi_cluster_activity_wide)])
  return(csi_cluster_activity_wide)
}

csi_cluster_activity_wide <- calc_csi_module_activity(clusters_df,
                                                      Regulon,
                                                      metadata_sub,
                                                      "group")

pheatmap(csi_cluster_activity_wide,
         show_colnames = TRUE,
         color = viridis(n = 10),
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         border_color = NA,
         angle_col = 0)
# save 10_Myeloid_Regulon_module_group_heatmap 4*6

# for subtype
metadata_sub = MB.Myeloid@meta.data

binary_regulons_trans <- as.matrix(binaryRegulonActivity)

# calculate regulon specificity score (RSS) based on Jensen-Shannon divergence
# need metadata has a column named cell_type
rrs_df = calculate_rrs(metadata_sub,
                       binary_regulons = binary_regulons_trans,
                       cell_type_column = "subtype") # the name of the type which you want to calculate the RSS for

# visualization
plot_rrs_ranking(rrs_df,
                 cell_type = "all",    # the cell type name
                 ggrepel_force = 1,
                 ggrepel_point_padding = 0.2,
                 top_genes = 40,
                 plot_extended = FALSE) # whether you'd like to plot the extend regulon which is low confidence regulon-relation or contains genes based on motif prediction
#save /ix/gkohanbash/SimonsLab/R_Share/fig/10_Myeloid_SCENIC_regulon_RSS_rank_for_subtypes 10*10

# visualize the distribution of RSS on all cell types
library(ggridges)

rrs_df_nona <- subset(rrs_df,RSS > 0)
ggplot(rrs_df_nona,aes(RSS,cell_type, fill = cell_type)) +
  geom_density_ridges(scale = 5, alpha = 0.75) +
  geom_vline(xintercept = 0.1) +
  theme(legend.position = "none")

# heatmap of the specific regulon on each cell type
rrs_df_wide <- rrs_df %>%
  spread(cell_type,RSS)
rownames(rrs_df_wide) <- rrs_df_wide$regulon 
rrs_df_wide <- rrs_df_wide[,2:ncol(rrs_df_wide)]

## Subset all regulons that don't have at least an RSS of 0.7 for one cell type
a = as.matrix(rrs_df_wide)
a = as.vector(a)
b = data.frame(RRS = a,
               group = c(rep("MAC",nrow(rrs_df_wide)),
                         rep("MG_1",nrow(rrs_df_wide)),
                         rep("MG_10",nrow(rrs_df_wide)),
                         rep("MG_2",nrow(rrs_df_wide)),
                         rep("MG_3",nrow(rrs_df_wide)),
                         rep("MG_4",nrow(rrs_df_wide)),
                         rep("MG_5",nrow(rrs_df_wide)),
                         rep("MG_6",nrow(rrs_df_wide)),
                         rep("MG_7",nrow(rrs_df_wide)),
                         rep("MG_8",nrow(rrs_df_wide)),
                         rep("MG_9",nrow(rrs_df_wide))),
               regulon = rep(rownames(rrs_df_wide),ncol(rrs_df_wide)))
b = b[-grep("extend",b$regulon),]
top5 <- b %>% group_by(group) %>% top_n(5, RRS)

rrs_df_wide_specific <- rrs_df_wide[unique(top5$regulon),]

library(heatmaply)
heatmaply(rrs_df_wide_specific)
pheatmap(rrs_df_wide_specific,color = viridis(n=256, alpha = 1, begin = 0, end = 1, option = "viridis"),border_color =NA,
         angle_col = 0)

# save 10_Myeloid_subtype_specific_regulon_heatmap 7*14

# plot heatmap for regulon module in each subtype

clusters_df <- data.frame("regulon" = names(clusters),
                          "csi_cluster" = clusters)

csi_cluster_activity_wide <- calc_csi_module_activity(clusters_df,
                                                      Regulon,
                                                      metadata_sub,
                                                      "subtype")

pheatmap(csi_cluster_activity_wide,
         show_colnames = TRUE,
         color = viridis(n = 10),
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         border_color = NA,
         angle_col = 90)
# save 10_Myeloid_Regulon_module_subtype_heatmap 4*6

write.csv(clusters_df,file = "/ix/gkohanbash/SimonsLab/R_Share/table/No_acetate_group/10_Myeloid_regulon_CSI_cluster.csv")

# check long-term TBI consequence caused by myeloid cells
load("/ix/gkohanbash/SimonsLab/R_Share/processing/10_Myeloid_Seurat_object.RData")
load("/ix/gkohanbash/SimonsLab/R_Share/processing/10_Myeloid_GO_AUC.RData")
Myeloid_GO_AUC = cells_AUC@assays@data@listData[["AUC"]]
Myeloid_GO_AUC = Myeloid_GO_AUC[,colnames(MB.Myeloid)]
G_item = as.data.frame(rownames(Myeloid_GO_AUC))


load("/ix/gkohanbash/SimonsLab/R_Share/processing/10_Myeloid_KEGG_AUC.RData")
Myeloid_KEGG_AUC = cells_AUC@assays@data@listData[["AUC"]]
Myeloid_KEGG_AUC = Myeloid_KEGG_AUC[,colnames(MB.Myeloid)]
K_item = as.data.frame(rownames(Myeloid_KEGG_AUC))

load("/ix/gkohanbash/SimonsLab/R_Share/processing/5_GO_AUC_genelist.RData")
load("/ix/gkohanbash/SimonsLab/R_Share/processing/5_KEGG_AUC_genelist.RData")

exprMatrix <- MB.Myeloid@assays$SCT@counts # counts
exprMatrix <- as.matrix(exprMatrix)
logMat <- log2(exprMatrix+2)

## Neurodegenerative disease
geneSetName = c("GOBP_MICROGLIAL_CELL_ACTIVATION","GOBP_POSITIVE_REGULATION_OF_MICROGLIAL_CELL_ACTIVATION","GOBP_WOUND_HEALING",
                "GOBP_MACROPHAGE_ACTIVATION","GOBP_POSITIVE_REGULATION_OF_MACROPHAGE_ACTIVATION")

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)

# violin plot in different groups - MAC
vln.df=as.data.frame(t(meanByGs))
# vln.df = as.data.frame(t(Myeloid_GO_AUC[geneSetName,]))
vln.df$group = MB.Myeloid$group
vln.df$subtype = MB.Myeloid$subtype
vln.df = vln.df[vln.df$subtype == "MAC",]

group=levels(factor(vln.df$group))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  
dodge <- position_dodge(width = 1)
vln.df%>%ggplot(aes(group,GOBP_MACROPHAGE_ACTIVATION))+geom_violin(aes(fill=group),scale = "width",position = dodge)+
  geom_boxplot(aes(fill=group),position = dodge,outlier.size = 0,notch = T,width = 0.3)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "right"
  )+
  stat_compare_means(comparisons = my_comparisons,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),label = "p.signif")+
  ggtitle("GOBP_MACROPHAGE_ACTIVATION")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/long-term_TBI_consquence/10_Macrophage_activation_groups_comparison.pdf",width = 5,height = 5)

vln.df%>%ggplot(aes(group,GOBP_WOUND_HEALING))+geom_violin(aes(fill=group),scale = "width",position = dodge)+
  geom_boxplot(aes(fill=group),position = dodge,outlier.size = 0,notch = T,width = 0.3)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "right"
  )+
  stat_compare_means(comparisons = my_comparisons,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),label = "p.signif")+
  ggtitle("GOBP_WOUND_HEALING")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/long-term_TBI_consquence/10_Macrophage_GOBP_WOUND_HEALING_groups_comparison.pdf",width = 5,height = 5)

## axonogenesis inhibition
geneSetName = c("REACTOME_REGULATION_OF_EXPRESSION_OF_SLITS_AND_ROBOS")

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)

# violin plot in different groups
vln.df=as.data.frame(t(meanByGs))
# vln.df = as.data.frame(t(Astrocyte_GO_AUC[geneSetName,]))
vln.df$group = MB.Myeloid$group
vln.df$subtype = MB.Myeloid$subtype
vln.df = vln.df[vln.df$subtype == "MAC",]

group=levels(factor(vln.df$group))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

dodge <- position_dodge(width = 1)
vln.df%>%ggplot(aes(group,REACTOME_REGULATION_OF_EXPRESSION_OF_SLITS_AND_ROBOS))+geom_violin(aes(fill=group),scale = "width",position = dodge)+
  geom_boxplot(aes(fill=group),position = dodge,outlier.size = 0,notch = T,width = 0.2)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "right"
  )+
  stat_compare_means(comparisons = my_comparisons,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),label = "p.signif")+
  ggtitle("REACTOME_REGULATION_OF_EXPRESSION_OF_SLITS_AND_ROBOS")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/long-term_TBI_consquence/10_MAC_REACTOME_REGULATION_OF_EXPRESSION_OF_SLITS_AND_ROBOS_groups_comparison.pdf",width = 5,height = 5)

# violin plot in different groups - MG
vln.df=as.data.frame(t(meanByGs))
# vln.df = as.data.frame(t(Myeloid_GO_AUC[geneSetName,]))
vln.df$group = MB.Myeloid$group
vln.df$subtype = MB.Myeloid$subtype
vln.df = vln.df[vln.df$subtype != "MAC",]

group=levels(factor(vln.df$group))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

dodge <- position_dodge(width = 1)
vln.df%>%ggplot(aes(group,GOBP_MICROGLIAL_CELL_ACTIVATION))+geom_violin(aes(fill=group),scale = "width",position = dodge)+
  geom_boxplot(aes(fill=group),position = dodge,outlier.size = 0,notch = T,width = 0.3)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "right"
  )+
  stat_compare_means(comparisons = my_comparisons,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),label = "p.signif")+
  ggtitle("GOBP_MICROGLIAL_CELL_ACTIVATION")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/long-term_TBI_consquence/10_Microglia_activation_groups_comparison.pdf",width = 5,height = 5)

## axonogenesis inhibition
geneSetName = c("REACTOME_REGULATION_OF_EXPRESSION_OF_SLITS_AND_ROBOS")

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)

# violin plot in different groups
vln.df=as.data.frame(t(meanByGs))
# vln.df = as.data.frame(t(Astrocyte_GO_AUC[geneSetName,]))
vln.df$group = MB.Myeloid$group
vln.df$subtype = MB.Myeloid$subtype
vln.df = vln.df[vln.df$subtype != "MAC",]

group=levels(factor(vln.df$group))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

dodge <- position_dodge(width = 1)
vln.df%>%ggplot(aes(group,REACTOME_REGULATION_OF_EXPRESSION_OF_SLITS_AND_ROBOS))+geom_violin(aes(fill=group),scale = "width",position = dodge)+
  geom_boxplot(aes(fill=group),position = dodge,outlier.size = 0,notch = T,width = 0.2)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "right"
  )+
  stat_compare_means(comparisons = my_comparisons,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),label = "p.signif")+
  ggtitle("REACTOME_REGULATION_OF_EXPRESSION_OF_SLITS_AND_ROBOS")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/long-term_TBI_consquence/10_MG_REACTOME_REGULATION_OF_EXPRESSION_OF_SLITS_AND_ROBOS_groups_comparison.pdf",width = 5,height = 5)

load("/ix/gkohanbash/SimonsLab/R_Share/processing/10_Myeloid_Seurat_object.RData")
MB.Myeloid$group = mapvalues(MB.Myeloid$group,"sham_control","1_sham_control")
MB.Myeloid$group = mapvalues(MB.Myeloid$group,"injury_control","2_injury_control")
MB.Myeloid$group = mapvalues(MB.Myeloid$group,"injury_SCFA","3_injury_SCFA")

load("/ix/gkohanbash/SimonsLab/R_Share/processing/10_Myeloid_GO_AUC.RData")
Myeloid_GO_AUC = cells_AUC@assays@data@listData[["AUC"]]
Myeloid_GO_AUC = Myeloid_GO_AUC[,colnames(MB.Myeloid)]
G_item = as.data.frame(rownames(Myeloid_GO_AUC))


load("/ix/gkohanbash/SimonsLab/R_Share/processing/10_Myeloid_KEGG_AUC.RData")
Myeloid_KEGG_AUC = cells_AUC@assays@data@listData[["AUC"]]
Myeloid_KEGG_AUC = Myeloid_KEGG_AUC[,colnames(MB.Myeloid)]
K_item = as.data.frame(rownames(Myeloid_KEGG_AUC))

load("/ix/gkohanbash/SimonsLab/R_Share/processing/5_GO_AUC_genelist.RData")

exprMatrix <- MB.Myeloid@assays$SCT@counts # counts
exprMatrix <- as.matrix(exprMatrix)
logMat <- log2(exprMatrix+2)

## Neurodegenerative disease
geneSetName = c("GOBP_INTERLEUKIN_12_PRODUCTION","GOBP_POSITIVE_REGULATION_OF_TRANSFORMING_GROWTH_FACTOR_BETA_PRODUCTION",
                "GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_PRODUCTION","GOBP_INTERLEUKIN_23_PRODUCTION","GOBP_INTERLEUKIN_6_PRODUCTION")

meanByGs <- data.frame()

for (i in 1:length(geneSetName)) {
  for (j in 1:ncol(logMat)) {
    gene = rownames(logMat)[rownames(logMat) %in% h.sets[[geneSetName[i]]]]
    meanByGs[i,j] = sum(logMat[gene,j])/length(h.sets[[geneSetName[i]]])
  }
}

rownames(meanByGs) <- geneSetName
colnames(meanByGs) = colnames(logMat)

# violin plot in different groups
vln.df=as.data.frame(t(meanByGs))
# vln.df = as.data.frame(t(Astrocyte_GO_AUC[geneSetName,]))
vln.df$group = MB.Myeloid$group
vln.df$subtype = MB.Myeloid$subtype
# vln.df = vln.df[vln.df$subtype != "MAC",]

group=levels(factor(vln.df$group))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

dodge <- position_dodge(width = 1)
vln.df%>%ggplot(aes(group,GOBP_INTERLEUKIN_12_PRODUCTION))+geom_violin(aes(fill=group),scale = "width",position = dodge)+
  geom_boxplot(aes(fill=group),position = dodge,outlier.size = 0,notch = T,width = 0.3)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "right"
  )+
  stat_compare_means(comparisons = my_comparisons,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),label = "p.signif")+
  ggtitle("GOBP_INTERLEUKIN_12_PRODUCTION")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/long-term_TBI_consquence/10_Myeloid_GOBP_INTERLEUKIN_12_PRODUCTION_groups_comparison.pdf",width = 5,height = 5)

vln.df%>%ggplot(aes(group,GOBP_POSITIVE_REGULATION_OF_TRANSFORMING_GROWTH_FACTOR_BETA_PRODUCTION))+geom_violin(aes(fill=group),scale = "width",position = dodge)+
  geom_boxplot(aes(fill=group),position = dodge,outlier.size = 0,notch = T,width = 0.3)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "right"
  )+
  stat_compare_means(comparisons = my_comparisons,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),label = "p.signif")+
  ggtitle("GOBP_POSITIVE_REGULATION_OF_TRANSFORMING_GROWTH_FACTOR_BETA_PRODUCTION")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/long-term_TBI_consquence/10_Myeloid_GOBP_POSITIVE_REGULATION_OF_TRANSFORMING_GROWTH_FACTOR_BETA_PRODUCTION_groups_comparison.pdf",width = 5,height = 5)

vln.df%>%ggplot(aes(group,GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_PRODUCTION))+geom_violin(aes(fill=group),scale = "width",position = dodge)+
  geom_boxplot(aes(fill=group),position = dodge,outlier.size = 0,notch = T,width = 0.3)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "right"
  )+
  stat_compare_means(comparisons = my_comparisons,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),label = "p.signif")+
  ggtitle("GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_PRODUCTION")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/long-term_TBI_consquence/10_Myeloid_GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_PRODUCTION_groups_comparison.pdf",width = 5,height = 5)

vln.df%>%ggplot(aes(group,GOBP_INTERLEUKIN_23_PRODUCTION))+geom_violin(aes(fill=group),scale = "width",position = dodge)+
  geom_boxplot(aes(fill=group),position = dodge,outlier.size = 0,notch = T,width = 0.3)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "right"
  )+
  stat_compare_means(comparisons = my_comparisons,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),label = "p.signif")+
  ggtitle("GOBP_INTERLEUKIN_23_PRODUCTION")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/long-term_TBI_consquence/10_Myeloid_GOBP_INTERLEUKIN_23_PRODUCTION_groups_comparison.pdf",width = 5,height = 5)

vln.df%>%ggplot(aes(group,GOBP_INTERLEUKIN_6_PRODUCTION))+geom_violin(aes(fill=group),scale = "width",position = dodge)+
  geom_boxplot(aes(fill=group),position = dodge,outlier.size = 0,notch = T,width = 0.3)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "right"
  )+
  stat_compare_means(comparisons = my_comparisons,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),label = "p.signif")+
  ggtitle("GOBP_INTERLEUKIN_6_PRODUCTION")
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/long-term_TBI_consquence/10_Myeloid_GOBP_INTERLEUKIN_6_PRODUCTION_groups_comparison.pdf",width = 5,height = 5)

