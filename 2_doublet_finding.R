#
##
###
library(DoubletFinder)
library(Seurat)
library(dplyr)

load("/bgfs/gkohanbash/SimonsLab/R_Share/processing/1_integration_input.RData")

MB_list = lapply(MB_list, FUN = function(x){
  # cluster for doublet identification
  x = RunPCA(x,npcs = 50, features = rownames(x))
  pc.num=1:30
  x = RunUMAP(x , reduction = "pca", dims = pc.num)
  # cluster cells
  x = FindNeighbors(x, verbose = FALSE, reduction = "pca",dims = pc.num)
  x = FindClusters(x,resolution = 0.5)
  
  # find doublets
  sweep.res.list <- paramSweep_v3(x, PCs = 1:30, sct = T) # sct = T for SCTransform normalization 
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
  
  # find optimal pK
  bcmvn <- find.pK(sweep.stats) 
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric() # extract optimal pK
  # 0.18
  
  # define homology doublet rate. Because DoubletFinder is weak for homology doublets
  DoubletRate = ncol(x)*8*1e-6 # every 1000 cells doublet rate increase 0.8% 
  homotypic.prop <- modelHomotypic(x$seurat_clusters) 
  nExp_poi <- round(DoubletRate*ncol(x)) 
  # adjust the doublet rate with the homology doublet rate 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  # identify doublet with optimal parameters defined previously
  x <- doubletFinder_v3(x, PCs = 1:30, pN = 0.25, pK = pK_bcmvn, 
                        nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
})


save(MB_list,file = "/bgfs/gkohanbash/SimonsLab/R_Share/processing/2_doublet_identification.RData") 

# annotate doublets
load("/bgfs/gkohanbash/SimonsLab/R_Share/processing/1_integrateData_PCA.RData")
# load("/bgfs/gkohanbash/SimonsLab/R_Share/processing/2_doublet_identification.RData")

doublet = c(MB_list[[1]]@meta.data$DF.classifications_0.25_0.22_121,
            MB_list[[2]]@meta.data$DF.classifications_0.25_0.22_80,
            MB_list[[3]]@meta.data$DF.classifications_0.25_0.11_66,
            MB_list[[4]]@meta.data$DF.classifications_0.25_0.17_124,
            MB_list[[5]]@meta.data$DF.classifications_0.25_0.12_891,
            MB_list[[6]]@meta.data$DF.classifications_0.25_0.27_450,
            MB_list[[7]]@meta.data$DF.classifications_0.25_0.05_912,
            MB_list[[8]]@meta.data$DF.classifications_0.25_0.005_1120,
            MB_list[[9]]@meta.data$DF.classifications_0.25_0.005_819,
            MB_list[[10]]@meta.data$DF.classifications_0.25_0.17_548,
            MB_list[[11]]@meta.data$DF.classifications_0.25_0.09_639,
            MB_list[[12]]@meta.data$DF.classifications_0.25_0.23_844,
            MB_list[[13]]@meta.data$DF.classifications_0.25_0.26_862,
            MB_list[[14]]@meta.data$DF.classifications_0.25_0.25_637,
            MB_list[[15]]@meta.data$DF.classifications_0.25_0.3_967,
            MB_list[[16]]@meta.data$DF.classifications_0.25_0.29_1186)

# > table(doublet)
# doublet
# Doublet Singlet 
# 10266  132024 

MB.combined$doublet = doublet
save(MB.combined,file = "/bgfs/gkohanbash/SimonsLab/R_Share/processing/2_integrateData_doublet_anno.RData") 

## plot singlet-doublet percentage bar plot
barplot = data.frame(cluster = MB.combined$seurat_clusters, Doublet = MB.combined$doublet)
ggplot(barplot,aes(x = cluster, fill = Doublet)) +
  geom_bar(position = position_fill()) +
  scale_color_brewer(palette = "Set1") +
  theme_classic() +
  labs(y = "Percentage")+
  geom_hline(yintercept = 0.75, linetype= 2)

# removing doublet percentage > 25% clusters (19,26,28) & all doublet
MB.combined = subset(MB.combined,subset = seurat_clusters %in% c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,
                                                                 20,21,22,23,24,25,27,29,30,31))
MB.combined = subset(MB.combined,subset = doublet == "Singlet")

save(MB.combined,file = "/bgfs/gkohanbash/SimonsLab/R_Share/processing/2_integrateData_doublet_removed.RData")

UMAPPlot(MB.combined,group.by = "seurat_clusters",label = T)+NoLegend()
