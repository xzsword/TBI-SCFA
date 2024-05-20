
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
