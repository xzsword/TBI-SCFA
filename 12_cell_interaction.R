#
##
###

# devtools::install_github("sqjin/CellChat")

library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 8000 * 1024^2)

# mark the subtype of each cell
load("/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/4_integrateData_with_group_info.RData")
load("/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/5_Neuron_Seurat_object.RData")
load("/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/6_Astrocyte_Seurat_object.RData")
load("/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/7_Oligodendrocyte_Seurat_object.RData")
load("/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/8_Ependymal_Seurat_object.RData")
load("/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/9_Vascular_Seurat_object.RData")
load("/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/10_Myeloid_Seurat_object.RData")
load("/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/11_TB_Seurat_object.RData")

subtype_n = data.frame(CB = colnames(MB.neuron),
                       subtype = as.character(MB.neuron$subtype))
subtype_a = data.frame(CB = colnames(MB.Astrocyte),
                       subtype = as.character(MB.Astrocyte$subtype))
subtype_o = data.frame(CB = colnames(MB.Oligodendrocyte),
                       subtype = as.character(MB.Oligodendrocyte$subtype))
subtype_e = data.frame(CB = colnames(MB.Ependymal),
                       subtype = as.character(MB.Ependymal$subtype))
subtype_v = data.frame(CB = colnames(MB.Vascular),
                       subtype = as.character(MB.Vascular$subtype))
subtype_m = data.frame(CB = colnames(MB.Myeloid),
                       subtype = as.character(MB.Myeloid$subtype))
subtype_t = data.frame(CB = colnames(MB.TB),
                       subtype = as.character(MB.TB$subtype))
subtype = rbind(subtype_n,subtype_a,subtype_o,subtype_e,subtype_v,subtype_m,subtype_t)
rownames(subtype) = subtype$CB
subtype = subtype[colnames(MB.combined),]
MB.combined$subtype = subtype$subtype

MB.combined = SCTransform(MB.combined,
                          variable.features.n = 5000,
                          ncells = 20000,
                          vars.to.regress = c("percent.mt"))

save(MB.combined,file = "/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/12_integrateData_subtype.RData")

# injury_control overview
MB.injury = MB.combined[,MB.combined$group == "injury_control"]

data.input = GetAssayData(MB.injury,slot = "data") # normalized data matrix
meta = MB.injury@meta.data # a dataframe with rownames containing cell mata data
cell.use = rownames(meta)# extract the cell names from disease data

data.input = as.matrix(data.input[VariableFeatures(MB.injury), cell.use])
meta = meta[cell.use, ]
unique(meta$celltype) # check the cell labels/ types

# create cellchat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")            
# if the object = seurat_object, then the meta will be automatically equals metadata in seurat and it changes
# cellchat <- createCellChat(object = seurat_object)

# # add new cell identity as the default cell type 
# cellchat <- addMeta(cellchat, meta = meta)
# cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
# levels(cellchat@idents) # show factor levels of the cell labels
# groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

# load cellchat database
CellChatDB <- CellChatDB.mouse # use  if running on mouse data
showDatabaseCategory(CellChatDB)

# use a subset of CellChatDB for cell-cell communication analysis
unique(CellChatDB[["interaction"]][["annotation"]])
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

# pre-processing for cellchat analysis
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multisession", workers = 4) # do parallel
#> Warning: [ONE-TIME WARNING] Forked processing ('multicore') is disabled
#> in future (>= 1.13.0) when running R from RStudio, because it is
#> considered unstable. Because of this, plan("multicore") will fall
#> back to plan("sequential"), and plan("multiprocess") will fall back to
#> plan("multisession") - not plan("multicore") as in the past. For more details,
#> how to control forked processing or not, and how to silence this warning in
#> future R sessions, see ?future::supportsMulticore
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

#calculate the cell-cell interaction
cellchat <- computeCommunProb(cellchat, type = "triMean") # if the expression percentage of any cell cluster < 0.25, the L-R pair will be remove
#type = "truncatedMean",trim = 0.1) # to reset the cutoff for the percentage

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

# calculate the interaction in pathway level
cellchat <- computeCommunProbPathway(cellchat)
# net : ligand-receptor pairs
# netP : LR pathways

cellchat@netP$pathways
head(cellchat@LR$LRsig)

# aggregate results for visualization
cellchat <- aggregateNet(cellchat)

# check the total interaction counts among cells and the total interaction strength among cells
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Interaction weights/strength")
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_celltype_interaction_overview.pdf 8*16

# check the cell type sending signal strength
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {  
  mat2 <- matrix(0, nrow = nrow(mat),
                 ncol = ncol(mat), 
                 dimnames = dimnames(mat))  
  mat2[i, ] <- mat[i, ]  
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), # line width = interaction strength
                   arrow.size = 0.1,
                   title.name = rownames(mat)[i])
}
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_celltype_interaction_separately 16*16

# injury_SCFA overview
MB.SCFA = MB.combined[,MB.combined$group == "injury_SCFA"]


data.input = GetAssayData(MB.SCFA,slot = "data") # normalized data matrix
meta = MB.SCFA@meta.data # a dataframe with rownames containing cell mata data
cell.use = rownames(meta)# extract the cell names from disease data

data.input = as.matrix(data.input[VariableFeatures(MB.SCFA), cell.use])
meta = meta[cell.use, ]
unique(meta$celltype) # check the cell labels/ types

# create cellchat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")            
# if the object = seurat_object, then the meta will be automatically equals metadata in seurat and it changes
# cellchat <- createCellChat(object = seurat_object)

# # add new cell identity as the default cell type 
# cellchat <- addMeta(cellchat, meta = meta)
# cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
# levels(cellchat@idents) # show factor levels of the cell labels
# groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

# load cellchat database
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# use a subset of CellChatDB for cell-cell communication analysis
unique(CellChatDB[["interaction"]][["annotation"]])
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

# pre-processing for cellchat analysis
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multisession", workers = 4) # do parallel
#> Warning: [ONE-TIME WARNING] Forked processing ('multicore') is disabled
#> in future (>= 1.13.0) when running R from RStudio, because it is
#> considered unstable. Because of this, plan("multicore") will fall
#> back to plan("sequential"), and plan("multiprocess") will fall back to
#> plan("multisession") - not plan("multicore") as in the past. For more details,
#> how to control forked processing or not, and how to silence this warning in
#> future R sessions, see ?future::supportsMulticore
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

#calculate the cell-cell interaction
cellchat <- computeCommunProb(cellchat, type = "triMean") # if the expression percentage of any cell cluster < 0.25, the L-R pair will be remove
#type = "truncatedMean",trim = 0.1) # to reset the cutoff for the percentage

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

# calculate the interaction in pathway level
cellchat <- computeCommunProbPathway(cellchat)
# net : ligand-receptor pairs
# netP : LR pathways

cellchat@netP$pathways
head(cellchat@LR$LRsig)

# aggregate results for visualization
cellchat <- aggregateNet(cellchat)

# check the total interaction counts among cells and the total interaction strength among cells
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Interaction weights/strength")
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_SCFA_celltype_interaction_overview.pdf 8*16

# check the cell type sending signal strength
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {  
  mat2 <- matrix(0, nrow = nrow(mat),
                 ncol = ncol(mat), 
                 dimnames = dimnames(mat))  
  mat2[i, ] <- mat[i, ]  
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), # line width = interaction strength
                   arrow.size = 0.1,
                   title.name = rownames(mat)[i])
}
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_SCFA_celltype_interaction_separately 16*16

# sham_control overview
MB.sham = MB.combined[,MB.combined$group == "sham_control"]

data.input = GetAssayData(MB.sham,slot = "data") # normalized data matrix
meta = MB.sham@meta.data # a dataframe with rownames containing cell mata data
cell.use = rownames(meta)# extract the cell names from disease data

data.input = as.matrix(data.input[VariableFeatures(MB.sham), cell.use])
meta = meta[cell.use, ]
unique(meta$celltype) # check the cell labels/ types

# create cellchat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")            
# if the object = seurat_object, then the meta will be automatically equals metadata in seurat and it changes
# cellchat <- createCellChat(object = seurat_object)

# # add new cell identity as the default cell type 
# cellchat <- addMeta(cellchat, meta = meta)
# cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
# levels(cellchat@idents) # show factor levels of the cell labels
# groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

# load cellchat database
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# use a subset of CellChatDB for cell-cell communication analysis
unique(CellChatDB[["interaction"]][["annotation"]])
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

# pre-processing for cellchat analysis
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multisession", workers = 4) # do parallel
#> Warning: [ONE-TIME WARNING] Forked processing ('multicore') is disabled
#> in future (>= 1.13.0) when running R from RStudio, because it is
#> considered unstable. Because of this, plan("multicore") will fall
#> back to plan("sequential"), and plan("multiprocess") will fall back to
#> plan("multisession") - not plan("multicore") as in the past. For more details,
#> how to control forked processing or not, and how to silence this warning in
#> future R sessions, see ?future::supportsMulticore
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

#calculate the cell-cell interaction
cellchat <- computeCommunProb(cellchat, type = "triMean") # if the expression percentage of any cell cluster < 0.25, the L-R pair will be remove
#type = "truncatedMean",trim = 0.1) # to reset the cutoff for the percentage

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

# calculate the interaction in pathway level
cellchat <- computeCommunProbPathway(cellchat)
# net : ligand-receptor pairs
# netP : LR pathways

cellchat@netP$pathways
head(cellchat@LR$LRsig)

# aggregate results for visualization
cellchat <- aggregateNet(cellchat)

# check the total interaction counts among cells and the total interaction strength among cells
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Interaction weights/strength")
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_sham_control_celltype_interaction_overview.pdf 8*16

# check the cell type sending signal strength
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {  
  mat2 <- matrix(0, nrow = nrow(mat),
                 ncol = ncol(mat), 
                 dimnames = dimnames(mat))  
  mat2[i, ] <- mat[i, ]  
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), # line width = interaction strength
                   arrow.size = 0.1,
                   title.name = rownames(mat)[i])
}
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_sham_control_celltype_interaction_separately 16*16

## cellchat comparison between injury_control & injury_SCFA

# diffType = setdiff(unique(MB.injury$subtype),unique(MB.SCFA$subtype))
# MB.injury = MB.injury[,MB.injury$subtype != diffType]

cco.injury <- createCellChat(MB.injury@assays$SCT@data, meta = MB.injury@meta.data, group.by = "celltype")
cco.SCFA <- createCellChat(MB.SCFA@assays$SCT@data, meta = MB.SCFA@meta.data, group.by = "celltype")
cco.sham <- createCellChat(MB.sham@assays$SCT@data, meta = MB.sham@meta.data, group.by = "celltype")

# cellchat analysis for injury_control object
cellchat <- cco.injury
cellchat@DB <- CellChatDB.mouse
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = T, population.size = T,type = "truncatedMean",trim = 0.1)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cco.injury <- cellchat

# ellchat analysis for injury_SCFA object
cellchat <- cco.SCFA
cellchat@DB <- CellChatDB.mouse
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = T, population.size = T,type = "truncatedMean",trim = 0.1)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cco.SCFA <- cellchat

# cellchat analysis for injury_sham object
cellchat <- cco.sham
cellchat@DB <- CellChatDB.mouse
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = T, population.size = T,type = "truncatedMean",trim = 0.1)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cco.sham <- cellchat

save(cco.injury, cco.SCFA, cco.sham,file = "/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/12_cellchat_injury_control_SCFA_sham_celltype_compare_input.RData")

# merge cellchat object
cco.list <- list(injury_control = cco.injury, # control group
                 injury_SCFA = cco.SCFA,# test group
                 injury_sham =cco.sham)     
cellchat <- mergeCellChat(cco.list, add.names = names(cco.list), cell.prefix = T)

# compare the total count/strength of the cell interaction
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "weight")
p <- gg1 + gg2
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_Cellchat_injury_control_SCFA_sham_celltype_comparison_Overview_number_strength.pdf", p, width = 6, height = 4)

# merge cellchat object - injury_SCFA - injury_control
cco.list <- list(injury_control = cco.injury, # control group
                 injury_SCFA = cco.SCFA)      
cellchat <- mergeCellChat(cco.list, add.names = names(cco.list), cell.prefix = T)

# compare number and strength of cellchat
par(mfrow = c(1,2))
weight.max <- getMaxWeight(cco.list, attribute = c("idents", "count"))
for(i in 1:length(cco.list)){
  netVisual_circle(cco.list[[i]]@net$count, weight.scale = T, label.edge = T, edge.weight.max = weight.max[2],
                   edge.width.max = 12, title.name = paste0("Number of interactions - ", names(cco.list)[i]))
}
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_celltype_interaction_number_comparison 8*16

weight.max <- getMaxWeight(cco.list, attribute = c("idents", "weight"))
for(i in 1:length(cco.list)){
  netVisual_circle(cco.list[[i]]@net$weight, weight.scale = T, label.edge = F, edge.weight.max = weight.max[2],
                   edge.width.max = 12, title.name = paste0("Strength of interactions - ", names(cco.list)[i]))
}

# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_celltype_interaction_strength_comparison 8*16

par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_celltype_interaction_comparison_difference 8*16
# red: SCFA   blue:injuery_control

# comparison of strength (conservative/specific)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = T)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = T)
p <- gg1 + gg2
p
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_celltype_interaction_Compare_pathway_strength.pdf",p,width = 10,height = 6)

# differential pathway strength
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
# netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
# netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2)
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
p <- rankSimilarity(cellchat, type = "structural") + ggtitle("structural similarity of pathway")
p
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_celltype_interaction_Pathway_Similarity.pdf", p, width = 8, height = 5)

save(cellchat, file = "/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/12_cellchat_injury_control_SCFA_celltype.RData")

# sender and receptor of pathway
library(ComplexHeatmap)

## total
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 <- netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "all", signaling = pathway.union,
                                         title = names(cco.list)[1], width = 8, height = 20)
ht2 <- netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "all", signaling = pathway.union,
                                         title = names(cco.list)[2], width = 8, height = 20)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

## outgoing signaling
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 <- netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "outgoing", signaling = pathway.union,
                                         title = names(cco.list)[1], width = 8, height = 22)
ht2 <- netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "outgoing", signaling = pathway.union,
                                         title = names(cco.list)[2], width = 8, height = 22)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_celltype_outgoing_signal_comparison 12*10

## matrix for outgoing signal
mat.o.c = ht1@matrix
mat.o.c[is.nan(mat.o.c)] <- NA
mat.o.c[is.na(mat.o.c)] <- 0

mat.o.s = ht2@matrix
mat.o.s[is.nan(mat.o.s)] <- NA
mat.o.s[is.na(mat.o.s)] <- 0

## heatmap
mat = mat.o.s - mat.o.c

color.heatmap = "BuGn"
color.use <- scPalette(length(colnames(mat)))
color.heatmap.use = colorRampPalette(colors = c("blue", "white", "red"))(100)

df <- data.frame(group = colnames(mat))
rownames(df) <- colnames(mat)
names(color.use) <- colnames(mat)
col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
                                    which = "column", show_legend = FALSE, show_annotation_name = FALSE, 
                                    simple_anno_size = grid::unit(0.2, "cm"))
ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(mat), 
                                                border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
                        show_annotation_name = FALSE)
pSum <- rowSums(mat)
pSum.original <- pSum
pSum <- -1/log(pSum)
pSum[is.na(pSum)] <- 0
idx1 <- which(is.infinite(pSum) | pSum < 0)
if (length(idx1) > 0) {
  values.assign <- seq(max(pSum) * 1.1, max(pSum) * 1.5, 
                       length.out = length(idx1))
  position <- sort(pSum.original[idx1], index.return = TRUE)$ix
  pSum[idx1] <- values.assign[match(1:length(idx1), position)]
}
ha1 = rowAnnotation(Strength = anno_barplot(pSum, border = FALSE), 
                    show_annotation_name = FALSE)
if (min(mat, na.rm = T) == max(mat, na.rm = T)) {
  legend.break <- max(mat, na.rm = T)
}else {
  legend.break <- c(round(min(mat, na.rm = T), digits = 1), 
                    round(max(mat, na.rm = T), digits = 1))
}

mat[mat == 0] <- NA

Heatmap(mat, col = color.heatmap.use, na_col = "white", 
        name = "Relative strength", bottom_annotation = col_annotation, 
        top_annotation = ha2, right_annotation = ha1, cluster_rows = F, 
        cluster_columns = F, row_names_side = "left", 
        row_names_rot = 0, row_names_gp = gpar(fontsize = 8), 
        column_names_gp = gpar(fontsize = 8), width = unit(10,"cm"), 
        height = unit(20, "cm"), column_title = "Different outgoing signal patterns", 
        column_title_gp = gpar(fontsize = 10), 
        column_names_rot = 90, heatmap_legend_param = list(title_gp = gpar(fontsize = 8,fontface = "plain"), title_position = "leftcenter-rot", 
                                                           border = NA, at = legend.break, legend_height = unit(20,"mm"), labels_gp = gpar(fontsize = 8), grid_width = unit(2,"mm")))
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_celltype_outgoing_signal_DE_strength 12*10
save(mat,file = "/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/12_cellchat_injury_control_SCFA_celltype_outgoing_signal_compare_result.RData")

## incoming signaling
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 <- netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "incoming", signaling = pathway.union,
                                         title = names(cco.list)[1], width = 8, height = 22)
ht2 <- netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "incoming", signaling = pathway.union,
                                         title = names(cco.list)[2], width = 8, height = 22)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_celltype_incoming_signal_comparison 12*10

## matrix for incoming signal
mat.o.c = ht1@matrix
mat.o.c[is.nan(mat.o.c)] <- NA
mat.o.c[is.na(mat.o.c)] <- 0

mat.o.s = ht2@matrix
mat.o.s[is.nan(mat.o.s)] <- NA
mat.o.s[is.na(mat.o.s)] <- 0

## heatmap
mat = mat.o.s - mat.o.c

color.heatmap = "BuGn"
color.use <- scPalette(length(colnames(mat)))
color.heatmap.use = colorRampPalette(colors = c("blue", "white", "red"))(100)

df <- data.frame(group = colnames(mat))
rownames(df) <- colnames(mat)
names(color.use) <- colnames(mat)
col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
                                    which = "column", show_legend = FALSE, show_annotation_name = FALSE, 
                                    simple_anno_size = grid::unit(0.2, "cm"))
ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(mat), 
                                                border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
                        show_annotation_name = FALSE)
pSum <- rowSums(mat)
pSum.original <- pSum
pSum <- -1/log(pSum)
pSum[is.na(pSum)] <- 0
idx1 <- which(is.infinite(pSum) | pSum < 0)
if (length(idx1) > 0) {
  values.assign <- seq(max(pSum) * 1.1, max(pSum) * 1.5, 
                       length.out = length(idx1))
  position <- sort(pSum.original[idx1], index.return = TRUE)$ix
  pSum[idx1] <- values.assign[match(1:length(idx1), position)]
}
ha1 = rowAnnotation(Strength = anno_barplot(pSum, border = FALSE), 
                    show_annotation_name = FALSE)
if (min(mat, na.rm = T) == max(mat, na.rm = T)) {
  legend.break <- max(mat, na.rm = T)
}else {
  legend.break <- c(round(min(mat, na.rm = T), digits = 1), 
                    round(max(mat, na.rm = T), digits = 1))
}

mat[mat == 0] <- NA

Heatmap(mat, col = color.heatmap.use, na_col = "white", 
        name = "Relative strength", bottom_annotation = col_annotation, 
        top_annotation = ha2, right_annotation = ha1, cluster_rows = F, 
        cluster_columns = F, row_names_side = "left", 
        row_names_rot = 0, row_names_gp = gpar(fontsize = 8), 
        column_names_gp = gpar(fontsize = 8), width = unit(10,"cm"), 
        height = unit(20, "cm"), column_title = "Different incoming signal patterns", 
        column_title_gp = gpar(fontsize = 10), 
        column_names_rot = 90, heatmap_legend_param = list(title_gp = gpar(fontsize = 8,fontface = "plain"), title_position = "leftcenter-rot", 
                                                           border = NA, at = legend.break, legend_height = unit(20,"mm"), labels_gp = gpar(fontsize = 8), grid_width = unit(2,"mm")))
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_celltype_incoming_signal_DE_strength 12*10
save(mat,file = "/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/12_cellchat_injury_control_SCFA_celltype_incoming_signal_compare_result.RData")

## show up / down L-R pairs in individual celltypes
cell = sort(unique(as.character(cco.list[[1]]@idents)))

for(i in 1:length(cell)){
p1 <- netVisual_bubble(cellchat, sources.use = c(i), targets.use = c(1:12), comparison = c(1,2),
                       max.dataset = 2, title.name = "Increased signaling in injury_SCFA", angle.x = 45, remove.isolate = T)
p2 <- netVisual_bubble(cellchat, sources.use = c(i), targets.use = c(1:12), comparison = c(1,2),
                       max.dataset = 1, title.name = "Decreased signaling in injury_SCFA", angle.x = 45, remove.isolate = T)
pc <- p1 + p2
pc
ggsave(paste0("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_",cell[i],"_diff_LR_pairs.pdf"), pc, width = 16, height = 10)
}

# comparison for specific signaling
## MK pathway: MK-deficiency reduced tissue infiltration of microglia/macrophages 
##             and altered their polarization status thereby reducing neuroinflammation, neuronal apoptosis, 
##             and tissue loss and improving neurological outcomes after TBI.
## citation: Takada, S., Sakakima, H., Matsuyama, T. et al. Disruption of Midkine gene reduces traumatic brain injury through the modulation of neuroinflammation. J Neuroinflammation 17, 40 (2020). https://doi.org/10.1186/s12974-020-1709-8
pathways.show <- c("MK")
weight.max <- getMaxWeight(cco.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd = TRUE)
for(i in 1:length(cco.list)){
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(cco.list)[i]))
}

# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_celltype_MK_signaling_pathway_comparison 6*12

# explor MK pathway in subtype
## astrocyte - myeloid
MB.s.injury = MB.injury[,MB.injury$celltype %in% c("Astrocyte","Myeloid_cell")]
MB.s.SCFA = MB.SCFA[,MB.SCFA$celltype %in% c("Astrocyte","Myeloid_cell")]

## generate cellchat object for subtype
cco.injury <- createCellChat(MB.s.injury@assays$SCT@data, meta = MB.s.injury@meta.data, group.by = "subtype")
cco.SCFA <- createCellChat(MB.s.SCFA@assays$SCT@data, meta = MB.s.SCFA@meta.data, group.by = "subtype")

# cellchat analysis for injury_control object
s.cellchat <- cco.injury
s.cellchat@DB <- CellChatDB.mouse
s.cellchat <- subsetData(s.cellchat)
s.cellchat <- identifyOverExpressedGenes(s.cellchat)
s.cellchat <- identifyOverExpressedInteractions(s.cellchat)
s.cellchat <- computeCommunProb(s.cellchat, raw.use = T, population.size = T)
s.cellchat <- computeCommunProbPathway(s.cellchat)
s.cellchat <- aggregateNet(s.cellchat)
s.cellchat <- netAnalysis_computeCentrality(s.cellchat, slot.name = "netP")
cco.injury <- s.cellchat

# ellchat analysis for injury_SCFA object
s.cellchat <- cco.SCFA
s.cellchat@DB <- CellChatDB.mouse
s.cellchat <- subsetData(s.cellchat)
s.cellchat <- identifyOverExpressedGenes(s.cellchat)
s.cellchat <- identifyOverExpressedInteractions(s.cellchat)
s.cellchat <- computeCommunProb(s.cellchat, raw.use = T, population.size = T)
s.cellchat <- computeCommunProbPathway(s.cellchat)
s.cellchat <- aggregateNet(s.cellchat)
s.cellchat <- netAnalysis_computeCentrality(s.cellchat, slot.name = "netP")
cco.SCFA <- s.cellchat

s.cco.list <- list(injury_control = cco.injury, # control group
                 injury_SCFA = cco.SCFA)      # test group
s.cellchat <- mergeCellChat(s.cco.list, add.names = names(s.cco.list), cell.prefix = T)

# visualize MK pathway for Astrocyte - Myeloid cell subtype
pathways.show <- c("MK")
weight.max <- getMaxWeight(s.cco.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd = TRUE)
for(i in 1:length(cco.list)){
  netVisual_aggregate(s.cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(s.cco.list)[i]))
}
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_Astrocyte_Myeloid_subtype_MK_signaling_pathway_comparison 6*12

## astrocyte - neuron
MB.s.injury = MB.injury[,MB.injury$celltype %in% c("Myeloid_cell","OPC")]
# MB.s.injury = MB.s.injury[,MB.s.injury$subtype != "GLUT_2"]
MB.s.SCFA = MB.SCFA[,MB.SCFA$celltype %in% c("Myeloid_cell","OPC")]

## generate cellchat object for subtype
cco.injury <- createCellChat(MB.s.injury@assays$SCT@data, meta = MB.s.injury@meta.data, group.by = "subtype")
cco.SCFA <- createCellChat(MB.s.SCFA@assays$SCT@data, meta = MB.s.SCFA@meta.data, group.by = "subtype")

# cellchat analysis for injury_control object
s.cellchat <- cco.injury
s.cellchat@DB <- CellChatDB.mouse
s.cellchat <- subsetData(s.cellchat)
s.cellchat <- identifyOverExpressedGenes(s.cellchat)
s.cellchat <- identifyOverExpressedInteractions(s.cellchat)
s.cellchat <- computeCommunProb(s.cellchat, raw.use = T, population.size = T)
s.cellchat <- computeCommunProbPathway(s.cellchat)
s.cellchat <- aggregateNet(s.cellchat)
s.cellchat <- netAnalysis_computeCentrality(s.cellchat, slot.name = "netP")
cco.injury <- s.cellchat

# ellchat analysis for injury_SCFA object
s.cellchat <- cco.SCFA
s.cellchat@DB <- CellChatDB.mouse
s.cellchat <- subsetData(s.cellchat)
s.cellchat <- identifyOverExpressedGenes(s.cellchat)
s.cellchat <- identifyOverExpressedInteractions(s.cellchat)
s.cellchat <- computeCommunProb(s.cellchat, raw.use = T, population.size = T)
s.cellchat <- computeCommunProbPathway(s.cellchat)
s.cellchat <- aggregateNet(s.cellchat)
s.cellchat <- netAnalysis_computeCentrality(s.cellchat, slot.name = "netP")
cco.SCFA <- s.cellchat

s.cco.list <- list(injury_control = cco.injury, # control group
                   injury_SCFA = cco.SCFA)      # test group
s.cellchat <- mergeCellChat(s.cco.list, add.names = names(s.cco.list), cell.prefix = T)

# visualize MK pathway for Astrocyte - Myeloid cell subtype
pathways.show <- c("MK")
weight.max <- getMaxWeight(s.cco.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd = TRUE)
for(i in 1:length(cco.list)){
  netVisual_aggregate(s.cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(s.cco.list)[i]))
}
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_OPC_Myeloid_subtype_MK_signaling_pathway_comparison 6*12

## ICAM signaling pathway
pathways.show <- c("ICAM")
weight.max <- getMaxWeight(cco.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd = TRUE)
for(i in 1:length(cco.list)){
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(cco.list)[i]))
}

# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_celltype_ICAM_signaling_pathway_comparison 6*12
## no diff
## myeloid - endothelia
MB.s.injury = MB.injury[,MB.injury$celltype %in% c("Myeloid_cell","Endothelial_cell")]
MB.s.injury = MB.s.injury[,-grep("PC",MB.s.injury$subtype)]
MB.s.injury = MB.s.injury[,-grep("VSMC",MB.s.injury$subtype)]

# MB.s.injury = MB.s.injury[,MB.s.injury$subtype != "GLUT_2"]
MB.s.SCFA = MB.SCFA[,MB.SCFA$celltype %in% c("Myeloid_cell","Endothelial_cell")]
MB.s.SCFA = MB.s.SCFA[,-grep("PC",MB.s.SCFA$subtype)]
MB.s.SCFA = MB.s.SCFA[,-grep("VSMC",MB.s.SCFA$subtype)]

## generate cellchat object for subtype
cco.injury <- createCellChat(MB.s.injury@assays$SCT@data, meta = MB.s.injury@meta.data, group.by = "subtype")
cco.SCFA <- createCellChat(MB.s.SCFA@assays$SCT@data, meta = MB.s.SCFA@meta.data, group.by = "subtype")

# cellchat analysis for injury_control object
s.cellchat <- cco.injury
s.cellchat@DB <- CellChatDB.mouse
s.cellchat <- subsetData(s.cellchat)
s.cellchat <- identifyOverExpressedGenes(s.cellchat)
s.cellchat <- identifyOverExpressedInteractions(s.cellchat)
s.cellchat <- computeCommunProb(s.cellchat, raw.use = T, population.size = T)
s.cellchat <- computeCommunProbPathway(s.cellchat)
s.cellchat <- aggregateNet(s.cellchat)
s.cellchat <- netAnalysis_computeCentrality(s.cellchat, slot.name = "netP")
cco.injury <- s.cellchat

# ellchat analysis for injury_SCFA object
s.cellchat <- cco.SCFA
s.cellchat@DB <- CellChatDB.mouse
s.cellchat <- subsetData(s.cellchat)
s.cellchat <- identifyOverExpressedGenes(s.cellchat)
s.cellchat <- identifyOverExpressedInteractions(s.cellchat)
s.cellchat <- computeCommunProb(s.cellchat, raw.use = T, population.size = T)
s.cellchat <- computeCommunProbPathway(s.cellchat)
s.cellchat <- aggregateNet(s.cellchat)
s.cellchat <- netAnalysis_computeCentrality(s.cellchat, slot.name = "netP")
cco.SCFA <- s.cellchat

s.cco.list <- list(injury_control = cco.injury, # control group
                   injury_SCFA = cco.SCFA)      # test group
s.cellchat <- mergeCellChat(s.cco.list, add.names = names(s.cco.list), cell.prefix = T)

# visualize ICAM pathway for Astrocyte - Myeloid cell subtype
pathways.show <- c("ICAM")
weight.max <- getMaxWeight(s.cco.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd = TRUE)
for(i in 1:length(cco.list)){
  netVisual_aggregate(s.cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(s.cco.list)[i]))
}
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_Endothelial_Myeloid_subtype_ICAM_signaling_pathway_comparison 7*14

netAnalysis_contribution(s.cco.list[[1]], signaling = pathways.show)
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_ICAM_signaling_pathway_contribution 4*6

## VCAM signaling pathway
pathways.show <- c("VCAM")
weight.max <- getMaxWeight(cco.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd = TRUE)
for(i in 1:length(cco.list)){
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(cco.list)[i]))
}

# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_celltype_VCAM_signaling_pathway_comparison 6*12
## no diff
## myeloid - endothelia
MB.s.injury = MB.injury[,MB.injury$celltype %in% c("Myeloid_cell","Endothelial_cell")]
MB.s.injury = MB.s.injury[,-grep("PC",MB.s.injury$subtype)]
MB.s.injury = MB.s.injury[,-grep("VSMC",MB.s.injury$subtype)]

# MB.s.injury = MB.s.injury[,MB.s.injury$subtype != "GLUT_2"]
MB.s.SCFA = MB.SCFA[,MB.SCFA$celltype %in% c("Myeloid_cell","Endothelial_cell")]
MB.s.SCFA = MB.s.SCFA[,-grep("PC",MB.s.SCFA$subtype)]
MB.s.SCFA = MB.s.SCFA[,-grep("VSMC",MB.s.SCFA$subtype)]

## generate cellchat object for subtype
cco.injury <- createCellChat(MB.s.injury@assays$SCT@data, meta = MB.s.injury@meta.data, group.by = "subtype")
cco.SCFA <- createCellChat(MB.s.SCFA@assays$SCT@data, meta = MB.s.SCFA@meta.data, group.by = "subtype")

# cellchat analysis for injury_control object
s.cellchat <- cco.injury
s.cellchat@DB <- CellChatDB.mouse
s.cellchat <- subsetData(s.cellchat)
s.cellchat <- identifyOverExpressedGenes(s.cellchat)
s.cellchat <- identifyOverExpressedInteractions(s.cellchat)
s.cellchat <- computeCommunProb(s.cellchat, raw.use = T, population.size = T)
s.cellchat <- computeCommunProbPathway(s.cellchat)
s.cellchat <- aggregateNet(s.cellchat)
s.cellchat <- netAnalysis_computeCentrality(s.cellchat, slot.name = "netP")
cco.injury <- s.cellchat

# ellchat analysis for injury_SCFA object
s.cellchat <- cco.SCFA
s.cellchat@DB <- CellChatDB.mouse
s.cellchat <- subsetData(s.cellchat)
s.cellchat <- identifyOverExpressedGenes(s.cellchat)
s.cellchat <- identifyOverExpressedInteractions(s.cellchat)
s.cellchat <- computeCommunProb(s.cellchat, raw.use = T, population.size = T)
s.cellchat <- computeCommunProbPathway(s.cellchat)
s.cellchat <- aggregateNet(s.cellchat)
s.cellchat <- netAnalysis_computeCentrality(s.cellchat, slot.name = "netP")
cco.SCFA <- s.cellchat

s.cco.list <- list(injury_control = cco.injury, # control group
                   injury_SCFA = cco.SCFA)      # test group
s.cellchat <- mergeCellChat(s.cco.list, add.names = names(s.cco.list), cell.prefix = T)

# visualize VCAM pathway for Astrocyte - Myeloid cell subtype
pathways.show <- c("VCAM")
weight.max <- getMaxWeight(s.cco.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd = TRUE)
for(i in 1:length(cco.list)){
  netVisual_aggregate(s.cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(s.cco.list)[i]))
}
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_Endothelial_Myeloid_subtype_VCAM_signaling_pathway_comparison 7*14

netAnalysis_contribution(s.cco.list[[1]], signaling = pathways.show)
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_VCAM_signaling_pathway_contribution 4*6

## VEGF signaling pathway
pathways.show <- c("VEGF")
weight.max <- getMaxWeight(cco.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd = TRUE)
for(i in 1:length(cco.list)){
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(cco.list)[i]))
}

# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_celltype_VEGF_signaling_pathway_comparison 6*12

## astrocyte - endothelia
MB.s.injury = MB.injury[,MB.injury$celltype %in% c("Astrocyte","Endothelial_cell")]
MB.s.injury = MB.s.injury[,-grep("PC",MB.s.injury$subtype)]
MB.s.injury = MB.s.injury[,-grep("VSMC",MB.s.injury$subtype)]

# MB.s.injury = MB.s.injury[,MB.s.injury$subtype != "GLUT_2"]
MB.s.SCFA = MB.SCFA[,MB.SCFA$celltype %in% c("Astrocyte","Endothelial_cell")]
MB.s.SCFA = MB.s.SCFA[,-grep("PC",MB.s.SCFA$subtype)]
MB.s.SCFA = MB.s.SCFA[,-grep("VSMC",MB.s.SCFA$subtype)]

## generate cellchat object for subtype
cco.injury <- createCellChat(MB.s.injury@assays$SCT@data, meta = MB.s.injury@meta.data, group.by = "subtype")
cco.SCFA <- createCellChat(MB.s.SCFA@assays$SCT@data, meta = MB.s.SCFA@meta.data, group.by = "subtype")

# cellchat analysis for injury_control object
s.cellchat <- cco.injury
s.cellchat@DB <- CellChatDB.mouse
s.cellchat <- subsetData(s.cellchat)
s.cellchat <- identifyOverExpressedGenes(s.cellchat)
s.cellchat <- identifyOverExpressedInteractions(s.cellchat)
s.cellchat <- computeCommunProb(s.cellchat, raw.use = T, population.size = T)
s.cellchat <- computeCommunProbPathway(s.cellchat)
s.cellchat <- aggregateNet(s.cellchat)
s.cellchat <- netAnalysis_computeCentrality(s.cellchat, slot.name = "netP")
cco.injury <- s.cellchat

# ellchat analysis for injury_SCFA object
s.cellchat <- cco.SCFA
s.cellchat@DB <- CellChatDB.mouse
s.cellchat <- subsetData(s.cellchat)
s.cellchat <- identifyOverExpressedGenes(s.cellchat)
s.cellchat <- identifyOverExpressedInteractions(s.cellchat)
s.cellchat <- computeCommunProb(s.cellchat, raw.use = T, population.size = T)
s.cellchat <- computeCommunProbPathway(s.cellchat)
s.cellchat <- aggregateNet(s.cellchat)
s.cellchat <- netAnalysis_computeCentrality(s.cellchat, slot.name = "netP")
cco.SCFA <- s.cellchat

s.cco.list <- list(injury_control = cco.injury, # control group
                   injury_SCFA = cco.SCFA)      # test group
s.cellchat <- mergeCellChat(s.cco.list, add.names = names(s.cco.list), cell.prefix = T)

# visualize VEGF pathway for Astrocyte - Myeloid cell subtype
pathways.show <- c("VEGF")
weight.max <- getMaxWeight(s.cco.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd = TRUE)
for(i in 1:length(cco.list)){
  netVisual_aggregate(s.cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(s.cco.list)[i]))
}
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_Endothelial_Astrocyte_subtype_VEGF_signaling_pathway_comparison 7*14

## MHC-I signaling pathway
pathways.show <- c("MHC-I")
weight.max <- getMaxWeight(cco.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd = TRUE)
for(i in 1:length(cco.list)){
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(cco.list)[i]))
}

# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_celltype_MHC-I_signaling_pathway_comparison 6*12

## Myeloid - T
MB.s.injury = MB.injury[,MB.injury$celltype %in% c("Myeloid_cell","T_cell")]
# MB.s.injury = MB.s.injury[,MB.s.injury$subtype != "GLUT_2"]
MB.s.SCFA = MB.SCFA[,MB.SCFA$celltype %in% c("Myeloid_cell","T_cell")]
MB.s.SCFA = MB.s.SCFA[,MB.s.SCFA$subtype != "Mature_B"]

## generate cellchat object for subtype
cco.injury <- createCellChat(MB.s.injury@assays$SCT@data, meta = MB.s.injury@meta.data, group.by = "subtype")
cco.SCFA <- createCellChat(MB.s.SCFA@assays$SCT@data, meta = MB.s.SCFA@meta.data, group.by = "subtype")

# cellchat analysis for injury_control object
s.cellchat <- cco.injury
s.cellchat@DB <- CellChatDB.mouse
s.cellchat <- subsetData(s.cellchat)
s.cellchat <- identifyOverExpressedGenes(s.cellchat)
s.cellchat <- identifyOverExpressedInteractions(s.cellchat)
s.cellchat <- computeCommunProb(s.cellchat, raw.use = T, population.size = T)
s.cellchat <- computeCommunProbPathway(s.cellchat)
s.cellchat <- aggregateNet(s.cellchat)
s.cellchat <- netAnalysis_computeCentrality(s.cellchat, slot.name = "netP")
cco.injury <- s.cellchat

# ellchat analysis for injury_SCFA object
s.cellchat <- cco.SCFA
s.cellchat@DB <- CellChatDB.mouse
s.cellchat <- subsetData(s.cellchat)
s.cellchat <- identifyOverExpressedGenes(s.cellchat)
s.cellchat <- identifyOverExpressedInteractions(s.cellchat)
s.cellchat <- computeCommunProb(s.cellchat, raw.use = T, population.size = T)
s.cellchat <- computeCommunProbPathway(s.cellchat)
s.cellchat <- aggregateNet(s.cellchat)
s.cellchat <- netAnalysis_computeCentrality(s.cellchat, slot.name = "netP")
cco.SCFA <- s.cellchat

s.cco.list <- list(injury_control = cco.injury, # control group
                   injury_SCFA = cco.SCFA)      # test group
s.cellchat <- mergeCellChat(s.cco.list, add.names = names(s.cco.list), cell.prefix = T)

# visualize MHC-I pathway for Astrocyte - Myeloid cell subtype
pathways.show <- c("MHC-I")
weight.max <- getMaxWeight(s.cco.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd = TRUE)
for(i in 1:length(cco.list)){
  netVisual_aggregate(s.cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(s.cco.list)[i]))
}
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_Myeloid_T_subtype_MHC-I_signaling_pathway_comparison 7*14

## CCL signaling pathway
pathways.show <- c("CCL")
weight.max <- getMaxWeight(cco.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd = TRUE)
for(i in 1:length(cco.list)){
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(cco.list)[i]))
}

# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_celltype_CCL_signaling_pathway_comparison 6*12

## Myeloid_cell - T_cell
MB.s.injury = MB.injury[,MB.injury$celltype %in% c("Myeloid_cell","T_cell")]
# MB.s.injury = MB.s.injury[,MB.s.injury$subtype != "GLUT_2"]
MB.s.SCFA = MB.SCFA[,MB.SCFA$celltype %in% c("Myeloid_cell","T_cell")]
MB.s.SCFA = MB.s.SCFA[,MB.s.SCFA$subtype != "Mature_B"]

## generate cellchat object for subtype
cco.injury <- createCellChat(MB.s.injury@assays$SCT@data, meta = MB.s.injury@meta.data, group.by = "subtype")
cco.SCFA <- createCellChat(MB.s.SCFA@assays$SCT@data, meta = MB.s.SCFA@meta.data, group.by = "subtype")

# cellchat analysis for injury_control object
s.cellchat <- cco.injury
s.cellchat@DB <- CellChatDB.mouse
s.cellchat <- subsetData(s.cellchat)
s.cellchat <- identifyOverExpressedGenes(s.cellchat)
s.cellchat <- identifyOverExpressedInteractions(s.cellchat)
s.cellchat <- computeCommunProb(s.cellchat, raw.use = T, population.size = T)
s.cellchat <- computeCommunProbPathway(s.cellchat)
s.cellchat <- aggregateNet(s.cellchat)
s.cellchat <- netAnalysis_computeCentrality(s.cellchat, slot.name = "netP")
cco.injury <- s.cellchat

# ellchat analysis for injury_SCFA object
s.cellchat <- cco.SCFA
s.cellchat@DB <- CellChatDB.mouse
s.cellchat <- subsetData(s.cellchat)
s.cellchat <- identifyOverExpressedGenes(s.cellchat)
s.cellchat <- identifyOverExpressedInteractions(s.cellchat)
s.cellchat <- computeCommunProb(s.cellchat, raw.use = T, population.size = T)
s.cellchat <- computeCommunProbPathway(s.cellchat)
s.cellchat <- aggregateNet(s.cellchat)
s.cellchat <- netAnalysis_computeCentrality(s.cellchat, slot.name = "netP")
cco.SCFA <- s.cellchat

s.cco.list <- list(injury_control = cco.injury, # control group
                   injury_SCFA = cco.SCFA)      # test group
s.cellchat <- mergeCellChat(s.cco.list, add.names = names(s.cco.list), cell.prefix = T)

# visualize CCL pathway for Astrocyte - Myeloid cell subtype
pathways.show <- c("CCL")
weight.max <- getMaxWeight(s.cco.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd = TRUE)
for(i in 1:length(cco.list)){
  netVisual_aggregate(s.cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(s.cco.list)[i]))
}
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_Myeloid_T_subtype_CCL_signaling_pathway_comparison 7*14

# 计算injury_control组中配体受体对选定信号通路的贡献值
netAnalysis_contribution(s.cco.list[[1]], signaling = pathways.show)
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_CCL_signaling_pathway_contribution 4*6

## TGFb signaling pathway - promote myeloid cell from M1 to M0/M2 phenotype
## cite: Divolis G, Stavropoulos A, Manioudaki M, Apostolidou A, Doulou A, Gavriil A, Dafnis I, Chroni A, Mummery C, Xilouri M, Sideras P. Activation of both transforming growth factor-β and bone morphogenetic protein signalling pathways upon traumatic brain injury restrains pro-inflammatory and boosts tissue reparatory responses of reactive astrocytes and microglia. Brain Commun. 2019 Oct 21;1(1):fcz028. doi: 10.1093/braincomms/fcz028. PMID: 32954268; PMCID: PMC7425383.
pathways.show <- c("TGFb")
weight.max <- getMaxWeight(cco.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd = TRUE)
for(i in 1:length(cco.list)){
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(cco.list)[i]))
}

# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_celltype_TGFb_signaling_pathway_comparison 6*12

netAnalysis_contribution(s.cco.list[[1]], signaling = pathways.show)
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_TGFb_signaling_pathway_contribution 4*6

## Myeloid_cell - Choroid_plexus_epithelial_cell
MB.s.injury = MB.injury[,MB.injury$celltype %in% c("Myeloid_cell","Choroid_plexus_epithelial_cell")]
MB.s.injury = MB.s.injury[,-grep("EPC",MB.s.injury$subtype)]

# MB.s.injury = MB.s.injury[,MB.s.injury$subtype != "GLUT_2"]
MB.s.SCFA = MB.SCFA[,MB.SCFA$celltype %in% c("Myeloid_cell","Choroid_plexus_epithelial_cell")]
MB.s.SCFA = MB.s.SCFA[,-grep("EPC",MB.s.SCFA$subtype)]

## generate cellchat object for subtype
cco.injury <- createCellChat(MB.s.injury@assays$SCT@data, meta = MB.s.injury@meta.data, group.by = "subtype")
cco.SCFA <- createCellChat(MB.s.SCFA@assays$SCT@data, meta = MB.s.SCFA@meta.data, group.by = "subtype")

# cellchat analysis for injury_control object
s.cellchat <- cco.injury
s.cellchat@DB <- CellChatDB.mouse
s.cellchat <- subsetData(s.cellchat)
s.cellchat <- identifyOverExpressedGenes(s.cellchat)
s.cellchat <- identifyOverExpressedInteractions(s.cellchat)
s.cellchat <- computeCommunProb(s.cellchat, raw.use = T, population.size = T)
s.cellchat <- computeCommunProbPathway(s.cellchat)
s.cellchat <- aggregateNet(s.cellchat)
s.cellchat <- netAnalysis_computeCentrality(s.cellchat, slot.name = "netP")
cco.injury <- s.cellchat

# ellchat analysis for injury_SCFA object
s.cellchat <- cco.SCFA
s.cellchat@DB <- CellChatDB.mouse
s.cellchat <- subsetData(s.cellchat)
s.cellchat <- identifyOverExpressedGenes(s.cellchat)
s.cellchat <- identifyOverExpressedInteractions(s.cellchat)
s.cellchat <- computeCommunProb(s.cellchat, raw.use = T, population.size = T)
s.cellchat <- computeCommunProbPathway(s.cellchat)
s.cellchat <- aggregateNet(s.cellchat)
s.cellchat <- netAnalysis_computeCentrality(s.cellchat, slot.name = "netP")
cco.SCFA <- s.cellchat

s.cco.list <- list(injury_control = cco.injury, # control group
                   injury_SCFA = cco.SCFA)      # test group
s.cellchat <- mergeCellChat(s.cco.list, add.names = names(s.cco.list), cell.prefix = T)

# visualize TGFb pathway for Astrocyte - Myeloid cell subtype
pathways.show <- c("TGFb")
weight.max <- getMaxWeight(s.cco.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd = TRUE)
for(i in 1:length(cco.list)){
  netVisual_aggregate(s.cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(s.cco.list)[i]))
}
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_Myeloid_CPC_subtype_TGFb_signaling_pathway_comparison 7*14

## PECAM1 signaling pathway - maintain the BBB integrity after injury
pathways.show <- c("PECAM1")
weight.max <- getMaxWeight(cco.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd = TRUE)
for(i in 1:length(cco.list)){
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(cco.list)[i]))
}

# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_celltype_PECAM1_signaling_pathway_comparison 6*12

netAnalysis_contribution(s.cco.list[[1]], signaling = pathways.show)
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_PECAM1_signaling_pathway_contribution 4*6

## Endothelial_cell - Pericyte
MB.s.injury = MB.injury[,MB.injury$celltype %in% c("Endothelial_cell","Pericyte")]
MB.s.injury = MB.s.injury[,-grep("VSMC",MB.s.injury$subtype)]

# MB.s.injury = MB.s.injury[,MB.s.injury$subtype != "GLUT_2"]
MB.s.SCFA = MB.SCFA[,MB.SCFA$celltype %in% c("Endothelial_cell","Pericyte")]
MB.s.SCFA = MB.s.SCFA[,-grep("VSMC",MB.s.SCFA$subtype)]

## generate cellchat object for subtype
cco.injury <- createCellChat(MB.s.injury@assays$SCT@data, meta = MB.s.injury@meta.data, group.by = "subtype")
cco.SCFA <- createCellChat(MB.s.SCFA@assays$SCT@data, meta = MB.s.SCFA@meta.data, group.by = "subtype")

# cellchat analysis for injury_control object
s.cellchat <- cco.injury
s.cellchat@DB <- CellChatDB.mouse
s.cellchat <- subsetData(s.cellchat)
s.cellchat <- identifyOverExpressedGenes(s.cellchat)
s.cellchat <- identifyOverExpressedInteractions(s.cellchat)
s.cellchat <- computeCommunProb(s.cellchat, raw.use = T, population.size = T)
s.cellchat <- computeCommunProbPathway(s.cellchat)
s.cellchat <- aggregateNet(s.cellchat)
s.cellchat <- netAnalysis_computeCentrality(s.cellchat, slot.name = "netP")
cco.injury <- s.cellchat

# ellchat analysis for injury_SCFA object
s.cellchat <- cco.SCFA
s.cellchat@DB <- CellChatDB.mouse
s.cellchat <- subsetData(s.cellchat)
s.cellchat <- identifyOverExpressedGenes(s.cellchat)
s.cellchat <- identifyOverExpressedInteractions(s.cellchat)
s.cellchat <- computeCommunProb(s.cellchat, raw.use = T, population.size = T)
s.cellchat <- computeCommunProbPathway(s.cellchat)
s.cellchat <- aggregateNet(s.cellchat)
s.cellchat <- netAnalysis_computeCentrality(s.cellchat, slot.name = "netP")
cco.SCFA <- s.cellchat

s.cco.list <- list(injury_control = cco.injury, # control group
                   injury_SCFA = cco.SCFA)      # test group
s.cellchat <- mergeCellChat(s.cco.list, add.names = names(s.cco.list), cell.prefix = T)

# visualize PECAM1 pathway for Astrocyte - Myeloid cell subtype
weight.max <- getMaxWeight(s.cco.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd = TRUE)
for(i in 1:length(cco.list)){
  netVisual_aggregate(s.cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(s.cco.list)[i]))
}
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_vascular_lineage_subtype_PECAM1_signaling_pathway_comparison 7*14

## CADM signaling pathway - increase the inflammatory cells interaction and may promote inflammation
# cite: Sona C, Yeh YT, Patsalos A, Halasz L, Yan X, Kononenko NL, Nagy L, Poy MN. Evidence of islet CADM1-mediated immune cell interactions during human type 1 diabetes. JCI Insight. 2022 Mar 22;7(6):e153136. doi: 10.1172/jci.insight.153136. PMID: 35133983; PMCID: PMC8986082.
pathways.show <- c("CADM")
weight.max <- getMaxWeight(cco.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd = TRUE)
for(i in 1:length(cco.list)){
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(cco.list)[i]))
}

# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_celltype_CADM_signaling_pathway_comparison 6*12

netAnalysis_contribution(s.cco.list[[1]], signaling = pathways.show)
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_CADM_signaling_pathway_contribution 4*6

## Myeloid_cell - Astrocyte
MB.s.injury = MB.injury[,MB.injury$celltype %in% c("Myeloid_cell","Astrocyte")]
# MB.s.injury = MB.s.injury[,MB.s.injury$subtype != "GLUT_2"]
MB.s.SCFA = MB.SCFA[,MB.SCFA$celltype %in% c("Myeloid_cell","Astrocyte")]

## generate cellchat object for subtype
cco.injury <- createCellChat(MB.s.injury@assays$SCT@data, meta = MB.s.injury@meta.data, group.by = "subtype")
cco.SCFA <- createCellChat(MB.s.SCFA@assays$SCT@data, meta = MB.s.SCFA@meta.data, group.by = "subtype")

# cellchat analysis for injury_control object
s.cellchat <- cco.injury
s.cellchat@DB <- CellChatDB.mouse
s.cellchat <- subsetData(s.cellchat)
s.cellchat <- identifyOverExpressedGenes(s.cellchat)
s.cellchat <- identifyOverExpressedInteractions(s.cellchat)
s.cellchat <- computeCommunProb(s.cellchat, raw.use = T, population.size = T)
s.cellchat <- computeCommunProbPathway(s.cellchat)
s.cellchat <- aggregateNet(s.cellchat)
s.cellchat <- netAnalysis_computeCentrality(s.cellchat, slot.name = "netP")
cco.injury <- s.cellchat

# ellchat analysis for injury_SCFA object
s.cellchat <- cco.SCFA
s.cellchat@DB <- CellChatDB.mouse
s.cellchat <- subsetData(s.cellchat)
s.cellchat <- identifyOverExpressedGenes(s.cellchat)
s.cellchat <- identifyOverExpressedInteractions(s.cellchat)
s.cellchat <- computeCommunProb(s.cellchat, raw.use = T, population.size = T)
s.cellchat <- computeCommunProbPathway(s.cellchat)
s.cellchat <- aggregateNet(s.cellchat)
s.cellchat <- netAnalysis_computeCentrality(s.cellchat, slot.name = "netP")
cco.SCFA <- s.cellchat

s.cco.list <- list(injury_control = cco.injury, # control group
                   injury_SCFA = cco.SCFA)      # test group
s.cellchat <- mergeCellChat(s.cco.list, add.names = names(s.cco.list), cell.prefix = T)

# visualize PECAM1 pathway for Astrocyte - Myeloid cell subtype
weight.max <- getMaxWeight(s.cco.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd = TRUE)
for(i in 1:length(cco.list)){
  netVisual_aggregate(s.cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(s.cco.list)[i]))
}
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_Myeloid_Astrocyte_subtype_CADM_signaling_pathway_comparison 7*14

## CDH5 signaling pathway - microvessel formation and vessel integrity
pathways.show <- c("CDH5")
weight.max <- getMaxWeight(cco.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd = TRUE)
for(i in 1:length(cco.list)){
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(cco.list)[i]))
}

# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_celltype_CDH5_signaling_pathway_comparison 6*12

## CDH2 signaling pathway - scar forming astrocyte marker
pathways.show <- c("CDH")
weight.max <- getMaxWeight(cco.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd = TRUE)
for(i in 1:length(cco.list)){
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(cco.list)[i]))
}

# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_celltype_CDH2_signaling_pathway_comparison 6*12

netAnalysis_contribution(s.cco.list[[1]], signaling = pathways.show)
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_CDH2_signaling_pathway_contribution 4*6

## CD86 signaling pathway 
pathways.show <- c("CD86")

weight.max <- getMaxWeight(cco.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd = TRUE)
for(i in 1:length(cco.list)){
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(cco.list)[i]))
}

# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_celltype_CD86_signaling_pathway_comparison 6*12

MB.s.injury = MB.injury[,MB.injury$celltype %in% c("Myeloid_cell","T_cell")]
# MB.s.injury = MB.s.injury[,MB.s.injury$subtype != "GLUT_2"]
MB.s.SCFA = MB.SCFA[,MB.SCFA$celltype %in% c("Myeloid_cell","T_cell")]
MB.s.SCFA = MB.s.SCFA[,MB.s.SCFA$subtype != "Mature_B"]

## generate cellchat object for subtype
cco.injury <- createCellChat(MB.s.injury@assays$SCT@data, meta = MB.s.injury@meta.data, group.by = "subtype")
cco.SCFA <- createCellChat(MB.s.SCFA@assays$SCT@data, meta = MB.s.SCFA@meta.data, group.by = "subtype")

# cellchat analysis for injury_control object
s.cellchat <- cco.injury
s.cellchat@DB <- CellChatDB.mouse
s.cellchat <- subsetData(s.cellchat)
s.cellchat <- identifyOverExpressedGenes(s.cellchat)
s.cellchat <- identifyOverExpressedInteractions(s.cellchat)
s.cellchat <- computeCommunProb(s.cellchat, raw.use = T, population.size = T)
s.cellchat <- computeCommunProbPathway(s.cellchat)
s.cellchat <- aggregateNet(s.cellchat)
s.cellchat <- netAnalysis_computeCentrality(s.cellchat, slot.name = "netP")
cco.injury <- s.cellchat

# ellchat analysis for injury_SCFA object
s.cellchat <- cco.SCFA
s.cellchat@DB <- CellChatDB.mouse
s.cellchat <- subsetData(s.cellchat)
s.cellchat <- identifyOverExpressedGenes(s.cellchat)
s.cellchat <- identifyOverExpressedInteractions(s.cellchat)
s.cellchat <- computeCommunProb(s.cellchat, raw.use = T, population.size = T)
s.cellchat <- computeCommunProbPathway(s.cellchat)
s.cellchat <- aggregateNet(s.cellchat)
s.cellchat <- netAnalysis_computeCentrality(s.cellchat, slot.name = "netP")
cco.SCFA <- s.cellchat

s.cco.list <- list(injury_control = cco.injury, # control group
                   injury_SCFA = cco.SCFA)      # test group
s.cellchat <- mergeCellChat(s.cco.list, add.names = names(s.cco.list), cell.prefix = T)

# visualize PECAM1 pathway for Astrocyte - Myeloid cell subtype
weight.max <- getMaxWeight(s.cco.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd = TRUE)
for(i in 1:length(s.cco.list)){
  netVisual_aggregate(s.cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(s.cco.list)[i]))
}
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_Myeloid_T_subtype_CD86_signaling_pathway_comparison 7*14


netAnalysis_contribution(s.cco.list[[1]], signaling = pathways.show)
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_CD86_signaling_pathway_contribution 4*6

## CD86 signaling pathway 
pathways.show <- c("COLLAGEN")
weight.max <- getMaxWeight(cco.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd = TRUE)
for(i in 1:length(cco.list)){
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(cco.list)[i]))
}

# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_celltype_COLLAGEN_signaling_pathway_comparison 6*12

MB.s.injury = MB.injury[,MB.injury$celltype %in% c("Endothelial_cell","Choroid_plexus_epithelial_cell","Vascular_smooth_muscle_cells","Pericyte")]
# MB.s.injury = MB.s.injury[,MB.s.injury$subtype != "GLUT_2"]
MB.s.SCFA = MB.SCFA[,MB.SCFA$celltype %in% c("Endothelial_cell","Choroid_plexus_epithelial_cell","Vascular_smooth_muscle_cells","Pericyte")]

## generate cellchat object for subtype
cco.injury <- createCellChat(MB.s.injury@assays$SCT@data, meta = MB.s.injury@meta.data, group.by = "subtype")
cco.SCFA <- createCellChat(MB.s.SCFA@assays$SCT@data, meta = MB.s.SCFA@meta.data, group.by = "subtype")

# cellchat analysis for injury_control object
s.cellchat <- cco.injury
s.cellchat@DB <- CellChatDB.mouse
s.cellchat <- subsetData(s.cellchat)
s.cellchat <- identifyOverExpressedGenes(s.cellchat)
s.cellchat <- identifyOverExpressedInteractions(s.cellchat)
s.cellchat <- computeCommunProb(s.cellchat, raw.use = T, population.size = T)
s.cellchat <- computeCommunProbPathway(s.cellchat)
s.cellchat <- aggregateNet(s.cellchat)
s.cellchat <- netAnalysis_computeCentrality(s.cellchat, slot.name = "netP")
cco.injury <- s.cellchat

# ellchat analysis for injury_SCFA object
s.cellchat <- cco.SCFA
s.cellchat@DB <- CellChatDB.mouse
s.cellchat <- subsetData(s.cellchat)
s.cellchat <- identifyOverExpressedGenes(s.cellchat)
s.cellchat <- identifyOverExpressedInteractions(s.cellchat)
s.cellchat <- computeCommunProb(s.cellchat, raw.use = T, population.size = T)
s.cellchat <- computeCommunProbPathway(s.cellchat)
s.cellchat <- aggregateNet(s.cellchat)
s.cellchat <- netAnalysis_computeCentrality(s.cellchat, slot.name = "netP")
cco.SCFA <- s.cellchat

s.cco.list <- list(injury_control = cco.injury, # control group
                   injury_SCFA = cco.SCFA)      # test group
s.cellchat <- mergeCellChat(s.cco.list, add.names = names(s.cco.list), cell.prefix = T)

# visualize COLLAGEN pathway for Astrocyte - Myeloid cell subtype
weight.max <- getMaxWeight(s.cco.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd = TRUE)
for(i in 1:length(cco.list)){
  netVisual_aggregate(s.cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(s.cco.list)[i]))
}
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_SCFA_CPC_Vascular_subtype_COLLAGEN_signaling_pathway_comparison 7*14


netAnalysis_contribution(s.cco.list[[2]], signaling = pathways.show)
# save /ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/12_MB_injury_control_COLLAGEN_signaling_pathway_contribution 4*6
