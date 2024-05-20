#
##
###

#cell type identification
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)

# load data
load("/bgfs/gkohanbash/SimonsLab/R_Share/processing/2_integrateData_doublet_removed.RData")

# cell type identification -- Mouse Cell Atlas
MCA_exp <- readRDS("/bgfs/gkohanbash/Zujuan_dhr11/TBI/ref/cellType_identification/MCA_merged_mat.rds")
MCA_celltypeinfo <- read.csv("/bgfs/gkohanbash/Zujuan_dhr11/TBI/ref/cellType_identification/MCA_CellAssignments.csv")
MCA_cellname <- read.csv("/bgfs/gkohanbash/Zujuan_dhr11/TBI/ref/cellType_identification/MCA_All-batch-removed-assignments.csv")
index <- match(MCA_cellname$Cell.name, MCA_celltypeinfo$Cell.name)
selected_cellname <- MCA_celltypeinfo[na.omit(index), ]
dim(selected_cellname)
brain_cellname <- selected_cellname[selected_cellname$Tissue %in% c("Brain","Peripheral_Blood"), ]
brain_cellname = brain_cellname[brain_cellname$Annotation != "Schwann cell(Brain)",]
brain_cellname = brain_cellname[brain_cellname$Annotation != "Erythroblast_Car2 high(Peripheral_Blood)",]
brain_cellname = brain_cellname[brain_cellname$Annotation != "Erythroblast_Hba-a2 high(Peripheral_Blood)",]
brain_cellname = brain_cellname[brain_cellname$Annotation != "Pan-GABAergic(Brain)",]
brain_cellname = brain_cellname[brain_cellname$Annotation != "Granulocyte_Il33 high(Brain)",]
brain_cellname = brain_cellname[brain_cellname$Annotation != "Granulocyte_Ngp high(Brain)",]

brain_exp <- MCA_exp[, as.character(brain_cellname$Cell.name)]

all.equal(as.character(brain_cellname$Cell.name), colnames(brain_exp))
MCA_brain_rds <- list()
MCA_brain_rds$exp <- brain_exp
MCA_brain_rds$meta.data <- brain_cellname
saveRDS(MCA_brain_rds, "/bgfs/gkohanbash/SimonsLab/R_Share/processing/3_MCA_brain_exp.rds")

# build seurat reference for MCA 
ref_exp <- brain_exp
ref.seurat <- CreateSeuratObject(counts = ref_exp)
ref.seurat$celltype <- as.character(brain_cellname$Annotation)
ref.seurat[["percent.mt"]] = PercentageFeatureSet(ref.seurat, pattern = "^mt-")
ref.seurat <- subset(ref.seurat, subset = percent.mt < 20)
ref.seurat <- SCTransform(ref.seurat,verbose = F, vars.to.regress = "percent.mt",return.only.var.genes = F)

ref.seurat <- RunPCA(object = ref.seurat, npcs = 30, verbose = FALSE)
saveRDS(ref.seurat, file = "/bgfs/gkohanbash/SimonsLab/R_Share/processing/3_MCA_NormalizeData_ref_seurat.rds")

test.anchors <- FindTransferAnchors(reference = ref.seurat, 
                                    query = MB.combined, 
                                    dims = 1:30)
predictions_MCA <- TransferData(anchorset = test.anchors, 
                                refdata = ref.seurat$celltype, 
                                dims = 1:30)
MB.combined$predictions_MCA <- predictions_MCA$predicted.id

# integrate the celltypes
MB.combined$celltype = MB.combined$predictions_MCA

MB.combined$celltype[grep("Astro",MB.combined$predictions_MCA)] = "Astrocyte"
MB.combined$celltype[grep("B cell",MB.combined$predictions_MCA)] = "B_cell"
MB.combined$celltype[grep("Basop",MB.combined$predictions_MCA)] = "Basophil"
MB.combined$celltype[grep("Dendritic",MB.combined$predictions_MCA)] = "Dendritic_cell"
MB.combined$celltype[grep("ependy",MB.combined$predictions_MCA)] = "Ependymal_cell"
MB.combined$celltype[grep("Macrophage",MB.combined$predictions_MCA)] = "Macrophage"
MB.combined$celltype[grep("Microglia",MB.combined$predictions_MCA)] = "Microglia"
MB.combined$celltype[grep("Monocyte",MB.combined$predictions_MCA)] = "Monocyte"
MB.combined$celltype[grep("Neutro",MB.combined$predictions_MCA)] = "Neutrophil"
MB.combined$celltype[grep("NK",MB.combined$predictions_MCA)] = "NK_cell"
MB.combined$celltype[grep("precursor",MB.combined$predictions_MCA)] = "Oligodendrocyte_precursor_cell"
MB.combined$celltype[grep("Myelina",MB.combined$predictions_MCA)] = "Oligodendrocyte"
MB.combined$celltype[grep("GABA",MB.combined$predictions_MCA)] = "Neuron"
MB.combined$celltype[grep("T cell",MB.combined$predictions_MCA)] = "T_cell"
MB.combined$celltype[grep("Plas",MB.combined$predictions_MCA)] = "Plasma_cell"
MB.combined$celltype[grep("Neuron",MB.combined$predictions_MCA)] = "Neuron"

celltype = MB.combined$celltype
save(predictions_MCA,celltype,file = "/bgfs/gkohanbash/SimonsLab/R_Share/processing/3_Mouse_Brain_cell_type_MCA_ref_identification.RData")

MB.combined$MCA_celltype = celltype

# cell type identify by cell marker
# get the marker gene for each clusters
MB.markers <- FindAllMarkers(object = MB.combined, only.pos = TRUE, min.pct = 0.25, 
                             thresh.use = 0.25)
save(MB.markers,file = "/bgfs/gkohanbash/SimonsLab/R_Share/processing/3_cluster_markers.RData")

# get top 10 marker genes in each cluster
top5 <- MB.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
top2 = MB.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)

# visulize the marker gene by heatmap
DoHeatmap(object = MB.combined, features = top5$gene,assay = "SCT")

# visulize the marker gene by dotplot
DefaultAssay(MB.combined) = "SCT"
cellMarker = c("Sall1","Tmem119","Dock2","Cx3cr1",
               "Ccr2", "F13a1","Itgam",
               "Flt3","Zbtb46",
               "Cd19","Cd3e","Ncr1",
               "Atp1b1","Rbfox3","Camk2a","Kif5c",
               "Aldoc","Slc1a3","Slc7a10","Gja1",
               "Mog","Olig1","Enpp2",
               "Calml4","Rarres2","Dnali1","Tmem212",
               "Vcan","Pdgfra",
               "Pecam1","Rgs5","Ly6c1",
               "Vtn","Tagln","Acta2",
               "Ttr",              # Choroid_plexus_epithelial_cells
               "Kcnj8",             # pericyte
               "Col1a1",           # fibroblasts, one of the mural cells
               "Cdk1",             # Neuron_restricted precursor
               "Sox11"             # immature neuron
               )
setdiff(cellMarker,rownames(MB.combined))

DotPlot(MB.combined, features = cellMarker)+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
# save 3_seurat_cluster_cell_marker_dotplot 8*12

# cell marker
# microglia/macrophage:  Sall1,Tmem119 Dock2 Cx3cr1
# monocyte: Ccr2, F13a1, Itgam(CD11b)
# dendritic cell marker: Flt3 Zbtb46
# B: Cd19
# T: CD3e
# NK: Ncr1
# Neutrophil: Ly6g 
# Neuron:  Atp1b1 Rbfox3(NeuN) "Camk2a","Kif5c",
# Astrocyte: Aldoc Slc1a3 "Slc7a10",Gja1
# Oligodendrocyte: Mog, Olig1,Enpp2
# Ependymal: Rarres2, Tmem212, Dnali1,Calml4
# OPC: Vcan Pdgfra
# Endothelial cell: Pecam1, Rgs5,Ly6c1
# mural cell: "Vtn","Tagln","Acta2"
# Vascular_smooth_muscle_cells: "Acta2"
# Choroid_plexus_epithelial_cells: Ttr
# pericyte: Kcnj8
# Neuronal_restricted_precursor: Cdk1
# immature_neuron: Sox11

# Marker source:
# 1. CellMarker web
# 2. Single cell molecular alterations reveal target cells and pathways of concussive brain injury (Nat. comm)
# 3. CCR2 deficiency alters activation of microglia subsets in traumatic brain injury (Cell. Rep)
# 4. Detecting Activated Cell Populations Using Single-Cell RNA-Seq (Neuron)
# 5. Ximerakis, M., Lipnick, S.L., Innes, B.T. et al. Single-cell transcriptomic profiling of the aging mouse brain. Nat Neurosci 22, 1696–1708 (2019). https://doi.org/10.1038/s41593-019-0491-3

# validate the identification results
# load cell type ref markers from paper ref: Single cell molecular alterations reveal target cells and pathways of concussive brain injury (Nat. comm) 
#                                            Detecting Activated Cell Populations Using Single-Cell RNA-Seq (Neuron)

# ref = top50 logFC in Nat.comm + all in Neuron paper
cell_marker_ref = readxl::read_xlsx("/bgfs/gkohanbash/SimonsLab/R_Share/raw_data/cellmarker_ref.xlsx",skip = 2)
cell_marker_ref = cell_marker_ref [,9:16]
cell_marker_ref = as.data.frame(cell_marker_ref)
cell_marker_ref = cell_marker_ref[1:50,]
cell_marker_ref = cell_marker_ref[,c("Neurons","Astrocytes","Microglia","Oligodendrocyte PCs",
                                     "Oligodendrocytes","Endothelial","Mural","Ependymal")]

cell_marker_1 = readxl::read_xlsx("/bgfs/gkohanbash/SimonsLab/R_Share/raw_data/mmc2.xlsx",skip = 4)
cell_marker_1 = cell_marker_1[,c(1,5,9,13,17,21,25)]
colnames(cell_marker_1) = c("Neurons","Astrocytes","Microglia","Oligodendrocyte PCs",
                            "Oligodendrocytes","Endothelial","Mural")
cell_marker_1 = as.data.frame(cell_marker_1)
cell_marker_1$Ependymal = NA

cell_marker_ref = rbind(cell_marker_ref,cell_marker_1)
# count cluster markers and ref markers overlap and evaluate by Fisher exact test

overlap_count = data.frame()
overlap_fdr = data.frame()
overlap_p = vector(mode = "numeric")

cluster = unique(MB.markers$cluster)

for (i in 1:ncol(cell_marker_ref)) {
  for (j in 1:length(cluster)) {
    ref = na.omit(cell_marker_ref[,i])
    test = MB.markers$gene[which(MB.markers$cluster == cluster[j])]
    o = length(intersect(ref,test))    # the number of marker genes that in ref gene list 
    k = length(ref)            # ref gene count
    n = nrow(MB.combined)      # total gene count
    m = length(test)           # the number of input DEGs
    overlap_p = c(overlap_p,1-phyper(o-1,k,n-k,m))  # Fisher exact test
    overlap_count[j,i] = o
  }
}

rownames(overlap_count) = cluster#paste0("cluster_",cluster)
colnames(overlap_count) = colnames(cell_marker_ref)

overlap_fdr = p.adjust(overlap_p, method = "BH", n = length(overlap_p))
# dim(overlap_fdr) = c(29,8)
# rownames(overlap_fdr) = paste0("cluster_",cluster)
# colnames(overlap_fdr) = colnames(cell_marker_ref)

overlap_count$cluster = rownames(overlap_count)
hm <- overlap_count %>% gather(Celltype, value, 1:8)
overlap_fdr = if_else(overlap_fdr > 0.1,0,
                      if_else(overlap_fdr > 0.001, 1,
                              if_else(overlap_fdr > 0.00001, 2,
                                      if_else(overlap_fdr > 1e-07,3,
                                              if_else(overlap_fdr > 1e-09,4,
                                                      if_else(overlap_fdr > 1e-11,5,6))))))

# buuild confusion matrix input
hm$fdr = overlap_fdr

# plot the confusion matrix
ggplot(hm, aes(x=cluster, y=Celltype, fill=fdr)) +
  geom_tile() + theme_bw() + coord_equal() +
  scale_fill_distiller(palette="Greens", direction=1) +
  #guides(fill=F) + # removing legend for `fill`
  labs(title ="Value distribution") + # using a title instead
  geom_text(aes(label=value), color="black") # printing values
# color = FDR value
# number in cell = number of cluster marker gene overlaped with ref cell markers
# FDR: 0 = fdr -> (0.1,1]
#      1 = fdr -> (0.001,0.1]
#      2 = fdr -> (1e-05,1e-03]
#      3 = fdr -> (1e-07,1e-05]
#      4 = fdr -> (1e-09,1e-07]
#      5 = fdr -> (1e-11,1e-09]
#      6 = fdr -> [0,1e-11]

# combining the results from MCA auto-annotation and marker genes + ref marker gene list validation to annotate the cell type
library(plyr)
celltype_final = as.character(MB.combined$seurat_clusters)

celltype_final = mapvalues(celltype_final,"0","Myeloid_cell")
celltype_final = mapvalues(celltype_final,"1","Endothelial_cell")
celltype_final = mapvalues(celltype_final,"2","Choroid_plexus_epithelial_cell")
celltype_final = mapvalues(celltype_final,"3","Neuron")
celltype_final = mapvalues(celltype_final,"4","Myeloid_cell")
celltype_final = mapvalues(celltype_final,"5","Myeloid_cell")
celltype_final = mapvalues(celltype_final,"6","Astrocyte")
celltype_final = mapvalues(celltype_final,"7","Choroid_plexus_epithelial_cell")
celltype_final = mapvalues(celltype_final,"8","Oligodendrocyte")
celltype_final = mapvalues(celltype_final,"9","Pericyte")
celltype_final = mapvalues(celltype_final,"10","T_cell")
celltype_final = mapvalues(celltype_final,"11","Astrocyte")
celltype_final = mapvalues(celltype_final,"12","OPC")
celltype_final = mapvalues(celltype_final,"13","Choroid_plexus_epithelial_cell")
celltype_final = mapvalues(celltype_final,"14","Choroid_plexus_epithelial_cell")
celltype_final = mapvalues(celltype_final,"15","Myeloid_cell")
celltype_final = mapvalues(celltype_final,"16","Endothelial_cell")
celltype_final = mapvalues(celltype_final,"17","Pericyte")
celltype_final = mapvalues(celltype_final,"18","Vascular_smooth_muscle_cells")
celltype_final = mapvalues(celltype_final,"20","Endothelial_cell")
celltype_final = mapvalues(celltype_final,"21","Ependymocyte")
celltype_final = mapvalues(celltype_final,"22","Endothelial_cell")
celltype_final = mapvalues(celltype_final,"23","Choroid_plexus_epithelial_cell")
celltype_final = mapvalues(celltype_final,"24","Neuron")
celltype_final = mapvalues(celltype_final,"25","B_cell")
celltype_final = mapvalues(celltype_final,"27","Myeloid_cell")
celltype_final = mapvalues(celltype_final,"29","Oligodendrocyte")
celltype_final = mapvalues(celltype_final,"30","T_cell")
celltype_final = mapvalues(celltype_final,"31","Myeloid_cell")

MB.combined$celltype = celltype_final
MB.combined = MB.combined[,MB.combined$group != "group1"]

UMAPPlot(MB.combined,group.by = "celltype",label = T) + NoLegend()
ggsave("/ix/gkohanbash/SimonsLab/R_Share/fig/No_acetate_group/3_Umap_with_final_cellType_annotation.pdf",height = 5,width = 5)
# save 3_Umap_with_final_cellType_annotation 4*6

save(MB.combined,file = "/ix/gkohanbash/SimonsLab/R_Share/processing/No_acetate_group/3_integrateData_with_cellType_annotation.RData")
