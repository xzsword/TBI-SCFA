#
##
###

library(Seurat)
library(ggplot2)
library(dplyr)
library(DoubletFinder)

# data loading
## metadata
metadata = readxl::read_xlsx("/bgfs/gkohanbash/SimonsLab/R_Share/raw_data/Simons_Lab_MetaData.xlsx")

# scRNA-seq data
MB757 = Read10X(data.dir = "/bgfs/gkohanbash/SimonsLab/usftp21.novogene.com/raw_data/MB757/MB757_CKDL210026110-1a-SI_TT_D1_HW7WLDSX2_S4_L004/outs/filtered_feature_bc_matrix")
MB758 = Read10X(data.dir = "/bgfs/gkohanbash/SimonsLab/usftp21.novogene.com/raw_data/MB758/MB758_CKDL210026111-1a-SI_TT_D2_HW7WLDSX2_S8_L004/outs/filtered_feature_bc_matrix")
MB759 = Read10X(data.dir = "/bgfs/gkohanbash/SimonsLab/usftp21.novogene.com/raw_data/MB759/MB759_CKDL210026112-1a-SI_TT_D3_HW7WLDSX2_S9_L004/outs/filtered_feature_bc_matrix")
MB760 = Read10X(data.dir = "/bgfs/gkohanbash/SimonsLab/usftp21.novogene.com/raw_data/MB760/MB760_CKDL210026113-1a-SI_TT_D4_HW7WLDSX2_S13_L004/outs/filtered_feature_bc_matrix")
MB761 = Read10X(data.dir = "/bgfs/gkohanbash/SimonsLab/usftp21.novogene.com/raw_data/MB761/MB761_CKDL210026114-1a-SI_TT_D9_HW7WLDSX2_S12_L004/outs/filtered_feature_bc_matrix")
MB762 = Read10X(data.dir = "/bgfs/gkohanbash/SimonsLab/usftp21.novogene.com/raw_data/MB762/MB762_CKDL210026115-1a-SI_TT_D10_HW7WLDSX2_S11_L004/outs/filtered_feature_bc_matrix")
MB763 = Read10X(data.dir = "/bgfs/gkohanbash/SimonsLab/usftp21.novogene.com/raw_data/MB763/MB763_CKDL210026116-1a-SI_TT_D11_HW7WLDSX2_S10_L004/outs/filtered_feature_bc_matrix")
MB764 = Read10X(data.dir = "/bgfs/gkohanbash/SimonsLab/usftp21.novogene.com/raw_data/MB764/MB764_CKDL210026117-1a-SI_TT_D12_HW7WLDSX2_S3_L004/outs/filtered_feature_bc_matrix")
MB765 = Read10X(data.dir = "/bgfs/gkohanbash/SimonsLab/usftp21.novogene.com/raw_data/MB765/MB765_CKDL210026118-1a-SI_TT_D5_HW7WLDSX2_S16_L004/outs/filtered_feature_bc_matrix")
MB766 = Read10X(data.dir = "/bgfs/gkohanbash/SimonsLab/usftp21.novogene.com/raw_data/MB766/MB766_CKDL210026119-1a-SI_TT_D6_HW7WLDSX2_S7_L004/outs/filtered_feature_bc_matrix")
MB767 = Read10X(data.dir = "/bgfs/gkohanbash/SimonsLab/usftp21.novogene.com/raw_data/MB767/MB767_CKDL210026120-1a-SI_TT_D7_HW7WLDSX2_S6_L004/outs/filtered_feature_bc_matrix")
MB768 = Read10X(data.dir = "/bgfs/gkohanbash/SimonsLab/usftp21.novogene.com/raw_data/MB768/MB768_CKDL210026121-1a-SI_TT_D8_HW7WLDSX2_S5_L004/outs/filtered_feature_bc_matrix")
MB785 = Read10X(data.dir = "/bgfs/gkohanbash/SimonsLab/usftp21.novogene.com/raw_data/MB785/MB785_CKDL210026122-1a-SI_TT_E1_HW7WLDSX2_S1_L004/outs/filtered_feature_bc_matrix")
MB786 = Read10X(data.dir = "/bgfs/gkohanbash/SimonsLab/usftp21.novogene.com/raw_data/MB786/MB786_CKDL210026123-1a-SI_TT_E2_HW7WLDSX2_S14_L004/outs/filtered_feature_bc_matrix")
MB787 = Read10X(data.dir = "/bgfs/gkohanbash/SimonsLab/usftp21.novogene.com/raw_data/MB787/MB787_CKDL210026124-1a-SI_TT_E3_HW7WLDSX2_S2_L004/outs/filtered_feature_bc_matrix")
MB788 = Read10X(data.dir = "/bgfs/gkohanbash/SimonsLab/usftp21.novogene.com/raw_data/MB788/MB788_CKDL210026125-1a-SI_TT_E4_HW7WLDSX2_S15_L004/outs/filtered_feature_bc_matrix")

MB_list = list(MB757,MB758,MB759,MB760,
               MB761,MB762,MB763,MB764,
               MB765,MB766,MB767,MB768,
               MB785,MB786,MB787,MB788)

# annotate the cells with sample info
for (i in 1:16) {
  colnames(MB_list[[i]]) = gsub("-1",paste0("-",i),colnames(MB_list[[i]]))
}

# normalize and identify variable features for each dataset independently: min feature = 500; mt.prep < 20; red.cell.prep < 1
# SCT normalization & doublet finding
MB_list = lapply(MB_list, FUN = function(x){
  x = CreateSeuratObject(counts = x, min.cells = 3, min.features = 500)
  x[["percent.mt"]] = PercentageFeatureSet(x, pattern = "^mt-")
  x <- subset(x, subset = percent.mt < 20)
  x <- SCTransform(x,verbose = F, vars.to.regress = "percent.mt",return.only.var.genes = F)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = MB_list,nfeatures = 2000) # combined the most variable gene to 3000 from all samples
MB_list <- PrepSCTIntegration(object.list = MB_list, anchor.features = features)# set the Most variable gene to former 3000 genes in all samples
save(MB_list,file = "/bgfs/gkohanbash/SimonsLab/R_Share/processing/1_integration_input.RData")

# performing integration
anchors <- FindIntegrationAnchors(object.list = MB_list, anchor.features = features,normalization.method = "SCT")
save(anchors,file = "/bgfs/gkohanbash/SimonsLab/R_Share/processing/1_integrate_anchors.RData")

MB.combined <- IntegrateData(anchorset = anchors,normalization.method = "SCT")
save(MB.combined,file = "/bgfs/gkohanbash/SimonsLab/R_Share/processing/1_integrateData.RData")

sample = unlist(lapply(strsplit(colnames(MB.combined),"-"),"[",2))
MB.combined$sample = sample
MB.combined$batch = if_else(sample %in% c(1:4), "batch1",
                            if_else(sample %in% c(5:8),"batch2",
                                    if_else(sample %in% c(9:12),"batch3","batch4")))

MB.combined$group = if_else(sample %in% c(1:4), "group1",
                            if_else(sample %in% c(5:8),"group2",
                                    if_else(sample %in% c(9:12),"group3","group4")))

# UMAP for integration data
DefaultAssay(MB.combined) <- "integrated"

MB.combined = RunPCA(MB.combined,npcs = 50, features = rownames(MB.combined))
save(MB.combined,file = "/bgfs/gkohanbash/SimonsLab/R_Share/processing/1_integrateData_PCA.RData")

ElbowPlot(MB.combined, ndims=50)
pc.num=1:30

MB.combined = RunUMAP(MB.combined , reduction = "pca", dims = pc.num)

# cluster cells
MB.combined = FindNeighbors(MB.combined, verbose = FALSE, reduction = "pca",dims = pc.num)
MB.combined = FindClusters(MB.combined,resolution = 0.5)

save(MB.combined,file = "/bgfs/gkohanbash/SimonsLab/R_Share/processing/1_integrateData_final.RData")

# plan("sequential") # back to single core calculation