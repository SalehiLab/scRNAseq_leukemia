# SCTransform
# PCA
# Harmony integration
# Find Neighbors
# Find Clusters
# Run UMAP
# save integrated object

library(dplyr)
library(Seurat)
library(tidyr)
library(ggplot2)
setwd("....")

####################Dimension Reduction & Clustering##################

leukemia_Integ_obj <- SCTransform(leukemia_Integ_obj, vst.flavor = "v2", verbose = FALSE)
leukemia_Integ_obj <- RunPCA(leukemia_Integ_obj, npcs = 30, verbose = FALSE)

leukemia_Integ_obj <- IntegrateLayers(object = leukemia_Integ_obj, method = HarmonyIntegration,
                                      orig.reduction = "pca", new.reduction = 'harmony',
                                      assay = "SCT", verbose = FALSE)
leukemia_Integ_obj <- FindNeighbors(leukemia_Integ_obj, reduction = "harmony", dims = 1:30)
leukemia_Integ_obj <- FindClusters(leukemia_Integ_obj, resolution = 0.2, cluster.name = "harmony_clusters")


clusdata <-table(Idents(leukemia_Integ_obj), leukemia_Integ_obj$DataType)
write.csv(clusdata, file= "./Results/ClusDataType.csv")

#Maximum modularity in 10 random starts: 0.9554
#Number of communities: 13

#Run UMAP

leukemia_Integ_obj <- RunUMAP(leukemia_Integ_obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

###################################cluster annotation##############################

DefaultAssay(leukemia_Integ_obj) <- "RNA"
pdf("Markers_label_cluster.pdf", width=15)
DotPlot(leukemia_Integ_obj, features = c("SOX4", "CD7","CD99", "ARMH1", "STMN1", "HES4",
                                         "CD68", "CD14", "FCGR3A", "ITGAM", "LYZ", "FCN1",
                                         "TIGIT", "GNLY", "CCL5", "KLRD1", "GZMA", "XCL1", "CD8A", "NKG7",
                                         "CD19", "MS4A1", "CD79A", "CD79B", "IGHM",
                                         "PTPRC", "CD3D", "CD3E", "IL7R", "TRAC",
                                         "GYPA", "HBB", "HBD", "CD34",
                                         "JCHAIN", "GZMB", "IRF8", "LILRA4")
        ,cols=c("blue","red")) + RotatedAxis()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()

current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
new.cluster.ids <- c("BLAST", "BLAST", "BLAST", "BLAST",
                     "Monocyte", "NK", "B-Cell", "T-Cell", "BLAST", "Erythrocyte","T-Cell", "DC", "B-Cell")
Idents(leukemia_Integ_obj) <- plyr::mapvalues(x = Idents(leukemia_Integ_obj), from = current.cluster.ids, to = new.cluster.ids)

cellTypes<-Idents(leukemia_Integ_obj)

saveRDS(leukemia_Integ_obj, "leukemia_Integ_obj.rds")
