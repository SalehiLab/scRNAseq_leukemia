library(dplyr)
library(Seurat)
library(tidyr)
library(ggplot2)
setwd("....")


######################################Initial Data######################################
###Loading data

###########################GSE154109##########################
data_dir <- "..."

gsm_ids <- c("GSM4664009_0064", "GSM4664010_3958","GSM4664011_5286","GSM4664012_7903","GSM4664013_1150", "GSM4664014_1197", "GSM4664015_1210", "GSM4664016_1258", 
             "GSM4664017_1323", "GSM4664018_1260", "GSM4664019_1249", "GSM4664020_1528", 
             "GSM4664021_1236", "GSM4664022_1199", "GSM4664023_1167", "GSM4664024_1308", 
             "GSM4664025_1281", "GSM4664026_1277", "GSM4664027_1209")
data_types <- c("Normal", "Normal", "Normal", "Normal", "AML", "AML", "AML", "AML", "AML",
                "AML", "AML", "AML", "B-ALL", "B-ALL", "B-ALL", "B-ALL", "B-ALL", "B-ALL", "B-ALL")

leukemia_objects <- list()

for (i in 1:length(gsm_ids)) {
  gsm <- gsm_ids[i]
  matrix_path <- file.path(data_dir, paste0(gsm, "_matrix.mtx.gz"))
  features_path <- file.path(data_dir, paste0(gsm, "_features.tsv.gz"))
  barcodes_path <- file.path(data_dir, paste0(gsm, "_barcodes.tsv.gz"))
  
  data <- ReadMtx(mtx = matrix_path, features = features_path, cells = barcodes_path)
  
  seurat_obj <- CreateSeuratObject(counts = data, project = gsm, min.cells = 3, min.features = 500)
  seurat_obj$DataType <- data_types[i]
  leukemia_objects[[gsm]] <- seurat_obj
}

###########################GSE227122##########################
data_dir <- "..."
gsm_ids2 <- c("GSM7091996_T1Dx","GSM7091997_T2Dx","GSM7091998_T3Dx","GSM7091999_T4Dx","GSM7092000_T5Dx",
              "GSM7092001_T6Dx","GSM7092003_T7Dx","GSM7092005_T8Dx","GSM7092007_T9Dx","GSM7092009_T10Dx")
data_types <- c("T-ALL", "T-ALL", "T-ALL", "T-ALL", "T-ALL", "T-ALL", "T-ALL", "T-ALL", "T-ALL", "T-ALL")

for (i in 1:length(gsm_ids2)) {
  gsm <- gsm_ids2[i]
  matrix_path <- file.path(data_dir, paste0(gsm, "_matrix.mtx.gz"))
  features_path <- file.path(data_dir, paste0(gsm, "_features.tsv.gz"))
  barcodes_path <- file.path(data_dir, paste0(gsm, "_barcodes.tsv.gz"))
  
  data <- ReadMtx(mtx = matrix_path, features = features_path, cells = barcodes_path)
  seurat_obj <- CreateSeuratObject(counts = data, project = gsm,  min.cells = 3, min.features = 500)
  seurat_obj$DataType <- data_types[i]
  leukemia_objects[[gsm]] <- seurat_obj
}


###########################Merging Data##########################
setwd("...")

leukemia_Integ_obj <- merge(x = leukemia_objects[[1]], y = leukemia_objects[2:29], 
                            add.cell.ids = c("Normal", "Normal", "Normal", "Normal", "AML", "AML", "AML", "AML", "AML",
                                             "AML", "AML", "AML", "B-ALL", "B-ALL", "B-ALL", "B-ALL", "B-ALL", "B-ALL", "B-ALL",
                                             "T-ALL", "T-ALL", "T-ALL", "T-ALL", "T-ALL", "T-ALL", "T-ALL", "T-ALL", "T-ALL", "T-ALL"),
                            project = "leukemia_Integ_obj")


####################Data Integration & Clustering##################

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
