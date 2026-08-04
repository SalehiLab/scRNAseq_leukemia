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


###########################Data Integration##########################
setwd("...")

leukemia_Integ_obj <- merge(x = leukemia_objects[[1]], y = leukemia_objects[2:29], 
                            add.cell.ids = c("Normal", "Normal", "Normal", "Normal", "AML", "AML", "AML", "AML", "AML",
                                             "AML", "AML", "AML", "B-ALL", "B-ALL", "B-ALL", "B-ALL", "B-ALL", "B-ALL", "B-ALL",
                                             "T-ALL", "T-ALL", "T-ALL", "T-ALL", "T-ALL", "T-ALL", "T-ALL", "T-ALL", "T-ALL", "T-ALL"),
                            project = "leukemia_Integ_obj")