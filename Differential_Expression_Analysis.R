# Differential expression analysis

library(dplyr)
library(Seurat)
library(tidyr)

setwd("....")


leukemia_Integ_obj <- readRDS("leukemia_Integ_obj.rds")

##################Differential Expression Analysis#############

MIN_LOGFOLD_CHANGE = 0.25 
MIN_PCT_CELLS_EXPR_GENE = 0.1
leukemia_Integ_obj <- PrepSCTFindMarkers(leukemia_Integ_obj)

# minimum percent of cells that must express gene in either cluster
all.markers = FindAllMarkers(leukemia_Integ_obj,
                             min.pct = MIN_PCT_CELLS_EXPR_GENE,
                             logfc.threshold = MIN_LOGFOLD_CHANGE,
                             slot = "data",
                             test.use = "wilcox",
                             return.thresh=0.05,
                             verbose = T,
                             only.pos = TRUE)
# sort all the markers by p-value
all.markers.sortedByPval = all.markers[order(all.markers$p_val),]

#the top most significant markers
head(all.markers.sortedByPval)

###
library(dplyr)
top500 <- all.markers.sortedByPval %>%  group_by(cluster)  %>% do(head(., n=500))
library(writexl)
write_xlsx(top500,"differentially_expressed_cluster_500.xlsx")

#######################

DefaultAssay(leukemia_Integ_obj) <- "RNA"

Idents(leukemia_Integ_obj) <- "cell_annotation"

# Differential expression analysis was performed independently for each annotated cell type.
# The same workflow was applied to all leukemic conditions (AML, B-ALL, and T-ALL).
# To compare each disease with its corresponding normal counterpart.

markers <- FindMarkers(
  leukemia_Integ_obj,
  assay = "RNA",
  slot = "data",
  ident.1 = "B-ALL",
  ident.2 = "Normal",
  group.by = "DataType",
  subset.ident = "B-Cell",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  test.use = "wilcox"
)

markers$gene <- rownames(markers)
markers <- markers[, c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]

#for each annotated cell type and datatype
#For instance for B-ALL 
library(writexl)
write_xlsx(markers, "markers_bcell_Normal_BALL.xlsx")

