
library(dplyr)
library(Seurat)
library(tidyr)
library(ggplot2)
#UMAP Plotting

pdf("Umap_integrated_harmony.pdf", width=10)
DimPlot(leukemia_Integ_obj,  reduction = "umap.harmony", 
        group.by = c("DataType","harmony_clusters"), combine = FALSE, label = TRUE)
dev.off() 

pdf("Umap_integrated_harmony_DataType.pdf", width=10)
DimPlot(leukemia_Integ_obj,  reduction = "umap.harmony", 
        group.by = c("DataType"), combine = FALSE, label = FALSE)
dev.off() 

pdf("Umap_integrated_DT.pdf", width=10)

DimPlot(leukemia_Integ_obj, reduction = "umap.harmony", group.by = 'DataType', repel = TRUE, 
        cols= c("gray","gray","gray","#88CEE9"),
        order = c("Normal", "AML", "B-ALL", "T-ALL"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
DimPlot(blast, reduction = "umap.harmony", group.by = 'DataType', repel = TRUE, 
        cols= c("gray","gray","gray","#FDDAB9"),
        order = c("AML", "B-ALL", "T-ALL","Normal"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
DimPlot(blast, reduction = "umap.harmony", group.by = 'DataType', repel = TRUE, 
        cols= c("gray","gray","gray","#F8B8C1"),
        order = c("B-ALL", "T-ALL","Normal", "AML"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
DimPlot(blast, reduction = "umap.harmony", group.by = 'DataType', repel = TRUE, 
        cols= c("gray","gray","gray","#DD9FDC"),
        order = c( "T-ALL","Normal", "AML", "B-ALL"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
dev.off()
##########
#Dot plot
DefaultAssay(leukemia_Integ_obj) <- "RNA"
pdf("Markers_Annotated_cluster.pdf", width=15)
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
###########

#Cell type specific expression of DRGs

library(readxl)
library(ggplot2)
library(dplyr)
library(stringr)

directory <- "..."

# Get a list of all xlsx files in the directory
files <- list.files(directory, pattern = "*.xlsx", full.names = TRUE)
#
#store combined data
combined_data <- data.frame()

# Loop through each file, read the data, and add a column for cell type

for (file in files) {
  data <- read_excel(file)
  cell_type <- str_extract(file, "(?<=disulfid_)[^_]+(?=_TALL_normal)")
  
  if ("p_val" %in% colnames(data)) {
    data <- data %>% mutate(p_val = as.numeric(p_val))
  }
  if ("p_val_adj" %in% colnames(data)) {
    data <- data %>% mutate(p_val_adj = as.numeric(p_val_adj))
  }
  
  data <- data %>%
    mutate(cell_type = cell_type) %>%
    filter(avg_log2FC > 1 | avg_log2FC < -1)
  
  combined_data <- bind_rows(combined_data, data)
}

#Dot plotting

pdf("TALL-Normal_disulfid.pdf", height = 6, width = 8) 

# Define breaks and set limits for the size scale
breaks <- c(0, 50, 100, 150, 200, 250)

# Add a column to flag genes with p_val_adj = 0
combined_data <- combined_data %>%
  mutate(is_zero = ifelse(p_val_adj == 0, "Zero p_val_adj", "Non-zero p_val_adj")) 

# Filter out rows where cell_type is NA
#combined_data <- combined_data %>%
  #filter(!is.na(cell_type))

############
# Create the dot plot with the triangle legend positioned under −log10(p_val_adj)
# Main points for non-zero p_val_adj
# Triangles for zero p_val_adj
# However, there was not any DRG in our results with zero p_val_adj

ggplot(combined_data, aes(x = cell_type, y = Marker)) +
  geom_point(aes(size = -log10(p_val_adj), color = avg_log2FC), na.rm = TRUE) +
  geom_point(data = combined_data %>% filter(is_zero == "Zero p_val_adj"),
             aes(x = cell_type, y = Marker, color = avg_log2FC, shape = is_zero), 
             size = 12, inherit.aes = FALSE) +
  # Color gradient for avg_log2FC
  scale_color_gradientn(colors = c("blue", "green", "orange", "red")) +
  scale_size(range = c(1, 10), breaks = breaks, limits = range(breaks)) +
  # Add shape scale for triangles
  scale_shape_manual(values = c("Zero p_val_adj" = 19), labels = c("p_val_adj = 0")) +
  theme_minimal() +
  labs(title = "TALL vs Normal",
       x = "Cell Type",
       y = "Gene",
       size = expression(-log[10](p_val_adj)),
       color = "avg_log2FC",
       shape = NULL) +
  guides(size = guide_legend(order = 1,  # −log10(p_val_adj) legend first
                             shape = guide_legend(order = 2,  # Triangle legend second
                                                  override.aes = list(size = 12, shape = 19, color = "black")), 
                             color = guide_colorbar(order = 3)))
dev.off()