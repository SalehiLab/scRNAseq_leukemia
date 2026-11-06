#monocle3 trajectory inference

library(Seurat)
library(monocle3)
library(Matrix)

obj <- readRDS("..")

DefaultAssay(obj) <- "RNA"

Layers(obj[["RNA"]])

rna_counts <- GetAssayData(
  obj,
  assay = "RNA",
  layer = "counts"
)

cell_metadata <- obj@meta.data

head(cell_metadata)

gene_metadata <- data.frame(
  gene_short_name = rownames(rna_counts),
  row.names = rownames(rna_counts)
)

cds <- new_cell_data_set(
  expression_data = rna_counts,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)


umap <- Embeddings(obj, "umap.harmony")
all(colnames(cds) == rownames(umap))
reducedDims(cds)$UMAP <- umap
reducedDims(cds)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(cds, alignment_group = "batch")
plot_cells(cds, color_cells_by = "subcluster")

cds <- learn_graph(cds)


p1 <- plot_cells(
  cds,
  color_cells_by = "subcluster",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
)

p1$layers[[3]]$aes_params$linewidth <- 1.5
p1

cds <- order_cells(cds)

p2 <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  graph_label_size = 2
)

p2$layers[[3]]$aes_params$linewidth <- 1.5
p2

saveRDS(cds, "blast_cds.rds")
############################################################################
#Blast subclusters pseudo time analysis

cds <- readRDS("blast_cds.rds")

df <- as.data.frame(colData(cds))
df$pseudotime <- pseudotime(cds)

ggplot(df, aes(x = subcluster, y = pseudotime, fill = subcluster)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 0.5) +
  theme_classic() +
  labs(
    x = "Blast subcluster",
    y = "Pseudotime"
  )
#########
#cell cycle score

ggplot(df, aes(x = pseudotime, y = S.Score))+
  geom_point(alpha = 0.5, size = 0.6) +
  geom_smooth(method = "loess", se = TRUE) +
  theme_classic() +
  labs(
    x = "Pseudotime",
    y = "S phase score"
  )

ggplot(df, aes(x = pseudotime, y = G2M.Score)) +
  geom_point(alpha = 0.5, size = 0.6) +
  geom_smooth(method = "loess", se = TRUE) +
  theme_classic() +
  labs(
    x = "Pseudotime",
    y = "G2/M score"
  )

ggplot(df, aes(x = pseudotime, y = QuiescenceLike1)) +
  geom_point(alpha = 0.5, size = 0.6) +
  geom_smooth(method = "loess", se = TRUE) +
  theme_classic() +
  labs(
    x = "Pseudotime",
    y = "Quiescence score"
  )

ggplot(df, aes(x = DataType, y = pseudotime, fill = DataType)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 0.5) +
  theme_classic() +
  labs(
    x = "Data_Type",
    y = "Pseudotime"
  )

ggplot(df, aes(x = Phase, y = pseudotime, fill = Phase)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 0.5) +
  theme_classic() +
  labs(
    x = "Phase",
    y = "Pseudotime"
  )

plot_cells(
  cds,
  color_cells_by = "Phase",
  label_cell_groups = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE
)

# Calculate Spearman correlations between pseudotime and cell-cycle/quiescence scores.
cor_results <- data.frame(
  variable = c(
    "S.Score",
    "G2M.Score",
    "QuiescenceLike1"
  ),
  rho = c(
    cor(df$pseudotime, df$S.Score, method = "spearman", use = "complete.obs"),
    cor(df$pseudotime, df$G2M.Score, method = "spearman", use = "complete.obs"),
    cor(df$pseudotime, df$QuiescenceLike1, method = "spearman", use = "complete.obs")
  )
)
##################################
cor_results