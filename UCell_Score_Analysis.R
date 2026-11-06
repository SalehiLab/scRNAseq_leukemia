#UCell score

library(UCell)

obj <- readRDS("blast_cds.rds")
#disulfidptosis UCell score

disulfidptosis_markers <- c(
  "SLC7A11", "SLC3A2", "RPN1", "NCKAP1", "NUBPL", "NDUFA11", "LRPPRC", "OXSM", "NDUFS1",
  "GYS1", "ACTN4", "ACTB", "CD2AP", "CAPZB", "DSTN", "FLNA", "FLNB", "INF2", "IQGAP1",
  "MYH10", "MYL6", "MYH9", "PDLIM1", "TLN1")
obj <- AddModuleScore_UCell(
  obj,
  features = list(
    Disulfidptosis = disulfidptosis_markers
  )
)
grep(
  "Disulfidptosis",
  colnames(obj@meta.data),
  value = TRUE
)
head(obj)

df <- as.data.frame(colData(cds))
df$pseudotime <- pseudotime(cds)[rownames(df)]
df$Disulfidptosis <- obj$Disulfidptosis_UCell[
  rownames(df)
]

library(ggplot2)

ggplot(
  df,
  aes(x = pseudotime, y = Disulfidptosis)
) +
  geom_point(
    alpha = 0.15,
    size = 0.5
  ) +
  geom_smooth(
    method = "loess",
    se = TRUE
  ) +
  theme_classic() +
  labs(
    x = "Pseudotime",
    y = "Disulfidptosis UCell score"
  )
cor.test(
  df$pseudotime,
  df$Disulfidptosis,
  method = "spearman"
)
ggplot(
  df,
  aes(
    x = subcluster,
    y = Disulfidptosis,
    fill = subcluster
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA
  ) +
  theme_classic() +
  labs(
    x = "Blast subcluster",
    y = "Disulfidptosis UCell score"
  )