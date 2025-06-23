# Load libraries
library(Seurat)
library(ggplot2)

#Load Data
level3.seurat <- readRDS("./PATH/TO/DIRECTORY/fibroblast-lean-seurat-20240308.RDS") # SCT transformed preprocessed fibroblast population

# Assigning new cluster identities
level3.seurat_L <- FindClusters(level3.seurat, resolution = 0.3)

# Assigning new cluster identities
new.cluster.ids <- c("Myofibroblasts", "PLIN2 High Fibroblasts", "Fibroblasts", "Fibroblasts", "Myofibroblasts", "HAS1 High Fibroblasts", "PLIN2 High Fibroblasts",
                     "Myofibroblasts", "Myofibroblasts")
names(new.cluster.ids) <- levels(level3.seurat_L)
level3.seurat_L <- RenameIdents(level3.seurat_L, new.cluster.ids)

# Define the new order of clusters
new_order <- c("HAS1 High Fibroblasts", "PLIN2 High Fibroblasts", "Myofibroblasts", "Fibroblasts")

# Reorder the levels of the Idents
level3.seurat_L@active.ident <- factor(level3.seurat_L@active.ident, levels = new_order)

# Simplify the label
DefaultAssay(level3.seurat_L) <- "RNA" # rds file is at SCT assay
new.cluster.ids <- c(
  "Fibroblasts" = "Fibroblasts",
  "Myofibroblasts" = "Myofibroblasts",
  "PLIN2 High Fibroblasts" = "PLIN2+",
  "HAS1 High Fibroblasts" = "HAS1 High"
)

level3.seurat_L <- RenameIdents(level3.seurat_L, new.cluster.ids)

# Organize Idents for the plot later (Default)
Idents(level3.seurat_L) <- factor(Idents(level3.seurat_L), levels = c("HAS1 High", "PLIN2+", "Myofibroblasts", "Fibroblasts"))

# Normalize RNA assay
level3.seurat_L <- NormalizeData(level3.seurat_L)
level3.seurat_L <- ScaleData(level3.seurat_L)
level3.seurat_L <- FindVariableFeatures(level3.seurat_L)

# Organize Idents for Coloring the plot
Idents(level3.seurat_L) <- factor(Idents(level3.seurat_L), levels = c("Myofibroblasts", "Fibroblasts","HAS1 High", "PLIN2+"))

# Figure 5E
scDimplot <- DimPlot(level3.seurat_L, split.by = "condition", cols = "alphabet")+ xlim(-8,8) + ylim(-5,5)
ggsave("./output/2025-05-09/scDimplot_level3_split_diff_cols.png",scDimplot, width = 25, height = 6, units = "cm")

# Figure S15A
scDimplot <- DimPlot(level3.seurat_L, split.by = "condition", cols = "alphabet") + NoLegend() + xlim(-8,8) + ylim(-5,5)
ggsave("./output/2025-05-12/scDimplot_level3_split_noLegend_diff_cols.png",scDimplot, width = 20, height = 6, units = "cm")

# Re-Organize Idents for the plot later (Default)
Idents(level3.seurat_L) <- factor(Idents(level3.seurat_L), levels = c("HAS1 High", "PLIN2+", "Myofibroblasts", "Fibroblasts"))
