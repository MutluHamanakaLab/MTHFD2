# Load Libraries
library(ggplot2)
library(dplyr)
library(magrittr)
library(patchwork)
library(viridis)
library(Seurat)
library(scCustomize)
library(qs)

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

# Organize Idents for Coloring of the plot
Idents(level3.seurat_L) <- factor(Idents(level3.seurat_L), levels = c("Myofibroblasts", "Fibroblasts","HAS1 High", "PLIN2+"))

# Figure 5E
scDimplot <- DimPlot(level3.seurat_L,cols = "alphabet") + NoLegend()
ggsave("./PATH/TO/DIRECTORY/scDimplot_level3_noLegend_diff_cols_NoLegend.png",scDimplot, width = 10, height = 6, units = "cm")

scDimplot <- DimPlot(level3.seurat_L,cols = "alphabet")
ggsave("./PATH/TO/DIRECTORY/scDimplot_level3_Legend_diff_cols_Legend.png",scDimplot, width = 15, height = 6, units = "cm")

# Figure S15A
scDimplot <- DimPlot(level3.seurat_L, split.by = "condition", cols = "alphabet")+ xlim(-8,8) + ylim(-5,5)
ggsave("./PATH/TO/DIRECTORY/scDimplot_level3_split_Legend_diff_cols.png",scDimplot, width = 25, height = 6, units = "cm")

scDimplot <- DimPlot(level3.seurat_L, split.by = "condition", cols = "alphabet") + NoLegend() + xlim(-8,8) + ylim(-5,5)
ggsave("./PATH/TO/DIRECTORY/scDimplot_level3_split_noLegend_diff_cols.png",scDimplot, width = 20, height = 6, units = "cm")

# Re-Organize Idents for the plot later (Default)
Idents(level3.seurat_L) <- factor(Idents(level3.seurat_L), levels = c("HAS1 High", "PLIN2+", "Myofibroblasts", "Fibroblasts"))

# ---
# Dot Plot
# Organize Condition factor for plot
level3.seurat_L$condition <- factor(level3.seurat_L$condition, levels = c("PF", "Control"))

# Define Gene list
MTHFD2_paper_h <- c("MTHFD2" ,"ALDH1L2", "MTHFD1L", "MTHFD1")

# Figure S15B
plot <- DotPlot(level3.seurat_L, features = MTHFD2_paper_h, group.by = "condition") + RotatedAxis()
W <- length(MTHFD2_paper_h) * 1.3
H <- length(unique(level3.seurat.split@active.ident)) * 0.6
ggsave("./PATH/TO/DIRECTORY/scDotPlot_level3_condition_MTHFD2_paper_TALL.png", plot, width = W, height = H)
H <- length(unique(level3.seurat.split@active.ident)) * 0.25
ggsave("./PATH/TO/DIRECTORY/scDotPlot_level3_condition_MTHFD2_paper_LOW.png", plot, width = W, height = H)

# Figure 5G
plot <- DotPlot(level3.seurat_L, features = MTHFD2_paper_h) + RotatedAxis()
W <- length(MTHFD2_paper_h) *1.3
H <- length(unique(level3.seurat_L@active.ident)) *1.1
ggsave("./PATH/TO/DIRECTORY/scDotPlot_Level3_MTHFD2_paper_RNA_Norm_TALL.png", width = W, height = H)
H <- length(unique(level3.seurat_L@active.ident)) *0.55
ggsave("./PATH/TO/DIRECTORY/scDotPlot_Level3_MTHFD2_paper_RNA_Norm_LOW.png", width = W, height = H)

# ---
# Feature Plot
# Figure 5F | ? 
for (gene in MTHFD2_paper_h) {
  plot <- FeaturePlot(level3.seurat_L, features = gene, order = T)
  ggsave(paste0("./PATH/TO/DIRECTORY/scFeaturePlot_Level3_",gene,"_RNA_Norm.png"), plot, width = 10, height = 6, units = "cm")
}

# Figure S15C
for (gene in MTHFD2_paper_h) {
  plot <- FeaturePlot_scCustom(level3.seurat_L, split.by = "condition", features = gene, order = T, colors_use = c("lightgrey", "blue"))
  ggsave(paste0("./PATH/TO/DIRECTORY/scFeaturePlot_level3_split_", gene,"_RNA_Norm.png"),
         width = 20,
         height = 6, units = "cm")
}
