# Load libraries
library(Seurat)
library(ggplot2)

# Load Data
level4.seurat_L <- readRDS(file = "./PATH/TO/DIRECTORY/level4.seurat_L.rds") # SCT transformed preprocessed fibroblast (MyoFibroblast, Fibroblast except PI16 high)

# Set default assay as "RNA" rds file default assay is SCT
DefaultAssay(level4.seurat_L) <- "RNA"

# Normalize RNA assay
level4.seurat_L <- NormalizeData(level4.seurat_L)
level4.seurat_L <- ScaleData(level4.seurat_L)
level4.seurat_L <- FindVariableFeatures(level4.seurat_L)

# Organize Idents for plot later (Default)
Idents(level4.seurat_L) <- factor(Idents(level4.seurat_L), levels = c("Fibrotic", "Inflammatory", "Alveolar"))

# Create a genes list
MTHFD2_paper_h <- c("MTHFD2" ,"ALDH1L2", "MTHFD1L", "MTHFD1")

# Figure 5J
scDotplot <- DotPlot(level4.seurat_L, features = MTHFD2_paper_h) + RotatedAxis()
W <- length(MTHFD2_paper_h) * 1.2
H <- length(unique(level4.seurat_L@active.ident)) *1.4
ggsave("./PATH/TO/DIRECTORY/scDotPlot_Level4_new_MTHFD2_RNA_TALL.png", scDotplot, width = W, height = H)
H <- length(unique(level4.seurat_L@active.ident)) *0.7
ggsave("./PATH/TO/DIRECTORY/scDotPlot_Level4_new_MTHFD2_RNA_LOW.png", scDotplot, width = W, height = H)

# ---------
# Dim Plot
# Figure 5H
dimplot <- DimPlot(level4.seurat_L, label = F, cols = "alphabet2") + NoLegend()
ggsave("./PATH/TO/DIRECTORY/scDimPlot_Level4_new_noLegend_diff_cols.png", dimplot, width = 3.5, height = 3.5)
dimplot <- DimPlot(level4.seurat_L, label = F)
ggsave("./PATH/TO/DIRECTORY/scDimPlot_Level4_new_diff_cols.png", dimplot, width = 5, height = 3.5)

# ---------
# Feature Plot
# Figure 5I / Figure S18B
DefaultAssay(level4.seurat_L) <- "RNA"
for (gene in MTHFD2_paper_h) {
  featureplots <- FeaturePlot(level4.seurat_L, features = gene, order = T)
  ggsave(paste0("./PATH/TO/DIRECTORY/scFeaturePlots_Level4_paper_RNA_", gene,".png"), featureplots, width = 3.5, height = 3.5)
}

# Figure S18A
# Create a marker list
markers <- c("SCN7A", "AOC3", "FMO2", "CCL2", "SFRP2", "GSN", "CTHRC1", "POSTN", "COL1A1")

for (gene in markers) {
  featureplots <- FeaturePlot(level4.seurat_L, features = gene, order = T)
  ggsave(paste0("./PATH/TO/DIRECTORY/scFeaturePlots_Level4_paper_SCT_", gene,".png"), featureplots, width = 3.5, height = 3.5)
}
