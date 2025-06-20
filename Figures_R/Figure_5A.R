# Load libraries
library(Seurat)
library(ggplot2)

# Import pre-processed Seurat Object
TdTomato.seurat <- readRDS("/PATH/TO/DIRECTORY/GSM6428697_seurat_v4_obj.rds")

# Subset tdTomato+ Alveolar lineage
TdTomato.pos <- subset(TdTomato.seurat, subset = lineage == "tdTomato+")

# Combine Alveolar1 and Alveolar2 as Alveolar
TdTomato.pos <- RenameIdents(TdTomato.pos, Alveolar1 = "Alveolar", Alveolar2 = "Alveolar")

# subset only Alveolar, Fibrotic, and Inflammatory labeled cells and organize for later plots
TdTomato.fibroblast <- subset(TdTomato.pos, idents = c("Alveolar", "Fibrotic", "Inflammatory"))
Idents(TdTomato.fibroblast) <- factor(Idents(TdTomato.fibroblast), levels = c("Fibrotic", "Inflammatory", "Alveolar"))
TdTomato.fibroblast$Day <- factor(TdTomato.fibroblast$Day, levels = c("Day21", "Day14", "Day7", "Day0"))

# Normalize subset Seurat object
TdTomato.fibroblast <- NormalizeData(TdTomato.fibroblast)
TdTomato.fibroblast <- ScaleData(TdTomato.fibroblast)
TdTomato.fibroblast <- FindVariableFeatures(TdTomato.fibroblast)

# DimPlot Figure 5A
plot <- DimPlot(TdTomato.fibroblast) + ylim(-5, 8) + xlim(NA, 5)
ggsave("./PATH/TO/DIRECTORY/DimPlot_AllDays_Legend_cut.png", plot, width = 5, height = plot_height)

plot <- DimPlot(TdTomato.fibroblast) + NoLegend() + ylim(-5, 8) + xlim(NA, 5)
ggsave("./PATH/TO/DIRECTORY/DimPlot_AllDays_NoLegend_cut.png", plot, width = plot_width, height = plot_height)
