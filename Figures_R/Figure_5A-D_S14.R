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

# Subset the Seurat object based on the Day metadata
TdTomato.D0 <- subset(TdTomato.fibroblast, subset = Day == "Day0")
TdTomato.D7 <- subset(TdTomato.fibroblast, subset = Day == "Day7")
TdTomato.D14 <- subset(TdTomato.fibroblast, subset = Day == "Day14")
TdTomato.D21 <- subset(TdTomato.fibroblast, subset = Day == "Day21")

days_list <- list(TdTomato.D0, TdTomato.D7, TdTomato.D14, TdTomato.D21)
days_labels <- c("Day0", "Day7", "Day14", "Day21")

# ---------------
# Dim Plots
# Visualize UMAP for each subset and save the plots
for (i in seq_along(days_list)) {
  days_list[[i]] <- NormalizeData(days_list[[i]])
  days_list[[i]] <- ScaleData(days_list[[i]])
  days_list[[i]] <- FindVariableFeatures(days_list[[i]])
  scDimPlot <- DimPlot(days_list[[i]])
  ggsave(paste0("./PATH/TO/DIRECTORY/scDimPlot_", days_labels[i], ".png"), scDimPlot)
}

# Collect the subsets into a list for iteration
days <- list(Day0 = TdTomato.D0, Day7 = TdTomato.D7, Day14 = TdTomato.D14, Day21 = TdTomato.D21)

# Set plot dimensions
plot_width <- 3.5
plot_height <- 3.5

# Loop over each Seurat objects
for (day in names(days)) {
  plot <- DimPlot(days[[day]]) + NoLegend() + ylim(-5, 8) + xlim(NA, 5)
  ggsave(paste0("./PATH/TO/DIRECTORY/DimPlot_", day, "_cut.png"), plot, width = plot_width, height = plot_height)
}

# Set uniform x-axis tic
plot <- DimPlot(TdTomato.D7) + NoLegend() + ylim(-5, 8) + xlim(NA, 5) + scale_x_continuous(breaks = c(-2.5, 0.0, 2.5, 5.0))
ggsave("./PATH/TO/DIRECTORY/DimPlot_Day7_cut.png", plot, width = plot_width, height = plot_height)

plot <- DimPlot(TdTomato.D0) + NoLegend() + ylim(-5, 8) + xlim(NA, 5.0) + scale_x_continuous(breaks = c(-2.5, 0.0, 2.5, 5.0)) + expand_limits(x = 5.0)
ggsave("./PATH/TO/DIRECTORY/DimPlot_Day0_cut.png", plot, width = plot_width, height = plot_height)

# ---------------
# Feature Plots
# Create a genes list
MTHFD2_paper_m <- c("Mthfd2" ,"Aldh1l2", "Mthfd1l", "Mthfd1")

# Set plot dimensions
plot_width <- 3.5
plot_height <- 3.5

for (day in names(days)) {
  for (gene in MTHFD2_paper_m) {
    plot <- FeaturePlot(days[[day]], features = gene, order = T) + ylim(-5, 8) + xlim(NA, 5)
  ggsave(paste0("./PATH/TO/DIRECTORY/scFeaturePlot_", day, "_" ,gene,"_cut.png"), plot, width = plot_width, height = plot_height)
  }
}

# Set uniform x-axis tic
for (gene in MTHFD2_paper_m) {
  plot <- FeaturePlot(TdTomato.D7, features = gene, order = T) + ylim(-5, 8) + xlim(NA, 5.0)
ggsave(paste0("./PATH/TO/DIRECTORY/scFeaturePlot_Day7_",gene,"_cut.png"), plot, width = plot_width, height = plot_height)
}

for (gene in MTHFD2_paper_m) {
  plot <- FeaturePlot(TdTomato.D0, features = gene, order = T) + ylim(-5, 8) + xlim(NA, 5.0) + scale_x_continuous(breaks = c(-2.5, 0.0, 2.5, 5.0)) + expand_limits(x = 5.0)
ggsave(paste0("./PATH/TO/DIRECTORY/scFeaturePlot_Day0_",gene,"_cut.png"), plot, width = plot_width, height = plot_height)
}

for (gene in MTHFD2_paper_m) {
  plot <- FeaturePlot(TdTomato.fibroblast, features = gene, order = T) + ylim(-5, 8) + xlim(NA, 5)
  ggsave(paste0("./PATH/TO/DIRECTORY/scFeaturePlot_AllDays_", gene, "_cut.png"), plot, width = plot_width, height = plot_height)
}

# ---------
# Dot Plot
# Figure 5C
plot <- DotPlot(TdTomato.fibroblast, features = MTHFD2_paper_m, group.by = "Day") + RotatedAxis()
W <- length(MTHFD2_paper_m) * 1
H <- length(unique(TdTomato.fibroblast@active.ident)) * 1.2
ggsave("./PATH/TO/DIRECTORY/scDotPlot_TdTomato_MTHFD2_paper_Days_TALL.png", plot, width = W, height = H)
H <- length(unique(TdTomato.fibroblast@active.ident)) * 0.6
ggsave("./PATH/TO/DIRECTORY/scDotPlot_TdTomato_MTHFD2_paper_Days_LOW.png", plot, width = W, height = H)

# Figure 5D
scDotplot <- DotPlot(TdTomato.fibroblast, features = MTHFD2_paper_m) + RotatedAxis()
H <- length(unique(TdTomato.fibroblast@active.ident)) *1.4
ggsave("./PATH/TO/DIRECTORY/scDotPlot_Tsukui_MTHFD2_TALL.png", scDotplot, width = 4.5, height = H)
H <- length(unique(level4.seurat_3@active.ident)) *0.7
ggsave("./PATH/TO/DIRECTORY/scDotPlot_Tsukui_MTHFD2_LOW.png", scDotplot, width = 4.5, height = H)

# Figure S14B
for (day in names(days)) {
  scDotPlot <- DotPlot(days[[day]], features = MTHFD2_paper_m) + RotatedAxis()
  H <- length(unique(TdTomato.fibroblast@active.ident)) * 1.2
  ggsave(paste0("./PATH/TO/DIRECTORY/DotPlots_TdTomato/scDotPlot_",day,"_MTHFD2_paper_TALL.png"),scDotPlot,width = W,height = H)
  H <- length(unique(TdTomato.fibroblast@active.ident)) * 0.6
  ggsave(paste0("./PATH/TO/DIRECTORY/DotPlots_TdTomato/scDotPlot_",day,"_MTHFD2_paper_LOW.png"),scDotPlot,width = W,height = H)
}
