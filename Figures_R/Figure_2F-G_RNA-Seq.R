# Load Libraries
library(readr)
library(readxl)
library(stringr)
library(dplyr)
library(data.table)
library(tidyverse)
library(dbplyr)
library(RColorBrewer)
library(data.table)
library(ggplot2)
library(BiocManager)
library(Biobase)
library(tximport)
library(tximportData)
library(S4Arrays)
library(DelayedArray)
library(tximeta)
library(Rhtslib)
library(AnnotationDbi)
library(GenomicFeatures)
library(GenomeInfoDb)
library(GenomicRanges)
library(BiocGenerics)
library(parallel)
library(S4Vectors)
library(Rsamtools)
library(stats4)
library(IRanges)
library(limma)
library(edgeR)
library(DESeq2)
library(pheatmap)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# ------------
# RNA Count
# Set Directory
dir <- system.file("extdata", package = "tximportData")
list.files(dir)

# Transcript to gene table
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
tx2gene <- na.omit(tx2gene)
tx2gene

txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)

head(txi$counts)
head(txi$abundance)
head(txi$countsFromAbundance)

txi.tx <- tximport(files, type = "kallisto", txOut = TRUE)
txi.sum <- summarizeToGene(txi.tx, tx2gene)
#transcripts missing from tx2gene: 201721
all.equal(txi$counts, txi.sum$counts)

attributes(txi)

# Write the counts to an object
All_siMTHFD2_data <- txi$counts %>% round() %>% data.frame()
head(All_siMTHFD2_data)

duplicates <- duplicated(colnames(All_siMTHFD2_data))
if(any(duplicates)) {
  print("There are duplicate column names.")
  print(colnames(All_siMTHFD2_data)[duplicates])  # Prints the names of duplicate columns
} else {
  print("All column names are unique.")
}

# Write the abundance to an object
read_abundance <- data.frame(txi$abundance)

duplicates <- duplicated(rownames(read_abundance))
if(any(duplicates)) {
  print("There are duplicate row names.")
  print(colnames(read_abundance)[duplicates])  # Prints the names of duplicate columns
} else {
  print("All row names are unique.")
}

oldcolnames <- colnames(read_abundance)
newcolnames <- samples$Name

setnames(read_abundance, old = oldcolnames, new = newcolnames)

read_abundance <- rownames_to_column(read_abundance, var = "entrezgene_id")
head(read_abundance)

# Load biomaRt
library(biomaRt)

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# for human dataset = "hsapiens_gene_ensembl"

id_map <- getBM(attributes=c('ensembl_gene_id',
                          'external_gene_name','entrezgene_id'),
             filters = 'entrezgene_id',
             values = read_abundance$entrezgene_id,
             mart = ensembl)


read_abundance <- read_abundance %>%  mutate(entrezgene_id = as.numeric(entrezgene_id))

read_abundance <- left_join(id_map, read_abundance, by = "entrezgene_id")

read_abundance <- read_abundance %>% rename("external_gene_name"="Symbol")

head(read_abundance)
rep_Symbol <- dim(read_abundance)[1] - length(unique(read_abundance$Symbol))
rep_Symbol #3856

read_abundance <- read_abundance[!duplicated(read_abundance$Symbol), ] #21157 left

write_csv(read_abundance, "./PATH/TO/DIRECTORY/siMTHFD2_abundance.csv")
read_abundance[newcolnames] <- sapply(read_abundance[newcolnames],as.numeric)#set data columns as numeric

sapply(read_abundance, class)
write_csv(read_abundance, "./PATH/TO/DIRECTORY/siMTHFD2_abundance.csv")
head(read_abundance)

TableOfCounts <- read_abundance[-c(1,3)] #or TableOfCounts <- select(read_abundance,-1,-3)
head(TableOfCounts)
write_csv(TableOfCounts, "./PATH/TO/DIRECTORY/TableOfCounts_siMTHFD2.txt")

# check if there is any NA
sum(is.na(TableOfCounts))

# ------------
# EdgeR
# read/use TableOfCounts
ALL_tpm <- TableOfCounts[,-1 ]  # Exclude the "Name" column
rownames(ALL_tpm) <- TableOfCounts$Symbol
colnames(ALL_tpm)

# DefineGroups 
group <- factor(c(samples$Group))

# Differentially Expressed Genes
DEG <- DGEList(counts=ALL_tpm, group=group)
keep <- filterByExpr(DEG)
DEG <- DEG[keep,,keep.lib.sizes=FALSE] 
DEG <- calcNormFactors(DEG)

# group design matrix
design <- model.matrix(~0+group)
design
Groupnames <- unique(c(samples$Treatment))
colnames(design) <- Groupnames
design

DEG <- estimateDisp(DEG,design) # sometimes DEG is called as NormData
fit <- glmQLFit(DEG,design)

# UT vs TGFb
qlf.siContr_TGF<- glmQLFTest(fit, contrast=c(-1,1,0,0,0,0))
topTags(qlf.siContr_TGF) #Coefficient: Coefficient:  -1*siContUT 1*siContTGF 
summary(dt<-decideTestsDGE(qlf.siContr_TGF, p.value = 0.01))

#       -1*siContUT 1*siContTGF
#Down              1426
#NotSig            6147
#Up                1231

write.csv(qlf.siContr_TGF, "./PATH/TO/DIRECTORY/CtrUT_vs_CtrTGFB.csv")

# siControl TGFb vs siMTHFD2 TGFb
qlf.ctrTGF_siTGF<- glmQLFTest(fit, contrast=c(0,-1,0,1))
topTags(qlf.ctrTGF_siTGF) #Coefficient: -1*Ctr_TGF 1*siKD_TGF 
summary(dt<-decideTestsDGE(qlf.ctrTGF_siTGF, p.value = 0.05))

#       -1*siContUT 1*siContTGF
#Down             1
#NotSig            8803
#Up              0

write.csv(qlf.ctrTGF_siTGF, "./PATH/TO/DIRECTORY/CtrTGF_vs_siTGFB.csv") 

# ---------
# Volcano Plot
# Load data
Ctr_TGF_vs_siMTHFD2_TGF <- read.csv("./PATH/TO/DIRECTORY/CtrTGF_vs_siTGFB.csv", row.names=1)

# Figure 2G
VP_Ctr_TGF_vs_siMTHFD2_TGF  <- EnhancedVolcano(Ctr_TGF_vs_siMTHFD2_TGF ,
    lab = rownames(Ctr_TGF_vs_siMTHFD2_TGF ),
    x = 'logFC',
    y = 'PValue',
    pCutoff = 0.01,
    FCcutoff =1.0,
    pointSize = 4.0,
    labSize = 16.0,
    col=c('black', 'blue', 'grey', 'red3'),
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    title = 'Ctr_TGF_vs_siMTHFD2_TGFβ', 
    axisLabSize = 24, 
    xlim = c(-10, 10),  # Set X-axis limits
    ylim = c(0.0001, 25))

# Set the desired width and height for the plot (in inches)
plot_width <- 20
plot_height <- 20

# Save the plot using ggsave
ggsave("./PATH/TO/DIRECTORY/VP_Ctr_TGF_vs_siMTHFD2_TGFB.png", plot = VP_Ctr_TGF_vs_siMTHFD2_TGF, width = plot_width, height = plot_height, limitsize = FALSE)

# ---------
# Heatmap
counts.in <- TableOfCounts[,-1]
rownames(counts.in) <- TableOfCounts$Symbol

head(counts.in)
colnames(counts.in)

# Define metadata
counts.metadata<- data.frame(
  dataset= c (colnames(counts.in)), 
  Treatment = c(Samples_siMTHFD2$Treatment),
   stringsAsFactors = FALSE)

# EdgeR Zscore Calc
group_HM <- counts.metadata$Treatment ## sample #
y_HM <- DGEList(counts=counts.in, genes=row.names.data.frame(counts.in), group=group_HM)
# Low count removal
keep_HM <- rowSums(cpm(y_HM)>1) >= 1
table(keep_HM)
y_HM <- y_HM[keep_HM, , keep.lib.sizes=FALSE]
y_HM <- calcNormFactors(y_HM, method="TMM")
counts_HM <- as.matrix(y_HM$counts)
logCPM_HM <- cpm(counts_HM, prior.count=1, log=TRUE) #normalization using cpm here
ZScore_HM <- t(scale(t(logCPM_HM)))
ZScore_HM <- as.data.frame(ZScore_HM)

# Define significant genes
Sig_HM_all  <- subset.data.frame(CtrUT_vs_CtrTGFB ,PValue<0.01)
Sig_HM_all  <- rownames_to_column(Sig_HM_all, "Symbol")

ctrSymbol_all <- c(Sig_HM_all$Symbol)

# Keep Only significant genes
sig.zscore <- ZScore_HM[ctrSymbol_all  ,]#change for each cluster
sig.zscore.mat <- as.matrix(sig.zscore)

sig.zscore.mat <- sig.zscore.mat[complete.cases(sig.zscore.mat),]

# Heatmap Annotation
heat.annotation <- data.frame(counts.metadata[ ,2])
colnames(heat.annotation) <- "Treatment"
row.names(heat.annotation) <- counts.metadata[ ,1]
nrow(sig.zscore.mat)

d.colnames <- c(counts.metadata[,1])
colnames(sig.zscore.mat) <- d.colnames

# Assign Color for each group
ann.colors <- list(Treatment= c(`siContUT` = "pink", `siCtr_TGFβ` = "red",`si9UT` = "gold", `si9_TGFβ` = "brown",`si11UT` = "cyan2", `si11_TGFβ` = "blue" ))

newnames <- lapply(
  rownames(sig.zscore.mat),
  function(x) bquote(italic(.(x)))) ## italicizes

# Figure 2F
heatmap_all <- pheatmap(sig.zscore.mat,
                         annotation_col = heat.annotation, 
                         cluster_cols = FALSE,
                         cluster_rows = TRUE,  # Perform row clustering
                         #treeheight_row = 0,  # Hide row dendrogram
                         main = "",
                         annotation_colors = ann.colors,
                         show_colnames = FALSE,
                         show_rownames = FALSE,,  # Keep row labels
                         fontsize_row = 6
)

# Define Plot Dimension
plot_width <- 20
plot_height <- 20

# Save the plot using ggsave
ggsave("./PATH/TO/DIRECTORY/heatmap_all.png", plot = heatmap_all, width = plot_width, height = plot_height, limitsize = FALSE)
