#install dplyr, patchwork and Seurat 
library(dplyr)
library(patchwork)
library(Seurat)
DataM3 <-Read10X(data.dir = "C:/Users/lhd1g21/OneDrive - University of Southampton/Dissertation/DataM3")
# Example work for Metastatic data (DataM3)

Meta3 <- CreateSeuratObject(counts = DataM3, project = "Meta3", min.cells = 3, min.features = 200)
Meta3[["percent.mt"]] <- PercentageFeatureSet(Meta3, pattern= "^MT-")
VlnPlot(Meta3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=5)
plot1 <- FeatureScatter(Meta3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Meta3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2
plot1
Meta3 <- subset(Meta3, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
Meta3 <- NormalizeData(Meta3, normalization.method = "LogNormalize", scale.factor = 10000)
Meta3 <- PercentageFeatureSet(Meta3, pattern = "^MT-", col.name = "percent.mt")
Meta3 <- FindVariableFeatures(Meta3, selection.method = "vst", nfeatures = 2000, loess.span = 1)
top10 <- head(VariableFeatures(Meta3), 10)
plot1 <- VariableFeaturePlot(Meta3)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot2
all.genes <- rownames(Meta3)
Meta3 <- ScaleData(Meta3, features = all.genes)
Meta3 <- RunPCA(Meta3, npcs = 20, features = VariableFeatures(object = Meta3))
print(Meta3[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Meta3, dims = 1:2, reduction = "pca")
DimPlot(Meta3, reduction = "pca") + NoLegend()
DimHeatmap(Meta3, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(Meta3)
DimHeatmap(Meta3, dims = 1:10, cells = 500, balanced = TRUE)
Meta3 <- FindNeighbors(Meta3, dims = 1:10)
Meta3 <- FindClusters(Meta3, resolution = 0.1)
Meta3 <- FindClusters(Meta3, resolution = 0.2)
Meta3 <- FindClusters(Meta3, resolution = 0.4)
Meta3 <- FindClusters(Meta3, resolution = 0.6)
Meta3 <- FindClusters(Meta3, resolution = 0.8)
library(clustree)
clustree(Meta3)
head(Idents(Meta3), 5)
Meta3 <- RunUMAP(Meta3, dims = 1:10)
DimPlot(Meta3, reduction = "umap")

BiocManager::install('celldex')
BiocManager::install('SingleR')
BiocManager::install('pheatmap')
BiocManager::install('tidyverse')
library(celldex)
library(SingleR)
library(pheatmap)
library(tidyverse)
ref <- celldex::HumanMeta3aryCellAtlasData()
Meta3_counts <- GetAssayData(Meta3, slot = 'counts')
pred <- SingleR(test = Meta3_counts,
        ref = ref,
        labels = ref$label.main)
Meta3$singleR.labels <- pred$labels[match(rownames(Meta3@Meta3.data),rownames(pred))]
DimPlot(Meta3, reduction = 'umap', group.by = 'singleR.labels')
plotScoreHeatmap(pred)
plotDeltaDistribution(pred)
tab <- table(Assigned=pred$labels, Clusters=Meta3$seurat_clusters)
pheatmap(log10(tab+10), color = colorRampPalette(c('white','blue'))(10))

saveRDS(Meta3, file = "C:/Users/leiah/OneDrive - University of Southampton/Dissertation/DataM3.rds")

cluster0.markers <- FindMarkers(Meta3, ident.1 = 0)
head(cluster0.markers, n = 5)
cluster2.markers <- FindMarkers(Meta3, ident.1 = 2)
head(cluster2.markers, n = 5)
cluster1.markers <- FindMarkers(Meta3, ident.1 = 1)
head(cluster1.markers, n = 5)
cluster3.markers <- FindMarkers(Meta3, ident.1 = 3)
head(cluster3.markers, n = 5)
cluster4.markers <- FindMarkers(Meta3, ident.1 = 4)
head(cluster4.markers, n = 5)
cluster5.markers <- FindMarkers(Meta3, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)
cluster6.markers <- FindMarkers(Meta3, ident.1 = 6)
head(cluster6.markers, n = 5)
cluster7.markers <- FindMarkers(Meta3, ident.1 = 7)
head(cluster7.markers, n = 5)
cluster8.markers <- FindMarkers(Meta3, ident.1 = 8)
head(cluster8.markers, n = 5)
cluster9.markers <- FindMarkers(Meta3, ident.1 = 9)
head(cluster9.markers, n = 5)
cluster10.markers <- FindMarkers(Meta3, ident.1 = 10)
head(cluster10.markers, n = 5)
cluster11.markers <- FindMarkers(Meta3, ident.1 = 11)
head(cluster11.markers, n = 5)
cluster12.markers <- FindMarkers(Meta3, ident.1 = 12)
head(cluster12.markers, n = 5)
Meta3.markers <- FindAllMarkers(Meta3, only.pos = TRUE)
Meta3.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
cluster0.markers <- FindMarkers(Meta3, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
Meta3.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 6) %>%
  ungroup() -> top6
DoHeatmap(Meta3, features = top6$gene) + NoLegend()

VlnPlot(Meta3, features = c("GNA11", "GNAQ"), y.max = 4)
VlnPlot(Meta3, features = c("EIF1AX", "SF3B1"), y.max = 4)
VlnPlot(Meta3, features = c("JUN", "JUNB"))
VlnPlot(Meta3, features = c("PRAME", "MDM2", "BAP1"))
VlnPlot(Meta3, features = c("CD44", "CD74", "CD47", "ALCAM", "NES", "PROM1", "MET"), y.max = 6.5)
VlnPlot(Meta3, features = c("MYC", "ARNT", "MET"))
VlnPlot(Meta3, features = c("ITGA4", "MMP9", "VEGFA", "ID3"))
VlnPlot(Meta3, features = c("SERTAD1", "CDK4", "IER2", "PLEKHB2","PARK2", "STAT3"))
VlnPlot(Meta3, features = c("FAS", "MMP3", "MMP7","CCR7", "IGF1", "CXCR4", "CXCL12"))
VlnPlot(Meta3, features = c("BIRC3", "IL1B", "BCL2A1","CXCL8", "TNFAIP3", "ICAM1", "PLAU"))
VlnPlot(Meta3, features = c("S100A6", "CITED6", "PAEP","QPCT", "CCL4", "STAB1", "PDE4B"))
VlnPlot(Meta3, features = c("RELA", "KRAS", "NRAS", "HRAS"))
VlnPlot(Meta3, features = c("ASAP1", "ARF6", "HGF","IGF1", "NFAT1"))
VlnPlot(Meta3, features = c("HIF1A", "TP53", "BMF","BAD", "CXCR2"))
VlnPlot(Meta3, features = c("NKG7", "PMEL"))
VlnPlot(Meta3, features = c("TYR", "PMEL", "MLANA"))
VlnPlot(Meta3, features = c("YAP1", "TNF"))
VlnPlot(Meta3, features = c("SPP1"))
VlnPlot(Meta3, features = c("PIK3CA"))
VlnPlot(Meta3, features = c("BRAF", "NRAS", "KIT", "NF1"))
VlnPlot(Prim, features = c("CD163", "CD68"))

FeaturePlot(Meta3, features = c("BIRC3", "IL1B", "BCL2A1","CXCL8", "TNFAIP3", "ICAM1", "PLAU"))
FeaturePlot(Meta3, features = c("XIAP", "apollon", "ML-IAP","Survivin", "Naip", "BIRC2", "BIRC8"))
FeaturePlot(Meta3, features = c("S100A6", "CITED6", "PAEP","QPCT", "CCL4", "STAB1", "PDE4B"))
FeaturePlot(Meta3, features = c("TYR", "PMEL", "MLANA"))
FeaturePlot(Meta3, features = c("MET", "NES", "ALCAM", "PROM1", "MCAM", "FAS")) 
FeaturePlot(Meta3, features = c("ITGA4", "MMP9", "VEGFA", "ID3"))
FeaturePlot(Meta3, features = c("CD44", "CD74", "CD47", "CCR7"))
FeaturePlot(Meta3, features = c("HLA-DRB1", "HLA-DQA1", "HLA-DPA1"))
FeaturePlot(Meta3, features = c("GNAQ", "GNA11"))
FeaturePlot(Meta3, features = c("YAP1", "PIK3CA", "TNF"))
FeaturePlot(Meta3, features = c("SPP1", "PRAME","MDM2"))
FeaturePlot(Meta3, features = c("HLA-A", "HLA-B","HLA-C"))
FeaturePlot(Meta3, features = c("SERTAD1", "ALDH1A1", "KDM5B", "CDK4"))
FeaturePlot(Meta3, features = c("RELA", "KRAS", "NRAS", "HRAS"))
FeaturePlot(Meta3, features = c("ASAP1", "ARF6", "HGF", "IGF1"))
FeaturePlot(Meta3, features = c("HIF1A", "TP53", "BMF","BAD", "CXCR2"))
FeaturePlot(Meta3, features = c("BAP1", "PLAU", "PLAUR", "PMEL"))
FeaturePlot(Meta3, features = c("SERTAD1", "IER2", "STAT3", "PLEKHB2"))
FeaturePlot(Meta3, features = c("SERTAD1", "ICE2", "MAFF", "TADA1"))
FeaturePlot(Prim, features = c("CD163", "CD68"))

# Chromosome plot 
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75
gene_symbols <- rownames(Meta3)
genes <- genes(edb, filter = GeneNameFilter(gene_symbols))
genes <- as.data.frame(genes)
chromosome_info <- data.frame(
  gene_symbol = genes$gene_name,
  chromosome = genes$seqnames)
chromosome_vector <- setNames(chromosome_info$chromosome, chromosome_info$gene_symbol)
expression_data <- GetAssayData(Meta3, slot = "data")
gene_chromosome_map <- data.frame(
  gene = rownames(expression_data),
  chromosome = chromosome_vector[rownames(expression_data)],
  stringsAsFactors = FALSE)
expression_and_chromosome <- expression_data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  inner_join(gene_chromosome_map, by = "gene")
average_expression_by_chromosome <- expression_and_chromosome %>%
  group_by(chromosome) %>%
  summarise(across(-gene, sum))
average_expression_by_chromosome_mtx <- as.data.frame(average_expression_by_chromosome)
average_expression_by_chromosome_mtx <- average_expression_by_chromosome_mtx[-72,]
rownames(average_expression_by_chromosome_mtx) <- average_expression_by_chromosome_mtx[,1]
average_expression_by_chromosome_mtx <- average_expression_by_chromosome_mtx[,-1]
column_sums <- colSums(average_expression_by_chromosome_mtx)
normalized_matrix <- sweep(average_expression_by_chromosome_mtx, 2, column_sums, `/`)
library(reshape2)
rows_to_keep <- !grepl("PATCH", rownames(normalized_matrix))
normalized_matrix <- normalized_matrix[rows_to_keep, ]
rows_to_keep <- !grepl("GL", rownames(normalized_matrix))
normalized_matrix <- normalized_matrix[rows_to_keep, ]
long_data <- as.data.frame(normalized_matrix)
long_data$Chromosome <- rownames(normalized_matrix)
long_data <- reshape2::melt(long_data, id.vars = "Chromosome")
names(long_data) <- c("Chromosome", "Cell", "Value")
selected_chromosomes_data <- long_data[long_data$Chromosome %in% c("1", "3"),]
plot<-ggplot(selected_chromosomes_data, aes(x = Value, fill = Chromosome)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +  # Adjust alpha for transparency
  scale_fill_manual(values = c("blue", "red")) +  # Manually set colors for each chromosome
  labs(title = "Distribution of Reads Mapping to Chromosome 1 and 3",
       x = "Proportion of Reads",
       y = "Count") +
  theme_minimal()







