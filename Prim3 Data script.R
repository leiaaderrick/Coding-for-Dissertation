library(dplyr)
library(patchwork)
library(Seurat)
DataP2 <-Read10X(data.dir = "C:/Users/leiah/OneDrive - University of Southampton/Dissertation/DataP2")
DataP2 <-Read10X(data.dir = "C:/Users/lhd1g21/OneDrive - University of Southampton/Dissertation/DataP2")
# work for Prim2ary data (DataP2)

Prim2 <- CreateSeuratObject(counts = DataP2, project = "Prim2", min.cells = 3, min.features = 200)
Prim2[["percent.mt"]] <- PercentageFeatureSet(Prim2, pattern= "^MT-")
VlnPlot(Prim2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=5)
plot1 <- FeatureScatter(Prim2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Prim2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2
plot1
Prim2 <- subset(Prim2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 15)
Prim2 <- NormalizeData(Prim2, normalization.method = "LogNormalize", scale.factor = 10000)
Prim2 <- PercentageFeatureSet(Prim2, pattern = "^MT-", col.name = "percent.mt")
Prim2 <- FindVariableFeatures(Prim2, selection.method = "vst", nfeatures = 2000, loess.span = 1)
top10 <- head(VariableFeatures(Prim2), 10)
plot1 <- VariableFeaturePlot(Prim2)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot2
all.genes <- rownames(Prim2)
Prim2 <- ScaleData(Prim2, features = all.genes)
Prim2 <- RunPCA(Prim2, npcs = 20, features = VariableFeatures(object = Prim2))
print(Prim2[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Prim2, dims = 1:2, reduction = "pca")
DimPlot(Prim2, reduction = "pca") + NoLegend()
DimHeatmap(Prim2, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(Prim2, dims = 1:14, cells = 500, balanced = TRUE)
ElbowPlot(Prim2)
Prim2 <- FindNeighbors(Prim2, dims = 1:8)
Prim2 <- FindClusters(Prim2, resolution = 0.2, 0.3,0.4, 0.5, 0.6 0.7)
library(clustree)
tree <- clustree(Prim2)
tree
Prim2 <FindClusters(Prim2, resolution = 0.4)
head(Idents(Prim2), 5)
Prim2 <- RunUMAP(Prim2, dims = 1:8)
DimPlot(Prim2, reduction = "umap")
library(celldex)
library(SingleR)
library(pheatmap)
library(tidyverse)
ref <- celldex::HumanPrimaryCellAtlasData()
Prim2_counts <- GetAssayData(Prim2, slot = 'counts')
pred <- SingleR(test = Prim2_counts,
        ref = ref,
        labels = ref$label.main)
Prim2$singleR.labels <- pred$labels[match(rownames(Prim2@meta.data),rownames(pred))]
DimPlot(Prim2, reduction = 'umap', group.by = 'singleR.labels')
plotScoreHeatmap(pred)
plotDeltaDistribution(pred)
tab <- table(Assigned=pred$labels, Clusters=Prim2$seurat_clusters)
pheatmap(log10(tab+10), color = colorRampPalette(c('white','blue'))(10))

saveRDS(Prim2, file = "C:/Users/leiah/OneDrive - University of Southampton/Dissertation/DataP2.rds")
cluster2.markers <- FindMarkers(Prim2, ident.1 = 1)
head(cluster2.markers, n = 5)
cluster5.markers <- FindMarkers(Prim2, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)
Prim2.markers <- FindAllMarkers(Prim2, only.pos = TRUE)
Prim2.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
cluster0.markers <- FindMarkers(Prim2, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(Prim2, features = c("GNA11", "GNAQ"))
VlnPlot(Prim2, features = c("EIF1AX", "SF3B1"))
VlnPlot(Prim2, features = c("BAP1", "PRAME"))
VlnPlot(Prim2, features = c("MET", "CD44", "CD74", "CD47", "ALCAM", "NES", "PROM1"))
VlnPlot(Prim2, features = c("MYC", "JUN", "ARNT"))
VlnPlot(Prim2, features = c("ITGA4", "MMP9", "VEGFA", "ID3"))
VlnPlot(Prim2, features = c("FAS", "MMP3", "MMP7","CCR7", "IGF1", "CXCR4", "CXCL12"))
VlnPlot(Prim2, features = c("NKG7", "EIF1AX"), slot = "counts", log = TRUE)
FeaturePlot(Prim2, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))
FeaturePlot(Prim2, features = c("MET", "NES", "ALCAM", "PROM1", "MCAM", "FAS")) 
FeaturePlot(Prim2, features = c("ITGA4", "MMP9", "VEGFA", "ID3"))
FeaturePlot(Prim2, features = c("CD44", "CD74", "CD47", "CCR7"))
FeaturePlot(Prim2, features = c("GNAQ", "GNA11"))
FeaturePlot(Prim2, features = c("EIF1AX", "SF3B1"))
FeaturePlot(Prim2, features = c("BAP1", "PRAME", "PLAUR", "PMEL"))
Prim2.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(Prim2, features = top10$gene) + NoLegend()
new_cluster_ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet", "Additional1", "Additional2", "Additional3", "Additional4", "Additional15")
Prim2$seurat_clusters <- factor(Prim2$seurat_clusters, levels = levels(Prim2$seurat_clusters), labels = new_cluster_ids)
DimPlot(Prim2, reduction = "umap", group.by ="seurat_clusters",label = TRUE, pt.size = 0.5) + NoLegend()
library(ggplot2)
plot <- DimPlot(Prim2, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
plot





