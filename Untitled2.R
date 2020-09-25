library(ggplot2)
library(cowplot)
library(patchwork)
BiocManager::install("devtools")
library(devtools)
install_github("immunogenomics/harmony")
library(harmony)
library(Seurat)
BiocManager::install("SeuratData")
install_github('satijalab/seurat-data')
library(SeuratData)

InstallData("panc8")
data("panc8")

# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(pancreas.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype", label = TRUE, 
              repel = TRUE) + NoLegend()
p1 + p2



InstallData("pbmcsca")
data("pbmcsca")

withharmony <- NormalizeData(pbmcsca) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
withharmony <- RunHarmony(withharmony, group.by.vars = "Method")
withharmony <- RunUMAP(withharmony, reduction = "harmony", dims = 1:30)
withharmony <- RunTSNE(withharmony, reduction = "harmony", dims = 1:30)
withharmony <- FindNeighbors(withharmony, reduction = "harmony", dims = 1:30) %>% FindClusters()
p1 = DimPlot(withharmony, group.by = c("Method", "ident", "CellType"), ncol = 3)
p2 = DimPlot(withharmony, reduction = "harmony", group.by = c("Method", "ident", "CellType"), ncol = 3)
p3 = DimPlot(withharmony, reduction = "tsne", group.by = c("Method", "ident", "CellType"), ncol = 3)

withoutharmony <- NormalizeData(pbmcsca) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
withoutharmony <- RunUMAP(withoutharmony, reduction = "pca", dims = 1:30)
withoutharmony <- RunTSNE(withoutharmony, reduction = "pca", dims = 1:30)
p4 <- DimPlot(withoutharmony, reduction = "umap", group.by = c("Method", "ident", "CellType"), ncol = 3)
p5 <- DimPlot(withoutharmony, reduction = "tsne", group.by = c("Method", "ident", "CellType"), ncol = 3)

library(gridExtra)
pdf(file="TRAVAIL/relecture_publication/TEST_harmony.pdf")
print(grid.arrange(p1,p4,nrow=2))
print(grid.arrange(p3,p5,nrow=2))
dev.off()