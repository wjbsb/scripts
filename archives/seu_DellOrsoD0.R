#!/usr/bin/env Rscript
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)

source(file="~/INMG_SingleCell/scripts/functions_stock.R",local=T)

prloc="~/INMG_SingleCell/"
setwd(prloc)

resu = "results/DellOrsoD0/"
dir.create(resu,recursive = T)
rdsdir = "rds/DellOrsoD0/"
dir.create(rdsdir,recursive = T)


dorso1.data <- Read10X(data.dir=paste0("data/DellOrsoD0/dorsowt1"))
dorso2.data <- Read10X(data.dir=paste0("data/DellOrsoD0/dorsowt2"))

dorso1 <- CreateSeuratObject(dorso1.data, project="DellOrso", min.cells=3, min.features=200)
dorso2 <- CreateSeuratObject(dorso2.data, project="DellOrso", min.cells=3, min.features=200)
rm(dorso1.data)
rm(dorso2.data)

dorso <- merge(dorso1,y=dorso2, add.cells.ids=c("wt1","wt2"),project="DellOrso")
rm(dorso1)
rm(dorso2)

# first steps seurat until PCA
dorso <- NormFeatScalePCA(dorso,2500,5) # function from functions_stock.R

pdf(paste0(resu,"preprocessSeu.pdf"))
VlnPlot(dorso, features =c("nFeature_RNA","nCount_RNA","percent.mt"), 
        ncol = 3) 
FeatureScatter(dorso, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
DimHeatmap(dorso, dims=1:15, cells = 500, balanced = TRUE) 
ElbowPlot(dorso)
dev.off()

# clustering and non-linear reduction for visuals
dorso <- KNNplusFITSNE(dorso, 15, 0.5)
# howmanyclusters we have found out:
nb.clus = max(as.integer(levels(dorso@meta.data$seurat_clusters)))+1
mycols <- colorRampPalette(brewer.pal(8, "Set2"))(nb.clus)

dorso.markers <- FindAllMarkers(dorso, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


topn <- dorso.markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_logFC)

# print results (plots)
pdf(paste0(resu,"FITSNE.pdf"))
DimPlot(dorso, reduction="tsne", label=TRUE, pt.size=0.5, 
        cols=mycols, shape.by="orig.ident") + ggtitle("DellOrso D0")
dev.off()

pdf(paste0(resu,"HEATMAP.pdf"),width=13)
DoHeatmap(dorso, features = topn$gene, group.colors=mycols) + NoLegend()
dev.off()

# save .rds object
saveRDS(dorso,file=paste0(rdsdir,"dorso_seu_fitsne.rds"))
#write table of all markers
write.table(dmizero.markers, paste0(resu,"ALLMARKERS_DellOrsoD0.txt"))
print("finished")
print("finished")




