#!/usr/bin/env Rscript
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)

source(file="~/INMG_SingleCell/scripts/functions_stock.R",local=T)

prloc="~/INMG_SingleCell/"
setwd(prloc)

resu = "results/DeMicheliD0/"
dir.create(resu,recursive = T)
rdsdir = "rds/DeMicheliD0/"
dir.create(rdsdir,recursive = T)


dmizerol <- read.table(paste0("data/DeMicheliD0/rawdataD0.txt"), sep="\t",
                       header=T, row.names=1)

# check tech replicates
unique(substring(colnames(dmizerol), first=1, last=5)) #"D0_A_" "D0_B_" "D0_Cv"
# remember GEO notes: "Cell 3' Reagent Kit (10X Genomics, v2 kit for all sample except D0_Cv3 with version 3) following the directions of the manufacturer"

dmizero <- CreateSeuratObject(dmizerol, project="DeMicheli", min.cells=3, min.features=200)
dmizero

# steps simple preproc until pca
dmizero <- NormFeatScalePCA(dmizero,2500,5) # function from functions_stock.R

pdf(paste0(resu,"preprocessSeu.pdf"))
VlnPlot(dmizero, features =c("nFeature_RNA","nCount_RNA","percent.mt"), 
        ncol = 3) 
FeatureScatter(dmizero, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
DimHeatmap(dmizero, dims=1:15, cells = 500, balanced = TRUE) 
ElbowPlot(dmizero)
dev.off()

# clustering and non-linear reduction for visuals
dmizero <- KNNplusFITSNE(dmizero, 15, 0.5)
# howmanyclusters we have found out:
nb.clus = max(as.integer(levels(dmizero@meta.data$seurat_clusters)))+1
mycols <- colorRampPalette(brewer.pal(8, "Set2"))(nb.clus)

dmizero.markers <- FindAllMarkers(dmizero, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

topn <- dmizero.markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_logFC)

# print results (plots)
pdf(paste0(resu,"FITSNE.pdf"))
DimPlot(dmizero, reduction="tsne", label=TRUE, pt.size=0.5, 
        cols=mycols, shape.by="orig.ident") + ggtitle("De Micheli D0")
dev.off()

pdf(paste0(resu,"HEATMAP.pdf"),width=13)
DoHeatmap(dmizero, features = topn$gene, group.colors=mycols) + NoLegend()
dev.off()

# save .rds object
saveRDS(dmizero,file=paste0(rdsdir,"dmizero_seu_fitsne.rds"))

#write table of all markers
write.table(dmizero.markers, paste0(resu,"ALLMARKERS_DeMicheliD0.txt"))
print("finished")
