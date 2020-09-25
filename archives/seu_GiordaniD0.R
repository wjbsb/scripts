#!/usr/bin/env Rscript
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
## VERY USEFUL: 
#https://satijalab.org/seurat/v3.0/merge_vignette.html

source(file="~/INMG_SingleCell/scripts/functions_stock.R",local=T)

prloc="~/INMG_SingleCell/"
setwd(prloc)

resu = "results/GiordaniD0/"
dir.create(resu,recursive = T)
rdsdir = "rds/GiordaniD0/"
dir.create(rdsdir,recursive = T)


wt1 <- read.csv(paste0("data/GiordaniD0/GSM3520458_20171018_uninjured_wt_filtered.csv"), 
                sep=",", header=TRUE, row.names=1)
wt2 <- read.csv(paste0("data/GiordaniD0/GSM3520459_20180917_uninjured_wt_filtered.csv"),
                sep=",", header=TRUE, row.names=1)


gio1 <- CreateSeuratObject(wt1, project="Giordani", min.cells=3, min.features=200 )
gio1
gio2 <- CreateSeuratObject(wt2, project="Giordani", min.cells=3, min.features=200)
gio2

# merge the two replicates:
gio <- merge(gio1, y=gio2, add.cell.ids=c("wt1","wt2"), project = "Giordani")
gio

dim(gio[["RNA"]]@counts) 

# first steps seurat until PCA
gio <- NormFeatScalePCA(gio,2500,5) # function from functions_stock.R

pdf(paste0(resu,"preprocessSeu.pdf"))
VlnPlot(gio, features =c("nFeature_RNA","nCount_RNA","percent.mt"), 
        ncol = 3) 
FeatureScatter(gio, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
DimHeatmap(gio, dims=1:15, cells = 500, balanced = TRUE) 
ElbowPlot(gio)
dev.off()

# clustering and non-linear reduction for visuals
gio <- KNNplusFITSNE(gio, 15, 0.5)
# howmanyclusters we have found out:
nb.clus = max(as.integer(levels(gio@meta.data$seurat_clusters)))+1
mycols <- colorRampPalette(brewer.pal(8, "Set2"))(nb.clus)

gio.markers <- FindAllMarkers(gio, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


topn <- gio.markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_logFC)

# print results (plots)
pdf(paste0(resu,"FITSNE.pdf"))
DimPlot(gio, reduction="tsne", label=TRUE, pt.size=0.5, 
        cols=mycols, shape.by="orig.ident") + ggtitle("Giordani D0")
dev.off()

pdf(paste0(resu,"HEATMAP.pdf"),width=13)
DoHeatmap(gio, features = topn$gene, group.colors=mycols) + NoLegend()
dev.off()

# save .rds object
saveRDS(gio,file=paste0(rdsdir,"gio_seu_fitsne.rds"))


#write table of all markers
write.table(gio.markers, paste0(resu,"ALLMARKERS_GiordaniD0.txt"))
print("finished")


