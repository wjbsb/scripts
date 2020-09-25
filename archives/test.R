suppressMessages(require(Seurat))
suppressMessages(require(scater))
suppressMessages(require(Matrix))

prloc="/Users/williamjarassier/TRAVAIL/snRNA_Johana/INMG_SingleCell/"
dorso1.data <- Read10X(data.dir=paste0(prloc,"data/DellOrsoD0/dorsowt1"))
dorso2.data <- Read10X(data.dir=paste0(prloc,"data/DellOrsoD0/dorsowt2"))
gio1.data <- read.csv(paste0(prloc,"data/GiordaniD0/GSM3520458_20171018_uninjured_wt_filtered.csv"), sep=",", header=TRUE, row.names=1)
gio2.data <- read.csv(paste0(prloc, "data/GiordaniD0/GSM3520459_20180917_uninjured_wt_filtered.csv"),sep=",", header=TRUE, row.names=1)
dmi.data <- read.table(paste0(prloc,"data/DeMicheliD0/rawdataD0.txt"), sep="\t",header=T, row.names=1)

dataraw=list(dorso1.data,dorso2.data,gio1.data,gio2.data,dmi.data)
rm(dorso1.data,dorso2.data,gio1.data,gio2.data,dmi.data)
project_name=c("dorso1","dorso2","giordani1","giodarni2","demicheli")


dataseurat=list()
for (i in 1:length(dataraw)) {
  dataseurat[i]=CreateSeuratObject(dataraw[[i]],project=project_name[i], min.cells=3, min.features=200)
}
names(dataseurat)=project_name
rm(dataraw)

alldata = merge(x=dataseurat$dorso1,y=c(dataseurat$dorso2,
                                        dataseurat$giordani1,dataseurat$giodarni2,
                                        dataseurat$demicheli),
                add.cell.ids = c("dorso1_wt1","dorso_wt2","gio_wt1","gio_wt2","demichili"))
alldata$orig.ident[alldata$orig.ident %in% "D0"]="demichelli"
alldata$orig.ident[alldata$orig.ident %in% "i1"]="giordani1"
alldata$orig.ident[alldata$orig.ident %in% "i2"]="giordani2"
#table(alldata@meta.data$orig.ident)
#QC
alldata[["percent.mt"]] <- PercentageFeatureSet(alldata, pattern = "^mt-")


library(grid)
library(gridExtra)

QCplot = function(data) {
  return(grid.arrange(
    grid.arrange(
      vlnqc(data@meta.data,"nCount_RNA"),
      vlnqc(data@meta.data,"nFeature_RNA"),
      vlnqc(data@meta.data,"percent.mt"),
      nrow=3),
    grid.arrange(
      corrqc(data@meta.data,"nCount_RNA","nFeature_RNA"),
      corrqc(data@meta.data,"nCount_RNA","percent.mt"),
      nrow=1),
    nrow=2))
}
vlnqc = function(data,var) {
  return(
    ggplot(data,aes(x=data[,1],y=data[,var])) + 
      geom_violin(fill="darkblue") + 
      geom_jitter(shape=16, position=position_jitter(0.2),alpha=1/5) + 
      theme_classic() +
      theme(plot.title = element_text( size=8, face="bold.italic")) +
      labs(x="",y=var))
}
corrqc = function(data,var1,var2) {
  corr=round(cor(x=data[,var1],y=data[,var2],method = "pearson"),2)
  return(
    ggplot(data,aes(x=data[,var1],y=data[,var2],color=data[,1])) + 
      geom_jitter(shape=16, position=position_jitter(0.2),alpha=1/5) + 
      theme_classic() +
      ggtitle(paste0("Pearson Corr=",corr)) +
      theme(plot.title = element_text( size=8, face="bold.italic")) +
      labs(x=var1,y=var2))
}
qcdata = alldata@meta.data
#grid.arrange(
#  corrqc(alldata@meta.data,"nCount_RNA","nFeature_RNA"),
#  corrqc(alldata@meta.data,"nCount_RNA","percent.mt"),
#  nrow=1)
rm(qcdata)
alldata_filtered = subset(alldata, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
QCplot(alldata)
QCplot(alldata_filtered)
save(alldata,file="/Users/williamjarassier/TRAVAIL/DATA/QC_single_cell/alldata.Rdata")
alldata_filtered <- NormalizeData(alldata_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
alldata_filtered <- FindVariableFeatures(alldata_filtered, selection.method = "vst", nfeatures = 2000)
#plot1 <- VariableFeaturePlot(alldata_filtered)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#plot1 + plot2
all.genes <- rownames(alldata_filtered)
alldata_filtered <- ScaleData(alldata_filtered, features = all.genes)
alldata_filtered <- RunPCA(alldata_filtered, features = VariableFeatures(object = alldata_filtered))
#VizDimLoadings(alldata_filtered, dims = 1:2, reduction = "pca")
#DimPlot(alldata_filtered, reduction = "pca")
alldata_filtered <- JackStraw(alldata_filtered, num.replicate = 100)
alldata_filtered <- ScoreJackStraw(alldata_filtered, dims = 1:20)
alldata_filtered <- FindNeighbors(alldata_filtered, dims = 1:10)
alldata_filtered <- FindClusters(alldata_filtered, resolution = 0.5)
alldata_filtered <- RunUMAP(alldata_filtered, dims = 1:10)
alldata_filtered <- RunTSNE(alldata_filtered, dims = 1:13, resolution = 0.5)
#DimPlot(alldata_filtered, reduction = "umap")
alldata_filtered.markers <- FindAllMarkers(alldata_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
alldata_filtered.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#VlnPlot(alldata_filtered, features = c("Pax7"), slot = "counts", log = TRUE)
#FeaturePlot(alldata_filtered, features = c("Pax7"))

top10 <- alldata_filtered.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#DoHeatmap(alldata_filtered, features = top10$gene) + NoLegend()

#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
sweep.res.alldata_filtered <- paramSweep_v3(alldata_filtered, PCs = 1:10, sct = FALSE)
sweep.stats_muscle <- summarizeSweep(sweep.res.alldata_filtered, GT = FALSE)
bcmvn_muscle <- find.pK(sweep.stats_muscle)

#table(Idents(alldata_filtered), alldata_filtered$clusters)
#homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*length(alldata_filtered@meta.data$orig.ident))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
#nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
alldata_filtered <- doubletFinder_v3(alldata_filtered, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
#alldata_filtered <- doubletFinder_v3(alldata_filtered, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)

save(alldata_filtered,file="/Users/williamjarassier/TRAVAIL/DATA/QC_single_cell/alldata_filtered.Rdata")

umapdata=data.frame(alldata_filtered@reductions$umap@cell.embeddings)
tsnedata=data.frame(alldata_filtered@reductions$tsne@cell.embeddings)
metadata=alldata_filtered@meta.data
umapdata$id=as.character(rownames(umapdata))
metadata$id=as.character(rownames(metadata))
tsnedata$id=as.character(rownames(tsnedata))
plotdoublet = inner_join(x=metadata,y=umapdata,by="id")
plotdoublet = inner_join(x=plotdoublet,y=tsnedata,by="id")

rm(umapdata,metadata)
colnames(plotdoublet) = c("groupe","nCount_RNA","nFeature_RNA","percent.mt",
                          "RNA_snn_res","cluster","pANN","DFclass","id",
                          "UMAP_1","UMAP_2","TSNE_1","TSNE_2")
plotdoublet$DFclass[plotdoublet$DFclass %in% "Doublet" & plotdoublet$pANN <= 0.5]="Doublets - Low confidence"
plotdoublet$DFclass[plotdoublet$DFclass %in% "Doublet" & plotdoublet$pANN > 0.5 ]="Doublets - High confidence"
plotdoublet$groupe[plotdoublet$groupe %in% "D0"]="demichelli"
plotdoublet$groupe[plotdoublet$groupe %in% "i1"]="giordani1"
plotdoublet$groupe[plotdoublet$groupe %in% "i2"]="giordani2"

save(plotdoublet,file="/Users/williamjarassier/TRAVAIL/DATA/QC_single_cell/plotdoublet.Rdata")

#grid.arrange(
#  ggplot(plotdoublet,aes(x=UMAP_1,y=UMAP_2,color=DFclass)) +
#               geom_point() + 
#               scale_color_manual(values=c("darkred", "darkblue", "gray")) +
#               theme_classic() +
#               theme(plot.title = element_text( size=8, face="bold.italic")),
#  ggplot(plotdoublet,aes(x=UMAP_1,y=UMAP_2,color=clusters)) +
#    geom_jitter(shape=16, position=position_jitter(0.2),alpha=1/10) + 
#    theme_classic() +
#    theme(plot.title = element_text( size=8, face="bold.italic")),
#  ggplot(plotdoublet,aes(x=UMAP_1,y=UMAP_2,color=groupe)) +
#    geom_jitter(shape=16, position=position_jitter(0.2),alpha=1/10) + 
#    theme_classic() +
#    theme(plot.title = element_text( size=8, face="bold.italic")),
#  FeaturePlot(alldata_filtered, features = c("Pax7")),
#  nrow=2)

library(scran)
score=doubletCells(alldata_filtered@assays$RNA@counts)
#plotdoublet2 = data.frame(
#  score=score,
#  cluster=plotdoublet$cluster,
#  groupe=plotdoublet$groupe,
#  DFclass=plotdoublet$DFclass,
#  id=plotdoublet$id
#)

#ggplot(plotdoublet2) + 
#  geom_boxplot(aes(y=log10(score),x=cluster,color=groupe)) +
#  facet_grid(~DFclass)
#table(plotdoublet2$DFclass, plotdoublet2$groupe)

ggplot(plotdoublet,aes(x=TSNE_1,y=TSNE_2,color=DFclass)) +
            geom_point() + 
            scale_color_manual(values=c("darkred", "darkblue", "gray")) +
            theme_classic() +
            theme(plot.title = element_text( size=8, face="bold.italic"))
TSNEPlot(alldata_filtered)
