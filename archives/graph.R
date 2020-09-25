##################### LIBRARIES
library(gridExtra)
library(ggplot2)
library(dplyr)
library(Seurat)
library(RColorBrewer)
library(Rsamtools)
library(stringr)
library(UpSetR)
library(ggforce)
##################### PATH
path="E:/TRAVAIL/William/Lyon/DATA/"
prloc=paste0(path,"INMG_SingleCell-master/")
source(file=paste0(prloc,"scripts/functions_stock.R"),local=T)
paths_output=c(paste0(prloc,"results/DellOrsoD0/"),
               paste0(prloc,"results/DeMicheliD0/"),
               paste0(prloc,"results/GiordaniD0/"))
##################### VARIABLE
project_name=c("DellOrso","DeMicheli","Giordani")
##################### FUNCTION
loading_data = function(path,Rdata,pj) {
  # loop on data folder
  data_list=list()
  for (i in 1:length(pj)) {
    # find Rdata in data folder
    if (file.exists(paste0(path,"data/",pj[i],"D0/",Rdata,".Rdata"))) {
      print(paste0("Loading ",path,"data/",pj[i],"D0/",Rdata,".Rdata"))
      load(paste0(path,"data/",pj[i],"D0/",Rdata,".Rdata"))
      data_list[[i]]=data
    } 
    else {
      data_list[[i]]=NA
    }
  }
  names(data_list)=pj
  return(data_list)
}
VFmodeUPSET=function(vst,mvp) {
  df=list(vst,mvp)
  names(df)=c("vst","mvp")
  
  return(
    upset(df %>% stack() %>% unique() %>% cbind(val = 1) %>% 
            reshape2::acast(values ~ ind, value.var = "val", fill = 0) %>% as.data.frame(),
          nsets = 2, matrix.color = "#DC267F",main.bar.color = "#648FFF", sets.bar.color = "#FE6100")
  )
}
FVF_plot=function(data,pj){
  
  tmp1=table(factor(data@assays$RNA@meta.features$vst.variable, levels = c(TRUE,FALSE)))
  tmp2=table(factor(data@assays$RNA@meta.features$mvp.variable, levels = c(TRUE,FALSE)))
  lab1=c(paste0("Variable Count : ",tmp1[1]),paste0("Non-variable Count : ",tmp1[2]))
  lab2=c(paste0("Variable Count : ",tmp2[1]),paste0("Non-variable Count : ",tmp2[2]))
  return(
    grid.arrange(
      ggplot(data@assays$RNA@meta.features,aes(x=log10(vst.mean),y=vst.variance.standardized,colour=vst.variable)) + 
        geom_count(alpha=0.5) +
        theme_classic() +
        ggtitle(pj,"") +
        theme(plot.title = element_text( size=8, face="bold.italic")) +
        labs(y="Standardized Variance",x="Log10(Average Expression)",
             size="Quantity",colour="") +
        scale_colour_hue(labels = lab1),
      ggplot(data@assays$RNA@meta.features,aes(x=mvp.mean,y=mvp.dispersion,colour=mvp.variable)) + 
        geom_point(alpha=0.5) +
        theme_classic() +
        ggtitle(pj,"") +
        theme(plot.title = element_text( size=8, face="bold.italic")) +
        labs(y="Dispersion",x="Average Expression",
             size="Quantity",colour="") +
        scale_colour_hue(labels = lab2),
      nrow=3
    )
  )
}

##################### LOADING
NormFeatScalePCA_list=loading_data(prloc,"NormFeatScalePCA",project_name)
KNNplusFITSNE_list=loading_data(prloc,"KNNplusFITSNE",project_name)
FindAllMarkers_list=loading_data(prloc,"FindAllMarkers",project_name)
FindVariableFeatures_list=loading_data(prloc,"FindVariableFeatures",project_name)
##################### TREATMENT
FindVariableFeatures_list=list()
for (i in 1:length(KNNplusFITSNE_list)) {
  FindVariableFeatures_list[i] = FindVariableFeatures(KNNplusFITSNE_list[[i]],selection.method = "vst",nfeatures=2000)
  vst=FindVariableFeatures_list[[i]]@assays$RNA@var.features
  FindVariableFeatures_list[i] = FindVariableFeatures(KNNplusFITSNE_list[[i]],selection.method = "mvp",nfeatures=2000)
  mvp=FindVariableFeatures_list[[i]]@assays$RNA@var.features

  save(FindVariableFeatures_list[i],file=paste0(prloc,"data/",project_name[i],"D0/FindVariableFeatures.Rdata"))
  #print(VFmodeUPSET(vst,mvp))
}


########BOUCLE A FAIRE
ScalingData_list=list()
PCAData_list=list()
FindAllMarkers_list=list()
top10_list=list()
for (i in 1:length(KNNplusFITSNE_list)){
  all.genes=rownames(FindVariableFeatures_list[[i]])
  ScalingData_list[i] <- ScaleData(KNNplusFITSNE_list[[i]], features = all.genes)
  
  PCAData_list[i] <- RunPCA(ScalingData_list[[i]], features = VariableFeatures(object = ScalingData_list[[i]]))
  
  #VizDimLoadings(PCAData_list[[i]], dims = 1:2, reduction = "pca")
  #DimPlot(PCAData_list[[i]], reduction = "pca")
  #DimHeatmap(PCAData_list[[i]], dims = 1, cells = 500, balanced = TRUE)
  #DimHeatmap(PCAData_list[[i]], dims = 1:15, cells = 500, balanced = TRUE)
  PCAData_list[i] <- JackStraw(PCAData_list[[i]], num.replicate = 100)
  PCAData_list[i] <- ScoreJackStraw(PCAData_list[[i]], dims = 1:20)
  #JackStrawPlot(PCAData_list[[i]], dims = 1:15)
  #ElbowPlot(PCAData_list[[i]])
  PCAData_list[i] <- FindNeighbors(PCAData_list[[i]], dims = 1:10)
  PCAData_list[i] <- FindClusters(PCAData_list[[i]], resolution = 0.5)
  PCAData_list[i] <- RunUMAP(PCAData_list[[i]], dims = 1:10)
  DimPlot(PCAData_list[[i]], reduction = "umap")
  data <- FindAllMarkers(PCAData_list[[i]], only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  save(data,file=paste0(prloc,"data/",project_name[i],"D0/FindAllMarkers.Rdata"))
  
}


PCAData_list[i] <- RunTSNE(PCAData_list[[i]], dims = 1:10)

cell_cluster_all=data.frame(cluster=integrated@active.ident)
cell_cluster_all$id=rownames(cell_cluster)
values=data.frame(integrated@reductions$tsne@cell.embeddings[,1:2])
values$id=rownames(values)

a=data.frame(cluster=PCAData_list[[1]]@active.ident)
a$id=rownames(a)
a$grp=project_name[1]
b=data.frame(cluster=PCAData_list[[2]]@active.ident)
b$id=rownames(b)
b$grp=project_name[2]
c=data.frame(cluster=PCAData_list[[3]]@active.ident)
c$id=rownames(c)
c$grp=project_name[3]

d=full_join(a,b)
d=full_join(d,c)
data_plot=inner_join(x=values,y=cell_cluster_all)
data_plot=inner_join(x=data_plot,y=d)
pdf(file=paste0("E:/TRAVAIL/William/Lyon/DATA/INMG_SingleCell-master/results/DimPlot.pdf"))
grid.arrange(
  ggplot(data_plot,aes(x=tSNE_1,y=tSNE_2,colour=grp)) +
    geom_jitter(shape=16, position=position_jitter(0.2),alpha=1/5) + 
    theme_classic() +
    theme(plot.title = element_text( size=8, face="bold.italic")),
  FeaturePlot(integrated,features="Pax7",reduction="tsne"),
  ggplot(data_plot,aes(x=tSNE_1,y=tSNE_2,colour=cluster)) +
    geom_point() + 
    theme_classic() +
    theme(plot.title = element_text( size=8, face="bold.italic")),
  nrow=2
)
grid.arrange(
  ggplot(data_plot,aes(x=tSNE_1,y=tSNE_2,colour=grp)) +
    geom_jitter(shape=16, position=position_jitter(0.2),alpha=1/5) + 
    theme_classic() +
    theme(plot.title = element_text( size=8, face="bold.italic")),
  FeaturePlot(integrated,features="Myod1",reduction="tsne"),
  ggplot(data_plot,aes(x=tSNE_1,y=tSNE_2,colour=cluster)) +
    geom_point() + 
    theme_classic() +
    theme(plot.title = element_text( size=8, face="bold.italic")),
  nrow=2
)
dev.off()


grid.arrange(
  ggplot(data_plot,aes(x=tSNE_1,y=tSNE_2,colour=grp)) +
    geom_jitter(shape=16, position=position_jitter(0.2),alpha=1/5) + 
    theme_classic() +
    theme(plot.title = element_text( size=8, face="bold.italic")),
  FeaturePlot(integrated,features="Myog",reduction="tsne"),
  ggplot(data_plot,aes(x=tSNE_1,y=tSNE_2,colour=cluster)) +
    geom_point() + 
    theme_classic() +
    theme(plot.title = element_text( size=8, face="bold.italic")),
  nrow=2
)

VlnPlot(integrated, features = c("Pax7", "Myog"), slot = "counts", log = TRUE)

pdf(file=paste0("E:/TRAVAIL/William/Lyon/DATA/INMG_SingleCell-master/results/DimPlot.pdf"))
  DimPlot(PCAData_list[[1]], reduction = "umap")
  DimPlot(PCAData_list[[2]], reduction = "umap")
  DimPlot(PCAData_list[[3]], reduction = "umap",group.by = "orig.ident")
  DimPlot(PCAData_list[[1]], reduction = "pca")
  DimPlot(PCAData_list[[2]], reduction = "pca")
  DimPlot(PCAData_list[[3]], reduction = "pca",group.by = "orig.ident")
  DimPlot(PCAData_list[[1]], reduction = "tsne")
  DimPlot(PCAData_list[[2]], reduction = "tsne")
  DimPlot(PCAData_list[[3]], reduction = "tsne",group.by = "orig.ident")
dev.off()
all_data <- SplitObject(KNNplusFITSNE_list, split.by = "project.name")


for (i in 1:length(KNNplusFITSNE_list)) {
  all_data[i] = NormalizeData(KNNplusFITSNE_list[[i]],verbose=F)
  all_data[i] = FindVariableFeatures(all_data[[i]],selection.method = "vst",nfeatures=2000)
}

reference = all_data[project_name]
anchors <- FindIntegrationAnchors(object.list = reference, dims = 1:30)
integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

library(cowplot)
library(patchwork)

integrated = ScaleData(integrated,verbose=F)
integrated = RunPCA(integrated,npcs=30,verbose=F)
integrated = RunUMAP(integrated, reduction="pca",dims=1:30)
integrated = RunTSNE(integrated, reduction="pca",dims=1:30)

p1 = DimPlot(integrated,reduction="umap",group.by = "tech")
p2 = DimPlot(integrated,reduction="umap",group.by = "celltype",label=T,repel=T) + NoLegend()

top10 <- FindAllMarkers_list[[1]] %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

biomarkers=read.csv(file="C:/Users/illir/Downloads/CellMarker.csv")
biomarkers_muscle=filter(biomarkers,Tissue=="Muscle")


pdf(file=paste0("E:/TRAVAIL/William/Lyon/DATA/INMG_SingleCell-master/results/Expression_Level.pdf"))
genetop=top10$gene
for (i in 1:length(genetop)) {
  print(paste0(genetop[i]," ",i,"/",length(genetop)))
  print(grid.arrange(
    VlnPlot(PCAData_list[[1]], features = genetop[i]),
    VlnPlot(PCAData_list[[2]], features = genetop[i]),
    VlnPlot(PCAData_list[[3]], features = genetop[i]),
    nrow=3))

}
dev.off()




#################################################################### PLOTTING


a=ExpressionLevelPlot(NormFeatScalePCA_list,"tomato1","Expression Level of nCount_RNA","nCount_RNA",project_name)
b=ExpressionLevelPlot(NormFeatScalePCA_list,"tan2","Expression Level of nFeature_RNA","nFeature_RNA",project_name)
c=ExpressionLevelPlot(NormFeatScalePCA_list,"steelblue2","Expression Level of %MT","percent.mt",project_name)
grid.arrange(a,b,c,nrow=3)
rm(a,b,c)
QCplot(NormFeatScalePCA_list,"nCount_RNA VS nFeature_RNA","nCount_RNA","nFeature_RNA",project_name)
QCplot(NormFeatScalePCA_list,"nCount_RNA VS %MT","nCount_RNA","percent.mt",project_name)






