



############## OTHERS SOURCES
source(file="~/INMG_SingleCell/scripts/functions_stock.R",local=T)





############## LIBRARIES
library(clusterProfiler)
library(AnnotationDbi)
library(cowplot)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyverse)
library(Seurat)
library(sctransform)
library(RColorBrewer)
library(stringr)
############## DATABASE
library(org.Mm.eg.db)
genome_annotation = org.Mm.eg.db

############## FUNCTION
# Differential expression + Gene Ontologies
# return DE & 2 graphs plots
GOtest_toy_1 = function(genome_annotation,genessymbols,genelistFC) {
  
  convert2entrezid=select(x = genome_annotation,
                          keys = genessymbols,
                          columns = c("ENTREZID","SYMBOL"),
                          keytype = "SYMBOL"
                          )
  
  DE = convert2entrezid$ENTREZID
  names(genelistFC) = DE
  
  rm(convert2entrezid)
  
  
  goALL <- enrichGO(DE, OrgDb = genome_annotation, ont="all", readable=TRUE)
  # goBP <- enrichGO(DE, OrgDb = genome_annotation, ont="BP", readable=TRUE)
  # hist(goBP@result[["p.adjust"]])  #simple view of FDR values in results 
  # selected <- subset(goBP, subset=goBP@result[["p.adjust"]]<=0.03)
  # dim(selected) 
  # topcat = dim(selected)[1] #rows number 
  # barplot(goBP,showCategory = topcat)
  
  p1=dotplot(goALL, split="ONTOLOGY") + facet_grid(ONTOLOGY~.,scale="free")
  p2=cnetplot(goALL, foldChange=genelistFC, circular=TRUE, colorEdge = TRUE)
  
  print(grid.arrange(p1,p2,nrow = 2,main = "Differential expression & Gene Ontology analysis"))
  
  return(DE)
}  
# Return graphs plots after clustering
# need cluster 
GOtest_toy_2 = function(clusteringDE,genome_annotation) {
  xx.formula <- compareCluster(Entrez~group,
                               data=clusteringDE,
                               fun="enrichGO", 
                               ont="all", 
                               OrgDb=genome_annotation)
  rm(mydf)
  
  print(dotplot(xx.formula, split="ONTOLOGY") + 
          facet_grid(ONTOLOGY~.,scale="free") +
          labs(x="cluster label (number of genes)"))
}
# genelistFC <- c(2.36354, 2.04829, 1.80765, 2.2143, 2.11138, 2.1023)
# genessymbols = c("Cxcl14", "Smoc2", "Hsd11b1", "Ugdh", "Fn1","Has1")
# DEtest = GOtest_toy_1(genome_annotation,genessymbols,genelistFC)
# clusteringDE <- data.frame(Entrez=DEtest, group= c(rep("0",3),rep("1",3)))
# GOtest_toy_2(clusteringDE,genome_annotation)


# Create a work repertory environment with rds & results subfolders
workEnvironment = function(pathListData) {
  #prloc="~/INMG_SingleCell/"
  #setwd(prloc)
  
  origintargetdf = data.frame()
  
  for (i in 1:length(pathListData)) {
    
    origintargetdf[i,1] = paste0("data/",pathListData[i])
    origintargetdf[i,2] = paste0(origintargetdf[i,1],"D0/")
    origintargetdf[i,3] = pathListData[i]

    #dir.create("./results"),recursive=T)
    #dir.create("./rds"),recursive=T)
  }
  colnames(origintargetdf) = c("origindir","targetdir","name_project")
  
  
  return(origintargetdf)
}
##pathListData = c("DeMicheli","DelOrso","Oprescu","Giordani")
##origintargetdf = workEnvironment(pathListData)

# Create mega SEURAT object which contains list of cso by data by project
# thresholdcell = 3
# thresholdfeaturemin = 200
# thresholdfeaturemax = 2500
# nbDIM = 13
# resolution = 0.5
CreateCSO = function(origintargetdf,thresholdcell,thresholdfeature) {

  listCSO = list()
  for (i in 1:dim(origintargetdf)[1]) {

    data_folder = list.dirs(origintargetdf[i,"origindir"])
    listCSOproject = list()
    for (j in 1:length(data_folder)) {
    
      if (str_detect(string = data_folder[j],pattern = ".txt")) {
        ############### txt files
        tmp_data = read.table(path = data_folder[j],
                              sep = "\t",
                              header = T, 
                              row.names = 1)
      } else if (str_detect(string = data_folder[j],pattern = ".csv")) {
        ############### csv files 
        tmp_data = read.csv(path = data_folder[j], 
                            sep=",", 
                            header = TRUE, 
                            row.names = 1)
      } else {
        ############### 10X files
        tmp_data = Read10X(data.dir = data_folder[j])
      }
      listCSOproject[[j]] = CreateSeuratObject(counts = tmp_data,
                                       project = origintargetdf[i,"name_project"],
                                       min.cells  = thresholdcell,
                                       min.features = thresholdfeature)
      rm(tmp_data)
    }
    rm(data_folder)
  }
  if (length(listCSOproject) > 1) {
    listCSO[[i]] = merge(x = listCSOproject[[1]],
                         y = listCSOproject[[2]],
                         add.cells.ids = c("wt1","wt2"),
                         project = origintargetdf[i,"name_project"]
    )
  } else {
    listCSO[[i]] = listCSOproject[[1]]
  }
  rm(listCSOproject)
  return(listCSO)
}
# Control QC of features & to do normalisation
QCseurat = function(data,thresholdfeaturemin,thresholdfeaturemax) {
  data[["percent.mt"]] = PercentageFeatureSet(object = data,
                                              pattern="^mt-")
  data = subset(x = data,
                subset = nFeature_RNA > thresholdfeaturemin & nFeature_RNA < thresholdfeaturemax & percent.mt <  5)
  data = NormalizeData(object = data,
                       verbose = F)
  data = FindVariableFeatures(object = data,
                              selection.method = "vst",
                              nfeatures = thresholdfeature,
                              verbose = F)

  pdf(paste0("./",targetdir,"preprocessSeu.pdf"))
  VlnPlot(data, features =c("nFeature_RNA","nCount_RNA","percent.mt"), 
          ncol = 3) 
  FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  DimHeatmap(data, dims=1:15, cells = 500, balanced = TRUE) 
  ElbowPlot(data)
  
  p1 = VariableFeaturePlot(data)
  LabelPoints(plot = p1,
              points = head(VariableFeaturePlot(data),20),
              repel = T)
  
  data <- ScaleData(data, 
                    features = rownames(data))
  
  data <- RunPCA(data, 
                 features = VariableFeatures(object = data))
  
  VizDimLoadings(data, dims = 1:2, reduction = "pca")
  
  data <- JackStraw(data, num.replicate = 100)
  JackStrawPlot(data, dims = 1:15)
  data <- ScoreJackStraw(data, dims = 1:20)
  ElbowPlot(data)
  
  dev.off()
  
  return(data)
}
# Clustering & non-linear reduction for visuals
clusterseurat = function(data,targetdir,project_name) {
  
  data = ScaleData(data, verbose = F)
  data = RunPCA(data, npcs = 30, verbose = F)
  data = FindNeighbors(data, dims = 1:10)
  data = FindClusters(data, resolution)
  data = RunUMAP(data, dims = 1:10)
  data = RunTSNE(data, dims = 1:10)

  nb.clus = max(as.integer(levels(data@meta.data$seurat_clusters)))+1
  mycols <- colorRampPalette(brewer.pal(8, "Set2"))(nb.clus)
  
  data.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  topn <- data.markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_logFC)
  
  # print results (plots)
  pdf(paste0("./",targetdir,"graph_reports.pdf"))
  p1 = DimPlot(data, 
          reduction="tsne", 
          label=TRUE, 
          pt.size=0.5, 
          cols=mycols, 
          shape.by="orig.ident") 
        + ggtitle(paste0(project_name))
  p2 = DimPlot(data, 
                reduction="umap", 
                label=TRUE, 
                pt.size=0.5, 
                cols=mycols, 
                shape.by="orig.ident") 
        + ggtitle(paste0(project_name))

  grid.arrange(p1,p2,ncol = 2)
  
  print(DoHeatmap(data, 
                  features = topn$gene, 
                  group.colors=mycols) + 
          NoLegend()
  )
  # save .rds object
  saveRDS(tmp,file=paste0("./",targetdir,"/rds/CSO.rds"))
  #write table of all markers
  write.table(tmp.markers, paste0("./",targetdir,"ALLMARKERS_",project_name,".txt"))
}



