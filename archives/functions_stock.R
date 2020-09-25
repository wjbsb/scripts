# Frequently used functions
# 
# --
# Joha 2020


NormFeatScalePCA <- function(seu, nFeatRNA, percentmit){
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern= "^mt-")
  seu <- subset(seu, subset = nFeature_RNA > 200 & 
                  nFeature_RNA < nFeatRNA & percent.mt < percentmit)
  seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor=10000)
  seu <- FindVariableFeatures(seu, selection.method="vst", nfeatures= 2000)
  all.genes <- rownames(seu)
  seu <- ScaleData(seu, features = all.genes)
  seu <- RunPCA(seu,features=VariableFeatures(object = seu))  
  return(seu)
}


KNNplusFITSNE <- function(seu, nbDIM, resolu, perplex){
  #note: more resolution more groups,range 0.4-1.12
  seu <- FindNeighbors(seu, dims=1:nbDIM)
  seu <- FindClusters(seu, resolution=resolu)
  # Run FIt-SNE
  source("~/INMG_SingleCell/programsSC/FIt-SNE/fast_tsne.R", chdir=T)
  nbCELLS.div100 <- round(dim(seu)[2]/100) #nombre de cellules divisé par 100
  seu <- RunTSNE(object = seu,
                         perplexity=nbCELLS.div100, 
                         reduction="pca",
                         dims=1:nbDIM,
                         tsne.method = "FIt-SNE",
                         nthreads=4,
                         reduction.key="FItSNE_",
                         fast_tsne_path="~/INMG_SingleCell/programsSC/FIt-SNE/bin/fast_tsne",
                         max_iter=1000)
   return(seu)
}

