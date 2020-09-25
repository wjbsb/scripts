####################### IMPORT
currentDirectory = "/Volumes/LaCie/InstitutNeuroMyoGene_DATA/SatCellAlert"
setwd(currentDirectory)

project = c("GSE47177","GSE55490")
setwd("./DATA")

createSDRF= function(TT,AD,CO,SN,CI,CP) {
  
  SDRF = data.frame(
    Tech.type = TT,
    Array.Data.File = AD,
    Characteristics.organism = CO,
    Source.Name = SN,
    Characteristics.individual = CI,
    Characteristics.phenotype = CP
  )
  rm(TT,AD,CO,SN,CI,CP)
  
  rownames(SDRF) = paste0(SDRF$Source.Name,"_",SDRF$Characteristics.phenotype)
  return(SDRF)
}

TT=c()
CO=c()
CI=c()
CP=c()
SN=c()
AD = dir(paste0("./",project[2]),pattern=".CEL",full.names = F)
for(i in 1:9) {
  CO[i]="Mus musculus"
  TT[i]="array assay"
  CI[i]=strsplit(x=AD[i],split="_")[[1]][3]
  CP[i]=ifelse((i >= 1 & i <= 5),"Control","Alert")
  SN[i] = gsub(pattern = ".CEL",replacement = "", x = AD[i])
}
ALLSDRF = createSDRF(TT,AD[1:9],CO,SN,CI,CP)

library(Biobase)
library(stringr)
library(ggplot2)
library(gridExtra)
plotPCA=function(data,title) {
  exp = exprs(data)
  PCA_raw = prcomp(t(exp),scale=F)
  percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
  sd_ratio <- sqrt(percentVar[2] / percentVar[1])
  dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                       Phenotype = pData(data)$Characteristics.phenotype,
                       Individual = pData(data)$Characteristics.individual)
  
  return(ggplot(dataGG, aes(PC1, PC2)) +
           geom_point(aes(shape = Individual, colour = Phenotype)) +
           ggtitle(title) +
           xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
           ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
           coord_fixed(ratio = sd_ratio) +
           scale_shape_manual(values = c(4,15)) + 
           scale_color_manual(values = c("darkorange2", "dodgerblue4")))
  
  rm(exp,PCA_raw,percentVar, sd_ratio, dataGG)
}
plotRLE=function(data,title) {
  row_medians_assayData = rowMedians(as.matrix(data))
  RLE_data = sweep(exprs(data),1,row_medians_assayData)
  RLE_data = as.data.frame(RLE_data)
  RLE_data_gathered = tidyr::gather(RLE_data, sample_array, log2_expression_deviation)
  
  #tmp = c("CTL-40PA","CTL-12PA","Wnt-12PA","Wnt-40PA")
  tmp = c("Control","Alert")
  for (i in 1:length(tmp)) {
    RLE_data_gathered$Sample[str_detect(string=RLE_data_gathered$sample_array,pattern = tmp[i])] = tmp[i]
  }
  
  return(ggplot(RLE_data_gathered, aes(x=sample_array,y=log2_expression_deviation,fill=Sample)) +
           geom_boxplot(outlier.shape = NA) + 
           ylim(c(-2, 2)) + 
           theme_bw() +
           ggtitle(title) +
           theme(axis.text.x = element_text(colour = "aquamarine4", 
                                            angle = 60, size = 6.5, hjust = 1 ,
                                            face = "bold"))
         
  )
}
plotIntensity = function(data1,data2) {
  
  logtmp = log(as.data.frame(data1@assayData$exprs))
  
  tmp = data.frame()
  for (i in 1:dim(logtmp)[2]) {
    tmp2 = data.frame(col = logtmp[,i])
    tmp2$grp = colnames(logtmp)[i]
    tmp = rbind(tmp,tmp2)
    rm(tmp2)
  }
  
  p1=ggplot(data.frame(tmp)) + 
    geom_density(aes(x=col,colour=grp)) +
    theme_bw() +
    labs(title="Density estimate for log-intensities before RMA",x="log-intensity", y = "Density")
  
  logtmp = log(as.data.frame(data2@assayData$exprs))
  
  tmp = data.frame()
  for (i in 1:dim(intensities)[2]) {
    tmp2 = data.frame(col = logtmp[,i])
    tmp2$grp = colnames(logtmp)[i]
    tmp = rbind(tmp,tmp2)
    rm(tmp2)
  }
  
  p2 = ggplot(data.frame(tmp)) + 
    geom_density(aes(x=col,colour=grp)) +
    theme_bw() +
    labs(title="Density estimate for intensities after RMA",x="intensity", y = "Density")
  
  return(grid.arrange(p1,p2,ncol = 2))
  rm(tmp)
}


#for (i in 1:length(SDRFlist)) {
#  SDRFlist[[i]] = AnnotatedDataFrame(SDRFlist[[i]])
#}
SDRFlist = list()
SDRFlist[[1]] = AnnotatedDataFrame(ALLSDRF)

library(oligo)
filescel = dir("./GSE55490",pattern=".CEL",full.names = T)
raw_data <- read.celfiles(filenames = filescel[1:9],
                          verbose = FALSE, 
                          phenoData = SDRFlist[[1]])
eset_raw=oligo::rma(raw_data, normalize = F)
p1=plotPCA(eset_raw,"PCA plot of the  raw expression data")
p2=plotRLE(eset_raw,"Relative Log Expression before calibration")
#grid.arrange(p1,p2,ncol=1)
####################### RELATIVE LOG EXPRESSION DATA QUALITY ANALYSIS
eset_norm=oligo::rma(raw_data, normalize = T)
p3=plotPCA(eset_norm,"PCA plot of the calibrated expression data")
p4=plotRLE(eset_norm,"Relative Log Expression after calibration")
p5 = plotIntensity(eset_raw,eset_norm)
#grid.arrange(p3,p4,ncol=1)

pdf("./resultsQC.pdf")
print(grid.arrange(p1,p3,p2,p4,ncol=2,nrow=2))
print(p5)
dev.off()
rm(p1,p2,p3,p4)

annotation_for_heatmap = SDRFlist[[i]]@data[,c("Characteristics.phenotype","Characteristics.individual")]
colnames(annotation_for_heatmap) = c("phenotype_names","condition")

dists = as.matrix(dist(t(exprs(eset_norm)),method = "manhattan"))
rownames(dists) = row.names(pData(eset_norm))
colnames(dists) = rownames(annotation_for_heatmap)
diag(dists) <- NA

library(pheatmap)
pheatmap(dists,
         annotation_col= annotation_for_heatmap,
         legend = TRUE, 
         main = "Clustering heatmap for the calibrated samples")

rm(dists)

####################### LINEAR MODEL
library(limma)

plottingVolcanoPlot = function(data_final,design,title,coef) {
  fit1 <- lmFit(data_final, design)
  fit1 <- eBayes(fit1)
  #lod <- -log10(ebayes[["p.value"]][,2]) 
  #mtstat<- ebayes[["t"]][,2]
  
  tmp =  topTable(fit1, coef, number = dim(fit1)[1])
  print(EnhancedVolcano(tmp,x="logFC",y="P.Value",lab=tmp$SYMBOL,title = title)) 
  return(tmp)
}
results = list()

library(EnhancedVolcano)
#FILTERING
data_medians = rowMedians(exprs(eset_norm))
hist(data_medians,100,col="cornsilk1",freq=F,
     main="Histogram of the median intensities",
     border="antiquewhite4",
     xlab="Median intensities")
abline(v = 5.5, col = "coral4", lwd = 2)
abline(v = 8, col = "coral4", lwd = 2)
no_of_samples = table(paste0(pData(eset_norm)$Characteristics.phenotype, "_", 
                             pData(eset_norm)$Characteristics.individual))
samples_cutoff = min(no_of_samples)
idx_man_threshold = apply(exprs(eset_norm),1,function(x){sum(x > 5.5 & x < 8) >= samples_cutoff})
table(idx_man_threshold)

data_manfiltered <- subset(eset_norm, idx_man_threshold)

rm(data_medians,no_of_samples,samples_cutoff,idx_man_threshold)
####################### 
#ANNOTATION
#BiocManager::install("mogene10st.db")
library(pd.mogene.1.0.st.v1)
#BiocManager::install("affycoretools")
library(affycoretools)
data_annotation = annotateEset(eset_norm, pd.mogene.1.0.st.v1)@featureData@data
data_annotation <- subset(data_annotation, !is.na(SYMBOL))

####################### EXCLUDE PROBE ID
probe_stats = data_annotation
nrow(probe_stats)
ids_to_exlude <- (featureNames(data_manfiltered) %in% probe_stats$PROBEID)
table(ids_to_exlude)
data_final <- subset(data_manfiltered, ids_to_exlude)

rm(probe_stats,ids_to_exlude)
# ####################### EXCLUDE FEATURE 
fData(data_final)$PROBEID <- rownames(fData(data_final))
fData(data_final) <- left_join(fData(data_final), data_annotation)
rownames(fData(data_final)) <- fData(data_final)$PROBEID  

Condition = factor(pData(data_final)$Characteristics.individual)
PH = factor(pData(data_final)$Characteristics.phenotype)

design = model.matrix(~Condition+PH)
results = plottingVolcanoPlot(data_final,design,"Samples QSC Control vs Samples CSC Alert",2)

datacaroQSC = read.csv(file = "/Volumes/LaCie/InstitutNeuroMyoGene_DATA/GSE144871/publication_caro/QSC_Gli3KO_vs_QSC_control.csv")

totalfilteredAQ = inner_join(datacaroQSC,results,by=c("mgi_symbol"="SYMBOL"))
totalfiltered = filter(totalfilteredAQ, padj < 0.05 & adj.P.Val < 0.05)
thresholdFC = 1
commun_filtered_1= filter(totalfilteredAQ,(log2FoldChange > thresholdFC & logFC > thresholdFC))
commun_filtered_2= filter(totalfilteredAQ,(log2FoldChange < thresholdFC & logFC < thresholdFC))
different_1 = filter(totalfilteredAQ,(log2FoldChange < thresholdFC & logFC > thresholdFC))
different_2 = filter(totalfilteredAQ,(log2FoldChange > thresholdFC & logFC < thresholdFC))

ggplot(totalfilteredAQ,aes(log2FoldChange,logFC)) + 
  geom_point() + 
  ggrepel::geom_text_repel(data=commun_filtered_1,
                           aes(x=log2FoldChange,
                               y=logFC,
                               label=mgi_symbol),
                           color="darkgreen") +
  ggrepel::geom_text_repel(data=commun_filtered_2,
                           aes(x=log2FoldChange,
                               y=logFC,
                               label=mgi_symbol),
                           color="darkred") +
  ggrepel::geom_text_repel(data=different_1,
                           aes(x=log2FoldChange,
                               y=logFC,
                               label=mgi_symbol),
                           color="darkblue") +
  ggrepel::geom_text_repel(data=different_2,
                           aes(x=log2FoldChange,
                               y=logFC,
                               label=mgi_symbol),
                           color="black") +
  geom_abline(intercept=0,slope = 1) + 
  geom_vline(xintercept = c(1,-1), colour = "darkblue") + 
  geom_hline(yintercept = c(1,-1), colour = "darkblue") + 
  theme_classic() +
  labs(x=paste0("LogFC QSC_Gli3KO_vs_QSC_control"),
       y=paste0("LogFC QSC Control vs CSC Alert"),
       title = "Differentially Expressed Genes in both contrasts")
