

dds=dds
dfid=IDwgetBM

genelist_mgisymbol = c("Gli1","Gli2","Gli3","Ptch1","Ptch2","Smo","Hhip","Hhat","Spop","Kif7",
                       "Kif3a","Ift88","Sufu","Prkaca","Gpr161","Prkar1a","Pde4b",
                       "Prkacb","Pkia","Prkar2b","Prkar1b","Prkar2a","Akap5","Cul3","Cdk6")
genelist_mgisymbol = genelist_mgisymbol[order(genelist_mgisymbol)]

a = filter(listeResults[[1]],mgi_symbol %in% genelist_mgisymbol)
b = filter(listeResults[[2]],mgi_symbol %in% genelist_mgisymbol)
tmp = filter(a,mgi_symbol %in% genelist_mgisymbol)
tmp = filter(b,mgi_symbol %in% genelist_mgisymbol)

genelist = a$ensembl_gene_id[a$mgi_symbol %in% genelist_mgisymbol]
genelist = b$ensembl_gene_id[b$mgi_symbol %in% genelist_mgisymbol]

  write.csv(listeResults[[1]],file = "./ASC_Gli3KO_vs_ASC_control.csv")
  write.csv(listeResults[[2]],file = "./QSC_Gli3KO_vs_QSC_control.csv")
###### ONLY FOR DFID WHICH NO KEGG
  print("OK create pheatmap")
  countNorm = data.frame(counts(dds, normalized=TRUE))
  countNorm$ensembl_gene_id = rownames(countNorm)
  dfid = unique(dfid[!colnames(dfid) %in% "entrezgene_id"])
  tmp = unique(inner_join(countNorm,dfid,by="ensembl_gene_id"))
  tmp = tmp[tmp$ensembl_gene_id %in% genelist,]


  rownames(tmp) = tmp$mgi_symbol

  
  tmp = tmp[,c(1:dim(countNorm)[2]-1)]
  
  tmp1 = tmp[,c("ASC.540","ASC.543","ASC.928","ASC.246","ASC.960","ASC.1818")]
  tmp2 = tmp[,c("Gli3.2","Gli3.1","Gli3.3","Pax7.CE1","Pax7.CE2","Pax7.CE3")]
  
  mat = as.matrix(tmp1)
  mat = as.matrix(tmp2)
  mat = as.matrix(tmp)
  
  data_subset_norm <- t(apply(mat, 1, cal_z_score))
  data_subset_norm[is.nan(data_subset_norm)] = 0
  
  
  col1 = c("ASC.1818","ASC.928","ASC.540","ASC.246","ASC.960","ASC.543")
  col2 = colnames(tmp2)
  col2 = col2[order(colnames(tmp2),decreasing = T)]
  
  pheatmap_plot = function(data,names_dataset,col) {

    mat = as.matrix(data)
    
    data_subset_norm <- t(apply(mat, 1, cal_z_score))
    data_subset_norm[is.nan(data_subset_norm)] = 0
    
    return(pheatmap::pheatmap(data_subset_norm,
                       cluster_rows = T,
                       cluster_cols = F,
                       annotation_col = mysample,
                       #labels_row = genelist_mgisymbol,
                       labels_col = col,
                       main=paste0("Scalling expressed matrix ",names_dataset)))
  }
  p1=pheatmap_plot(tmp1,names_dataset[1],col1)
  p2=pheatmap_plot(tmp2,names_dataset[2],col2)

  pheatmap::pheatmap(data_subset_norm,
                              cluster_rows = T,
                              cluster_cols = T,
                              annotation_col = mysample,
                              #labels_row = genelist_mgisymbol,
                              #labels_col = col,
                              main=paste0("Scalling expressed matrix ",names_dataset))

  pdf(file="./Heatmap4Caro.pdf")
  print(p1)
  grid.newpage()
  print(p2)
  dev.off()
  
  pdf(file="./GOtermBP4Caro.pdf",width=20,height=10)
  plottingGOSEQ(listeResults[[1]],names_dataset[1],1,"mm9")
  plottingGOSEQ(listeResults[[2]],names_dataset[2],1,"mm9")
  dev.off()
  
  
  garethQSC = read.delim(file = "/Volumes/LaCie/InstitutNeuroMyoGene_DATA/GSE144871/DATA/DATApretreated/Canada/GSE144871_QSC_Gli3_Ctr_vs_QSC_Gli3_KO_190926.txt")
  garethASC = read.delim(file = "/Volumes/LaCie/InstitutNeuroMyoGene_DATA/GSE144871/DATA/DATApretreated/Canada/GSE144871_ASC_Gli3_Ctr_vs_ASC_Gli3_KO_190926.txt")
  
  pdf(file="./GOtermBPCanada4Caro.pdf",width=20,height=10)
  plottingGOSEQ(garethQSC,"QSC_Gli3_Ctr_vs_QSC_Gli3_KO",1,"mm9")
  plottingGOSEQ(garethASC,"ASC_Gli3_Ctr_vs_ASC_Gli3_KO",1,"mm9")
  dev.off()
  
  
  
  
  for (i in 1:length(listeResults)) {
    ########DATABASE : REACTOME
    geneList = listeResults[[i]]$padj
    names(geneList) = listeResults[[i]]$entrezgene_id
    
    tmp = filterdata2pathwayfoundplotting(listeResults[[i]],0.05,1,"mouse")
    print(names_dataset[i])
    
    typeData = c("WITH ALL GENES","WITH UP GENES","WITH DOWN GENES")
    
    pdf(file=paste0("./ANR4Caro_",names_dataset[i],".pdf"))
    for (j in 1:3) {
      if (dim(tmp[[j]])[1] != 0) {
        if (dim(tmp[[j]])[1] > 1) {  
          print(emapplot(tmp[[j]]) + ggtitle(typeData[j])) 
          #print(upsetplot(tmp[[j]]) + ggtitle(typeData[j]))
          #print(dotplot(tmp[[j]]) + ggtitle(typeData[j]))
          #print(heatplot(tmp[[j]], foldChange=geneList) + ggtitle(typeData[j]))
          print(cnetplot(tmp[[j]], node_label="all"))
        }
        #p1 <- cnetplot(tmp[[j]], node_label="category") 
        #p2 <- cnetplot(tmp[[j]], node_label="gene") 
        #print(cowplot::plot_grid(p1, p2, p3, ncol=2, labels=LETTERS[1:3]) + ggtitle(typeData[j]) )
      }
    }
    dev.off()
    
    rm(geneList,tmp)
  }
  
setwd(currentDirectory)
save.image("./Graphs4Caro.Rdata")
load.image("./Graphs4Caro.Rdata")
######################################
library(pathview)
library(gage)
library(gageData)
KEGGPATHWAY=readRDS(file = "/Users/williamjarassier/TRAVAIL/references/Mus_musculus/Ensembl/GRCm38/GMT/Mm.c2.cp.kegg.v7.1.entrez.rds")

KEGGPATHWAYfiltered = KEGGPATHWAY[c("KEGG_INSULIN_SIGNALING_PATHWAY",
                                    "KEGG_CELL_CYCLE",
                                    "KEGG_MTOR_SIGNALING_PATHWAY",
                                    "KEGG_HEDGEHOG_SIGNALING_PATHWAY")]
rm(KEGGPATHWAY)
foldchangesContrast1 = listeResults[[1]]$log2FoldChange
foldchangesContrast2 = listeResults[[2]]$log2FoldChange

names(foldchangesContrast1) = listeResults[[1]]$entrezgene_id
names(foldchangesContrast2) = listeResults[[2]]$entrezgene_id
head(foldchangesContrast1)
# Get the results
keggres1 = gage(foldchangesContrast1, gsets=KEGGPATHWAYfiltered, same.dir=TRUE)
keggres2 = gage(foldchangesContrast2, gsets=KEGGPATHWAYfiltered, same.dir=TRUE)
# Look at both up (greater), down (less), and statatistics.
lapply(keggres, head)

keggrespathways = c("04911","04340","04150","04110")

# Define plotting function for applying later
plot_pathway = function(pid) pathview(gene.data=foldchangesContrast1, pathway.id=pid, species="mmu", new.signature=FALSE)
# plot multiple pathways (plots saved to disk and returns a throwaway list object)
tmp = sapply(keggrespathways, function(pid) pathview(gene.data=foldchangesContrast1, pathway.id=pid, species="mmu"))


# Define plotting function for applying later
plot_pathway = function(pid) pathview(gene.data=foldchangesContrast2, pathway.id=pid, species="mmu", new.signature=FALSE)
# plot multiple pathways (plots saved to disk and returns a throwaway list object)
tmp = sapply(keggrespathways, function(pid) pathview(gene.data=foldchangesContrast2, pathway.id=pid, species="mmu"))


pathview(gene.data=foldchangesContrast2, pathway.id="04340", species="mmu") 
######################################

printVolcanoPlotall_selectGenes = function(data,genelist_mgisymbol,contrastlist) {
res_modified = data.frame(data)

xlimit = c(min(res_modified$log2FoldChange, na.rm=TRUE) - 1,max(res_modified$log2FoldChange, na.rm=TRUE) + 1)
ylimit = c(0, max(-log10(res_modified$padj), na.rm=TRUE) + 5)
library(EnhancedVolcano)
p1=EnhancedVolcano(res_modified,
                      lab = res_modified$mgi_symbol,
                      x = 'log2FoldChange',
                      y = 'padj',
                      FCcutoff = 1,
                      pCutoff = 0.05,
                      xlim = xlimit,
                      ylim = ylimit,
                      subtitle = "All differential expression genes",title ="")
tmp = filter(res_modified,mgi_symbol %in% genelist_mgisymbol)
p2=EnhancedVolcano(tmp,
                lab = tmp$mgi_symbol,
                selectLab = genelist_mgisymbol,
                x = 'log2FoldChange',
                y = 'padj',
                FCcutoff = 1,
                pCutoff = 0.05,
                xlim = xlimit,
                ylim = ylimit,
                subtitle = "Genes selection of Hedgehog signaling pathway",
                title = "")
print(p1)
print(p2)
#print(grid.arrange(grobs = list(p1,p2) ,ncol = 2, top = contrastlist))
}
pdf(file = "./VolcanoPlot_allgenes_and_hedgehog_pathway.pdf")
printVolcanoPlotall_selectGenes(listeResults[[1]],genelist_mgisymbol,"ASC Gli3KO vs ASC Control")
printVolcanoPlotall_selectGenes(listeResults[[2]],genelist_mgisymbol,"QSC Gli3KO vs QSC Control")
dev.off()
################################

load("/Users/williamjarassier/TRAVAIL/references/Mus_musculus/Ensembl/GRCm38/GMT/mouse_H_v5p2.rdata")
HALLMARKgenes = data.frame(extractGeneFromPathway(Mm.H))
colnames(HALLMARKgenes) = "entrezgene_id"
IDwgetBM$entrezgene_id = as.character(IDwgetBM$entrezgene_id)

IDw_HALLMARK = inner_join(IDwgetBM,HALLMARKgenes,by="entrezgene_id")
IDw_HALLMARK$entrezgene_id = as.numeric(IDw_HALLMARK$entrezgene_id)

HALLMARK_genes = extractGeneFromPathway(Mm.H)
IDw_HALLMARK = IDw_pathway(IDwgetBM,HALLMARK_genes,"entrezgene_id")

HALLMARK = list()
HALLMARK[[1]] = plottingFGSEA(listeResults[[1]],IDw_HALLMARK,Mm.H,"noKEGG")
HALLMARK[[2]] = plottingFGSEA(listeResults[[2]],IDw_HALLMARK,Mm.H,"noKEGG")

plotfgsea = function(data) {
  gseaDat = inner_join(data,IDw_HALLMARK,by="ensembl_gene_id")
  gseaDat = filter(gseaDat,!is.na(padj))
  gseaDat$signed = ifelse(gseaDat$log2FoldChange < 0,-1,1)
  #https://github.com/ctlab/fgsea/issues/50
  gseaDat$ranks_signed_padj = - gseaDat$signed * log10(gseaDat$padj)
  gseaDat = gseaDat[order(gseaDat[,"ranks_signed_padj"],decreasing=T),]
  ranks = gseaDat$ranks_signed_padj
  names(ranks) = as.character(gseaDat$entrezgene_id.x)
  
  fgseaRes <- fgsea(Mm.H, ranks, minSize=15, maxSize = 500, nperm=1000)
  topUp <- fgseaRes %>% 
    filter(ES > 0) %>% 
    top_n(10, wt=-padj)
  topDown <- fgseaRes %>% 
    filter(ES < 0) %>% 
    top_n(10, wt=-padj)
  topPathways <- bind_rows(topUp, topDown) %>% 
    arrange(-ES)
  return(plotGseaTable(Mm.H[topPathways$pathway], 
                   ranks, 
                   fgseaRes, 
                   gseaParam = 0.5,
                   render=T))
}

pdf(file = "./fgseaResults.pdf",width = 12)
print(plotfgsea(listeResults[[1]])) #ASC
grid.newpage()
print(plotfgsea(listeResults[[2]])) #QSC
dev.off()



rld<- rlog(dds,blind=TRUE)

config=na.omit(read.delim(file="/Volumes/LaCie/InstitutNeuroMyoGene_DATA/GSE144871/design.tsv",header=T,sep="\t",na.strings=NA,row.names=1))
col="grp"

ntop <- 500
rv <- rowVars(assay(rld))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat <- t( assay(rld)[select, ] )
pc <- prcomp(mat)
eig <- (pc$sdev)^2
variance <- eig*100/sum(eig)
PCAdata<-as.data.frame(pc$x[,1:3])
PCAdata$condition <- config$grp

library(corrplot)
source("http://www.sthda.com/upload/rquery_cormat.r")
rquery.cormat()
rquery.cormat(t(mat), type="upper")
a = rquery.cormat(t(mat), type="flatten", graph=T)[[1]]
PCAdata$row = rownames(PCAdata)
library(dplyr)
b = data.frame(condition1=PCAdata$condition,row=rownames(PCAdata))
c = data.frame(condition2=PCAdata$condition,column=rownames(PCAdata))
d = inner_join(a,b,by="row")
d = inner_join(d,c,by="column")

save = d %>% group_by(condition1,condition2) %>% summarise(Mean=mean(cor), Median=median(cor))

gd <- PCAdata %>% 
  group_by(condition) %>% 
  summarise(PC1 = mean(PC1),
            PC2 = mean(PC2))

n <- nrow(gd) - 1
gd1 <- data.frame(X = c(rep(gd$PC1[1], n), gd$PC1[-1]) ,
                  Y = c(rep(gd$PC2[1], n), gd$PC2[-1]),
                  condition = c(rep(gd$condition[1], n), gd$condition[-1]))
gd1$grp <- as.factor(rep(1:n, times = 1))

gd2 <- data.frame(X = c(rep(gd$PC1[2], n), gd$PC1[-2]) ,
                  Y = c(rep(gd$PC2[2], n), gd$PC2[-2]),
                  condition = c(rep(gd$condition[2], n), gd$condition[-2]))
gd2$grp <- as.factor(rep(1:n, times = 1))

gd3 <- data.frame(X = c(rep(gd$PC1[3], n), gd$PC1[-3]) ,
                  Y = c(rep(gd$PC2[3], n), gd$PC2[-3]),
                  condition = c(rep(gd$condition[3], n), gd$condition[-3]))
gd3$grp <- as.factor(rep(1:n, times = 1))

nameACP1 <- paste(projectName,"_deseq2_Unsupervised_PCA_PC1vsPC2_top500varGenes.jpg",sep="")
ggplot(PCAdata,aes(x=PC1,y=PC2,color = condition)) +
  geom_label_repel(label = rownames(PCAdata)) +
  geom_point(aes(shape=condition, color=condition)) +
  geom_point(data = gd, size = 4) + 
  theme_bw() + 
  guides(color = guide_legend("Condition"),  shape = guide_legend("Condition")) +
  labs(
    x = paste0("PC1: ",round(variance[1],1),"% variance"),
    y = paste0("PC2: ",round(variance[2],1),"% variance") 
  ) +
  geom_line(data = gd1 ,aes(X, Y, group = grp), color = "red", linetype = "dashed") +
  geom_line(data = gd2 ,aes(X, Y, group = grp), color = "blue", linetype = "dashed") +
  geom_line(data = gd3 ,aes(X, Y, group = grp), color = "gray", linetype = "dashed")



pdf("./heatmap4significativesgenes.pdf",width = 10,height = 10)
for (i in 1:2) {
  listeResults[[i]] = resData(dds,unlist(contrastList[i,]),IDwgetBM)
  
  a = filteredDE(listeResults[[i]],0.05,1)

  countNorm = data.frame(counts(dds, normalized=TRUE))
  countNorm$ensembl_gene_id = rownames(countNorm)
  dfid = unique(dfid[!colnames(dfid) %in% "entrezgene_id"])
  tmp = unique(inner_join(countNorm,dfid,by="ensembl_gene_id"))
  tmp = tmp[tmp$mgi_symbol %in% a$mgi_symbol,]
  
  rownames(tmp) = tmp$mgi_symbol
  tmp = tmp[,c(1:dim(countNorm)[2]-1)]
  mat = as.matrix(tmp)
  
  data_subset_norm <- t(apply(mat, 1, cal_z_score))
  data_subset_norm[is.nan(data_subset_norm)] = 0
  
  print(pheatmap::pheatmap(data_subset_norm,
                     annotation_col = mysample,
                     cutree_rows = 2,
                     cutree_cols = 2,
                     main=paste0("Scalling expressed matrix ",names_dataset[i])))
}
dev.off()
