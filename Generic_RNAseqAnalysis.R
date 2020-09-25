####### V2 WILLIAM JARASSIER AOUT 2020
rm(list = ls())
args = commandArgs(trailingOnly=TRUE)

currentDirectory = args[1]
singleAnalysis = args[2]
findPathways = args[3]
genome = args[4]
thresholdFC = args[5]
chemincount = args[6]
pathresults = args[7]
modeCanada = args[8]

source("/Users/williamjarassier/TRAVAIL/scripts/function4RNAseqAnalysis.R")

#currentDirectory = "/Volumes/LaCie/InstitutNeuroMyoGene_DATA/DATA_benedicteDurand"
#singleAnalysis = "yes"
#findPathways = "yes"
#genome = "drosophila"          
#thresholdFC = 2
#chemincount = "Durand_0216-rawcounts.txt"
#pathresults = "resultsRNAseq"
#modeCanada = "no"
currentDirectory = "/Volumes/LaCie/InstitutNeuroMyoGene_DATA/GSE144871"
singleAnalysis = "yes"
findPathways = "yes"
genome = "mouse"          
thresholdFC = 1
chemincount = "/DATA/DATApretreated/Canada/GSE144871_filtered_counts_190926.txt"
pathresults = "resultsRNAseq_vtest"
modeCanada  = "yes"
design = ~grp

setwd(currentDirectory)
dir.create(paste0("./",pathresults))

# METHODE DESEQ2
#####################################################

library(dplyr)
library(biomaRt)
library(DESeq2)
library(ggplot2)
library(VennDiagram)
library(VennDetail)
library(EnhancedVolcano)
library(reshape2)
library(gridExtra)
library(fgsea)
library(grid)
library(pheatmap)
library(stringr)
library(goseq)
library(ReactomePA)
library(clusterProfiler)
library(enrichplot)
library(org.Dm.eg.db)

  
if (singleAnalysis == "yes") {
  #### 
  # rawcountW : merge count gene stringtie from nf-core pipeline
  # design : configuration of this study
  #### 
  
  configuration = na.omit(read.delim(file="./design.tsv",header=T,sep="\t",na.strings=NA,row.names=1))
  namessamples = rownames(configuration)
  rawcountW = readCount(currentDirectory,chemincount,modeCanada,namessamples)

  #### 
  # IDw : table of genes in rawcountW
  # IDwgetBM : table of interested genes from nf-core pipeline
  ###

  IDwgetBM = IDgetBM(rawcountW$rowid,"mouse") 

  ###
  # FCw : filtered count from rawcountW 
  #
  
  FC = filter(rawcountW, rownames(rawcountW) %in% IDwgetBM$ensembl_gene_id) %>% dplyr::select(namessamples)
  
  contrastList = read.delim("./contrastlist.tsv",header=T,sep="\t",na.strings=NA,row.names=1)

  pdf(file = paste0(pathresults,"/dataanalysis4_with_DESEQ2.pdf"),width = 15, height = 10)
  AnalyzingDESeq2(FC,configuration,~grp,"grp",IDwgetBM,contrastList,genome)
  dev.off()
}

if (findPathways == "yes") {
  
  ###
  # Calcul du DDS + génération des results ASC et QSC
  
  #dds = createDDS(FC,configuration,~grp) #15457

  listeResultsFiltered = list()
  for (i in 1:dim(contrastList)[1]) {
    listeResultsFiltered[[i]] = filteredDE(listeResults[[i]],0.05,1)
  }
  names(listeResultsFiltered) = paste0(contrastList[,2], " VS ", contrastList[,3])
  
  #combinaison
  combinaison = combn(dim(contrastList)[1],2)
  
  pdf(paste0(pathresults,"/compare_contrasts_by_FC.pdf"),width = 10, height = 10)
  for (i in 1:dim(combinaison)[2]) {
    x = listeResultsFiltered[[combinaison[1,i]]]
    y = listeResultsFiltered[[combinaison[2,i]]]
  
    communDE(x,y,genome,names(listeResultsFiltered),thresholdFC)
    rm(x,y)
  }
  dev.off()
  
  rm(combinaison)
  
  pdf(paste0(pathresults,"/DE4goodGenes.pdf"),width = 10, height = 70)
  for (i in 1:length(listeResultsFiltered)) {
    genelist = listeResultsFiltered[[i]]$ensembl_gene_id
    pheatmap4genesOK(dds,IDwgetBM,genelist,configuration,genome,names(listeResultsFiltered)[i])
  }
  dev.off()
  
  data = inner_join(listeResults[[1]],listeResults[[2]],by="ensembl_gene_id") %>%
    filter(padj.x < 0.05 & padj.y < 0.05) %>%
    dplyr::select(c("log2FoldChange.x","log2FoldChange.y","mgi_symbol.x"))
  rownames(data) = data$mgi_symbol.x
  colnames(data) = names(listeResultsFiltered)
  data = data[,c(1:2)]
  
  mat = as.matrix(data)

  pheatmap::pheatmap(mat)
  
  ###HEDGEHOG
  #genelist_mgisymbol = c("Gli1","Gli2","Gli3","Ptch1","Ptch2","Smo","Hhip","Hhat","Spop","Kif7",
  #             "Kif3a","Ift88","Sufu","Prkaca","Gpr161","Prkar1a","Adcy8","Pde4b","Prkaca",
  #             "Prkacb","Pkia","Prkar2b","Prkar1b","Prkar2a","Akap5","Cul3","Cdk6")
  #genelist = listeResults[[1]]$ensembl_gene_id[listeResults[[1]]$mgi_symbol %in% genelist_mgisymbol]
  #
  #pheatmap4genesOK(dds,IDwgetBM,genelist,configuration,genome,names_dataset[i])

  ###
  # LOADING GMT DATABASE
  #
  #
  
  if (genome == "mouse") {
    
    print("LOADING OF GMT FILES")
    
    tmp <- PGSEA::readGmt("/Users/williamjarassier/TRAVAIL/references/Mus_musculus/Ensembl/GRCm38/GMT/MousePath_Pathway.gmt")
    pathwayM_GSKB = list()
    names = c()
    for (i in 2:length(tmp)) {
      pathwayM_GSKB[i-1] = list(stringr::str_to_upper(tmp[[i]]@ids))
      names[i-1] = tmp[[i]]@reference
    }
    names(pathwayM_GSKB) = names
    
    pathwayM_GSKB_MOUSECYC = pathwayM_GSKB[grep("^MOUSECYC_MM",names(pathwayM_GSKB))]
    pathwayM_GSKB_INOH = pathwayM_GSKB[grep("^INOH",names(pathwayM_GSKB))]
    pathwayM_GKSB_WIKIPATHWAYS = pathwayM_GSKB[grep("^WIKIPATHWAY",names(pathwayM_GSKB))]
    #ONLY KEGG
    KEGGPATHWAY=readRDS(file = "/Users/williamjarassier/TRAVAIL/references/Mus_musculus/Ensembl/GRCm38/GMT/Mm.c2.cp.kegg.v7.1.entrez.rds")
    
    rm(pathwayM_GSKB,names,tmp)
    
    MOUSECYC_genes = extractGeneFromPathway(pathwayM_GSKB_MOUSECYC)
    INOH_genes = extractGeneFromPathway(pathwayM_GSKB_INOH)
    WIKIPATHWAY_genes = extractGeneFromPathway(pathwayM_GKSB_WIKIPATHWAYS)
    KEGG_genes = extractGeneFromPathway(KEGGPATHWAY)
    
    #########GLOSSAIRE AVEC IDw
    IDw_KEGG = IDw_pathway(IDwgetBM,KEGG_genes,"entrezgene_id")
    IDw_MOUSECYC = IDw_pathway(IDwgetBM,MOUSECYC_genes,"mgi_symbol")
    IDw_INOH = IDw_pathway(IDwgetBM,INOH_genes,"mgi_symbol")
    IDw_WIKIPATHWAY = IDw_pathway(IDwgetBM,WIKIPATHWAY_genes,"mgi_symbol")
    
    rm(KEGG_genes,MOUSECYC_genes,INOH_genes,WIKIPATHWAY_genes)

    KEGG = lapply(listeResults,function(x) {plottingFGSEA(x,IDw_KEGG,KEGGPATHWAY,"KEGG")})
    MOUSECYC = lapply(listeResults,function(x) {plottingFGSEA(x,IDw_MOUSECYC,pathwayM_GSKB_MOUSECYC,"noKEGG")})
    INOH = lapply(listeResults,function(x) {plottingFGSEA(x,IDw_INOH,pathwayM_GSKB_INOH,"noKEGG")})
    WIKIPATHWAY = lapply(listeResults,function(x) {plottingFGSEA(x,IDw_WIKIPATHWAY,pathwayM_GKSB_WIKIPATHWAYS,"noKEGG")})
    
    ####PLOT FGSEA PLOT DESACTIVE
    #fgseaPDF = lapply(names(listeResults), function(x) paste0(paste0(pathresults,
    #                                                                 "/FGSEA_",
    #                                                                 stringr::str_replace_all(string = x ,pattern = " ",replacement = "_"),
    #                                                                 ".pdf")))
    #filespdfPATHWAY = lapply(names(listeResults), function(x) paste0(paste0("/ALLPATHWAY_FGSEA_",
    #                                                                        stringr::str_replace_all(string = x ,pattern = " ",replacement = "_"),
    #                                                                        ".pdf")))
  
    rm(fgseaPDF)
    
    tmp = data.frame(t(contrastList[,c(2:3)]))
    listeSamples = lapply(tmp, function(x) {rownames(configuration)[configuration$grp %in% unlist(t(x))]}) 
    rm(tmp)
    names(listeSamples) = names(listeResults)
    
    listDir = c(paste0(pathresults,"/DEmatrix2Pathway_KEGG"),
                paste0(pathresults,"/DEmatrix2Pathway_INOH"),
                paste0(pathresults,"/DEmatrix2Pathway_MOUSECYC"),
                paste0(pathresults,"/DEmatrix2Pathway_WIKIPATHWAY"))
    
    fgseaResultsList = list(KEGG,INOH,MOUSECYC,WIKIPATHWAY)
    dbIDList = list(IDw_KEGG,IDw_INOH,IDw_MOUSECYC,IDw_WIKIPATHWAY)
    
    for (i in 1:length(listDir)) {
      dir.create(listDir[i])
      
      fgseaResults = fgseaResultsList[[i]]
      
      typeDB = ifelse(i==1,"KEGG","noKEGG")
      dbID = dbIDList[[i]]
      
      for (j in 1:length(filespdfPATHWAY)) {
        path = paste0(listDir[i],filespdfPATHWAY[j])  
        Samples = unlist(listeSamples[j])
        pdf(path)
        pheatmapPathwayGenesExpressed(fgseaResults[[j]],dds,Samples,typeDB,dbID)
        dev.off()
      }
    }
    
    rm(listDir,filespdfPATHWAY,fgseaResultsList,dbIDList)
    
    #####################################################
    pdf(file=paste0(pathresults,"/GOTERM_ANALYSIS_RESULTS.pdf"),width=20,height=10)
    for (i in 1:length(listeResults)) {
      plottingGOSEQ(listeResults[[i]],names_dataset[i],thresholdFC,"mm9")
    }
    dev.off()
    #####################################################
    for (i in 1:length(listeResults)) {
      ########DATABASE : REACTOME
      geneList = listeResults[[i]]$padj
      names(geneList) = listeResults[[i]]$entrezgene_id
      
      tmp = filterdata2pathwayfoundplotting(listeResults[[i]],0.05,1,"mouse")

      typeData = c("WITH ALL GENES","WITH UP GENES","WITH DOWN GENES")
      
      pdf(file=paste0(pathresults,"/ANALYSIS_NETWORK_REACTOME_",names_dataset[i],".pdf"),width=20,height=10)
        for (j in 1:3) {
          if (dim(tmp[[j]])[1] != 0) {
            if (dim(tmp[[j]])[1] > 1) {  
              print(emapplot(tmp[[j]]) + ggtitle(typeData[j])) 
              print(upsetplot(tmp[[j]]) + ggtitle(typeData[j]))
              print(dotplot(tmp[[j]]) + ggtitle(typeData[j]))
              print(heatplot(tmp[[j]], foldChange=geneList) + ggtitle(typeData[j]))
            }
            p1 <- cnetplot(tmp[[j]], node_label="category") 
            p2 <- cnetplot(tmp[[j]], node_label="gene") 
            p3 <- cnetplot(tmp[[j]], node_label="all") 
            print(cowplot::plot_grid(p1, p2, p3, ncol=2, labels=LETTERS[1:3]) + ggtitle(typeData[j]) )
            rm(p1,p2,p3)
          }
        }
      dev.off()
      
      rm(geneList,tmp)
    }
  } 
  else {
    #library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
    #####################################################
    pdf(file=paste0(pathresults,"/GOTERM_ANALYSIS_RESULTS.pdf"),width=20,height=10)
    for (i in 1:length(listeResults)) {
      plottingGOSEQ(listeResults[[i]],names_dataset[i],thresholdFC,"dm3")
    }
    dev.off()
    #####################################################
    for (i in 1:length(listeResults)) {
      ########DATABASE : REACTOME
      geneList = listeResults[[i]]$padj
      names(geneList) = listeResults[[i]]$entrezgene_id
      
      tmp = filterdata2pathwayfoundplotting(listeResults[[i]],0.05,1,"fly")
      print(names_dataset[i])
      
      typeData = c("WITH ALL GENES","WITH UP GENES","WITH DOWN GENES")
      
      pdf(file=paste0(pathresults,"/ANALYSIS_NETWORK_REACTOME_",names_dataset[i],".pdf"),width=20,height=10)
      for (j in 1:3) {
        if (dim(tmp[[j]])[1] != 0) {
          if (dim(tmp[[j]])[1] > 1) {  
            print(emapplot(tmp[[j]]) + ggtitle(typeData[j])) 
            print(upsetplot(tmp[[j]]) + ggtitle(typeData[j]))
            print(dotplot(tmp[[j]]) + ggtitle(typeData[j]))
            print(heatplot(tmp[[j]], foldChange=geneList) + ggtitle(typeData[j]))
          }
          p1 <- cnetplot(tmp[[j]], node_label="category") 
          p2 <- cnetplot(tmp[[j]], node_label="gene") 
          p3 <- cnetplot(tmp[[j]], node_label="all") 
          print(cowplot::plot_grid(p1, p2, p3, ncol=2, labels=LETTERS[1:3]) + ggtitle(typeData[j]) )
          rm(p1,p2,p3)
        }
      }
      dev.off()
      
      rm(geneList,tmp)
    }
  }
}
