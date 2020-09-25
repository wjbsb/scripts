################FUNCTIONS



readCount = function(racine,chemincount,modeCanada,namessamples) {
  setwd(racine)
  chemincount = paste0("./",chemincount)
  
  tmp = read.csv(chemincount,
                 header=T,
                 dec=".",
                 sep="\t",
                 na.strings=NA,
                 row.names=1) %>% dplyr::select(namessamples)

  if (modeCanada == "no") {
    tmp = tmp %>% mutate(rowid = rownames(tmp))
  } else {
    tmp = data.frame(mapply(round,tmp)) %>% mutate(rowid = stringr::word(rownames(tmp),1,sep = "\\."))
  }
  
  rownames(tmp) = tmp$rowid
  return(tmp)
}
modifID = function(dfid) {
  dfid$mgi_symbol[dfid$ensembl_gene_id == "ENSMUSG00000095675"] = "Ccl21b"
  dfid$mgi_symbol[dfid$ensembl_gene_id == "ENSMUSG00000094121"] = "LOC100041504"
  dfid$mgi_symbol[dfid$ensembl_gene_id == "ENSMUSG00000096271"] = "LOC100041593"
  dfid$mgi_symbol[dfid$ensembl_gene_id == "ENSMUSG00000096873"] = "Ccl21c"
  
  
  dfid$mgi_symbol[dfid$ensembl_gene_id == "ENSMUSG00000038729"] = "Pakap_201_210"
  dfid$mgi_symbol[dfid$ensembl_gene_id == "ENSMUSG00000089945"] = "Pakap_211_213"
  dfid$mgi_symbol[dfid$ensembl_gene_id == "ENSMUSG00000090053"] = "Pakap_214_216"
  
  
  return(dfid)
}
IDgetBM = function(IDs,genome) {
  
  ensembl <- useMart("ensembl")
  
  if (genome == "mouse") {
    ensembl <- useDataset("mmusculus_gene_ensembl",mart=ensembl)
    keep_attributes = c('ensembl_gene_id','mgi_symbol','gene_biotype','entrezgene_id')
  } else if (genome == "drosophila") {
    ensembl = useDataset("dmelanogaster_gene_ensembl",mart=ensembl)
    keep_attributes = c('ensembl_gene_id','flybasename_gene','gene_biotype','external_gene_name')
  }
  
  tmp = getBM(attributes = keep_attributes, 
              filters = 'ensembl_gene_id', 
              values = IDs, 
              mart = ensembl)
  tmp = modifID(tmp)
  return(filter(tmp,gene_biotype == "protein_coding"))
}
make_nice_clusters <- function(rld,configuration,col,projectName) {

  matDist <- as.matrix(dist(t(assay(rld))))
  
  Conditions <- data.frame(condition = configuration[,col],row.names=rownames(configuration))
  nameClustering <- paste(projectName,"deseq2_Unsupervised_clustering_euclidean-complete.png",sep="_")
  return(pheatmap::pheatmap(matDist, 
                            col=colorRampPalette(RColorBrewer::brewer.pal(9, 'GnBu'))(100), 
                            annotation_col = Conditions,
                            show_rownames=T,
                            show_colnames=T,
                            main=paste0("Unsupervised Clustering Euclidean (Meth.DESeq2)",projectName)))
}
make_nice_pca <- function(rld,configuration,ntop,thresholdimension) {

  if (is.na(ntop) | is.na(thresholdimension)) {
    ntop <- 500
    thresholdimension = 3
  }

  rv <- rowVars(assay(rld))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  mat <- t(assay(rld)[select, ])
  res.pca <- prcomp(mat,scale=F)

  eig.val <- get_eigenvalue(res.pca)[c(1:thresholdimension),]
  #Montre le pourcentage de variances expliquées par chaque axe principal.
  contributionDimension = fviz_eig(res.pca,addlabels = TRUE, ylim = c(0, 50))
  #Coloration en fonction du cos2 (qualité de représentation). Les individus similaires sont groupés ensemble.
  listDimension = list(c(1,2),c(1,3),c(2,3))
  
  PCAcos2 = lapply(listDimension, function(x) {fviz_pca_ind(res.pca, axes = x,
                                                        col.ind = "cos2", # Colorer par le cos2
                                                        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                        repel = TRUE,
                                                        title = ""
                                                        )})
  names(PCAcos2) = c("Dim1 vs Dim2","Dim1 vs Dim3","Dim2 vs Dim3")
  PCAcos2[[4]] = contributionDimension
  
  listDimensionUnique = c(1:3)
  contributionIndividu2Dimension = lapply(listDimensionUnique, function(x) {fviz_contrib(res.pca, choice = "ind", axes = x)})
  names(contributionIndividu2Dimension) = c("Dim1","Dim2","Dim3")
  
  # Contribution totale sur PC1 et PC2
  contributionDimensionTotal = fviz_contrib(res.pca, choice = "ind", axes = 1:thresholdimension)
  
  contributionIndividu2Dimension[[4]] = contributionDimensionTotal
  
  PCAgrp = lapply(listDimension, function(x) {fviz_pca_ind(res.pca, axe = x,
                                                           geom.var = c("point","text"),
                                                           col.ind = configuration$grp,
                                                           legend.title = "Groups",
                                                           title = ""
                                              )})
  names(PCAgrp) = c("Dim1 vs Dim2","Dim1 vs Dim3","Dim2 vs Dim3")
  

  library(ggpubr)
  
  print(ggarrange(
    plotlist = PCAcos2, labels = c("A", "B","C","D"),
    common.legend = TRUE, legend = "bottom"
  ))
  print(ggarrange(
    plotlist = PCAgrp, labels = c("A", "B","C"),
    common.legend = TRUE, legend = "bottom"
  ))
  print(ggarrange(
    plotlist = contributionIndividu2Dimension, labels = c("A", "B","C","D"),
    common.legend = TRUE, legend = "bottom"
  ))
}
AnalyzingDESeq2 = function(count,configuration,design,col,dfID,contrastList,genome) {
  
  ### CREATE DESEQ OBJECT
  dds = createDDS(count,configuration,design) #15457

  print("Dispersion")
  plotDispEsts(dds, ylim = c(1e-6, 1e1))
  print("QC")
  
  configuration$Name = rownames(configuration)
  
  countNorm<-counts(dds, normalized=TRUE)
  pseudoNormCount <- data.frame(log2(countNorm+1),
                                row.names=row.names(countNorm),
                                Ids = row.names(countNorm))
  datNorm <- melt(pseudoNormCount, 
                  id.vars = "Ids", 
                  variable.name = "Samples", 
                  value.name = "count") %>% merge(configuration, by.x = "Samples", by.y = "Name")
  
  countRaw<-counts(dds, normalized=FALSE)
  pseudoRawCount <- data.frame(log2(countRaw+1),
                               row.names=row.names(countRaw),
                               Ids = row.names(countRaw))
  
  datRaw <- melt(pseudoRawCount, 
                 id.vars = "Ids", 
                 variable.name = "Samples", 
                 value.name = "count") %>% merge(configuration, by.x = "Samples", by.y = "Name")

  gDens <- ggplot(datNorm, aes(x = count, colour = Samples, fill = Samples)) + 
    geom_density(alpha = 0.2, size = 1.25) + 
    facet_wrap(~ Condition) + theme(legend.position = "top") + 
    xlab(expression(log[2](count + 1)))
  gNorm <- ggplot(data = datNorm, aes(x = Samples, y = count, fill = datNorm[,col])) + 
    geom_boxplot() + xlab("Samples") + ylab("log2(Count+1)") + 
    ggtitle("Normalized log2(counts) per sample") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(fill = col)
  gRaw <- ggplot(data = datRaw, aes(x = Samples, y = count, fill = datNorm[,col])) + 
    geom_boxplot() + xlab("Samples") + ylab("log2(Count+1)") + 
    ggtitle("Raw log2(counts) per sample") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(fill = col)

  ggarrange(gDens,gNorm,gRaw, labels = c("A", "B","C"), common.legend = TRUE, legend = "bottom")
  
  rm(countNorm,countRaw,pseudoNormCount,pseudoRawCount,datNorm,datRaw,gDens,gNorm,gRaw)
  
  rld<- rlog(dds,blind=TRUE)
  print("MAKE NICE PLOT")
  make_nice_clusters(rld,configuration,col,"GSE144871")
  print("MAKE RESULTS + plotMA")
  make_nice_pca(rld,configuration,NA,NA)

  
  listeResults = list()
  for (i in 1:dim(contrastList)[1]) {
    listeResults[[i]] = resData(dds,unlist(contrastList[i,]),IDwgetBM)
    names(listeResults)[i] = paste0(contrastList[,2], " VS ", contrastList[,3]) 
   if (genome == "mouse") {
     print(EnhancedVolcano(listeResults[[i]],
                           lab = listeResults[[i]]$mgi_symbol,
                           x = 'log2FoldChange',
                           y = 'padj',
                           FCcutoff = 2,
                           pCutoff = 0.05,
                           title = names(listeResults)[i]))
    } else {
      print(EnhancedVolcano(listeResults[[i]],
                           lab = listeResults[[i]]$external_gene_name,
                           x = 'log2FoldChange',
                           y = 'padj',
                           FCcutoff = 2,
                           pCutoff = 0.05,
                           title = names(listeResults)[i]))
    }
  }
  
  
  
  save(dds,listeResults, file = paste0("./dds_resDESEQ2.Rdata"))
  
  rm(rld)
}
communDE = function(data1,data2,genome,names_xy,thresholdFC) {
  
  totalfilteredAQ = full_join(data1,data2,by="ensembl_gene_id")
  
  if (genome == "mouse") {
    datalist = list(na.omit(totalfilteredAQ$mgi_symbol.x),na.omit(totalfilteredAQ$mgi_symbol.y))
    names(datalist) = names_xy
    #ven <- venndetail(datalist)
    commun = filter(totalfilteredAQ,!is.na(mgi_symbol.x) & !is.na(mgi_symbol.y)) #2  proteincoding
  } else if (genome == "drosophila") {
    datalist = list(na.omit(totalfilteredAQ$external_gene_name.x),na.omit(totalfilteredAQ$external_gene_name.y))
    names(datalist) = names_xy
    #ven <- venndetail(datalist)
    commun = filter(totalfilteredAQ,!is.na(external_gene_name.x) & !is.na(external_gene_name.y)) #2  proteincoding
  }
  ######ERROR ARGUMENT INUTILISE
  #plot(ven, type ="ven")
  
  commun_filtered_1= filter(commun,(log2FoldChange.x > thresholdFC & log2FoldChange.y > thresholdFC))
  commun_filtered_2= filter(commun,(log2FoldChange.x < thresholdFC & log2FoldChange.y < thresholdFC))
  different_1 = filter(commun,(log2FoldChange.x < thresholdFC & log2FoldChange.y > thresholdFC))
  different_2 = filter(commun,(log2FoldChange.x > thresholdFC & log2FoldChange.y < thresholdFC))
  
  if (genome == "mouse") {
    print(ggplot(commun,aes(log2FoldChange.x,log2FoldChange.y)) + 
            geom_point() + 
            ggrepel::geom_text_repel(data=commun_filtered_1,
                                     aes(x=log2FoldChange.x,
                                         y=log2FoldChange.y,
                                         label=mgi_symbol.x),
                                     color="darkgreen") +
            ggrepel::geom_text_repel(data=commun_filtered_2,
                                     aes(x=log2FoldChange.x,
                                         y=log2FoldChange.y,
                                         label=mgi_symbol.x),
                                     color="darkred") +
            ggrepel::geom_text_repel(data=different_1,
                                     aes(x=log2FoldChange.x,
                                         y=log2FoldChange.y,
                                         label=mgi_symbol.x),
                                     color="darkblue") +
            ggrepel::geom_text_repel(data=different_2,
                                     aes(x=log2FoldChange.x,
                                         y=log2FoldChange.y,
                                         label=mgi_symbol.x),
                                     color="black") +
            geom_abline(intercept=0,slope = 1) + 
            geom_vline(xintercept = c(1,-1), colour = "darkblue") + 
            geom_hline(yintercept = c(1,-1), colour = "darkblue") + 
            theme_classic() +
            labs(x=paste0("LogFC ", names_xy[1]),
                 y=paste0("LogFC ", names_xy[2]),
                 title = "Differentially Expressed Genes in both contrasts"))
  } else if (genome == "drosophila") {
    print(ggplot(commun,aes(log2FoldChange.x,log2FoldChange.y)) + 
            geom_point() + 
            ggrepel::geom_text_repel(data=commun_filtered_1,
                                     aes(x=log2FoldChange.x,
                                         y=log2FoldChange.y,
                                         label=external_gene_name.x),
                                     color="darkgreen") +
            ggrepel::geom_text_repel(data=commun_filtered_2,
                                     aes(x=log2FoldChange.x,
                                         y=log2FoldChange.y,
                                         label=external_gene_name.x),
                                     color="darkred") +
            ggrepel::geom_text_repel(data=different_1,
                                     aes(x=log2FoldChange.x,
                                         y=log2FoldChange.y,
                                         label=external_gene_name.x),
                                     color="darkblue") +
            ggrepel::geom_text_repel(data=different_2,
                                     aes(x=log2FoldChange.x,
                                         y=log2FoldChange.y,
                                         label=external_gene_name.x),
                                     color="black") +
            geom_abline(intercept=0,slope = 1) + 
            geom_vline(xintercept = c(1,-1), colour = "darkblue") + 
            geom_hline(yintercept = c(1,-1), colour = "darkblue") + 
            theme_classic() +
            labs(x=paste0("LogFC ", names_xy[1]),
                 y=paste0("LogFC ", names_xy[2]),
                 title = "Differentially Expressed Genes in both contrasts"))
  }
  rm(totalfilteredAQ,commun,commun_filtered_1,commun_filtered_2)
}
createDDS = function(count,configuration,design) {
  dds = DESeqDataSetFromMatrix(countData = count,
                               colData = configuration,
                               design = design) %>% estimateSizeFactors()
  ### FILTER DESEQ OBJECT
  idx <- rowSums(counts(dds, normalized=TRUE) >= 10 ) >= 3
  return(dds[idx,] %>% DESeq())
}
resData = function(dds,contrastList,dfID) {
  res = results(dds,contrast = contrastList)

  identifiant = filter(dfID, ensembl_gene_id %in% rownames(res))

  return(data.frame(res@listData) %>% 
           mutate(ensembl_gene_id = res@rownames) %>%
           inner_join(identifiant,by="ensembl_gene_id"))
}
filteredDE = function(data,padjvallue,logFCvalue) {
  return(
    filter(data,padj < padjvallue & abs(log2FoldChange) > logFCvalue)
  )
}
pheatmap4genesOK = function(dds,dfid,genelist,configuration,genome,names_dataset) {
  
  ###### ONLY FOR DFID WHICH NO KEGG
  print("OK create pheatmap")
  countNorm = data.frame(counts(dds, normalized=TRUE)) %>% mutate(ensembl_gene_id = rownames(countNorm))
  
  dfid = unique(dfid[!colnames(dfid) %in% "entrezgene_id"])
  tmp = unique(inner_join(countNorm,dfid,by="ensembl_gene_id")) %>% filter(ensembl_gene_id %in% genelist)

  print(paste0("DATA GRAPH = ",dim(tmp)[1]))
  if (dim(tmp)[1] > 0) {
    if (genome == "mouse") {
      rownames(tmp) = tmp$mgi_symbol
    } else if (genome == "drosophila") {
      rownames(tmp) = tmp$external_gene_name
    }
    tmp = tmp[,c(1:dim(countNorm)[2]-1)]
    
    mat = as.matrix(tmp)
    
    mysample = data.frame(configuration[,c("Condition","Genotype")])
    rownames(mysample) = rownames(configuration)
  
    pheatmap::pheatmap(mat,
                       annotation_col = mysample,
                       main=paste0("Expressed Matrix from Normalized Count ",names_dataset),
                       cutree_rows = 2,
                       cutree_cols = 2)
    
    log2mat = log2(mat)
    log2mat[is.nan(log2mat)] = 0
    log2mat[is.infinite(log2mat)] = 0
    
    pheatmap::pheatmap(log2mat,
                       annotation_col = mysample,
                       main=paste0("Expressed Matrix with Log2 regression of Normalized Count ",names_dataset),
                       cutree_rows = 2,
                       cutree_cols = 2)
    data_subset_norm <- t(apply(mat, 1, cal_z_score))
    data_subset_norm[is.nan(data_subset_norm)] = 0
    
    pheatmap::pheatmap(data_subset_norm,
                       annotation_col = mysample,
                       cutree_rows = 2,
                       cutree_cols = 2,
                       main=paste0("Scalling expressed matrix ",names_dataset))
  } else {
    print("NO DATA")
  }
}
extractGeneFromPathway = function(dfpathway_gmt) {
  genelist = data.frame()
  for (i in 1:length(dfpathway_gmt)) {
    if (i == 1) {
      print("0%")
    } else if ((i/length(dfpathway_gmt)*100) %% 10 ==0) {
      print(paste0(i/length(dfpathway_gmt)*100,"%"))
    }

    tmp = data.frame(genes=dfpathway_gmt[[i]])
    genelist = rbind(genelist,tmp)
  }
  return(unique(genelist))
}
IDw_pathway = function(dfid,dfpathway,ID) {
  
  dfid$mgi_symbol = stringr::str_to_upper(dfid$mgi_symbol)
  if (ID == "entrezgene_id") {
    return(dfpathway %>% 
             mutate(genes = as.numeric(genes)) %>% 
             inner_join(dfid,by = c("genes"="entrezgene_id")))
    
  } else if (ID == "external_gene_name") {
    return(dfpathway %>% 
             mutate(genes = stringr::str_to_upper(as.character(genes))) %>% 
             inner_join(dfid,by = c("genes"="external_gene_name")))
  } else  {
    return(dfpathway %>% 
             mutate(genes = as.character(genes)) %>% 
             inner_join(dfid,by = c("genes"="mgi_symbol")))
  }
}
plottingFGSEA = function(data,genelist,databasepathway,typeDB) {
  #IDw_MOUSECYC,pathwayM_GSKB_MOUSECYC
  gseaDat = inner_join(data,genelist,by="ensembl_gene_id") %>% 
    filter(!is.na(padj)) %>%
    mutate(signed = ifelse(log2FoldChange < 0,-1,1),
           ranks_signed_padj = - signed * log10(padj))

  #https://github.com/ctlab/fgsea/issues/50
  
  gseaDat = gseaDat[order(gseaDat[,"ranks_signed_padj"],decreasing=T),]
  
  ranks = gseaDat$ranks_signed_padj
  if (typeDB == "noKEGG") {
    names(ranks) = stringr::str_to_upper(gseaDat$genes)
  } else if (typeDB == "fly") {
    names(ranks) = as.character(gseaDat$external_gene_name.x)
  } else {
    names(ranks) = as.character(gseaDat$entrezgene_id)
  }
  
  #fgseaRes <- fgsea(databasepathway, ranks, minSize=15, maxSize = 500, nperm=1000)
  fgseaRes <- fgsea(databasepathway, ranks, minSize=15, maxSize = 500)
  
  topUp <- fgseaRes %>% 
    filter(ES > 0) %>% 
    top_n(10, wt=-padj)
  topDown <- fgseaRes %>% 
    filter(ES < 0) %>% 
    top_n(10, wt=-padj)
  topPathways <- bind_rows(topUp, topDown) %>% 
    arrange(-ES)
  
  rm(gseaDat,topUp,topDown)
  
  #PGT = list(databasepathway[topPathways$pathway],
  #           ranks,
  #           fgseaRes,
  #           0.5,
  #           F)
  #names(PGT)= c("pathways","stats","fgseaRes","gseaParam","render")
  
  p = plotGseaTable(databasepathway[topPathways$pathway], 
                      ranks, 
                      fgseaRes, 
                      gseaParam = 0.5,
                      render=F)
  d = filter(fgseaRes,pathway %in% topPathways$pathway)
  
  object = list(p,d)
  names(object) = c("plotGseaTable","DATA")
  return(object)
}
pheatmapPathwayGenesExpressed = function(fgseaResults,dds,Samples,typeDB,dbID) {
  
  data = fgseaResults[,c("pathway","leadingEdge")]
  countNorm = counts(dds, normalized=TRUE)
  
  for (i in 1:dim(data)[1]) {
    
    tmp = data.frame(sapply(data[i,"leadingEdge"],c))
    if (dim(tmp)[1] >= 2) {
      print(as.character(data[i,"pathway"]))
      #FOR OTHERS DATABASE
      #colnames(tmp) = "mgi_symbol"
      #IDsWithNamesDesc$mgi_symbol = stringr::str_to_upper(IDsWithNamesDesc$mgi_symbol)
      #tmp = inner_join(tmp,IDsWithNamesDesc,by="mgi_symbol")
      #suppression de entrezgene_id pour eviter les doublons de genes
      #tmp = unique(tmp[,1:3])
      #pour eviter les ensembl_gene_id multiples et parce que la table des counts est ensembl ID like
      #rownames(tmp) = tmp$ensembl_gene_id
      #mat = as.array(CountNormMatrix[rownames(CountNormMatrix) %in% rownames(tmp) ,Samples])
      #rownames(mat) = tmp[rownames(tmp) %in% rownames(mat),"mgi_symbol"]
      
      #FOR KEGGPATHWAY
      if (typeDB == "KEGG") { 
        colnames(tmp) = "entrezgene_id"
        dbID$entrezgene_id = as.character(dbID$genes)
        tmp = inner_join(tmp,dbID,by="entrezgene_id")
        
        tmp = unique(tmp[,c("ensembl_gene_id","mgi_symbol")])
      } else {
        colnames(tmp) = "mgi_symbol"
        dbID$entrezgene_id = stringr::str_to_upper(as.character(dbID$entrezgene_id))
        tmp = inner_join(tmp,dbID,by="mgi_symbol")
        tmp = unique(tmp[,c("ensembl_gene_id","mgi_symbol")])
      }
      
      if (dim(tmp)[1] >= 2) {
        
        #la table des counts est ensembl ID like
        rownames(tmp) = tmp$ensembl_gene_id
        
        mat = countNorm[rownames(countNorm) %in% rownames(tmp),Samples]
        rownames(mat) = tmp[rownames(tmp) %in% rownames(mat),"mgi_symbol"]
        
        #afficher matrice expression
        mysample = data.frame(sample=configuration[Samples,"Condition"])
        rownames(mysample) = rownames(configuration[Samples,])
        
        #pheatmap::pheatmap(mat,
        #                   annotation_col = mysample,
        #                   main=paste0("Expressed Matrix from Normalized for genes in \n",data[i,"pathway"]),
        #                   cutree_rows = 2,
        #                   cutree_cols = 2)
        log2mat = log2(mat)
        log2mat[is.nan(log2mat)] = 0
        log2mat[is.infinite(log2mat)] = 0
        
        pheatmap::pheatmap(log2mat,
                         annotation_col = mysample,
                           main=paste0("Expressed Matrix for genes active in \n",data[i,"pathway"]),
                           cutree_rows = 2,
                           cutree_cols = 2)
        #data_subset_norm <- t(apply(mat, 1, cal_z_score))
        #data_subset_norm[is.nan(data_subset_norm)] = 0
        #
        #pheatmap::pheatmap(data_subset_norm,
        #                   annotation_col = mysample,
        #                   cutree_rows = 2,
        #                   cutree_cols = 2,
        #                   main=paste0("Scalling expressed matrix for genes active in \n",data[i,"pathway"]))
      } else {print(paste0(data[i,"pathway"]," ===> ",unlist(data[i,"leadingEdge"])))} 
    } else {print(paste0(data[i,"pathway"]," ===> ",unlist(data[i,"leadingEdge"])))}
  }
}
plotgoseq = function(data,filtervalue) {
  return(data %>% 
    filter(ontology %in% filtervalue) %>%
    top_n(10, wt=-over_represented_pvalue) %>%
    ggplot(aes(x=hitsPerc, 
               y=term, 
               colour=over_represented_pvalue, 
               size=numDEInCat)) +
    geom_point() +
    expand_limits(x=0) +
    labs(x="Hits (%)", y="GO term", colour="p value", size="Count", title=filtervalue))
}
calculateGOSEQ = function(data,TITLE,genome) {
  #mouse = mm9
  #drosophila = TxDb.Dmelanogaster.UCSC.dm6.ensGene
  idnullp = ifelse(genome == "mm9","ensGene","geneSymbol")
  
  #txdb.ucsc <- makeTxDbFromUCSC("dm6","ensGene")
  #genes(txdb.ucsc)
  
  pwf <- nullp(DEgenes = data,
               genome = genome,
               id = "ensGene",
               bias.data = NULL,
               plot.fit = F)
  goResults <- goseq(pwf, genome,"ensGene", test.cats=c("GO:CC", "GO:BP", "GO:MF"))
  #table(goResults$over_represented_pvalue < 0.05)
  tmp = filter(goResults,over_represented_pvalue < 0.05) %>% mutate(hitsPerc = numDEInCat*100/numInCat)
  p1=plotgoseq(tmp,"BP")
  #p2=plotgoseq(tmp,"CC")
  #p3=plotgoseq(tmp,"MF")
  #print(p1)
  #print(grid.arrange(grobs = list(p1,p2,p3),nrow=3,top=TITLE))
  #rm(pwf,goResults,tmp,p1,p2,p3,data,TITLE)
}
plottingGOSEQ = function(data,TITLE,thresholdFC,genome) {
  goseqDat <- dplyr::filter(data, gene_biotype == "protein_coding")
  goseqDat = unique(goseqDat[,c(1:9)]) #14831
  #rownames(goseqDat) = goseqDat$entrezgene_id
  
  gseaDatUP = filter(goseqDat,log2FoldChange > thresholdFC)
  gseaDatDown = filter(goseqDat,log2FoldChange < thresholdFC)
  
  genes = as.integer(goseqDat$padj < 0.05 & !is.na(goseqDat$padj))
  
  names(genes)=goseqDat$ensembl_gene_id
  table(genes)
  
  genesUP = as.integer(gseaDatUP$padj < 0.05 & !is.na(gseaDatUP$padj))
  names(genesUP)=gseaDatUP$ensembl_gene_id
  table(genesUP)
  
  genesDOWN = as.integer(gseaDatDown$padj < 0.05 & !is.na(gseaDatDown$padj))
  names(genesDOWN)=gseaDatDown$ensembl_gene_id
  table(genesDOWN)
  
  print("ALL GENES")
  p1=calculateGOSEQ(genes,paste0(TITLE,"\nGOTERM for all genes (pval < 0.05)"),genome)
  print("UP GENES")
  p2=calculateGOSEQ(genesUP,paste0(TITLE,"\nGOTERM for UP regulated genes (pval < 0.05 & |logFC| > ",thresholdFC,")"),genome)
  print("DOWN GENES")
  p3=calculateGOSEQ(genesDOWN,paste0(TITLE,"\nGOTERM for DOWN regulated genes (pval < 0.05 & |logFC| > ",thresholdFC,")"),genome)
  
  ggarrange(p1,p2,p3, labels = c("A", "B","C"), common.legend = TRUE, legend = "bottom",nrow=3)
  
}
filterdata2pathwayfoundplotting = function(data,thresholdpadj,thresholdFC,genome) {
  
  dataf = filter(data,padj < thresholdpadj)
  up = filter(dataf,log2FoldChange > thresholdFC)
  down = filter(dataf,log2FoldChange < -thresholdFC)
  
  p1=pathwayfound(dataf,genome)
  p2=pathwayfound(up,genome)
  p3=pathwayfound(down,genome)
  
  tmp = list(p1,p2,p3)
  names(tmp) = c("ALL","UP","DOWN")
  return(tmp)
}
plotFGSEA = function(data,index) {
  return(plotGseaTable(pathways = data[[index]]$plotGseaTable$pathways,
                       stats = data[[index]]$plotGseaTable$stats,
                       fgseaRes = data[[index]]$plotGseaTable$fgseaRes,
                       gseaParam = data[[index]]$plotGseaTable$gseaParam,
                       render = data[[index]]$plotGseaTable$render))
}
