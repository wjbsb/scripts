



for (i in 1:length(paths_data)) {
  load(paths_data[i])
  pdf(paste0(paths_output[i],"preprocessSeu_V1.pdf"))
  VlnPlot(data1, 
          features =c("nFeature_RNA","nCount_RNA","percent.mt"), 
          ncol = 3,
          pt.size=0.1) 
  FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  DimHeatmap(data, dims=1:15, cells = 500, balanced = TRUE) 
  ElbowPlot(data)
  dev.off()
  
  pdf(paste0(paths_output[i],"preprocessSeu_V2.pdf"))
  vlnplotdata=data@meta.data
  col=colnames(vlnplotdata)
  colour=c("tomato1","tan2","steelblue2")
  # In replacement to VlnPlot
  ifelse(length(unique(as.factor(vlnplotdata[,1])))==1,
         grid.arrange(plottingmetadata1("tomato1","nCount_RNA"),
                      plottingmetadata1("tan2","nFeature_RNA"),
                      plottingmetadata1("steelblue2","percent.mt"),
                      nrow=1,
                      top=as.character(project_name[i])),
         grid.arrange(plottingmetadata1bis("tomato1","nCount_RNA"),
                      plottingmetadata1bis("tan2","nFeature_RNA"),
                      plottingmetadata1bis("steelblue2","percent.mt"),
                      nrow=1,
                      top=as.character(project_name[i]))
         
  )
  # In replacement to FeatureScatter
  plottingmetadata2("chartreuse4","nCount_RNA","nFeature_RNA")
  DimHeatmap(data, dims=1:15, cells = 500, balanced = TRUE) 
  ElbowPlot(data)       
  dev.off()
  rm(data,vlnplotdata,col,colour)
}


plottingmetadata1 = function(color,var) {
  return(ggplot(vlnplotdata,aes(x=vlnplotdata[,1],y=vlnplotdata[,var])) + 
           geom_violin(fill=color) + 
           geom_jitter(shape=16, position=position_jitter(0.2),alpha=1/5) + 
           theme_classic() +
           labs(y="Count",x=var) +
           theme(axis.text.x = element_blank())
  )
}
plottingmetadata1bis = function(color,var) {
  return(ggplot(vlnplotdata,aes(x=vlnplotdata[,1],y=vlnplotdata[,var])) + 
           geom_violin(fill=color) + 
           geom_jitter(shape=16, position=position_jitter(0.2),alpha=1/5) + 
           theme_classic() +
           labs(y="Count",x=var)
  )
}
plottingmetadata2 = function(color,var1,var2) {
  corr=round(cor(x=vlnplotdata[,var1],y=vlnplotdata[,var2],method = "pearson"),2)
  return(
    grid.arrange(
      ggplot(vlnplotdata,aes(x=vlnplotdata[,var1],y=vlnplotdata[,var2])) +
        geom_smooth() +
        geom_jitter(shape=16, position=position_jitter(0.2),alpha=1/5,colour=color) +
        theme_classic() +
        labs(x=var1,y=var2),
      ggplot(vlnplotdata,aes(y=vlnplotdata[,var1],x=vlnplotdata[,var2])) +
        geom_boxplot() +
        theme_classic() +
        coord_flip() +
        labs(x=var2,y=var1),
      nrow=2,
      top=paste0(unique(as.character(vlnplotdata[,1]))," (Pearson Corr = ",corr,")")
    )
  )
}



corr1=round(cor(x=vlnplotdata3[,"nCount_RNA"],y=vlnplotdata3[,"nFeature_RNA"],method = "pearson"),2)
corr2=round(cor(x=vlnplotdata3[,"nCount_RNA"],y=vlnplotdata3[,"nFeature_RNA"],method = "pearson"),2)


abis=grid.arrange(      
  ggplot(vlnplotdata3[as.character(vlnplotdata3$orig.ident)=="i1",],aes(x=nCount_RNA,y=nFeature_RNA)) +
    geom_smooth() +
    geom_jitter(shape=16, position=position_jitter(0.2),alpha=1/5,colour="darkred") +
    theme_classic() +
    ggtitle(paste0(project_name[3]," group : i1\n(Pearson Corr=",corr1,")")) +
    theme(plot.title = element_text( size=8, face="bold.italic")) +
    labs(x="",y=""),
  ggplot(vlnplotdata3[as.character(vlnplotdata3$orig.ident)=="i2",],aes(x=nCount_RNA,y=nFeature_RNA)) +
    geom_smooth() +
    geom_jitter(shape=16, position=position_jitter(0.2),alpha=1/5,colour="darkblue") +
    theme_classic() +
    ggtitle(paste0(project_name[3]," group : i2\n(Pearson Corr=",corr1,")")) +
    theme(plot.title = element_text( size=8, face="bold.italic")) +
    labs(x="",y=""),
  nrow=1)
rm(corr1,corr2)
ater=grid.arrange(abis,
                  ggplot(vlnplotdata3,aes(y=nCount_RNA,x=nFeature_RNA,fill=orig.ident)) +
                    geom_boxplot() +
                    theme_classic() +
                    coord_flip() + 
                    labs(x="",y="",fill="group"),
                  nrow=2)


for (i in 1:length(KNNplusFITSNE_list)) {
  tmp=KNNplusFITSNE_list[[i]]
  nb.clus = max(as.integer(levels(tmp@meta.data$seurat_clusters)))+1
  print(nb.clus)
  mycols <- colorRampPalette(brewer.pal(8, "Set2"))(nb.clus)
  
  data <- FindAllMarkers(tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  save(data,file=paste0(prloc,"data/",project_name[i],"D0/FindAllMarkers.Rdata"))
}