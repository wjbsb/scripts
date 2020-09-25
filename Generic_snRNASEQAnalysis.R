






listCSO = CreateCSO(origintargetdf,thresholdcell,thresholdfeaturemin)
listCSO = sapply(listCSO, function(QC) {tmp = QCseurat(data=QC,
                                                       thresholdfeaturemin,
                                                       thresholdfeaturemax)
return(tmp)})

listAnchors = FindIntegrationAnchors(object.list = listCSO, dims=1:30)
listIntegrated = IntegrateData(anchorset = listAnchors, dims = 1:30)
DefaultAssay(object = listIntegrated) <- "integrated"


listIntegrated <- KNNplusFITSNE(listIntegrated, nbDIM, resolution)

saveRDS(listIntegrated,file="./rds/integrated_seu_fitsne.rds")

# thresholdlogFCmin = 0.25
listIntegratedMarkers <- FindAllMarkers(listIntegrated, 
                                        only.pos = F, 
                                        min.pct = 0.25, 
                                        logfc.threshold = thresholdlogFCmin)
# thresholdlogFC = 2
positive <- listIntegratedMarkers %>% filter(avg_logFC > thresholdlogFC)
negative <- listIntegratedMarkers %>% filter(avg_logFC < -thresholdlogFC)
most.neg <- negative %>% group_by(cluster) %>% top_n(n=10, wt= -avg_logFC)
most.pos <- positive %>% group_by(cluster) %>% top_n(n=10, wt= avg_logFC)

write.table(positive, file = "./ALLMARKERS_Pos_integratedD0.tsv", sep="\t")
write.table(negative, file = "./ALLMARKERS_Neg_integratedD0.tsv", sep="\t")


pdf("./grap_reports_integratedD0.pdf")
DimPlot(listIntegrated, 
        reduction="tsne", 
        label=T,
        cols=definecolors(listIntegrated@active.ident)) + 
  ggtitle("Integrated D0")
DoHeatmap(listIntegrated, features = most.pos$gene, 
          group.colors=definecolors(listIntegrated@active.ident)) 
dev.off()
