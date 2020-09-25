library(limma)
library(edgeR)
library(RColorBrewer)
library(dplyr)
library(topGO)
library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)

currentDirectory = "/Users/williamjarassier/Downloads/NGS19-109_Resultats/"
a = c(paste0(currentDirectory,"ProjetTGF/Analyses_supervisees/NGS19-109_ensemblGRCm38r96_ProjetTGF__deseq2_results_contrast_TAM_vs_OIL.tsv"),
          paste0(currentDirectory,"ProjetTGF_2v2/Analyses_supervisees/NGS19-109_ensemblGRCm38r96_ProjetTGF_2v2__deseq2_results_contrast_TAM_vs_OIL_MAplot.tsv"),
          paste0(currentDirectory,"ProjetSET/Analyses_supervisees/NGS19-109_ensemblGRCm38r96_ProjetSET__deseq2_results_contrast_TAM_vs_OIL.tsv"))


b = read.delim(a[1],header=T,sep="\t",dec=".")
c = read.delim(a[2],header=T,sep="\t",dec=".")
d = read.delim(a[3],header=T,sep="\t",dec=".")

library(DESeq2)
b=DESeq(b)
res = results(b,contrast=c("TAM","OIL"))

geneList = b[,3]
names(geneList) = as.character(b[,1])
geneList = sort(geneList,decreasing=T)
gene <- names(geneList)[abs(geneList) > 2]
gene2 <- bitr(gene, fromType = "ENSEMBL",
           toType = "ENTREZID",
           OrgDb = org.Mm.eg.db)
egoBP = enrichGO(gene=as.vector(gene2$ENTREZID),OrgDb = org.Mm.eg.db,
               ont="BP",minGSSize = 1, pvalueCutoff = 0.05, qvalueCutoff = 1,
               readable = T)
egoMF = enrichGO(gene=as.vector(gene2$ENTREZID),OrgDb = org.Mm.eg.db,
                 ont="MF",minGSSize = 1, pvalueCutoff = 0.05, qvalueCutoff = 1,
                 readable = T)
egoCC = enrichGO(gene=as.vector(gene2$ENTREZID),OrgDb = org.Mm.eg.db,
                 ont="CC",minGSSize = 1, pvalueCutoff = 0.05, qvalueCutoff = 1,
                 readable = T)
barplot(egoBP, showCategory=10,title="GO terms enriched in group_4")
barplot(egoMF, showCategory=10,title="GO terms enriched in group_4")
barplot(egoCC, showCategory=10,title="GO terms enriched in group_4")



row.names(analyses_super_data)=analyses_super_data$Row.names


ggo <- groupGO(gene     = ,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)


#ggplot(analyses_super_data) +
#  geom_point(aes(x=log2FoldChange,y=-log(pvalue,10)))
#table(analyses_super_data$padj < 0.05)
#subset=analyses_super_data[analyses_super_data$padj < 0.05,]
#subset=analyses_super_data[order(subset$pvalue),]
#over = subset[subset$stat>=0,]
#under = subset[subset$stat<0,]
install.packages("vroom")






cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Mouse_cell_markers.txt') %>%
  tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
  dplyr::select(cellMarker, geneID) %>%
  dplyr::mutate(geneID = strsplit(geneID, ', '))



tmp<- bitr(analyses_super_data$ensembl_gene_id, fromType = "ENSEMBL",
           toType = c("ENTREZID", "SYMBOL"),
           OrgDb = org.Mm.eg.db)
colnames(tmp)=c("ensembl_gene_id","entrezid","mgi_symbol")
analyses_super_data2 = inner_join(analyses_super_data,tmp,by=c("ensembl_gene_id","mgi_symbol"))
gene = analyses_super_data2$entrezid
res =  enricher(as.character(analyses_super_data2$entrezid), TERM2GENE=cell_markers, minGSSize=1)
x = filter(res@result,p.adjust<.05,qvalue<0.2)
y = mutate(x, geneRatio = parse_ratio(GeneRatio)) %>%
  arrange(desc(geneRatio))
z <- mutate(y, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
library(ggplot2)
library(forcats)
library(enrichplot)

ggplot(z, showCategory = 20, 
       aes(richFactor, fct_reorder(Description, richFactor))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(2, 10)) +
  theme_minimal() + 
  xlab("rich factor") +
  ylab(NULL) + 
  ggtitle("Enriched Ontology")

w=mutate(z, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))

library(ReactomePA)
x <- gsePathway(analyses_super_data$ensembl_gene_id,organism="mouse")


