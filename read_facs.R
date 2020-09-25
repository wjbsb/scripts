

library(flowCore)

a = read.FCS(filename = "/Volumes/LaCie/InstitutNeuroMyoGene_DATA/snRNAseq_FACS/Specimen_001_Unst_001.fcs")

b = data.frame(a@exprs)
colnames(b)

#FSC = evalue la taille des cellules
#SSC = evalue la granularité

library(ggplot2)

hist(b$FSC.A)
hist(b$SSC.A)


ggplot(b,aes(x=log2(FSC.A),y=SSC.A)) + geom_point()
