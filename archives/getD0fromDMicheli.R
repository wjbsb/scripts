## Get DeMichelli D0 data
# From GEO .txt matrix

library(dplyr)
#library(tibble)
library(Seurat)

prloc="~/INMG_SingleCell/"
setwd(prloc)
origdir="data/DeMicheli/"
targetdir="data/DeMicheliD0/"
dir.create(targetdir,recursive=T)

#°# ~~~~ Data from DeMicheli ~~~~
# first attempt yield error cause no colname for 1st column, so in terminal: cd data/DeMicheli
#sed -e 1's/.*/genenames &/' GSE143437_DeMicheli_MuSCatlas_rawdata.txt > DeMicheli_rawdata.txt

D0to7 <- read.table(paste0(origdir,"DeMicheli_rawdata.txt"), 
                    sep="\t", header=TRUE, row.names=1) # correctly take rownames

# extract data in homeostasis (day zero):
dmi0 <- D0to7 %>% select(starts_with("D0")) 
svgdmi0 <- cbind(rownames(dmi0), dmi0)
colnames(svgdmi0)[1] <- 'genenames'
write.table( svgdmi0, paste0(targetdir,"rawdataD0.txt"), sep="\t", 
             quote=F, row.names = F)
ho <- data.frame(be=c(1,2))
write.table(ho, paste0(prloc,targetdir,"rawdataD0.txt"))
