
BiocManager::install("Rgb")
library(Rgb)


path1 =  "/Volumes/LaCie/InstitutNeuroMyoGene_DATA/GSE144871/results/stringtieFPKM/transcripts"
path2 = "/Volumes/LaCie/InstitutNeuroMyoGene_DATA/autre_data_caroline/results_exon/results/stringtieFPKM/transcripts"

create_dataframe_transcript = function(path) {
  files = dir(path = path,
            pattern = ".gtf$",
            full.names = T
            )
  names = unlist(lapply(files, function(x) {
    tmp = stringr::str_replace(string = x , pattern = paste0(path,"/"), replacement = "")
    stringr::str_replace(string = tmp , pattern = ".sorted_transcripts.gtf", replacement = "")
    }))
  
  print("reading gtf")
  tmp_gft = lapply(files, function(x) {
    print(x)
    read.gtf(file = x, attr = "intact", features = NULL)
    })
  print("create dataframe")
  alldata = data.frame()
  for (i in 1:length(tmp_gft)) {
    
    print(names[i])
    tmp = stringr::str_replace_all(tmp_gft[[i]]$attributes,";"," ")
    tmp = stringr::str_split(tmp," ")
    
    gene_id = do.call(rbind,lapply(tmp, function(x) {return(x[2])}))
    transcript_id = do.call(rbind,lapply(tmp, function(x) {return(x[5])}))
    ref_gene_name = do.call(rbind,lapply(tmp, function(x) {
      if (length(x) == 18) {
        return(x[8])
      } else {
        return(x[11])
      }
    }))
    exon_number = do.call(rbind,lapply(tmp, function(x) {  if (length(x) == 18) {
      return(NA)
    } else {
      return(x[8])
    }
    }))
    cov = do.call(rbind,lapply(tmp, function(x) {  if (length(x) == 18) {
      return(x[11])
    } else {
      return(x[14])
    }
    }))
    FPKM = do.call(rbind,lapply(tmp, function(x) {  if (length(x) == 18) {
      return(x[14])
    } else {
      return(NA)
    }
    }))
    TPM = do.call(rbind,lapply(tmp, function(x) {  if (length(x) == 18) {
      return(x[17])
    } else {
      return(NA)
    }
    }))
    
    d = data.frame(gene_id,transcript_id,ref_gene_name,exon_number,cov,FPKM,TPM)
    d$sample = names[i]
    
    alldata = rbind(alldata,d)
  }
  saveRDS(alldata,file = paste0(path,"/globaltranscriptgtf.rds"))
  return(alldata)
}

GSE144871transcript = create_dataframe_transcript(path1)
autre_data_carolinetranscript = create_dataframe_transcript(path2)

reading_design = function(path) {
  tmp = read.table(path)
  rownames(tmp) = stringr::str_replace_all(rownames(tmp),"\\.","-")
  return(tmp)
}
design1 = reading_design("/Volumes/LaCie/InstitutNeuroMyoGene_DATA/GSE144871/design.tsv")
design2 = reading_design("/Volumes/LaCie/InstitutNeuroMyoGene_DATA/autre_data_caroline/design.tsv")


library(dplyr)

subset_data = function(listegene,df,design) {
  
  tmp = list()
  
  tmp[[1]] = df %>% filter(ref_gene_name %in% listegene  & !is.na(exon_number))
  tmp[[2]] = df %>% filter(ref_gene_name %in% listegene  & is.na(exon_number))
  
  design$sample = rownames(design)
  
  tmp[[1]] = inner_join(tmp[[1]],design,by="sample")
  tmp[[2]] = inner_join(tmp[[2]],design,by="sample")
  
  names(tmp) = c("exons","transcripts")
  return(tmp)
}

gli3 = subset_data("Gli3",GSE144871transcript,design1)
dmd = subset_data("Dmd",autre_data_carolinetranscript,design2)


names = c("mdx-ASC-78","mdx-ASC-2","mdx-ASC-3",
          "WT-ASC-1","WT-ASC-2","WT-ASC-3",
          "mdx-QSC-50","mdx-QSC-76","mdx-QSC-77",
          "WT-QSC-782","WT-QSC-A","WT-QSC-B", 
          "mdx1-2DM","mdx2-2DM","mdx3-2DM", 
          "WT1-2DM","WT2-2DM","WT3-2DM",
          "mdx1-myob","mdx2-myob","mdx3-myob",
          "WT1-myob","WT2-myob","WT3-myob")

graph_export = function(df,names) {
  


  print("plot1")
  p1 = df[[2]] %>% 
    mutate(sample = factor(sample, levels=names)) %>%
    ggplot(aes(x=sample,y=as.numeric(TPM),fill=Condition)) + 
    geom_bar(stat = "identity") + 
    coord_flip() +
    labs(y ="TPM", x = "Sample") + 
    theme_bw() 
  print("plot2")
  p2 = df[[2]] %>%
    mutate(sample = factor(sample, levels=names)) %>%
    ggplot(aes(x=sample,y=as.numeric(TPM),fill=Condition)) + 
    geom_bar(stat = "identity") + 
    facet_wrap(~transcript_id) + 
    coord_flip() +
    labs(y ="TPM", x = "Sample") + 
    theme_bw()
  print("plot3")
  p3 = df[[2]] %>%
    ggplot(aes(x=transcript_id,y=as.numeric(TPM))) + 
    geom_boxplot() + 
    geom_point(data=df[[2]],aes(x=transcript_id,y=as.numeric(TPM),colour=Condition)) + 
    coord_flip() +
    labs(y ="TPM", x = "Sample") + 
    theme_bw()
  print("plot4")
  p4 = df[[1]] %>%
    ggplot(aes(x=as.numeric(exon_number),y=as.numeric(cov),fill=transcript_id)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    facet_wrap(~transcript_id) + 
    theme_bw() +
    labs(x = "exon number", y = "coverage")
  print("plot5")
  p5 = df[[1]] %>%
    ggplot(aes(x=as.numeric(exon_number),y=as.numeric(cov),fill=Condition)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_bw() +
    labs(x = "exon number", y = "coverage")
  print("plot6")
  p6 = df[[1]] %>%
    ggplot(aes(x=as.numeric(exon_number),y=as.numeric(cov),fill=Condition)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    facet_wrap(~transcript_id) + 
    theme_bw() +
    labs(x = "exon number", y = "coverage")
  
  return(list(p1,p2,p3,p4,p5,p6))
}
 
library(ggplot2)
pdf(file="/Volumes/LaCie/InstitutNeuroMyoGene_DATA/autre_data_caroline/TAR.pdf",width = 20,height = 15)
graph_export(dmd,names)
dev.off()

pdf(file="/Volumes/LaCie/InstitutNeuroMyoGene_DATA/GSE144871/TAR.pdf",width = 20,height = 15)
graph_export(gli3)
dev.off()

library(tidyverse)
library(rstatix)
library(ggpubr)
dmd[[2]]$TPM = as.numeric(dmd[[2]]$TPM)


p1=dmd[[2]] %>%
  ggplot(aes(x=transcript_id,y=log2(as.numeric(TPM)))) + 
  geom_boxplot() + 
  geom_point(data=dmd[[2]],aes(x=transcript_id,y=log2(as.numeric(TPM)),colour=Condition)) + 
  coord_flip() +
  labs(y ="log2 TPM", x = "Sample") + 
  theme_bw()
p2=dmd[[2]] %>%
  ggplot(aes(x=sample,y=log2(as.numeric(TPM)))) + 
  geom_boxplot() + 
  geom_point(data=dmd[[2]],aes(x=sample,y=log2(as.numeric(TPM)),colour=transcript_id)) + 
  coord_flip() +
  labs(y ="log2 TPM", x = "Sample") + 
  theme_bw()
library(gridExtra)
grid.arrange(p1,p2,ncol = 2)

tmp = dmd[[2]]
tmp2 = tmp[tmp$Condition == "QSC",]

a = ggpaired(dmd[[2]], x = "transcript_id", y = "TPM", fill="Condition", ylab = "TPM", xlab = "transcripts",orientation="horizontal")
b = ggpaired(tmp2, x = "transcript_id", y = "TPM", fill = "grp",ylab = "TPM", xlab = "transcripts",orientation="horizontal")
grid.arrange(a,b,ncol = 2)



############################################
design2$sample = rownames(design2)
data2 = inner_join(autre_data_carolinetranscript,design2,by="sample")

MYOB = filter(data2,!is.na(FPKM) & stringr::str_detect(string = data2$grp,pattern = "MYOB"))
DM = filter(data2,!is.na(FPKM) & stringr::str_detect(string = data2$grp,pattern = "DM"))
QSC = filter(data2,!is.na(FPKM) & stringr::str_detect(string = data2$grp,pattern = "QSC"))
ASC = filter(data2,!is.na(FPKM) & stringr::str_detect(string = data2$grp,pattern = "ASC"))

mdx = filter(MYOB,stringr::str_detect(string=MYOB$grp,pattern="MDX"))
wt = filter(MYOB,stringr::str_detect(string=MYOB$grp,pattern="WT"))

dim(mdx)[1] - dim(wt)[1]

uniquetranscript = unique(mdx$transcript_id)
wt2 = filter(wt,transcript_id %in% mdx)

wilcox.test(x=as.numeric(mdx$FPKM), y=as.numeric(wt$FPKM), paired = F)

wilcox_results = dmd[[2]] %>% rstatix::wilcox_test(TPM ~ transcript_id) %>% add_significance()
install.packages("coin")
library(coin)
effectsize = dmd[[2]] %>% wilcox_effsize(TPM+transcript_id ~ transcript_id)
