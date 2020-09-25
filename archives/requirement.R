


library(gridExtra)
library(ggplot2)
library(dplyr)
library(Seurat)
library(RColorBrewer)
library(Rsamtools)
library(stringr)
library(UpSetR)
library(ggforce)
library(sctransform)
library("viridis")

path="E:/TRAVAIL/William/Lyon/DATA/"
prloc=paste0(path,"INMG_SingleCell-master/")
source(file=paste0(prloc,"scripts/functions_stock.R"),local=T)
paths_output=c(paste0(prloc,"results/DellOrsoD0/"),
               paste0(prloc,"results/DeMicheliD0/"),
               paste0(prloc,"results/GiordaniD0/"))
##################### VARIABLE
project_name=c("DellOrso","DeMicheli","Giordani")
