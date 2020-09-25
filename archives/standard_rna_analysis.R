args = commandArgs(trailingOnly=TRUE)

currentDirectory = args[1]
datafiles = args[2]
output= args[3]

source(file="/Users/williamjarassier/TRAVAIL/scripts/functions_stock.R",local=T)