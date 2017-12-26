

rm(list=ls())

source('./installpackages.R')
source('./fc_analysis_func.R')

netmats = read.table('/home/share/LeftHand/LHconnectivity/Berlin/melodic.gica/netmats_ridge.csv', header=F, sep = ',')
netmats.table = read.table('/home/share/LeftHand/LHconnectivity/Berlin/ICA_table.csv', header=T, sep = ',')

netmats.table$SUBJECT = sapply(as.character(netmats.table$SUBJECT), function(x) strsplit(x, '-')[[1]][2])

NNETS = ncol(netmats)
NROIS = (1 + sqrt(1 + 8*NNETS))/2
results = lapply(seq(NNETS), function(x) dotest_Berlin(x, netmats, netmats.table, ylim = c(-10, 10))) 


