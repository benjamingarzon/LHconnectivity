
rm(list=ls())

source('~/Software/LeftHand/VBM/get_session_data.R')
source('./fc_analysis_func.R')

ROI.lh = read.table('/home//share/LeftHand/LHconnectivity/Stockholm/zstat_lh_all_ROI.csv', header=F, sep = '')
ROI.rh = read.table('/home//share/LeftHand/LHconnectivity/Stockholm/zstat_rh_all_ROI.csv', header=F, sep = '')
ROI.lh_rh = read.table('/home//share/LeftHand/LHconnectivity/Stockholm/zstat_lh-rh_all_ROI.csv', header=F, sep = '')
ROI.rh_lh = read.table('/home//share/LeftHand/LHconnectivity/Stockholm/zstat_rh-lh_all_ROI.csv', header=F, sep = '')
colnames(ROI.lh) = colnames(ROI.rh) = colnames(ROI.lh_rh) = colnames(ROI.rh_lh) = ROInames

data.activation = read.table('/home//share/LeftHand/LHconnectivity/Stockholm/zstat_lh_all.txt', header=T, sep = ',')
data.activation$SUBJECT = sapply(as.character(data.activation$SUBJECT), function(x) strsplit(x, '-')[[1]][2])

data.activation = merge(data.activation, data, by = c("SUBJECT", "SESSION"))

fc_cols = 1:ncol(ROI.lh)
results.lh = lapply(seq(length(fc_cols)), function(x) dotest_Stockholm(x, ROI.lh, data.activation, ylim = c(-10, 10) ))
results.rh = lapply(seq(length(fc_cols)), function(x) dotest_Stockholm(x, ROI.rh, data.activation, ylim = c(-10, 10) ))
results.lh_rh = lapply(seq(length(fc_cols)), function(x) dotest_Stockholm(x, ROI.lh_rh, data.activation, ylim = c(-10, 10) ))
results.rh_lh = lapply(seq(length(fc_cols)), function(x) dotest_Stockholm(x, ROI.rh_lh, data.activation, ylim = c(-10, 10) ))

p.values = sapply(results, function(x) x$p)
coefs = sapply(results, function(x) x$coef)

hist(p.values, 20, xlab = "p value")

p.values.fdr = p.adjust(p.values, method = "fdr")

