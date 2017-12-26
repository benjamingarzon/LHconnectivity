
# definitions

selected.Stockholm = c(1, 1, 1, 1, 1, 1, 1,
             1, 1, 1, 1,
             1, 1, 1, 1)

selected.Berlin = c(0, 0, 1, 1, 1, 1, 1,
                       1, 1, 1, 1,
                       1, 1, 1, 1)

selected.Berlin = c(1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1,
                    1, 1, 1, 1)


#selected = c(1, 1, 1, 1, 1, 1, 1,
#             0, 0, 0, 0,
#             1, 1, 1, 1)

Berlin_experimental = c('1001', '1003', '1004', '1005', '1006', '1007', '1008', '1009', 
                   '1010', '1011', '1012', '1013', '1014', '1015', '2008')

Berlin_control = c('1002', '2001', '2002', '2003', '2004', '2006', '2007', '2009', '2011', 
                   '2013', '2014', '2015', '2016', '2017', '2019', '2020')



ROInames = unlist(read.table('/home/ALDRECENTRUM/benjamin.garzon/Software/LeftHand/LHconnectivity/LHconnectivity/data/atlas_lefthand8.txt'))

# 
# ROInames2 = c(
#   'CerebVIIIa_Left',
#   'CerebVIIIa_Right',
#   'CerebV_Left',
#   'CerebV_Right',
#   'GMPrimarymotorcortexBA4a_Left',
#   'GMPrimarymotorcortexBA4a_Right',
#   'GMSMA',
#   'GMVisualcortexV3V_Left',
#   'GMVisualcortexV3V_Right',
#   'GMVisualcortexV5_Left',
#   'GMVisualcortexV5_Right',
#   'Putamen_Left',
#   'Putamen_Right',
#   'Thalamus_Left',
#   'Thalamus_Right')[selected == 1]


ROIcolors = c('red', 'red', 'red', 'red', 'green', 'green', 'gray',               
              'blue', 'blue', 'cyan', 'cyan', 'yellow', 'yellow', 'orange', 'orange')




#-functions ---------------------------------------------

se = function(x){return(sd(x)/sqrt(length(x)))}

detrend = function(X, id.vars){
  # use lmerTest
  fc_cols = grep( "fc_", colnames(X))
  X.melt = melt(X, id.vars = id.vars)
  model = lm(value ~ 1 + REALWEEK, data = X.melt)
#  y = rowMeans(X[fc_cols])   
#  model = lm(y ~ 1 + X$REALWEEK)
  a = coef(model)[1]
  b = coef(model)[2]
  for (i in fc_cols) {
     X[, i] = scale(X[, i] - a - b*X$REALWEEK, center = T, scale = F)
   }
    
  print(a)
  print(b)
  print(summary(model))
  # detrend by the average
#  for (i in fc_cols) {
#    y = unlist(X[, i])
#    model = lm(y ~ 1 + X$SESSION)
#    X[, i] = scale(model$residuals, center = T, scale = T)
#  }
  return(X)
}

find_outliers = function(myfc, fc_cols){
  
  for (sess in unique(myfc$SESSION)){
    for (fc_col in fc_cols){
      y = myfc[myfc$SESSION == sess, fc_col]
      m = mean(y)
      s = sd(y)
      outlier = abs(y - m) > 3*s
      y[outlier] = NA
      myfc[myfc$SESSION == sess, fc_col] = y
    } 
  }  
  
  print(paste("Total outliers: ", sum(is.na(myfc[, fc_cols] )), "of", prod(dim(myfc[, fc_cols]))))
  return(myfc)
}


dotest_Stockholm = function(index, X, myfc, ylim = c(-.4, .4), name = NULL, plotit = F){
  
  y = unlist(X[index])
  myfc$y = y
  
  #  myfc.12 = subset(myfc, WEEK.MEAN %in% c(1, 2) & PHASE == 1)  
  #  model = lmer(y ~ 1 + WEEK.MEAN*GROUP + FD + (1|SUBJECT), data = myfc.12) #/WEEK.MEAN
  #p = summary(model)$coefficients["WEEK.MEAN:GROUPearly", "Pr(>|t|)"]
  
# ACHTUNG!!
  myfc.1 = subset(myfc, PHASE == 1)
  myfc = myfc.1

  model = lmer(y ~ 1 + TRAINING*GROUP + FD + (1|SUBJECT), data = myfc.1) #/WEEK.MEAN + + REALWEEK
  p = summary(model)$coefficients["TRAINING:GROUPearly", "Pr(>|t|)"]
  
  FFX = fixef(model)
  RFX = ranef(model)

  myfc$y.dem = 0
  
  for (subject in unique(myfc$SUBJECT)){
    myfc$y.dem[myfc$SUBJECT == subject] = myfc$y[myfc$SUBJECT == subject] - RFX$SUBJECT[subject, ]
    #      mean(myfc$y[myfc$SUBJECT == subject], na.rm = T)
  }
  myfc$y.dem = myfc$y.dem - FFX["FD"]*myfc$FD# - FFX["(Intercept)"] #- FFX["GROUPearly"]*(myfc$GROUP == 'early') #-
  #  FFX["REALWEEK"]*myfc$REALWEEK
  
  if (p < 0.05 & plotit)
  {
    if (!is.null(name)) {
      print(name)  
    } else {    
      print(paste('Connection', index))
    }
    
    print("-----------SIGNIFICANT ------------------")
    
    
    print(summary(model))#$coefficients)
    
    fc.mean = myfc %>% group_by(SESSION.MEAN, GROUP) %>% summarise(
      fc_mean = mean(y.dem, na.rm=T),
      fc_se = sqrt(var(y.dem, na.rm=T)/sum(!is.na(y.dem)))
    )
    fc.run = myfc %>% group_by(SUBJECT, TRAINING, PHASE, GROUP) %>% summarise(
      fc_mean = mean(y.dem, na.rm=T)
    )    
    
    plot.p1 = ggplot(data = fc.mean, aes(x = SESSION.MEAN, y = fc_mean, col = GROUP)) + geom_line(lwd = 2) + 
      geom_point() +  geom_errorbar(aes(ymin=fc_mean - fc_se, ymax=fc_mean + fc_se), width=.1) + ylim(ylim)
    
    #plot.p2 = ggplot(data = fc.run, aes(x = TRAINING, y = fc_mean, col = SUBJECT, group = SUBJECT)) + geom_line() + 
    #  geom_point() + facet_grid(. ~  PHASE)
    
    print(plot.p1)
    #print(plot.p2)
  }
  
  return(list(p = p, coef = FFX["TRAINING:GROUPearly"], y.dem = myfc$y.dem))
  

} 


doplot_Stockholm = function(index, X, myfc, ylim = c(-.4, .4), name = NULL, plotit = T){

  y = unlist(X[, index])
  

  if (plotit[index])
  {
    if (!is.null(name)) {
      print(name)  
    } else {    
      print(paste('Connection', index))
    }
    
    # ACHTUNG!!!
    myfc = subset(myfc, PHASE == 1)
    myfc$y.dem = y
    
    fc.mean = myfc %>% group_by(SESSION.MEAN, GROUP, CONTRAST) %>% summarise(
      fc_mean = mean(y.dem, na.rm=T),
      fc_se = sqrt(var(y.dem, na.rm=T)/sum(!is.na(y.dem)))
    )
    fc.run = myfc %>% group_by(SUBJECT, TRAINING, PHASE, GROUP, CONTRAST) %>% summarise(
      fc_mean = mean(y.dem, na.rm=T)
    )    
    
    plot.p0 = ggplot(data = myfc, aes(x = SUBJECT, y = y.dem, col = GROUP)) + 
      geom_point() +  xlab('WEEK') + ylab('a. u.') +
      facet_grid(CONTRAST ~ SESSION.MEAN) + theme(legend.position="none")
    
    
    print(plot.p0)
    plot.p1 = ggplot(data = fc.mean, aes(x = SESSION.MEAN, y = fc_mean, col = GROUP)) + geom_line(lwd = 2) + 
      geom_point() +  geom_errorbar(aes(ymin=fc_mean - fc_se, ymax=fc_mean + fc_se), width=.1) + xlab('WEEK') + ylab('a. u.') +
      facet_grid(. ~ CONTRAST) + theme(legend.position="none") #+ ylim(ylim)
    
    #plot.p2 = ggplot(data = fc.run, aes(x = TRAINING, y = fc_mean, col = SUBJECT, group = SUBJECT)) + geom_line() + 
    #  geom_point() + facet_grid(. ~  PHASE)
    
    print(plot.p1)
    #print(plot.p2)
  } 
 return()
}

dotest_Berlin = function(index, X, myfc, ylim = c(-.4, .4)){
  y = unlist(X[index])
  myfc$y = y
  
  myfc$GROUP = 'none' 
  myfc$GROUP[myfc$SUBJECT %in% Berlin_experimental] = 'experimental'
  myfc$GROUP[myfc$SUBJECT %in% Berlin_control] = 'control'
  
#  myfc.1 = subset(myfc, SESSION < 15 & GROUP != 'none')# & RUN == 1 
#  myfc.1 = subset(myfc, GROUP != 'none' & RUN == 1 ) #& SESSION < 11
  myfc.1 = subset(myfc, GROUP != 'none' & SESSION < 11)
  model = lmer(y ~ 1 + SESSION*GROUP + FD  + (1|SUBJECT/RUN), data = myfc.1) 

#  myfc.1 = subset(myfc, GROUP != 'none' & RUN == 1 )
#  model = lmer(y ~ 1 + SESSION*GROUP + FD  + (1|SUBJECT), data = myfc.1)

  myfc = myfc.1
  p = summary(model)$coefficients["SESSION:GROUPexperimental", "Pr(>|t|)"]
  
  FFX = fixef(model)
  RFX = ranef(model)

#myfc$y.dem = 0

for (subject in unique(myfc$SUBJECT)){
  myfc$y.dem[myfc$SUBJECT == subject] = myfc$y[myfc$SUBJECT == subject] - RFX$SUBJECT[subject, ]
}
#myfc$y.dem = myfc$y.dem - FFX["FD"]*myfc$FD #- FFX["GROUPexperimental"]*(myfc$GROUP == 'experimental') - FFX["(Intercept)"] #-

myfc$y.dem = myfc$y - FFX["FD"]*myfc$FD



  
  if (p < 0.05)
  {  
    print(paste('Connection', index))

    print("-----------SIGNIFICANT ------------------")
    
    
    print(summary(model))#$coefficients)
    
    
    fc.mean = myfc %>% group_by(SESSION, GROUP) %>% summarise(
      fc_mean = mean(y.dem, na.rm=T),
      fc_se = sqrt(var(y.dem, na.rm=T)/sum(!is.na(y.dem)))
    )
    fc.run = myfc %>% group_by(SUBJECT, SESSION, GROUP) %>% summarise(
      fc_mean = mean(y.dem, na.rm=T)
    )    
    
    plot.p1 = ggplot(data = fc.mean, aes(x = SESSION, y = fc_mean, col = GROUP)) + geom_line(lwd = 2) + 
      geom_point() +  geom_errorbar(aes(ymin=fc_mean - fc_se, ymax=fc_mean + fc_se), width=.1) + ylim(ylim)
    
    #plot.p2 = ggplot(data = fc.run, aes(x = TRAINING, y = fc_mean, col = SUBJECT, group = SUBJECT)) + geom_line() + 
    #  geom_point() + facet_grid(. ~  PHASE)
    
    print(plot.p1)
    #print(plot.p2)
  }
  
  return(list(p = p, coef = FFX["SESSION:GROUPexperimental"], y.dem = myfc$y.dem)) 
  
} 

dotest_Berlin_difference = function(index, X, myfc, ylim = c(-.4, .4)){
  y = unlist(X[index])
  myfc$y = y
  
  myfc$GROUP = 'none' 
  myfc$GROUP[myfc$SUBJECT %in% Berlin_experimental] = 'experimental'
  myfc$GROUP[myfc$SUBJECT %in% Berlin_control] = 'control'
  
  myfc.1 = subset(myfc, GROUP != 'none')# & RUN == 1 
  # nest the session
  model = lmer(y ~ 1 + RUN*GROUP + FD  + (1|SUBJECT), data = myfc.1) 

  p = summary(model)$coefficients["RUN:GROUPexperimental", "Pr(>|t|)"]
  
  FFX = fixef(model)
  RFX = ranef(model)
  
  if (p < 0.05)
  {
    print(paste('Connection', index))
    
    print("-----------SIGNIFICANT ------------------")
    
    myfc$y.dem = 0
    
    print(summary(model))#$coefficients)
    for (subject in unique(myfc$SUBJECT)){
      myfc$y.dem[myfc$SUBJECT == subject] = myfc$y[myfc$SUBJECT == subject] - RFX$SUBJECT[subject, ]
    }
    myfc$y.dem = myfc$y.dem - FFX["FD"]*myfc$FD - FFX["GROUPexperimental"]*(myfc$GROUP == 'experimental') - FFX["(Intercept)"] #-
    
    #fc.mean = myfc %>% group_by(SESSION, GROUP) %>% summarise(
    #  fc_mean = mean(y.dem, na.rm=T),
    #  fc_se = sqrt(var(y.dem, na.rm=T)/sum(!is.na(y.dem)))
    #)
    #fc.run = myfc %>% group_by(SUBJECT, SESSION, GROUP) %>% summarise(
    #  fc_mean = mean(y.dem, na.rm=T)
    #)    
    
    #plot.p1 = ggplot(data = fc.mean, aes(x = SESSION, y = fc_mean, col = GROUP)) + geom_line(lwd = 2) + 
    #  geom_point() +  geom_errorbar(aes(ymin=fc_mean - fc_se, ymax=fc_mean + fc_se), width=.1) + ylim(-0.4, 0.4)
    
    #plot.p2 = ggplot(data = fc.run, aes(x = TRAINING, y = fc_mean, col = SUBJECT, group = SUBJECT)) + geom_line() + 
    #  geom_point() + facet_grid(. ~  PHASE)
    
    #print(plot.p1)
    #print(plot.p2)
  }
  
  return(list(p = p, coef = FFX["RUN:GROUPexperimental"]))
  
} 


select_dotest =function(dataset) {
  
  if (dataset == 'Stockholm')
    return(dotest_Stockholm)
  
  if (dataset == 'Berlin')
    return(dotest_Berlin)
  
  if (dataset == 'Berlin_difference')
    return(dotest_Berlin_difference)  
  
}

dotest.2 = function(y, myfc){
  myfc$y = y
  
  #  myfc.12 = subset(myfc, WEEK.MEAN %in% c(1, 2) & PHASE == 1)  
  #  model = lmer(y ~ 1 + WEEK.MEAN*GROUP + FD + (1|SUBJECT), data = myfc.12) #/WEEK.MEAN
  #p = summary(model)$coefficients["WEEK.MEAN:GROUPearly", "Pr(>|t|)"]
  
  #  myfc.12 = subset(myfc, SESSION.MEAN < 5)  
  myfc = subset(myfc,  PHASE == 1)
  model = lmer(y ~ 1 + t.TRAINING.1*GROUP + t2.TRAINING.1*GROUP + FD + (1|SUBJECT), data = myfc) #/WEEK.MEANREALWEEK + 
  p = summary(model)$coefficients["t.TRAINING.1:GROUPearly", "Pr(>|t|)"]
  p2 = summary(model)$coefficients["GROUPearly:t2.TRAINING.1", "Pr(>|t|)"]  
  
  FFX = fixef(model)
  RFX = ranef(model)
  
  
  if (p < 0.05)
  {
    print("-----------SIGNIFICANT ------------------")
    
    myfc$y.dem = 0
    
    print(summary(model))#$coefficients)
    for (subject in unique(myfc$SUBJECT)){
      myfc$y.dem[myfc$SUBJECT == subject] = myfc$y[myfc$SUBJECT == subject] - RFX$SUBJECT[subject, ]
      #      mean(myfc$y[myfc$SUBJECT == subject], na.rm = T)
    }
    myfc$y.dem = myfc$y.dem - FFX["FD"]*myfc$FD - FFX["GROUPearly"]*(myfc$GROUP == 'early') - FFX["(Intercept)"] #-
    #    FFX["REALWEEK"]*myfc$REALWEEK
    
    fc.mean = myfc %>% group_by(SESSION.MEAN, GROUP) %>% summarise(
      fc_mean = mean(y.dem, na.rm=T),
      fc_se = sqrt(var(y.dem, na.rm=T)/sum(!is.na(y.dem)))
    )
    fc.run = myfc %>% group_by(SUBJECT, TRAINING, PHASE, GROUP) %>% summarise(
      fc_mean = mean(y.dem, na.rm=T)
    )    
    
    plot.p1 = ggplot(data = fc.mean, aes(x = SESSION.MEAN, y = fc_mean, col = GROUP)) + geom_line(lwd = 2) + 
      geom_point() +  geom_errorbar(aes(ymin=fc_mean - fc_se, ymax=fc_mean + fc_se), width=.1) + ylim(-0.2, 0.2)
    
    #plot.p2 = ggplot(data = fc.run, aes(x = TRAINING, y = fc_mean, col = SUBJECT, group = SUBJECT)) + geom_line() + 
    #  geom_point() + facet_grid(. ~  PHASE)
    
    print(plot.p1)
    #print(plot.p2)
  }
  
  return(rbind(p, p2))
  
} 



dotest_global = function(y, myfc){
  myfc$y = y
  var.names = c("Intercept", 
                "SESSION.MEAN.4","SESSION.MEAN.8","SESSION.MEAN.13","SESSION.MEAN.17","SESSION.MEAN.21",
                "GROUPearly", "FD", 
                "SESSION.MEAN.4_x_GROUPearly","SESSION.MEAN.8_x_GROUPearly","SESSION.MEAN.13_x_GROUPearly",
                "SESSION.MEAN.17_x_GROUPearly","SESSION.MEAN.21_x_GROUPearly")
  
  contrast.names = c("GROUPearly", 
                     "SESSION.MEAN.4_x_GROUPearly", "SESSION.MEAN.8_x_GROUPearly","SESSION.MEAN.13_x_GROUPearly",
                     "SESSION.MEAN.17_x_GROUPearly","SESSION.MEAN.21_x_GROUPearly", "SESSION.MEAN.4.PEAK", 
                     "SESSION.MEAN.8.PEAK")
  
  c.1       = c(0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0)
  c.2       = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0)
  c.3       = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0)
  c.4       = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0)
  c.5       = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)
  c.6       = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)
  c.7       = c(0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0)
  c.8       = c(0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0)
  
  #cont.mat = rbind(c.1, c.2, c.3, c.4, c.5, c.6)
  cont.mat = rbind(c.2, c.3)  
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names[c(2, 3)]
  
  cont.mat2 = rbind(c.7, c.8)  
  colnames(cont.mat2) = var.names
  rownames(cont.mat2) = contrast.names[c(7, 8)]
  
  model = lmer(y ~ 1 + as.factor(SESSION.MEAN)*as.factor(GROUP) + FD + (1 |SUBJECT), data = myfc, REML = T)
  FFX = fixef(model)
  RFX = ranef(model)
  
  glh = glht(model, linfct = cont.mat, alternative="greater")
  contrast.pvalues = summary(glh, test = Ftest())$test$pvalue[[1]]
  
  glh2 = glht(model, linfct = cont.mat2, alternative="greater")
  contrast.pvalues2 = summary(glh2)$test$pvalues  
  
  #contrast.coefs = coef(glh)
  if (contrast.pvalues < 0.05)
  {
    print(contrast.pvalues)
    print("-----------SIGNIFICANT ------------------")
    myfc$y.dem = 0
    
    print(summary(model))#$coefficients)
    for (subject in unique(myfc$SUBJECT)){
      myfc$y.dem[myfc$SUBJECT == subject] = myfc$y[myfc$SUBJECT == subject] - RFX$SUBJECT[subject, ]
    }
    myfc$y.dem = myfc$y.dem - FFX["FD"]*myfc$FD - FFX["as.factor(GROUP)early"]*(myfc$GROUP == 'early') - FFX["(Intercept)"] #-
    #      FFX["as.factor(SESSION.MEAN)4"]*(myfc$SESSION.MEAN == '4') -
    #      FFX["as.factor(SESSION.MEAN)8"]*(myfc$SESSION.MEAN == '8') -
    #      FFX["as.factor(SESSION.MEAN)13"]*(myfc$SESSION.MEAN == '13') -
    #      FFX["as.factor(SESSION.MEAN)17"]*(myfc$SESSION.MEAN == '17') -
    #      FFX["as.factor(SESSION.MEAN)21"]*(myfc$SESSION.MEAN == '21')
    
    
    fc.mean = myfc %>% group_by(SESSION.MEAN, GROUP) %>% summarise(
      fc_mean = mean(y.dem, na.rm=T),
      fc_se = sqrt(var(y.dem, na.rm=T)/sum(!is.na(y.dem)))
    )
    fc.run = myfc %>% group_by(SUBJECT, TRAINING, PHASE, GROUP) %>% summarise(
      fc_mean = mean(y.dem, na.rm=T)
    )    
    
    plot.p1 = ggplot(data = fc.mean, aes(x = SESSION.MEAN, y = fc_mean, col = GROUP)) + geom_line(lwd = 2) + 
      geom_point() +  geom_errorbar(aes(ymin=fc_mean - fc_se, ymax=fc_mean + fc_se), width=.1) #+ ylim(-0.2, 0.2)
    
    
    print(plot.p1)
    
  }   
  
  return(c(contrast.pvalues, contrast.pvalues2))
}

get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  diag(cormat) = 1
  return(cormat)
}

BWRES=1200
CRES=1200
POINTSIZE=35

save_fig = function(figname=NULL, width=6.5, height=6.5, res=600, jpg=F){
  #size in inches
  if (dev.cur()!=1) dev.off()
  if (is.null(figname)) {
    figname = paste0('plot', fig_count)
    
    fig_count <<- fig_count+1
  } 
  print(FIGS_DIR)
  
  print(paste("Generating figure: ", figname))
  if (jpg){
    figname = paste0(FIGS_DIR, figname, '.jpg')
    jpeg(figname, width = floor(res/2.54*width), height = floor(res/2.54*height), pointsize=POINTSIZE)
  } else {
    figname = paste0(FIGS_DIR, figname, '.png')
    png(figname, width = floor(res/2.54*width), height = floor(res/2.54*height), pointsize=POINTSIZE)
  }
  
  
}
