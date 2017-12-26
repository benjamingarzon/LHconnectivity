rm(list = ls())

# to do: impute data, if not make it hierarchical but only fitting wherever data are valid


setwd("~/Software/LeftHand/LHconnectivity/scripts/stats")
source('./installpackages.R')
source('~/Software/LeftHand/VBM/get_session_data.R')
source('./fc_analysis_func.R')
source('./init.R')


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
##################################
# load behavioural data
##################################

load('~/Software/LeftHand//Behaviour/analysis/speed_accuracy/speed_accuracy.RData')


##################################
# load connectivity data
##################################

fc_evol = read.table('/home//share/LeftHand/LHconnectivity/Stockholm/lw_fc_rest_nlp.csv', header=T, sep = ',')
fc_evol$SUBJECT = sapply(as.character(fc_evol$SUBJECT), function(x) strsplit(x, '-')[[1]][2])


##################################
# fit exponentials to behavioural data
##################################

normData <- function (dv, sess, nSess = 5) {
   dv / mean(dv[sess<=nSess], na.rm = T)
}

myvar = 'dev'
dAll.orig = subset(training.daily, (GROUP == 'early' & PHASE == 1) | (GROUP == 'late' & PHASE == 2))


mydata = dAll.orig[c("ID", "session", "trnum", "taskname", "rep", "PHASE", "GROUP", myvar)]
mydata = mydata[complete.cases(mydata), ]
##################################
# fit stan model
##################################

stan_data <- list(
   NOBS = nrow(mydata),
   NSUBS = length(unique(mydata$ID)),  
   NTASKS = length(unique(mydata$taskname)), 
   NSESSIONS = length(unique(mydata$session)),
   y = unlist(mydata[myvar]),
   ID = as.numeric(mydata$ID),
   session = mydata$session,
   task = as.numeric(as.factor(mydata$taskname))
)


stan_control = list(adapt_delta=0.99,  stepsize = 0.01, max_treedepth = 15)
MODEL_FILE = 'learning_model.stan'
USE_VB = T
mypars = c("mu_gp","sigma_gp", "sigma", "lp__")
parlist = list()


if (!USE_VB){
   fit <- stan(file = MODEL_FILE, data = stan_data, 
               iter = 5000, chains = 4, control = stan_control, 
               warmup = 4000)
   
   pairs(fit, pars = mypars, las = 1)
   
   sampler_params <- get_sampler_params(fit, inc_warmup = TRUE)
   summary(do.call(rbind, sampler_params), digits = 2)
   Rhat  = summary(fit)$summary[, "Rhat"] 
   n_eff = summary(fit)$summary[, "n_eff"] 
   
   print(max(Rhat))
   print(min(n_eff))
   
   lapply(sampler_params, summary, digits = 2)
} else {  
   fit <- vb(stan_model(MODEL_FILE), data = stan_data, output_samples = 1000)
}

pars <- rstan::extract(fit, permuted = TRUE)

A = colMeans(pars$A) + mean(pars$A_task)
tau = colMeans(exp(pars$log_tau + mean(pars$log_tau_task)))
C = colMeans(pars$C) + mean(pars$C_task)


parlist[[1]] = list(A = A, tau = tau, C = C)

ids <- unique(mydata$ID)
id_index <- as.numeric(ids)

for (i in id_index) {
  plot_trace_pars(dAll.orig, id = ids[i], id_index = i, pars = parlist, doprint =T , rm.out = 5, vars = 'dev')
}

mystats = data.frame(ID = ids, ID_index = id_index, A = A, tau = tau, C = C)
mystats = mystats %>% mutate(SUBJECT = as.character(paste0('LH', ID)))
mystats$var = myvar

dAll <- dAll.orig %>% group_by (ID, taskname) %>%
   mutate (dev = log(normData (dev, session)), 
           meanpace = log(normData(1/meanspeed, session))
           ) %>%
   group_by(ID, taskname, session) %>%
   summarize (dev = mean(dev, na.rm = T), 
              meanpace = mean(meanpace, na.rm =T),
              c_session = median(c_session, na.rm = T)) 


# plot mean across all tasks and participants
plot_trace(dAll, doprint = T, rm.out = 5)

# plot mean across tasks for each participant
subjstats <- c()
for (id in ids) {
   stats <- plot_trace(dAll, id = id, doprint =T , rm.out = 5)
   subjstats <- rbind(subjstats, stats)
}

# plot mean across participants, for each task
#tasks <- levels(dAll$taskname)
#for (task in tasks) {
#   plot_trace(dAll, task = task, doprint =T, rm.out = 5 )
#}

# do not print to file but determine exponential fit for all ids and tasks
ids <- levels(dAll$ID)
#tasks <- levels(dAll$taskname)
allstats <- c()
#for (task in tasks) {
   for (id in ids) {
      print(paste(id))
#      print(paste(id, task))
      stats <- plot_trace(dAll, id = id, task = NULL, doprint = F, rm.out = 5)
      allstats <- rbind(allstats, stats)
   }
#}

summary(allstats$R2 )
allstats = allstats %>% mutate(SUBJECT = as.character(paste0('LH', ID)))


##################################
# get connectivity patterns
##################################

fc_evol = merge(fc_evol, data, by = c("SUBJECT", "SESSION"))

fc_evol.mean = fc_evol %>% group_by(SUBJECT, WEEK.MEAN, GROUP, PHASE) %>%
      summarise_at(vars(contains("fc_")), funs(mean(., na.rm = TRUE)))  
fc_evol.mean = fc_evol.mean[
      with(fc_evol.mean, order(SUBJECT, GROUP, WEEK.MEAN)),
      ]
   
fc_training = subset(fc_evol.mean, (GROUP == 'early' & PHASE == 1) | (GROUP == 'late' & PHASE == 2)) 
fc_training.1 = subset(fc_training, WEEK.MEAN ==1) 
fc_training.2 = subset(fc_training, WEEK.MEAN ==2) 
fc_training.3 = subset(fc_training, WEEK.MEAN ==3)
colnames(fc_training.1)[-seq(4)] = paste(colnames(fc_training.1)[-seq(4)], 1, sep = '.')
colnames(fc_training.2)[-seq(4)] = paste(colnames(fc_training.2)[-seq(4)], 2, sep = '.')
colnames(fc_training.3)[-seq(4)] = paste(colnames(fc_training.3)[-seq(4)], 3, sep = '.')

fc = merge(fc_training.1[-2], fc_training.2[-2], by = c('SUBJECT', 'GROUP', 'PHASE'))
fc = merge(fc, fc_training.3[-2], by = c('SUBJECT', 'GROUP','PHASE'))

##################################
# sync datasets
##################################

plot(subset(allstats, var == myvar)$tau, subset(mystats, var == myvar)$tau)
plot(subset(mystats, var == myvar)$A, subset(mystats, var == myvar)$tau)

myvar = 'dev'
mydata = merge(subset(allstats, var == myvar), fc, by = 'SUBJECT' )
valid = mydata$R2 > .2

##################################
# try prediction
##################################
fc_cols = grep( "fc_*", colnames(mydata))#[seq(105)]
y = mydata$tau
#y = mydata$asymp[valid]
#y[ y < 0] = 0
# do LOOCV
X = mydata[, fc_cols]

maxcomp = 10
steps = seq(0.1, 0.9, 0.2)

X.1 = X[mydata$GROUP == 'early', ]
X.2 = X[mydata$GROUP == 'late', ]
y.1 = y[mydata$GROUP == 'early']
y.2 = y[mydata$GROUP == 'late']

cv <- cv.spls( X.1, y.1, eta = steps, K = c(1:maxcomp), plot.it = T )
mypls <- spls( X.1, y.1, eta = cv$eta.opt, K = cv$K.opt  )
y.pred = predict(mypls, X.2)
cor.test(y.2, y.pred)
plot(y, y.pred)


cv <- cv.spls( X, y, eta = steps, K = c(1:maxcomp), plot.it = T )
mypls <- spls( X, y, eta = cv$eta.opt, K = cv$K.opt  )
y.pred = predict(mypls, X)

NITER = 10
y.pred = y.train.pred = c()
for (i in seq(length(y))){
   print(i)
   
   X.train = X[-i, ]
   X.test = X[i, ]
   y.train = y[-i]
   
   coefs.iter = y.pred.iter = y.train.pred.iter = NULL    
   for (iter in seq(NITER)){
      mysample = sort(sample(seq(length(y.train)), replace = T))
      cv <- cv.spls( X.train[mysample, ], y.train[mysample], eta = steps, K = c(1:maxcomp), plot.it = T )
      mypls <- spls( X.train[mysample, ], y.train[mysample], eta = cv$eta.opt, K = cv$K.opt  )
      
      coefs.iter = cbind(coefs.iter, coef.spls(mypls))
      y.pred.iter = cbind(y.pred.iter, predict(mypls, X.test))
#      y.train.pred.iter = cbind(y.train.pred.iter, predict(mypls, X.train))
   }
   
   coefs = rowMeans(coefs.iter)
   y.pred[i] = rowMeans(y.pred.iter)
#   y.train.pred[, i] = rowMeans(y.train.pred.iter)
      
}

cor.test(y, y.pred)
plot(y, y.pred)

