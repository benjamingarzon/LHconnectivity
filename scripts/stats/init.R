# init.R
# ======================
#
# load libraries, define aux functions, analysis settings
#


############################
# import libraries
############################

library(dplyr)

#source('../../useful/myanovaplot.r');
#source('../../useful/epsilon.r');
#source('../../useful/latex_IO.r');
#source('../../useful/julstats.r');
#source('../../useful/julgraph.r');



############################
# auxiliary functions
############################

# compute zscore 
zscore <- function (a) {(a-mean(a, na.rm =T))/sd(a, na.rm=T) }

# counts number of no-NaN elements
count<-function(x) length(which(!is.nan(x)))

# wrapper for paste function
ps <- function (...) {paste(...,sep='')}


# determine outliers (separately for each group, if indicated), 
# based on mad criterion (default = 3 mad)
outliers_mad <- function (d,groups=NULL,nmad = 3) {  
   out <- rep(TRUE,length(d))
   if (is.null(groups)) { groups <- factor(rep('a',length(d)))}
   
   for (l in levels(groups)) {
      temp <- (groups==l);
      m <- median(d[temp],na.rm = T);
      s <- mad(d[temp],na.rm = T);
      out[temp] <- abs(d[temp]-m)>nmad*s;
   }	
   out
}


# plot dvs for trace data
#
# arguments
#   d = data frame
#   id = id (or NULL --> take mean across all ids)
#   task = task (or NULL --> take mean across all tasks)
#   doprint = TRUE --> create plot file
#   rm.out = FALSE -> don't remove outliers, value n > 0 --> remove outliers with n*mad criterion
#   ylim = limits for vertical axis(default: c(-0.5, 0.1))
#   subdir = subdir to which to plot (if doprint == T); defaults to 'plots'
#   vars = dvs to be plotted (should be 4, to fit in 2x2 plot)
#
# returns
#   data frame with statistics (fit info, e.g., R2, for each of the dvs)
plot_trace <- function(d, id = NULL, task = NULL, doprint = F, rm.out = F, ylim = c(-0.5, 0.3),
                       subdir = 'plots', 
                       vars = c('dev', 'meanpace', 'c_session')
) {
   
   # select (if necessary)
   if  (!is.null(task)) {
       d <- filter(d, taskname == task)
   } else {
      task = 'meanTask'
   }
      
   if  (!is.null(id)) {
      d <- filter(d, ID == id)
   }  else {
      id <- 'meanID'
   }
   
   # aggregate (if necessary)
   d <- d %>% group_by(session) %>% summarise_at(vars,  function(d) {mean(d, na.rm = T)})
   
   if (doprint) {
      png(ps(subdir, '/trace_',id,'_',task,'.png'), height=20, width=20,res=300,units='cm')
   }
   
   par(mfrow=c(2,2))

   stats <- c()
   t <- as.numeric(d$session)
   for (k in 1:length(vars)) {
      vn <- vars[k]
      y <- d[[vn]]
      
      if (rm.out) { # remove outliers?
         out <- outliers_mad(y, nmad = rm.out)
         y[out] <- NaN
      }
      fm <- fit_explearn(y,t, doplot=F)
      plot(y~t,ylab='a.u.', main=vars[k],ylim = ylim, xlab='session')
      lines(fm$yy~fm$t)
      text(mean(t),ylim[2],adj = c(0.5, 1),sprintf('tau=%.2f, asymp=%.2f\nR2 = %.2f', fm$cc, fm$aa ,fm$R2))
      
      # collect stats info
      newline <- data.frame(ID = id, task = task, var = vars[k], R2 = fm$R2, asymp = fm$aa, tau = fm$cc, fac = fm$bb)
      stats <- rbind(stats, newline)
   }
   if (doprint) {
      dev.off()
   }
   stats
}  

plot_trace_pars <- function(d, id = NULL, id_index = NULL, task = NULL, doprint = F, rm.out = F, ylim = c(-0.5, 0.3),
                       subdir = 'plots', pars = NULL,
                       vars = c('dev', 'meanpace', 'c_session')
) {
   
   # select (if necessary)
   if  (!is.null(task)) {
      d <- filter(d, taskname == task)
   } else {
      task = 'meanTask'
   }
   
   if  (!is.null(id)) {
      d <- filter(d, ID == id)
   }  else {
      id <- 'meanID'
   }
   
   # aggregate (if necessary)
   d <- d %>% group_by(session) %>% summarise_at(vars,  function(d) {mean(d, na.rm = T)})
   
   if (doprint) {
      png(ps(subdir, '/trace_',id,'_',task,'_1.png'), height=20, width=20, res=300, units='cm')
   }
   
   par(mfrow=c(2,2))
   
   t <- as.numeric(d$session)
   for (k in 1:length(vars)) {
      vn <- vars[k]
      y <- d[[vn]]
      
      if (rm.out) { # remove outliers?
         out <- outliers_mad(y, nmad = rm.out)
         y[out] <- NaN
      }

      plot(y ~ t,ylab='a.u.', main=vars[k], xlab='session') #ylim = ylim,

      A = pars[[k]]$A[id_index]
      tau = pars[[k]]$tau[id_index]
      C = pars[[k]]$C[id_index]
      tt = (seq(min(t), max(t)))
      ff = A*exp(-(tt - 1)/tau) + C
      lines( tt,  ff)
      text(mean(t), -diff(range(y))*.1 + max(y), adj = c(0.5, 1),sprintf('A=%.2f, tau=%.2f, asymp=%.2f', A, tau ,C))
       
   }
   if (doprint) {
      dev.off()
   }
}  


# fit exponential ("learning") curve to data
#
# arguments
#   y = data
#   t = time/session info
#
# returns 
#   list with fit statistics (R2, asymptote, time constant tau)
fit_explearn <- function(y,t, doplot=F) {
   
   minls <- function(cc) {
      temp = lm (y~exp(-(t-1)/cc))
      sum(temp$residuals^2)
   }
   
   opt = optimize(minls, lower=0, upper=100)
   cc = opt$minimum
   fm = lm (y~exp(-(t-1)/cc))
   aa = fm$coefficients[[1]]
   bb = fm$coefficients[[2]]
   fmsumm = summary(fm)
   yy = aa+bb*exp(-(t-1)/cc)
   if (doplot) {
      plot(y~t, type='b', pch = 20)
      lines(yy~t,col='green')
   }
   res=c()
   res$t = t
   res$y = y
   res$yy = yy
   res$aa = aa
   res$bb = bb
   res$cc = cc
   res$res = opt$objective
   res$R2 = fmsumm$r.squared
   res$fn = 'y ~ aa+bb*exp(-(t-1)/cc)'
   res
}

