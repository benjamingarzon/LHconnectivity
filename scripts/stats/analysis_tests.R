
do_tests_LME_SESSION_Berlin = function(y, X, ALTERNATIVE = "greater")
{
   var.names = c("Intercept", "SESSION", "GROUPtraining", 
                 "SESSION_x_GROUPtraining")
   
   contrast.names = c("contrast_SESSION_x_GROUPtraining")
   
   c.1       = c(0, 1, 0, 0)
   c.2       = c(0, 0, 0, 1)
   
   #   cont.mat = rbind(c.1, c.2)
   cont.mat = rbind(c.2)
   colnames(cont.mat) = var.names
   rownames(cont.mat) = contrast.names
   
   tags = c(paste0(var.names, '_coef'),
            paste0(contrast.names, '_coef'),
            paste0(var.names, '_p'),
            paste0(contrast.names, '_p'))
   
   X$y = y
   X = subset(X, SESSION < 20 & GROUP != 'none' )
   
   if( sum(y == 0) == 0){
      model = lmer(y ~ 1 + SESSION*GROUP + (1|SUBJECT), data= X)
      coefs = fixef(model)
      
      pvalues = summary(model)$coefficients[-5 , "Pr(>|t|)"]
      
      glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE)
      contrast.pvalues = summary(glh)$test$pvalues
      contrast.coefs = coef(glh)
      val = c(coefs, contrast.coefs, pvalues,  contrast.pvalues)
      
   } else {
      val = c(rep(0, ncol(cont.mat)), rep(0, nrow(cont.mat)), 
              rep(1, ncol(cont.mat)), rep(1, nrow(cont.mat)))
      
   }
   
   print(val)
   print(tags)
   names(val) = tags
   return(val)
   
}

do_tests_LME_WEEK = function(y, X, ALTERNATIVE = "greater")
{
   var.names = c("Intercept", "WEEK","GROUPearly", 
                 "WEEK_x_GROUPearly")
   
   #contrast.names = c("contrast_WEEK", "contrast_WEEK_x_GROUPearly")
   contrast.names = c("contrast_WEEK_x_GROUPearly")
   
   c.1       = c(0, 1, 0, 0)
   c.2       = c(0, 0, 0, 1)
   
   #   cont.mat = rbind(c.1, c.2)
   cont.mat = rbind(c.2)
   
   colnames(cont.mat) = var.names
   rownames(cont.mat) = contrast.names
   
   tags = c(paste0(var.names, '_coef'),
            paste0(contrast.names, '_coef'),
            paste0(var.names, '_p'),
            paste0(contrast.names, '_p'))
   
   X$y = y
   
   X = subset(X, PHASE == 1 & WEEK < 8)
   if( sum(y == 0) == 0){
      
      model = lmer(y ~ 1 + WEEK*GROUP + (1|SUBJECT), data= X)
      coefs = fixef(model)
      
      pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
      
      glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE)
      contrast.pvalues = summary(glh)$test$pvalues
      contrast.coefs = coef(glh)
      val = c(coefs, contrast.coefs, pvalues,  contrast.pvalues)
   } else {
      val = c(rep(0, ncol(cont.mat)), rep(0, nrow(cont.mat)), 
              rep(1, ncol(cont.mat)), rep(1, nrow(cont.mat)))
      
   }
   
   print(val)
   print(tags)
   names(val) = tags
   return(val)
   
}


do_tests_LME_TRAINING = function(y, X, ALTERNATIVE = "greater")
{
  var.names = c("Intercept", "TRAINING","GROUPearly", 
                "TRAINING_x_GROUPearly", "FD")
  
  #contrast.names = c("contrast_WEEK", "contrast_WEEK_x_GROUPearly")
  contrast.names = c("contrast_TRAINING_x_GROUPearly")
  
  c.1       = c(0, 1, 0, 0, 0)
  c.2       = c(0, 0, 0, 1, 0)
  
  #   cont.mat = rbind(c.1, c.2)
  cont.mat = rbind(c.2)
  
  colnames(cont.mat) = var.names
  rownames(cont.mat) = contrast.names
  
  tags = c(paste0(var.names, '_coef'),
           paste0(contrast.names, '_coef'),
           paste0(var.names, '_p'),
           paste0(contrast.names, '_p'))
  
  X$y = y
  
  X = subset(X, PHASE == 1)# & WEEK < 8)
  if( sum(y == 0) == 0){
    
    model = lmer(y ~ 1 + TRAINING*GROUP + FD + (1|SUBJECT), data= X)
    coefs = fixef(model)
    
    pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
    
    glh = glht(model, linfct = cont.mat, alternative=ALTERNATIVE)
    contrast.pvalues = summary(glh)$test$pvalues
    contrast.coefs = coef(glh)
    val = c(coefs, contrast.coefs, pvalues,  contrast.pvalues)
  } else {
    val = c(rep(0, ncol(cont.mat)), rep(0, nrow(cont.mat)), 
            rep(1, ncol(cont.mat)), rep(1, nrow(cont.mat)))
    
  }
  
  print(val)
  print(tags)
  names(val) = tags
  return(val)
  
}
