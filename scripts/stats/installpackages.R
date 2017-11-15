
libs = c('lmerTest', 'dplyr', 'reshape2', 'ggplot2', 'lme4', 'lmerTest', 'multcomp', 'GGally', 'pracma', 'geomnet', 'ggnetwork', 'network',
         'tsne', 'factoextra', 'cluster')

for (i in libs){
  if( !is.element(i, .packages(all.available = TRUE)) ) {
    install.packages(i)
  }
  library(i,character.only = TRUE)
}

c('lmerTest', 'dplyr', 'reshape2', 'ggplot2', 'lme4', 'lmerTest', 'multcomp', 'GGally', 'pracma', 'geomnet', 'ggnetwork', 'network',
'tsne', 'factoextra', 'cluster')