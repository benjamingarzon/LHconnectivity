
# drift?
# do 1-sided tests!!!!!!!!!!!!
# do only for rest...

setwd("~/Software/LeftHand/LHconnectivity/scripts/stats")
source('./installpackages.R')
source('~/Software/LeftHand/VBM/get_session_data.R')
source('./fc_analysis_func.R')

FIGS_DIR = '/home/share/LeftHand/LHconnectivity/figs/'

alpha = 0.05
dotest = select_dotest(dataset)

#cor.test(unlist(fc_evol.orig[-seq(4)]), unlist(fc_evol.denoised[-seq(4)]))
#plot(fc_evol.orig$fc_001, fc_evol.denoised$fc_001)


fc_evol$SUBJECT = sapply(as.character(fc_evol$SUBJECT), function(x) strsplit(x, '-')[[1]][2])

if (dataset == 'Stockholm') {  
  fc_evol = merge(fc_evol, data, by = c("SUBJECT", "SESSION"))
  selected = selected.Stockholm
  } else {
  selected = selected.Berlin    
}

ROInames = ROInames[selected == 1]
ROIcolors = ROIcolors[selected == 1]

selected = selected %*% t(selected)
diag(selected) = 0
selected.conn = squareform(selected)


fc_cols = grep( "fc_", colnames(fc_evol))[selected.conn == 1]

#plot(fc_evol$SESSION, rowMeans(fc_evol[fc_cols]))

fc_conn = squareform(colMeans(fc_evol[fc_cols], na.rm = T))
colnames(fc_conn) = rownames(fc_conn) = ROInames


##################
# univariate analysis
##################


# plot similarity between subjects/scans
image(cor(fc_evol[fc_cols]))


fc_cor = cor(t(fc_evol[fc_cols]))

fc_cor.melt = melt(fc_cor)

plot.cor <- ggplot(fc_cor.melt, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name="FC") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 5, hjust = 1)) +
  theme(axis.text.y = element_text(size = 5, hjust = 1)) +
  coord_fixed() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()
#    legend.justification = c(1, 0),
#    legend.position = c(0.6, 0.7),
#    legend.direction = "horizontal")+
#  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
#                               title.position = "top", title.hjust = 0.5))
)
#print(plot.cor)

print('Connection layout')
fc_names = squareform(as.numeric(substr(colnames(fc_evol)[fc_cols], 4, 6)))
colnames(fc_names) = rownames(fc_names) = ROInames  
# plot location of the connections
plot.names <- ggplot(melt(fc_names), aes(Var2, Var1, fill = value))+
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 7, hjust = 1)) +
  theme(axis.text.y = element_text(size = 7, hjust = 1)) +
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank())

print(plot.names)

#results = apply(fc_evol[fc_cols], 2, dotest, fc_evol[-fc_cols])
results = lapply(seq(length(fc_cols)), function(x) dotest(x, fc_evol[fc_cols], fc_evol[-fc_cols]))

p.values = sapply(results, function(x) x$p)
coefs = sapply(results, function(x) x$coef)

hist(p.values, 20, xlab = "p value")

p.values.fdr = p.adjust(p.values, method = "fdr")


# represent connectivity matrix 
fc_conn.melt = melt(get_upper_tri(round(fc_conn, 2)), na.rm = T)

plot.conn <- ggplot(fc_conn.melt, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1.5,1.5), space = "Lab", 
                       name="FC") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1)) +
  theme(axis.text.y = element_text(size = 12, hjust = 1)) +
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

# plot pvalues
p.values.mat = squareform(p.values)
coefs.mat = squareform(coefs)

diag(p.values.mat) = NA
colnames(p.values.mat) = rownames(p.values.mat) = ROInames

p.values.melt = melt(p.values.mat)
coefs.melt = melt(coefs.mat)

p.values.melt$sig = sign(coefs.melt$value)*(p.values.melt$value < alpha)

plot.pvals = ggplot(p.values.melt, aes(Var2, Var1, fill = sig)) + geom_tile() +  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) +
  theme(axis.text.y = element_text(size = 12)) +
  scale_fill_gradientn(values = c(1, 0), colours = c("red", "white", "cyan")) + 
  geom_text(aes(Var2, Var1, label = round(value, 2)), color = "black", size = 4)


# save plots
save_fig(paste("FC", suffix, sep = "_"))
print(plot.conn)
save_fig(paste("FC_pvals", suffix, sep = "_"))
print(plot.pvals)
dev.off()

##################
# clustering analysis
##################

if (dataset == 'Stockholm'){     #Berlin_difference' ){
  
fc.dem = sapply(results, function(x) x$y.dem)
fc_evol.dem = fc_evol  
fc_evol.dem[fc_cols] = fc.dem 


# cluster connections by trajectory
# aggregate by controls / patient
id.vars = c("SESSION", "GROUP", "REALWEEK", "PHASE")
NVARS = length(id.vars)
fc_evol.mean = fc_evol.dem %>% group_by(SESSION, REALWEEK, GROUP, PHASE) %>%
  summarise_at(vars(contains("fc_")), funs(mean(., na.rm = TRUE)))  
fc_evol.mean = fc_evol.mean[
  with(fc_evol.mean, order(GROUP, SESSION)),
  ]
# plot(fc_evol.mean$REALWEEK, fc_evol.mean$fc_002, pch = 20, ylim = c(-1, 1))
# #fc_evol.mean = detrend(fc_evol.mean, id.vars)
# #points(fc_evol.mean$REALWEEK, fc_evol.mean$fc_002, pch = 20, col = 'red')
# 
# #matplot(fc_evol.mean$REALWEEK, fc_evol.mean[-seq(NVARS)])
cluster_colors = c('red', 'blue', 'green', 'yellow', 'grey', 'cyan','black','orange','magenta', 'brown')
X = t(subset(fc_evol.mean, SESSION <= 5)[-seq(NVARS)])
 
mydist = as.dist(1 - cor(t(X)))
clustmethod = "complete" # ward.D2

#mydist = dist(X, method = "euclidean")

heatmap(X, distfun = function(x) as.dist(1 - cor(t(X))), #dist(x, method = 'euclidean'),
        hclustfun = function(d) hclust(d, method = clustmethod), Colv = NA)

res.hc <- hclust(mydist, method = clustmethod)     # Compute hierachical clustering  

K = 20
fviz_dend(res.hc, k = K, # Cut in four groups
          cex = 0.5, # label size
          k_colors = c(cluster_colors),
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE # Add rectangle around groups
)

# # analyze different cluster numbers
fc_evol.clust = fc_evol.mean[seq(NVARS)]
memb <- cutree(res.hc, k = K)
for (i in 1:K){
  fc_evol.clust[ paste0("cluster", i)] = rowMeans(fc_evol.mean[-seq(NVARS)][, memb==i])
  
}

fc_evol.melt = melt(fc_evol.clust, id.vars = id.vars)
plot.clus <- ggplot(fc_evol.melt, aes(x = REALWEEK, y = value, group = GROUP, col = GROUP)) + geom_line() + facet_grid( . ~ variable)
print(plot.clus)

table(memb)

}
#fc_evol.melt = melt(fc_evol.detrended, id.vars = id.vars)
#fc_evol.melt$cluster = memb[fc_evol.melt$variable]
#plot.clus <- ggplot(fc_evol.melt, aes(x = SESSION, y = value, col = GROUP)) + geom_point() + facet_grid( . ~ cluster)
#table(memb)
#print(plot.clus)
# 
# library(igraph)
# 
# #  diag(p.values.mat) = 1
# #  colnames(p.values.mat) = rownames(p.values.mat) = ROInames
# memb.mat = squareform(memb)
# memb.mat[memb.mat %in% c(0, 2, 5)  ] = NA
# colnames(memb.mat) = rownames(memb.mat) = ROInames
# edge.color = c('red', 'black', 'blue', 'green', 'cyan')[memb.mat]
# #net <- graph_from_adjacency_matrix(!is.na(memb.mat))
# 
# net.mat = exp(fc_conn.all)
# net.mat[is.na(memb.mat)] = NA
# net <- graph_from_adjacency_matrix(net.mat)
# 
# l = layout.fruchterman.reingold(net)
# 
# plot(net, edge.color = edge.color, edge.width = 3, vertex.size = 4, layout = l, edge.curved = 0.3, arrow.mode = 0, edg.lty = 2)
# 
# plot.net <- ggplot(melt(edge.color), aes(Var2, Var1, fill = value)) + geom_tile() +  theme_minimal()+ # minimal theme
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, 
#                                    size = 12, hjust = 1))
# 
# 
# print(plot.net)
# 
# #  p.values.all.mat = squareform(p.values.all)
# 
# #  p.values.rest.fdr = p.adjust(p.values.rest, method = "fdr")
# 
# #  p.values.mat = squareform(p.values.rest)
# #  diag(p.values.mat) = 1
# #  colnames(p.values.mat) = rownames(p.values.mat) = ROInames
# #  p.values.net = network::network(p.values.mat, directed = F)
# 
# 
# #colnames(p.values.rest.net) = c("from_id", "to_id", "value")
# #  fc_conn.net = network::network(exp(fc_conn.rest), directed = F)
# 
# #  alpha = 0.05
# #set.edge.attribute(p.values.net, "size", ifelse(p.values.mat < alpha), 1, 0.001)
# 
# #plot.net = ggnet2(net = p.values.net, label = T,  color = ROIcolors, 
# #                  mode = "fruchtermanreingold", label.size = 5,  edge.size = "size",
# #                  label.color = "black", layout.exp = 0.5)
# 
# #plot.net <- ggplot(melt(p.values.mat < alpha), aes(Var2, Var1, fill = value)) + geom_tile() +  theme_minimal()+ # minimal theme
# #theme(axis.text.x = element_text(angle = 45, vjust = 1, 
# #                                 size = 12, hjust = 1))
# 
# #print(plot.net)
# 
# 
# get_upper_tri <- function(cormat){
#   cormat[lower.tri(cormat)]<- NA
#   diag(cormat) = 1
#   return(cormat)
# }
# 
# fc_conn.melt = melt(get_upper_tri(round(fc_conn.rest, 2)), na.rm = T)
# 
# plot.conn <- ggplot(fc_conn.melt, aes(Var2, Var1, fill = value))+
#   geom_tile(color = "white")+
#   scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
#                        midpoint = 0, limit = c(-1,1), space = "Lab", 
#                        name="FC") +
#   theme_minimal()+ # minimal theme
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, 
#                                    size = 12, hjust = 1)) +
#   coord_fixed() +
#   geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) +
#   theme(
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.border = element_blank(),
#     panel.background = element_blank(),
#     axis.ticks = element_blank(),
#     legend.justification = c(1, 0),
#     legend.position = c(0.6, 0.7),
#     legend.direction = "horizontal")+
#   guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
#                                title.position = "top", title.hjust = 0.5))
# 
# 
# print(plot.conn)
# 
# 
# 
# if (F){
#   
#   dotest_global = function(y, myfc){
#     myfc$y = y
#     var.names = c("Intercept", 
#                   "SESSION.MEAN.4","SESSION.MEAN.8","SESSION.MEAN.13","SESSION.MEAN.17","SESSION.MEAN.21",
#                   "GROUPearly",
#                   "SESSION.MEAN.4_x_GROUPearly","SESSION.MEAN.8_x_GROUPearly","SESSION.MEAN.13_x_GROUPearly",
#                   "SESSION.MEAN.17_x_GROUPearly","SESSION.MEAN.21_x_GROUPearly", "FD")
#     
#     contrast.names = c("GROUPearly", 
#                        "SESSION.MEAN.4_x_GROUPearly", "SESSION.MEAN.8_x_GROUPearly","SESSION.MEAN.13_x_GROUPearly",
#                        "SESSION.MEAN.17_x_GROUPearly","SESSION.MEAN.21_x_GROUPearly")
#     
#     c.1       = c(0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0)
#     c.2       = c(0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)
#     c.3       = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0)
#     c.4       = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0)
#     c.5       = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0)
#     c.6       = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)
#     
#     cont.mat = rbind(c.1, c.2, c.3, c.4, c.5, c.6)
#     colnames(cont.mat) = var.names
#     rownames(cont.mat) = contrast.names
#     
#     model = lmer(y ~ 1 + as.factor(SESSION.MEAN)*as.factor(GROUP) + FD + (1 |SUBJECT), data = myfc, REML = T)
#     #print(model)
#     
#     coefs = fixef(model)
#     pvalues = summary(model)$coefficients[ , "Pr(>|t|)"]
#     
#     glh = glht(model, linfct = cont.mat, alternative="greater")
#     contrast.pvalues = summary(glh, test = Ftest())$test$pvalue[[1]]
#     
#     contrast.coefs = coef(glh)
#     
#     return(contrast.pvalues)
#   }
#   
#   
#   p.values.left = apply(fc_evol_left[fc_cols], 2, dotest_global, fc_evol_left[-fc_cols])
#   
#   which(p.values.left < 0.01)
#   
#   fc.run = fc %>% group_by(SUBJECT, WEEK.MEAN, PHASE, block, GROUP) %>% summarise(
#     fc_025 = mean(fc_025, na.rm=T)
#   )
#   
#   fc.mean = fc %>% group_by(WEEK.MEAN, PHASE, block, GROUP) %>% summarise(
#     fc_025 = mean(fc_025, na.rm=T)
#   )   
#   
#   
#   plot.p = ggplot(data = fc.mean, aes(x = WEEK.MEAN, y = fc_025, col = GROUP)) + geom_line(lwd = 2) + 
#     geom_line(data = fc.run, aes(x = WEEK.MEAN, y = fc_025, group = SUBJECT, col = GROUP), lwd = 0.5, lty = 2) +   
#     facet_grid(block ~  PHASE)
#   
#   print(plot.p)
#   
#   #p.values.left = apply(fc13.left[fc_cols], 2, dotest, fc13.left[-fc_cols])
#   #p.values.right = apply(fc13.right[fc_cols], 2, dotest, fc13.right[-fc_cols])
#   #p.values.rest = apply(fc13.rest[fc_cols], 2, dotest, fc13.rest[-fc_cols])
#   
#   par(mfrow=c(3, 1))
#   fc13.left = subset(fc, SESSION %in% c(1, 3) & block=="left")
#   fc13.right = subset(fc, SESSION %in% c(1, 3) & block=="right")
#   fc13.rest = subset(fc, SESSION %in% c(1, 3) & block=="rest")
#   fc13.rest = subset(fc, SESSION %in% c(1, 3) & block=="rest")
#   
#   
# }

