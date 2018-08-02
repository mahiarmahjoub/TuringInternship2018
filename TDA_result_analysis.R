## ======================= analyse & visualise links ======================= ##


dir <- "C://Users/mahia/Google Drive/MPhil/CompBio/Internship"
require(stringi, quietly = TRUE)
require(GGally, quietly = TRUE)
require(network, quietly = TRUE)
require(sna, quietly = TRUE)
require(ggplot2, quietly = TRUE)
require(RColorBrewer, quietly = TRUE)
require(intergraph, quietly = TRUE)
require(ndtv, quietly = TRUE)
require(igraph, quietly = TRUE)
require(gplots, quietly = TRUE)
require(randomcoloR, quietly = TRUE)
require(gridExtra, quietly = TRUE)
source(file.path(dir,"TuringInternship2018/network_stats_functions.R"))



## ==== load expression matrix, wanted gene list etc. ---------------------- ##
gene.list.GO <- read.csv(file.path(dir,"changing_genes_GO_list_v3_180702.csv"), 
                         skip=1, stringsAsFactors = FALSE, sep = ",", 
                         header=FALSE)
gene.list <- read.csv(file.path(dir,"changing_genes_list_v3_180702.csv"), 
                      header = TRUE, stringsAsFactors = FALSE, sep = ",")
TF.list <- read.csv(file.path(dir,"changing_genes_list_v3_180702_TF.csv"), 
                    header = TRUE, stringsAsFactors = FALSE, sep = ",")
target.list <- read.csv(file.path(dir,"changing_genes_list_v3_180702_targets.csv"), 
                        header = TRUE, stringsAsFactors = FALSE, sep = ",")
gene.alias <- read.delim(file.path(dir,"TAIR_genealiases.txt"), header=TRUE, 
                         stringsAsFactors = TRUE)
expr.all <- read.table(file = file.path(dir, "dawnburst_slcu", 
                                        "rawcounts", 
                                        "dawnburst_masterdf.txt"), 
                       header = TRUE, stringsAsFactors = FALSE
)
edges.aranet <- read.table(file.path(dir,"AraNet-network.txt"), header=FALSE, 
                           stringsAsFactors = FALSE, 
                           col.names = c("regulator", "target", "weight"))
genes.of.interest <- read.csv(file.path(dir,"CircadianClock_genelist_literature.csv"), 
                              header=TRUE, stringsAsFactors = FALSE)





## choose the right expression set for analysis 
temp <- "22"; mut <- "Col"
time <- c(-10, 10, 24, 30, 45, 60, 105, 120)
# pick selected experimental conditions
expr <- expr.all[,grepl(temp,colnames(expr.all)) & grepl(mut,colnames(expr.all))]

if(mut=="Ler" & temp=="22"){time <- time[-2]; expr <- expr[,-ncol(expr)]}
if(mut=="prr579" & temp=="22"){expr <- expr[,-ncol(expr)]}
if(mut=="prr579" & temp=="27"){time <- time[1:6]}
if(mut=="phyA" & temp=="22"){expr <- expr[,1:8]}
if(mut=="phyA" & temp=="27"){expr <- expr[,1:8]}

colnames(expr) <- as.character(time)
rownames(expr) <- expr.all$tracking_id
expr.selected <- as.matrix(expr[gene.list$x,])  # pick selected genes



## ==== dynGENIE3 ---------------------------------------------------------- ##

res.genereg <- read.csv(file.path(dir,"TDA_output",
                                  paste0("geneReg_links_180702_",mut,"_",temp,".csv")), 
                        header = TRUE, stringsAsFactors = FALSE)[,2:6]
res.dynGENIE <- read.csv(file.path(dir, "TDA_output",
                                   paste0("dynGENIE_links_180702_v3_",mut,"_",temp,".csv")),
                         header = TRUE, stringsAsFactors = FALSE)

# plot the histogram of weights to assess the selection of a threshold
par(mfrow=c(2,2))
barplot(log10(hist(res.dynGENIE$weight, plot = FALSE)$counts+1),
        main="dynGENIE weight histogram (log)",
        xlab="weight", ylab="log(frequency)", las=2,
        names.arg = hist(res.dynGENIE$weight, plot = FALSE)$breaks[-1]
)
hist(res.genereg$coef,  breaks = 30, main="GeneReg coef histogram",
     xlab="coefficient")
hist(res.genereg$delay,  breaks = 30, main="GeneReg delay histogram",
     xlab="delay (mins)")
hist(res.genereg$adj.r.squared,  breaks = 150,
     main="GeneReg adjusted R^2 histogram", xlab="adjusted R squared")
par(mfrow=c(1,1))


threshold.dynGENIE <- 0.02
threshold.genereg <- 0.9999999

# exclude edges below cutoff
res.dynGENIE.filt <- res.dynGENIE[res.dynGENIE$weight > threshold.dynGENIE,]
res.genereg.filt <- res.genereg[abs(res.genereg$coef) > threshold.genereg,]


# pdf(file.path(dir,"ExprPlots_dynGENIE_edges_filt_v3_splined.pdf"))
# for(i in 1:nrow(res.dynGENIE.filt)){
#   y.range <- range(c(expr.selected[res.dynGENIE.filt[i,1],],
#                expr.selected[res.dynGENIE.filt[i,2],]))
#   plot(time, expr.selected[res.dynGENIE.filt[i,1],], type="b", lwd=2, col="blue",
#        ylim=y.range, ylab= "expression (CPM)", xlab="time (mins)",
#        main=paste0(res.dynGENIE.filt[i,1],"->",res.dynGENIE.filt[i,1]))
#   lines(time, expr.selected[res.dynGENIE.filt[i,2],], type="b", lwd=2,
#         col="red", lty=2)
#   legend("topright", legend = c("from","to"), lty = c(1,2), lwd=c(2,2),
#          col = c("blue","red"), horiz = TRUE)
# }
# dev.off()






## ==== set up IGRAPH object ----------------------------------------------- ##
## IGRAPH object 
vertices.available.dynGENIE <- unique(c(res.dynGENIE.filt[,1],
                                        res.dynGENIE.filt[,2]))
igraph.dynGENIE <- graph.data.frame(res.dynGENIE.filt, 
                                    vertices.available.dynGENIE, 
                                    directed=T)
igraph.dynGENIE <- induced_subgraph(igraph.dynGENIE, 
                                    components(igraph.dynGENIE)$membership==1)


vertices.available.genereg <- unique(c(res.genereg.filt[,1],
                                       res.genereg.filt[,2]))
igraph.genereg <- graph.data.frame(res.genereg.filt, 
                                   vertices.available.genereg, 
                                   directed=T)







## ==== (1) analyse NETWORK ------------------------------------------------ ##

## ---- obtain quantitative information first
# DEGREE 
V(igraph.dynGENIE)$"degree" <- igraph::degree(igraph.dynGENIE, 
                                              V(igraph.dynGENIE)$name)
V(igraph.dynGENIE)$"degree.in" <- igraph::degree(igraph.dynGENIE, 
                                                 V(igraph.dynGENIE)$name,
                                                 mode = "in")
V(igraph.dynGENIE)$"degree.out" <- igraph::degree(igraph.dynGENIE, 
                                                  V(igraph.dynGENIE)$name,
                                                  mode = "out")
# EDGE WEIGHT
E(igraph.dynGENIE)$"edgeweight" <- E(igraph.dynGENIE)$weight
# assign TF or other to each gene 
V(igraph.dynGENIE)$"TF" <- ifelse(V(igraph.dynGENIE)$name %in% TF.list$x, "TF","other")
# COMMUNITIES
# directed (walktrap)
community.walktrap <- igraph::cluster_walktrap(igraph.dynGENIE)
communities.dir <- communities(community.walktrap)
V(igraph.dynGENIE)$"membership.dir" <- membership(community.walktrap)
# undirected (louvin)
community.undir <- igraph::cluster_louvain(igraph::as.undirected(igraph.dynGENIE))
communities.undir <- communities(community.undir)
V(igraph.dynGENIE)$"membership.undir" <- membership(community.undir)
# BETWEENNESS
V(igraph.dynGENIE)$"betweenness" <- igraph::betweenness(igraph.dynGENIE) # vertex
E(igraph.dynGENIE)$"betweenness" <- igraph::edge.betweenness(igraph.dynGENIE) # edge
# CLUSTERING COEFFICIENT
V(igraph.dynGENIE)$"clustCoeff" <- igraph::transitivity(as.undirected(igraph.dynGENIE), 
                                                        type = "local")
# CLOSENESS CENTRALITY
V(igraph.dynGENIE)$"closnessCentrality" <- igraph::centralization.closeness(igraph.dynGENIE)$res
# EIGENVECTOR CENTRALITY
V(igraph.dynGENIE)$"EigenCentrality" <- igraph::eigen_centrality(igraph.dynGENIE)$vector
# HUB SCORE
V(igraph.dynGENIE)$"hubscore" <- igraph::hub.score(igraph.dynGENIE)$vector
# PEAK EXPRESSION TIME 
V(igraph.dynGENIE)$"peak.expr.time" <- temporal.order.genes(expr.selected[V(igraph.dynGENIE)$name,])$peak.time
peak.expr.time <- temporal.order.genes(expr.selected[V(igraph.dynGENIE)$name,])



## ---- assign colours/sizes based on quantitative measures 
degree.threshold <- 20

# change the ALPHA of NODES according to their degree
# degree > threhold --> alpha=0.95
V(igraph.dynGENIE)$"node.alpha.degree" <- ifelse(V(igraph.dynGENIE)$"degree" > degree.threshold,
                                                 0.9,0)
V(igraph.dynGENIE)$"node.alpha.degree"[
  V(igraph.dynGENIE)$"degree" <= degree.threshold
  ] <- seq(0.01,0.8,length.out = degree.threshold)[
    V(igraph.dynGENIE)$"degree"[
      V(igraph.dynGENIE)$"degree" <= degree.threshold]
    ]

# change the SIZE of NODES according to their degree
V(igraph.dynGENIE)$"node.size.degree" <- ifelse(V(igraph.dynGENIE)$"degree" > degree.threshold,
                                                (V(igraph.dynGENIE)$"degree"+40)/5,
                                                V(igraph.dynGENIE)$"degree"/5)

# set the NODE COLOUR to BLUE if TF
V(igraph.dynGENIE)$"node.color.TF" <- ifelse(V(igraph.dynGENIE)$name %in% TF.list$x, 
                                             brew_colors(col = "blue"),
                                             brew_colors(col = "red"))

# set the NODE SHAPE to 3 if TF
V(igraph.dynGENIE)$"node.shape.TF" <- ifelse(V(igraph.dynGENIE)$name %in% TF.list$x,15,16) 


# set the NODE COLOUR according to PEAK EXPRESSION TIME
V(igraph.dynGENIE)$"node.color.exptime" <- colorRampPalette(
  c("blue","lightgreen","yellow"))(length(time))[peak.expr.time$peak.time]

# set the NODE COLOUR to their respective COMMUNITY/MEMBERSHIP
set.seed(3246)
V(igraph.dynGENIE)$"node.col.membership.dir" <- distinctColorPalette(
  length(unique(V(igraph.dynGENIE)$"membership.dir"))
)[V(igraph.dynGENIE)$"membership.dir"]
set.seed(3246)
V(igraph.dynGENIE)$"node.col.membership.undir" <- distinctColorPalette(
  length(unique(V(igraph.dynGENIE)$"membership.undir"))
)[V(igraph.dynGENIE)$"membership.undir"]

# set NODE COLOUR to BETWEENNESS
V(igraph.dynGENIE)$"node.color.betweenness" <- ifelse(V(igraph.dynGENIE)$"betweenness"==0,
                                                      brew_colors(col = "orange"),
                                                      brew_colors(col = "blue"))
# set NODE SIZE to BETWEENNESS
V(igraph.dynGENIE)$"node.size.betweenness" <- ifelse(log(V(igraph.dynGENIE)$"betweenness") > 1,
                                                     (log(V(igraph.dynGENIE)$"betweenness"))+10,
                                                     1)

# set NODE SHAPE to BETWEENNESS + DEGREE
V(igraph.dynGENIE)$"node.shape.betweenDegree" <- ifelse(log(V(igraph.dynGENIE)$"betweenness")> 1,18,16)
V(igraph.dynGENIE)$"node.shape.betweenDegree"[V(igraph.dynGENIE)$"node.shape.betweenDegree"==16] <- 
  ifelse(V(igraph.dynGENIE)$"degree"[V(igraph.dynGENIE)$"node.shape.betweenDegree"==16]> degree.threshold,15,16)

# set NODE SIZE to BETWEENNESS + DEGREE
V(igraph.dynGENIE)$"node.size.betweenDegree" <- ifelse(log(V(igraph.dynGENIE)$"betweenness") > 1,
                                                       (log(V(igraph.dynGENIE)$"betweenness"))+10,
                                                       1)
V(igraph.dynGENIE)$"node.size.betweenDegree"[V(igraph.dynGENIE)$"node.size.betweenDegree"==1] <- 
  ifelse(V(igraph.dynGENIE)$"degree"[V(igraph.dynGENIE)$"node.size.betweenDegree"==1] > degree.threshold,
         (V(igraph.dynGENIE)$"degree"[V(igraph.dynGENIE)$"node.size.betweenDegree"==1]+40)/5,
         V(igraph.dynGENIE)$"degree"[V(igraph.dynGENIE)$"node.size.betweenDegree"==1]/5)

# set NODE SIZE basde on clustCoeff*100
V(igraph.dynGENIE)$"node.size.clust" <- sapply(V(igraph.dynGENIE)$"clustCoeff",
                                               function(x){
                                                 if(is.nan(x)){
                                                   5
                                                 } else if (x > 0){
                                                   x*500
                                                 } else if (x==0){
                                                   10
                                                 }
                                               })

# set NODE COLOUR on CLUST.COEFF
V(igraph.dynGENIE)$"node.color.clust" <- sapply(V(igraph.dynGENIE)$"clustCoeff",
                                                function(x){
                                                  if(is.nan(x)){
                                                    brew_colors("red")
                                                  } else if (x > 0){
                                                    brew_colors("orange")
                                                  } else if (x==0){
                                                    brew_colors("blue")
                                                  }
                                                })

# set NODE ALPHA on CLUST.COEFF
V(igraph.dynGENIE)$"node.alpha.clust" <- unlist(sapply(V(igraph.dynGENIE)$"clustCoeff",
                                                       function(x){
                                                         if(is.nan(x)){
                                                           0.05 # nan
                                                         } else if (x >= 0.3){
                                                           1 # orange
                                                         } else if (x < 0.3){
                                                           0.5 # blue
                                                         }
                                                       }
))

# set NODE ALPHA to log(EigenCentrality)
V(igraph.dynGENIE)$"node.alpha.EigCent" <- sapply(log10(V(igraph.dynGENIE)$"EigenCentrality"),
                                                  function(x){
                                                    if(x==0){
                                                      1
                                                    } else if (x<0 & x>=-1){
                                                      0.9
                                                    } else if (x<-1 & x>=-2){
                                                      0.5
                                                    } else if (x<-2 & x>=-4){
                                                      0.2
                                                    } else if (x<-4 & x>=-6){
                                                      0.05
                                                    }
                                                  }
)






# set the edge colour proportional to its weight 
E(igraph.dynGENIE)$"edgecolor" <- colorRampPalette(c("black", 
                                                     "grey80"))(length(E(
                                                       igraph.dynGENIE)$weight))

# set EDGE THICKNESS higher for COMMON EDGES + COLOUR black
edges.dynGENIE.igraph <- as_edgelist(igraph.dynGENIE)
edges.dynGENIE.igraph <- data.frame(reg=edges.dynGENIE.igraph[,1],
                                    target=edges.dynGENIE.igraph[,2],
                                    weight=rep(0.1,nrow(edges.dynGENIE.igraph)),
                                    stringsAsFactors = FALSE
)
# vertices.dynGENIE.igraph <- unique(c(edges.dynGENIE.igraph$reg, edges.dynGENIE.igraph$target))
# edges.aranet.dynGENIE <- edges.aranet[edges.aranet$regulator %in% vertices.dynGENIE.igraph & 
#                                         edges.aranet$target %in% vertices.dynGENIE.igraph,]
# common.edges.dynGENIE <- network.pairwise.edge.comparison(sample.list = edges.dynGENIE.igraph,
#                                                           reference.list = edges.aranet.dynGENIE,
#                                                           node.names = vertices.dynGENIE.igraph
# )
# E(igraph.dynGENIE)$"common.edges.width" <- ifelse(common.edges.dynGENIE$common.edges.index,1,0.1)


## SCALAR STATISTICS
#closeness.global <- igraph::centralization.closeness(igraph.dynGENIE)$centralization
global.clustering <- igraph::transitivity(as.undirected(igraph.dynGENIE), 
                                          type = "global")
shortest.path.matrix <- igraph::shortest.paths(igraph.dynGENIE, weights = NA)
diameter <- diameter(igraph.dynGENIE)
mean.shortest.path <- mean(shortest.path.matrix[shortest.path.matrix!=Inf])





## ---- plot GGNET2 (network object)
coords <- read.csv(file.path(dir,"ggnet2_graph_coords_dynGENIE3_splined_250718.csv"),header = TRUE,
                   stringsAsFactors = FALSE)
V(igraph.dynGENIE)$"x" <- coords$x; V(igraph.dynGENIE)$"y" <- coords$y
ggnet.layout.par <- list(niter=500,
                         repulse.rad=1000*log(vcount(igraph.dynGENIE))*vcount(igraph.dynGENIE)^2,
                         cool.exp=10,
                         ncell=vcount(igraph.dynGENIE)^0.5
)

set.seed(3246)
ggnet2(igraph.dynGENIE, 
       # nodes
       node.color = "node.col.membership.dir",
       node.size = "node.size.betweenDegree", 
       node.alpha = "node.alpha.degree",
       node.shape = "node.shape.betweenDegree",
       # edges
       #edge.size = "common.edges.width",
       edge.color = "edgecolor", 
       
       # layout coords
       #mode = c("x","y"), 
       # layout pars
       layout.par = ggnet.layout.par
       # neighbour
       #na.rm = "BLH3_order1"
) +
  guides(size=FALSE) + 
  ggplot2::ggtitle("col=membership.dir, size=betweenDegree, alpha=degree, shape=betweenDegree,
                   edgesize=common.edges.width, edgecol=weight")






# V(igraph.dynGENIE)$"x" <- x$data$x; V(igraph.dynGENIE)$"y" <- x$data$y
# write.csv(cbind(x=V(igraph.dynGENIE)$"x",y=V(igraph.dynGENIE)$"y"),
#           file.path(dir,"ggnet2_graph_coords_dynGENIE3_170718.csv"),
#           quote = FALSE,row.names = FALSE)




## ==== (2) Plots ---------------------------------------------------------- ##
## Degree distribution 
par(mfrow=c(1,2))
barplot(degree.distribution(igraph.dynGENIE),
        names.arg = as.character(1:length(degree.distribution(igraph.dynGENIE))), 
        xlab = "degree (k)", ylab="frequency (p_k)",
        col = c(rep("red",degree.threshold),
                    rep("blue",(length(degree.distribution(igraph.dynGENIE))+1-degree.threshold)
                        )),
        main = "degree distribution")

ckk.col <- ifelse(V(igraph.dynGENIE)$"degree" > degree.threshold, "blue","red")
## k vs C(k) - clustering coefficient
plot(V(igraph.dynGENIE)$"degree", V(igraph.dynGENIE)$"clustCoeff",
     pch=16, xlab="degree (k)", ylab="Clustering coefficient C(k)",
     main="clustering k vs C(k)", col=ckk.col)

par(mfrow=c(1,1))



par(mfrow=c(1,2))
dd.col <- c(rep("red", degree.threshold),
            rep("blue",
                length(degree.distribution(igraph.dynGENIE))-degree.threshold))
dd.col[degree.distribution(igraph.dynGENIE)==0] <- "grey20"
plot(c(0.9, 1:(length(degree.distribution(igraph.dynGENIE))-1)),
     degree.distribution(igraph.dynGENIE)+0.001, log = "xy",
     main = "degree distribution (logscale)",
     xlab = "degree (k)", ylab="frequency (p_k) + 0.001", pch=16, cex=1.5,
     col = dd.col
)
## distance histogram
hist(shortest.path.matrix, xlab="shortest path",
     main=paste0("shortest path histogram \n diameter = ",
                 round(diameter,3),", mean = ", round(mean.shortest.path,3)))
par(mfrow=c(1,1))








## ==== (3) COMPARISON ----------------------------------------------------- ##
edges.aranet.selected <- edges.aranet[edges.aranet$regulator %in% gene.list$x &
                                        edges.aranet$target %in% gene.list$x,]
dynGENIE.cutoff <- seq(0,max(res.dynGENIE[,3])*0.9, length.out = 15)
precision <- recall <- specificity <- 
  sensitivity <- rep(NA, length(dynGENIE.cutoff))


# for AraNet
progBar <- txtProgressBar(0, length(dynGENIE.cutoff), style= 3)
for (i in 1:length(dynGENIE.cutoff)){
  compare <- network.pairwise.edge.comparison(
    sample.list = res.dynGENIE[res.dynGENIE[,3]>dynGENIE.cutoff[i],], 
    reference.list = edges.aranet.selected, #, 
    node.names = gene.list$x)
  
  precision[i] <- compare$metrics["prec"]
  recall[i] <- compare$metrics["recall"]
  specificity[i] <- compare$metrics["specificity"]
  sensitivity[i] <- compare$metrics["sensitivity"]
  setTxtProgressBar(progBar, i)
}
close(progBar)

par(mfrow=c(1,2))
plot(1-specificity, sensitivity, type="b",
     main="ROC Curve \n -Pairwise AraNet-", xlim=c(0,1),ylim=c(0,1),
     col=colorRampPalette(c("red","blue"))(length(dynGENIE.cutoff)),
     pch=16)
lines(1:100/100,1:100/100)
plot(recall, precision, type="b", 
     main="Precision-Recall curve \n -Pairwise AraNet-",
     col=colorRampPalette(c("red","blue"))(length(dynGENIE.cutoff)), pch=16)
legend("topright", 
       legend = as.character(round(dynGENIE.cutoff[!is.na(sensitivity)],4)),
       fill = colorRampPalette(c("red","blue"))(length(dynGENIE.cutoff)))
par(mfrow=c(1,1))





## ==== (4) TF link HEATMAP ------------------------------------------------ ##
adj.mat.dynGENIE.TFs <- get.adjmatrix.from.list(res.dynGENIE.filt, gene.list$x)
adj.mat.dynGENIE.TFs <- adj.mat.dynGENIE.TFs[TF.list$x, TF.list$x]
connected.TFs <- rowSums(adj.mat.dynGENIE.TFs) > 0.02 | 
  colSums(adj.mat.dynGENIE.TFs) > 0.02
# heatmap of all TFs with h-clustering
TF.heatmap <- heatmap.2(adj.mat.dynGENIE.TFs,
                        trace = 'none', labRow = NA, labCol = NA,
                        col=brewer.pal(5,"OrRd"), #dendrogram = 'none',
                        main = "Heatmap of TF connectivity", key.title = "edge weight"
)


TF.heatmap.cluster <- heatmap.2(TF.heatmap$carpet[1:30,1:30],
                                trace = 'none', labRow = NA, labCol = NA, 
                                col=c("lightyellow",colorRampPalette(c("lightyellow","red3"))(50)),
                                main = "Selected TF connectivity", key.title = "edge weight"
)

TF.heatmap.cluster.names <- (unlist(list(from=rownames(TF.heatmap.cluster$carpet),
                                         to=colnames(TF.heatmap.cluster$carpet))))
TF.heatmap.cluster.names.col <- ifelse(grepl("from",names(TF.heatmap.cluster.names)),'blue','red')
names(TF.heatmap.cluster.names.col) <- TF.heatmap.cluster.names

# plot expr of the TF clusters 
plot(time,rep(NA,length(time)),ylim=c(-2.5,2.5), type='b', main="TF cluster expression",
     xlab="time (mins)", ylab="CPM Z-score")
for(i in TF.heatmap.cluster.names){
  lines(time,Zscore.normalise.expr(expr.selected[i,]),type='b',
        col=TF.heatmap.cluster.names.col[i], lwd=2, pch=16)
}
legend("bottom",legend = c("from","to"),col = c("blue","red"),
       lty = c(1,1),lwd=c(2,2), horiz = TRUE, bty = 'n')


V(igraph.dynGENIE)$"node.color.TF.heatmapcluster" <- 
  V(igraph.dynGENIE)$"node.color.TF"
V(igraph.dynGENIE)$"node.color.TF.heatmapcluster"[
  V(igraph.dynGENIE)$name %in% unlist(TF.heatmap.cluster.names)] <- 
  brew_colors(col = "yellow")





## ==== (5) COMMUNITY analysis --------------------------------------------- ##
pdf(file.path(dir,"ExprPlots_communities_Col_22_splined.pdf"),width = 15, height = 10)
# plot for communities 
communities.sorted.population <- sort(table(V(igraph.dynGENIE)$"membership.dir"), decreasing = TRUE)
for (i in 1:length(communities.sorted.population)){
  community.chosen <- as.numeric(names(communities.sorted.population))[i]
  communities.sorted.population.nodenames <- V(igraph.dynGENIE)$name[V(igraph.dynGENIE)$"membership.dir"==community.chosen]
  
  plot(time,rep(NA,length(time)),ylim=c(-2.5,2.5), type='b', 
       main=paste0("expression for community ", community.chosen),
       xlab="time (mins)", ylab="CPM Z-score")
  for(j in communities.sorted.population.nodenames){
    lines(time,Zscore.normalise.expr(expr.selected[j,]),type='b',
          #col=TF.heatmap.cluster.names.col[i], 
          #col=exprtime.communities.col[i],
          lwd=2, pch=16)
  }
}
dev.off()






## ==== (6) FUNCTIONAL CONNECTIVITY HEATMAP -------------------------------- ##
## -- ALL genes --
GO.terms.interest <- c("circadian","light", "stress", "heat", "cold", "water",
                          "auxin", "abscisic", "jasmonic","ethylene",
                          "flower", "stimulus", "metabolic","signaling",
                          "biosynthetic","phosphorylation", "growth",
                          "temperature", "photosynthesis", "kinase",
                          "transferase", "ion binding",
                       "DNA binding", "other")  
gene.list.GO.modified <- matrix(NA, length(gene.list$x), 2,
                                dimnames = list(NULL, c("gene","GO")))
gene.list.GO.modified[,"gene"] <- gene.list$x
gene.list.GO.modified[,"GO"] <- apply(gene.list.GO[gene.list.GO$V1 %in% 
                                                        gene.list.GO.modified[,1], 3:7], 
                                         1, 
                                         function(x){
                                           x <- paste(x,collapse = "_")
                                           x <- unique(unlist(strsplit(x,split = "_")))
                                           paste(x,collapse = "_")
                                         })

GO.adj.matrix <- matrix(0, length(GO.terms.interest), length(GO.terms.interest),
                        dimnames = list(GO.terms.interest,GO.terms.interest))
edge.list.GO.link <- igraph::get.edgelist(igraph.dynGENIE)
gene.list.GO.modified[gene.list.GO.modified[,2]=="",2] <- sub("","other",gene.list.GO.modified[gene.list.GO.modified[,2]=="",2])

# raw
for (i in 1:nrow(edge.list.GO.link)){
  # find genes 
  regulator.index <- which(gene.list.GO.modified[,1]==edge.list.GO.link[i,1])
  target.index <- which(gene.list.GO.modified[,1]==edge.list.GO.link[i,2])
  # extract functions of each gene 
  reg.fun <- as.vector(unlist(strsplit(gene.list.GO.modified[regulator.index,2],
                                       split = "_")))
  target.fun <- as.vector(unlist(strsplit(gene.list.GO.modified[target.index,2],
                                          split = "_")))
  # see whether our GO terms of interest are in each gene 
  reg.fun <- (sapply(GO.terms.interest, function(x){grepl(x,reg.fun)}))
  target.fun <- (sapply(GO.terms.interest, function(x){grepl(x,target.fun)}))
  
  if(is.matrix(reg.fun)){
    reg.fun <- apply(reg.fun,2,sum)
  } else {
    reg.fun <- as.numeric(reg.fun)
  }
  if(is.matrix(target.fun)){
    target.fun <- apply(target.fun,2,sum)
  } else {
    target.fun <- as.numeric(target.fun)
  }
  
  # exdtract the GO terms of interest for each gene 
  if (sum(reg.fun)==0){  # regulator
    reg.fun <- "other"
    # reg.fun <- as.vector(unlist(strsplit(gene.list.GO.modified[regulator.index,2],
    #                                      split = "_")))
    # if (sum(grepl("DNA binding", reg.fun))>0){
    #   reg.fun <- "DNA binding"
    # } else {
    #   reg.fun <- "other"
    # }
  } else {
    reg.fun <- names(reg.fun[reg.fun>0])
  }
  if (sum(target.fun)==0){  # target
    target.fun <- "other"
    # target.fun <- as.vector(unlist(strsplit(gene.list.GO.modified.TF[target.index,2],
    #                                         split = "_")))
    # if (sum(grepl("DNA binding", target.fun))>0){
    #   target.fun <- "DNA binding"
    # } else {
    #   target.fun <- "other"
    # }
  } else {
    target.fun <- names(target.fun[target.fun>0])
  }
  
  # fill the connectivity matrix 
  for (j in reg.fun){
    for (k in target.fun){
      GO.adj.matrix[j,k] <- GO.adj.matrix[j,k] + 1
    }
  }
}

# random
gene.list.GO.modified.random <- gene.list.GO.modified
GO.adj.matrix.random <- GO.adj.matrix
GO.adj.matrix.random <- GO.adj.matrix.random-GO.adj.matrix.random
GO.adj.matrix.random.total <- GO.adj.matrix.random
GO.adj.matrix.pvalue <- GO.adj.matrix.random
n.shuffles <- 200
progBar <- txtProgressBar(0,n.shuffles,style= 3)
set.seed(123)
for (m in 1:n.shuffles){
  gene.list.GO.modified.random[,"gene"] <- sample(gene.list.GO.modified[,"gene"])
  GO.adj.matrix.random <- GO.adj.matrix.random-GO.adj.matrix.random
  for (i in 1:nrow(edge.list.GO.link)){
    # find genes 
    regulator.index <- which(gene.list.GO.modified.random[,1]==edge.list.GO.link[i,1])
    target.index <- which(gene.list.GO.modified.random[,1]==edge.list.GO.link[i,2])
    # extract functions of each gene 
    reg.fun <- as.vector(unlist(strsplit(gene.list.GO.modified.random[regulator.index,2],
                                         split = "_")))
    target.fun <- as.vector(unlist(strsplit(gene.list.GO.modified.random[target.index,2],
                                            split = "_")))
    # see whether our GO terms of interest are in each gene 
    reg.fun <- (sapply(GO.terms.interest, function(x){grepl(x,reg.fun)}))
    target.fun <- (sapply(GO.terms.interest, function(x){grepl(x,target.fun)}))
    
    if(is.matrix(reg.fun)){
      reg.fun <- apply(reg.fun,2,sum)
    } else {
      reg.fun <- as.numeric(reg.fun)
    }
    if(is.matrix(target.fun)){
      target.fun <- apply(target.fun,2,sum)
    } else {
      target.fun <- as.numeric(target.fun)
    }
    
    # exdtract the GO terms of interest for each gene 
    if (sum(reg.fun)==0){  # regulator
      reg.fun <- "other"
      # reg.fun <- as.vector(unlist(strsplit(gene.list.GO.modified[regulator.index,2],
      #                                      split = "_")))
      # if (sum(grepl("DNA binding", reg.fun))>0){
      #   reg.fun <- "DNA binding"
      # } else {
      #   reg.fun <- "other"
      # }
    } else {
      reg.fun <- names(reg.fun[reg.fun>0])
    }
    if (sum(target.fun)==0){  # target
      target.fun <- "other"
      # target.fun <- as.vector(unlist(strsplit(gene.list.GO.modified.TF[target.index,2],
      #                                         split = "_")))
      # if (sum(grepl("DNA binding", target.fun))>0){
      #   target.fun <- "DNA binding"
      # } else {
      #   target.fun <- "other"
      # }
    } else {
      target.fun <- names(target.fun[target.fun>0])
    }
    
    # fill the connectivity matrix 
    for (j in reg.fun){
      for (k in target.fun){
        GO.adj.matrix.random[j,k] <- GO.adj.matrix.random[j,k] + 1
      }
    }
  }
  
  GO.adj.matrix.random.total <- GO.adj.matrix.random.total + GO.adj.matrix.random
  GO.adj.matrix.pvalue[GO.adj.matrix>GO.adj.matrix.random] <- GO.adj.matrix.pvalue[GO.adj.matrix>GO.adj.matrix.random]+1
  
  setTxtProgressBar(progBar, m)
  
}
close(progBar)

GO.adj.matrix.norm <- log(GO.adj.matrix+1) - log((GO.adj.matrix.random.total/n.shuffles)+1)
GO.adj.matrix.pvalue <- GO.adj.matrix.pvalue/n.shuffles





## -- TFs only ----------------------------------------------------------------
# pick TFs (DNA binding) only 
gene.list.GO.modified.TF <- gene.list.GO.modified[gene.list.GO.modified[,1] %in% TF.list$x,]
# pick edges containing DNA-binding genes only 
edge.list.GO.link.TF <- edge.list.GO.link[edge.list.GO.link[,1] %in% 
                                            gene.list.GO.modified.TF[,1] & 
                                            edge.list.GO.link[,2] %in% 
                                            gene.list.GO.modified.TF[,1],
                                          ]
GO.adj.matrix.TF <- GO.adj.matrix.norm
GO.adj.matrix.TF <- GO.adj.matrix.TF-GO.adj.matrix.TF

# fill the adjacency matrix for the TF roles connectivity (REAL NETWORK)                                           
for (i in 1:nrow(edge.list.GO.link.TF)){
  # find genes 
  regulator.index <- which(gene.list.GO.modified.TF[,1]==edge.list.GO.link.TF[i,1])
  target.index <- which(gene.list.GO.modified.TF[,1]==edge.list.GO.link.TF[i,2])
  # extract functions of each gene 
  reg.fun <- as.vector(unlist(strsplit(gene.list.GO.modified.TF[regulator.index,2],
                                       split = "_")))
  target.fun <- as.vector(unlist(strsplit(gene.list.GO.modified.TF[target.index,2],
                                          split = "_")))
  # see whether our GO terms of interest are in each gene 
  reg.fun <- (sapply(GO.terms.interest, function(x){grepl(x,reg.fun)}))
  target.fun <- (sapply(GO.terms.interest, function(x){grepl(x,target.fun)}))
  
  if(is.matrix(reg.fun)){
    reg.fun <- apply(reg.fun,2,sum)
  } else {
    reg.fun <- as.numeric(reg.fun)
  }
  if(is.matrix(target.fun)){
    target.fun <- apply(target.fun,2,sum)
  } else {
    target.fun <- as.numeric(target.fun)
  }
  
  # exdtract the GO terms of interest for each gene 
  if (sum(reg.fun)==0){  # regulator
    reg.fun <- "other"
    # reg.fun <- as.vector(unlist(strsplit(gene.list.GO.modified[regulator.index,2],
    #                                      split = "_")))
    # if (sum(grepl("DNA binding", reg.fun))>0){
    #   reg.fun <- "DNA binding"
    # } else {
    #   reg.fun <- "other"
    # }
  } else {
    reg.fun <- names(reg.fun[reg.fun>0])
  }
  if (sum(target.fun)==0){  # target
    target.fun <- "other"
    # target.fun <- as.vector(unlist(strsplit(gene.list.GO.modified.TF[target.index,2],
    #                                         split = "_")))
    # if (sum(grepl("DNA binding", target.fun))>0){
    #   target.fun <- "DNA binding"
    # } else {
    #   target.fun <- "other"
    # }
  } else {
    target.fun <- names(target.fun[target.fun>0])
  }
  
  # fill the connectivity matrix 
  for (j in reg.fun){
    for (k in target.fun){
      GO.adj.matrix.TF[j,k] <- GO.adj.matrix.TF[j,k] + 1
    }
  }
}


# fill the adjacency matrix for the TF roles connectivity (RANDOM NETWORK)  
gene.list.GO.modified.TF.random <- gene.list.GO.modified.TF
GO.adj.matrix.TF.random <- GO.adj.matrix.TF
GO.adj.matrix.TF.random <- GO.adj.matrix.TF.random-GO.adj.matrix.TF.random
GO.adj.matrix.TF.random.total <- GO.adj.matrix.TF.random
GO.adj.matrix.TF.pvalue <- GO.adj.matrix.TF.random

progBar <- txtProgressBar(0,n.shuffles,style= 3)
set.seed(123)
for(m in 1:n.shuffles){
  gene.list.GO.modified.TF.random[,"gene"]<-sample(gene.list.GO.modified.TF[,"gene"])
  GO.adj.matrix.TF.random <- GO.adj.matrix.TF.random-GO.adj.matrix.TF.random
  for (i in 1:nrow(edge.list.GO.link.TF)){
    # find genes 
    regulator.index <- which(gene.list.GO.modified.TF.random[,1]==edge.list.GO.link.TF[i,1])
    target.index <- which(gene.list.GO.modified.TF.random[,1]==edge.list.GO.link.TF[i,2])
    # extract functions of each gene 
    reg.fun <- as.vector(unlist(strsplit(gene.list.GO.modified.TF.random[regulator.index,2],
                                         split = "_")))
    target.fun <- as.vector(unlist(strsplit(gene.list.GO.modified.TF.random[target.index,2],
                                            split = "_")))
    # see whether our GO terms of interest are in each gene 
    reg.fun <- (sapply(GO.terms.interest, function(x){grepl(x,reg.fun)}))
    target.fun <- (sapply(GO.terms.interest, function(x){grepl(x,target.fun)}))
    
    if(is.matrix(reg.fun)){
      reg.fun <- apply(reg.fun,2,sum)
    } else {
      reg.fun <- as.numeric(reg.fun)
    }
    if(is.matrix(target.fun)){
      target.fun <- apply(target.fun,2,sum)
    } else {
      target.fun <- as.numeric(target.fun)
    }
    
    # exdtract the GO terms of interest for each gene 
    if (sum(reg.fun)==0){  # regulator
      reg.fun <- "other"
      # reg.fun <- as.vector(unlist(strsplit(gene.list.GO.modified[regulator.index,2],
      #                                      split = "_")))
      # if (sum(grepl("DNA binding", reg.fun))>0){
      #   reg.fun <- "DNA binding"
      # } else {
      #   reg.fun <- "other"
      # }
    } else {
      reg.fun <- names(reg.fun[reg.fun>0])
    }
    if (sum(target.fun)==0){  # target
      target.fun <- "other"
      # target.fun <- as.vector(unlist(strsplit(gene.list.GO.modified.TF[target.index,2],
      #                                         split = "_")))
      # if (sum(grepl("DNA binding", target.fun))>0){
      #   target.fun <- "DNA binding"
      # } else {
      #   target.fun <- "other"
      # }
    } else {
      target.fun <- names(target.fun[target.fun>0])
    }
    
    # fill the connectivity matrix 
    for (j in reg.fun){
      for (k in target.fun){
        GO.adj.matrix.TF.random[j,k] <- GO.adj.matrix.TF.random[j,k] + 1
      }
    }
  }
  
  GO.adj.matrix.TF.random.total <- GO.adj.matrix.TF.random.total+GO.adj.matrix.TF.random
  GO.adj.matrix.TF.pvalue[GO.adj.matrix.TF>GO.adj.matrix.TF.random] <- GO.adj.matrix.TF.pvalue[GO.adj.matrix.TF>GO.adj.matrix.TF.random] + 1
  
  setTxtProgressBar(progBar, m)
}
close(progBar)

GO.adj.matrix.TF.norm <- log(GO.adj.matrix.TF+1) - log((GO.adj.matrix.TF.random.total/n.shuffles)+1)
GO.adj.matrix.TF.pvalue <- GO.adj.matrix.TF.pvalue/n.shuffles




#pdf(file.path(dir,"network_prelim_GO_connectivity.pdf"))
heatmap.2(GO.adj.matrix.pvalue, trace = 'none', key.title = "p-value" ,
          col=c(colorRampPalette(c("blue3","lightblue"))(100),
                colorRampPalette(c("white"))(40),
                colorRampPalette(c("pink","red3"))(100)
          ),
          #col=c("lightyellow",colorRampPalette(c("yellow","orange","red","red2",
          #                                       "red3"))(200),"red4"),
          main = "Enriched functional connectivity \n ALL genes (based on GO terms)")
heatmap.2(GO.adj.matrix.norm, trace = 'none', key.title = "log norm counts" ,
          col=c(colorRampPalette(c("blue4","blue","lightblue"))(100),
                colorRampPalette(c("white"))(10),
                colorRampPalette(c("pink","red","red4"))(100)),
          #col=c("lightyellow",colorRampPalette(c("yellow","orange","red","red2",
          #                                       "red3"))(200),"red4"),
          main = "Normalised functional connectivity \n ALL genes (based on GO terms)")
heatmap.2(GO.adj.matrix.TF.pvalue,
          trace = 'none', key.title = "p-value" , 
          col=c(colorRampPalette(c("blue3","lightblue"))(100),
                colorRampPalette(c("white"))(40),
                colorRampPalette(c("pink","red3"))(100)),
          #col=c("lightyellow",colorRampPalette(c("yellow","orange","red","red2",
          #                                       "red3"))(200),"red4"), 
          main = "Normalised functional connectivity \n TFs (based on GO terms)")
heatmap.2(GO.adj.matrix.TF.norm,
          trace = 'none', key.title = "log norm counts" , 
          col=c(colorRampPalette(c("blue4","blue","lightblue"))(100),
                colorRampPalette(c("white"))(10),
                colorRampPalette(c("pink","red","red4"))(100)),
          #col=c("lightyellow",colorRampPalette(c("yellow","orange","red","red2",
          #                                       "red3"))(200),"red4"), 
          main = "Enriched functional connectivity \n TFs (based on GO terms)")


# graphs of GO connectivities via adj matrices 
GO.adj.matrix.norm.network <- GO.adj.matrix.pvalue
cutoff.GO.allgenes <- 0.8
GO.adj.matrix.norm.network[GO.adj.matrix.pvalue >= cutoff.GO.allgenes] <- 1
GO.adj.matrix.norm.network[!(GO.adj.matrix.pvalue >= cutoff.GO.allgenes)] <- 0

igraph.GO.allgenes <- igraph::graph.adjacency(GO.adj.matrix.norm.network, mode = "directed")
igraph.GO.allgenes <- igraph::simplify(igraph.GO.allgenes)
E(igraph.GO.allgenes)$"col.circadian.path" <- ifelse(get.edgelist(igraph.GO.allgenes)[,1]=="circadian","blue","grey")
#E(igraph.GO.allgenes)$"col.circadian.path"[get.edgelist(igraph.GO.allgenes)[,1]=="abscisic"] <- "pink"
#V(igraph.GO.allgenes)$"node.col.degree" <- igraph::degree(igraph.GO.allgenes, mode = "in")

set.seed(1238)
ggnet2(igraph.GO.allgenes, label = TRUE,arrow.size = 7, arrow.gap = 0.025, 
       edge.color = "col.circadian.path") + ggtitle("all genes")
  

# graphs of GO connectivities via adj matrices - TFs 
GO.adj.matrix.TF.norm.network <- GO.adj.matrix.TF.pvalue
cutoff.GO.TF <- 0.8
GO.adj.matrix.TF.norm.network[GO.adj.matrix.TF.pvalue >= cutoff.GO.TF] <- 1
GO.adj.matrix.TF.norm.network[!(GO.adj.matrix.TF.pvalue >= cutoff.GO.TF)] <- 0

igraph.GO.TF <- igraph::graph.adjacency(GO.adj.matrix.TF.norm.network, mode = "directed")
igraph.GO.TF <- igraph::simplify(igraph.GO.TF)
E(igraph.GO.TF)$"col.circadian.path" <- ifelse(get.edgelist(igraph.GO.TF)[,1]=="circadian","blue","grey")
#E(igraph.GO.allgenes)$"col.circadian.path"[get.edgelist(igraph.GO.allgenes)[,1]=="abscisic"] <- "pink"
#V(igraph.GO.allgenes)$"node.col.degree" <- igraph::degree(igraph.GO.allgenes, mode = "in")

set.seed(1237)
ggnet2(igraph.GO.TF, label = TRUE,arrow.size = 7, arrow.gap = 0.025, 
       edge.color = "col.circadian.path", mode="kamadakawai") + ggtitle("TF")
#dev.off()





## ==== (7) EXPRESSION TIME HEATMAP ---------------------------------------- ##
# extract edge list
edge.list.peakExprTime <- igraph::get.edgelist(igraph.dynGENIE)
edge.list.peakExprTime.times <- edge.list.peakExprTime
# assign peak expression time for each gene and each edge (time list)
edge.list.peakExprTime.times[,1] <- (peak.expr.time$peak.time[edge.list.peakExprTime.times[,1]])
edge.list.peakExprTime.times[,2] <- (peak.expr.time$peak.time[edge.list.peakExprTime.times[,2]])
edge.list.peakExprTime.times <- apply(edge.list.peakExprTime.times,2,as.numeric)
peakExpr.adj.matrix <- matrix(0,length(time),length(time),
                              dimnames=list(as.character(time), 
                                            as.character(time)))
# find the time interval of max expression for each gene 
# find # links that traverse between each time point 
for (i in 1:nrow(edge.list.peakExprTime)){
  regulator <- edge.list.peakExprTime[i,1]
  target <- edge.list.peakExprTime[i,2]
  tp.reg <- peak.expr.time$peak.time[regulator]
  tp.target <- peak.expr.time$peak.time[target]
  peakExpr.adj.matrix[tp.reg,tp.target] <- peakExpr.adj.matrix[tp.reg,tp.target] + 1
}

# random connectivity
peak.expr.time.shuffled <- peak.expr.time$peak.time
names(peak.expr.time.shuffled) <- names(peak.expr.time$peak.time)
edge.list.peakExprTime.times.random <- edge.list.peakExprTime.times
edge.list.peakExprTime.times.random[,1] <- (peak.expr.time.shuffled[edge.list.peakExprTime[,1]])
edge.list.peakExprTime.times.random[,2] <- (peak.expr.time.shuffled[edge.list.peakExprTime[,2]])
edge.list.peakExprTime.times.random <- apply(edge.list.peakExprTime.times.random,2,as.numeric)
peakExpr.adj.matrix.random <- matrix(0,length(time),length(time),
                                     dimnames=list(as.character(time), 
                                                   as.character(time)))
# find the time interval of max expression for each gene 
# find # links that traverse between each time point 
for (i in 1:nrow(edge.list.peakExprTime)){
  regulator <- edge.list.peakExprTime[i,1]
  target <- edge.list.peakExprTime[i,2]
  tp.reg <- peak.expr.time.shuffled[regulator]
  tp.target <- peak.expr.time.shuffled[target]
  peakExpr.adj.matrix.random[tp.reg,tp.target] <- peakExpr.adj.matrix.random[tp.reg,tp.target] + 1
}


peakExpr.adj.matrix.norm <- (peakExpr.adj.matrix - peakExpr.adj.matrix.random)/peakExpr.adj.matrix.random


# determining repressive or activating effect of TFs 
# time.delay.connectivity <- apply(edge.list.peakExprTime.times,1,diff)
# edge.list.peakExprTime.noTD <- edge.list.peakExprTime[time.delay.connectivity==0,]
# colnames(edge.list.peakExprTime.noTD) <- c("from.max","to.min")
# edge.list.peakExprTime.noTD.times <- edge.list.peakExprTime.noTD
# edge.list.peakExprTime.noTD.times[,1] <- peak.expr.time$peak.time[edge.list.peakExprTime.noTD.times[,1]]
# edge.list.peakExprTime.noTD.times[,2] <- temporal.order.genes(expr.selected[edge.list.peakExprTime.noTD.times[,2],],sorting.order = "min")$peak.time
# edge.list.peakExprTime.noTD.times <- apply(edge.list.peakExprTime.noTD.times,2,as.numeric)
# hist(time.delay.connectivity, breaks=100, 
#      main="time delay \n regulator -> target",
#      xlab="time steps")



# heatmap
heatmap.2(peakExpr.adj.matrix, Rowv = NA, Colv = NA, dendrogram = 'none', 
          trace = 'none', 
          col=c("lightyellow",colorRampPalette(c("yellow","orange","red1",
                                                 "red2","red3","red4"))(100)), 
          main="Peak expression time connectivity")
heatmap.2(peakExpr.adj.matrix, 
          trace = 'none', #col=brewer.pal(9,"Blues"), 
          col=c("lightyellow",colorRampPalette(c("yellow","orange","red1",
                                                 "red2","red3","red4"))(100)), 
          main="Peak expression time connectivity (clustered)")
heatmap.2((abs(peakExpr.adj.matrix.norm)), 
          Rowv = NA, Colv = NA, dendrogram = 'none', 
          trace = 'none', 
          col=c("lightyellow",colorRampPalette(c("yellow","orange","red1",
                                                 "red2","red3","red4"))(100)), 
          main="Normalised peak expression time connectivity")
heatmap.2((abs(peakExpr.adj.matrix.norm)), 
          #Rowv = NA, Colv = NA, dendrogram = 'none', 
          trace = 'none', 
          col=c("lightyellow",colorRampPalette(c("yellow","orange","red1",
                                                 "red2","red3","red4"))(100)), 
          main="Normalised peak expression time connectivity \n (clustered)")



weird.edges.times <- edge.list.peakExprTime[which(edge.list.peakExprTime.times[,1] > 
                                                    edge.list.peakExprTime.times[,2]),]


# pdf(file.path(dir,"regulator_expressed_later_exprPlots_splined_col22.pdf"))
# for (i in 1:nrow(weird.edges.times)){
# y.range <- range(c(expr.selected[weird.edges.times[i,1],],
#                   expr.selected[weird.edges.times[i,2],]))
# y.range[1] <- y.range[1]*0.9; y.range[2] <- y.range[2]*1.1
# plot(time,expr.selected[weird.edges.times[i,1],],
#      ylim=y.range, col="blue",type="b", lwd=2,
#      main=paste0(weird.edges.times[i,1],"->",weird.edges.times[i,2]))
# lines(time,expr.selected[weird.edges.times[i,2],],col="red",
#       type="b",lty=2,lwd=2)
# legend("topright",c("from","to"),col=c("blue","red"),lty = c(1,2),
#        lwd = c(2,2))}
# dev.off()






## ==== (7) TEMPORAL NETWORK VIZ ------------------------------------------- ## 

## --- NEIGHBOURHOODS -----------------------------
# # for BLH3
# test.neighbor <- V(igraph.dynGENIE)$name[
#     as.vector(neighborhood(igraph.dynGENIE,order = 4, 
#                            nodes = "AT1G75410",mode = "out")[[1]])]
#   V(igraph.dynGENIE)$"BLH3_order4" <- ifelse(V(igraph.dynGENIE)$name %in% 
#                                                         test.neighbor,1,NA)





## grid - timelapse 
# t1 <- ggnet2(igraph.dynGENIE, 
#        # nodes
#        node.color = node.color,
#        node.size = node.size,
#        node.alpha = node.alpha,
#        node.shape = node.shape,
#        # edges
#        edge.size = edge.size,
#        edge.color = edge.color, 
#        
#        # layout coords
#        mode = c("x","y"), 
#        # layout pars
#        layout.par = ggnet.layout.par,
#        # neighbour
#        na.rm = "BLH3_order1"
# )
# t2 <- ggnet2(igraph.dynGENIE, 
#              # nodes
#              node.color = node.color,
#              node.size = node.size,
#              node.alpha = node.alpha,
#              node.shape = node.shape,
#              # edges
#              edge.size = edge.size,
#              edge.color = edge.color, 
#              
#              # layout coords
#              mode = c("x","y"), 
#              # layout pars
#              layout.par = ggnet.layout.par,
#              # neighbour
#              na.rm = "BLH3_order2"
# )
# t3 <- ggnet2(igraph.dynGENIE, 
#              # nodes
#              node.color = node.color,
#              node.size = node.size,
#              node.alpha = node.alpha,
#              node.shape = node.shape,
#              # edges
#              edge.size = edge.size,
#              edge.color = edge.color, 
#              
#              # layout coords
#              mode = c("x","y"), 
#              # layout pars
#              layout.par = ggnet.layout.par,
#              # neighbour
#              na.rm = "BLH3_order3"
# )
# t4 <- ggnet2(igraph.dynGENIE, 
#              # nodes
#              node.color = node.color,
#              node.size = node.size,
#              node.alpha = node.alpha,
#              node.shape = node.shape,
#              # edges
#              edge.size = edge.size,
#              edge.color = edge.color, 
#              
#              # layout coords
#              mode = c("x","y"), 
#              # layout pars
#              layout.par = ggnet.layout.par,
#              # neighbour
#              na.rm = "BLH3_order4"
# )
# 
# g <- guides(size=FALSE,alpha=FALSE,color=FALSE) 
# y.coord <- scale_y_continuous(limits = range(V(igraph.dynGENIE)$y * 1.1), 
#                               breaks = NULL)
# x.coord <- scale_x_continuous(limits = range(V(igraph.dynGENIE)$x * 1.1), 
#                               breaks = NULL)
# gridExtra::grid.arrange(t1 + x.coord + y.coord + g + 
#                           ggtitle(paste0("col=",node.color, ", size=", node.size, 
#                                                                     ", \n alpha=", node.alpha, ", shape=", node.shape,
#                                                                     ", \n edgesize=", edge.size, ", \n edgecol=",edge.color, 
#                                                                     ", na.rm=",na.rm)
# ) ,
#                         t2 + x.coord + y.coord + g + ggtitle("order 2"),
#                         t3 + x.coord + y.coord + g + ggtitle("order 3"),
#                         t4 + x.coord + y.coord + g + ggtitle("order 4"),
#                         nrow=2,ncol=2)











## == PLOT EVERYTHING ====================================================== ##

pdf(file.path(dir,paste0("network_summary_v3_",mut,"_",temp,"_splined.pdf")), 
    width = 15, height = 10)
## -- (1) TDA scores -----------
par(mfrow=c(2,2))
barplot(log10(hist(res.dynGENIE$weight, plot = FALSE)$counts+1),
        main="dynGENIE weight histogram (log)",
        xlab="weight", ylab="log(frequency)", las=2,
        names.arg = hist(res.dynGENIE$weight, plot = FALSE)$breaks[-1]
)
hist(res.genereg$coef,  breaks = 30, main="GeneReg coef histogram",
     xlab="coefficient")
hist(res.genereg$delay,  breaks = 30, main="GeneReg delay histogram",
     xlab="delay (mins)")
hist(res.genereg$adj.r.squared,  breaks = 150,
     main="GeneReg adjusted R^2 histogram", xlab="adjusted R squared")
par(mfrow=c(1,1))



## -- (2) DEGREE DISTRIBUTION & CLUSTERING --------------------
## Degree distribution 
par(mfrow=c(1,2))
barplot(degree.distribution(igraph.dynGENIE),
        names.arg = as.character(1:length(degree.distribution(igraph.dynGENIE))), 
        xlab = "degree (k)", ylab="frequency (p_k)",
        col = c(rep("red",degree.threshold),
                rep("blue",(length(degree.distribution(igraph.dynGENIE))+1-degree.threshold)
                )),
        main = "degree distribution")

ckk.col <- ifelse(V(igraph.dynGENIE)$"degree" > degree.threshold, "blue","red")
## k vs C(k) - clustering coefficient
plot(V(igraph.dynGENIE)$"degree", V(igraph.dynGENIE)$"clustCoeff",
     pch=16, xlab="degree (k)", ylab="Clustering coefficient C(k)",
     main="clustering k vs C(k)", col=ckk.col)

par(mfrow=c(1,1))



par(mfrow=c(1,2))
dd.col <- c(rep("red", degree.threshold),
            rep("blue",
                length(degree.distribution(igraph.dynGENIE))-degree.threshold))
dd.col[degree.distribution(igraph.dynGENIE)==0] <- "grey20"
plot(c(0.9, 1:(length(degree.distribution(igraph.dynGENIE))-1)),
     degree.distribution(igraph.dynGENIE)+0.001, log = "xy",
     main = "degree distribution (logscale)",
     xlab = "degree (k)", ylab="frequency (p_k) + 0.001", pch=16, cex=1.5,
     col = dd.col
)
## distance histogram
hist(shortest.path.matrix, xlab="shortest path",
     main=paste0("shortest path histogram \n diameter = ",
                 round(diameter,3),", mean = ", round(mean.shortest.path,3)))
par(mfrow=c(1,1))






## -- (3) COMPARISON: Pre-Recall + AROC: GeneReg vs dynGENIE -------------------
par(mfrow=c(1,2))
plot(1-specificity, sensitivity, type="b",
     main="ROC Curve \n -Pairwise AraNet-", xlim=c(0,1),ylim=c(0,1),
     col=colorRampPalette(c("red","blue"))(length(dynGENIE.cutoff)),
     pch=16)
lines(1:100/100,1:100/100)
plot(recall, precision, type="b", 
     main="Precision-Recall curve \n -Pairwise AraNet-",
     col=colorRampPalette(c("red","blue"))(length(dynGENIE.cutoff)), pch=16)
legend("topright", 
       legend = as.character(round(dynGENIE.cutoff[!is.na(sensitivity)],4)),
       fill = colorRampPalette(c("red","blue"))(length(dynGENIE.cutoff)))
par(mfrow=c(1,1))





## -- (4) TF LINK HEATMAP -----------------------------------------------
heatmap.2(adj.mat.dynGENIE.TFs,
          trace = 'none', labRow = NA, labCol = NA,
          col=brewer.pal(5,"OrRd"), #dendrogram = 'none',
          main = "Heatmap of TF connectivity", key.title = "edge weight"
)

heatmap.2(TF.heatmap$carpet[1:20,1:20],
          trace = 'none', labRow = NA, labCol = NA, 
          col=c("lightyellow",
                colorRampPalette(c("lightyellow",
                                   "red3"))(50)),
          main = "Selected TF connectivity", 
          key.title = "edge weight"
)

plot(time,rep(NA,length(time)),ylim=c(-2.5,2.5), type='b', main="TF cluster expression",
     xlab="time (mins)", ylab="CPM Z-score")
for(i in TF.heatmap.cluster.names){
  lines(time,Zscore.normalise.expr(expr.selected[i,]),type='b',
        col=TF.heatmap.cluster.names.col[i], lwd=2, pch=16)
}
legend("bottom",legend = c("from","to"),col = c("blue","red"),
       lty = c(1,1),lwd=c(2,2), horiz = TRUE, bty = 'n')





## -- (5) FUNCTIONAL CONNECTIVITY HEATMAP -------------------------------
heatmap.2(GO.adj.matrix.pvalue, trace = 'none', key.title = "p-value" ,
          col=c(colorRampPalette(c("blue3","lightblue"))(100),
                colorRampPalette(c("white"))(40),
                colorRampPalette(c("pink","red3"))(100)
          ),
          #col=c("lightyellow",colorRampPalette(c("yellow","orange","red","red2",
          #                                       "red3"))(200),"red4"),
          main = "Enriched functional connectivity \n ALL genes (based on GO terms)")
heatmap.2(GO.adj.matrix.norm, trace = 'none', key.title = "log norm counts" ,
          col=c(colorRampPalette(c("blue4","blue","lightblue"))(100),
                colorRampPalette(c("white"))(10),
                colorRampPalette(c("pink","red","red4"))(100)),
          #col=c("lightyellow",colorRampPalette(c("yellow","orange","red","red2",
          #                                       "red3"))(200),"red4"),
          main = "Normalised functional connectivity \n ALL genes (based on GO terms)")
heatmap.2(GO.adj.matrix.TF.pvalue,
          trace = 'none', key.title = "p-value" , 
          col=c(colorRampPalette(c("blue3","lightblue"))(100),
                colorRampPalette(c("white"))(40),
                colorRampPalette(c("pink","red3"))(100)),
          #col=c("lightyellow",colorRampPalette(c("yellow","orange","red","red2",
          #                                       "red3"))(200),"red4"), 
          main = "Normalised functional connectivity \n TFs (based on GO terms)")
heatmap.2(GO.adj.matrix.TF.norm,
          trace = 'none', key.title = "log norm counts" , 
          col=c(colorRampPalette(c("blue4","blue","lightblue"))(100),
                colorRampPalette(c("white"))(10),
                colorRampPalette(c("pink","red","red4"))(100)),
          #col=c("lightyellow",colorRampPalette(c("yellow","orange","red","red2",
          #                                       "red3"))(200),"red4"), 
          main = "Enriched functional connectivity \n TFs (based on GO terms)")






## -- (6) EXPRESSION TIME HEATMAP ---------------------------------------
heatmap.2(peakExpr.adj.matrix, Rowv = NA, Colv = NA, dendrogram = 'none', 
          trace = 'none', key.title = "counts",
          col=c("lightyellow",colorRampPalette(c("yellow","orange","red1",
                                                 "red2","red3","red4"))(100)), 
          main="Peak expression time connectivity")
heatmap.2(peakExpr.adj.matrix, key.title = "counts",
          trace = 'none', #col=brewer.pal(9,"Blues"), 
          col=c("lightyellow",colorRampPalette(c("yellow","orange","red1",
                                                 "red2","red3","red4"))(100)), 
          main="Peak expression time connectivity (clustered)")
heatmap.2((abs(peakExpr.adj.matrix.norm)), 
          Rowv = NA, Colv = NA, dendrogram = 'none', 
          trace = 'none', key.title = "%change from random",
          col=c("lightyellow",colorRampPalette(c("yellow","orange","red1",
                                                 "red2","red3","red4"))(100)), 
          main="Normalised peak expression time connectivity")
heatmap.2((abs(peakExpr.adj.matrix.norm)), 
          #Rowv = NA, Colv = NA, dendrogram = 'none', 
          trace = 'none', key.title = "%change from random",
          col=c("lightyellow",colorRampPalette(c("yellow","orange","red1",
                                                 "red2","red3","red4"))(100)), 
          main="Normalised peak expression time connectivity \n (clustered)")


## -- (7) Networks ----------------------------

ggnet2(igraph.dynGENIE, 
       # nodes
       node.color = "node.color.exptime",
       node.size = "node.size.betweenDegree",
       node.alpha = "node.alpha.degree",
       node.shape = "node.shape.TF",
       # edges
       edge.size = "common.edges.width",
       edge.color = "edgecolor", 
       
       # layout coords
       mode = c("x","y"), 
       # layout pars
       layout.par = ggnet.layout.par
       # neighbour
       #na.rm = na.rm
) +
  guides(size=FALSE,alpha=FALSE)  + 
  ggplot2::ggtitle("node.color.exptime, size=betweenDegree, node.alpha.degree, node.shape.TF,
                   edgesize=common.edges.width, edgecol=weight")

ggnet2(igraph.dynGENIE, 
       # nodes
       node.color = "node.color.exptime",
       node.size = "node.size.betweenDegree",
       node.alpha = "node.alpha.degree",
       node.shape = "node.shape.betweenDegree",
       # edges
       edge.size = "common.edges.width",
       edge.color = "edgecolor", 
       
       # layout coords
       mode = c("x","y"), 
       # layout pars
       layout.par = ggnet.layout.par
       # neighbour
       #na.rm = na.rm
) +
  guides(size=FALSE,alpha=FALSE)  + 
  ggplot2::ggtitle("node.color.exptime, node.size.betweenDegree, node.alpha.degree, node.shape.betweenDegree,
                   edgesize=common.edges.width, edgecol=weight")

ggnet2(igraph.dynGENIE, 
       # nodes
       node.color = "node.color.TF",
       node.size = "node.size.degree",
       node.alpha = "node.alpha.degree",
       node.shape = "node.shape.TF",
       # edges
       edge.size = "common.edges.width",
       edge.color = "edgecolor", 
       
       # layout coords
       mode = c("x","y"), 
       # layout pars
       layout.par = ggnet.layout.par
       # neighbour
       #na.rm = na.rm
) +
  guides(size=FALSE,alpha=FALSE)  + 
  ggplot2::ggtitle("node.color.TF, node.size.degree, node.alpha.degree, node.shape.TF,
                   edgesize=common.edges.width, edgecol=weight")


ggnet2(igraph.dynGENIE, 
       # nodes
       node.color = "membership.dir",
       node.size = "node.size.betweenDegree",
       node.alpha = "node.alpha.degree",
       node.shape = "node.shape.betweenDegree",
       # edges
       edge.size = "common.edges.width",
       edge.color = "edgecolor", 
       
       # layout coords
       mode = c("x","y"), 
       # layout pars
       layout.par = ggnet.layout.par
       # neighbour
       #na.rm = na.rm
) +
  guides(size=FALSE,alpha=FALSE)  + 
  ggplot2::ggtitle("col=membership.dir, node.size.betweenDegree, node.alpha.degree, node.shape.betweenDegree,
                   edgesize=common.edges.width, edgecol=weight")

ggnet2(igraph.dynGENIE, 
       # nodes
       node.color = "membership.dir",
       node.size = "node.size.betweenDegree",
       node.alpha = "node.alpha.degree",
       node.shape = "node.shape.betweenDegree",
       # edges
       edge.size = "common.edges.width",
       edge.color = "edgecolor", 
       
       # layout coords
       mode = c("x","y"), 
       # layout pars
       layout.par = ggnet.layout.par
       # neighbour
       #na.rm = na.rm
) +
  guides(size=FALSE,alpha=FALSE)  + 
  ggplot2::ggtitle("col=membership.dir, node.size.betweenDegree, node.alpha.degree, node.shape.betweenDegree,
                   edgesize=common.edges.width, edgecol=weight")


ggnet2(igraph.dynGENIE, 
       # nodes
       node.color = "membership.dir",
       node.size = "node.size.betweenDegree",
       node.alpha = "node.alpha.EigCent",
       node.shape = "node.shape.betweenDegree",
       # edges
       edge.size = "common.edges.width",
       edge.color = "edgecolor", 
       
       # layout coords
       mode = c("x","y"), 
       # layout pars
       layout.par = ggnet.layout.par
       # neighbour
       #na.rm = na.rm
) +
  guides(size=FALSE,alpha=FALSE)  + 
  ggplot2::ggtitle("col=membership.dir, node.size.betweenDegree, node.alpha.EigCent, node.shape.betweenDegree,
                   edgesize=common.edges.width, edgecol=weight")

ggnet2(igraph.dynGENIE, 
       # nodes
       node.color = "node.color.TF.heatmapcluster",
       node.size = "node.size.betweenDegree",
       node.alpha = "node.alpha.degree",
       node.shape = "node.shape.betweenDegree",
       # edges
       edge.size = "common.edges.width",
       edge.color = "edgecolor", 
       
       # layout coords
       mode = c("x","y"), 
       # layout pars
       layout.par = ggnet.layout.par
       # neighbour
       #na.rm = na.rm
) +
  guides(size=FALSE,alpha=FALSE)  + 
  ggplot2::ggtitle("node.color.TF.heatmapcluster, node.size.betweenDegree, node.alpha.degree, 
                   node.shape.betweenDegree, edgesize=common.edges.width, edgecol=weight")

ggnet2(igraph.dynGENIE, 
       # nodes
       node.color = "node.color.TF.heatmapcluster",
       node.size = "node.size.betweenDegree",
       node.alpha = "node.alpha.EigCent",
       node.shape = "node.shape.betweenDegree",
       # edges
       edge.size = "common.edges.width",
       edge.color = "edgecolor", 
       
       # layout coords
       mode = c("x","y"), 
       # layout pars
       layout.par = ggnet.layout.par
       # neighbour
       #na.rm = na.rm
) +
  guides(size=FALSE,alpha=FALSE)  + 
  ggplot2::ggtitle("node.color.TF.heatmapcluster, node.size.betweenDegree, node.alpha.EigCent, 
                   node.shape.betweenDegree, edgesize=common.edges.width, edgecol=weight")



dev.off()




























