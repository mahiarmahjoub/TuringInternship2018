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
colnames(expr) <- as.character(time)
rownames(expr) <- expr.all$tracking_id
expr.selected <- as.matrix(expr[gene.list$x,])  # pick selected genes
expr.selected.log <- log2(expr.selected+1)  # log2 normalise



## ==== dynGENIE3 ---------------------------------------------------------- ##

res.genereg <- read.csv(file.path(dir, "TDA_output",
                                  "geneReg_links_180702_Col_22.csv"), 
                        header = TRUE, stringsAsFactors = FALSE)[,2:6]
res.dynGENIE <- read.csv(file.path(dir, "TDA_output",
                                   "dynGENIE_links_180702_v3_Col_22.csv"), 
                         header = TRUE, stringsAsFactors = FALSE)

# plot the histogram of weights to assess the selection of a threshold
# par(mfrow=c(2,2))
# hist(res.dynGENIE$weight,  breaks = 300, main="dynGENIE weight histogram",
#      xlab="weight")
# hist(res.genereg$coef,  breaks = 30, main="GeneReg coef histogram", 
#      xlab="coefficient")
# hist(res.genereg$delay,  breaks = 30, main="GeneReg delay histogram",
#      xlab="delay (mins)")
# hist(res.genereg$adj.r.squared,  breaks = 150, 
#      main="GeneReg adjusted R^2 histogram", xlab="adjusted R squared")
# par(mfrow=c(1,1))


threshold.dynGENIE <- 0.02
threshold.genereg <- 0.9999999

# exclude edges below cutoff
res.dynGENIE.filt <- res.dynGENIE[res.dynGENIE$weight > threshold.dynGENIE,]
res.genereg.filt <- res.genereg[abs(res.genereg$coef) > threshold.genereg,]







## ==== set up IGRAPH object ----------------------------------------------- ##
## IGRAPH object 
vertices.available <- unique(c(res.dynGENIE.filt[,1],res.dynGENIE.filt[,2]))
igraph.dynGENIE <- graph.data.frame(res.dynGENIE.filt, vertices.available, 
                                    directed=T)




## ==== (1) plot NETWORK --------------------------------------------------- ##

## ---- obtain quantitative information first
# DEGREE 
V(igraph.dynGENIE)$"degree" <- igraph::degree(igraph.dynGENIE, 
                                              V(igraph.dynGENIE)$name)
# EDGE WEIGHT
E(igraph.dynGENIE)$"edgeweight" <- E(igraph.dynGENIE)$weight*5
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
V(igraph.dynGENIE)$"betweennes" <- igraph::betweenness(igraph.dynGENIE) # vertex
E(igraph.dynGENIE)$"betweennes" <- igraph::edge.betweenness(igraph.dynGENIE) # edge
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
# change the transparency of nodes according to their degree
V(igraph.dynGENIE)$"nodealpha" <- seq(0.1,1,
                                      length.out = max(igraph::degree(igraph.dynGENIE, 
                                                                      vertices.available)
                                      ))[igraph::degree(igraph.dynGENIE, V(igraph.dynGENIE)$name)]
# set the node colour to BLUE if TF
V(igraph.dynGENIE)$"nodecolor" <- ifelse(V(igraph.dynGENIE)$name %in% TF.list$x, 
                                         brew_colors(col = "blue"),
                                         brew_colors(col = "red"))
# set the edge colour proportional to its weight 
E(igraph.dynGENIE)$"edgecolor" <- colorRampPalette(c("black", 
                                                     "grey80"))(length(E(
                                                       igraph.dynGENIE)$weight))


## SCALAR STATISTICS
closeness.global <- igraph::centralization.closeness(igraph.dynGENIE)$centralization
global.clustering <- igraph::transitivity(as.undirected(igraph.dynGENIE), 
                                          type = "global")
shortest.path.matrix <- igraph::shortest.paths(igraph.dynGENIE, weights = NA)
diameter <- diameter(igraph.dynGENIE)
mean.shortest.path <- mean(shortest.path.matrix[shortest.path.matrix!=Inf])



## ---- plot GGNET2 (network object)
# ggnet2(igraph.dynGENIE, node.color = "TF", palette = "Set1",
#        node.size = "degree", edge.size = "edgeweight", 
#        edge.color = "edgecolor", node.alpha = "nodealpha", 
#        mode = "fruchtermanreingold", color.legend = "TF", 
#        layout.par = list(cell.jitter = 0.1)) +
#   guides(size=FALSE)




## ==== (2) Plots ---------------------------------------------------------- ##
## Degree distribution 
# par(mfrow=c(1,2))
# barplot(degree.distribution(igraph.dynGENIE), 
#         names.arg = as.character(1:length(degree.distribution(igraph.dynGENIE)
#         )
#         ), xlab = "degree (k)", ylab="frequency (p_k)", 
#         main = "degree distribution")
# dd.col <- c(rep("red", degree.threshold), 
#             rep("blue", 
#                 length(degree.distribution(igraph.dynGENIE))-degree.threshold)
#             )
# dd.col[degree.distribution(igraph.dynGENIE)==0] <- "grey20"
# plot(c(0.9, 1:(length(degree.distribution(igraph.dynGENIE))-1)), 
#      degree.distribution(igraph.dynGENIE)+0.001, log = "xy",
#      main = "degree distribution (logscale)", 
#      xlab = "degree (k)", ylab="frequency (p_k) + 0.001", pch=16, cex=2, 
#      col = dd.col
#      )
# par(mfrow=c(1,1))


# par(mfrow=c(1,2))
# ckk.col <- ifelse(V(igraph.dynGENIE)$"degree" > degree.threshold, "blue","red")
# ## k vs C(k) - clustering coefficient 
# plot(V(igraph.dynGENIE)$"degree", V(igraph.dynGENIE)$"clustCoeff", 
#      pch=16, xlab="degree (k)", ylab="Clustering coefficient C(k)", 
#      main="clustering k vs C(k)", col=ckk.col)
# 
# ## distance histogram
# hist(shortest.path.matrix, xlab="shortest path", 
#      main=paste0("shortest path histogram \n diameter = ",
#                  round(diameter,3),", mean = ", round(mean.shortest.path,3)))
# par(mfrow=c(1,1))






## ==== (3) COMPARISON ----------------------------------------------------- ##
dynGENIE.cutoff <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.015, 0.02, 
                     0.025, 0.03)
precision <- recall <- specificity <- 
  sensitivity <- rep(NA, length(dynGENIE.cutoff))

for (i in 1:length(dynGENIE.cutoff)){
  compare <- network.pairwise.edge.comparison(
    sample.list = res.dynGENIE[res.dynGENIE[,3]>dynGENIE.cutoff[i],], 
    reference.list = res.genereg, node.names = gene.list$x)
  
  precision[i] <- compare$metrics["prec"]
  recall[i] <- compare$metrics["recall"]
  specificity[i] <- compare$metrics["specificity"]
  sensitivity[i] <- compare$metrics["sensitivity"]
}

# par(mfrow=c(1,2))
# plot(specificity, sensitivity, type="b", 
#      main="ROC Curve", 
#      col=colorRampPalette(c("red","blue"))(length(dynGENIE.cutoff)), 
#      pch=16)
# plot(recall, precision, type="b", main="Precision-Recall curve",
#      col=colorRampPalette(c("red","blue"))(length(dynGENIE.cutoff)), pch=16)
# legend("topright", legend = as.character(dynGENIE.cutoff[!is.na(sensitivity)]), 
#        fill = colorRampPalette(c("red","blue"))(length(dynGENIE.cutoff)))
# par(mfrow=c(1,1))







## ==== (4) TF link HEATMAP ------------------------------------------------ ##
adj.mat.dynGENIE.TFs <- get.adjmatrix.from.list(res.dynGENIE.filt, gene.list$x)
adj.mat.dynGENIE.TFs <- adj.mat.dynGENIE.TFs[TF.list$x, TF.list$x]
connected.TFs <- rowSums(adj.mat.dynGENIE.TFs) > 0.02 | 
  colSums(adj.mat.dynGENIE.TFs) > 0.02
# heatmap of all TFs with h-clustering
# heatmap.2(adj.mat.dynGENIE.TFs, 
#           trace = 'none', labRow = NA, labCol = NA,
#           col=brewer.pal(5,"OrRd"), dendrogram = 'none', 
#           main = "Heatmap of TF connectivity", key.title = "edge weight"
#           )
# # heatmap of TFs edges with weight>0 only 
# heatmap.2(adj.mat.dynGENIE.TFs[connected.TFs,connected.TFs], 
#           trace = 'none', labRow = NA, labCol = NA,
#           col=brewer.pal(5,"OrRd"), 
#           main = "Selected TF connectivity", key.title = "edge weight"
#           )





## ==== (5) FUNCTIONAL CONNECTIVITY HEATMAP -------------------------------- ##
# assign a GO term to each vertex 
gene.list.GO.modified <- matrix(NA, length(gene.list$x), 2, 
                                dimnames = list(NULL, c("gene","GO"))
)
gene.list.GO.modified[,"gene"] <- gene.list$x
gene.list.GO.modified[,"GO"] <- apply(gene.list.GO[gene.list.GO$V1 %in% gene.list$x, 3:7], 
                                      1, 
                                      function(x){
                                        if(grepl("DNA binding",paste(x,collapse = "_"))){
                                          "DNA binding"
                                        } else if (grepl("biosynthetic",paste(x,collapse = "_"))){
                                          "biosynthetic"
                                        } else if (grepl("signaling",paste(x,collapse = "_"))){
                                          "signalling"
                                        } else if (grepl("light",paste(x,collapse = "_"))){
                                          "light"
                                        } else if (grepl("stress",paste(x,collapse = "_"))){
                                          "stress"
                                        } else if (grepl("phosphorylation",paste(x,collapse = "_"))){
                                          "phosphorylation"
                                        } else if (grepl("metabolic",paste(x,collapse = "_"))){
                                          "metabolic"
                                        } else if (grepl("heat",paste(x,collapse = "_")) | 
                                                   grepl("cold",paste(x,collapse = "_"))){
                                          "temperature"
                                        } else if (grepl("circadian",paste(x,collapse = "_"))){
                                          "circadian"
                                        } else if (grepl("water",paste(x,collapse = "_"))){
                                          "water"
                                        } else if (grepl("salt",paste(x,collapse = "_"))){
                                          "salt"
                                        } else if (grepl("auxin",paste(x,collapse = "_")) | 
                                                   grepl("jasmonic",paste(x,collapse = "_")) | 
                                                   grepl("absiscic",paste(x,collapse = "_")) |
                                                   grepl("ethylene",paste(x,collapse = "_"))){
                                          "hormone"
                                        } else if (grepl("kinase",paste(x,collapse = "_"))){
                                          "kinase"
                                        } else if (grepl("transferase",paste(x,collapse = "_"))){
                                          "transferase"
                                        } else if (grepl("transporter",paste(x,collapse = "_"))){
                                          "transporter"
                                        } else if (grepl("protein binding",paste(x,collapse = "_"))){
                                          "protein binding"
                                        }
                                        else {"other"}
                                      }
)


# create a adjacency matrix of GO terms 
GO.terms.interest <- unique(gene.list.GO.modified[,"GO"])
GO.adj.matrix <- matrix(0, length(GO.terms.interest), length(GO.terms.interest), 
                        dimnames = list(GO.terms.interest,GO.terms.interest))
edge.list.GO.link <- igraph::get.edgelist(igraph.dynGENIE)

for (i in 1:nrow(edge.list.GO.link)){
  target <- gene.list.GO.modified[
    gene.list.GO.modified[,"gene"]==edge.list.GO.link[i,1],"GO"]
  regulator <- gene.list.GO.modified[
    gene.list.GO.modified[,"gene"]==edge.list.GO.link[i,2],"GO"]
  GO.adj.matrix[target,regulator] <- GO.adj.matrix[target,regulator] + 1
}

# normalise: Gij/(Gi*Gj)
# GO.adj.matrix.norm <- GO.adj.matrix
# for (i in 1:nrow(GO.adj.matrix)){
#   for (j in 1:ncol(GO.adj.matrix)){
#     GO.adj.matrix.norm[i,j] <- GO.adj.matrix[i,j]/(sum(
#       GO.adj.matrix[i,])*sum(GO.adj.matrix[,j]))
#   }
# }
# GO.adj.matrix.norm[is.nan(GO.adj.matrix.norm)] <-0 

heatmap.2(GO.adj.matrix, Rowv = NA, Colv = NA, dendrogram = 'none', 
          trace = 'none', 
          col=c("lightyellow",colorRampPalette(c("orange","red","red3"))(100)), 
          main = "Functional connectivity \n (based on GO terms)")



## ---- GO connectivities of TFs only based on their regulatory job -------- ##
# pick TFs (DNA binding) only 
gene.list.GO.modified.TF <- gene.list.GO.modified[gene.list.GO.modified[,2]=="DNA binding",]
# pick edges containing DNA-binding genes only 
edge.list.GO.link.TF <- edge.list.GO.link[edge.list.GO.link[,1] %in% 
                                            gene.list.GO.modified.TF[,1] & 
                                            edge.list.GO.link[,2] %in% 
                                            gene.list.GO.modified.TF[,1],
                                          ]
GO.terms.interest.TF <- c("circadian","light","metabolic","signaling",
                          "stress", "heat", "cold", "water", "flower", "stimulus",
                          "biosynthetic","phosphorylation",
                          "temperature", "photosynthesis", "kinase","transporter",
                          "transferase", "ion binding")   
# replace 'DNA binding' with their regulatory job e.g. temp, salt etc. regulation
# collate all GO together for each gene with "DNA binding" GO
GO.adj.matrix.TF <- matrix(0,length(GO.terms.interest.TF), 
                           length(GO.terms.interest.TF), 
                           dimnames = list(GO.terms.interest.TF,
                                           GO.terms.interest.TF)
)
GO.adj.matrix.TF <- rbind(GO.adj.matrix.TF,other=rep(0,ncol(GO.adj.matrix.TF)))
GO.adj.matrix.TF <- cbind(GO.adj.matrix.TF,other=rep(0,nrow(GO.adj.matrix.TF)))
GO.adj.matrix.TF <- rbind(GO.adj.matrix.TF,`DNA binding`=rep(0,ncol(GO.adj.matrix.TF)))
GO.adj.matrix.TF <- cbind(GO.adj.matrix.TF,`DNA binding`=rep(0,nrow(GO.adj.matrix.TF)))


gene.list.GO.modified.TF[,"GO"] <- apply(gene.list.GO[gene.list.GO$V1 %in% 
                                                        gene.list.GO.modified.TF[,1], 3:7], 
                                         1, 
                                         function(x){
                                           x <- paste(x,collapse = "_")
                                           x <- unique(unlist(strsplit(x,split = "_")))
                                           paste(x,collapse = "_")
                                         })

# fill the adjacency matrix for the TF roles connectivity                                            
for (i in 1:nrow(edge.list.GO.link.TF)){
  # find genes 
  regulator.index <- which(gene.list.GO.modified.TF[,1]==edge.list.GO.link.TF[i,1])
  target.index <- which(gene.list.GO.modified.TF[,1]==edge.list.GO.link.TF[i,2])
  # extract functions of each gene 
  reg.fun <- as.vector(unlist(strsplit(gene.list.GO.modified.TF[regulator.index,2],
                                       split = "_")))
  target.fun <- as.vector(unlist(strsplit(gene.list.GO.modified.TF[target.index,2],
                                          split = "_")))
  
  reg.fun <- sapply(GO.terms.interest.TF, function(x){grepl(x,reg.fun)})
  target.fun <- sapply(GO.terms.interest.TF, function(x){grepl(x,target.fun)})
  
  reg.fun <- apply(reg.fun,2,sum); target.fun <- apply(target.fun,2,sum)
  
  if (sum(reg.fun)==0){
    reg.fun <- as.vector(unlist(strsplit(gene.list.GO.modified.TF[regulator.index,2],
                                         split = "_")))
    if (sum(grepl("DNA binding", reg.fun))>0){
      reg.fun <- "DNA binding"
    } else {
      reg.fun <- "other"
    }
  } else {
    reg.fun <- names(reg.fun[reg.fun>0])
  }
  
  
  
  if (sum(target.fun)==0){
    target.fun <- as.vector(unlist(strsplit(gene.list.GO.modified.TF[target.index,2],
                                            split = "_")))
    if (sum(grepl("DNA binding", target.fun))>0){
      target.fun <- "DNA binding"
    } else {
      target.fun <- "other"
    }
  } else {
    target.fun <- names(target.fun[target.fun>0])
  }
  
  for (i in reg.fun){
    for (j in target.fun){
      GO.adj.matrix.TF[i,j] <- GO.adj.matrix.TF[i,j] + 1
    }
  }
  
}



heatmap.2(GO.adj.matrix.TF,
          trace = 'none', key.title = "edge count" ,
          col=c("lightyellow",colorRampPalette(c("yellow",
                                                 "orange","red",
                                                 "red3"))(100)), 
          main = "Functional connectivity TFs \n (based on GO terms)")








## ==== (6) EXPRESSION TIME HEATMAP ---------------------------------------- ##

edge.list.peakExprTime <- igraph::get.edgelist(igraph.dynGENIE)
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


heatmap.2(peakExpr.adj.matrix, Rowv = NA, Colv = NA, dendrogram = 'none', 
          trace = 'none', col=brewer.pal(5,"Blues"), 
          main="Peak expression time connectivity")










## == PLOT EVERYTHING ====================================================== ##

pdf(file.path(dir,"network_summary_Col_22.pdf"), width = 7, height = 5)
## -- (0) TDA scores -----------
par(mfrow=c(2,2))
hist(res.dynGENIE$weight,  breaks = 300, main="dynGENIE weight histogram",
     xlab="weight")
hist(res.genereg$coef,  breaks = 30, main="GeneReg coef histogram", 
     xlab="coefficient")
hist(res.genereg$delay,  breaks = 30, main="GeneReg delay histogram",
     xlab="delay (mins)")
hist(res.genereg$adj.r.squared,  breaks = 150, 
     main="GeneReg adjusted R^2 histogram", xlab="adjusted R squared")
par(mfrow=c(1,1))

## -- (1) NETWORK ---------------
ggnet2(igraph.dynGENIE, node.color = "TF", palette = "Set1",
       node.size = "degree", edge.size = "edgeweight", 
       edge.color = "edgecolor", node.alpha = "nodealpha", 
       mode = "fruchtermanreingold", color.legend = "TF", 
       layout.par = list(cell.jitter = 0.1)) +
  guides(size=FALSE)


## -- (2.1) DEGREE DISTRIBUTION --------------------
par(mfrow=c(1,2))
barplot(degree.distribution(igraph.dynGENIE), 
        names.arg = as.character(1:length(degree.distribution(igraph.dynGENIE)
        )
        ), xlab = "degree (k)", ylab="frequency (p_k)", 
        main = "degree distribution")
dd.col <- c(rep("red", degree.threshold), 
            rep("blue", 
                length(degree.distribution(igraph.dynGENIE))-degree.threshold)
)
dd.col[degree.distribution(igraph.dynGENIE)==0] <- "grey20"
plot(c(0.9, 1:(length(degree.distribution(igraph.dynGENIE))-1)), 
     degree.distribution(igraph.dynGENIE)+0.001, log = "xy",
     main = "degree distribution (logscale)", 
     xlab = "degree (k)", ylab="frequency (p_k) + 0.001", pch=16, cex=2, 
     col = dd.col
)
par(mfrow=c(1,1))



## -- (2.2) CLUSTERING & PATH distribution -----------------------------
par(mfrow=c(1,2))
ckk.col <- ifelse(V(igraph.dynGENIE)$"degree" > degree.threshold, "blue","red")
## k vs C(k) - clustering coefficient 
plot(V(igraph.dynGENIE)$"degree", V(igraph.dynGENIE)$"clustCoeff", 
     pch=16, xlab="degree (k)", ylab="Clustering coefficient C(k)", 
     main="clustering k vs C(k)", col=ckk.col)

## distance histogram
hist(shortest.path.matrix, xlab="shortest path", 
     main=paste0("shortest path histogram \n diameter = ",
                 round(diameter,3),", mean = ", round(mean.shortest.path,3)))
par(mfrow=c(1,1))




## -- (3) COMPARISON: Pre-Recall + AROC: GeneReg vs dynGENIE -------------------
par(mfrow=c(1,2))
plot(specificity, sensitivity, type="b", 
     main=paste0("ROC Curve \n GR.adjR.cutoff = ",threshold.genereg), 
     col=colorRampPalette(c("red","blue"))(length(dynGENIE.cutoff)), 
     pch=16)
plot(recall, precision, type="b", main="Precision-Recall curve",
     col=colorRampPalette(c("red","blue"))(length(dynGENIE.cutoff)), pch=16)
legend("topright", legend = as.character(dynGENIE.cutoff[!is.na(sensitivity)]), 
       fill = colorRampPalette(c("red","blue"))(length(dynGENIE.cutoff)))
par(mfrow=c(1,1))




## -- (4) TF LINK HEATMAP -----------------------------------------------
heatmap.2(adj.mat.dynGENIE.TFs, 
          trace = 'none', labRow = NA, labCol = NA,
          #col=brewer.pal(7,"OrRd"), 
          col = c("lightyellow",colorRampPalette(c("orange","red",
                                                   "red3"))(100)),
          dendrogram = 'none', 
          main = "Heatmap of TF connectivity", key.title = "edge weight"
)
# heatmap of TFs edges with weight>0 only 
heatmap.2(adj.mat.dynGENIE.TFs[connected.TFs,connected.TFs], 
          trace = 'none', labRow = NA, labCol = NA,
          #col=brewer.pal(5,"OrRd"), 
          col = c("lightyellow",colorRampPalette(c("orange","red",
                                                   "red3"))(100)),
          main = "Selected TF connectivity", key.title = "edge weight"
)






## -- (5) FUNCTIONAL CONNECTIVITY HEATMAP -------------------------------
heatmap.2(GO.adj.matrix, Rowv = NA, Colv = NA, dendrogram = 'none', 
          trace = 'none', 
          col=c("lightyellow",colorRampPalette(c("orange","red","red3"))(100)), 
          main = "Functional connectivity \n (based on GO terms)")

# TFs only 
heatmap.2(GO.adj.matrix.TF,
          trace = 'none', key.title = "edge count" ,
          col=c("lightyellow",colorRampPalette(c("yellow",
                                                 "orange","red",
                                                 "red3"))(100)), 
          main = "Functional connectivity TFs \n (based on GO terms)")



## -- (6) EXPRESSION TIME HEATMAP ---------------------------------------
heatmap.2(peakExpr.adj.matrix, Rowv = NA, Colv = NA, dendrogram = 'none', 
          trace = 'none', col=brewer.pal(5,"Blues"), 
          main="Peak expression time connectivity")

dev.off()















## ==== (7) investigate CC activity  --------------------------------------- ##







# # finds genes with links going out only 
# out.only <- V(igraph.dynGENIE)$name[centr_degree(igraph.dynGENIE, 
#                                                  mode="in")$res==0]
# 
# # pick genes that have non-NaN and non-zero clustering coefficient 
# x <- V(igraph.dynGENIE)$name[igraph::transitivity(igraph.dynGENIE, 
#                                                   type = "local", vids = V(igraph.dynGENIE)) > 0 &
#                                !is.nan(igraph::transitivity(igraph.dynGENIE, type = "local", 
#                                                             vids = V(igraph.dynGENIE)))]
# 
# x <- unlist(strsplit(paste(c(gene.list.GO[,3], gene.list.GO[,4], 
#                              gene.list.GO[,5], gene.list.GO[,6], 
#                              gene.list.GO[,7]), collapse = "_"), split = "_"))
# table(x)

