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
source(file.path(dir,"TuringInternship2018/network_stats_functions.R"))


## ==== load expression matrix, wanted gene list etc. ---------------------- ##
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
# temp <- "22"; mut <- "Col"
# time <- c(-10, 10, 24, 30, 45, 60, 105, 120)
# # pick selected experimental conditions 
# expr <- expr.all[,grepl(temp,colnames(expr.all)) & grepl(mut,colnames(expr.all))]
# colnames(expr) <- as.character(time)
# rownames(expr) <- expr.all$tracking_id
# expr.selected <- as.matrix(expr[gene.list$x,])  # pick selected genes 
# expr.selected.log <- log2(expr.selected+1)  # log2 normalise



## ==== dynGENIE3 ---------------------------------------------------------- ##

res.dynGENIE <- read.csv(file.path(dir, "TDA_output",
                                   "dynGENIE_links_180702_v3_Col_22.csv"), 
                         header = TRUE, stringsAsFactors = FALSE)

# plot the histogram of weights to assess the selection of a threshold
hist(res.dynGENIE$weight,  breaks = 300)
threshold <- 0.02

# exclude edges below cutoff
res.dynGENIE.filt <- res.dynGENIE[res.dynGENIE$weight > threshold,]

# NETWORK object 
# net.dynGENIE <- network(res.dynGENIE.filt,
#                         matrix.type="edgelist",
#                         loops=F, multiple=F, ignore.eval = F)
# network.vertex.names(net.dynGENIE) <- gene.list$x
# net.dynGENIE %v% "degree" <- 
# ggnet2(net.dynGENIE, color = "degree")





## ==== set up IGRAPH object ----------------------------------------------- ##
## IGRAPH object 
vertices.available <- unique(c(res.dynGENIE.filt[,1],res.dynGENIE.filt[,2]))
igraph.dynGENIE <- graph.data.frame(res.dynGENIE.filt, vertices.available, 
                                    directed=T)



## ==== (1) plot NETWORK --------------------------------------------------- ##
# set size of each node to its degree/2
V(igraph.dynGENIE)$"degree" <- igraph::degree(igraph.dynGENIE, 
                                              V(igraph.dynGENIE)$name)
# set colour of each node to its degree
# set the edge width proportional to weight
E(igraph.dynGENIE)$"edgeweight" <- E(igraph.dynGENIE)$weight*5
# set the edge colour proportional to its weight 
E(igraph.dynGENIE)$"edgecolor" <- colorRampPalette(c("black", 
                                                     "grey80"))(length(E(
                                                       igraph.dynGENIE)$weight))
# set the node colour to BLUE if TF
V(igraph.dynGENIE)$"TF" <- ifelse(V(igraph.dynGENIE)$name %in% TF.list$x, "TF","other")
V(igraph.dynGENIE)$"nodecolor" <- ifelse(V(igraph.dynGENIE)$name %in% TF.list$x, 
                                         brew_colors(col = "blue"),
                                         brew_colors(col = "red"))
# change the transparency of nodes according to their degree
V(igraph.dynGENIE)$"nodealpha" <- seq(0.1,1,
                                      length.out = max(igraph::degree(igraph.dynGENIE, 
                                                                      vertices.available)
                                      ))[igraph::degree(igraph.dynGENIE, V(igraph.dynGENIE)$name)]
# place gene names on the high degree TF nodes 
V(igraph.dynGENIE)$"nodelabel" <- rep(NA, vcount(igraph.dynGENIE))



# plot GGNET2 (network object)
#pdf(file.path(dir,"test_network.pdf"), width = 7, height = 5)
ggnet2(igraph.dynGENIE, node.color = "TF", palette = "Set1",
       node.size = "degree", edge.size = "edgeweight", 
       edge.color = "edgecolor", node.alpha = "nodealpha", 
       mode = "fruchtermanreingold", color.legend = "TF") +
  guides(size=FALSE)
#dev.off()




## ==== (2) Degree distribution -------------------------------------------- ##
barplot(degree.distribution(igraph.dynGENIE), 
        names.arg = as.character(1:length(degree.distribution(igraph.dynGENIE)
        )
        ), xlab = "degree (k)", ylab="frequency (p_k)", 
        main = "degree distribution")



## ==== (3) COMPARISON ----------------------------------------------------- ##
edges.aranet.filt <- edges.aranet[(edges.aranet$regulator %in% gene.list$x) & 
                                    (edges.aranet$target %in% gene.list$x),]

compare <- network.pairwise.edge.comparison(sample.list = res.dynGENIE.filt, 
                                            reference.list = edges.aranet.filt, 
                                            node.names = gene.list$x)
# plot the common links
# calculate and plot ROC/AUROC/APROC etc.






## ==== (4) TF link HEATMAP ------------------------------------------------ ##








## ==== (5) METRICS -------------------------------------------------------- ##





## ==== (6) FUNCTIONAL CONNECTIVITY ---------------------------------------- ##






## ==== (7) COMPARISON ----------------------------------------------------- ##














