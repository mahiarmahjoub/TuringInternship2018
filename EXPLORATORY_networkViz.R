## =================== Network visualisation exploratory =================== ##

## -- simulate a random network -------------------------------------------- ## 
n.nodes <- 20
node.names <- LETTERS[1:n.nodes]
adj.matrix <- matrix(sample(c(0,1),n.nodes^2, c(0.8,0.2),replace = TRUE),
                     n.nodes,
                     n.nodes,dimnames = list(node.names, node.names))
adj.matrix[,"A"] <- 0
adj.matrix[,"S"] <- 0
adj.matrix["S",] <- 0



## ==== VISUALISATION ------------------------------------------------------ ##
## ---- ggnet -------------------------------------------------------------- ##
# find ways to draw the networks in very pretty ways 
# https://briatte.github.io/ggnet/

library(GGally)
library(network)
library(sna)
library(ggplot2)
library(RColorBrewer)
library(intergraph)
library(ndtv)

# use adjacency matrix to construct 'network' graph object
net <- rgraph(10, mode = "graph", tprob = 0.5)
net <- network(net, directed = FALSE)

# vertex names
network.vertex.names(net) = letters[1:10]

ggnet2(net, mode = "target", edge.color = "black", label=TRUE)

# vertex attribute 
net %v% "phono" <- ifelse(letters[1:10] %in% c("a", "e", "i"), 
                          "vowel", "consonant")
net %v% "color" = ifelse(net %v% "phono" == "vowel", 
                         "steelblue", "tomato")
ggnet2(net, color = "color")
ggnet2(net, color = "phono", palette = "Set1")
ggnet2(net, size = "phono", size.palette = c("vowel"=10, "consonant"=1))
ggnet2(net, size = "degree")  # outdegree, indegree, degree
ggnet2(net, size = sample(0:2, 10, replace = TRUE), size.zero = TRUE) # plot zero sized nodes 
ggnet2(net, size = "degree", size.min = 5)  # remove smaller sized nodes

# remove legends 
ggnet2(net, color = "phono", size = "degree") +
  guides(color = FALSE, size = FALSE)
ggnet2(net, alpha = "phono", alpha.legend = "Phonetics")
ggnet2(net, shape = "phono", shape.legend = "Phonetics")
ggnet2(net, color = "phono", color.legend = "Phonetics")
ggnet2(net, size = "degree", size.legend = "Centrality")

ggnet2(net, color = "phono", legend.size = 12, 
       legend.position = "bottom", node.label = TRUE, 
       label.color = "white"
       ) +
  theme(panel.background = element_rect(color = "grey"))

ggnet2(net, alpha = "phono", 
       alpha.palette = c("vowel" = 0.2, "consonant" = 1), label = TRUE)



## bipartite network 
# weighted adjacency matrix
bip = data.frame(event1 = c(1, 2, 1, 0),
                 event2 = c(0, 0, 3, 0),
                 event3 = c(1, 1, 0, 4),
                 row.names = letters[1:4])

# weighted bipartite network
bip = network(bip,
              matrix.type = "bipartite",
              ignore.eval = FALSE,
              names.eval = "weights")
# set colors for each mode
col = c("actor" = "grey", "event" = "gold")
# detect and color the mode
ggnet2(bip, color = "mode", palette = col, label = TRUE, 
       edge.label = "weights")

# change edge size 
ggnet2(bip, color = "mode", palette = col, edge.size = "weights")
bip %e% "color" <- ifelse(bip %e% "weights" > 1, "black", "grey75")
ggnet2(bip, color = "mode", palette = col, edge.size = "weights", 
       edge.color = "color")
bip %e% "lty" <- ifelse(bip %e% "weights" >1, 1, 2)
ggnet2(bip, color = "mode", palette = col, edge.size = "weights", 
       edge.lty = "lty")

# icelandic legal code 
source("https://goo.gl/q1JFih")
x = cut_number(as.integer(net %v% "year"), 4)
col = c("#E1AF00", "#EBCC2A", "#78B7C5", "#3B9AB2")
names(col) = levels(x)
ggnet2(net, color = x, color.legend = "period", palette = col,
       edge.alpha = 1/4, edge.size = "weight",
       size = "degree", max_size = 4, size.cut = 3,
       legend.size = 12, legend.position = "bottom") +
  coord_equal()






## ==== GRAPH INFORMATION -------------------------------------------------- ##
## ---- igraph ------------------------------------------------------------- ##
# to extract the following information from the graph:
# (1) number of nodes and vertices
# (2) number of connections per node
# (3) number of parent (upstream) and children (downstream) nodes 
# (4) shortest path between two nodes 

library(igraph)

nodes <- read.csv("Dataset1-Media-Example-NODES.csv", header=T, as.is=T)
links <- read.csv("Dataset1-Media-Example-EDGES.csv", header=T, as.is=T)
nodes2 <- read.csv("Dataset2-Media-User-Example-NODES.csv", header=T, as.is=T)
links2 <- read.csv("Dataset2-Media-User-Example-EDGES.csv", header=T, row.names=1)

# if non-unique links present, then aggregate (their weight) into one link
links <- aggregate(links[,3], links[,-3], sum)
links <- links[order(links$from, links$to),]
colnames(links)[4] <- "weight"
rownames(links) <- NULL

# create igraph object 
net <- graph.data.frame(links, nodes, directed=T)
# access edges and vertices 
E(net)
V(net)
# manipulate the matrix directly 
net[1,]
# plot
plot(net)
# improve plot 
net <- simplify(net, remove.multiple = F, remove.loops = T) 
plot(net, edge.arrow.size=0.4, vertex.label=NA)
# colours 
pal2 <- rainbow(5, alpha=.5); pal1 <- heat.colors(5, alpha=1)
grDevices::adjustcolor("557799", alpha=0.7)
colorRampPalette(c("gray80", "dark red"))(10)
colorRampPalette(c(rgb(1,1,1, .2),rgb(.8,0,0, .7)), alpha=TRUE)(10)
library(RColorBrewer) # <<<<<-----------------!!!!!!!!!!!! best
display.brewer.pal(8, "Set3") # Spectral, Blues


# you can set graphing parameters directly
plot(net, edge.arrow.size=.2, edge.color="orange",
     vertex.color="orange", vertex.frame.color="#ffffff",
     vertex.label=V(net)$media, vertex.label.color="black") 
# OR change them in the igraph object 
colrs <- c("gray50", "tomato", "gold")
V(net)$color <- colrs[V(net)$media.type] # colours based on media.type
deg <- degree(net, mode="all")
V(net)$size <- deg*3
E(net)$width <- E(net)$weight/6  # Set edge width based on weight
#change arrow size and edge color:
E(net)$arrow.size <- .2
E(net)$edge.color <- "gray80"
E(net)$width <- 1+E(net)$weight/12
plot(net) 
legend(x=-1.5, y=-1.1, c("Newspaper","Television", "Online News"), pch=21,
       col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)

# color edge based on source node 
edge.start <- get.edges(net, 1:ecount(net))[,1]
edge.col <- V(net)$color[edge.start]
plot(net, edge.color=edge.col, edge.curved=.1)

net.bg <- barabasi.game(80) 
V(net.bg)$frame.color <- "white"
V(net.bg)$color <- "orange"
V(net.bg)$label <- "" 
V(net.bg)$size <- 10
E(net.bg)$arrow.mode <- 0
plot(net.bg, layout=layout.circle)

# force-directed plot 
l <- layout.fruchterman.reingold(net.bg, repulserad=vcount(net.bg)^3, 
                                 area=vcount(net.bg)^2.4)
plot(net.bg, layout=layout.fruchterman.reingold)
plot(net.bg, layout=l)
# another force-directed plot 
l <- layout.kamada.kawai(net.bg)
plot(net.bg, layout=l)
# LGL useful for large connected networks 
plot(net.bg, layout=layout.lgl)

# all layout options 
grep("^layout\\.", ls("package:igraph"), value=TRUE)

net - E(net)[E(net)$type=="hyperlink"] # another way to delete edges

# place nodes in communities
V(net)$community <- optimal.community(net)$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)
plot(net, vertex.color=colrs[V(net)$community])

# colour nodes depending on distance from a chosen node 
# SHORTEST.PATH  ....
# NEIGHBOURS finds first degree connected nodes 
# INCIDENT for edges 
dist.from.NYT <- shortest.paths(net, algorithm="unweighted")[1,]
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.NYT)+1)[dist.from.NYT+1]
plot(net, vertex.color=col, vertex.label=dist.from.NYT, edge.arrow.size=.6, 
     vertex.label.color="white")

# mark a group of ndoes 
# Mark multiple groups:
plot(net, mark.groups=list(c(1,4,5,8), c(15:17)), 
     mark.col=c("#C5E5E7","#ECD89A"), mark.border=NA)

V(net)$label <- NA
# highlight a path in network 
news.path <- get.shortest.paths(net, V(net)[media=="MSNBC"], 
                                V(net)[media=="New York Post"],
                                mode="all", output="both")
ecol <- rep("gray80", ecount(net))
ecol[unlist(news.path$epath)] <- "orange"
# Generate edge width variable:
ew <- rep(2, ecount(net))
ew[unlist(news.path$epath)] <- 4
# Generate node color variable:
vcol <- rep("gray40", vcount(net))
vcol[unlist(news.path$vpath)] <- "gold"
plot(net, vertex.color=vcol, edge.color=ecol, 
     edge.width=ew, edge.arrow.mode=0)


## DEGREE DISTRIBUTION
dd <- degree.distribution(net, cumulative=T, mode="all")
plot(dd, pch=19, cex=1, col="orange", xlab="Degree", 
     ylab="Cumulative Frequency")





## ------- work with NETWORK package --------------------------------------- ## 
library(network)
net3 <- network(links,  vertex.attr=nodes, matrix.type="edgelist", 
                loops=F, multiple=F, ignore.eval = F)
net3[,]  # access adjacency matrix 
net3 %v% "col" <- c("gray70", "tomato", "gold")[net3 %v% "media.type"]
plot(net3, vertex.cex=(net3 %v% "audience.size")/7, vertex.col="col")
# allows you to store the coordinates of a plot and use later on 
l <- plot(net3, vertex.cex=(net3 %v% "audience.size")/7, vertex.col="col")
plot(net3, vertex.cex=(net3 %v% "audience.size")/7, vertex.col="col", coord=l)





## ------- work with NETWORKD3 package ------------------------------------- ## 
library(networkD3)
# need to give numeric IDs to all nodes, starting from 0
el <- data.frame(from=as.numeric(factor(links$from))-1, 
                 to=as.numeric(factor(links$to))-1 )
nl <- cbind(idn=factor(nodes$media, levels=nodes$media), nodes) 

D3.file <- forceNetwork(Links = el, Nodes = nl, Source="from", Target="to",
             NodeID = "idn", Group = "type.label",linkWidth = 1,
             linkColour = "#afafaf", fontSize=12, zoom=T, legend=T,
             Nodesize=6, opacity = 0.8, charge=-600, 
             width = 600, height = 400)

networkD3::saveNetwork(D3.file, "D3_test.html")












