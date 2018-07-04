## ===================== network visualisation packages ==================== ##


## ==== netbiov ------------------------------------------------------------ ##
# nice package but I think it is more useful for thesis/paper
biocLite("netbiov")
library("netbiov")
library("igraph")

data("artificial2.graph")
xx <- plot.modules(g1,mod.lab=TRUE, color.random=TRUE, mod.edge.col="grey",
                   ed.color="gold", sf=15, v.size=.5,
                   layout.function=layout.fruchterman.reingold,
                   lab.color="grey", modules.name.num=FALSE, 
                   lab.cex=1, lab.dist=5)


## ==== Rgraphviz ---------------------------------------------------------- ##
# https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/r/rgraphviz/
library("Rgraphviz")
library("graph")

test.matrix<-matrix(rep(c(0,1,0,0), 9), ncol=6, nrow=6)
rownames(test.matrix)<-c("a", "b", "c", "d", "e", "f")
colnames(test.matrix)<-c("a", "b", "c", "d", "e", "f")
am.graph<-new("graphAM", adjMat=test.matrix, edgemode="directed")
nel.graph <- new("graphNEL", nodes=rownames(test.matrix), edgemode="directed")
nel.graph <- addEdge(sample(rownames(test.matrix), 10, replace = TRUE), 
                     sample(rownames(test.matrix), 10, replace = TRUE), 
                     graph = nel.graph, weights = rnorm(10))

plot(nel.graph, attrs = list(node = list(fillcolor = "lightblue"),
                              edge = list(arrowsize=0.5)))



