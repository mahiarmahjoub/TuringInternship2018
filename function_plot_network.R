## ========================== plot gene networks =========================== ##

## can use a (1) link list or an (2) adjacency matrix  
## using RgraphViz 
require(Rgraphviz, quietly = TRUE)


construct.network.graph <- function(x, input.mode=c("matrix","list"), 
                         edge.mode=c("directed","undirected")){
  
  # collate all the genes in the list
  genes <- unique(c(as.character(as.matrix(x[,1])), 
                    as.character(as.matrix(x[,2]))
                    )
                  )
  
  # setup graph objects based on the list/matrix provided 
  if(input.mode=="list"){
    x.graph <- new("graphNEL", nodes=genes, edgemode=edge.mode)
    x.graph <- addEdge(from = as.character(as.matrix(x[,1])), 
                       to = as.character(as.matrix(x[,2])), 
                       graph = x.graph, 
                       weights = as.numeric(x[,3])
                       )
  } else if (input.mode=="matrix"){
    x.graph <- new("graphAM", adjMat=x, edgemode=edge.mode)
  } else {print("Object to be drawn not compatible!")}
  
  return(graph.object=x.graph)
  
  
}
