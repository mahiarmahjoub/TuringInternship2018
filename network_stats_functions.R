## ===================== functions to assess networks ====================== ##

# FUNCTIONS:
# get.adjmatrix.from.list
# get.list.from.adjmatrix
# network.pairwise.edge.comparison



## -- simulate a random network -------------------------------------------- ## 
n.nodes <- 20
node.names <- LETTERS[1:n.nodes]
adj.matrix <- matrix(sample(c(0,1),n.nodes^2, c(0.8,0.2),replace = TRUE),
                     n.nodes,
                     n.nodes,dimnames = list(node.names, node.names))
adj.matrix[,"A"] <- 0
adj.matrix[,"S"] <- 0
adj.matrix["S",] <- 0





## -- FUNCTION: extract adjacency matrix from list ------------------------- ##

get.adjmatrix.from.list <- function(link.list, node.names){
  n <- length(node.names)
  adjMat <- matrix(0, n, n, dimnames = list(node.names, node.names))
  
  adjMat[as.matrix(link.list[link.list[,3]>0,1:2])] <- link.list[link.list[,3]>0,3]
  
  return(as.matrix(adjMat))
}

#get.adjmatrix.from.list(link.list, node.names)






## -- FUNCTION: extract list from adjacency matrix ------------------------- ##

get.list.from.adjmatrix <- function(matrix, threshold=list(FALSE,0)){
  # extract number of nodes 
  n <- nrow(matrix)
  names <- rownames(matrix)
  # set up a list matrix 
  temp.list <- matrix(NA,n,3, 
                      dimnames = list(NULL,c("regulator","target","weight"))
  )
  connect.list <- NULL
  
  # fill the connectivity list matrix 
  for (i in 1:n){
    temp.list[,1] <- names[i]
    for (j in 1:n){
      temp.list[j,2] <- names[j]
      temp.list[j,3] <- matrix[names[i],names[j]]
    }
    connect.list <- rbind(connect.list, temp.list)
  }
  
  # change matrix into a data.frame with numeric weights
  connect.list <- as.data.frame(connect.list, stringsAsFactors=FALSE)
  connect.list$weight <- as.numeric(connect.list$weight)
  
  # remove links with weight less then the threshold
  if (isTRUE(threshold[[1]])){
    connect.list <- connect.list[connect.list[,3] > threshold[[2]],]
  }
  
  return(connect.list)
}

#list.links <- get.list.from.adjmatrix(adj.matrix, threshold = list(TRUE,0.1))







## -- FUNCTION: construct graph object for plotting network ---------------- ##
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

# graph <- construct.network.graph(x = list.links, input.mode = "list",
#                                  edge.mode = "directed")
# 
# Attr <- getDefaultAttrs()
# Attr$node$fixedsize <- FALSE
# Attr$node$fillcolor <- "lightblue"
# Attr$node$height <- "100"
# plot(graph, attrs = Attr)







## -- FUNCTION: find connections in/out of each node ----------------------- ##
get.first.connection.information <- function(links, link.type="list", node.names){
  
  if(link.type=="adj.matrix"){
    links <- get.list.from.adjmatrix(links, threshold = list(TRUE,0))
  } 
  
  # number of nodes
  n <- length(node.names)
  # initialise matrix holding connectivity information 
  connection.info <- matrix(NA, n, 5, 
                            dimnames = list(NULL,
                                            c("node","n.out",
                                              "n.in","wsum.out", 
                                              "wsum.in")
                            )
  )
  # fill in the node names 
  connection.info[,1] <- node.names
  for (i in 1:n){
    # pick connections with the wanted node as regulator 
    select.connections <- links[links[,1]==node.names[i],]
    # rows equal to the number of edges connected to our wanted node 
    connection.info[i,2] <- nrow(select.connections)
    connection.info[i,4] <- sum(select.connections[,3]) # sum of edge weights
    
    # pick connections with the wanted node as target
    select.connections <- links[links[,2]==node.names[i],]
    # rows equal to the number of edges connected to our wanted node 
    connection.info[i,3] <- nrow(select.connections)
    connection.info[i,5] <- sum(select.connections[,3]) # sum of edge weights
    
  }
  
  # convert matrix to data.frame and approrpriate data into numerics 
  connection.info <- as.data.frame(connection.info, stringsAsFactors=FALSE)
  class(connection.info$n.out) <- "numeric"
  class(connection.info$n.in) <- "numeric"
  class(connection.info$wsum.out) <- "numeric"
  class(connection.info$wsum.in) <- "numeric"
  
  return(connection.info)
}


# connections.info <- get.first.connection.information(list.links, 
#                                                      node.names = node.names)







## -- FUNCTIONS: find genes upstream and downstream of each node ----------- ## 

find.connected.genes <- function(link.list, node.wanted, collapsed=FALSE){
  ## find the downstream target genes 
  regulator <- node.wanted
  all.targets <- c()
  targets <- "A"
  while(length(targets) > 0){
    # find gene targets of the regulator node 
    targets <- link.list[link.list[,1] %in% regulator,2]
    # only choose targets not in the all.downstream list already
    # avoids issues with cyclic regulation 
    targets <- targets[!(targets %in% all.targets)]
    # add downstream target to the list 
    all.targets <- c(all.targets, unique(targets))
    # make the immediate downstream targets the reguator in the next round 
    regulator <- unique(targets)
  }
  
  ## find the upstream (regulon) genes 
  targets <- node.wanted
  all.regulators <- c()
  regulator <- "A"
  while(length(regulator) > 0){
    # find gene regulons of the target node 
    regulator <- link.list[link.list[,2] %in% targets,1]
    # only choose targets not in the all.regulator list already
    # avoids issues with cyclic regulation 
    regulator <- regulator[!(regulator %in% all.regulators)]
    # add upstream regulators to the list 
    all.regulators <- c(all.regulators, unique(regulator))
    # make the immediate upstream regulators the target in the next round 
    targets <- unique(regulator)
  }
  
  if(length(all.targets)==0){all.targets <- NULL}
  if(length(all.regulators)==0){all.regulators <- NULL}
  
  if(collapsed==TRUE){
    all.targets <- paste(all.targets, collapse = "_")
    all.regulators <- paste(all.regulators, collapse = "_")
  }
  
  return(list(downstream=all.targets, upstream=all.regulators))
}

# connected.genes.to.A <- find.connected.genes(link.list = list.links, 
#                                              node.wanted = "A", collapsed = FALSE)






## -- FUNCTION: find number and name of upstream/downstream genes from all nodes 

extract.network.connections <- function(link.list, node.names){
  
  # set up table holding info
  connected.genes <- matrix(NA, length(node.names), 5,
                            dimnames = list(NULL,
                                            c("node", "upstream", 
                                              "downstream", "n.up","n.down")))
  # add node names 
  connected.genes[,"node"] <- node.names
  
  # run function "find.connected.genes" for each node 
  # store both connected nodes' names and their counts
  for (i in 1:length(node.names)){
    x <- find.connected.genes(link.list = link.list, 
                              node.wanted = node.names[i], 
                              collapsed = FALSE)
    connected.genes[i,"n.up"] <- length(x$upstream)
    connected.genes[i,"n.down"] <- length(x$downstream)
    connected.genes[i,"upstream"] <- paste(x$upstream, collapse = "_")
    connected.genes[i,"downstream"] <- paste(x$downstream, collapse = "_")
  }
  
  connected.genes <- as.data.frame(connected.genes, stringsAsFactors=FALSE)
  connected.genes$n.up <- as.numeric(connected.genes$n.up)
  connected.genes$n.down <- as.numeric(connected.genes$n.down)
  
  return(connected.genes)
  
}

# network.node.connections.stats <- extract.network.connections(link.list = list.links, 
#                                                               node.names = node.names)









## -- FUNCTION: find pairwise common edges --------------------------------- ##

network.pairwise.edge.comparison <- function(sample.list, reference.list, 
                                             node.names, threshold=0){
  n <- length(node.names)
  # extract adjacency matrix from lists 
  sample.matrix <- get.adjmatrix.from.list(link.list = sample.list, 
                                           node.names = node.names)
  ref.matrix <- get.adjmatrix.from.list(link.list = reference.list, 
                                           node.names = node.names)
  
  # turn all entries >threshold into 1
  sample.matrix[sample.matrix > threshold] <- 1
  ref.matrix[ref.matrix > threshold] <- 1
  
  # simliarity table to calculate TP, FP, TN, FN etc. 
  sim.table <- table(sample.matrix, ref.matrix) 
  # obtain the metrics
  TN <- sim.table["0","0"]
  FN <- sim.table["0","1"]
  TP <- sim.table["1","1"]
  FP <- sim.table["1","0"]
  precision <- TP/(TP+FP)
  recall <- TP/(TP+FN)
  TNR <- TN/(TN+FP)
  accuracy <- (TP+TN)/sum(sim.table)
  F1 <- 2*precision*recall/(precision+recall)
  
  # extract the edges names (direction sensitive) from list

  sample.edges <- paste(sample.list$regulator[sample.list$weight>0], 
                        sample.list$target[sample.list$weight>0], 
                        sep = "~")
  ref.edges <- paste(reference.list$regulator[reference.list$weight>0], 
                     reference.list$target[sample.list$weight>0], 
                        sep = "~")
  # find common egdes
  common.edges <- sample.edges[sample.edges %in% ref.edges]
  
  list.common.edges <- matrix(NA,length(common.edges), 2, 
                              dimnames = list(NULL,c("regulator","target")))
  for (i in 1:length(common.edges)){
    list.common.edges[i,] <- unlist(strsplit(common.edges[i], split = "~"))
  }
  
  return(list(common.edges=list.common.edges, 
              metrics=c(TN=TN, FN=FN, TP=TP, FP=FP, 
                        prec=precision, recall=recall, 
                        true.neg.rate=TNR, accuracy=accuracy, 
                        F1=F1)
              )
         )
  
}



















