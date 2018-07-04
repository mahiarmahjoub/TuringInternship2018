## ================== gene expression ordering source code ================= ##

# order the genes based on the time of peak expression 
# future plan: order genes OR clusters of genes based of peak expression time


temporal.order.genes <- function(x, order.fun=which.max){
  
  gene.order <- order(apply(x, 1, order.fun))
  
  return(list(peak.order=gene.order))
  
}
