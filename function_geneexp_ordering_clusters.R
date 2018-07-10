## ================== gene expression ordering source code ================= ##

# order the genes based on the time of peak expression 
# future plan: order genes OR clusters of genes based of peak expression time


temporal.order.genes <- function(x){
  
  peak.time <- apply(x, 1, which.max)
  gene.order <- order(peak.time, decreasing = FALSE)
  gene.order.names <- rownames(x[gene.order,])
  
  return(list(peak.time=peak.time, peak.order=gene.order, 
              names=gene.order.names))
  
}


# gene.clusters <- hclust(dist((expr.selected)))
# plot(gene.clusters)
# gene.clusterCut <- cutree(gene.clusters, 500)
# 
# mean.time.peak.clust <- rep(NA,8)
# for (i in 1:length(mean.time.peak.clust)){
#   clust.members <- names(gene.clusterCut[gene.clusterCut==i])
#   mean.time.peak.clust[i] <- mean(apply(expr.selected[clust.members,],1,which.max))
# }