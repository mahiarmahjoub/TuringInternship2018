## ================== gene expression heatmap source code ================== ##

# code based on 'gplot' or 'ComplexHeatmap' (Bioconductor) packages 
# library(c('gplots', 'ComplexHeatmaps'))

require(gplots, quietly = TRUE)
require(RColorBrewer, quietly = TRUE)

# MAIN input variables: expression matrix, gene.list=NULL
# other input variables: scaled=TRUE/FALSE, order=c('none','peak.time','given')

expression.heatmap <- function(x, y, gene.list=NULL, gene.order=NULL,
                               order=c('none','given','peak.time'),
                               scaled=c('none','col','row'), 
                               value='raw', tile.col=NULL, ncol=2,
                               col.clust=NA, row.clust=NA, dend='none',
                               ...){
  
  # if gene.list is not provided, rownames are taken in instead 
  if(is.null(gene.list)){gene.list <- rownames(x)}
  # generate log expression matrix if requested 
  if(value=='log2'){x.log <- as.matrix(log2(x + 1))}
  # determine the order of gene list 
  if(order=='peak.time'){
    gene.order <- order(apply(x,1,which.max))
  } else if (order=='given'){
    gene.order <- gene.order
  } else if (order=='none') {
    gene.order <- 1:nrow(x)
  }
  # tile colours
  if (is.null(tile.col)){tile.col <- c("red","blue")}

  
  # plot - column (gene) scaled + non-log
  heatmap.2(t(x[gene.order,]), 
            Rowv = row.clust, Colv=col.clust, dendrogram = dend,
            trace = 'none', scale = scaled,
            col=colorRampPalette(colors = tile.col)(ncol),
            labCol = gene.list[gene.order], 
            labRow = as.character(y), key.title = NA, key.ylab = NA,
            cexCol = 0.55, srtCol = 45, ... 
  )

  
  
}


# example of use:
# expression.heatmap(x=exp, y=time, order = 'given', 
#                    gene.order = exp.temporal.max.sorted, scaled = 'col',
#                    ncol=4, gene.list = genes.selected.name, main="lsfgb", 
#                    ylab="oaerty")

