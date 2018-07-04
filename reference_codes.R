## ============================ reference codes ============================ ##



## PCA analysis 
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r
# -practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/
library(factoextra)
data(decathlon2)
decathlon2.active <- decathlon2[1:23, 1:10]

res.pca <- prcomp(decathlon2.active, scale = TRUE)
fviz_eig(res.pca)
fviz_pca_ind(res.pca, col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)     # Avoid text overlapping

# Eigenvalues
eig.val <- get_eigenvalue(res.pca)
eig.val

# Results for Variables
res.var <- get_pca_var(res.pca)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 
# Results for individuals
res.ind <- get_pca_ind(res.pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 


# Data for the supplementary individuals - for prediction
ind.sup <- decathlon2[24:27, 1:10]
ind.sup.coord <- predict(res.pca, newdata = ind.sup)
# Plot of active individuals
p <- fviz_pca_ind(res.pca, repel = TRUE)
# Add supplementary individuals
fviz_add(p, ind.sup.coord, color ="blue")







## heirarchical clustering 
clusters <- hclust(dist(iris[, 3:4]))
plot(clusters)
clusterCut <- cutree(clusters, 4)
table(clusterCut, iris$Species)

clusters <- hclust(dist(iris[, 3:4]), method = 'average')
plot(clusters)
clusterCut <- cutree(clusters, 3)
table(clusterCut, iris$Species)

#install.packages("rpuHclust") # supposed to be faster than hclust
#library("rpuHclust")



## h-clustering of genes 
gene.clusters <- hclust(dist((exp)))
plot(gene.clusters)
gene.clusterCut <- cutree(gene.clusters, 9)
names(gene.clusterCut) <-  genes.selected.name
gene.clusterCut.ordered <- sort(gene.clusterCut)
col.cluster <- rep('black', length(gene.clusterCut))
col.cluster[gene.clusterCut.ordered==2] <- 'red'
col.cluster[gene.clusterCut.ordered==3] <- 'blue'
col.cluster[gene.clusterCut.ordered==4] <- 'green'
col.cluster[gene.clusterCut.ordered==5] <- 'magenta'
col.cluster[gene.clusterCut.ordered==6] <- 'orange'
col.cluster[gene.clusterCut.ordered==7] <- 'grey'
col.cluster[gene.clusterCut.ordered==8] <- 'brown'
expression.heatmap(x = exp, y = time, gene.list = names(gene.clusterCut.ordered), 
                   scaled = 'col', 
                   order = 'none', ncol = 5, colCol=col.cluster)




### ----------- TopGo ------------------------------------------------------ ##


# library(topGO, quietly = TRUE)
# library(ALL)
# data(geneList)
# data(ALL)
# affyLib <- paste(annotation(ALL), "db", sep = ".")
# library(package = affyLib, character.only = TRUE)
# sum(topDiffGenes(geneList))
# sampleGOdata <- new("topGOdata",
#                     description = "Simple session", ontology = "BP",
#                     allGenes = geneList, geneSel = topDiffGenes,
#                     nodeSize = 10,
#                     annot = annFUN.db, affyLib = affyLib)
# resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
# resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
# resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")
# allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
#                    classicKS = resultKS, elimKS = resultKS.elim,
#                    orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)
# pValue.classic <- score(resultKS)
# pValue.elim <- score(resultKS.elim)[names(pValue.classic)]
# gstat <- termStat(sampleGOdata, names(pValue.classic))
# gSize <- gstat$Annotated / max(gstat$Annotated) * 4
# colMap <- function(x) {
#   .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
#   return(.col[match(1:length(x), order(x))])
# }
# gCol <- colMap(gstat$Significant)
# plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim",
#        pch = 19, cex = gSize, col = gCol)
# sel.go <- names(pValue.classic)[pValue.elim < pValue.classic]
# cbind(termStat(sampleGOdata, sel.go),
#       elim = pValue.elim[sel.go],
#       classic = pValue.classic[sel.go])
# showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = 'all')
# 
