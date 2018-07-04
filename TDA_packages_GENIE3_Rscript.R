## ==== TDCor -------------------------------------------------------------- ##

setwd("/home/mahiar.mahjoub/TD_algorithms")
require(stringi, quietly = TRUE)

## ==== load expression matrix AND wanted gene list ------------------------ ##
gene.list <- read.csv("changing_genes_list_v1_180622.csv", 
                      header = TRUE, stringsAsFactors = FALSE, sep = ",")
TF.list <- read.csv("changing_genes_list_v1_180622_TF.csv", 
                    header = TRUE, stringsAsFactors = FALSE, sep = ",")
target.list <- read.csv("changing_genes_list_v1_180622_targets.csv", 
                        header = TRUE, stringsAsFactors = FALSE, sep = ",")
expr.all <- read.table("dawnburst_masterdf.txt", 
                       header = TRUE, stringsAsFactors = FALSE)


## choose the right expression set for analysis 
temp <- "22"; mut <- "Col"
time <- c(-10, 10, 24, 30, 45, 60, 105, 120)
# pick selected experimental conditions 
expr <- expr.all[,grepl(temp,colnames(expr.all)) & grepl(mut,colnames(expr.all))]
colnames(expr) <- as.character(time)
rownames(expr) <- expr.all$tracking_id
expr.selected <- expr[as.matrix(gene.list),]  # pick selected genes 


## ------------------------------------------------------------------------- ##
require("GENIE3", quietly = TRUE)
res.GENIE3 <- GENIE3(as.matrix(expr.selected[,1:3]))
#weightMat <- GENIE3(exprMatr, treeMethod="ET", K=7, nTrees=50)
res.linkList.GENIE <- getLinkList(res.GENIE3)

write.csv(res.linkList.GENIE, "GENIE3_links_180626.csv", quote = FALSE, 
          row.names = FALSE)