args <- commandArgs(trailingOnly=TRUE)
if(length(args)==0){stop("Please input mutant type and temperature.", 
			 call.=FALSE)}

## ==== TD-ARACNE ---------------------------------------------------------- ##

setwd("/home/mahiar.mahjoub/TD_algorithms")
source("dynGENIE3.R")
require(TDARACNE, quietly = TRUE)
require(stringi, quietly = TRUE)

## ==== load expression matrix AND wanted gene list ------------------------ ##
gene.list <- read.csv("changing_genes_list_v3_180702.csv", 
                      header = TRUE, stringsAsFactors = FALSE, sep = ",")
TF.list <- read.csv("changing_genes_list_v3_180702_TF.csv", 
                    header = TRUE, stringsAsFactors = FALSE, sep = ",")
target.list <- read.csv("changing_genes_list_v3_180702_targets.csv", 
                        header = TRUE, stringsAsFactors = FALSE, sep = ",")
expr.all <- read.table("dawnburst_masterdf.txt", 
                       header = TRUE, stringsAsFactors = FALSE)


## choose the right expression set for analysis 
temp <- args[2]; mut <- args[1]
#temp <- "22"; mut <- "Col"
time <- c(-10, 10, 24, 30, 45, 60, 105, 120)
# pick selected experimental conditions 
expr <- expr.all[,grepl(temp,colnames(expr.all)) & grepl(mut,colnames(expr.all))]
if(mut=="Ler" & temp=="22"){time <- time[-2]; expr <- expr[,-ncol(expr)]}
if(mut=="prr579" & temp=="22"){expr <- expr[,-ncol(expr)]}
if(mut=="prr579" & temp=="27"){time <- time[1:6]}
if(mut=="phyA" & temp=="22"){expr <- expr[,1:8]}
if(mut=="phyA" & temp=="27"){expr <- expr[,1:8]}


colnames(expr) <- as.character(time)
rownames(expr) <- expr.all$tracking_id
expr.selected <- expr[as.matrix(gene.list),]  # pick selected genes 


## ------------------------------------------------------------------------- ##
require("Biobase", quietly = TRUE)

expr.selected.ESobject <- ExpressionSet(assayData = 
                                              as.matrix(expr.selected)
                                            )

res.TDARACNE <- TDARACNE(expr.selected.ESobject, 11, adj = TRUE)


write.csv(res.TDARACNE, 
	  paste0("TDARACNE_adjmatrix_180702_v3_",mut,"_",temp,".csv"), 
	  quote = FALSE, 
          row.names = TRUE)
write.csv(get.link.list(res.TDARACNE, threshold = 0.01), 
          paste0("TDARACNE_list_180702_v3_",mut,"_",temp,".csv"),
	  quote = FALSE, 
          row.names = FALSE)
