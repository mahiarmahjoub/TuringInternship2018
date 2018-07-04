args <- commandArgs(trailingOnly=TRUE)
if(length(args)==0){stop("Please specify mutant type and temperature.", 
		         call.=FALSE)}

## ==== dynGENIE3 ---------------------------------------------------------- ##

setwd("/home/mahiar.mahjoub/TD_algorithms")
source("dynGENIE3.R")
require(reshape2); require(doRNG); require(doParallel)
require(stringi, quietly = TRUE)

## ==== load expression matrix AND wanted gene list ------------------------ ##
gene.list <- read.csv("changing_genes_list_v3_180702.csv", 
                      header = TRUE, stringsAsFactors = FALSE, sep = ",")
TF.list <- read.csv("changing_genes_list_v3_180702_TF.csv", 
                    header = TRUE, stringsAsFactors = FALSE, sep = ",")
TF.list <- TF.list$x
target.list <- read.csv("changing_genes_list_v3_180702_targets.csv", 
                        header = TRUE, stringsAsFactors = FALSE, sep = ",")
expr.all <- read.table("dawnburst_masterdf.txt", 
                       header = TRUE, stringsAsFactors = FALSE)


## choose the right expression set for analysis 
#temp <- "22"; mut <- "Col"
temp <- args[2]; mut <- args[1]
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
expr.selected <- as.matrix(expr[gene.list$x,])  # pick selected genes

# splining expression data
require("GeneReg", quietly = TRUE)
expr.bspline.data <- ts.bspline(expr.selected, 
                                ts.point = as.numeric(colnames(expr.selected)), 
                                data.predict = 50)


# remove genes with 0 expression time points  
#if(mut=="Ler" & temp=="27"){
expr.bspline.data <- expr.bspline.data[apply(expr.bspline.data,1,
					function(x){sum(x==0)})==0,]
TF.list <- TF.list[TF.list %in% rownames(expr.bspline.data)]
#}


## setup inputs of dynGENIE3 
TS.data <- list(as.matrix(expr.bspline.data))
time.points <- list(as.numeric(colnames(expr.bspline.data)))


## ==== run dynGENIE3 ------------------------------------------------------ ##
res.dynGENIE3 <- dynGENIE3(TS.data, time.points, 
		           regulators = TF.list
			  )

link.list.dynGENIE3 <- get.link.list(res.dynGENIE3$weight.matrix)

write.csv(link.list.dynGENIE3, 
	  paste0("dynGENIE_links_180702_v3_",mut,"_",temp,"_splined.csv"), 
          quote = FALSE, row.names = FALSE)



