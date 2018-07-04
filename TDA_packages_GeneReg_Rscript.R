args <- commandArgs(trailingOnly=TRUE)
if(length(args)==0){stop("Please input mutant type and temperature.", 
			 call.=FALSE)}

## ==== GeneReg ------------------------------------------------------------ ##

new.folder <- paste0(args[1], "_",args[2],"_",args[3])
setwd("/home/mahiar.mahjoub/TD_algorithms")
require(stringi, quietly = TRUE)

print(new.folder)

## ==== load expression matrix AND wanted gene list 
gene.list <- read.csv("changing_genes_list_v3_180702.csv", 
                      header = TRUE, stringsAsFactors = FALSE, sep = ",")
TF.list <- read.csv("changing_genes_list_v3_180702_TF.csv", 
                    header = TRUE, stringsAsFactors = FALSE, sep = ",")
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
expr.selected <- expr[as.matrix(gene.list),]  # pick selected genes 






## ==== GeneReg
require("GeneReg", quietly = TRUE)
expr.bspline.data <- ts.bspline(expr.selected, 
                                ts.point = as.numeric(colnames(expr.selected)), 
                                data.predict = 50)
setwd(file.path("/home/mahiar.mahjoub/TD_algorithms/genereg_output",new.folder))
res.genereg <- timedelay.lm.batch(bspline.data=expr.bspline.data,
                                  expr.data=expr.selected, 
                                  regulator.list=as.character(as.matrix(TF.list)),
                                  target.list=rownames(expr.bspline.data),
                                  single.adj.r.squared=0.9, 
                                  multiple.adj.r.squared=0.9,
                                  maxdelay=ncol(expr.bspline.data)*0.1, 
                                  min.coef=0.5, max.coef=4,
                                  output=T, topdf=T, 
                                  xlab='Time point (lifeline)',
                                  ylab='Relative expression level (in log ratio)'
                                  )
setwd("../..")



write.csv(res.genereg, paste0("geneReg_links_180702_",mut,"_",temp,".csv"), 
	  row.names = TRUE, quote = FALSE)
