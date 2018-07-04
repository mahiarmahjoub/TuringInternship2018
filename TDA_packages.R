## ================== priliminary testing out various TDAs ================= ##


dir <- "C://Users/mahia/Google Drive/MPhil/CompBio/Internship"
require(stringi, quietly = TRUE)
source(file.path(dir,"TuringInternship2018/function_plot_network.R"))


## ==== load expression matrix AND wanted gene list 
gene.list <- read.csv(file.path(dir,"changing_genes_list_v1_180622.csv"), 
                      header = TRUE, stringsAsFactors = FALSE, sep = ",")
TF.list <- read.csv(file.path(dir,"changing_genes_list_v1_180622_TF.csv"), 
                    header = TRUE, stringsAsFactors = FALSE, sep = ",")
target.list <- read.csv(file.path(dir,"changing_genes_list_v1_180622_targets.csv"), 
                    header = TRUE, stringsAsFactors = FALSE, sep = ",")
gene.alias <- read.delim(file.path(dir,"TAIR_genealiases.txt"), header=TRUE, 
                         stringsAsFactors = TRUE)
expr.all <- read.table(file = file.path(dir, "dawnburst_slcu", 
                                        "rawcounts", 
                                        "dawnburst_masterdf.txt"), 
                       header = TRUE, stringsAsFactors = FALSE
                       )


## choose the right expression set for analysis 
temp <- "27"; mut <- "Col"
time <- c(-10, 10, 24, 30, 45, 60, 105, 120)
# pick selected experimental conditions 
expr <- expr.all[,grepl(temp,colnames(expr.all)) & grepl(mut,colnames(expr.all))]
colnames(expr) <- as.character(time)
rownames(expr) <- expr.all$tracking_id
expr.selected <- expr[as.matrix(gene.list),]  # pick selected genes 
expr.selected.log <- log2(expr.selected+1)  # log2 normalise


## ==== GeneReg
require("GeneReg", quietly = TRUE)
expr.bspline.data <- ts.bspline(expr.selected, 
                              ts.point = as.numeric(colnames(expr.selected)), 
                              data.predict = 50)
setwd(file.path(dir,"genereg_output"))
res.genereg <- timedelay.lm.batch(bspline.data=expr.bspline.data[1:100,],
                                expr.data=expr.selected[1:100,], 
                                regulator.list=as.matrix(TF.list),
                                target.list=rownames(expr.bspline.data[1:100,]),
                                single.adj.r.squared=0.9, 
                                multiple.adj.r.squared=0.9,
                                maxdelay=ncol(expr.bspline.data)*0.1, 
                                min.coef=0.9, max.coef=4,
                                output=T, topdf=T, 
                                xlab='Time point (lifeline)',
                                ylab='Relative expression level (in log ratio)'
                                )
setwd("..")


# setup graph
graph.genereg <- construct.network.graph(x = res.genereg[,1:3], 
                                         input.mode = "list", 
                                         edge.mode = "directed")

# extract gene names based on ID
for (i in 1:length(nodes(graph.genereg))){
  if (any(gene.alias$locus_name == nodes(graph.genereg)[i])){
    nodes(graph.genereg)[i] <- as.character(gene.alias$symbol[which(gene.alias$locus_name == 
                                                         nodes(graph.genereg)[i])][1])
  }
}

# plot graph
plot(graph.genereg, attrs = list(node = list(fillcolor = "lightblue"),
                           edge = list(arrowsize=0.5)))



## ==== TD-ARACNE ---------------------------------------------------------- ##
require("TDARACNE", quietly = TRUE)
require("Biobase")
data("dataSOSmean")

expr.bspline.data.ESobject <- ExpressionSet(assayData = as.matrix(
  expr.selected[1:100,])
  )
res.TDARACNE <- TDARACNE(expr.bspline.data.ESobject, 11)
res.TDARACNE.test <- TDARACNE(dataSOSmean, 11, adj = TRUE)
plot(res.TDARACNE, attrs = list(node = list(fillcolor = "lightblue"),
                                edge = list(arrowsize=0.5)))




## ==== GENIE3 ------------------------------------------------------------- ##
require("GENIE3", quietly = TRUE)
res.GENIE3 <- GENIE3(as.matrix(expr.selected), 
                     regulators = as.character(as.matrix(TF.list)),
                     targets = rownames(expr.selected))
#weightMat <- GENIE3(exprMatr, treeMethod="ET", K=7, nTrees=50)
res.linkList.GENIE <- getLinkList(res.GENIE3, threshold = 0.013)

# setup graph
graph.GENIE <- construct.network.graph(x = res.linkList.GENIE, 
                                         input.mode = "list", 
                                         edge.mode = "directed")

# extract gene names based on ID
for (i in 1:length(nodes(graph.GENIE))){
  if (any(gene.alias$locus_name == nodes(graph.GENIE)[i])){
    nodes(graph.GENIE)[i] <- as.character(gene.alias$symbol[which(gene.alias$locus_name == 
                                                                      nodes(graph.GENIE)[i])][1])
  }
}

# plot graph
plot(graph.GENIE, attrs = list(node = list(fillcolor = "lightblue"),
                                 edge = list(arrowsize=0.5)))


## ==== TDCor
require("TDCor", quietly = TRUE)
data(LR_dataset)
data(l_genes)
data(l_names)
data(l_prior)
data(times)


TPI10 <- CalculateTPI(dataset=as.matrix(expr.selected[1:50,]), 
                      l_genes=rownames(expr.selected)[1:50], l_prior=rep(2,50),
                      times=time, time_step=1, N=10000, ks_int=c(0.5,3),
                      kd_int=c(0.5,3), delta_int=c(0.5,3), noise=0.1, delay=3)
DPI15 <- CalculateDPI(dataset=as.matrix(expr.selected[1:50,]), 
                      l_genes=rownames(expr.selected)[1:50], l_prior=rep(2,50),
                     times=time, time_step=1, N=10000, ks_int=c(0.5,3),
                     kd_int=c(0.5,3), delta_int=c(0.5,3), noise=0.15, delay=3)

ptime_step=1
ptol=0.13
pdelayspan=12
pthr_cor=c(0.65,0.8)
pdelaymax=c(2.5,3.5)
pdelaymin=0
pdelay=3
pthrpTPI=c(0.55,0.8)
pthrpDPI=c(0.65,0.8)
pthr_overlap=c(0.4,0.6)
pthr_ind1=0.65
pthr_ind2=3.5
pn0=1000
pn1=10
pregmax=5
pthr_isr=c(4,6)
pTPI=TPI10
pDPI=DPI15
pMinTarNumber=5
pMinProp=0.6
poutfile_name="TDCor_output.txt"


tdcor_out= TDCOR(dataset=as.matrix(expr.selected[1:50,]), 
                 l_genes=rownames(expr.selected)[1:50], l_names=rep(2,50),
                 n0=pn0, n1=pn1, l_prior=l_prior, thr_ind1=pthr_ind1, 
                 thr_ind2=pthr_ind2, regmax=pregmax, thr_cor=pthr_cor,
                 delayspan=pdelayspan, delaymax=pdelaymax, delaymin=pdelaymin,
                 delay=pdelay, thrpTPI=pthrpTPI, thrpDPI=pthrpDPI, TPI=pTPI,
                 DPI=pDPI, thr_isr=pthr_isr, time_step=ptime_step, 
                 thr_overlap=pthr_overlap, tol=ptol, MinProp=pMinProp,
                 MinTarNumber=pMinTarNumber, outfile_name=poutfile_name)


res.tdcor <- read.delim(file.path(dir, "TuringInternship2018/TDCor_output.txt"), 
                          header = FALSE, stringsAsFactors = FALSE)












