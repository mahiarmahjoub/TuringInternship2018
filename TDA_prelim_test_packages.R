## ================== priliminary testing out various TDAs ================= ##

# source("https://bioconductor.org/biocLite.R")
# install.packages("GeneReg")
# install.packages("TDCor")
# biocLite("TDARACNE")
# biocLite("networkBMA")
# install.packages("G1DBN")
# install.packages("iRafNet")
# biocLite("GENIE3")

dir <- "C://Users/mahia/Google Drive/MPhil/CompBio/Internship"


## ==== GeneReg
require("GeneReg", quietly = TRUE)
data(wt.expr.data)
wt.bspline.data <- ts.bspline(wt.expr.data, 
                              ts.point = as.numeric(colnames(wt.expr.data)), 
                              data.predict = 100)
data("tf.list")
wt.models <- timedelay.lm.batch(bspline.data=wt.bspline.data[1:100,],
                                expr.data=wt.expr.data[1:100,], 
                                regulator.list=tf.list,
                                target.list=rownames(wt.bspline.data[1:100,]),
                                single.adj.r.squared=0.8, 
                                multiple.adj.r.squared=0.9,
                                maxdelay=ncol(wt.bspline.data)*0.1, 
                                min.coef=0.25, max.coef=4,
                                output=T, topdf=T, 
                                xlab='Time point (lifeline)',
                                ylab='Relative expression level (in log ratio)'
                                )

res.genereg <- wt.models
#read.csv(file.path(dir, "TuringInternship2018/algorithms/genereg_result_test.csv"), 
#                 header = TRUE, stringsAsFactors = FALSE)

#plot.GeneReg(as.matrix(wt.models),vertex.size=2,layout=layout.fruchterman.reingold)





## ==== TDCor
require("TDCor", quietly = TRUE)
data(LR_dataset)
data(l_genes)
data(l_names)
data(l_prior)
data(times)


TPI10 <- CalculateTPI(dataset=LR_dataset, l_genes=l_genes, l_prior=l_prior,
                      times=times, time_step=1, N=10000, ks_int=c(0.5,3),
                      kd_int=c(0.5,3), delta_int=c(0.5,3), noise=0.1, delay=3)
DPI15 <- CalculateDPI(dataset=LR_dataset, l_genes=l_genes, l_prior=l_prior,
                     times=times, time_step=1, N=10000, ks_int=c(0.5,3),
                     kd_int=c(0.5,3), delta_int=c(0.5,3), noise=0.15, delay=3)
TPI10 <- UpdateTPI(TPI10,LR_dataset,l_genes,l_prior)
DPI15 <- UpdateDPI(DPI15,LR_dataset,l_genes,l_prior)


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


tdcor_out= TDCOR(dataset=LR_dataset, l_genes=l_genes, l_names=l_names,
                 n0=pn0, n1=pn1, l_prior=l_prior, thr_ind1=pthr_ind1, 
                 thr_ind2=pthr_ind2, regmax=pregmax, thr_cor=pthr_cor,
                 delayspan=pdelayspan, delaymax=pdelaymax, delaymin=pdelaymin,
                 delay=pdelay, thrpTPI=pthrpTPI, thrpDPI=pthrpDPI, TPI=pTPI,
                 DPI=pDPI, thr_isr=pthr_isr, time_step=ptime_step, 
                 thr_overlap=pthr_overlap, tol=ptol, MinProp=pMinProp,
                 MinTarNumber=pMinTarNumber, outfile_name=poutfile_name)

res.tdcor <- read.delim(file.path(dir, "TuringInternship2018/TDCor_output.txt"), 
                          header = FALSE, stringsAsFactors = FALSE)




## ==== TD-ARACNE ---------------------------------------------------------- ##
require("TDARACNE", quietly = TRUE)
data(dataYeast)
data(threshYeast)
data(dataSOSmean)
data(threshSOSmean)
data(dataIRMAon)
data(threshIRMAon)

res.TDARACNE <- TDARACNE(dataSOSmean, 11)
plot(res.TDARACNE, attrs = list(node = list(fillcolor = "lightblue"),
                            edge = list(arrowsize=0.5)))




## ==== GENIE3 ------------------------------------------------------------- ##
require("GENIE3", quietly = TRUE)
exprMatr <- matrix(sample(1:10, 100, replace=TRUE), nrow=20)
rownames(exprMatr) <- paste("Gene", 1:20, sep="")
colnames(exprMatr) <- paste("Sample", 1:5, sep="")
weightMat <- GENIE3(exprMatr)
weightMat <- GENIE3(exprMatr, treeMethod="ET", K=7, nTrees=50)
res.linkList.GENIE <- getLinkList(weightMat)






## ==== dynGENIE3 ---------------------------------------------------------- ##
setwd("C://Users/mahia/Google Drive/MPhil/CompBio/Internship/TuringInternship2018/algorithms/dynGENIE3_R_C_wrapper/")
source("dynGENIE3.R")
# install.packages("reshape2")
# install.packages("doRNG")
# install.packages("doParallel")
require(reshape2); require(doRNG); require(doParallel)
TS1 <- read.expr.matrix("example_data/time_series_1.txt", form = "rows.are.samples")
TS2 <- read.expr.matrix("example_data/time_series_2.txt", form = "rows.are.samples")
TS3 <- read.expr.matrix("example_data/time_series_3.txt", form = "rows.are.samples")

time.points <- list(TS1[1,], TS2[1,], TS3[1,])
TS.data <- list(TS1[2:nrow(TS1),], TS2[2:nrow(TS2),], TS3[2:nrow(TS3),])

res <- dynGENIE3(TS.data, time.points)

link.list <- get.link.list(res$weight.matrix)

res.dynGENIE <- read.csv(file.path(dir,"TuringInternship2018/dynGENIE_test_links_results.csv"), 
                          header=TRUE, stringsAsFactors = FALSE)

