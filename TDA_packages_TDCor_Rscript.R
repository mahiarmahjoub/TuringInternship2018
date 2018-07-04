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
require("TDCor", quietly = TRUE)


TPI10 <- CalculateTPI(dataset=as.matrix(expr.selected), 
                      l_genes=rownames(expr.selected), 
                      l_prior=rep(2,nrow(expr.selected)),
                      times=time, time_step=1, N=10000, ks_int=c(0.5,3),
                      kd_int=c(0.5,3), delta_int=c(0.5,3), noise=0.1, delay=3)
DPI15 <- CalculateDPI(dataset=as.matrix(expr.selected), 
                      l_genes=rownames(expr.selected), 
                      l_prior=rep(2,nrow(expr.selected)),
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
poutfile_name="TDCor_links_180626.txt"


tdcor_out= TDCOR(dataset=as.matrix(expr.selected), 
                 l_genes=rownames(expr.selected), 
                 l_names=rownames(expr.selected),
                 n0=pn0, n1=pn1, 
                 l_prior=rep(2,nrow(expr.selected)), 
                 thr_ind1=pthr_ind1, 
                 thr_ind2=pthr_ind2, regmax=pregmax, thr_cor=pthr_cor,
                 delayspan=pdelayspan, delaymax=pdelaymax, delaymin=pdelaymin,
                 delay=pdelay, thrpTPI=pthrpTPI, thrpDPI=pthrpDPI, TPI=pTPI,
                 DPI=pDPI, thr_isr=pthr_isr, time_step=ptime_step, 
                 thr_overlap=pthr_overlap, tol=ptol, MinProp=pMinProp,
                 MinTarNumber=pMinTarNumber, outfile_name=poutfile_name)


