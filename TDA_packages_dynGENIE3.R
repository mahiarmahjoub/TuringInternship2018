## ==== dynGENIE3 ---------------------------------------------------------- ##

setwd("C://Users/mahia/Google Drive/MPhil/CompBio/Internship/TuringInternship2018/algorithms/dynGENIE3_R_C_wrapper/")
source("dynGENIE3.R")
# install.packages("reshape2")
# install.packages("doRNG")
# install.packages("doParallel")
require(reshape2); require(doRNG); require(doParallel)
require(stringi, quietly = TRUE)

## ==== load expression matrix AND wanted gene list ------------------------ ##
gene.list <- read.csv(file.path(dir,"changing_genes_list_v1_180622.csv"), 
                      header = TRUE, stringsAsFactors = FALSE, sep = ",")
TF.list <- read.csv(file.path(dir,"changing_genes_list_v1_180622_TF.csv"), 
                    header = TRUE, stringsAsFactors = FALSE, sep = ",")
target.list <- read.csv(file.path(dir,"changing_genes_list_v1_180622_targets.csv"), 
                        header = TRUE, stringsAsFactors = FALSE, sep = ",")
expr.all <- read.table(file = file.path(dir, "dawnburst_slcu", 
                                        "rawcounts", 
                                        "dawnburst_masterdf.txt"), 
                       header = TRUE, stringsAsFactors = FALSE
)


## choose the right expression set for analysis 
temp <- "22"; mut <- "Col"
time <- c(-10, 10, 24, 30, 45, 60, 105, 120)
# pick selected experimental conditions 
expr <- expr.all[,grepl(temp,colnames(expr.all)) & grepl(mut,colnames(expr.all))]
colnames(expr) <- as.character(time)
rownames(expr) <- expr.all$tracking_id
expr.selected <- expr[as.matrix(gene.list),]  # pick selected genes 


## ==== setup the data for analysis ---------------------------------------- ##
TS1 <- read.expr.matrix("example_data/time_series_1.txt", form = "rows.are.samples")
TS2 <- read.expr.matrix("example_data/time_series_2.txt", form="rows.are.samples")
time.points <- list(TS1[1,], TS2[2,])
TS.data <- list(TS1[2:nrow(TS1),], TS2[2:nrow(TS2),])


## ==== run dynGENIE3 ------------------------------------------------------ ##
res.dynGENIE3 <- dynGENIE3(TS.data, time.points)

link.list.dynGENIE3 <- get.link.list(res$weight.matrix)

res.dynGENIE <- read.csv(file.path(dir,"TuringInternship2018/dynGENIE_test_links_results.csv"), 
                         header=TRUE, stringsAsFactors = FALSE)



