## ================== finalising gene list for TD analysis ================= ##

# finding the right genes for TD analysis through the following ways: 
# (1) GO to include CC-related and TF-encoding genes 
# (2) remove unchanging and non-expressing genes (Z-score)
# (3) remove functionally irrelevant genes (DE analysis WT vs prr579 mutant)
# (4) ChIPseq PRR5/7/9

dir <- "C:/Users/mahia/Google Drive/MPhil/CompBio/Internship"
# import comprehensive dawn burst expression data 
expr.all <- read.table(file = file.path(dir, "dawnburst_slcu", 
                                       "rawcounts", 
                                       "dawnburst_masterdf.txt"), 
                       header = TRUE, stringsAsFactors = FALSE
                       )
time <- c(-10, 10, 24, 30, 45, 60, 105, 120)


## ---- GO analysis -------------------------------------------------------- ##

# import file containing ALL genes of arabidopsis and their GO annotation
GO2genes <- read.delim(file = file.path(dir,"tair.gaf"), skip = 12, 
                       stringsAsFactors = FALSE, header=FALSE
                       )
# import file og GO terms and functional information & db source (B, F, C)
GO2Info <- read.delim(file = file.path(dir,"GOTERM_info_arabidopsis.txt"), 
                      skip = 11, 
                       stringsAsFactors = FALSE, header = FALSE, 
                      col.names = c("GO","otherGO", "info", "db", "other")
)

# choose GO terms: "transcription factor" from Cellular Components 
TF.GOterms <- grepl("transcription factor", GO2Info$info) & GO2Info$db=="C"
# choose GO terms: "circadian clock/rhythm/behavior" from Biological Processes 
CC.GOterms <- (grepl("circadian clock", GO2Info$info) |
  grepl("circadian rhythm", GO2Info$info) |
  grepl("circadian behavior", GO2Info$info)) & GO2Info$db=="P" 
# put all GO terms from above categories together (union)
all.GOterms <- GO2Info$GO[TF.GOterms | CC.GOterms]

# obtain the corresponding GO terms for each wanted gene
genes.GO.terms <- GO2genes$V5[GO2genes$V5 %in% all.GOterms]
# find genes matching the wanted GO terms (in terms of IDs)
genes.GO.selected <- substring(GO2genes$V11[GO2genes$V5 %in% all.GOterms],1,9)
genes.GO.selected[genes.GO.selected=="UniProtKB"] <- 
  genes.GO.selected[genes.GO.selected=="UniProtKB"]
# obtain the gene names 
genes.GO.selected.name <- GO2genes$V3[GO2genes$V5 %in% all.GOterms]

# put all info together - remove non-AT genes 
# note that some genes were excluded with "At" beginnings or no AT notation
selected.GO.genes.df <- cbind(geneID=genes.GO.selected[grepl("AT",genes.GO.selected)], 
                              gene.name=genes.GO.selected.name[grepl("AT",genes.GO.selected)], 
                              GO=genes.GO.terms[grepl("AT",genes.GO.selected)])
# write.csv(selected.GO.genes.df, file = file.path(dir,"GO_selected_gene.csv"), quote = FALSE,
#             row.names = FALSE)



## ---- ChIPseq ------------------------------------------------------------ ##

chipseq.targets <- read.csv(file = file.path(dir,"boundgenes_500bp_MACS2filt.csv"), 
                            header = TRUE,stringsAsFactors = FALSE)
prr1.targets <- chipseq.targets$PRR1[!is.na(chipseq.targets$PRR1)]
prr5.targets <- chipseq.targets$PRR5[!is.na(chipseq.targets$PRR5)]
prr7.targets <- chipseq.targets$PRR7[!is.na(chipseq.targets$PRR7)]
prr7HA.targets <- chipseq.targets$PRR7HA[!is.na(chipseq.targets$PRR7HA)]




## ---- literature --------------------------------------------------------- ##
lit.genes <- read.csv(file = file.path(dir,"CircadianClock_genelist_literature.csv"), 
                      header=TRUE, stringsAsFactors = FALSE)
lit.genes.ID <- lit.genes$Gene.ID



## ==== put ChIPseq + lit + GO together ==================================== ##
# put all the genes together 
prelim.gene.ID.list <- c(selected.GO.genes.df[,"geneID"], 
                         lit.genes.ID, 
                         prr1.targets, prr5.targets, prr7.targets, prr7HA.targets
                         )
prelim.gene.ID.list <- unique(prelim.gene.ID.list)

# remove spaces in the character strings 
prelim.gene.ID.list <- unique(sapply(prelim.gene.ID.list, function(x){substring(x,1,9)}, 
                              USE.NAMES = FALSE)
                              )




## ---- calculate Z-score -------------------------------------------------- ##
# use the WT data for Col-0 and Ler
Col0.22 <- expr.all[,grepl("Col",colnames(expr.all)) & 
                      grepl("22",colnames(expr.all))]
Col0.27 <- expr.all[,grepl("Col",colnames(expr.all)) & 
                      grepl("27",colnames(expr.all))]
Ler.22 <- expr.all[,grepl("Ler",colnames(expr.all)) & 
                     grepl("22",colnames(expr.all))]
Ler.27 <- expr.all[,grepl("Ler",colnames(expr.all)) & 
                     grepl("27",colnames(expr.all))]

# choose one of the datasets 
selected.expr <- Col0.22
rownames(selected.expr) <- expr.all$tracking_id

# only pick genes selected from ChIPseq, GO and lit
selected.expr <- selected.expr[expr.all$tracking_id %in% prelim.gene.ID.list,]


# count the number of <1 counts per gene (unexpressing genes)
n.zero.expr <- apply(selected.expr, 1, function(x){sum(x<1)})
# calculate the Z-score for each gene and timepoint
selected.expr.Zscore <- (selected.expr - 
                           apply(selected.expr,1,mean))/apply(
                             selected.expr,1,sd)
# count the number of timepoints with NaN per gene (unchanging genes)
n.nan.Zscore <- apply(selected.expr.Zscore, 1, function(x){sum(is.nan(x))})

# NaN Z-score + more than 0 zero counts + low spread
# combination to remove unexpressing or unchanging genes 
changing.genes <- rownames(selected.expr[(n.zero.expr < 3) & (n.nan.Zscore == 0),])
changing.genes.GOinfo <- changing.genes.GOID <- changing.genes.db <- changing.genes

for (i in 1:length(changing.genes)){
  GO.term <- GO2genes$V5[which(GO2genes$V10 == changing.genes[i])]
  GO.info <- GO2Info$info[GO2Info$GO %in% GO.term]
  GO.process <- GO2Info$db[GO2Info$GO %in% GO.term]
  changing.genes.GOID[i] <- paste(GO.term, collapse = ";")
  changing.genes.GOinfo[i] <- paste(GO.info, collapse = ";")
  changing.genes.db[i] <- paste(GO.process, collapse = ";")
}

changing.genes.df <- cbind(GeneID=changing.genes,
                           GOID=changing.genes.GOID ,
                           GOinfo=changing.genes.GOinfo,
                           GOcategory=changing.genes.db
                           )
n.TFs <- sum(grepl("transcription factor", changing.genes.GOinfo))
write.csv(changing.genes.df, 
          file.path(dir,"changing_genes_list_v1_180622.csv"), 
          quote = FALSE, row.names = FALSE)








## ------------------------------------------------------------------------- ##
## ------------------------------------------------------------------------- ##
## ------------------------------------------------------------------------- ##
## ---- DE analysis -------------------------------------------------------- ##
# do DE analysis at each timepoint between WT and mutant PRR579?
library(DESeq2)

# grab Ler and mutant PRR579 counts (only for -10 till 60 mins)
Ler.27 <- expr.all[,grepl("Ler",colnames(expr.all)) & 
                     grepl("27",colnames(expr.all))]
Ler.27 <- Ler.27[,c(-7,-8)]
Ler.prr579.27 <- expr.all[,grepl("prr579",colnames(expr.all)) & 
                     grepl("27",colnames(expr.all))]

# put WT and mutant together in a matrix 
Ler.WT.prr579 <- cbind(Ler.27, Ler.prr579.27)
# standardise column names to reflect WT/mutant and time point 
colnames(Ler.WT.prr579) <- samples <- c("WTm10","WT10","WT24","WT30","WT45","WT60",
                             "MUTm10","MUT10","MUT24","MUT30","MUT45","MUT60"
                             )
# put all time combinations into one matrix
all.times <- cbind(time1=as.character(c(1,2,2,2,2,1)),
                   time2=as.character(c(1,1,2,2,2,2)),
                   time3=as.character(c(2,2,2,2,1,1))
                   )

# assign gene names as row names 
rownames(Ler.WT.prr579) <- expr.all$tracking_id

# specify WT/mutant files and place mutant+WT count data into DESeq object
coldata.mut <- cbind(samples, 
                 condition=c(rep("untreated",6),rep("treated",6))
                 )
dds.mut <- DESeqDataSetFromMatrix(countData = as.matrix(round(Ler.WT.prr579)),
                                  colData = coldata.mut,
                                  design= ~ condition 
                                  )
featureData <- data.frame(gene=rownames(Ler.WT.prr579))
mcols(dds.time) <- DataFrame(mcols(dds.time), featureData)

# find changing genes for the 3 different time structures 
temporal.changing.genes <- c()
for (i in 1:3){
  # place WT count data into DESeq object with the time design matrix 
  coldata.time <- cbind(samples, time=all.times[,i])
  dds.time <- DESeqDataSetFromMatrix(countData = as.matrix(round(Ler.WT.prr579[,1:6])),
                                     colData = coldata.time[1:6,],
                                     design= ~ time
                                     )
  mcols(dds.mut) <- DataFrame(mcols(dds.mut), featureData)

  dds <- dds.time
  # filter out unexpressed genes 
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]

  # set conditions for DE pairwise comparison 
  #dds$condition <- factor(dds$condition, levels = c("untreated","treated"))
  dds$time <- factor(dds$time, levels = as.character(c(2,1)))

  # run DESeq2
  dds <- DESeq(dds)
  resultsNames(dds) # lists the coefficients
  res <- results(dds, name=resultsNames(dds)[2], alpha = 0.05)
  resOrdered <- res[order(res$pvalue),]


  temporal.changing.genes <- c(temporal.changing.genes, 
                               rownames(resOrdered)[resOrdered$pvalue < 0.01])

}

temporal.changing.genes <- unique(temporal.changing.genes)

write.csv(temporal.changing.genes, 
          file.path(dir,"DE_temporal_gene_list.csv"), 
          quote = FALSE, row.names = FALSE)





