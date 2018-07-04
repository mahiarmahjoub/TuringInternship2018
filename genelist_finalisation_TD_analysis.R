## ================== finalising gene list for TD analysis ================= ##

# finding the right genes for TD analysis through the following ways: 
# (1) GO to include CC-related and TF-encoding genes 
# (2) remove unchanging and non-expressing genes (Z-score)
# (3) remove functionally irrelevant genes (DE analysis WT vs prr579 mutant)
# (4) ChIPseq PRR5/7/9

dir <- "C:/Users/mahia/Google Drive/MPhil/CompBio/Internship"
require(VennDiagram)
# import comprehensive dawn burst expression data 
expr.all <- read.table(file = file.path(dir, "dawnburst_slcu", 
                                       "rawcounts", 
                                       "dawnburst_masterdf.txt"), 
                       header = TRUE, stringsAsFactors = FALSE
                       )
time <- c(-10, 10, 24, 30, 45, 60, 105, 120)


## ==== GO analysis ======================================================== ##

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


## choose GO terms: "transcription factor" from Cellular Components 
# TF.GOterms <- grepl("transcription factor", GO2Info$info) & GO2Info$db=="C"
# TF.GOterms <- grepl("DNA binding", GO2Info$info) & 
#   grepl("transcription",GO2Info$info) & GO2Info$db=="F"
TF.GOterms <- GO2Info$GO[grepl("DNA binding", GO2Info$info) & 
  grepl("transcription",GO2Info$info) & 
    (GO2Info$db=="F" | GO2Info$db=="P")]

## choose GO terms: "light signaling" from Biological Processes
light.GOterms <- GO2Info$GO[grepl("light", GO2Info$info) & 
                              (GO2Info$db=="P")]

## choose GO terms: "abiotic response" from Biological Processes
stress.GOterms <- GO2Info$GO[grepl("response to stress", 
                                   GO2Info$info) & (GO2Info$db=="P")]

## choose GO terms: "circadian clock/rhythm/behavior" from Biological Processes
CC.GOterms <- GO2Info$GO[(grepl("circadian clock", GO2Info$info) |
  grepl("circadian rhythm", GO2Info$info) |
  grepl("circadian behavior", GO2Info$info)) & GO2Info$db=="P" ]

## obtain genes for each GO term 
# remove non-AT genes 
# note that some genes were excluded with "At" beginnings or no AT notation
TF.GO.genes <- substring(GO2genes$V10[GO2genes$V5 %in% 
                                        TF.GOterms],1,9)
TF.GO.genes <- unique(TF.GO.genes[grepl("^AT", TF.GO.genes)])
light.GO.genes <- substring(GO2genes$V10[GO2genes$V5 %in% 
                                           light.GOterms],1,9)
light.GO.genes <- unique(light.GO.genes[grepl("^AT", light.GO.genes)])
stress.GO.genes <- substring(GO2genes$V10[GO2genes$V5 %in% 
                                            stress.GOterms],1,9)
stress.GO.genes <- unique(stress.GO.genes[grepl("^AT", stress.GO.genes)])
CC.GO.genes <- substring(GO2genes$V10[GO2genes$V5 %in% 
                                        CC.GOterms],1,9)
CC.GO.genes <- unique(CC.GO.genes[grepl("^AT", CC.GO.genes)])

## ---- plot Venn Diagram of overlapping genes ----
# venn.diagram(x = list(TF=TF.GO.genes, light=light.GO.genes, 
#                       stress=stress.GO.genes, CC=CC.GO.genes), 
#              filename = file.path(dir, "GO_genes_Venn_v3.tiff"), 
#              fill=c("red","blue","green","orange"), 
#              main = "GO selected genes", main.fontfamily = "arial", 
#              main.cex = 2, fontfamily = "arial")


## put all GO genes from different categories together 
all.GO.genes <- unique(c(TF.GO.genes, light.GO.genes, 
                         stress.GO.genes, CC.GO.genes))
## extract GO term and gene name for each gene ID
all.GO.genes.name <- all.GO.genes.GOterm <- all.GO.genes
for (i in 1:length(all.GO.genes)){
  if(is.null(which(GO2genes$V10 == all.GO.genes[i]))){next}
  all.GO.genes.name[i] <- GO2genes$V3[which(GO2genes$V10 == all.GO.genes[i])][1]
  all.GO.genes.GOterm[i] <- paste(unique(GO2genes$V5[which(GO2genes$V10 == all.GO.genes[i])]), 
                                  collapse = "_")
}



selected.GO.genes.df <- cbind(geneID=all.GO.genes,
                              gene.name=all.GO.genes.name, 
                              GO=all.GO.genes.GOterm)
# write.csv(selected.GO.genes.df,
#           file = file.path(dir,"GO_selected_gene_180702_v3.csv"),
#           quote = FALSE,
#           row.names = FALSE)



## ==== ChIPseq ============================================================ ##

chipseq.targets <- read.csv(file = file.path(dir,"boundgenes_500bp_MACS2filt.csv"), 
                            header = TRUE,stringsAsFactors = FALSE)
prr1.targets <- chipseq.targets$PRR1[!is.na(chipseq.targets$PRR1)]
prr5.targets <- chipseq.targets$PRR5[!is.na(chipseq.targets$PRR5)]
prr7.targets <- chipseq.targets$PRR7[!is.na(chipseq.targets$PRR7)]
prr7HA.targets <- chipseq.targets$PRR7HA[!is.na(chipseq.targets$PRR7HA)]

## ---- Venn Diagram of ChIPseq targets ----
# venn.diagram(x = list(PRR1=prr1.targets, PRR5=prr5.targets, 
#                       PRR7=prr7.targets, PRR7HA=prr7HA.targets), 
#              filename = file.path(dir, "ChIPseq_genes_Venn_v3.tiff"), 
#              fill=c("red","blue","green","orange"), 
#              main = "ChIPseq selected genes", main.fontfamily = "arial", 
#              main.cex = 2, fontfamily = "arial")



## ==== literature ========================================================= ##
lit.genes <- read.csv(file = file.path(dir,"CircadianClock_genelist_literature.csv"), 
                      header=TRUE, stringsAsFactors = FALSE)
lit.genes.ID <- lit.genes$Gene.ID



## ---- Venn Diagrams ------------------------------------------------------ ##
## all GO, ChIPseq and literature genes 
# venn.diagram(x = list(literature=lit.genes.ID,
#                       ChIPseq=c(prr1.targets,prr5.targets,prr7.targets,
#                                 prr7HA.targets),
#                       GO=all.GO.genes),
#              filename = file.path(dir, "genes_GO_ChIPseq_lit_Venn_v3.tiff"),
#              fill=c("red","green","orange"),
#              main = "all pre-filtered selected genes",
#              main.fontfamily = "arial",
#              main.cex = 2, fontfamily = "arial", force.unique = TRUE)





## ==== put ChIPseq + lit + GO together ==================================== ##
# put all the genes together 
prelim.gene.ID.list <- c(selected.GO.genes.df[,"geneID"], 
                         lit.genes.ID, 
                         prr1.targets, prr5.targets, prr7.targets, prr7HA.targets
                         )
prelim.gene.ID.list <- unique(prelim.gene.ID.list)

# remove spaces in the character strings 
prelim.gene.ID.list <- unique(sapply(prelim.gene.ID.list, 
                                     function(x){substring(x,1,9)}, 
                              USE.NAMES = FALSE)
                              )




## ==== calculate Z-score, find unexpressed & unchanging genes ============= ##

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
selected.expr <- selected.expr[rownames(selected.expr) %in% prelim.gene.ID.list,]

# log2(CPM+1) transform of the non-filtered data 
selected.expr.log <- log2(selected.expr + 1)
# rank genes based on sd to choose for most fluctuating ones 
selected.expr.log.sd <- apply(selected.expr.log, 1, sd)

# count the number of <1 counts per gene (unexpressing genes)
n.zero.expr <- apply(selected.expr, 1, function(x){sum(x>5)})

# calculate the CV=sd/mean for each gene 
selected.expr.FC <- apply(selected.expr,1,function(x){diff(range(x))/min(x)})

# calculate the Z-score for each gene and timepoint
selected.expr.Zscore <- (selected.expr - 
                           apply(selected.expr,1,mean))/apply(
                             selected.expr,1,sd)
# count the number of timepoints with NaN per gene (unchanging genes)
n.nan.Zscore <- apply(selected.expr.Zscore, 1, function(x){sum(is.nan(x))})



## ---- FINAL FILTERING ---------------------------------------------------- ##
# NaN Z-score + more than 0 zero counts + low spread
# combination to remove unexpressing or unchanging genes 
selected.expr.filt <- selected.expr[(n.zero.expr > 6) & 
                                           (n.nan.Zscore == 0),]

# log2(CPM+1) transform of the initial filtered data 
selected.expr.filt.log <- log2(selected.expr.filt + 1)
# rank genes based on sd to choose for most fluctuating ones 
selected.expr.filt.log.sd <- apply(selected.expr.filt.log, 1, sd)
selected.expr.filt.log.sd.rank <- rank(selected.expr.filt.log.sd)

# choosing the genes based on the lowest ranks, showing highest sd
n.genes.wanted <- 1000
changing.genes <- names(
  selected.expr.filt.log.sd.rank[
    selected.expr.filt.log.sd.rank > 
      (length(selected.expr.filt.log.sd.rank)-n.genes.wanted)])


# calculate the fold-change between max and min expression for each gene 
selected.expr.filt.FC <- apply(selected.expr.filt[changing.genes,],1,
                               function(x){diff(range(x))/min(x)})

# select expression rows corresponding to the final listed genes 
selected.expr.filt <- selected.expr.filt[changing.genes,]


# selected.expr.filt.NA <- apply(selected.expr.filt, 1, 
#                                function(x){sum(is.na(x))})




# plot the histograms of all the metrics used for filtering
pdf(file.path(dir,"genelist_filtering_stats_180702_v3.pdf"))
hist(n.zero.expr, xlab="# timepoints with < 3 CPM",
     main="Lowly expressed genes")
hist(n.nan.Zscore, xlab="# timepoints with NaN Z-score",
     main="Unexpressed genes")
hist(selected.expr.FC, xlim = c(0,15), breaks=5000,
     xlab="FC (diff(range(CPM))/CPM_min)",
     main="Increase in temporal gene expression (pre-filtered)")
hist(selected.expr.filt.FC, xlim = c(0,15), breaks=5000,
     xlab="FC (diff(range(CPM))/CPM_min)",
     main="Increase in temporal gene expression (filtered)")
hist(log(selected.expr.log.sd), breaks=50,
     xlab="log(standard deviation)",
     main="SD of log transformed CPM log2(CPM+1) PRE-FILTER")
hist(log(selected.expr.filt.log.sd), breaks=50,
     xlab="log(standard deviation)",
     main="SD of log transformed CPM log2(CPM+1) FILTERED")
dev.off()



changing.genes.GOinfo <- changing.genes.GOID <- changing.genes.db <- changing.genes

for (i in 1:length(changing.genes)){
  GO.term <- GO2genes$V5[which(GO2genes$V10 == changing.genes[i])]
  GO.info <- GO2Info$info[GO2Info$GO %in% GO.term]
  GO.process <- GO2Info$db[GO2Info$GO %in% GO.term]
  changing.genes.GOID[i] <- paste(GO.term, collapse = "_")
  changing.genes.GOinfo[i] <- paste(GO.info, collapse = "_")
  changing.genes.db[i] <- paste(GO.process, collapse = "_")
}

changing.genes.df <- cbind(GeneID=changing.genes,
                           GOID=changing.genes.GOID ,
                           GOinfo=changing.genes.GOinfo,
                           GOcategory=changing.genes.db
                           )

TF.conditions <- grepl("transcription factor activity", changing.genes.GOinfo) & 
  grepl("DNA binding", changing.genes.GOinfo) 

TF.list <- changing.genes.df[TF.conditions,1]
target.list <- changing.genes.df[!TF.conditions,1]

write.csv(changing.genes.df, 
          file.path(dir,"changing_genes_GO_list_v3_180702.csv"), 
          quote = FALSE, row.names = FALSE)
write.csv(changing.genes.df[,1], 
          file.path(dir,"changing_genes_list_v3_180702.csv"), 
          quote = FALSE, row.names = FALSE)

write.csv(target.list, 
          file.path(dir,"changing_genes_list_v3_180702_targets.csv"), 
          quote = FALSE, row.names = FALSE)
write.csv(TF.list, 
          file.path(dir,"changing_genes_list_v3_180702_TF.csv"), 
          quote = FALSE, row.names = FALSE)


venn.diagram(x = list(literature=lit.genes.ID,
                      ChIPseq=c(prr1.targets,prr5.targets,prr7.targets,
                                prr7HA.targets),
                      GO=all.GO.genes, 
                      FINAL=changing.genes),
             filename = file.path(dir, "selectedGenes_GO_ChIPseq_lit_Venn_v3.tiff"),
             fill=c("red","green","orange", "blue"),
             main = "Final filtered selected genes",
             main.fontfamily = "arial",
             main.cex = 2, fontfamily = "arial", force.unique = TRUE)

venn.diagram(x = list(ChIPseq=c(prr1.targets,prr5.targets,prr7.targets,
                                prr7HA.targets),
                      GO=all.GO.genes, 
                      TF=TF.list, target=target.list),
             filename = file.path(dir, "selectedTFsTargets_GO_ChIPseq_Venn_v3.tiff"),
             fill=c("red","green","orange", "blue"),
             main = "Final filtered selected TFs & targets",
             main.fontfamily = "arial",
             main.cex = 2, fontfamily = "arial", force.unique = TRUE)


pdf(file.path(dir,"genelist_filtering_180702_v3_expressionPlots.pdf"))
for(i in 1:length(changing.genes)){
  plot(time,selected.expr.filt[i,], 
       main=paste0(rownames(selected.expr.filt[changing.genes[i],]),"\n",
                   selected.expr.filt.log.sd.rank[changing.genes[i]],
                   " (",round(selected.expr.filt.log.sd[changing.genes[i]],4),")"
                   ), 
       xlab="time (min)", ylab="CPM", type='b')
  }
dev.off()









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





