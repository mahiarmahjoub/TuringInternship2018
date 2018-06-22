## ================== analysis of PRR expression dynamics ================== ##

## ---- import the raw counts AND gene aliases ----------------------------- ##

# counts in normalised in CPM based on library size 
# contains selected timepoints based on NitPicker
# lowly expressed genes excluded 
dir <- "C:/Users/mahia/Google Drive/MPhil/CompBio/Internship"
cpm.raw <- read.table(file = file.path(dir, "dawnburst_slcu", 
                                       "rawcounts", 
                                       "dawnburst_masterdf.txt"
), header = TRUE, 
stringsAsFactors = FALSE
)
gene.alias <- read.delim(file.path(dir,"TAIR_genealiases.txt"), header = TRUE, 
                         stringsAsFactors = FALSE)




## ---- assign downstream target genes of PRR ------------------------------ ##
# from (1) GO, (2) ChIPseq (3) CLUST and (4) literature 

# genes of intereset from LITERATURE
# includes master regulators + downstream targets
genes.interest <- read.csv(file = file.path(dir,
                                            "CircadianClock_genelist_literature.csv"), 
                           header = TRUE, stringsAsFactors = FALSE)

# time points in minutes 
time <- c(-10, 10, 24, 30, 45, 60, 105, 120)




## ---- visualisation of gene expression ----------------------------------- ##
# consulted: https://rstudio-pubs-static.s3.amazonaws.com/123938_b1ce6ecfb4c342fca62cc7d4703b8dcd.html
require(gplots, quietly = TRUE)
require(RColorBrewer, quietly = TRUE)

## ==== literature based list ============================================== ##
# select the desired genes from the expression matrix 
ratio.genes.present <- sum(cpm.raw$tracking_id %in% genes.interest$Gene.ID)/nrow(genes.interest)
expression.selected.genes <- cpm.raw[cpm.raw$tracking_id %in% 
                                       genes.interest$Gene.ID,]
genes.selected <- expression.selected.genes$tracking_id
genes.selected.name <- genes.interest$Gene.Name[genes.interest$Gene.ID %in% genes.selected]
col.fun <- rep('black',length(genes.selected))
col.fun[genes.interest$Role == "flowering"] <- 'green'
col.fun[genes.interest$Role == "CC"] <- 'blue'
col.fun[genes.interest$Role == "light"] <- 'orange'
col.fun[genes.interest$Role == "temp"] <- 'red'
col.fun[genes.interest$Role %in% c("biosynthesis","ethylene response")] <- 'purple'

## -- choose condition of interest 
# Col Ler prr579 phyABcry12 hsfQK (22-27)
mutant <- "Col"
temp <- "22"
# extract relevant columns 
exp.cols <- grepl(mutant, colnames(expression.selected.genes)) & 
  grepl(temp, colnames(expression.selected.genes))
exp <- as.matrix(expression.selected.genes[,exp.cols])
# log normalise 
exp.log <- as.matrix(log2(exp + 1)) 
# order genes based on the order of time each gene reaches peak expression
exp.temporal.max.sorted <- order(apply(exp.log,1,which.max))
# log mean normalise 
exp.log.norm <- exp.log - apply(exp.log,1,mean)

# plot - log + unscaled
#pdf(file = file.path(dir,paste0("heatmap_",mutant,"_",temp,".pdf")), width = 7, height = 5)
heatmap.2(t(exp.log[exp.temporal.max.sorted,]), 
          Rowv = NA, Colv=NA, trace = 'none', scale = 'none',
          #col=colorRampPalette(colors = c("black","white"))(20),
          col=(brewer.pal(9,"OrRd")), dendrogram = 'none',
          labCol = genes.selected.name[exp.temporal.max.sorted], 
          labRow = as.character(time), key.title = NA, key.ylab = NA, 
          key.xlab = "log2(CPM+1)", ylab = "time (mins)", cexCol = 0.55, 
          colCol = col.fun[exp.temporal.max.sorted],
          srtCol = 45, main = paste0(mutant," ",temp," ","nonscaled log(CPM)")
)
# plot - column (gene) scaled + non-log
heatmap.2(t(exp[exp.temporal.max.sorted,]), 
          #Rowv = NA, 
          Colv=NA, trace = 'none', scale = 'col',
          #col=colorRampPalette(colors = c("black","white"))(20),
          col=(brewer.pal(4,"OrRd")), 
          labCol = genes.selected.name[exp.temporal.max.sorted], 
          labRow = as.character(time), key.title = NA, key.ylab = NA,
          ylab = "time (mins)", cexCol = 0.55, srtCol = 45, 
          colCol = col.fun[exp.temporal.max.sorted],
          main = paste0(mutant," ",temp," ","col.scaled CPM")
)
# plot - column (gene) scaled + log
heatmap.2(t(exp.log[exp.temporal.max.sorted,]), 
          Rowv = NA, Colv=NA, trace = 'none', scale = 'col',
          #col=colorRampPalette(colors = c("black","white"))(20),
          col=(brewer.pal(4,"OrRd")), dendrogram = 'none',
          labCol = genes.selected.name[exp.temporal.max.sorted], 
          labRow = as.character(time), key.title = NA, key.ylab = NA,
          ylab = "time (mins)",  cexCol = 0.55, srtCol = 45, 
          colCol = col.fun[exp.temporal.max.sorted],
          main = paste0(mutant," ",temp," ","col.scaled log(CPM)")
)
#dev.off()




## ==== ChIPseq based list ================================================= ##
# select the desired genes from the expression matrix 
chipseq.protein <- "PRR7"
genes.interest.chipseq <- read.csv(file.path(dir,"PRR_chipseq_analysis/boundgenes_500bp.csv"), 
                                   header = TRUE, stringsAsFactors = FALSE)
expression.selected.genes <- cpm.raw[cpm.raw$tracking_id %in% 
                                       genes.interest.chipseq[,chipseq.protein],]
genes.selected <- genes.selected.name <- expression.selected.genes$tracking_id
# genes.selected.name[genes.selected 
#                     %in% 
#                       gene.alias$locus_name] <- gene.alias$symbol[gene.alias$locus_name %in% 
#                                                                     genes.selected 
#                                                                     ]


## -- choose condition of interest 
# Col Ler prr579 phyABcry12 hsfQK (22-27)
mutant <- "Ler"
temp <- "22"
# extract relevant columns 
exp.cols <- grepl(mutant, colnames(expression.selected.genes)) & 
  grepl(temp, colnames(expression.selected.genes))
exp <- as.matrix(expression.selected.genes[,exp.cols])
# log normalise 
exp.log <- as.matrix(log2(exp + 1)) 
# order genes based on the order of time each gene reaches peak expression
exp.temporal.max.sorted <- sort(apply(exp.log,1,which.max), index.return=TRUE)$ix
# log mean normalise 
exp.log.norm <- exp.log - apply(exp.log,1,mean)

# plot - log + unscaled
pdf(file = file.path(dir,paste0("heatmap_chipseq_",chipseq.protein,
                                "_",mutant,"_",temp,".pdf")), width = 7, 
    height = 5)
heatmap.2(t(exp.log[exp.temporal.max.sorted,]), 
          Rowv = NA, Colv=NA, trace = 'none', scale = 'none',
          #col=colorRampPalette(colors = c("black","white"))(20),
          col=(brewer.pal(9,"OrRd")), dendrogram = 'none',
          #labCol = genes.selected[exp.temporal.max.sorted], 
          labCol = NA,
          labRow = as.character(time), key.title = NA, key.ylab = NA, 
          key.xlab = "log2(CPM+1)", ylab = "time (mins)", cexCol = 0.55, 
          srtCol = 45, main = paste0(mutant," ",temp," ",chipseq.protein,
                                     " ChIPseq \n nonscaled log(CPM)")
)
# plot - column (gene) scaled + non-log
heatmap.2(t(exp[exp.temporal.max.sorted,]), 
          Rowv = NA, Colv=NA, trace = 'none', scale = 'col',
          #col=colorRampPalette(colors = c("black","white"))(20),
          col=(brewer.pal(4,"OrRd")), dendrogram = 'none',
          #labCol = genes.selected[exp.temporal.max.sorted], 
          labCol = NA,
          labRow = as.character(time), key.title = NA, key.ylab = NA,
          ylab = "time (mins)", cexCol = 0.55, srtCol = 45, 
          main = paste0(mutant," ",temp," ",chipseq.protein,
                        " ChIPseq \n col.scaled CPM")
)
# plot - column (gene) scaled + log
heatmap.2(t(exp.log[exp.temporal.max.sorted,]), 
          Rowv = NA, Colv=NA, trace = 'none', scale = 'col',
          #col=colorRampPalette(colors = c("black","white"))(20),
          col=(brewer.pal(4,"OrRd")), dendrogram = 'none',
          #labCol = genes.selected[exp.temporal.max.sorted], 
          labCol = NA,
          labRow = as.character(time), key.title = NA, key.ylab = NA,
          ylab = "time (mins)",  cexCol = 0.55, srtCol = 45, 
          main = paste0(mutant," ",temp," ",chipseq.protein,
                        " ChIPseq \n col.scaled log(CPM)")
)
dev.off()







