## ====================== analyse chipseq PRR targets ====================== ##

setwd(dir = "C://Users/mahia/Google Drive/MPhil/CompBio/Internship/")


#gene.aliases <- read.delim("TAIR_genealiases.txt", header=FALSE, 
#stringsAsFactors = FALSE)

gtf.raw <- read.delim("PRR_chipseq_analysis/genes.gtf", header=FALSE, 
                      stringsAsFactors = FALSE)
gtf.genes <- unique(as.vector(sapply(gtf.raw$V9, 
                                     function(x){strsplit(
                                       unlist(strsplit(x, split = ";"))[2], 
                                       split=" ")[[1]][3]}
)))
  


# -- 1kbp upstream 
# import bound gene gtf (from bedtools intersect)
bg.raw.prr1.1kbp <- read.delim("PRR_chipseq_analysis/PRR1_1kbp_boundgenes.gtf",
                              header = FALSE,
                              stringsAsFactors = FALSE)
bg.raw.prr5.1kbp <- read.delim("PRR_chipseq_analysis/PRR5_1kbp_boundgenes.gtf", 
                               header = FALSE,
                               stringsAsFactors = FALSE)
bg.raw.prr7.1kbp <- read.delim("PRR_chipseq_analysis/PRR7_1kbp_boundgenes.gtf", 
                               header = FALSE,
                               stringsAsFactors = FALSE)
bg.raw.prr7ha.1kbp <- read.delim("PRR_chipseq_analysis/PRR7HA_1kbp_boundgenes.gtf", 
                               header = FALSE,
                               stringsAsFactors = FALSE)
# extract gene IDs
prr5.1kbp.bound.genes <- unique(as.vector(sapply(bg.raw.prr5.1kbp$V9, 
                           function(x){
                             strsplit(unlist(strsplit(x, split = ";"))[2], 
                                      split=" ")[[1]][3]}
                           )))

prr1.1kbp.bound.genes <- unique(as.vector(sapply(bg.raw.prr1.1kbp$V9, 
                            function(x){strsplit(unlist(
                              strsplit(x, split = ";"))[2], 
                                                 split=" ")[[1]][3]}
                            )))

prr7.1kbp.bound.genes <- unique(as.vector(sapply(bg.raw.prr7.1kbp$V9, 
                                                 function(x){strsplit(unlist(
                                                   strsplit(x, split = ";"))[2], 
                                                   split=" ")[[1]][3]}
)))
prr7ha.1kbp.bound.genes <- unique(as.vector(sapply(bg.raw.prr7ha.1kbp$V9, 
                                                 function(x){strsplit(unlist(
                                                   strsplit(x, split = ";"))[2], 
                                                   split=" ")[[1]][3]}
)))


# -- 500 bp upstream 
# import bound gene gtf (from bedtools intersect)
bg.raw.prr1.500bp <- read.delim("PRR_chipseq_analysis/chipseq_boundgenes_040618/PRR1_500bp_boundgenes.gtf",
                               header = FALSE,
                               stringsAsFactors = FALSE)
bg.raw.prr5.500bp <- read.delim("PRR_chipseq_analysis/chipseq_boundgenes_180618/PRR5_500bp_FCge10_boundgenes.gtf",
                               header = FALSE,
                               stringsAsFactors = FALSE)
bg.raw.prr7.500bp <- read.delim("PRR_chipseq_analysis/chipseq_boundgenes_040618/PRR7_500bp_boundgenes.gtf",
                                header = FALSE,
                                stringsAsFactors = FALSE)
bg.raw.prr7ha.500bp <- read.delim("PRR_chipseq_analysis/chipseq_boundgenes_180618/PRR7HA_500bp_FCge5_boundgenes.gtf",
                                header = FALSE,
                                stringsAsFactors = FALSE)
# extract gene IDs
prr5.500bp.bound.genes <- unique(as.vector(sapply(bg.raw.prr5.500bp$V9, 
                                                 function(x){strsplit(unlist(
                                                   strsplit(x, split = ";"))[2],
                                                   split=" ")[[1]][3]}
)))

prr1.500bp.bound.genes <- unique(as.vector(sapply(bg.raw.prr1.500bp$V9, 
                                                 function(x){strsplit(unlist(
                                                   strsplit(x, split = ";"))[2], 
                                                   split=" ")[[1]][3]}
)))
prr7.500bp.bound.genes <- unique(as.vector(sapply(bg.raw.prr7.500bp$V9, 
                                                  function(x){strsplit(unlist(
                                                    strsplit(x, split = ";"))[2], 
                                                    split=" ")[[1]][3]}
)))
prr7ha.500bp.bound.genes <- unique(as.vector(sapply(bg.raw.prr7ha.500bp$V9, 
                                                  function(x){strsplit(unlist(
                                                    strsplit(x, split = ";"))[2], 
                                                    split=" ")[[1]][3]}
)))


# -- analyse common (intersect) genes 
library(VennDiagram)
venn.diagram(x = list(prr5=prr5.500bp.bound.genes, 
                      prr1=prr1.500bp.bound.genes, 
                      prr7=prr7.500bp.bound.genes,
                      prr7HA=prr7ha.500bp.bound.genes), 
             filename = "PRR_chipseq_analysis/chipseq_Prr1_5_7_7HA_venn_500bp_filt.tiff", 
             fill = c('red','blue','orange','green')
             )
venn.diagram(x = list(prr7=prr7.500bp.bound.genes,
                      prr7HA=prr7ha.500bp.bound.genes), 
             filename = "PRR_chipseq_analysis/chipseq_Prr7_7HA_venn_500bp_filt.tiff", 
             fill = c('red','blue')
)
venn.diagram(x = list(prr5=prr5.1kbp.bound.genes, 
                      prr1=prr1.1kbp.bound.genes, 
                      prr7=prr7.1kbp.bound.genes,
                      prr7HA=prr7ha.1kbp.bound.genes), 
             filename = "PRR_chipseq_analysis/chipseq_Prr1_5_7_7HA_venn_1kbp.tiff", 
             fill = c('red','blue','orange','green')
)
venn.diagram(x = list(prr7=prr7.1kbp.bound.genes,
                      prr7HA=prr7ha.1kbp.bound.genes), 
             filename = "PRR_chipseq_analysis/chipseq_Prr7_7HA_venn_1kbp.tiff", 
             fill = c('red','blue')
)
venn.diagram(x = list(prr7.1kbp=prr7.1kbp.bound.genes,
                      prr7.500bp=prr7.500bp.bound.genes), 
             filename = "PRR_chipseq_analysis/chipseq_Prr7_1kbpvs500bp_venn_500bp.tiff", 
             fill = c('red','blue')
)


# -- collate and export lists of bound genes for all experiments 
ngenes.max <- length(prr5.500bp.bound.genes)
boundgenes.list.500bp <- as.data.frame(matrix(NA,ngenes.max,4))
boundgenes.list.1kbp <- as.data.frame(matrix(NA,ngenes.max,4))
colnames(boundgenes.list.1kbp) <- 
  colnames(boundgenes.list.500bp) <- 
  c("PRR1", "PRR5", "PRR7", "PRR7HA")

boundgenes.list.500bp$PRR1 <- c(prr1.500bp.bound.genes, rep(NA,ngenes.max-length(prr1.500bp.bound.genes)))
boundgenes.list.500bp$PRR5 <- c(prr5.500bp.bound.genes, rep(NA,ngenes.max-length(prr5.500bp.bound.genes)))
boundgenes.list.500bp$PRR7 <- c(prr7.500bp.bound.genes, rep(NA,ngenes.max-length(prr7.500bp.bound.genes)))
boundgenes.list.500bp$PRR7HA <- c(prr7ha.500bp.bound.genes, rep(NA,ngenes.max-length(prr7ha.500bp.bound.genes)))
boundgenes.list.1kbp$PRR1 <- c(prr1.1kbp.bound.genes, rep(NA,ngenes.max-length(prr1.1kbp.bound.genes)))
boundgenes.list.1kbp$PRR5 <- c(prr5.1kbp.bound.genes, rep(NA,ngenes.max-length(prr5.1kbp.bound.genes)))
boundgenes.list.1kbp$PRR7 <- c(prr7.1kbp.bound.genes, rep(NA,ngenes.max-length(prr7.1kbp.bound.genes)))
boundgenes.list.1kbp$PRR7HA <- c(prr7ha.1kbp.bound.genes, rep(NA,ngenes.max-length(prr7ha.1kbp.bound.genes)))

write.csv(x = boundgenes.list.1kbp, 
          file = "PRR_chipseq_analysis/boundgenes_1kbp.csv", 
          quote = FALSE, row.names = FALSE)
write.csv(x = boundgenes.list.500bp, 
          file = "PRR_chipseq_analysis/boundgenes_500bp_MACS2filt.csv", 
          quote = FALSE, row.names = FALSE)


