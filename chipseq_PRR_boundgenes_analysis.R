## ====================== analyse chipseq PRR targets ====================== ##

setwd(dir = "C://Users/mahia/Google Drive/MPhil/CompBio/Internship/")


#gene.aliases <- read.delim("TAIR_genealiases.txt", header=FALSE, 
stringsAsFactors = FALSE)

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
bg.raw.prr7.1kbp <- read.delim("PRR_chipseq_analysis/PRR5_1kbp_boundgenes.gtf", 
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


# -- 500 bp upstream 
# import bound gene gtf (from bedtools intersect)
bg.raw.prr1.500bp <- read.delim("PRR_chipseq_analysis/PRR1_500bp_boundgenes.gtf",
                               header = FALSE,
                               stringsAsFactors = FALSE)
bg.raw.prr5.500bp <- read.delim("PRR_chipseq_analysis/PRR5_500bp_boundgenes.gtf",
                               header = FALSE,
                               stringsAsFactors = FALSE)
bg.raw.prr7.500bp <- read.delim("PRR_chipseq_analysis/PRR5_500bp_boundgenes.gtf",
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
prr7.500bp.bound.genes <- unique(as.vector(sapply(bg.raw.prr1.500bp$V9, 
                                                  function(x){strsplit(unlist(
                                                    strsplit(x, split = ";"))[2], 
                                                    split=" ")[[1]][3]}
)))


# -- export lists of bound genes for all experiments 



