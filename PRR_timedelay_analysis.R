## ================== analysis of PRR expression dynamics ================== ##

## -- import the raw counts AND gene aliases

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


# import the gene aliases to convert ATG format to gene names 
gene.aliases <- read.delim(file.path(dir,"TAIR_genealiases.txt"), 
                           stringsAsFactors=FALSE
                           )


## -- assign downstream target genes of PRR 
# from (1) GO, (2) ChIPseq (3) CLUST and (4) literature Grundy et al 2015

# upstream regulatory genes of interest
master.regulators <- c("PRR7", "PRR8", "PRR9")
master.regulators.ID <- c("AT5G02810", "AT4G00760", "AT2G46790")

# some of the potential downstream targets
genes.interest <- c("PRR1", "PRR5",
                    "LHY", "CCA1", "ZTL", "GI", "LKP2", "ELF3",
                    "CBF1", "CBF2", "CBF3", "MYC2"
                    )
genes.interest.ID <- gene.aliases$locus_name[gene.aliases$symbol 
                                             %in% genes.interest]

## -- expression of wanted genes AND separate based on condition

time <- c(-10, 10, 24, 30, 45, 60, 105, 120)

# select the desired genes 
expression.all <- cpm.raw[cpm.raw$tracking_id %in% 
                            c(master.regulators.ID, genes.interest.ID),]

## place each separate condition in an object (maybe put in array later)
exp.27.col0 <- expression.all[,
                              colnames(expression.all)[
                                grepl("Col",colnames(expression.all)) & 
                                  grepl("27",colnames(expression.all))
                                ]
                              ]

exp.22.col0 <- expression.all[,
                              colnames(expression.all)[
                                grepl("Col",colnames(expression.all)) & 
                                  grepl("22",colnames(expression.all))
                                ]
                              ]

exp.27.Ler0 <- expression.all[,
                              colnames(expression.all)[
                                grepl("Ler",colnames(expression.all)) & 
                                  grepl("27",colnames(expression.all))
                                ]
                              ]

exp.27.PRR <- expression.all[,
                              colnames(expression.all)[
                                grepl("prr",colnames(expression.all)) & 
                                  grepl("27",colnames(expression.all))
                                ]
                              ]

exp.22.PRR <- expression.all[,
                              colnames(expression.all)[
                                grepl("prr",colnames(expression.all)) & 
                                  grepl("22",colnames(expression.all))
                                ]
                              ]

exp.27.phy <- expression.all[,
                             colnames(expression.all)[
                               grepl("phyAB",colnames(expression.all)) & 
                                 grepl("27",colnames(expression.all))
                               ]
                             ]

exp.22.phy <- expression.all[,
                             colnames(expression.all)[
                               grepl("phyAB",colnames(expression.all)) & 
                                 grepl("22",colnames(expression.all))
                               ]
                             ]

exp.22.hsf <- expression.all[,
                             colnames(expression.all)[
                               grepl("hsfQK",colnames(expression.all)) & 
                                 grepl("22",colnames(expression.all))
                               ]
                             ]


## -- plot results 

par(mfrow=c(2,2))
for (i in 1:4){
  plot(time, exp.22.col0[1,], xlim = c(time[1],time[length(time)]), 
       ylim=c(0,max(as.numeric(expression.all))))
  for (j in 1:nrow(expression.all)){
    
  }
  
  }
