## ======================= MACS2 narrowPeak analysis ======================= ##

dir <- "C://Users/mahia/Google Drive/MPhil/CompBio/Internship/PRR_chipseq_analysis"

# narrowPeak headings 
col.labs <- c("chr", "start", "end", "source", "score", ".", 
              "FoldChange", "logP", "logQ", "pos")

# import narrowPeak files 
prr1.peakcall <- read.delim(file = file.path(dir, 
                                             "chipseq_MACS2_PRR1_filt_peaks.narrowPeak"), 
                            header = FALSE, stringsAsFactors = FALSE, 
                            col.names = col.labs)
prr5.peakcall <- read.delim(file = file.path(dir, 
                                             "chipseq_MACS2_PRR5_filt_peaks.narrowPeak"), 
                            header = FALSE, stringsAsFactors = FALSE, 
                            col.names = col.labs)
prr7.peakcall <- read.delim(file = file.path(dir, 
                                             "chipseq_MACS2_PRR7_filt_peaks.narrowPeak"), 
                            header = FALSE, stringsAsFactors = FALSE, 
                            col.names = col.labs)
prr7ha.peakcall <- read.delim(file = file.path(dir, 
                                             "chipseq_MACS2_PRR7HA_filt_peaks.narrowPeak"), 
                            header = FALSE, stringsAsFactors = FALSE, 
                            col.names = col.labs)

## == distributions 
pdf(file = file.path(dir,"peakcall_stats_PRR1577HA.pdf"))
par(mfrow=c(2,2))
## fold change distributions 
hist(prr1.peakcall$FoldChange, breaks=30, xlab="fold change", main="PRR1", xlim = c(0,20))
hist(prr5.peakcall$FoldChange, breaks=30, xlab="fold change", main="PRR5", xlim = c(0,20))
hist(prr7.peakcall$FoldChange, breaks=30, xlab="fold change", main="PRR7", xlim = c(0,20))
hist(prr7ha.peakcall$FoldChange, breaks=30, xlab="fold change", main="PRR7-HA", xlim = c(0,20))

## score distributions 
hist(prr1.peakcall$score, breaks=40, xlab="score", main="PRR1", xlim = c(0,900))
hist(prr5.peakcall$score, breaks=40, xlab="score", main="PRR5", xlim = c(0,900))
hist(prr7.peakcall$score, breaks=40, xlab="score", main="PRR7", xlim = c(0,900))
hist(prr7ha.peakcall$score, breaks=40, xlab="score", main="PRR7-HA", xlim = c(0,900))

## -log(FDR) distributions 
hist(prr1.peakcall$logQ, breaks=50, xlab="-log(FDR)", main="PRR1")
hist(prr5.peakcall$logQ, breaks=50, xlab="-log(FDR)", main="PRR5")
hist(prr7.peakcall$logQ, breaks=50, xlab="-log(FDR)", main="PRR7")
hist(prr7ha.peakcall$logQ, breaks=50, xlab="-log(FDR)", main="PRR7-HA")
dev.off()

## == export files after filtering with new cutoffs 
# following options for threshold can be used: 
# (1) FC > 5 and/or (2) -log(q) >= 4 
# let's choose a threshold of FC > 5
write.table(prr5.peakcall[prr5.peakcall$FoldChange >= 10,], 
            file = file.path(dir,"chipseq_MACS2_PRR5_filt_peaks_filter_FCge10.narrowPeak"), 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t"
            )
write.table(prr7ha.peakcall[prr7ha.peakcall$FoldChange >= 5,], 
            file = file.path(dir,"chipseq_MACS2_PRR7HA_filt_peaks_filter_FCge5.narrowPeak"), 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t"
)
write.table(prr5.peakcall[prr5.peakcall$logQ >= 4,], 
            file = file.path(dir,"chipseq_MACS2_PRR5_filt_peaks_filter_q0.0001.narrowPeak"), 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t"
)
write.table(prr7ha.peakcall[prr7ha.peakcall$logQ >= 4,], 
            file = file.path(dir,"chipseq_MACS2_PRR7HA_filt_peaks_filter_q0.0001.narrowPeak"), 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t"
)


