## ===================== configure SRR download links ===================== ##

library(stringi)
dir <- "C://Users/mahia/Google Drive/MPhil/CompBio/Internship"

# csv containing the SRR and GEO IDs of all ChIPseq samples 
filenames <- read.csv(file = file.path(dir,"PRR_chipseq_files.csv"), 
                      header = TRUE, stringsAsFactors = FALSE)

# extract the first six letters of the SRR ID (using stringi)
srr.six <- sapply(filenames$SRA.Accession, function(x){stri_sub(x,1,6)})

# template download link (first bit)
url.template <- "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR"

# put SRR six digit, SRR ID and .sra at the end of the link
urls <- paste0(file.path(url.template,srr.six,filenames$SRA.Accession,
                         filenames$SRA.Accession),".sra")

# export the file containing all the download links for use by 'wget -i'
write.table(urls, file.path(dir,"chipseq_SRR_download_links.txt"), 
            quote = FALSE, row.names = FALSE, col.names = FALSE
            )
write.table(filenames$SRA.Accession, file.path(dir,"chipseq_SRR_ID.txt"), 
            quote = FALSE, row.names = FALSE, col.names = FALSE
)
