#!/bin/bash

awk 'BEGIN{FS=OFS="\t"}{ if($7=="+") { start = $4 - 1000; end = $5 + 500 }
\
else { start = $4 - 500; end=$5 + 1000 }; \
print $1,$2,$3,start,end,$6,$7,$8,$9}' \
genes.gtf > gtf_geneProm_1kbp500bp.gtf


#/local/data/public/genome_informatics_2017/programs/bedtools2/bin/bedtools \
#intersect -a genePromoter.gtf \
#-b PAX5_macs2.narrowPeak -u > PAX5_bound.genes.gtf
