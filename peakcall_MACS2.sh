#!/bin/bash 

# MACS2 peak calling for the ChIPseq files 
# pool all the replicates together - for both control and experiments
# value of -g 1.118e8 for a.thaliana obtained from: https://github.com/iamciera/chipSeqTutorial
# files with replicates were merged (as recommended per MACS manual) by the following way:
# samtools merge out.bam in1.bam in2.bam in3.bam ... 


# PRR1 - done
#macs2 callpeak -t /media/pw_synology3/MattM/SRR411113_trimmed_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bam \
#-c /media/pw_synology3/MattM/SRR411114_trimmed_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bam \
#-f BAM -g 1.118e8 -n /media/pw_synology3/MattM/chipseq_MACS2_PRR1_filt -B -q 0.01

# PRR5 - done 
#macs2 callpeak -t /media/pw_synology3/MattM/SRR442100_trimmed_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bam \
#-c /media/pw_synology3/MattM/SRR442101_trimmed_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bam \
#-f BAM -g 1.118e8 -n /media/pw_synology3/MattM/chipseq_MACS2_PRR5_filt -B -q 0.01

# PRR7
macs2 callpeak -t /media/pw_synology3/MattM/SRR943787_91_trimmed_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard_merged.bam \
-c /media/pw_synology3/MattM/SRR943786_trimmed_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bam -f BAM \
-g 1.118e8 -n /media/pw_synology3/MattM/chipseq_MACS2_PRR7_filt -B -q 0.01

# PRR7::HA-PRR7
macs2 callpeak -t /media/pw_synology3/MattM/SRR943789_90_trimmed_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard_merged.bam \
-c /media/pw_synology3/MattM/SRR943788_trimmed_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bam -f BAM \
-g 1.118e8 -n /media/pw_synology3/MattM/chipseq_MACS2_PRR7HA_filt -B -q 0.01

# PRR9::HA-PRR9 CCR2::LUC
#macs2 callpeak \
#-t /media/pw_synology3/MattM/SRR2131057_trimmed_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bam \
#-c /media/pw_synology3/MattM/SRR2131056_trimmed_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bam \ 
#-f BAM -g 1.118e8 -n /media/pw_synology3/MattM/chipseq_MACS2_PRR9HA_filt -B -q 0.01

# PRR9 CCR2::LUC
#macs2 callpeak \
#-t /media/pw_synology3/MattM/SRR2131059_trimmed_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bam \
#-c /media/pw_synology3/MattM/SRR2131058_trimmed_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bam \ 
#-f BAM -g 1.118e8 -n /media/pw_synology3/MattM/chipseq_MACS2_PRR9_filt -B -q 0.01

