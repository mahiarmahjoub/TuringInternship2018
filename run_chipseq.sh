#!/bin/bash 

# chipseq analysis 
# trimming, alignment and normalising chipseq SE reads 
# for the analysis of PRR1/5/7/9 chipseq pulldowns 

filename=$1

# remove adapters and low quality SE reads 
echo "... removing adapters" > "$filename"_run.log
java -jar /home/Program_NGS_sl-pw-srv01/Trimmomatic-0.32/trimmomatic-0.32.jar SE -threads 4 \
-trimlog "$filename"_trimmolog1.txt \
"$filename"_1.fastq "$filename"_trimmed.fastq \
ILLUMINACLIP:/home/Program_NGS_sl-pw-srv01/Trimmomatic-0.32/adapters/TruSeq3-SE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35


# map with bowtie2
echo "... mapping with bowtie2" >> "$filename"_run.log
bowtie2 -x /home/Reference_genomes/Arabidopsis_thaliana_Ensembl_TAIR10/Ensembl/TAIR10/Sequence/Bowtie2Index/genome \
-p 4 -U "$filename"_trimmed.fastq -S "$filename"_trimmed_bowtie2_TAIR10_ensembl_nomixed.sam \
--maxins 800 --no-mixed --no-discordant --no-unal 2>> "$filename"_run.log 


# convert to bam and sort 
echo "... converting sam to bam" >> "$filename"_run.log
samtools view -bS "$filename"_trimmed_bowtie2_TAIR10_ensembl_nomixed.sam \
> "$filename"_trimmed_bowtie2_TAIR10_ensembl_nomixed.bam 2>> "$filename"_run.log

echo "... sorting mapped reads" >> "$filename"_run.log
samtools sort -o "$filename"_trimmed_bowtie2_TAIR10_ensembl_nomixed_sorted.bam \
"$filename"_trimmed_bowtie2_TAIR10_ensembl_nomixed.bam 2>> "$filename"_run.log 


# mark and remove duplicates 
echo "...removing duplicates" >> "$filename"_run.log
java -Xmx4g -jar /home/Program_NGS_sl-pw-srv01/picard-tools-1.103/MarkDuplicates.jar \
INPUT="$filename"_trimmed_bowtie2_TAIR10_ensembl_nomixed_sorted.bam \
OUTPUT="$filename"_trimmed_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bam \
METRICS_FILE=dup.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true 2> markup_stderr.txt


# indexing 
echo "... samtools indexing bam" >> "$filename"_run.log
samtools index "$filename"_trimmed_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bam 


# word count raw data 
raw_line=$(wc -l < "$filename"_1.fastq)
raw_count=$(echo "$raw_line/4" | bc -l)
echo "number of raw reads: $raw_count" >> "$filename"_run.log


# get flagstat 
echo "number of clean reads mapped:" >> "$filename"_run.log
samtools flagstat "$filename"_trimmed_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bam >> "$filename"_run.log


# estimate genome average and normalise 
echo "... normalising reads" >> "$filename"_run.log
genomeCoverageBed -split -bg -ibam "$filename"_trimmed_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bam \
-g /home/Reference_genomes/Arabidopsis_thaliana_Ensembl_TAIR10/Ensembl/TAIR10/Annotation/Genes/ChromInfo.txt \
> "$filename"_trimmed_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bedgraph


sum=$(samtools depth "$filename"_trimmed_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bam | awk '{sum+=$3;cnt++}END{printf "%.0f", sum}')
sum_norm=$(echo "$sum/119667750" | bc -l)
echo "genome normalised coverage: $sum_norm" >> "$filename"_run.log

export MYVAR=$sum_norm
perl -e 'print $ENV{MYVAR}."\n"'


# normalise read counts by genome-wide coverage 
echo "... normalising read counts by genome-wide coverage" >> "$filename"_run.log
perl -ne 'chomp($_); @a=split(/\t/,$_);print $a[0]."\t".$a[1]."\t".$a[2]."\t".$a[3]/$ENV{MYVAR}."\n";' \
"$filename"_trimmed_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bedgraph \
> "$filename"_trimmed_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard_genomenorm.bedgraph


# convert bedgraph to bigwig
echo "... converting to bedgrap to BigWig" >> "$filename"_run.log
bedGraphToBigWig "$filename"_trimmed_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard_genomenorm.bedgraph \
/home/Reference_genomes/Arabidopsis_thaliana_Ensembl_TAIR10/Ensembl/TAIR10/Annotation/Genes/ChromInfo.txt \
"$filename"_trimmed_bowtie2_TAIR10_ensembl_nomixed_sorted_rmdup_picard_genomenorm.bw
