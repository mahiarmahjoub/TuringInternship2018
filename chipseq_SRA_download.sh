#!/bin/bash 

# download all the required SRA files from the server using the SRAtoolkit
# files downloaded in fastq; automatically converted from sra to fastq


#read name < file_containing _the_answer

#while read LINE; do 
#	echo "Downloading $LINE ..."; 
#	fastq-dump --split-files $LINE
#done < /home/mahiar.mahjoub/chipseq_SRR_ID.txt

#echo "Processing $f now"
#fastq-dump --split-files $f 

file="/home/mahiar.mahjoub/chipseq_SRR_ID.txt"

read -d $'\x04' name < "$file"

#echo $name 

for name in $(cat $file); do 
	echo "Downloading $name ..."
	fastq-dump --split-files $name
done
