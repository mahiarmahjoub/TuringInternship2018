#!/bin/bash 

# make sure you are in /home/mahiar.mahjoub


# obtain file containing the SRA IDs of all chipseq files 
file="/home/mahiar.mahjoub/chipseq_SRR_ID.txt"
read -d $'\x04' name < "$file"

echo "Starting running all chipseq files ..." > chipseq_multi.log

for name in $(cat $file); do 
	echo "Running run_chipseq.sh on $name ..." >> chipseq_multi.log
	/home/mahiar.mahjoub/run_chipseq.sh /media/pw_synology3/MattM/"$name"
done
