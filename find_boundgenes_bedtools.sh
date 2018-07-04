#!/bin/bash 

# find the TF bound genes from the MACS2 peak calls 
# use the GFF file with the appropraite upstream region size included 

# PRR1

for name in 'PRR7' 'PRR7HA'; do
	for type in '1kbp' '500bp' '1kbp500bp' '2kbp' '2kbp500bp'; do 
		bedtools intersect -a gtf_geneProm_"$type".gtf -b \
		/media/pw_synology3/MattM/chipseq_MACS2_"$name"_filt_peaks.narrowPeak \
		-u > "$name"_"$type"_boundgenes.gtf
	done
done  
