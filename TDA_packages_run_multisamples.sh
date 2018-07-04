#!/bin/bash

# script to run the TDAs for different conditions 
# e.g. WT vs mutants and the different temperatures 

for temp in "22" "27"; do 
	for mut in "Col" "Ler" "prr579" "phyA"; do
		#mkdir /home/mahiar.mahjoub/TD_algorithms/genereg_output/"$mut"_"$temp"_180702_v3
		Rscript TDA_packages_GeneReg_Rscript.R $mut $temp "180702_v3"
	done
done 
