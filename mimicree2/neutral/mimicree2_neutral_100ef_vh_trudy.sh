#!/bin/bash

##evolution of phenotypic variance 
##MimicrEE2
##WYLai 20191101

echo “name”
echo $1
echo “number of contributing loci”
echo $2
echo “repicate”
echo $3
echo “effectsize_distribution”
echo $4
date


mkdir $1_$2
cd $1_$2
mkdir tmp
for	((i=1;i<=$3;i=i+1)); do
	hr=`gshuf -n 1 /Volumes/cluster/Wei-Yun/simulation/heri_trudy_micro.txt`
	Rscript ../pick-SNPs-QTL_100.R -I ../mimihap_effectsize_0.5_100.rds -L $2 -O ./ -R $i -D $4
	java -jar /Volumes/cluster/Wei-Yun/simulation/mim2-v208.jar qff --haplotypes-g0 selectedloci_${2}_$i.mimhap --sex ../sexinfo.txt --recombination-rate ../dsim.rr_X_LDjump-LOESS-0.1.txt --population-size ../population-size_300.txt --effect-size effectsize_${4}_${2}_$i.txt --fitness-function ../fitness_function_tmp.txt --heritability $hr --snapshots 1 --replicate-runs 1 --output-sync tmp.gz --output-gpf tmp.gpf --output-dir tmp/ --threads 2
	Rscript ../generate_fitness_function_neutral.R -I tmp.gpf -O ./ -L $2 -R $i
	mkdir $1_$2_$i
	java -jar /Volumes/cluster/Wei-Yun/simulation/mim2-v208.jar qff --haplotypes-g0 selectedloci_${2}_$i.mimhap --sex ../sexinfo.txt --recombination-rate ../dsim.rr_X_LDjump-LOESS-0.1.txt --population-size ../population-size_300.txt --effect-size effectsize_${4}_${2}_$i.txt --fitness-function fitnessfunction_neutral_$2_$i.txt --heritability $hr --snapshots 1,5,10,15,20,25,35,50,75,100,200 --replicate-runs 1 --output-sync $1_$2_$i.gz --output-gpf $1_$2_$i.gpf --output-dir $1_$2_$i/ --threads 2
done

date
