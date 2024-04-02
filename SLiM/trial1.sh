#!/bin/bash

mkdir /Volumes/Temp1/shengkai/evolVar_pleiotropy/cov$1_prop$2
cd /Volumes/Temp1/shengkai/evolVar_pleiotropy/cov$1_prop$2


d=""
for ((i=1;i<=100;i++));
do
	echo $i
	mkdir simulation_$i
	mkdir simulation_$i/freq
	d="$d dir='/Volumes/Temp1/shengkai/evolVar_pleiotropy/cov$1_prop$2/simulation_$i'"
done

#echo $d
parallel -j 20 slim -d N=20 -d prop=$2 -d cov=$1 -d {} ../correlate_trait_uncorrelate_fitness ::: $d
#parallel -j 20 slim  ::: $d
