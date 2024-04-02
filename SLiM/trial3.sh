#!/bin/bash

mkdir /Volumes/Temp1/shengkai/evolVar_pleiotropy/cov$1-$2_prop$3
mkdir /Volumes/Temp1/shengkai/evolVar_pleiotropy/cov$1-$2_prop$3/freq
cd /Volumes/Temp1/shengkai/evolVar_pleiotropy/cov$1-$2_prop$3

slim -d N=1000 -d n=20 -d prop=$3 -d cov1=$1 -d cov2=$2 -d "dir='/Volumes/Temp1/shengkai/evolVar_pleiotropy/cov$1-$2_prop$3'" ../correlate_trait_uncorrelate_fitness_v2

