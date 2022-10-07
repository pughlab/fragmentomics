#!/bin/bash

## Usage splitBedpe.sh input.bedpe

input=$1
name=${1::-6}

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22;
do
echo $chr
grep -w $chr $input > ${name}_${chr}.bedpe
done
