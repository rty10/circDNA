#!/bin/bash

declare -a ARRAY=(*.dat)

for ((i=0; i<${#ARRAY[@]}; i++));
do
 #echo ${ARRAY[$i]%.*}
 emDNA_parser \
 --x3DNA-bp-input=${ARRAY[$i]%.*}.dat \
 --get-x3DNA-params>${ARRAY[$i]%.*}.par
 
done

