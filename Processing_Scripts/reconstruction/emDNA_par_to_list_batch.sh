#!/bin/bash

declare -a ARRAY=(*.par)

for ((i=0; i<${#ARRAY[@]}; i++));
do
 #echo ${ARRAY[$i]%.*}
 emDNA_parser \
 --x3DNA-bp-step-params-input=${ARRAY[$i]%.*}.par \
 --get-bp-list>${ARRAY[$i]%.*}.txt
 
done

