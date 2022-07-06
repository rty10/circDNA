#!/bin/bash

declare -a ARRAY=(*.txt)

for ((i=0; i<${#ARRAY[@]}; i++));
do
 LENGTH=$(wc -l < "${ARRAY[$i]%.*}.txt")
 for ((k=1; k<="$LENGTH"; k++));
 do
  SEQ="${SEQ}A"
 done
 
 #echo ${ARRAY[$i]%.*}
 emDNA_parser \
 --bp-list-input=${ARRAY[$i]%.*}.txt \
 --bp-list-sequence=$SEQ \
 --get-x3DNA-params>${ARRAY[$i]%.*}.par
 
 unset LENGTH
 unset SEQ
 
done

