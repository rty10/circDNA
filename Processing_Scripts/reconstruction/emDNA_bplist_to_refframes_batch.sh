#!/bin/bash

declare -a ARRAY=(*.txt)

for ((i=0; i<${#ARRAY[@]}; i++));
do

 echo "" >> ${ARRAY[$i]%.*}.txt
 
 LENGTH=$(wc -l < "${ARRAY[$i]%.*}.txt")
 for ((k=1; k<"$LENGTH"; k++));
 do
  SEQ="${SEQ}A"
 done
 
 #echo ${ARRAY[$i]%.*}
 emDNA_parser \
 --bp-list-input=${ARRAY[$i]%.*}.txt \
 --bp-list-sequence=$SEQ \
 --get-x3DNA-bp>${ARRAY[$i]%.*}.dat
 
 unset LENGTH
 unset SEQ
 
done

