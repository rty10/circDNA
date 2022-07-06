#!/bin/bash

LENGTH=336
LOOPINTS=$((LENGTH-(140+140+1)))
MODEL="nuc2_n${LENGTH}"

NUC1="1:140"
NUC2A=142
NUC2B=281

cp ${MODEL}.txt ${MODEL}_slide0.txt

#for ((i=0; i<=${LOOPINTS}; i++));
for ((i=0; i<=11; i++));
do
 NAME=${MODEL}_slide${i}
 
 emDNA \
 --bp-list-input=${NAME}.txt \
 --frozen-steps=${NUC1},$((NUC2A+i)):$((NUC2B+i)) \
 --hold-last-bp \
 --energy-progress \
 --DNA-seqdep-model=IdealDNA \
 --output-name=${NAME}
 mv ${NAME}.log ${NAME}_opt.log
 
 #if [ $i -lt $LOOPINTS ]
 if [ $i -lt 11 ]
 then
  python /home/rty10/Documents/scripts/optimization/dinuc_sliding_bplist.py \
  -in ${NAME}_opt.txt \
  -l ${LENGTH} \
  -ti 141 \
  -td ${LENGTH} \
  -out ${MODEL}_slide$((i+1)).txt
 fi
 
done
