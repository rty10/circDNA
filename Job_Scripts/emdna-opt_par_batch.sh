#!/bin/bash

declare -a ARRAY=(*.par)

for ((i=0; i<${#ARRAY[@]}; i++));
do
 echo "------------------------------------------------------"
 echo "--- ${ARRAY[$i]%.*} "
 echo "--- start time: $(/bin/date) ---"
 emDNA \
 --x3DNA-bp-step-params-input=${ARRAY[$i]%.*}.par \
 --hold-last-bp \
 --DNA-seqdep-model=IdealDNA \
 --output-name=${ARRAY[$i]%.*}
 
 mv ${ARRAY[$i]%.*}_opt.txt ${ARRAY[$i]%.*}_opt.par
 
 emDNA_topology --x3DNA-bp-step-params-input=${ARRAY[$i]%.*}.par --virtual-last-bp

 echo "--- end time : $(/bin/date) ---" 
 echo "------------------------------------------------------"

 
done
