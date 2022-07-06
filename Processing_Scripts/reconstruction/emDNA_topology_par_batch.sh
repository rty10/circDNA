#!/bin/bash

declare -a ARRAY=(*.par)
for ((i=0; i<${#ARRAY[@]}; i++));
do
 emDNA_topology \
 --x3DNA-bp-step-params-input=${ARRAY[$i]%.*}.par>topo_${ARRAY[$i]%.*}.txt \
 --virtual-last-bp
done


