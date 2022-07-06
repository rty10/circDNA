#!/bin/bash

declare -a ARRAY=(*.txt)
for ((i=0; i<${#ARRAY[@]}; i++));
do
 emDNA_topology \
 --bp-list-input=${ARRAY[$i]%.*}.txt>topo_${ARRAY[$i]%.*}.txt \
 --virtual-last-bp
done


