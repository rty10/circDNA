#!/bin/bash

declare -a LIST

for FILE in *.dat;
do
 LIST=("${LIST[@]}" "${FILE%.*}")
done

for ITEM in "${LIST[@]}";
do
 python ~/research/Scripts/reconstruction/circle_rf-construct.py ${ITEM}.dat
 
 rm ${ITEM}.dat
 mv ${ITEM}_circ.dat ${ITEM}.dat
 
 emDNA-parser \
 --x3DNA-bp-input=${ITEM}.dat \
 --get-x3DNA-params>${ITEM}.par

 emDNA-topology \
 --x3DNA-bp-step-params-input=${ITEM}.par>topo_${ITEM}.txt \
 --virtual-last-bp

done

