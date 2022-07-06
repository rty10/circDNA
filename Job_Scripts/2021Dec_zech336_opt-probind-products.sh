#!/bin/bash

MODEL="zech336"
STATE="st33"
POS="bbx"

for ((i=0; i<=328; i++));
do

 SITE2START=$((5+i))
 SITE2END=$((SITE2START+3))
 MODELNAME=${MODEL}-${STATE}_${POS}_pos${SITE2START}
 
 if [ -f "${MODELNAME}.par" ];
 then
  echo "$MODELNAME"
 
  emDNA \
  --x3DNA-bp-step-params-input=${MODELNAME}.par \
  --frozen-steps=1:4,${SITE2START}:${SITE2END} \
  --hold-last-bp \
  --DNA-seqdep-model=IdealDNA \
  --output-name=${MODELNAME}
   
  mv ${MODELNAME}_opt.txt ${MODELNAME}_opt.par
  mv ${MODELNAME}.log ${MODELNAME}_opt.log
  
  emDNA_topology \
  --x3DNA-bp-step-params-input=${MODELNAME}_opt.par>topo_${MODELNAME}_opt.txt
  
  emDNA_parser \
  --x3DNA-bp-step-params-input=${MODELNAME}_opt.par \
  --get-x3DNA-bp>${MODELNAME}_opt.dat
 
 fi
 
done
