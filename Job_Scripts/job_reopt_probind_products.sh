#!/bin/bash

# 2 Dec 2021
# start with parameter files of the emDNA_probind product with zech336 and specific bent/kinked/melted models
# optimize with IdealDNA for "idt105" models

# site1 = 1:5
# site2 = 60:64
# site3 = 159:163

BASENAME="zech336"

declare -A MODELS=( [bbb]=1:5,60:64,159:163 [mmm]=1:5,60:64,159:163 [kkk]=1:5,60:64,159:163 \
[bbx]=1:5,60:64 [mmx]=1:5,60:64 [kkx]=1:5,60:64 \
[bxb]=1:5,159:163 [mxm]=1:5,159:163 [kxk]=1:5,159:163 \
[xbb]=60:64,159:163 [xmm]=60:64,159:163 [xkk]=60:64,159:163 \
[bxx]=1:5 [kxx]=1:5 [mxx]=1:5)

for KEY in "${!MODELS[@]}";

do

# echo "$KEY - ${MODELS[$KEY]}"

 OPTNAME=${BASENAME}_${KEY}_idt105
 echo ${OPTNAME}
 
 emDNA \
 --x3DNA-bp-step-params-input=${BASENAME}_${KEY}.par \
 --frozen-steps=${MODELS[$KEY]} \
 --hold-last-bp \
 --DNA-seqdep-model=IdealDNA \
 --output-name=${OPTNAME}
  
 mv ${OPTNAME}_opt.txt ${OPTNAME}.par
 
 emDNA_topology \
 --x3DNA-bp-step-params-input=${OPTNAME}.par>topo_${OPTNAME}.txt
 
 emDNA_parser \
 --x3DNA-bp-step-params-input=${OPTNAME}.par \
 --get-x3DNA-bp>${OPTNAME}.dat
 
 echo "------------------------------------------------------------"
done



