#!/bin/bash

EMDNAPATH="/home/rty10/emDNA-install/bin/Release"


declare -a MODELS=(IdealDNA Olson1998 Young2019)


for MODEL in ${MODELS[@]};
do

 ${EMDNAPATH}/emDNA-ff-packager \
 --force-constants-input=ForceConstants_${MODEL}_tet.txt \
 --intrinsic-steps-input=StepParameters_${MODEL}_tet.txt \
 --model-name=${MODEL}-tet

done



