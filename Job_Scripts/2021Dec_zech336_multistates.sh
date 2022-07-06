#!/bin/bash

MODEL="n336"

for ((i=100; i<111; i++));
do 
 emDNA \
 --x3DNA-bp-step-params-input=${MODEL}_${i}.par \
 --hold-last-bp \
 --DNA-seqdep-model=IdealDNA \
 --output-name=${MODEL}_${i}
 mv ${MODEL}_${i}_opt.txt ${MODEL}_${i}_opt.par
 mv ${MODEL}_${i}.log ${MODEL}_${i}_opt.log
 
done
