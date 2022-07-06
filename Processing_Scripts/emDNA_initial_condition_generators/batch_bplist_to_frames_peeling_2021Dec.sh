#!/bin/bash
 
mkdir backup
mkdir bplist_peeled

LENGTH=341
SEQ="A"
declare -a FILENAMES=(*n${LENGTH}*.txt)
for ((k=1; k<="$LENGTH"; k++));
do
SEQ="${SEQ}A"
done

for MODEL in "${FILENAMES[@]}";
do
 python /home/rty10/Documents/scripts/reconstruction/modeling_move_refframes.py \
  -in ${MODEL%.*}.txt \
  -type bp-list \
  -circ yes \
  -range 1:70 \
  -loc upstream \
  -out ${MODEL%.*}_peel

 emDNA_parser \
 --bp-list-input=${MODEL%.*}_peel.txt \
 --bp-list-sequence=$SEQ \
 --get-x3DNA-bp>${MODEL%.*}.dat

 mv ${MODEL%.*}_peel.txt bplist_peeled
 mv ${MODEL%.*}.txt backup
done

LENGTH=359
SEQ="A"
declare -a FILENAMES=(*n${LENGTH}*.txt)
for ((k=1; k<="$LENGTH"; k++));
do
SEQ="${SEQ}A"
done

for MODEL in "${FILENAMES[@]}";
do
 python /home/rty10/Documents/scripts/reconstruction/modeling_move_refframes.py \
  -in ${MODEL%.*}.txt \
  -type bp-list \
  -circ yes \
  -range 1:70 \
  -loc upstream \
  -out ${MODEL%.*}_peel

 emDNA_parser \
 --bp-list-input=${MODEL%.*}_peel.txt \
 --bp-list-sequence=$SEQ \
 --get-x3DNA-bp>${MODEL%.*}.dat

 mv ${MODEL%.*}_peel.txt bplist_peeled
 mv ${MODEL%.*}.txt backup
done

### -- END --
