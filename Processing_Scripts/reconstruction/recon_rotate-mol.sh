#!/bin/bash

x3dna_utils cp_std BDNA
NAME=test_conf

for i in {001..100} 
do
 cp store_models_v01/${NAME}-${i}.pdb .

 rotate_mol -c ${NAME}-${i}.pdb ${NAME}-${i}_rm.pdb
 
 mv ${NAME}-${i}_rm.pdb store_models_rotatemol-c
 rm ${NAME}-${i}.pdb
 
done

rm Atomic* rotmat.dat pmiview*
