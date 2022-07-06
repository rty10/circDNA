#!/bin/bash

x3dna_utils cp_std BDNA
NAME=test_conf

for i in {001..100} 
do
 cp store_models_v01/${NAME}-${i}.pdb .
 find_pair ${NAME}-${i}.pdb | analyze
 cp ref_frames.dat ${NAME}-${i}.dat
 frame_mol -121 ${NAME}-${i}.dat ${NAME}-${i}.pdb ${NAME}-${i}_fm.pdb
 rm ${NAME}-${i}.out \
 ${NAME}-${i}.dat \
 auxiliary.par \
 bestpairs.pdb \
 bp_helical.par \
 bp_order.dat \
 bp_step.par \
 cf_7methods.par \
 col_chains* \
 col_helices* \
 hel_regions.pdb \
 hstacking.pdb \
 stacking.pdb \
 
 mv ${NAME}-${i}_fm.pdb store_models_framemol
 rm ${NAME}-${i}.pdb
done
rm Atomic*
rm rotmat.dat
rm ref_frames.dat