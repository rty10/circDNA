#! /usr/bin/python

import os, sys, subprocess, shutil
path = os.getcwd()

# ----------------------------------------------------------------------------------------------------------------

rfs = []
for filename in os.listdir(path):
    if filename.endswith('.dat'):
        rfs.append( filename.split('.')[0] )

subprocess.call(['x3dna_utils', 'cp_std', 'BDNA'])

for i in range(0, len(sorted(rfs))):
    print(rfs[i])
    subprocess.call("emDNA_parser --x3DNA-bp-input="+rfs[i]+".dat --get-x3DNA-params>"+rfs[i]+".par", shell=True)
    subprocess.call("rebuild -atomic "+rfs[i]+".par "+rfs[i]+".pdb", shell=True)

for filename in os.listdir('.'):
    if 'Atomic' in filename:
        os.remove(filename)
emDNA_files = ['auxiliary.par', 'bestpairs.pdb', 'bp_helical.par', 
                'bp_order.dat', 'bp_step.par', 'cf_7methods.par',
                'col_chains.scr', 'hel_regions.pdb', 'hstacking.pdb', 
                'stacking.pdb','ref_frames.dat', 'col_helices.scr', 
                'poc_haxis.r3d', 'ref_frames.dat', 'ref_frame.dat']
for i in range(0, len(emDNA_files)):
    if emDNA_files[i] in os.listdir('.'):
        os.remove(emDNA_files[i])

