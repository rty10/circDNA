#! /usr/bin/python

import os, sys, subprocess, shutil
path = os.getcwd()

# start from central optimized directory
# -!' Command line: >python /home/rty10/Documents/scripts/reconstruction/pdb+out_reconstruct.py

path_list = path.split('/')
path_list = path_list[:-1]
path_org  = '/'.join(path_list)
path_org = '/'+path_org

# ---------------------------------------------------------------------------------------------------------------

pars = []
for filename in os.listdir(path):
    if filename.endswith('.par'):
        pars.append( filename.split('.')[0] )


subprocess.call(['x3dna_utils', 'cp_std', 'BDNA'])

for i in range(0, len(sorted(pars))):   
    subprocess.call("rebuild -atomic "+pars[i]+".par "+pars[i]+".pdb", shell=True)
    
    subprocess.call("find_pair "+pars[i]+".pdb | analyze", shell=True)
    
    #os.remove(pars[i]+".pdb")

for filename in os.listdir('.'):
    if 'Atomic' in filename:
        os.remove(filename)
emDNA_files = ['auxiliary.par', 'bestpairs.pdb', 'bp_helical.par', 
                'bp_order.dat', 'bp_step.par', 'cf_7methods.par',
                'col_chains.scr', 'hel_regions.pdb', 'hstacking.pdb', 
                'stacking.pdb','ref_frames.dat', 'col_helices.scr', 
                'poc_haxis.r3d', 'ref_frames.dat']
for i in range(0, len(emDNA_files)):
    if emDNA_files[i] in os.listdir('.'):
        os.remove(emDNA_files[i])

os.chdir(os.pardir)
