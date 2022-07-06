#! /usr/bin/python

import os, sys, subprocess, shutil
path = os.getcwd()

# -!' Command line: >python /home/rty10/Documents/scripts/reconstruction/rf_reconstruct.py

def par_format(new_par, old_par):
    outfile = open(new_par, 'w')
    infile  = open(old_par, 'r')
    indata  = infile.readlines()
    infile.close()
    
    CIRCLELENGTH = int( indata[0].split()[0] )
    indata[0] = indata[0].replace(str(CIRCLELENGTH), str(CIRCLELENGTH-1))
    
    for line in indata:
        outfile.write(line)
    outfile.close()
    os.remove(path+'/'+old_par)
    os.rename(path+'/'+new_par, path+'/'+old_par)
    return


refframes = []
for filename in os.listdir('.'):
    if filename.endswith('.dat'):
        refframes.append( filename.split('.')[0] )

subprocess.call(['x3dna_utils', 'cp_std', 'BDNA'])

for i in range(0, len(sorted(refframes))):
    print(refframes[i])
    
    subprocess.call("emDNA_parser --x3DNA-bp-input="+refframes[i]+".dat --get-x3DNA-params>"+refframes[i]+".par", shell=True)
    
    par_format("nuc.par", refframes[i]+".par")
    
    #subprocess.call("rebuild -atomic "+refframes[i]+".par "+refframes[i]+".pdb", shell=True)
    
    #os.remove(refframes[i]+".dat")


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

