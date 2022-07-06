#! /usr/bin/python

import os, sys, subprocess, shutil
path = os.getcwd()

# -!' Command line: >python /home/rty10/Documents/scripts/reconstruction/rf_reconstruct.py

def par_format(new_par, old_par):
    with open(new_par, 'w') as new_file:
        with open(old_par, 'r') as f:
            lines = f.readlines()
            N_seq = int(lines[0].split()[0])
            seq2  = str(N_seq - 1)
            for line in lines:
                if line.startswith(seq2 + " # ***local"):
                    new_file.write(line.replace(line.split('#')[0], '  0 '))
                else:
                    new_file.write(line)
    new_file.close()
    f.close()
    subprocess.call(['rm', old_par])
    subprocess.call(['mv', new_par, old_par])
    return


refframes = []
for filename in os.listdir('.'):
    if filename.endswith('.dat'):
        refframes.append( filename.split('.')[0] )

subprocess.call(['x3dna_utils', 'cp_std', 'BDNA'])

for i in range(0, len(sorted(refframes))):
    print refframes[i]
    
    subprocess.call("emDNA_parser --x3DNA-bp-input="+refframes[i]+".dat --get-x3DNA-params>"+refframes[i]+".par", shell=True)
    
    par_format("nuc.par", refframes[i]+".par")
    
    subprocess.call("rebuild -atomic "+refframes[i]+".par "+refframes[i]+".pdb", shell=True)
    
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
        
        data_dirs = ['par_files','pdb_files','output_files']
for i in range(0, len(data_dirs)):
    if not os.path.exists(data_dirs[i]):
        os.mkdir(data_dirs[i])
for filename in os.listdir('.'):
    if filename.endswith('.pdb'):
        shutil.move(filename, 'pdb_files')
    elif filename.endswith('.par'):
        shutil.move(filename, 'par_files')
    elif filename.endswith('.out'):
        shutil.move(filename, 'output_files')


