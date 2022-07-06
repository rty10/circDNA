#! /usr/bin/python

import os, sys, subprocess, shutil
path = os.getcwd()

# -!' Command line: >python /home/rty10/Documents/scripts/reconstruction/pdb_reconstruct.py

path_list = path.split('/')
path_list = path_list[:-1]
path_org  = '/'.join(path_list)
path_org = '/'+path_org

# --- FUNCTIONS ---
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


# ----------------------------------------------------------------------------------------------------------------

pdbs = []
for filename in os.listdir(path):
	if filename.endswith('.pdb'):
		pdbs.append( filename.split('.')[0] )

subprocess.call(['x3dna_utils', 'cp_std', 'BDNA'])

for i in range(0, len(sorted(pdbs))):
	print(pdbs[i])
	subprocess.call("find_pair "+pdbs[i]+".pdb | analyze", shell=True)
	os.rename("bp_step.par", pdbs[i]+".par")
	os.rename("ref_frames.dat", pdbs[i]+".dat")

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

