#! /usr/bin/python
# Circular deconstructor for post-optimization using --hold-last-bp functionality

import os, sys
path=os.getcwd()

# Input: circular step_pars file
# Arguments Input:
#     python SCRIPT.py step_pars.par
# Output: new reference frame file where first frame != new last frame

STRUCTURE = sys.argv[1]
STRUCTURE = STRUCTURE.split('.')[0]

infile = open(path+"/"+STRUCTURE+".par", 'r')
indat  = infile.readlines()
infile.close()

CIRCLELENGTH = int( indat[0].split()[0] )
NEWLENGTH    = CIRCLELENGTH - 1
indat[0].replace(str(CIRCLELENGTH), str(NEWLENGTH))

outfile = open(path+"/"+STRUCTURE+"_new.par", 'w')
for line in indat:
    outfile.write(line)
outfile.close()

os.remove(path+"/"+STRUCTURE+".par")
os.rename(path+"/"+STRUCTURE+"_new.par", path+"/"+STRUCTURE+".par")
