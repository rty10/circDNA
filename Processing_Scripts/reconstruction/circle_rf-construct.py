#! /usr/bin/python
# Circular constructor

import os, sys
path=os.getcwd()

# Input: reference frame file
# Arguments Input:
#     python SCRIPT.py incon.dat
# Output: new reference frame file where first frame == new last frame

STRUCTURE = sys.argv[1]
STRUCTURE = STRUCTURE.split('.')[0]

infile = open(path+"/"+STRUCTURE+".dat", 'r')
indat  = infile.readlines()
infile.close()

SEQUENCELENGTH = int( indat[0].split()[0] )

if len(indat) == 5*SEQUENCELENGTH+1:
    circle = indat + indat[1:6]

outfile = open(path+"/"+STRUCTURE+"_circ.dat", 'w')
for line in circle:
    outfile.write(line)
outfile.close()
