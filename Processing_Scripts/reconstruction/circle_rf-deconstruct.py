#! /usr/bin/python
# Circular deconstructor
import os, sys
path=os.getcwd()

# Input: circular reference frame file
# Arguments Input:
#     python SCRIPT.py incon.dat
# Output: new reference frame file where first frame =/= new last frame

STRUCTURE = sys.argv[1]
STRUCTURE = STRUCTURE.split('.')[0]

infile = open(path+"/"+STRUCTURE+".dat", 'r')
indat  = infile.readlines()
uncircle = indat[:-5]
infile.close()

SEQUENCELENGTH = int( indat[0].split()[0] )
indat[0].replace( str(SEQUENCELENGTH), str(SEQUENCELENGTH - 1) )

outfile = open(path+"/"+STRUCTURE+"_new.dat", 'w')
for line in circle:
    outfile.write(line)
outfile.close()

os.remove(path+"/"+STRUCTURE+".dat")
os.rename(path+"/"+STRUCTURE+"_new.dat", path+"/"+STRUCTURE+".dat")