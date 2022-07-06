#! /usr/bin/python

import os, sys, argparse

path=os.getcwd()

# cli example: python modeling_move_refframes.py -in <FILENAME> -type {bp-list, x3dna-bp-list} -circ {yes, no} -range x:m -loc x -out <OUTNAME>
parser = argparse.ArgumentParser()
parser.add_argument("-in",  "--input", dest="input",required=True, help="This arguement loads file that contains the list of temp configs.")

parser.add_argument("-type","--format-type", choices=["bp-list", "x3dna-bp-list"], dest="type", required=True, help="Argument selects reference frame file type.")
parser.add_argument("-circ","--circular", choices=["yes", "no"], dest="circ", required=True, help="Argument stating whether the reference frames are for a circle.")

parser.add_argument("-range","--frames-range", dest="range", required=True, help="Argument for the inclusive slice of base pair frames that will be moved.")
parser.add_argument("-loc",  "--location", choices=["upstream", "downstream"], dest="loc", required=True, help="The argument that states which base pair to move the chunk of frames BEFORE.")

parser.add_argument("-out", "--outname", dest="output", required=True, help="The name of the output files")

args = parser.parse_args()
# -----------------------------------------------------------------------------------------

FILE_TYPES = {"bp-list":"txt", "x3dna-bp-list":"dat"}

STRUCTURE = args.input.split('.')[0]
EXT = FILE_TYPES[args.type]
infile = open(path+"/"+STRUCTURE+"."+EXT, 'r')
indat  = infile.readlines()
indat  = [i.rstrip('\n') for i in indat]
infile.close()


if args.circ=="yes":
    if args.type=="bp-list":
        indat = indat[:-1]
    elif args.type=="x3dna-bp-list":
        HEADER = indat[0:1]
        indat = indat[1:-5]

RANGESTART, RANGEEND = (args.range).split(':')
RANGESTART=int(RANGESTART)
RANGEEND=int(RANGEEND)

if args.type=="bp-list":
    if args.loc == "upstream":
        new_indat = indat[RANGEEND:] + indat[RANGESTART-1:RANGEEND]
    elif args.loc == "downstream":
        new_indat = indat[RANGESTART-1:RANGEEND] + indat[:RANGESTART-1] 
    
elif args.type=="x3dna-bp-list":
    new_indat = indat[:-5]

if args.circ=="yes":
    if args.type=="bp-list":
        new_indat = new_indat + new_indat[0:1]
    elif args.type=="x3dna-bp-list":
        new_indat = HEADER + new_indat + new_indat[0:5]
        del HEADER

NEWSTRUCTURE = args.output
outfile = open(path+"/"+NEWSTRUCTURE+"."+EXT, 'w')
for line in new_indat:
    outfile.write(line+"\n")
outfile.close()

del indat, new_indat
