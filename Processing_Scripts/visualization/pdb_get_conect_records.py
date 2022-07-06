#! /usr/bin/python

import os, sys, argparse
path = os.getcwd()
parser = argparse.ArgumentParser()
parser.add_argument("-in", "--input", help="PDB file that has CONECT records")
parser.add_argument("-out", "--output", help="TXT file that contains CONECT records")
args = parser.parse_args()

infile1=open(args.input, "r")
data_conect = infile1.readlines()
infile1.close()
data_conect = [i for i in data_conect if 'CONECT ' in i]


COMBINE = data_conect + ['END']

infile2=open(args.output, "w")
[infile2.write(line) for line in COMBINE]
infile2.close()
del data_conect, COMBINE

