#! /usr/bin/python
import argparse
'''
Load a bp reference frames list file, isolate a base pair frame in loop, and copy that right beside it.
Recall: Each line is a full reference frame {{ox, oy, ox},{{d1x, d1y, d1z},{d2x, d2y, d2z},{d3x, d3y, d3z}}}
'''

# cli example: python structure_length_change_bplist.py -l 336 -in a1_n336.txt -ti 141 -td 336 -out a1_n337.txt
parser = argparse.ArgumentParser()
parser.add_argument("-l",   "--length",  dest="length", type=int,               required=True, help="Size of structure.")
parser.add_argument("-in",  "--input",   dest="input",                          required=True, help="This arguement loads file that contains the list of temp configs.")
parser.add_argument("-out", "--outname", dest="output",                         required=True, help="The name of the output files")
parser.add_argument("-ti",   "--repeat-target",  dest="target_inc", type=int,   required=True, help="The base location to be replicated for fragment increase")
parser.add_argument("-td",   "--remove-target",  dest="target_dec", type=int,   required=True, help="The base location to be remove for fragment decrease")
args = parser.parse_args()
# -----------------------------------------------------------------------------------------

# ---CODE--- 

infile = open( args.input, "r" )
lst = infile.readlines()
lst = [i.rstrip('\n') for i in lst]
infile.close()

# Increase fragement first
# recall: targets are base pair number, not index number
inc_idx = (args.target_inc)-1
dec_idx = (args.target_dec)-1
new_loop = lst[: inc_idx+1] + lst[inc_idx:inc_idx+1] + lst[inc_idx+1:]
# Decrease other fragment
if dec_idx > inc_idx:
    new_loop = new_loop[: (1+dec_idx)] + new_loop[(1+dec_idx)+1:]
else:
    new_loop = new_loop[: dec_idx] + new_loop[dec_idx+1:]
del inc_idx, dec_idx


outfile = open(args.output, "w")
[outfile.write(i+'\n') for i in new_loop]
outfile.write('\n')
outfile.close()

del lst, new_loop
# -----------------------------------------------------------------------------------------
