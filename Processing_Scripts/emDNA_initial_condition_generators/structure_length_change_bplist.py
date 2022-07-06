#! /usr/bin/python
import argparse
'''
Load a bp reference frames list file, isolate a base pair frame in loop, and copy that right beside it.
Recall: Each line is a full reference frame {{ox, oy, ox},{{d1x, d1y, d1z},{d2x, d2y, d2z},{d3x, d3y, d3z}}}
'''

# cli example: python structure_length_change_bplist.py -c inc -in a1_n336.txt -t 300 -out a1_n337.txt
parser = argparse.ArgumentParser()
parser.add_argument("-c",   "--change",  choices=["inc", "dec"], dest="change", required=True, help="Argument to increase (inc) or decrease (dec) configuration size.")
parser.add_argument("-in",  "--input",   dest="input",                          required=True, help="This arguement loads file that contains the list of temp configs.")
parser.add_argument("-out", "--outname", dest="output",                         required=True, help="The name of the output files")
parser.add_argument("-t",   "--target",  dest="target", type=int,               required=True, help="The base location to be replicated for path increase")
args = parser.parse_args()
# -----------------------------------------------------------------------------------------

# ---CODE--- 

infile = open( args.input, "r" )
lst = infile.readlines()
lst = [i.rstrip('\n') for i in lst]
infile.close()
if args.change == "inc":
    new_loop = lst[: args.target+1] + lst[args.target:args.target+1] + lst[args.target+1:]
elif args.change == "dec":
    new_loop = lst[: args.target] + lst[args.target+1:]
  
outfile = open(args.output, "w")
[outfile.write(i+'\n') for i in new_loop]
outfile.close()
del lst, new_loop
# -----------------------------------------------------------------------------------------
