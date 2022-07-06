#! /usr/bin/python
import argparse
'''
Load a reference frames file, isolate a base pair frame in loop, and copy that right beside it.
Recall: For every n bp frame,
# - 1 = bp number and base-pair sequence [5*(n-1)+1]
# - 2 = origin position (x, y, z)        [5*(n-1)+2]
# - 3 = d1 vector                        [5*(n-1)+3]
# - 4 = d2 vector                        [5*(n-1)+4]
# - 5 = d3 vector                        [5*(n-1)+5]
'''
parser = argparse.ArgumentParser()
parser.add_argument("-in", "--input", help="This arguement loads file that contains the list of temp configs.")
parser.add_argument("-out", "--outname", help="The name of the output files")
parser.add_argument("-t", "--target", type=int, help="The base location to be replicated for path increase")
args = parser.parse_args()

infile = open( args.input, "r" )
lst = infile.readlines()
lst = [i.rstrip('\n') for i in lst]
infile.close()

inheader     = lst[0]
loop_len     = int(inheader.split()[0])
out_loop_len = int(loop_len+1)
outheader    = " "+(str(out_loop_len))+" base pairs"

new_loop = lst[ 1 : 5*(args.target)+1] + lst[5*((args.target)-1)+1:5*(args.target)+1] + lst[5*(args.target)+1:]

outfile = open(args.outname, "w")
outfile.write(outheader+'\n')
[outfile.write(i+'\n') for i in new_loop]
outfile.close()

del lst, inheader, loop_len, out_loop_len, outheader, new_loop
